"""
HiTE 工具函数模块|HiTE Utility Functions Module

提供Singularity容器管理和通用工具函数
Provides Singularity container management and utility functions
"""

import os
import subprocess
import shutil
import re
from typing import List, Union, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name"""
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment"""
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, command] + args
    else:
        return [command] + args


class SingularityContainerManager:
    """
    Singularity 容器管理器|Singularity Container Manager

    管理Singularity容器的运行和挂载
    Manages Singularity container execution and mounting
    """

    def __init__(self, config, logger):
        """
        初始化容器管理器|Initialize container manager

        Args:
            config: 配置对象(HiteConfig或PanHiteConfig)|Configuration object
            logger: 日志器对象|Logger object
        """
        self.config = config
        self.logger = logger
        self.singularity_cmd = config.singularity_cmd
        self.sif_file = config.sif_file

    def build_mount_options(self, config) -> List[str]:
        """
        构建Singularity挂载选项|Build Singularity mount options

        根据配置对象类型(HiteConfig或PanHiteConfig)构建相应的挂载选项
        Build mount options based on configuration type (HiteConfig or PanHiteConfig)

        注意：不挂载输出目录，避免HiTE清理时冲突
        Note: Don't mount output directory to avoid conflicts during HiTE cleanup

        Args:
            config: 配置对象|Configuration object

        Returns:
            List[str]: 挂载选项列表|List of mount options
        """
        mounts = []

        # 不挂载输出目录，让容器内处理输出，最后再复制出来
        # Don't mount output directory, let container handle it internally

        # HiTE 单基因组: 挂载基因组目录|HiTE single-genome: mount genome directory
        if hasattr(config, 'genome'):
            genome_dir = os.path.dirname(config.genome)
            mounts.extend(['-B', f'{genome_dir}:{genome_dir}'])

        # panHiTE: 挂载泛基因组相关目录|panHiTE: mount pan-genome related directories
        if hasattr(config, 'pan_genomes_dir'):
            pan_dir = config.pan_genomes_dir
            mounts.extend(['-B', f'{pan_dir}:{pan_dir}'])

            # 挂载基因组列表文件所在目录|Mount genome list file directory
            list_dir = os.path.dirname(config.genome_list)
            mounts.extend(['-B', f'{list_dir}:{list_dir}'])

            # 挂载基因注释目录|Mount genes annotation directory
            if config.genes_dir:
                mounts.extend(['-B', f'{config.genes_dir}:{config.genes_dir}'])

            # 挂载RNA-seq目录|Mount RNA-seq directory
            if config.rna_dir:
                mounts.extend(['-B', f'{config.rna_dir}:{config.rna_dir}'])

        # 挂载用户主目录|Mount user home directory
        if config.mount_home:
            home = os.path.expanduser("~")
            mounts.extend(['-B', f'{home}:{home}'])

        return mounts

    def run_container_command(self, config, internal_cmd: List[str],
                             cmd_type: str = "HiTE") -> bool:
        """
        运行容器内的命令|Run command inside container

        Args:
            config: 配置对象|Configuration object
            internal_cmd: 容器内执行的命令列表|Command list to execute inside container
            cmd_type: 命令类型描述|Command type description (default: "HiTE")

        Returns:
            bool: 执行成功返回True|Returns True if execution succeeds

        Raises:
            subprocess.CalledProcessError: 如果命令执行失败|If command execution fails
        """
        # 构建Singularity命令|Build Singularity command
        singularity_cmd_name = os.path.basename(self.singularity_cmd)
        singularity_args = ['run'] + self.build_mount_options(config) + [self.sif_file] + internal_cmd

        # 使用conda wrapper包装singularity命令|Use conda wrapper to wrap singularity command
        cmd = build_conda_command(singularity_cmd_name, singularity_args)

        self.logger.info(
            f"执行{cmd_type}命令|Executing {cmd_type} command: {' '.join(cmd)}"
        )

        # 执行命令并实时捕获输出|Execute command and capture output in real-time
        try:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )

            # 实时输出日志|Output logs in real-time
            for line in process.stdout:
                line = line.strip()
                if line:  # 跳过空行|Skip empty lines
                    self.logger.info(line)

            # 等待进程结束|Wait for process to complete
            process.wait()

            # 检查返回码|Check return code
            if process.returncode != 0:
                raise subprocess.CalledProcessError(
                    process.returncode, cmd
                )

            self.logger.info(
                f"{cmd_type}命令执行成功|{cmd_type} command executed successfully"
            )
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"{cmd_type}命令执行失败|{cmd_type} command execution failed: "
                f"返回码|returncode: {e.returncode}"
            )
            raise
        except Exception as e:
            self.logger.error(
                f"{cmd_type}命令执行出错|{cmd_type} command execution error: {str(e)}"
            )
            raise

    def copy_results_from_container(self, container_output_dir: str,
                                    host_output_dir: str) -> bool:
        """
        从容器复制结果文件到主机输出目录
        Copy result files from container to host output directory

        通过挂载输出目录的父目录到特殊挂载点，实现文件复制
        Copy files by mounting the parent directory of output to a special mount point

        Args:
            container_output_dir: 容器内输出目录路径|Container output directory path
            host_output_dir: 主机输出目录路径|Host output directory path

        Returns:
            bool: 复制成功返回True|Returns True if copy succeeds
        """
        import pathlib

        self.logger.info(f"容器输出目录|Container output dir: {container_output_dir}")
        self.logger.info(f"主机输出目录|Host output dir: {host_output_dir}")

        # 确保主机输出目录存在
        pathlib.Path(host_output_dir).mkdir(parents=True, exist_ok=True)

        # 获取输出目录的父目录|Get parent directory of output
        host_parent = os.path.dirname(host_output_dir)
        output_basename = os.path.basename(host_output_dir)

        # 挂载父目录到特殊挂载点|Mount parent directory to special mount point
        mount_point = '/host_output_mount'

        copy_cmd = [
            self.singularity_cmd,
            'exec',
            '-B', f'{host_parent}:{mount_point}',
            self.sif_file,
            'bash',
            '-c',
            f'''
            if [ -d "{container_output_dir}" ]; then
                cp -r {container_output_dir}/* {mount_point}/{output_basename}/ 2>/dev/null || true
                echo "Copy completed"
            else
                echo "Warning: Container output directory not found"
                exit 1
            fi
            '''
        ]

        self.logger.info(f"执行复制命令|Executing copy command...")
        self.logger.debug(f"复制命令|Copy command: {' '.join(copy_cmd[:6])} ...")

        try:
            result = subprocess.run(
                copy_cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5分钟超时|5 minute timeout
            )

            # 输出复制命令的结果|Output copy command result
            if result.stdout:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        self.logger.info(f"  {line}")

            if result.returncode != 0:
                self.logger.warning(f"复制结果时出现警告|Warning during copy: {result.stderr}")
            else:
                self.logger.info("复制命令执行成功|Copy command executed successfully")

            # 列出复制的文件|List copied files
            output_files = list(pathlib.Path(host_output_dir).iterdir())
            if output_files:
                self.logger.info(
                    f"成功复制|Successfully copied {len(output_files)} 个文件|files:"
                )
                for f in list(output_files)[:10]:
                    self.logger.info(f"  - {f.name}")
                if len(output_files) > 10:
                    self.logger.info(
                        f"  ... 还有 {len(output_files) - 10} 个文件|and {len(output_files) - 10} more files"
                    )
            else:
                self.logger.warning("输出目录为空|Output directory is empty")

            return True

        except subprocess.TimeoutExpired:
            self.logger.error("复制命令超时|Copy command timeout")
            return False
        except Exception as e:
            self.logger.error(f"复制结果时出错|Error during copy: {str(e)}")
            return False


def find_singularity_executable(custom_path: str = None) -> str:
    """
    查找Singularity可执行文件|Find Singularity executable

    Args:
        custom_path: 自定义Singularity路径|Custom Singularity path (optional)

    Returns:
        str: Singularity可执行文件的绝对路径|Absolute path to Singularity executable

    Raises:
        RuntimeError: 如果找不到Singularity|If Singularity is not found
    """
    # 如果指定了自定义路径|If custom path is specified
    if custom_path:
        expanded = os.path.expanduser(custom_path)
        if os.path.exists(expanded):
            return expanded

    # 候选路径列表|Candidate path list
    candidates = [
        '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
        'singularity',
        '/usr/bin/singularity',
        '/usr/local/bin/singularity'
    ]

    # 尝试候选路径|Try candidate paths
    for path in candidates:
        expanded = os.path.expanduser(path)
        if os.path.exists(expanded):
            return expanded

    raise RuntimeError(
        "Singularity未找到|Singularity not found. "
        "请安装Singularity或指定正确路径|Please install Singularity or specify correct path."
    )


def find_hite_sif_file(custom_path: str = None) -> str:
    """
    查找HiTE SIF文件|Find HiTE SIF file

    Args:
        custom_path: 自定义SIF文件路径|Custom SIF file path (optional)

    Returns:
        str: HiTE SIF文件的绝对路径|Absolute path to HiTE SIF file

    Raises:
        FileNotFoundError: 如果找不到SIF文件|If SIF file is not found
    """
    # 如果指定了自定义路径|If custom path is specified
    if custom_path and os.path.exists(custom_path):
        return os.path.abspath(custom_path)

    # 默认路径|Default path
    default_path = '~/software/singularity/hite_3.3.3.sif'
    if os.path.exists(default_path):
        return default_path

    # 检查环境变量|Check environment variable
    env_path = os.getenv('HITE_SIF')
    if env_path and os.path.exists(env_path):
        return env_path

    raise FileNotFoundError(
        "HiTE SIF文件未找到|HiTE SIF file not found. "
        "请下载SIF文件或指定正确路径|Please download SIF file or specify correct path."
    )


def validate_file_exists(file_path: str, file_desc: str = "文件|File") -> None:
    """
    验证文件是否存在|Validate if file exists

    Args:
        file_path: 文件路径|File path
        file_desc: 文件描述|File description (default: "文件|File")

    Raises:
        FileNotFoundError: 如果文件不存在|If file does not exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            f"{file_desc}不存在|{file_desc} not found: {file_path}"
        )


def validate_directory_exists(dir_path: str, dir_desc: str = "目录|Directory") -> None:
    """
    验证目录是否存在|Validate if directory exists

    Args:
        dir_path: 目录路径|Directory path
        dir_desc: 目录描述|Directory description (default: "目录|Directory")

    Raises:
        FileNotFoundError: 如果目录不存在|If directory does not exist
    """
    if not os.path.exists(dir_path):
        raise FileNotFoundError(
            f"{dir_desc}不存在|{dir_desc} not found: {dir_path}"
        )
