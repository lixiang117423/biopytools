"""
PanMAN工具函数模块|PanMAN Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, Tuple


class PanMANLogger:
    """PanMAN日志管理器|PanMAN Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "panman_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除旧日志文件|Delete old log file if exists
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class PanMANCommandRunner:
    """PanMAN命令执行器|PanMAN Command Runner"""

    def __init__(self, logger, conda_env: str, backend: str = "conda", threads: int = 8,
                 conda_base: str = None, sif_image: str = None, singularity_path: str = None):
        self.logger = logger
        self.conda_env = conda_env
        self.backend = backend
        self.threads = threads
        # 支持自定义conda路径或从环境变量读取
        import os as os_module
        self.conda_base = conda_base or os_module.environ.get(
            'CONDA_BASE',
            os_module.path.join(os_module.environ.get('HOME', ''), 'miniforge3')
        )
        # Singularity 相关参数|Singularity related parameters
        self.sif_image = sif_image
        self.singularity_path = singularity_path

    def run_panman_command(self, args: list, description: str = "") -> Tuple[bool, str]:
        """
        执行panmanUtils命令|Execute panmanUtils command

        Args:
            args: 命令参数列表|Command argument list
            description: 命令描述|Command description

        Returns:
            (成功状态, 输出内容)|(Success status, Output content)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 构建完整命令|Build full command
        if self.backend == "conda":
            cmd = self._build_conda_command(args)
        elif self.backend == "singularity":
            cmd = self._build_singularity_command(args)
        else:  # docker
            cmd = self._build_docker_command(args)

        # 添加线程数参数|Add threads parameter (总是添加以确保线程参数生效|Always add to ensure thread parameter takes effect)
        cmd.extend(["--threads", str(self.threads)])

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            # 使用实时输出避免缓冲区死锁（Singularity产生大量输出）|Use real-time output to avoid buffer deadlock (Singularity generates large output)
            # 对于conda后端，需要在输入文件目录运行以避免工作目录问题|For conda backend, run in input file dir to avoid working directory issues
            # 同时使用相对路径参数|Also use relative path arguments
            cwd = None
            modified_args = args.copy()

            if self.backend == "conda" and len(args) >= 2:
                # 找到输入文件路径，切换到该目录，并使用相对路径|Find input file path, change to that dir, use relative paths
                input_file = None
                input_dir = None

                for i, arg in enumerate(args):
                    if arg in ['-I', '--input-panman'] and i + 1 < len(args):
                        # Extract命令的输入文件|Input file for extract command
                        input_file = args[i + 1]
                        input_dir = os.path.dirname(input_file)
                        if input_dir and os.path.isdir(input_dir):
                            modified_args[i + 1] = os.path.basename(input_file)
                        break
                    elif arg in ['-P', '--input-pangraph'] and i + 1 < len(args):
                        input_file = args[i + 1]
                        input_dir = os.path.dirname(input_file)
                        if input_dir and os.path.isdir(input_dir):
                            # 改为相对路径|Change to relative path
                            modified_args[i + 1] = os.path.basename(input_file)
                        break
                    elif arg in ['-G', '--input-gfa'] and i + 1 < len(args):
                        input_file = args[i + 1]
                        input_dir = os.path.dirname(input_file)
                        if input_dir and os.path.isdir(input_dir):
                            modified_args[i + 1] = os.path.basename(input_file)
                        break
                    elif arg in ['-M', '--input-msa'] and i + 1 < len(args):
                        input_file = args[i + 1]
                        input_dir = os.path.dirname(input_file)
                        if input_dir and os.path.isdir(input_dir):
                            modified_args[i + 1] = os.path.basename(input_file)
                        break

                if input_dir and os.path.isdir(input_dir):
                    cwd = input_dir
                    # Newick文件也改为相对路径|Also change Newick to relative path
                    for i, arg in enumerate(modified_args):
                        if arg in ['-N', '--input-newick'] and i + 1 < len(modified_args):
                            newick_file = modified_args[i + 1]
                            if os.path.isabs(newick_file):
                                modified_args[i + 1] = os.path.basename(newick_file)
                            break

                    # 输出路径改为相对路径|Change output path to relative
                    for i, arg in enumerate(modified_args):
                        if arg in ['-o', '--output-file'] and i + 1 < len(modified_args):
                            output_path = modified_args[i + 1]
                            if os.path.isabs(output_path):
                                modified_args[i + 1] = os.path.basename(output_path)
                            break

                    self.logger.info(f"设置工作目录|Setting working directory: {cwd}")
                    self.logger.debug(f"修改后的参数|Modified args: {modified_args}")

            # 重新构建命令|Rebuild command with modified args
            final_cmd = [cmd[0]]  # 可执行文件
            final_cmd.extend(modified_args)

            result = subprocess.run(
                final_cmd,
                capture_output=False,  # 实时输出到终端|Output directly to terminal
                check=True,
                env=os.environ.copy(),
                cwd=cwd  # 设置工作目录|Set working directory
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return (True, "")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            error_msg = f"Return code: {e.returncode}"
            return (False, error_msg)
        except Exception as e:
            self.logger.error(f"执行出错|Execution error: {str(e)}")
            return (False, str(e))

    def _build_conda_command(self, args: list) -> list:
        """构建Conda命令|Build Conda command"""
        # 直接调用conda环境中的可执行文件，避免conda run的问题
        # Directly call executable in conda env to avoid conda run issues
        panman_utils_path = os.path.join(
            self.conda_base,
            "envs",
            self.conda_env,
            "bin",
            "panmanUtils"
        )
        cmd = [panman_utils_path]
        cmd.extend(args)
        return cmd

    def _build_docker_command(self, args: list) -> list:
        """构建Docker命令|Build Docker command"""
        cmd = [
            "docker", "run",
            "-i",  # 交互模式
            "--rm",  # 容器退出后自动删除
            "swalia14/panman:latest",
            "panmanUtils"
        ]
        cmd.extend(args)
        return cmd

    def _build_singularity_command(self, args: list) -> list:
        """构建Singularity命令|Build Singularity command"""
        cmd = [
            self.singularity_path,
            "exec",
            self.sif_image,
            "panmanUtils"
        ]
        cmd.extend(args)
        return cmd

    def check_conda_env(self) -> bool:
        """检查Conda环境是否存在|Check if Conda environment exists"""
        try:
            cmd = [
                os.path.join(self.conda_base, "bin", "conda"),
                "env", "list"
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            envs = result.stdout
            return self.conda_env in envs
        except Exception as e:
            self.logger.warning(f"无法检查Conda环境|Cannot check Conda environment: {e}")
            return False


class PanMANValidator:
    """PanMAN文件验证器|PanMAN File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def validate_json_file(self, json_file: str) -> bool:
        """验证JSON文件格式|Validate JSON file format"""
        try:
            import json
            with open(json_file, 'r', encoding='utf-8') as f:
                json.load(f)
            self.logger.debug(f"JSON文件格式有效|JSON file format valid: {json_file}")
            return True
        except json.JSONDecodeError as e:
            self.logger.error(f"无效的JSON文件|Invalid JSON file: {json_file} - {e}")
            return False
        except Exception as e:
            self.logger.error(f"读取JSON文件失败|Failed to read JSON file: {e}")
            return False

    def validate_fasta_file(self, fasta_file: str) -> bool:
        """验证FASTA文件格式|Validate FASTA file format"""
        try:
            with open(fasta_file, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"无效的FASTA文件|Invalid FASTA file: {fasta_file}")
                    return False
            self.logger.debug(f"FASTA文件格式有效|FASTA file format valid: {fasta_file}")
            return True
        except Exception as e:
            self.logger.error(f"读取FASTA文件失败|Failed to read FASTA file: {e}")
            return False

    def validate_newick_file(self, newick_file: str) -> bool:
        """验证Newick文件格式|Validate Newick file format"""
        try:
            with open(newick_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
                # Newick文件应以分号结尾|Newick file should end with semicolon
                if not content.endswith(';'):
                    self.logger.warning(
                        f"Newick文件可能格式不正确（缺少分号）|"
                        f"Newick file may be malformed (missing semicolon): {newick_file}"
                    )
                    return False
            self.logger.debug(f"Newick文件格式有效|Newick file format valid: {newick_file}")
            return True
        except Exception as e:
            self.logger.error(f"读取Newick文件失败|Failed to read Newick file: {e}")
            return False

    def validate_gfa_file(self, gfa_file: str) -> bool:
        """验证GFA文件格式|Validate GFA file format"""
        try:
            with open(gfa_file, 'r', encoding='utf-8') as f:
                # GFA文件应该有H或S行|GFA file should have H or S lines
                has_valid_lines = False
                for line in f:
                    if line.startswith('H\t') or line.startswith('S\t'):
                        has_valid_lines = True
                        break

                if not has_valid_lines:
                    self.logger.error(f"无效的GFA文件|Invalid GFA file: {gfa_file}")
                    return False

            self.logger.debug(f"GFA文件格式有效|GFA file format valid: {gfa_file}")
            return True
        except Exception as e:
            self.logger.error(f"读取GFA文件失败|Failed to read GFA file: {e}")
            return False


class PanMANParser:
    """PanMAN结果解析器|PanMAN Result Parser"""

    def __init__(self, logger):
        self.logger = logger

    def parse_summary_output(self, output: str) -> dict:
        """解析摘要输出|Parse summary output"""
        summary = {}
        for line in output.split('\n'):
            if '|' in line:
                key, value = line.split('|', 1)
                summary[key.strip()] = value.strip()
        return summary

    def format_file_size(self, size_bytes: int) -> str:
        """格式化文件大小|Format file size"""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.2f}{unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.2f}TB"

    def get_file_size(self, file_path: str) -> str:
        """获取文件大小|Get file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            return self.format_file_size(size_bytes)
        except Exception as e:
            self.logger.warning(f"无法获取文件大小|Cannot get file size: {e}")
            return "Unknown"
