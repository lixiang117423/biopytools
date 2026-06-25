"""
MEME运行器模块|MEME Runner Module
"""

import subprocess
import os
import shutil
import re
from pathlib import Path
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
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
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class MemeRunner:
    """MEME运行器|MEME Runner"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def run(self) -> bool:
        """
        运行MEME命令|Run MEME command

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("="*60)
        self.logger.info("运行MEME|Running MEME")
        self.logger.info("="*60)

        # 构建命令|Build command
        cmd = self._build_command()

        # 自动包装conda环境的命令|Auto-wrap conda environment commands
        if cmd:
            cmd_name = os.path.basename(cmd[0])
            wrapped_cmd = build_conda_command(cmd_name, cmd[1:])
        else:
            wrapped_cmd = cmd

        self.logger.info(f"命令|Command: {' '.join(wrapped_cmd)}")

        # 处理输出目录|Handle output directory
        output_dir = self.config.get_meme_output_dir()

        # 检查输出是否已存在|Check if output already exists
        xml_file = Path(output_dir) / "meme.xml"
        if xml_file.exists():
            self.logger.info(f"MEME输出已存在，跳过运行|MEME output already exists, skipping run: {xml_file}")
            return True

        # 如果输出目录已存在（但没有XML文件），删除它
        if Path(output_dir).exists():
            self.logger.info(f"删除不完整的输出目录|Removing incomplete output directory: {output_dir}")
            import shutil
            shutil.rmtree(output_dir)

        # 运行命令|Run command
        try:
            self.logger.info("开始运行MEME...|Starting MEME...")

            # 使用subprocess运行，不捕获输出以显示进度
            # 在当前目录运行，MEME会自己创建输出目录
            result = subprocess.run(
                wrapped_cmd,
                check=True
            )

            self.logger.info("MEME运行完成|MEME completed successfully")

            # 检查输出文件|Check output files
            xml_file = Path(output_dir) / "meme.xml"
            if xml_file.exists():
                self.logger.info(f"输出文件已生成|Output file generated: {xml_file}")
                return True
            else:
                self.logger.warning(f"MEME运行完成但未找到输出文件|MEME completed but output file not found: {xml_file}")
                return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"MEME运行失败|MEME run failed with exit code {e.returncode}")
            return False
        except Exception as e:
            self.logger.error(f"MEME运行出错|Error running MEME: {e}")
            return False

    def _build_command(self) -> list:
        """
        构建MEME命令|Build MEME command

        Returns:
            list: 命令列表|Command list
        """
        cmd = []

        # MEME可执行文件|MEME executable
        cmd.append(self.config.meme_path)

        # 输入文件|Input file
        cmd.append(self.config.input_file)

        # 输出目录|Output directory
        output_dir = self.config.get_meme_output_dir()
        cmd.extend(['-o', output_dir])

        # 序列类型|Sequence type
        if self.config.protein:
            cmd.append('-protein')
        elif self.config.dna:
            cmd.append('-dna')

        # Motif分布模式|Motif distribution mode
        cmd.extend(['-mod', self.config.mod])

        # Motif数量|Number of motifs
        cmd.extend(['-nmotifs', str(self.config.nmotifs)])

        # Motif宽度范围|Motif width range
        cmd.extend(['-minw', str(self.config.minw)])
        cmd.extend(['-maxw', str(self.config.maxw)])

        # 目标函数|Objective function
        cmd.extend(['-objfun', self.config.objfun])

        # Markov链阶数|Markov order
        cmd.extend(['-markov_order', str(self.config.markov_order)])

        return cmd

    def check_output_exists(self) -> bool:
        """
        检查MEME输出文件是否存在|Check if MEME output files exist

        Returns:
            bool: 是否存在|Whether exists
        """
        xml_file = Path(self.config.get_meme_output_dir()) / "meme.xml"
        txt_file = Path(self.config.get_meme_output_dir()) / "meme.txt"

        return xml_file.exists() or txt_file.exists()
