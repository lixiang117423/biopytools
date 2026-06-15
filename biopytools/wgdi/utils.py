"""
WGDI工具函数模块|WGDI Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional


class WGDILogger:
    """WGDI日志管理器|WGDI Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)
        handlers = [logging.StreamHandler(sys.stdout)]

        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = "."):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description

        Returns:
            bool: 是否成功|Whether succeeded
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class WGDIConfGenerator:
    """WGDI配置文件生成器|WGDI Configuration File Generator"""

    @staticmethod
    def generate_dotplot_conf(config) -> str:
        """生成DotPlot配置文件|Generate DotPlot configuration file"""
        conf_content = f"""[dotplot]
blast = {config.blast_file}
blast_reverse = {str(config.blast_reverse).lower()}
gff1 = {config.gff1_file}
gff2 = {config.gff2_file}
lens1 = {config.lens1_file}
lens2 = {config.lens2_file}
genome1_name = {config.genome1_name}
genome2_name = {config.genome2_name}
multiple = {config.multiple}
score = {config.score}
evalue = {config.evalue}
repeat_number = {config.repeat_number}
position = {config.position}
ancestor_left = {config.ancestor_left or 'none'}
ancestor_top = {config.ancestor_top or 'none'}
markersize = {config.markersize}
figsize = {config.figsize}
savefig = {config.savefig}
"""
        return conf_content

    @staticmethod
    def generate_collinearity_conf(config) -> str:
        """生成Collinearity配置文件|Generate Collinearity configuration file"""
        conf_content = f"""[collinearity]
gff1 = {config.gff1_file}
gff2 = {config.gff2_file}
lens1 = {config.lens1_file}
lens2 = {config.lens2_file}
blast = {config.blast_file}
blast_reverse = {str(config.blast_reverse).lower()}
comparison = {config.comparison}
multiple = {config.multiple}
process = {config.process}
evalue = {config.evalue}
score = {config.score}
grading = {config.grading}
mg = {config.mg}
pvalue = {config.pvalue}
repeat_number = {config.repeat_number}
position = {config.position}
savefile = {config.savefile}
"""
        return conf_content

    @staticmethod
    def generate_calks_conf(config) -> str:
        """生成CalKs配置文件|Generate CalKs configuration file"""
        conf_content = f"""[calks]
collinearity = {config.collinearity_file}
fasta1 = {config.fasta1_file}
fasta2 = {config.fasta2_file}
savefile = {config.savefile}
"""
        return conf_content
