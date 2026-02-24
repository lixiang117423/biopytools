"""
BUSCO分析工具函数模块|BUSCO Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import glob
import os
import re
from pathlib import Path
from typing import List, Tuple

class BUSCOLogger:
    """BUSCO分析日志管理器|BUSCO Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "busco_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器|Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> Tuple[bool, str, str]:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

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

            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False, e.stdout, e.stderr

class FileManager:
    """文件管理器|File Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def get_input_files(self) -> List[Tuple[str, str]]:
        """获取输入文件列表|Get input file list

        Returns:
            List[Tuple[str, str]]: [(file_path, sample_name), ...]
        """
        input_path = Path(self.config.input_path)
        files = []

        if input_path.is_file():
            # 单个文件|Single file
            sample_name = self.extract_sample_name(input_path.name)
            files.append((str(input_path), sample_name))
            self.logger.info(f"发现单个输入文件|Found single input file: {input_path.name}")

        elif input_path.is_dir():
            # 目录批处理|Directory batch processing
            pattern = self.config.sample_suffix.replace('*', '*')
            search_pattern = input_path / pattern

            matching_files = glob.glob(str(search_pattern))
            if not matching_files:
                raise ValueError(f"目录中未找到匹配模式 '{pattern}' 的文件|No files matching pattern '{pattern}' found in directory")

            for file_path in sorted(matching_files):
                file_name = Path(file_path).name
                sample_name = self.extract_sample_name(file_name)
                files.append((file_path, sample_name))

            self.logger.info(f"发现批处理文件|Found batch files: {len(files)} 个文件|files")

        else:
            raise ValueError(f"输入路径无效|Invalid input path: {input_path}")

        return files
    
    def extract_sample_name(self, filename: str) -> str:
        """从文件名提取样本名|Extract sample name from filename"""
        suffix_pattern = self.config.sample_suffix
        
        # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
        # 转义特殊字符，但保留*|Escape special chars but keep *
        escaped_pattern = re.escape(suffix_pattern)
        # 将转义的\*替换为捕获组|Replace escaped \* with capture group
        regex_pattern = escaped_pattern.replace(r'\*', r'(.*)')
        
        match = re.match(regex_pattern, filename)
        if match:
            return match.group(1)
        else:
            # 如果模式不匹配，使用文件名去除扩展名|If pattern doesn't match, use filename without extension
            return Path(filename).stem

def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        result = subprocess.run([config.busco_path, "--version"],
                              capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            version_info = result.stdout.strip()
            logger.info(f"BUSCO 可用|BUSCO available: {version_info}")
        else:
            raise RuntimeError("BUSCO版本检查失败|BUSCO version check failed")
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"BUSCO不可用|BUSCO not available: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True
