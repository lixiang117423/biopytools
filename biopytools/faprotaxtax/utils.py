"""
FAPROTAX工具函数模块|FAPROTAX Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple


class FaprotaxtaxLogger:
    """FAPROTAX日志管理器|FAPROTAX Logger Manager"""

    def __init__(self, log_file: str):
        self.log_file = Path(log_file)
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 创建日志目录|Create log directory
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 删除旧日志|Delete old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 配置根日志|Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        root_logger.handlers.clear()
        root_logger.propagate = False

        root_logger.addHandler(file_handler)
        root_logger.addHandler(stdout_handler)
        root_logger.addHandler(stderr_handler)

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str):
        self.logger = logger
        self.working_dir = working_dir

    def run_command(self, cmd: List[str], description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（shell=False安全模式）|Command list (shell=False safe mode)
            description: 步骤描述|Step description

        Returns:
            (success, stdout, stderr): 执行结果|Execution result
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 记录完整命令到INFO级别|Log complete command at INFO level
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                self.logger.error(f"错误代码|Error code: {result.returncode}")
                if result.stderr:
                    # 记录最后2000字符|Log last 2000 chars
                    stderr_tail = result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr
                    self.logger.error(f"错误信息|Error message: {stderr_tail}")

                    # 保存完整stderr到文件|Save full stderr to file
                    stderr_file = os.path.join(
                        self.working_dir, "99_logs",
                        f"stderr_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
                    )
                    os.makedirs(os.path.dirname(stderr_file), exist_ok=True)
                    with open(stderr_file, 'w') as f:
                        f.write(result.stderr)
                    self.logger.error(f"完整错误日志已保存|Full error log saved: {stderr_file}")
                return False, result.stdout, result.stderr

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True, result.stdout, result.stderr

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {description}")
            self.logger.error(f"异常信息|Exception: {str(e)}")
            return False, "", str(e)


def detect_input_format(filepath: str) -> str:
    """
    检测输入文件格式|Detect input file format

    Args:
        filepath: 输入文件路径|Input file path

    Returns:
        'biom' 或 'tsv'|'biom' or 'tsv'
    """
    lower = filepath.lower()
    if any(lower.endswith(ext) for ext in ['.biom', '.biom.gz', '.hbiom', '.hbiom.gz', '.jbiom', '.jbiom.gz']):
        return 'biom'
    return 'tsv'


def extract_faprotaxtax_version(script_path: str) -> str:
    """
    从collapse_table.py头部提取版本号|Extract version from collapse_table.py header

    Args:
        script_path: collapse_table.py路径|Path to collapse_table.py

    Returns:
        版本字符串或'unknown'|Version string or 'unknown'
    """
    import re
    try:
        with open(script_path, 'r') as f:
            for _ in range(60):
                line = f.readline()
                if not line:
                    break
                match = re.search(r'Script version:\s*([\d.]+)', line)
                if match:
                    return match.group(1)
    except Exception:
        pass
    return "unknown"


def generate_software_versions_yml(
    output_dir: str,
    config,
    start_time: datetime
) -> None:
    """
    生成software_versions.yml文件|Generate software_versions.yml file

    Args:
        output_dir: 输出目录|Output directory
        config: FaprotaxtaxConfig实例|FaprotaxtaxConfig instance
        start_time: 开始时间|Start time
    """
    try:
        import yaml
    except ImportError:
        logging.getLogger(__name__).warning(
            "PyYAML未安装，跳过生成software_versions.yml|PyYAML not installed, skipping software_versions.yml"
        )
        return

    end_time = datetime.now()
    runtime_seconds = int((end_time - start_time).total_seconds())

    # 获取collapse_table.py版本|Get collapse_table.py version
    faprotaxtax_version = extract_faprotaxtax_version(config.collapse_table_path)

    # 获取数据库信息|Get database info
    db_info = {}
    if os.path.exists(config.groups_file):
        db_info['path'] = config.groups_file
        db_info['size_mb'] = round(os.path.getsize(config.groups_file) / (1024 * 1024), 2)
        try:
            with open(config.groups_file, 'r') as f:
                # 读取第2行的版本信息|Read version from line 2
                first_line = f.readline()
                second_line = f.readline()
                if 'Version:' in second_line:
                    db_info['version'] = second_line.split('Version:')[1].strip()
        except Exception:
            pass

    info = {
        'pipeline': {
            'name': 'biopytools faprotaxtax',
            'version': '1.0.0'
        },
        'tools': {
            'collapse_table': {
                'version': faprotaxtax_version,
                'path': config.collapse_table_path
            },
            'python_interpreter': {
                'version': sys.version.split()[0],
                'path': config.python_interpreter
            }
        },
        'database': db_info,
        'parameters': {
            'input_table': config.input_table,
            'output_dir': config.output_dir,
            'groups_file': config.groups_file,
            'collapse_by_metadata': config.collapse_by_metadata,
            'group_leftovers_as': config.group_leftovers_as,
            'normalize': config.normalize,
            'average': config.average,
            'row_names_are_in_column': config.row_names_are_in_column,
            'output_format': config.output_format,
            'threads': config.threads
        },
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': runtime_seconds
        }
    }

    output_file = Path(output_dir) / '00_pipeline_info' / 'software_versions.yml'
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
