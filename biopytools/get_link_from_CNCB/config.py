"""
CNCB链接提取器配置管理模块 | CNCB Link Extractor Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class CNCBConfig:
    """CNCB链接提取器配置类 | CNCB Link Extractor Configuration Class"""

    # 必需输入 | Required inputs
    input_file: str

    # 输出文件配置 | Output file configuration
    output_file: Optional[str] = None
    failed_file: Optional[str] = None
    download_script: str = "download.sh"

    # FTP服务器配置 | FTP server configuration
    ftp_host: str = "download2.cncb.ac.cn"
    ftp_timeout: int = 60
    base_dirs: list = None

    # 搜索配置 | Search configuration
    max_threads: int = 4
    retry_attempts: int = 3
    retry_delay: float = 1.0

    # 日志配置 | Logging configuration
    verbose: bool = False
    log_file: Optional[str] = None

    # 输出格式配置 | Output format configuration
    generate_download_script: bool = True
    script_executable: bool = True
    include_timestamp: bool = True

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 设置默认base_dirs
        if self.base_dirs is None:
            self.base_dirs = ['INSDC', 'INSDC2', 'INSDC3', 'INSDC4', 'INSDC5']

        # 生成默认输出文件名
        if not self.output_file:
            base_name = os.path.splitext(os.path.basename(self.input_file))[0]
            self.output_file = f"{base_name}_links.txt"

        if not self.failed_file:
            base_name = os.path.splitext(os.path.basename(self.input_file))[0]
            self.failed_file = f"{base_name}_failed.txt"

        # 验证输入文件路径
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))

        # 创建输出目录
        output_dir = os.path.dirname(os.path.abspath(self.output_file))
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        failed_dir = os.path.dirname(os.path.abspath(self.failed_file))
        Path(failed_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查输入文件 | Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在 | Input file does not exist: {self.input_file}")

        if not os.path.isfile(self.input_file):
            errors.append(f"输入路径不是文件 | Input path is not a file: {self.input_file}")

        # 检查文件可读性 | Check file readability
        try:
            with open(self.input_file, 'r', encoding='utf-8') as f:
                f.read(1)  # 尝试读取第一个字符
        except (UnicodeDecodeError, PermissionError, IOError) as e:
            errors.append(f"输入文件无法读取 | Input file cannot be read: {e}")

        # 检查输出目录可写性 | Check output directory writability
        try:
            with open(self.output_file, 'a', encoding='utf-8'):
                pass
        except (PermissionError, IOError) as e:
            errors.append(f"输出文件无法写入 | Output file cannot be written: {e}")

        try:
            with open(self.failed_file, 'a', encoding='utf-8'):
                pass
        except (PermissionError, IOError) as e:
            errors.append(f"失败记录文件无法写入 | Failed file cannot be written: {e}")

        # 检查配置值 | Check configuration values
        if self.ftp_timeout <= 0:
            errors.append("FTP超时时间必须大于0 | FTP timeout must be greater than 0")

        if self.max_threads < 1:
            errors.append("最大线程数必须大于等于1 | Max threads must be >= 1")

        if self.retry_attempts < 0:
            errors.append("重试次数不能为负数 | Retry attempts cannot be negative")

        if self.retry_delay < 0:
            errors.append("重试延迟不能为负数 | Retry delay cannot be negative")

        # 检查脚本文件名 | Check script filename
        if self.download_script:
            if not self.download_script.strip():
                errors.append("下载脚本文件名不能为空 | Download script filename cannot be empty")

            # 检查是否包含路径
            if os.path.dirname(self.download_script):
                script_dir = os.path.dirname(os.path.abspath(self.download_script))
                try:
                    Path(script_dir).mkdir(parents=True, exist_ok=True)
                except (PermissionError, OSError) as e:
                    errors.append(f"下载脚本目录无法创建 | Download script directory cannot be created: {e}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_archive_map(self):
        """获取数据库映射 | Get archive mapping"""
        return {
            "SRR": "SRA",
            "ERR": "ERA",
            "DRR": "DRA"
        }

    def get_filename_templates(self):
        """获取文件名模板 | Get filename templates"""
        return [
            "{run_id}.sra",
            "{run_id}_1.fastq.gz",
            "{run_id}_2.fastq.gz",
            "{run_id}.fastq.gz",
            "{run_id}"
        ]

    def get_output_info(self):
        """获取输出文件信息 | Get output file information"""
        return {
            "links_file": self.output_file,
            "failed_file": self.failed_file,
            "download_script": self.download_script if self.generate_download_script else None,
            "log_file": self.log_file
        }