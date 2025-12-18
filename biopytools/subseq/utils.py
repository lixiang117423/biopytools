"""
序列子集提取工具函数模块 | Sequence Subsequence Extraction Utility Functions Module
"""

import logging
import os
import sys
from pathlib import Path


class SubseqLogger:
    """序列子集提取日志管理器 | Sequence Subsequence Extraction Logger Manager"""

    def __init__(self, output_dir, log_name: str = "subseq_extraction.log"):
        # 确保output_dir是Path对象 | Ensure output_dir is Path object
        self.output_dir = Path(output_dir) if not isinstance(output_dir, Path) else output_dir
        self.log_file = self.output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger


def validate_file_exists(file_path: str, description: str = "文件"):
    """验证文件是否存在 | Validate file exists"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"📂❌ {description}不存在 | {description} does not exist: {file_path}")
    return True


def validate_id_list(id_list_file: str, logger):
    """验证ID列表文件格式 | Validate ID list file format"""
    try:
        with open(id_list_file, 'r') as f:
            id_list = [line.strip() for line in f if line.strip()]

        if not id_list:
            logger.warning("⚠️ ID列表文件为空 | ID list file is empty")
            return False

        logger.info(f"📋 从文件中读取到 {len(id_list)} 个ID | Read {len(id_list)} IDs from file")
        return id_list

    except Exception as e:
        logger.error(f"❌ 读取ID列表文件失败 | Failed to read ID list file: {e}")
        return False