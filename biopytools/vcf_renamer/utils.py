"""
VCF重命名工具类 | VCF Renamer Utility Classes
包含日志和辅助函数 | Contains logging and helper functions
"""

import logging
import sys
import os
import subprocess
from typing import List
from pathlib import Path


class VCFRenamerLogger:
    """VCF重命名日志管理器 | VCF Renamer Logger Manager"""

    def __init__(self, output_dir: str = "."):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            output_dir: 输出目录 | Output directory
        """
        self.output_dir = output_dir
        log_file = os.path.join(output_dir, 'vcf_renamer.log') if output_dir else None
        self.setup_logging(log_file)

    def setup_logging(self, log_file: str = None):
        """设置日志 | Setup logging"""
        self.logger = logging.getLogger("VCFRenamer")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式 | Log format
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler(如果指定) | File handler (if specified)
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class VCFRenamerChecker:
    """VCF重命名检查器 | VCF Renamer Checker"""

    def __init__(self, logger):
        """
        初始化检查器 | Initialize checker

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def check_bcftools(self) -> bool:
        """
        检查bcftools是否可用 | Check if bcftools is available

        Returns:
            bcftools是否可用 | Whether bcftools is available
        """
        if subprocess.run(['which', 'bcftools'], capture_output=True).returncode != 0:
            self.logger.error("❌ 未找到bcftools命令 | bcftools command not found")
            self.logger.error("请先安装bcftools | Please install bcftools first:")
            self.logger.error("  conda install -c bioconda bcftools")
            return False
        return True

    def extract_samples(self, vcf_file: str) -> List[str]:
        """
        提取VCF文件中的样品名称 | Extract sample names from VCF file

        Args:
            vcf_file: VCF文件路径 | VCF file path

        Returns:
            样品名称列表 | List of sample names
        """
        self.logger.info("📋 正在提取样品名称... | Extracting sample names...")

        try:
            result = subprocess.run(
                ['bcftools', 'query', '-l', vcf_file],
                capture_output=True,
                text=True,
                check=True
            )
            samples = [s.strip() for s in result.stdout.strip().split('\n') if s.strip()]
            self.logger.info(f"✅ 检测到 {len(samples)} 个样品 | Found {len(samples)} samples")
            return samples
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 提取样品名称失败 | Failed to extract samples: {e.stderr}")
            return []

    def validate_vcf(self, vcf_file: str) -> bool:
        """
        验证VCF文件 | Validate VCF file

        Args:
            vcf_file: VCF文件路径 | VCF file path

        Returns:
            是否有效 | Whether valid
        """
        try:
            result = subprocess.run(
                ['bcftools', 'index', '-i', vcf_file],
                capture_output=True,
                text=True
            )
            return True
        except Exception as e:
            self.logger.error(f"❌ VCF文件验证失败 | VCF validation failed: {e}")
            return False


class VCFRenamerProcessor:
    """VCF重命名处理器 | VCF Renamer Processor"""

    def __init__(self, logger):
        """
        初始化处理器 | Initialize processor

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def generate_mapping(self, samples: List[str], prefix: str, mapping_file: str) -> bool:
        """
        生成新旧样品名称映射表 | Generate new-old sample name mapping

        Args:
            samples: 原始样品名称列表 | Original sample name list
            prefix: 新样品名前缀 | New sample name prefix
            mapping_file: 映射文件路径 | Mapping file path

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info(f"📝 正在生成映射表 (前缀: {prefix})... | Generating mapping (prefix: {prefix})...")

        try:
            with open(mapping_file, 'w') as f:
                for i, sample in enumerate(samples, 1):
                    new_name = f"{prefix}{i}"
                    f.write(f"{sample}\t{new_name}\n")

            self.logger.info(f"✅ 映射文件已生成 | Mapping file generated: {mapping_file}")
            return True
        except Exception as e:
            self.logger.error(f"❌ 生成映射文件失败 | Failed to generate mapping: {e}")
            return False

    def rename_samples(self, input_vcf: str, output_vcf: str, mapping_file: str) -> bool:
        """
        使用bcftools重命名样品 | Rename samples using bcftools

        Args:
            input_vcf: 输入VCF文件 | Input VCF file
            output_vcf: 输出VCF文件 | Output VCF file
            mapping_file: 映射文件 | Mapping file

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("🔄 正在重命名样品... | Renaming samples...")

        try:
            subprocess.run(
                ['bcftools', 'reheader', '-s', mapping_file, '-o', output_vcf, input_vcf],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"✅ 样品重命名完成 | Sample renaming completed: {output_vcf}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 重命名失败 | Renaming failed: {e.stderr}")
            return False

    def index_vcf(self, vcf_file: str) -> bool:
        """
        为VCF文件创建索引 | Create index for VCF file

        Args:
            vcf_file: VCF文件路径 | VCF file path

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("📇 正在创建VCF索引... | Creating VCF index...")

        try:
            subprocess.run(
                ['bcftools', 'index', '-t', vcf_file],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"✅ 索引创建完成 | Index created: {vcf_file}.tbi")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 创建索引失败 | Failed to create index: {e.stderr}")
            return False
