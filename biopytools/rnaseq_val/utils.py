"""
转录组验证注释工具函数模块|Transcriptome Validation Utility Functions Module
"""

import os
import glob
import logging
import subprocess
import sys
import time
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple


# ============================================================
# 日志管理器|Logger Manager
# ============================================================

class RnaseqValLogger:
    """转录组验证日志管理器|Transcriptome Validation Logger Manager"""

    def __init__(self, log_dir: Path, verbose: bool = False, quiet: bool = False):
        """初始化日志管理器|Initialize logger

        Args:
            log_dir: 日志目录|Log directory
            verbose: 详细模式|Verbose mode
            quiet: 静默模式|Quiet mode
        """
        self.log_dir = log_dir
        self.log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = log_dir / f"rnaseq_val_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"rnaseq_val_{timestamp}")
        self.logger.propagate = False

        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        # 清除已有 handler|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件 handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger实例|Get logger instance

        Returns:
            logging.Logger: logger 实例
        """
        return self.logger

    def step(self, message: str):
        """记录步骤分隔|Log step separator

        Args:
            message: 步骤描述|Step description
        """
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


# ============================================================
# 命令执行器|Command Runner
# ============================================================

class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        """初始化命令执行器|Initialize command runner

        Args:
            logger: logger 实例
        """
        self.logger = logger

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令字符串|Command string to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            bool: 执行成功返回 True|True if command succeeded
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        if timeout:
            self.logger.info(
                f"超时设置|Timeout: {timeout}秒|seconds ({timeout / 3600:.1f}小时|hours)"
            )

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            self.logger.info(f"{description} 完成|completed")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error(
                f"{description} 超时|timed out after {timeout}秒|seconds "
                f"({timeout / 3600:.1f}小时|hours)"
            )
            return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


# ============================================================
# 文件验证器|File Validator
# ============================================================

class FileValidator:
    """文件验证器|File Validator"""

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在|Check if file exists

        Args:
            file_path: 文件路径|File path
            description: 文件描述|File description

        Returns:
            bool: 文件存在返回 True|True if file exists
        """
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False


# ============================================================
# 样本解析器|Sample Parser
# ============================================================

class SampleParser:
    """样本解析器，支持二代配对和三代单端自动检测|Sample parser for paired SR and single-end LR"""

    # 二代配对 fastq 的 read 标识符|Read indicators for paired-end SR fastq
    READ_INDICATORS = [
        ("R1", "R2"),
        ("_1", "_2"),
        (".1", ".2"),
        ("_f1", "_r2"),
        ("_f1", "_f2"),
        ("_F1", "_R2"),
        ("_F1", "_F2"),
        ("_read1", "_read2"),
        ("_pair1", "_pair2"),
    ]

    # 三代单端 fastq 的扩展名模式|Extension patterns for single-end LR fastq
    LR_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    # ----- 二代样本解析 ----- SR sample parsing -----

    def parse_sr_samples(self, sr_dir: str, pattern: Optional[str] = None) -> List[Dict]:
        """解析二代配对样本|Parse paired-end short-read samples

        Args:
            sr_dir: 二代 reads 目录|Short-read directory
            pattern: 自定义命名模式|Custom naming pattern (e.g. *_1.clean.fq.gz)

        Returns:
            List[Dict]: 样本列表，每项包含 name/fastq1/fastq2
        """
        if pattern:
            self.logger.info(f"使用指定模式|Using specified pattern: {pattern}")
            return self._parse_sr_with_pattern(sr_dir, pattern)

        self.logger.info("自动检测二代配对样本|Auto-detecting paired-end SR samples")
        return self._parse_sr_default(sr_dir)

    def _parse_sr_with_pattern(self, sr_dir: str, pattern: str) -> List[Dict]:
        """使用指定模式解析二代样本|Parse SR samples with specified pattern"""
        samples = []

        if "*" not in pattern:
            self.logger.error("文件模式必须包含 * 作为样本名占位符|Pattern must contain * placeholder")
            return samples

        parts = pattern.split("*")
        if len(parts) != 2:
            self.logger.error("文件模式只能包含一个 * 占位符|Pattern can only contain one * placeholder")
            return samples

        prefix, suffix = parts
        r1_indicator = r2_indicator = None
        for r1, r2 in self.READ_INDICATORS:
            if r1 in suffix:
                r1_indicator, r2_indicator = r1, r2
                break

        if not r1_indicator:
            self.logger.error(
                f"无法从模式中识别 read 标识符|Cannot identify read indicator in pattern: {pattern}"
            )
            return samples

        search_pattern = os.path.join(sr_dir, pattern)
        fastq_files = sorted(glob.glob(search_pattern))

        self.logger.info(f"搜索模式|Search pattern: {search_pattern}")
        self.logger.info(f"找到 {len(fastq_files)} 个 read1 文件|Found {len(fastq_files)} read1 files")

        for fq1 in fastq_files:
            basename = os.path.basename(fq1)
            sample_name = basename.replace(prefix, "", 1).replace(suffix, "", 1)
            r2_suffix = suffix.replace(r1_indicator, r2_indicator)
            fq2 = os.path.join(sr_dir, prefix + sample_name + r2_suffix)

            if os.path.exists(fq2):
                samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
            else:
                self.logger.warning(f"找不到配对的 read2 文件|Missing read2: {fq2}")

        return samples

    def _parse_sr_default(self, sr_dir: str) -> List[Dict]:
        """使用默认模式自动检测二代配对样本|Auto-detect SR paired samples with default patterns"""
        samples = []

        patterns = [
            ("*_1.fq.gz", "*_2.fq.gz"),
            ("*_R1.fq.gz", "*_R2.fq.gz"),
            ("*.R1.fastq.gz", "*.R2.fastq.gz"),
            ("*_1.fastq.gz", "*_2.fastq.gz"),
            ("*.1.fq.gz", "*.2.fq.gz"),
            ("*_f1.fq.gz", "*_r2.fq.gz"),
            ("*_f1.fq.gz", "*_f2.fq.gz"),
            ("*_F1.fq.gz", "*_R2.fq.gz"),
            ("*_read1.fq.gz", "*_read2.fq.gz"),
            ("*_pair1.fq.gz", "*_pair2.fq.gz"),
        ]

        for p1, p2 in patterns:
            if "*" not in p1:
                continue
            prefix, suffix = p1.split("*", 1)
            r1_ind = r2_ind = None
            for r1, r2 in self.READ_INDICATORS:
                if r1 in suffix:
                    r1_ind, r2_ind = r1, r2
                    break
            if not r1_ind:
                continue

            search_path = os.path.join(sr_dir, p1)
            fq1_files = sorted(glob.glob(search_path))

            if not fq1_files:
                continue

            self.logger.info(f"尝试模式|Trying pattern: {p1} / {p2}")

            for fq1 in fq1_files:
                basename = os.path.basename(fq1)
                sample_name = basename
                if prefix:
                    sample_name = sample_name.replace(prefix, "", 1)
                if suffix:
                    sample_name = sample_name.replace(suffix, "", 1)

                r2_suffix = suffix.replace(r1_ind, r2_ind)
                fq2 = os.path.join(sr_dir, prefix + sample_name + r2_suffix)

                if os.path.exists(fq2):
                    samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})

            if samples:
                self.logger.info(f"使用模式|Using pattern: {p1} / {p2}")
                break

        return samples

    # ----- 三代样本解析 ----- LR sample parsing -----

    def parse_lr_samples(self, lr_dir: str) -> List[Dict]:
        """解析三代单端样本|Parse single-end long-read samples

        从目录中自动检测 fastq/fq.gz 文件（排除可能的配对 read2 文件）

        Args:
            lr_dir: 三代 reads 目录|Long-read directory

        Returns:
            List[Dict]: 样本列表，每项包含 name/reads
        """
        samples = []
        all_files = sorted(os.listdir(lr_dir))

        # 收集所有 fastq 文件|Collect all fastq files
        fq_files = []
        for f in all_files:
            if f.endswith(self.LR_EXTENSIONS) and not f.startswith("."):
                fq_files.append(f)

        # 过滤掉可能的 read2 文件|Filter out potential read2 files
        r2_indicators = {"_R2", "_2.", ".r2", "_r2", "_2.", "_read2", "_pair2"}
        sr_samples = []

        for f in fq_files:
            is_r2 = False
            for r2 in r2_indicators:
                if r2 in f.upper() or r2 in f.lower():
                    is_r2 = True
                    break
            if not is_r2:
                sr_samples.append(f)

        # 去重：如果有同名不同扩展名只保留一个|Dedup: keep one per unique base name
        seen_names = set()
        for f in sr_samples:
            # 取扩展名前的部分作为样本名|Base name before extension
            base = f
            for ext in self.LR_EXTENSIONS:
                if base.endswith(ext):
                    base = base[: -len(ext)]
                    break
            if base not in seen_names:
                seen_names.add(base)
                full_path = os.path.join(lr_dir, f)
                samples.append({"name": base, "reads": full_path})

        self.logger.info(f"在目录中检测到 {len(samples)} 个三代样本|Detected {len(samples)} LR samples in {lr_dir}")
        for s in samples:
            self.logger.info(f"  - {s['name']}: {os.path.basename(s['reads'])}")

        return samples


# ============================================================
# 辅助函数|Utility Functions
# ============================================================

def build_conda_command(cmd: str, env_name: str = "rnaseq_val") -> str:
    """构建 conda run 包装命令|Build conda run wrapped command

    Args:
        cmd: 原始命令字符串|Original command string
        env_name: conda 环境名称|Conda environment name

    Returns:
        str: 包装后的命令字符串|Wrapped command string
    """
    return f"conda run -n {env_name} --no-capture-output {cmd}"


def format_number(num: int) -> str:
    """格式化大数字|Format large numbers

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
