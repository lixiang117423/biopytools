"""
RNA-Bloom组装工具函数模块|RNA-Bloom Assembly Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, List


class RNABloomLogger:
    """RNA-Bloom组装日志管理器|RNA-Bloom Assembly Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "rnabloom_assembly.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 如果日志文件存在，先删除|Remove existing log file
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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Optional[Path] = None):
        self.logger = logger
        self.working_dir = working_dir or Path.cwd()

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 命令描述|Command description

        Returns:
            bool: 成功返回True，失败返回False|True if successful, False otherwise
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
            if e.stdout:
                self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False

    def run_java_jar(self, jar_path: str, args: List[str], description: str = "",
                     java_options: str = "") -> bool:
        """执行Java JAR文件|Execute Java JAR file

        Args:
            jar_path: JAR文件路径|JAR file path
            args: 参数列表|Argument list
            description: 命令描述|Command description
            java_options: Java选项|Java options (e.g., "-Xmx4g")

        Returns:
            bool: 成功返回True，失败返回False|True if successful, False otherwise
        """
        # 构建Java命令|Build Java command
        cmd_parts = ["java"]
        if java_options:
            cmd_parts.append(java_options)
        cmd_parts.extend(["-jar", jar_path])
        cmd_parts.extend(args)

        cmd = " ".join(cmd_parts)
        return self.run(cmd, description)


class SequenceValidator:
    """序列文件验证器|Sequence File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def validate_fastq(self, file_path: str) -> bool:
        """验证FASTQ文件格式|Validate FASTQ file format

        Args:
            file_path: FASTQ文件路径|FASTQ file path

        Returns:
            bool: 有效返回True，无效返回False|True if valid, False otherwise
        """
        if not os.path.exists(file_path):
            self.logger.error(f"文件不存在|File does not exist: {file_path}")
            return False

        # 检查文件是否为空|Check if file is empty
        if os.path.getsize(file_path) == 0:
            self.logger.error(f"文件为空|File is empty: {file_path}")
            return False

        # 检查文件是否为压缩文件|Check if file is compressed
        if file_path.endswith('.gz'):
            try:
                result = subprocess.run(
                    ["zcat", file_path],
                    capture_output=True,
                    text=True,
                    check=True
                )
                # 检查前几行|Check first few lines
                lines = result.stdout.split('\n')[:4]
            except subprocess.CalledProcessError:
                self.logger.error(f"无法读取压缩文件|Cannot read compressed file: {file_path}")
                return False
        else:
            try:
                with open(file_path, 'r') as f:
                    lines = [f.readline().strip() for _ in range(4)]
            except Exception as e:
                self.logger.error(f"无法读取文件|Cannot read file: {file_path}: {e}")
                return False

        # 验证FASTQ格式|Validate FASTQ format
        if len(lines) < 4:
            self.logger.error(f"FASTQ文件行数不足|FASTQ file has insufficient lines: {file_path}")
            return False

        if not lines[0].startswith('@'):
            self.logger.error(f"FASTQ文件第一行应以@开头|FASTQ file first line should start with @: {file_path}")
            return False

        try:
            quality_length = len(lines[3])
            sequence_length = len(lines[1])
            if quality_length != sequence_length:
                self.logger.warning(
                    f"FASTQ质量分数长度与序列长度不匹配|"
                    f"Quality score length doesn't match sequence length: {file_path}"
                )
        except Exception:
            pass

        self.logger.info(f"FASTQ文件验证通过|FASTQ file validation passed: {file_path}")
        return True

    def validate_fasta(self, file_path: str) -> bool:
        """验证FASTA文件格式|Validate FASTA file format

        Args:
            file_path: FASTA文件路径|FASTA file path

        Returns:
            bool: 有效返回True，无效返回False|True if valid, False otherwise
        """
        if not os.path.exists(file_path):
            self.logger.error(f"文件不存在|File does not exist: {file_path}")
            return False

        if os.path.getsize(file_path) == 0:
            self.logger.error(f"文件为空|File is empty: {file_path}")
            return False

        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"FASTA文件第一行应以>开头|FASTA file first line should start with >: {file_path}")
                    return False
        except Exception as e:
            self.logger.error(f"无法读取文件|Cannot read file: {file_path}: {e}")
            return False

        self.logger.info(f"FASTA文件验证通过|FASTA file validation passed: {file_path}")
        return True

    def count_sequences(self, file_path: str) -> int:
        """统计序列数量|Count number of sequences

        Args:
            file_path: 序列文件路径|Sequence file path

        Returns:
            int: 序列数量|Number of sequences
        """
        try:
            if file_path.endswith('.gz'):
                result = subprocess.run(
                    ["zcat", file_path, "|", "grep", "-c", "^@"],
                    capture_output=True,
                    text=True,
                    shell=True
                )
            else:
                result = subprocess.run(
                    ["grep", "-c", "^>", file_path],
                    capture_output=True,
                    text=True,
                    check=False
                )

            if result.returncode == 0 and result.stdout.strip():
                return int(result.stdout.strip())
            else:
                return 0
        except Exception as e:
            self.logger.warning(f"无法统计序列数量|Cannot count sequences: {file_path}: {e}")
            return 0


class DependencyChecker:
    """依赖检查器|Dependency Checker"""

    def __init__(self, logger):
        self.logger = logger

    def check_java(self) -> bool:
        """检查Java是否可用|Check if Java is available"""
        try:
            result = subprocess.run(
                ["java", "-version"],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"Java已安装|Java is available: {result.stderr.split()[2]}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError, IndexError):
            self.logger.error("Java未安装或不在PATH中|Java is not installed or not in PATH")
            return False

    def check_minimap2(self) -> bool:
        """检查minimap2是否可用|Check if minimap2 is available"""
        try:
            result = subprocess.run(
                ["minimap2", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"minimap2已安装|minimap2 is available: {result.stdout.strip()}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("minimap2未安装或不在PATH中|minimap2 is not installed or not in PATH")
            return False

    def check_ntcard(self) -> bool:
        """检查ntCard是否可用|Check if ntCard is available"""
        try:
            result = subprocess.run(
                ["ntCard", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"ntCard已安装|ntCard is available")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.warning("ntCard未安装或不在PATH中|ntCard is not installed or not in PATH")
            return False

    def check_rnabloom(self, rnabloom_path: str = "rnabloom") -> bool:
        """检查RNA-Bloom是否可用|Check if RNA-Bloom is available

        Args:
            rnabloom_path: RNA-Bloom路径|RNA-Bloom path

        Returns:
            bool: 可用返回True|True if available
        """
        try:
            # 如果是JAR文件路径|If it's a JAR file path
            if rnabloom_path.endswith('.jar'):
                return os.path.exists(rnabloom_path)

            # 尝试运行rnabloom命令|Try to run rnabloom command
            result = subprocess.run(
                [rnabloom_path, "--help"],
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"RNA-Bloom可用|RNA-Bloom is available")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error(
                f"RNA-Bloom不可用|RNA-Bloom is not available: {rnabloom_path}. "
                f"请安装RNA-Bloom: conda install -c bioconda rnabloom"
            )
            return False
