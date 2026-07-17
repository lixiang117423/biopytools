"""
Fastq到VCF (GTX) 工具模块|Fastq to VCF (GTX) Utilities Module
"""

import os
import subprocess
import time
import logging
import threading
from pathlib import Path
from typing import List, Optional, Dict, Any


class Fastq2VcfGTXLogger:
    """Fastq到VCF (GTX) 日志管理器|Fastq to VCF (GTX) Logger Manager"""

    def __init__(self, log_dir: str, verbose: bool = False, quiet: bool = False):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        self.quiet = quiet

        # 创建日志文件|Create log files
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.log_dir / f"pipeline_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"fastq2vcf_gtx_{timestamp}")
        self.logger.propagate = False  # 避免日志重复|Avoid duplicate logs

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        import sys
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

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤|Log step"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str, dry_run: bool = False):
        self.logger = logger
        self.output_dir = output_dir
        self.dry_run = dry_run

    def run(self, command, description: str = "") -> bool:
        """
        运行命令|Run command

        Args:
            command: 要执行的命令（字符串或列表）|Command to execute (string or list)
            description: 命令描述|Command description

        Returns:
            bool: 是否成功|Success status
        """
        if description:
            self.logger.info(f"{description}")

        # 格式化命令显示|Format command display
        if isinstance(command, list):
            command_str = ' '.join(command)
        else:
            command_str = command

        self.logger.info(f"执行命令|Running command: {command_str}")

        if self.dry_run:
            self.logger.info(f"[DRY RUN] 跳过实际执行|Skip actual execution")
            return True

        try:
            # 执行命令|Execute command
            if isinstance(command, list):
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=86400,  # 24小时超时|24 hours timeout
                    cwd=self.output_dir
                )
            else:
                result = subprocess.run(
                    command,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=86400,  # 24小时超时|24 hours timeout
                    cwd=self.output_dir
                )

            if result.returncode == 0:
                if result.stdout:
                    self.logger.info(f"标准输出|Stdout: {result.stdout}")
                else:
                    self.logger.info("标准输出|Stdout: (空输出)")

                if result.stderr:
                    self.logger.info(f"标准错误|Stderr: {result.stderr}")

                # 检查并读取fastp日志文件|Check and read fastp log files
                self._read_fastp_logs(command_str, command, log_level=logging.INFO)

                self.logger.info("命令执行成功|Command executed successfully")
                return True
            else:
                self.logger.error(f"命令执行失败|Command execution failed")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                if result.stdout:
                    self.logger.error(f"标准输出|Stdout: {result.stdout}")
                if result.stderr:
                    self.logger.error(f"错误输出|Stderr: {result.stderr}")

                # 失败时也尝试读取fastp日志文件|Also try to read fastp log file on failure
                self._read_fastp_logs(command_str, command, log_level=logging.ERROR)

                return False

        except subprocess.TimeoutExpired:
            self.logger.error("命令执行超时|Command execution timed out")
            return False
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {str(e)}")
            return False

    def _read_fastp_logs(self, command_str: str, command, log_level=logging.INFO):
        """查找并读取fastp日志文件|Find and read fastp log files"""
        # command_str 已统一为字符串,直接判断是否含 fastp|command_str is normalized to str, check fastp directly
        if "fastp" not in command_str:
            return

        if log_level >= logging.ERROR:
            self.logger.info("尝试读取FastP日志文件以获取更多错误信息|Trying to read FastP log file for more error information")

        fastp_log_files = []

        # 检查输出目录中的fastp报告目录|Check fastp reports directory in output directory
        reports_dir = os.path.join(self.output_dir, "01_fastp", "fastp_reports") if self.output_dir else None
        if reports_dir and os.path.exists(reports_dir):
            for log_file in os.listdir(reports_dir):
                if log_file.endswith(".log"):
                    fastp_log_files.append(os.path.join(reports_dir, log_file))

        # 检查输出目录中的日志文件|Check log files in output directory
        if self.output_dir and os.path.exists(self.output_dir):
            for log_file in os.listdir(self.output_dir):
                if log_file.endswith("_processing.log") or log_file.endswith("fastp.log"):
                    fastp_log_files.append(os.path.join(self.output_dir, log_file))

        # 显示找到的日志文件内容|Show contents of found log files
        for log_file in fastp_log_files:
            if os.path.exists(log_file):
                label = "FastP日志文件内容|FastP Log File Content"
                if log_level >= logging.ERROR:
                    label = "FastP错误日志文件内容|FastP Error Log File Content"
                self.logger.log(log_level, f"{label}: {log_file}")
                try:
                    with open(log_file, 'r', encoding='utf-8') as f:
                        log_content = f.read()
                        if log_content:
                            for line in log_content.split('\n'):
                                if line.strip():
                                    self.logger.log(log_level, f"  {line}")
                        else:
                            self.logger.log(log_level, "  (空日志文件)|(Empty log file)")
                except Exception as e:
                    self.logger.error(f"读取FastP日志文件失败|Failed to read FastP log file {log_file}: {str(e)}")
            else:
                self.logger.warning(f"FastP日志文件不存在|FastP log file does not exist: {log_file}")

    def run_with_progress(self, command: str, description: str = "",
                         timeout: int = 86400) -> bool:
        """
        运行命令并流式显示进度|Run command with streamed progress

        使用独立 reader 线程持续排空 stdout,避免:
        - 管道缓冲区写满导致子进程阻塞死锁
        - readline() 阻塞在无换行的进度输出上
        - 循环结束后再次 communicate() 的反模式
        Uses a dedicated reader thread to drain stdout continuously, avoiding:
        - pipe-buffer-full deadlock blocking the subprocess
        - readline() blocking on newline-free progress output
        - the communicate()-after-read anti-pattern

        Args:
            command: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间（秒）|Timeout in seconds

        Returns:
            bool: 是否成功|Success status
        """
        if description:
            self.logger.info(f"{description}")

        self.logger.info(f"执行命令|Running command: {command}")

        if self.dry_run:
            self.logger.info(f"[DRY RUN] 跳过实际执行|Skip actual execution")
            return True

        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,  # 行缓冲,便于实时读取|Line-buffered for real-time reading
                cwd=self.output_dir
            )
        except Exception as e:
            self.logger.error(f"启动命令失败|Failed to start command: {str(e)}")
            return False

        # reader 线程:持续排空 stdout 到日志,防止管道写满死锁
        # Reader thread: continuously drain stdout to log, preventing buffer-full deadlock
        def _drain_output():
            try:
                for line in iter(process.stdout.readline, ''):
                    line = line.rstrip('\n')
                    if line:
                        # 实时输出到 INFO 级别,以便在 .out 文件中显示进度
                        # Stream to INFO so progress appears in the .out file
                        self.logger.info(line)
            except Exception:
                pass
            finally:
                try:
                    process.stdout.close()
                except Exception:
                    pass

        reader = threading.Thread(target=_drain_output, daemon=True)
        reader.start()

        start_time = time.time()
        timed_out = False

        try:
            while True:
                try:
                    # 每秒唤醒一次检查超时,不阻塞读取|Wake every second to check timeout
                    process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    if time.time() - start_time > timeout:
                        timed_out = True
                        self.logger.error("命令执行超时|Command execution timed out")
                        process.terminate()
                        try:
                            process.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait()
                        break

            # 等待 reader 线程排空剩余输出|Wait for reader to flush remaining output
            reader.join(timeout=5)
        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution exception: {str(e)}")
            return False

        if timed_out:
            return False

        if process.returncode == 0:
            self.logger.info("命令执行成功|Command executed successfully")
            return True
        else:
            self.logger.error(f"命令执行失败|Command execution failed")
            self.logger.error(f"返回码|Return code: {process.returncode}")
            return False


class CheckpointManager:
    """检查点管理器|Checkpoint Manager"""

    def __init__(self, checkpoint_dir: str, logger: logging.Logger):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger

    def exists(self, step: str) -> bool:
        """检查检查点是否存在|Check if checkpoint exists"""
        return (self.checkpoint_dir / f"{step}.done").exists()

    def create(self, step: str):
        """创建检查点|Create checkpoint"""
        checkpoint_file = self.checkpoint_dir / f"{step}.done"
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        checkpoint_file.write_text(timestamp)
        self.logger.info(f"创建检查点|Create checkpoint: {step}")

    def remove(self, step: str):
        """移除检查点|Remove checkpoint"""
        checkpoint_file = self.checkpoint_dir / f"{step}.done"
        if checkpoint_file.exists():
            checkpoint_file.unlink()
            self.logger.info(f"移除检查点|Remove checkpoint: {step}")

    def list_completed(self) -> List[str]:
        """列出已完成的检查点|List completed checkpoints"""
        completed = []
        for checkpoint_file in self.checkpoint_dir.glob("*.done"):
            step = checkpoint_file.stem
            timestamp = checkpoint_file.read_text()
            completed.append(f"{step} ({timestamp})")
        return completed


class FileManager:
    """文件管理器|File Manager"""

    @staticmethod
    def count_files(directory: str, pattern: str = "*") -> int:
        """统计目录中文件数量|Count files in directory"""
        if not os.path.exists(directory):
            return 0
        return len(list(Path(directory).glob(pattern)))

    @staticmethod
    def find_files(directory: str, pattern: str) -> List[str]:
        """查找文件|Find files"""
        if not os.path.exists(directory):
            return []
        return [str(f) for f in Path(directory).glob(pattern)]

    @staticmethod
    def ensure_directory(directory: str):
        """确保目录存在|Ensure directory exists"""
        Path(directory).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_file_size(file_path: str) -> str:
        """获取文件大小|Get file size"""
        if not os.path.exists(file_path):
            return "0 B"
        size = os.path.getsize(file_path)
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size < 1024:
                return f"{size:.1f} {unit}"
            size /= 1024
        return f"{size:.1f} PB"


class SystemChecker:
    """系统检查器|System Checker"""

    @staticmethod
    def check_command_exists(command: str, logger: logging.Logger) -> bool:
        """检查命令是否存在|Check if command exists"""
        try:
            result = subprocess.run(
                ["which", command],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                logger.info(f"命令存在|Command exists: {command}")
                return True
            else:
                logger.error(f"命令不存在|Command does not exist: {command}")
                return False
        except Exception as e:
            logger.error(f"检查命令失败|Failed to check command {command}: {str(e)}")
            return False

    @staticmethod
    def check_disk_space(path: str, required_gb: int, logger: logging.Logger) -> bool:
        """检查磁盘空间|Check disk space"""
        try:
            result = subprocess.run(
                ["df", "-k", path],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                available_kb = int(result.stdout.split('\n')[1].split()[3])
                available_gb = available_kb // 1024 // 1024

                if available_gb < required_gb:
                    logger.warning(f"磁盘空间不足|Insufficient disk space: {available_gb}GB available, {required_gb}GB required")
                    return False
                else:
                    logger.info(f"磁盘空间充足|Sufficient disk space: {available_gb}GB available")
                    return True
            else:
                logger.error("无法检查磁盘空间|Cannot check disk space")
                return False
        except Exception as e:
            logger.error(f"检查磁盘空间失败|Failed to check disk space: {str(e)}")
            return False

    @staticmethod
    def check_memory(required_gb: int, logger: logging.Logger) -> bool:
        """检查系统内存|Check system memory"""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        total_mem_kb = int(line.split()[1])
                        total_mem_gb = total_mem_kb // 1024 // 1024

                        if total_mem_gb < required_gb:
                            logger.warning(f"系统内存可能不足|System memory may be insufficient: {total_mem_gb}GB total, {required_gb}GB recommended")
                            return True  # 警告但不阻止执行|Warning but don't block execution
                        else:
                            logger.info(f"系统内存充足|Sufficient system memory: {total_mem_gb}GB total")
                            return True
            logger.error("无法读取内存信息|Cannot read memory information")
            return False
        except Exception as e:
            logger.error(f"检查内存失败|Failed to check memory: {str(e)}")
            return False