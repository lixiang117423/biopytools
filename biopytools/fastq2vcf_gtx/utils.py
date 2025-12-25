"""
Fastq到VCF (GTX) 工具模块 | Fastq to VCF (GTX) Utilities Module
"""

import os
import subprocess
import time
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any
import signal


class Fastq2VcfGTXLogger:
    """Fastq到VCF (GTX) 日志管理器 | Fastq to VCF (GTX) Logger Manager"""

    def __init__(self, log_dir: str, verbose: bool = False):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose

        # 创建日志文件 | Create log files
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.log_dir / f"pipeline_{timestamp}.log"
        self.error_log = self.log_dir / f"error_{timestamp}.log"

        # 配置logger | Configure logger
        self.logger = logging.getLogger(f"fastq2vcf_gtx_{timestamp}")
        self.logger.setLevel(logging.DEBUG if self.verbose else logging.INFO)

        # 清除现有的处理器 | Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 文件处理器 | File handler
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter('[%(levelname)s] %(asctime)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)

        # 控制台处理器 - 确保始终输出到stdout | Console handler - ensure always output to stdout
        import sys
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)  # 使用DEBUG级别以捕获所有输出
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)

        # 错误日志处理器 | Error log handler
        error_handler = logging.FileHandler(self.error_log)
        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(file_formatter)
        self.logger.addHandler(error_handler)

    def get_logger(self):
        """获取logger实例 | Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤 | Log step"""
        self.logger.info("\n" + "="*50)
        self.logger.info(message)
        self.logger.info("="*50)


class CommandRunner:
    """命令执行器 | Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str, dry_run: bool = False):
        self.logger = logger
        self.output_dir = output_dir
        self.dry_run = dry_run

    def run(self, command, description: str = "") -> bool:
        """
        运行命令 | Run command

        Args:
            command: 要执行的命令（字符串或列表）| Command to execute (string or list)
            description: 命令描述 | Command description

        Returns:
            bool: 是否成功 | Success status
        """
        if description:
            self.logger.info(f"{description}")

        # 格式化命令显示 | Format command display
        if isinstance(command, list):
            command_str = ' '.join(command)
        else:
            command_str = command

        self.logger.info(f"🔄 执行命令 | Running command: {command_str}")

        if self.dry_run:
            self.logger.info(f"[DRY RUN] 跳过实际执行 | Skip actual execution")
            return True

        try:
            # 执行命令 | Execute command
            if isinstance(command, list):
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=7200,  # 2小时超时 | 2 hours timeout
                    cwd=self.output_dir
                )
            else:
                result = subprocess.run(
                    command,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=7200,  # 2小时超时 | 2 hours timeout
                    cwd=self.output_dir
                )

            if result.returncode == 0:
                # 总是输出标准输出，即使为空 | Always output stdout, even if empty
                if result.stdout:
                    # 输出标准内容到INFO级别，以便在.out文件中显示
                    self.logger.info(f"标准输出 | Stdout: {result.stdout}")
                else:
                    self.logger.info("标准输出 | Stdout: (空输出)")

                # 输出stderr信息（如果有）| Output stderr if any
                if result.stderr:
                    self.logger.info(f"标准错误 | Stderr: {result.stderr}")

                # 检查并读取biopytools fastp的日志文件 | Check and read biopytools fastp log file
                if "biopytools fastp" in command_str if isinstance(command, str) else "fastp" in " ".join(command) if isinstance(command, list) else "":
                    # 尝试找到fastp日志文件 | Try to find fastp log file
                    fastp_log_files = []

                    # 检查输出目录中的fastp报告目录 | Check fastp reports directory in output directory
                    reports_dir = os.path.join(self.output_dir, "fastp_reports") if self.output_dir else None
                    if reports_dir and os.path.exists(reports_dir):
                        for log_file in os.listdir(reports_dir):
                            if log_file.endswith(".log"):
                                fastp_log_files.append(os.path.join(reports_dir, log_file))

                    # 检查输出目录中的日志文件 | Check log files in output directory
                    if self.output_dir and os.path.exists(self.output_dir):
                        for log_file in os.listdir(self.output_dir):
                            if log_file.endswith("_processing.log") or log_file.endswith("fastp.log"):
                                fastp_log_files.append(os.path.join(self.output_dir, log_file))

                    # 显示找到的日志文件内容 | Show contents of found log files
                    for log_file in fastp_log_files:
                        if os.path.exists(log_file):
                            self.logger.info(f"📄 FastP日志文件内容 | FastP Log File Content: {log_file}")
                            try:
                                with open(log_file, 'r', encoding='utf-8') as f:
                                    log_content = f.read()
                                    if log_content:
                                        # 分行显示日志内容 | Display log content line by line
                                        for line in log_content.split('\n'):
                                            if line.strip():
                                                self.logger.info(f"  {line}")
                                    else:
                                        self.logger.info("  📄 (空日志文件) | (Empty log file)")
                            except Exception as e:
                                self.logger.error(f"❌ 读取FastP日志文件失败 | Failed to read FastP log file {log_file}: {str(e)}")
                        else:
                            self.logger.warning(f"⚠️ FastP日志文件不存在 | FastP log file does not exist: {log_file}")

                self.logger.info("✅ 命令执行成功 | Command executed successfully")
                return True
            else:
                self.logger.error(f"❌ 命令执行失败 | Command execution failed")
                self.logger.error(f"返回码 | Return code: {result.returncode}")
                if result.stdout:
                    self.logger.error(f"标准输出 | Stdout: {result.stdout}")
                if result.stderr:
                    self.logger.error(f"错误输出 | Stderr: {result.stderr}")

                # 失败时也尝试读取fastp日志文件 | Also try to read fastp log file on failure
                if "biopytools fastp" in command_str if isinstance(command, str) else "fastp" in " ".join(command) if isinstance(command, list) else "":
                    self.logger.info("🔍 尝试读取FastP日志文件以获取更多错误信息 | Trying to read FastP log file for more error information")
                    fastp_log_files = []

                    # 检查输出目录中的fastp报告目录 | Check fastp reports directory in output directory
                    reports_dir = os.path.join(self.output_dir, "fastp_reports") if self.output_dir else None
                    if reports_dir and os.path.exists(reports_dir):
                        for log_file in os.listdir(reports_dir):
                            if log_file.endswith(".log"):
                                fastp_log_files.append(os.path.join(reports_dir, log_file))

                    # 检查输出目录中的日志文件 | Check log files in output directory
                    if self.output_dir and os.path.exists(self.output_dir):
                        for log_file in os.listdir(self.output_dir):
                            if log_file.endswith("_processing.log") or log_file.endswith("fastp.log"):
                                fastp_log_files.append(os.path.join(self.output_dir, log_file))

                    # 显示找到的日志文件内容 | Show contents of found log files
                    for log_file in fastp_log_files:
                        if os.path.exists(log_file):
                            self.logger.error(f"📄 FastP错误日志文件内容 | FastP Error Log File Content: {log_file}")
                            try:
                                with open(log_file, 'r', encoding='utf-8') as f:
                                    log_content = f.read()
                                    if log_content:
                                        # 分行显示日志内容 | Display log content line by line
                                        for line in log_content.split('\n'):
                                            if line.strip():
                                                self.logger.error(f"  {line}")
                                    else:
                                        self.logger.error("  📄 (空日志文件) | (Empty log file)")
                            except Exception as e:
                                self.logger.error(f"❌ 读取FastP日志文件失败 | Failed to read FastP log file {log_file}: {str(e)}")
                        else:
                            self.logger.warning(f"⚠️ FastP日志文件不存在 | FastP log file does not exist: {log_file}")

                return False

        except subprocess.TimeoutExpired:
            self.logger.error("❌ 命令执行超时 | Command execution timed out")
            return False
        except Exception as e:
            self.logger.error(f"❌ 命令执行异常 | Command execution exception: {str(e)}")
            return False

    def run_with_progress(self, command: str, description: str = "",
                         timeout: int = 7200) -> bool:
        """
        运行命令并显示进度 | Run command with progress indication

        Args:
            command: 要执行的命令 | Command to execute
            description: 命令描述 | Command description
            timeout: 超时时间（秒）| Timeout in seconds

        Returns:
            bool: 是否成功 | Success status
        """
        if description:
            self.logger.info(f"{description}")

        self.logger.info(f"🔄 执行命令 | Running command: {command}")

        if self.dry_run:
            self.logger.info(f"[DRY RUN] 跳过实际执行 | Skip actual execution")
            return True

        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=self.output_dir
            )

            start_time = time.time()

            while True:
                # 检查进程是否结束 | Check if process has finished
                if process.poll() is not None:
                    break

                # 检查超时 | Check timeout
                if time.time() - start_time > timeout:
                    process.terminate()
                    process.wait()
                    self.logger.error("❌ 命令执行超时 | Command execution timed out")
                    return False

                # 读取输出 | Read output
                try:
                    output = process.stdout.readline()
                    if output:
                        # 实时输出到INFO级别，以便在.out文件中显示进度
                        self.logger.info(output.strip())
                except:
                    pass

                time.sleep(1)

            # 获取最终输出 | Get final output
            remaining_output, _ = process.communicate()
            if remaining_output:
                # 输出最终结果到INFO级别
                self.logger.info(remaining_output.strip())

            if process.returncode == 0:
                self.logger.info("✅ 命令执行成功 | Command executed successfully")
                return True
            else:
                self.logger.error(f"❌ 命令执行失败 | Command execution failed")
                self.logger.error(f"返回码 | Return code: {process.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"❌ 命令执行异常 | Command execution exception: {str(e)}")
            return False


class CheckpointManager:
    """检查点管理器 | Checkpoint Manager"""

    def __init__(self, checkpoint_dir: str, logger: logging.Logger):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger

    def exists(self, step: str) -> bool:
        """检查检查点是否存在 | Check if checkpoint exists"""
        return (self.checkpoint_dir / f"{step}.done").exists()

    def create(self, step: str):
        """创建检查点 | Create checkpoint"""
        checkpoint_file = self.checkpoint_dir / f"{step}.done"
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        checkpoint_file.write_text(timestamp)
        self.logger.info(f"📍 创建检查点 | Create checkpoint: {step}")

    def remove(self, step: str):
        """移除检查点 | Remove checkpoint"""
        checkpoint_file = self.checkpoint_dir / f"{step}.done"
        if checkpoint_file.exists():
            checkpoint_file.unlink()
            self.logger.info(f"🗑️ 移除检查点 | Remove checkpoint: {step}")

    def list_completed(self) -> List[str]:
        """列出已完成的检查点 | List completed checkpoints"""
        completed = []
        for checkpoint_file in self.checkpoint_dir.glob("*.done"):
            step = checkpoint_file.stem
            timestamp = checkpoint_file.read_text()
            completed.append(f"{step} ({timestamp})")
        return completed


class FileManager:
    """文件管理器 | File Manager"""

    @staticmethod
    def count_files(directory: str, pattern: str = "*") -> int:
        """统计目录中文件数量 | Count files in directory"""
        if not os.path.exists(directory):
            return 0
        return len(list(Path(directory).glob(pattern)))

    @staticmethod
    def find_files(directory: str, pattern: str) -> List[str]:
        """查找文件 | Find files"""
        if not os.path.exists(directory):
            return []
        return [str(f) for f in Path(directory).glob(pattern)]

    @staticmethod
    def ensure_directory(directory: str):
        """确保目录存在 | Ensure directory exists"""
        Path(directory).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_file_size(file_path: str) -> str:
        """获取文件大小 | Get file size"""
        if not os.path.exists(file_path):
            return "0 B"
        size = os.path.getsize(file_path)
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size < 1024:
                return f"{size:.1f} {unit}"
            size /= 1024
        return f"{size:.1f} PB"


class SystemChecker:
    """系统检查器 | System Checker"""

    @staticmethod
    def check_command_exists(command: str, logger: logging.Logger) -> bool:
        """检查命令是否存在 | Check if command exists"""
        try:
            result = subprocess.run(
                ["which", command],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                logger.info(f"✓ 命令存在 | Command exists: {command}")
                return True
            else:
                logger.error(f"❌ 命令不存在 | Command does not exist: {command}")
                return False
        except Exception as e:
            logger.error(f"❌ 检查命令失败 | Failed to check command {command}: {str(e)}")
            return False

    @staticmethod
    def check_disk_space(path: str, required_gb: int, logger: logging.Logger) -> bool:
        """检查磁盘空间 | Check disk space"""
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
                    logger.warning(f"⚠️ 磁盘空间不足 | Insufficient disk space: {available_gb}GB available, {required_gb}GB required")
                    return False
                else:
                    logger.info(f"✓ 磁盘空间充足 | Sufficient disk space: {available_gb}GB available")
                    return True
            else:
                logger.error("❌ 无法检查磁盘空间 | Cannot check disk space")
                return False
        except Exception as e:
            logger.error(f"❌ 检查磁盘空间失败 | Failed to check disk space: {str(e)}")
            return False

    @staticmethod
    def check_memory(required_gb: int, logger: logging.Logger) -> bool:
        """检查系统内存 | Check system memory"""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        total_mem_kb = int(line.split()[1])
                        total_mem_gb = total_mem_kb // 1024 // 1024

                        if total_mem_gb < required_gb:
                            logger.warning(f"⚠️ 系统内存可能不足 | System memory may be insufficient: {total_mem_gb}GB total, {required_gb}GB recommended")
                            return True  # 警告但不阻止执行 | Warning but don't block execution
                        else:
                            logger.info(f"✓ 系统内存充足 | Sufficient system memory: {total_mem_gb}GB total")
                            return True
            logger.error("❌ 无法读取内存信息 | Cannot read memory information")
            return False
        except Exception as e:
            logger.error(f"❌ 检查内存失败 | Failed to check memory: {str(e)}")
            return False