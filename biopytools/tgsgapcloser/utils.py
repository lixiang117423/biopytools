"""
TGS-GapCloser工具函数模块|TGS-GapCloser Utility Functions Module
"""

import logging
import sys
import subprocess
from pathlib import Path


class TGSGapCloserLogger:
    """TGS-GapCloser日志管理器|TGS-GapCloser Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志系统（符合规范2.3）|Setup logging system (spec 2.3 compliant)"""
        # 删除旧日志|Delete old log
        if self.log_file and Path(self.log_file).exists():
            Path(self.log_file).unlink()

        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = []

        # 文件handler - DEBUG级别|File handler - DEBUG level
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
            handlers.append(file_handler)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
        handlers.append(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
        handlers.append(stderr_handler)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=level,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, output_dir):
        """初始化命令执行器|Initialize command runner"""
        self.logger = logger
        self.output_dir = Path(output_dir)

    def run_command(self, command, check=True, timeout=None, cwd=None):
        """
        执行命令|Execute command

        Args:
            command: 命令列表|Command list
            check: 是否检查返回码|Whether to check return code
            timeout: 超时时间（秒）|Timeout in seconds
            cwd: 工作目录|Working directory

        Returns:
            subprocess.CompletedProcess: 命令执行结果|Command execution result
        """
        self.logger.info(f"执行命令|Executing command: {' '.join(str(c) for c in command)}")

        try:
            result = subprocess.run(
                command,
                check=check,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=cwd
            )

            if result.stdout:
                # 只在DEBUG级别记录stdout|Only log stdout at DEBUG level
                for line in result.stdout.splitlines():
                    self.logger.debug(f"stdout| {line}")

            if result.stderr:
                # stderr记录为WARNING级别|Log stderr as WARNING level
                for line in result.stderr.splitlines():
                    self.logger.warning(f"stderr| {line}")

            self.logger.info(f"命令执行完成|Command execution completed with return code: {result.returncode}")
            return result

        except subprocess.TimeoutExpired as e:
            self.logger.error(f"命令执行超时|Command execution timeout: {e}")
            raise

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed with return code: {e.returncode}")
            if e.stdout:
                self.logger.debug(f"stdout| {e.stdout}")
            if e.stderr:
                self.logger.error(f"stderr| {e.stderr}")
            raise

    def build_tgsgapcloser_command(self, config):
        """
        构建TGS-GapCloser命令|Build TGS-GapCloser command

        Args:
            config: TGSGapCloserConfig对象|TGSGapCloserConfig object

        Returns:
            list: 命令列表|Command list
        """
        cmd = [config.tgsgapcloser_path]

        # 必需参数|Required parameters (使用长参数格式|Use long option format)
        cmd.extend(['--scaff', config.scaff_file])
        cmd.extend(['--reads', config.reads_file])
        cmd.extend(['--output', config.output_prefix])

        # TGS类型|TGS type
        cmd.extend(['--tgstype', config.tgstype])

        # 纠错模式|Error correction mode
        if config.mode == 'none':
            cmd.append('--ne')  # no error correction
        elif config.mode == 'racon':
            if config.racon_path:
                cmd.extend(['--racon', config.racon_path])
        elif config.mode == 'pilon':
            if config.ngs_file:
                cmd.extend(['--ngs', config.ngs_file])
            if config.pilon_path:
                cmd.extend(['--pilon', config.pilon_path])
            if config.samtools_path:
                cmd.extend(['--samtools', config.samtools_path])
            if config.java_path:
                cmd.extend(['--java', config.java_path])

        # 过滤参数|Filter parameters
        if config.min_idy is not None:
            cmd.extend(['--min_idy', str(config.min_idy)])
        if config.min_match is not None:
            cmd.extend(['--min_match', str(config.min_match)])

        # 线程数|Threads
        cmd.extend(['--thread', str(config.threads)])

        # 分块数量|Chunk count
        cmd.extend(['--chunk', str(config.chunk)])

        # Gap大小差异检查|Gap size difference check
        if config.g_check:
            cmd.append('--g_check')

        # 最小/最大reads数量|Min/Max read count
        cmd.extend(['--min_nread', str(config.min_nread)])
        cmd.extend(['--max_nread', str(config.max_nread)])

        # 最大候选数|Max candidates
        cmd.extend(['--max_candidate', str(config.max_candidate)])

        # Racon参数|Racon parameters
        if config.mode == 'racon':
            cmd.extend(['--r_round', str(config.racon_round)])

        # Pilon参数|Pilon parameters
        if config.mode == 'pilon':
            cmd.extend(['--pilon_mem', config.pilon_mem])
            cmd.extend(['--p_round', str(config.pilon_round)])

        # 自定义minimap2参数|Custom minimap2 parameters
        if config.minmap_arg:
            cmd.extend(['--minmap_arg', config.minmap_arg])

        return cmd
