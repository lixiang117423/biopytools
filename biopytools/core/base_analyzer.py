"""
基础分析器类 | Base Analyzer Class
标准化分析器基类，包含通用的分析流程和错误处理
"""

import time
import logging
import signal
import sys
from typing import Optional
from .logger import setup_logger
from .config import BaseConfig


class BaseAnalyzer:
    """基础分析器类 | Base Analyzer Class"""

    def __init__(self, config: BaseConfig):
        """
        初始化分析器

        Args:
            config: 配置对象
        """
        self.config = config
        self.start_time = time.time()
        self.interrupted = False

        # 设置日志
        self.logger = setup_logger(
            self.__class__.__name__,
            level=config.get_log_level()
        )

        # 设置信号处理
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _signal_handler(self, signum, frame):
        """信号处理器"""
        self.interrupted = True
        self.logger.warning("进程被信号中断 | Process interrupted by signal")
        sys.exit(130)

    def log_step_start(self, step_name: str, step_number: Optional[int] = None):
        """记录步骤开始"""
        self.logger.info("=" * 60)
        if step_number:
            self.logger.info(f"STEP {step_number}: {step_name}")
        else:
            self.logger.info(f"STEP: {step_name}")
        self.logger.info("=" * 60)

    def log_step_end(self, step_name: str, success: bool = True):
        """记录步骤结束"""
        if success:
            self.logger.info(f"✅ {step_name} completed successfully")
        else:
            self.logger.error(f"❌ {step_name} failed")

    def log_program_info(self, program_name: str):
        """记录程序信息"""
        self.logger.info("=" * 60)
        self.logger.info(f"Program: {program_name}")
        self.logger.info(f"Version: {self.config.version}")
        self.logger.info("=" * 60)

    def log_parameters(self, **kwargs):
        """记录参数信息"""
        for key, value in kwargs.items():
            self.logger.info(f"{key}: {value}")

    def validate_inputs(self) -> bool:
        """
        验证输入参数 - 子类应该重写此方法

        Returns:
            bool: 验证是否通过
        """
        try:
            self.config.validate()
            return True
        except Exception as e:
            self.logger.error(f"配置验证失败 | Configuration validation failed: {e}")
            return False

    def run_analysis(self) -> bool:
        """
        运行分析流程 - 子类必须重写此方法

        Returns:
            bool: 分析是否成功
        """
        raise NotImplementedError("子类必须实现run_analysis方法 | Subclass must implement run_analysis method")

    def run_pipeline(self) -> bool:
        """
        运行完整的分析流程

        Returns:
            bool: 分析是否成功
        """
        try:
            # 记录开始信息
            self.logger.info("Pipeline started")

            # 验证输入
            if not self.validate_inputs():
                self.logger.critical("Input validation failed")
                return False

            # 如果是dry run模式，只显示参数
            if self.config.dry_run:
                self.logger.info("DRY RUN mode - no files will be modified")
                return True

            # 运行分析
            success = self.run_analysis()

            # 记录完成信息
            self._log_completion_summary(success)
            return success

        except KeyboardInterrupt:
            self.logger.warning("Process interrupted by user")
            return False
        except Exception as e:
            self.logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
            return False

    def _log_completion_summary(self, success: bool):
        """记录完成摘要"""
        elapsed_time = time.time() - self.start_time

        self.logger.info("=" * 60)
        if success:
            self.logger.info("Pipeline completed successfully")
        else:
            self.logger.info("Pipeline failed")
        self.logger.info("=" * 60)
        self.logger.info(f"Total runtime: {elapsed_time:.2f} seconds")

        # 子类可以重写此方法添加更多信息
        self.log_additional_summary()

    def log_additional_summary(self):
        """
        记录额外的摘要信息 - 子类可以重写此方法
        """
        pass

    def log_statistics(self, stats: dict):
        """记录统计信息"""
        self.logger.info("Statistics:")
        for key, value in stats.items():
            if isinstance(value, (int, float)):
                # 对于数字，使用千位分隔符
                if isinstance(value, int):
                    formatted_value = f"{value:,}"
                else:
                    formatted_value = f"{value:.3f}"
                self.logger.info(f"  {key}: {formatted_value}")
            else:
                self.logger.info(f"  {key}: {value}")