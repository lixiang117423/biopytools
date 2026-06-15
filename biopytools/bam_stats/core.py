"""
BAM统计分析核心模块|BAM Statistics Analysis Core Module
"""

import sys
from pathlib import Path

from .config import BAMStatsConfig
from .batch_processor import BAMBatchProcessor
from .utils import setup_logger, check_dependencies


class BAMAnalyzer:
    """BAM统计分析器|BAM Statistics Analyzer"""

    def __init__(self, config: BAMStatsConfig):
        """初始化BAM分析器|Initialize BAM analyzer"""
        config.validate()
        self.config = config
        self.logger = setup_logger(config.output_dir)

    def run_analysis(self) -> bool:
        """运行BAM统计分析|Run BAM statistics analysis"""
        try:
            self.logger.info(
                "开始BAM统计分析流程|Starting BAM statistics analysis pipeline"
            )

            if not check_dependencies(self.logger):
                return False

            bam_files = self.config.bam_files
            if not bam_files:
                self.logger.error(
                    f"未找到BAM文件|No BAM files found in: "
                    f"{self.config.input_path}"
                )
                return False

            processor = BAMBatchProcessor(self.config, self.logger)
            results = processor.process_bam_files(bam_files)
            if not results:
                return False

            success = processor.generate_report(
                results, self.config.output_file,
                self.config.output_dir, self.config.prefix
            )
            if success:
                self.logger.info(
                    f"分析完成，结果保存至|"
                    f"Analysis completed, results saved to: "
                    f"{self.config.output_dir}"
                )
                return True

            return False

        except KeyboardInterrupt:
            self.logger.info("用户中断操作|User interrupted operation")
            return False
        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            return False
