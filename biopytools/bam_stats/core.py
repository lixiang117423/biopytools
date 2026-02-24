"""
BAM统计分析核心模块|BAM Statistics Analysis Core Module
"""

import os
import subprocess
import sys
from pathlib import Path
from multiprocessing import cpu_count

from .batch_processor import BAMBatchProcessor
from .utils import BAMStatsLogger, check_dependencies


class BAMAnalyzerConfig:
    """BAM分析器配置类|BAM Analyzer Configuration Class"""

    def __init__(self, **kwargs):
        self.input_dir = Path(kwargs.get('input_dir', '')).expanduser()
        self.output_file = Path(kwargs.get('output_file', '')).expanduser()
        self.processes = kwargs.get('processes', cpu_count())
        self.log_dir = Path(kwargs.get('log_dir', '.')).expanduser()

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.input_dir:
            raise ValueError("输入路径不能为空|Input path cannot be empty")

        if not self.output_file:
            raise ValueError("输出文件不能为空|Output file cannot be empty")

        input_path = self.input_dir
        if not input_path.exists():
            raise FileNotFoundError(f"输入路径不存在|Input path does not exist: {input_path}")

        # 允许文件或目录|Allow file or directory
        if not (input_path.is_dir() or input_path.is_file()):
            raise ValueError(f"输入路径必须是文件或目录|Input path must be a file or directory: {input_path}")

        # 确保输出目录存在|Ensure output directory exists
        output_path = self.output_file
        output_path.parent.mkdir(parents=True, exist_ok=True)


class BAMAnalyzer:
    """BAM统计分析器|BAM Statistics Analyzer"""

    def __init__(self, **kwargs):
        """初始化BAM分析器|Initialize BAM analyzer"""
        self.config = BAMAnalyzerConfig(**kwargs)
        self.config.validate()

        # 设置日志|Setup logging
        self.logger_manager = BAMStatsLogger(self.config.log_dir, "bam_stats_analysis.log")
        self.logger = self.logger_manager.get_logger()

        # 初始化处理器|Initialize processor
        self.processor = BAMBatchProcessor(self.config, self.logger)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.logger)

    def find_bam_files(self):
        """查找输入目录或单个BAM文件|Find BAM files in input directory or single BAM file"""
        input_path = Path(self.config.input_dir)

        # 如果是单个BAM文件，直接返回|If single BAM file, return directly
        if input_path.is_file():
            if input_path.suffix.lower() != '.bam':
                self.logger.error(f"输入文件不是.bam格式|Input file is not .bam format: {input_path}")
                return []
            return [input_path]

        # 如果是目录，查找所有.bam文件|If directory, find all .bam files
        bam_files = sorted(list(input_path.glob('*.bam')))
        return bam_files

    def validate_bam_files(self, bam_files):
        """验证BAM文件|Validate BAM files"""
        valid_files = []
        for bam_file in bam_files:
            if self._validate_single_bam_file(bam_file):
                valid_files.append(bam_file)

        return valid_files

    def _validate_single_bam_file(self, bam_file: Path):
        """验证单个BAM文件|Validate single BAM file"""
        if not bam_file.exists():
            self.logger.error(f"文件不存在|File does not exist: {bam_file}")
            return False

        if not bam_file.suffix.lower() in ['.bam']:
            self.logger.error(f"文件格式错误，应为.bam文件|File format error, should be .bam file: {bam_file}")
            return False

        # 检查BAM文件是否完整|Check if BAM file is complete
        try:
            result = subprocess.run(['samtools', 'quickcheck', str(bam_file)],
                                  capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.warning(f"BAM文件可能不完整|BAM file may be incomplete: {bam_file}")
                self.logger.warning(f"   {result.stderr}")
                # 不返回False，继续尝试处理|Don't return False, try to process anyway
        except FileNotFoundError:
            self.logger.warning("samtools quickcheck不可用，跳过BAM文件完整性检查|samtools quickcheck not available, skipping BAM file integrity check")

        return True

    def run_analysis(self):
        """运行BAM统计分析|Run BAM statistics analysis"""
        try:
            self.logger.info("开始BAM统计分析流程|Starting BAM statistics analysis pipeline")

            # 检查依赖|Check dependencies
            if not self.check_dependencies():
                return False

            # 查找并验证BAM文件|Find and validate BAM files
            bam_files = self.find_bam_files()
            if not bam_files:
                self.logger.error(f"在目录中没有找到任何.bam文件|No .bam files found in directory: '{self.config.input_dir}'")
                return False

            valid_bam_files = self.validate_bam_files(bam_files)
            if not valid_bam_files:
                self.logger.error("没有有效的BAM文件可处理|No valid BAM files to process")
                return False

            # 批量处理BAM文件|Batch process BAM files
            results = self.processor.process_bam_files(valid_bam_files, self.config.processes)
            if not results:
                return False

            # 生成报告|Generate report
            success = self.processor.generate_report(results, self.config.output_file)
            if success:
                self.logger.info(f"报告生成完毕|Report generation completed! 结果已保存至|Results saved to: {self.config.output_file}")
                return True
            else:
                return False

        except KeyboardInterrupt:
            self.logger.info("\n用户中断操作|User interrupted operation")
            return False
        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error occurred during analysis: {str(e)}")
            return False
