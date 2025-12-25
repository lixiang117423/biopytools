# """
# FASTP质控主程序模块 | FASTP Quality Control Main Module
# """

# import argparse
# import sys
# from .config import FastpConfig
# from .utils import FastpLogger, CommandRunner
# from .data_processing import SampleFinder
# from .processing import FastpCore
# from .results import SummaryGenerator

# class FastpProcessor:
#     """FASTP质控主类 | Main FASTP Quality Control Class"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置 | Initialize configuration
#         self.config = FastpConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志 | Initialize logging
#         self.logger_manager = FastpLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger)
        
#         # 初始化各个处理器 | Initialize processors
#         self.sample_finder = SampleFinder(self.config, self.logger)
#         self.fastp_core = FastpCore(self.config, self.logger, self.cmd_runner)
#         self.summary_generator = SummaryGenerator(self.config, self.logger)
    
#     def run_batch_processing(self):
#         """运行批处理 | Run batch processing"""
        
#         self.logger.info("="*60)
#         self.logger.info("开始FASTQ数据质控批处理 | Starting FASTQ data quality control batch processing")
#         self.logger.info("="*60)
#         self.logger.info(f"输入目录 | Input directory: {self.config.input_dir}")
#         self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
#         self.logger.info(f"文件模式 | File pattern: *{self.config.read1_suffix}, *{self.config.read2_suffix}")
#         self.logger.info("="*60)
        
#         # 验证fastp可执行性 | Validate fastp executable
#         if not self.fastp_core.validate_fastp():
#             sys.exit(1)
        
#         # 创建输出目录 | Create output directories
#         self.fastp_core.create_output_directories()
        
#         # 查找样本配对 | Find sample pairs
#         sample_pairs = self.sample_finder.find_sample_pairs()
        
#         # 验证样本配对 | Validate sample pairs
#         if not self.sample_finder.validate_sample_pairs(sample_pairs):
#             sys.exit(1)
        
#         self.logger.info(f"找到 {len(sample_pairs)} 个有效样本配对 | Found {len(sample_pairs)} valid sample pairs:")
#         for sample_name, _, _ in sample_pairs:
#             self.logger.info(f"  - {sample_name}")
        
#         # 处理所有样本 | Process all samples
#         successful_count = 0
#         failed_count = 0
        
#         for sample_name, read1_file, read2_file in sample_pairs:
#             if self.fastp_core.process_sample(sample_name, read1_file, read2_file):
#                 successful_count += 1
#             else:
#                 failed_count += 1
        
#         # 生成总结报告 | Generate summary report
#         self.summary_generator.generate_summary_report(
#             successful_count, failed_count, len(sample_pairs)
#         )
        
#         # 输出最终统计 | Output final statistics
#         self.logger.info("="*60)
#         self.logger.info("FASTQ质控批处理完成 | FASTQ quality control batch processing completed!")
#         self.logger.info(f"总样本数 | Total samples: {len(sample_pairs)}")
#         self.logger.info(f"成功处理 | Successfully processed: {successful_count}")
#         self.logger.info(f"失败样本 | Failed samples: {failed_count}")
#         self.logger.info(f"成功率 | Success rate: {(successful_count/len(sample_pairs))*100:.1f}%")
#         self.logger.info(f"质控后的清洁数据位于 | Clean data location: {self.config.output_dir}")
#         self.logger.info(f"质控报告位于 | QC reports location: {self.config.report_path}")
#         self.logger.info("="*60)

# def main():
#     """主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description="FASTQ数据质控批处理脚本 | FASTQ Data Quality Control Batch Processing Script",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter
#     )
    
#     # 必需参数 | Required arguments
#     parser.add_argument("-i", "--input-dir", required=True,
#                        help="输入原始FASTQ数据目录 | Input raw FASTQ data directory")
#     parser.add_argument("-o", "--output-dir", required=True,
#                        help="输出清洁FASTQ数据目录 | Output clean FASTQ data directory")
    
#     # 可选参数 | Optional arguments
#     parser.add_argument("--fastp-path", default="fastp",
#                        help="fastp可执行文件路径 | fastp executable path")
#     parser.add_argument("-t", "--threads", type=int, default=12,
#                        help="线程数 | Number of threads")
#     parser.add_argument("-q", "--quality-threshold", type=int, default=30,
#                        help="质量阈值 | Quality threshold")
#     parser.add_argument("-l", "--min-length", type=int, default=50,
#                        help="最小长度 | Minimum length")
#     parser.add_argument("-u", "--unqualified-percent", type=int, default=40,
#                        help="不合格碱基百分比阈值 | Unqualified base percentage threshold")
#     parser.add_argument("-n", "--n-base-limit", type=int, default=10,
#                        help="N碱基数量限制 | N base count limit")
#     parser.add_argument("--read1-suffix", default="_1.fq.gz",
#                        help="Read1文件后缀 | Read1 file suffix")
#     parser.add_argument("--read2-suffix", default="_2.fq.gz",
#                        help="Read2文件后缀 | Read2 file suffix")
    
#     args = parser.parse_args()
    
#     # 创建处理器并运行 | Create processor and run
#     processor = FastpProcessor(
#         input_dir=args.input_dir,
#         output_dir=args.output_dir,
#         fastp_path=args.fastp_path,
#         threads=args.threads,
#         quality_threshold=args.quality_threshold,
#         min_length=args.min_length,
#         unqualified_percent=args.unqualified_percent,
#         n_base_limit=args.n_base_limit,
#         read1_suffix=args.read1_suffix,
#         read2_suffix=args.read2_suffix
#     )
    
#     processor.run_batch_processing()

# if __name__ == "__main__":
#     main()

"""
🚀 FASTP质控主程序模块 | FASTP Quality Control Main Module
"""

import argparse
import sys
from .config import FastpConfig
from .utils import FastpLogger, CommandRunner
from .data_processing import SampleFinder
from .processing import FastpCore
from .results import SummaryGenerator

class FastpProcessor:
    """🎯 FASTP质控主类 | Main FASTP Quality Control Class"""
    
    def __init__(self, **kwargs):
        # ⚙️ 初始化配置 | Initialize configuration
        self.config = FastpConfig(**kwargs)
        self.config.validate()
        
        # 📝 初始化日志 | Initialize logging
        self.logger_manager = FastpLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 💻 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)
        
        # 🔧 初始化各个处理器 | Initialize processors
        self.sample_finder = SampleFinder(self.config, self.logger)
        self.fastp_core = FastpCore(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_batch_processing(self):
        """🚀 运行批处理 | Run batch processing"""
        
        self.logger.info("="*60)
        self.logger.info("🚀 开始FASTQ数据质控批处理 | Starting FASTQ data quality control batch processing")
        self.logger.info("="*60)
        self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
        self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"📄 文件模式 | File pattern: *{self.config.read1_suffix}, *{self.config.read2_suffix}")
        self.logger.info("="*60)
        
        # ✅ 验证fastp可执行性 | Validate fastp executable
        if not self.fastp_core.validate_fastp():
            sys.exit(1)
        
        # 📁 创建输出目录 | Create output directories
        self.fastp_core.create_output_directories()
        
        # 🔍 查找样本配对 | Find sample pairs
        sample_pairs = self.sample_finder.find_sample_pairs()
        
        # ✅ 验证样本配对 | Validate sample pairs
        if not self.sample_finder.validate_sample_pairs(sample_pairs):
            sys.exit(1)
        
        self.logger.info(f"🧬 找到 {len(sample_pairs)} 个有效样本配对 | Found {len(sample_pairs)} valid sample pairs:")
        for sample_name, _, _ in sample_pairs:
            self.logger.info(f"  - 📋 {sample_name}")
        
        # ⚡ 处理所有样本 | Process all samples
        successful_count = 0
        failed_count = 0
        
        for sample_name, read1_file, read2_file in sample_pairs:
            if self.fastp_core.process_sample(sample_name, read1_file, read2_file):
                successful_count += 1
            else:
                failed_count += 1
        
        # 📊 生成总结报告 | Generate summary report
        self.summary_generator.generate_summary_report(
            successful_count, failed_count, len(sample_pairs)
        )
        
        # 📈 输出最终统计 | Output final statistics
        self.logger.info("="*60)
        self.logger.info("🎉 FASTQ质控批处理完成 | FASTQ quality control batch processing completed!")
        self.logger.info(f"🧮 总样本数 | Total samples: {len(sample_pairs)}")
        self.logger.info(f"✅ 成功处理 | Successfully processed: {successful_count}")
        self.logger.info(f"❌ 失败样本 | Failed samples: {failed_count}")
        self.logger.info(f"📊 成功率 | Success rate: {(successful_count/len(sample_pairs))*100:.1f}%")
        self.logger.info(f"🧬 质控后的清洁数据位于 | Clean data location: {self.config.output_dir}")
        self.logger.info(f"📋 质控报告位于 | QC reports location: {self.config.report_path}")
        self.logger.info("="*60)

def main():
    """🚀 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🧬 FASTQ数据质控批处理脚本 | FASTQ Data Quality Control Batch Processing Script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 📋 必需参数 | Required arguments
    parser.add_argument("-i", "--input-dir", required=True,
                       help="📂 输入原始FASTQ数据目录 | Input raw FASTQ data directory")
    parser.add_argument("-o", "--output-dir", required=True,
                       help="📁 输出清洁FASTQ数据目录 | Output clean FASTQ data directory")
    
    # ⚙️ 可选参数 | Optional arguments
    parser.add_argument("--fastp-path", default="fastp",
                       help="🛠️ fastp可执行文件路径 | fastp executable path")
    parser.add_argument("-t", "--threads", type=int, default=12,
                       help="🧵 线程数 | Number of threads")
    parser.add_argument("-q", "--quality-threshold", type=int, default=30,
                       help="🎯 质量阈值 | Quality threshold")
    parser.add_argument("-l", "--min-length", type=int, default=50,
                       help="📏 最小长度 | Minimum length")
    parser.add_argument("-u", "--unqualified-percent", type=int, default=40,
                       help="📊 不合格碱基百分比阈值 | Unqualified base percentage threshold")
    parser.add_argument("-n", "--n-base-limit", type=int, default=10,
                       help="🧬 N碱基数量限制 | N base count limit")
    parser.add_argument("--read1-suffix", default="_1.fq.gz",
                       help="📄 Read1文件后缀 | Read1 file suffix")
    parser.add_argument("--read2-suffix", default="_2.fq.gz",
                       help="📄 Read2文件后缀 | Read2 file suffix")
    
    args = parser.parse_args()
    
    # 🏗️ 创建处理器并运行 | Create processor and run
    processor = FastpProcessor(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        fastp_path=args.fastp_path,
        threads=args.threads,
        quality_threshold=args.quality_threshold,
        min_length=args.min_length,
        unqualified_percent=args.unqualified_percent,
        n_base_limit=args.n_base_limit,
        read1_suffix=args.read1_suffix,
        read2_suffix=args.read2_suffix
    )
    
    processor.run_batch_processing()

if __name__ == "__main__":
    main()