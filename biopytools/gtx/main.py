# """
# 🚀 GTX WGS分析主程序模块 | GTX WGS Analysis Main Module 🚀
# """

# import argparse
# import sys
# from .config import GTXConfig
# from .utils import GTXLogger, CommandRunner, check_dependencies, FileProcessor
# from .processor import GTXProcessor
# from .results import SummaryGenerator

# class GTXAnalyzer:
#     """🧬 GTX WGS分析主类 | Main GTX WGS Analyzer Class"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置 ⚙️ | Initialize configuration
#         self.config = GTXConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志 📝 | Initialize logging
#         self.logger_manager = GTXLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 ▶️ | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#         # 初始化各个处理器 🛠️ | Initialize processors
#         self.file_processor = FileProcessor(self.config, self.logger)
#         self.gtx_processor = GTXProcessor(self.config, self.logger, self.cmd_runner)
#         self.summary_generator = SummaryGenerator(self.config, self.logger)
    
#     def check_dependencies(self):
#         """检查依赖软件 🧩 | Check dependencies"""
#         return check_dependencies(self.config, self.logger)
    
#     def run_analysis(self):
#         """运行完整的GTX WGS分析流程 🚀 | Run complete GTX WGS analysis pipeline"""
#         try:
#             self.logger.info("🚀 开始GTX WGS批处理分析流程 | Starting GTX WGS batch processing analysis pipeline")
#             self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
#             self.logger.info(f"📤 输出目录 | Output directory: {self.config.output_dir}")
#             self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
#             self.logger.info(f"📁 R1文件模式 | R1 file pattern: {self.config.read1_pattern}")
#             self.logger.info(f"📁 R2文件模式 | R2 file pattern: {self.config.read2_pattern}")
            
#             # 检查依赖 🧩 | Check dependencies
#             self.check_dependencies()
            
#             # 查找输入文件 🔍 | Find input files
#             file_pairs = self.file_processor.find_fastq_files()
            
#             if not file_pairs:
#                 raise RuntimeError("❌ 没有找到需要处理的FASTQ文件对 | No FASTQ file pairs found to process")
            
#             # 处理每个样品 🔬 | Process each sample
#             self.logger.info("=" * 60)
#             self.logger.info("🔬 开始批量处理样品 | Starting batch processing of samples")
            
#             total_samples = len(file_pairs)
#             successful_samples = 0
#             failed_samples = []
            
#             for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
#                 self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing sample: {sample_name}")
                
#                 success = self.gtx_processor.process_sample(sample_name, r1_file, r2_file)
                
#                 if success:
#                     successful_samples += 1
#                     self.logger.info(f"  ✅ 样品 {sample_name} 分析完成 ({current}/{total_samples}) | Sample {sample_name} analysis completed ({current}/{total_samples})")
#                 else:
#                     failed_samples.append(sample_name)
                
#                 self.logger.info("-" * 40)
            
#             # 输出处理结果统计 📊 | Output processing results statistics
#             self.logger.info("=" * 60)
#             self.logger.info("🎉 批量处理完成 | Batch processing completed")
#             self.logger.info(f"✅ 成功处理样品数 | Successfully processed samples: {successful_samples}")
#             self.logger.info(f"❌ 失败样品数 | Failed samples: {len(failed_samples)}")
            
#             if failed_samples:
#                 self.logger.warning(f"⚠️ 失败样品列表 | Failed samples list: {', '.join(failed_samples)}")
            
#             # 生成总结报告 📄 | Generate summary report
#             self.summary_generator.generate_summary_report()
            
#             # 显示最终目录结构 🌳 | Display final directory structure
#             self.logger.info("")
#             self.logger.info("🌳 最终输出目录结构 | Final output directory structure:")
#             self.logger.info(f"  {self.config.output_dir}/")
#             self.logger.info(f"  ├── bam/     (💾 BAM文件: {len(list(self.config.bam_output_dir.glob('*.sorted.bam')))} 个)")
#             self.logger.info(f"  ├── vcf/     (📜 VCF文件: {len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))} 个)")
#             self.logger.info(f"  └── tmp/     (🗑️ 临时文件目录)")
            
#             self.logger.info("=" * 60)
#             self.logger.info("🎉 GTX WGS批处理分析全部完成 | GTX WGS batch processing analysis completed")
#             self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
#         except Exception as e:
#             self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
#             sys.exit(1)

# def main():
#     """主函数 🚀 | Main function"""
#     parser = argparse.ArgumentParser(
#         description='🧬 GTX WGS批处理分析脚本 (模块化版本) | GTX WGS Batch Processing Analysis Script (Modular Version)',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         epilog="""
# 💡 示例 | Examples:
#   %(prog)s -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
#   %(prog)s -i ./clean_data -o ./gtx_results -r ./reference.fa -t 64
#   %(prog)s -i /data/clean -o /results -r /genome/ref.fa --gtx-path /custom/gtx/path
#   %(prog)s -i ./data -o ./results -r ./ref.fa --min-confidence 25 --min-base-quality 15
#   %(prog)s -i ./data -o ./results -r ./ref.fa --read1-pattern "*_1.fq.gz" --read2-pattern "*_2.fq.gz"
#         """
#     )
    
#     # 必需参数 📌 | Required arguments
#     parser.add_argument('-i', '--input-dir', required=True, 
#                        help='输入目录路径 (包含clean FASTQ文件) 📂 | Input directory path (containing clean FASTQ files)')
#     parser.add_argument('-o', '--output-dir', required=True, 
#                        help='输出目录路径 📤 | Output directory path')
#     parser.add_argument('-r', '--reference', required=True, 
#                        help='参考基因组文件路径 🧬 | Reference genome file path')
    
#     # 可选参数 ✨ | Optional arguments
#     parser.add_argument('-t', '--threads', type=int, default=88, 
#                        help='线程数 🧵 | Number of threads')
#     parser.add_argument('--gtx-path', default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx', 
#                        help='GTX程序路径 💻 | GTX program path')
#     parser.add_argument('--tmp-dir', 
#                        help='临时目录路径 (默认使用输出目录下的tmp) 🗑️ | Temporary directory path (default: tmp under output directory)')
    
#     # 质量控制参数 🔬 | Quality control parameters
#     parser.add_argument('--min-confidence', type=int, default=30, 
#                        help='最小置信度阈值 🎯 | Minimum confidence threshold')
#     parser.add_argument('--min-base-quality', type=int, default=20, 
#                        help='最小碱基质量阈值 ✨ | Minimum base quality threshold')
#     parser.add_argument('--ploidy', type=int, default=2, 
#                        help='倍性 🧬 | Ploidy')
#     parser.add_argument('--pcr-indel-model', default='CONSERVATIVE', 
#                        help='PCR indel模型 🔬 | PCR indel model')
    
#     # 文件模式参数 📁 | File pattern parameters
#     parser.add_argument('--read1-pattern', default='*_1.fq.gz', 
#                        help='R1文件匹配模式 📄 | R1 file pattern')
#     parser.add_argument('--read2-pattern', default='*_2.fq.gz', 
#                        help='R2文件匹配模式 📄 | R2 file pattern')
    
#     args = parser.parse_args()
    
#     # 创建分析器并运行 🚀 | Create analyzer and run
#     analyzer = GTXAnalyzer(
#         input_dir=args.input_dir,
#         output_dir=args.output_dir,
#         reference=args.reference,
#         gtx_path=args.gtx_path,
#         threads=args.threads,
#         tmp_dir=args.tmp_dir,
#         min_confidence_threshold=args.min_confidence,
#         min_base_quality=args.min_base_quality,
#         ploidy=args.ploidy,
#         pcr_indel_model=args.pcr_indel_model,
#         read1_pattern=args.read1_pattern,
#         read2_pattern=args.read2_pattern
#     )
    
#     analyzer.run_analysis()

# if __name__ == "__main__":
#     main()


"""
🚀 GTX WGS分析主程序模块 | GTX WGS Analysis Main Module 🚀
"""

import argparse
import sys
from .config import GTXConfig
from .utils import GTXLogger, CommandRunner, check_dependencies, FileProcessor
from .processor import GTXProcessor, JointProcessor
from .results import SummaryGenerator
from .utils import ReferenceIndexManager

class GTXAnalyzer:
    """🧬 GTX WGS分析主类 | Main GTX WGS Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 ⚙️ | Initialize configuration
        self.config = GTXConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 📝 | Initialize logging
        self.logger_manager = GTXLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 ▶️ | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 🛠️ | Initialize processors
        self.file_processor = FileProcessor(self.config, self.logger)
        self.gtx_processor = GTXProcessor(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

        # 构建索引
        self.reference_manager = ReferenceIndexManager(self.config, self.logger, self.cmd_runner)
        
        # 如果启用joint calling，初始化joint处理器
        if self.config.enable_joint_calling:
            self.joint_processor = JointProcessor(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """检查依赖软件 🧩 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def check_and_build_reference_index(self):
        """检查并构建参考基因组索引 🧬 | Check and build reference genome index"""
        return self.reference_manager.ensure_reference_index()

    def run_individual_analysis(self) -> tuple:
        """运行单样品分析 🔬 | Run individual sample analysis"""
        self.logger.info("🔬 开始单样品WGS分析阶段 | Starting individual WGS analysis phase")
        
        # 查找输入文件 🔍 | Find input files
        file_pairs = self.file_processor.find_fastq_files()
        
        if not file_pairs:
            raise RuntimeError("❌ 没有找到需要处理的FASTQ文件对 | No FASTQ file pairs found to process")
        
        # 处理每个样品 🔬 | Process each sample
        self.logger.info("=" * 60)
        self.logger.info("🔬 开始批量处理样品 | Starting batch processing of samples")
        
        total_samples = len(file_pairs)
        successful_samples = 0
        failed_samples = []
        processed_sample_names = []
        
        for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
            self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing sample: {sample_name}")
            
            success = self.gtx_processor.process_sample(sample_name, r1_file, r2_file)
            
            if success:
                successful_samples += 1
                processed_sample_names.append(sample_name)
                self.logger.info(f"  ✅ 样品 {sample_name} 分析完成 ({current}/{total_samples}) | Sample {sample_name} analysis completed ({current}/{total_samples})")
            else:
                failed_samples.append(sample_name)
            
            self.logger.info("-" * 40)
        
        # 输出处理结果统计 📊 | Output processing results statistics
        self.logger.info("=" * 60)
        self.logger.info("🎉 单样品分析阶段完成 | Individual analysis phase completed")
        self.logger.info(f"✅ 成功处理样品数 | Successfully processed samples: {successful_samples}")
        self.logger.info(f"❌ 失败样品数 | Failed samples: {len(failed_samples)}")
        
        if failed_samples:
            self.logger.warning(f"⚠️ 失败样品列表 | Failed samples list: {', '.join(failed_samples)}")
        
        return successful_samples, failed_samples, processed_sample_names
    
    def run_joint_analysis(self, processed_sample_names) -> bool:
        """运行联合分析 🤝 | Run joint analysis"""
        if not self.config.enable_joint_calling:
            self.logger.info("ℹ️ Joint calling未启用，跳过 | Joint calling not enabled, skipping")
            return True
        
        if len(processed_sample_names) < 2:
            self.logger.warning("⚠️ 需要至少2个成功处理的样品才能进行Joint calling | Need at least 2 successfully processed samples for joint calling")
            return False
        
        self.logger.info("=" * 60)
        self.logger.info("🤝 开始Joint Calling阶段 | Starting Joint Calling phase")
        
        try:
            # 获取所有已处理的样品（包括之前可能已经处理过的）
            all_processed_samples = self.joint_processor.get_processed_samples()
            
            if not all_processed_samples:
                self.logger.error("❌ 没有找到任何已处理的VCF文件 | No processed VCF files found")
                return False
            
            # 生成样品映射文件
            if not self.joint_processor.generate_sample_map(all_processed_samples):
                self.logger.error("❌ 生成样品映射文件失败 | Failed to generate sample mapping file")
                return False
            
            # 执行joint calling
            success = self.joint_processor.run_joint_calling()
            
            if success:
                self.logger.info("✅ Joint Calling阶段完成 | Joint Calling phase completed")
                return True
            else:
                self.logger.error("❌ Joint Calling阶段失败 | Joint Calling phase failed")
                return False
                
        except Exception as e:
            self.logger.error(f"❌ Joint Calling过程中出现异常 | Exception during Joint Calling: {e}")
            return False
    
    def run_analysis(self):
        """运行完整的GTX WGS分析流程 🚀 | Run complete GTX WGS analysis pipeline"""
        try:
            analysis_mode = "Joint Calling模式" if self.config.enable_joint_calling else "单样品分析模式"
            
            self.logger.info(f"🚀 开始GTX WGS批处理分析流程 ({analysis_mode}) | Starting GTX WGS batch processing analysis pipeline ({analysis_mode})")
            self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
            self.logger.info(f"📤 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
            self.logger.info(f"📁 R1文件模式 | R1 file pattern: {self.config.read1_pattern}")
            self.logger.info(f"📁 R2文件模式 | R2 file pattern: {self.config.read2_pattern}")
            
            if self.config.enable_joint_calling:
                self.logger.info(f"🤝 Joint calling已启用 | Joint calling enabled")
                self.logger.info(f"  📜 Joint输出文件 | Joint output file: {self.config.joint_output_name}")
                self.logger.info(f"  🧵 Joint calling线程数 | Joint calling threads: {self.config.joint_threads}")
            
            # 检查依赖 🧩 | Check dependencies
            self.check_dependencies()

            # 🆕 检查并构建参考基因组索引
            if not self.check_and_build_reference_index():
                raise RuntimeError("❌ 参考基因组索引检查/构建失败")
            
            # 阶段1: 单样品分析
            successful_samples, failed_samples, processed_sample_names = self.run_individual_analysis()
            
            # 阶段2: Joint calling (如果启用)
            joint_success = True
            if self.config.enable_joint_calling:
                joint_success = self.run_joint_analysis(processed_sample_names)
            
            # 生成总结报告 📄 | Generate summary report
            self.summary_generator.generate_summary_report()
            
            # 显示最终目录结构 🌳 | Display final directory structure
            self._display_final_structure(successful_samples)
            
            # 最终状态汇总
            self.logger.info("=" * 60)
            if joint_success:
                self.logger.info("🎉 GTX WGS批处理分析全部完成 | GTX WGS batch processing analysis completed")
            else:
                self.logger.warning("⚠️ GTX WGS批处理分析部分完成 (Joint calling失败) | GTX WGS batch processing analysis partially completed (Joint calling failed)")
            
            self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)
    
    def _display_final_structure(self, successful_samples):
        """显示最终目录结构 🌳 | Display final directory structure"""
        self.logger.info("")
        self.logger.info("🌳 最终输出目录结构 | Final output directory structure:")
        self.logger.info(f"  {self.config.output_dir}/")
        
        bam_count = len(list(self.config.bam_output_dir.glob('*.sorted.bam')))
        vcf_count = len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))
        
        self.logger.info(f"  ├── bam/     (💾 BAM文件: {bam_count} 个)")
        self.logger.info(f"  ├── vcf/     (📜 VCF文件: {vcf_count} 个)")
        
        if self.config.enable_joint_calling:
            joint_files = list(self.config.joint_output_dir.glob('*'))
            self.logger.info(f"  ├── joint/   (🤝 Joint calling结果: {len(joint_files)} 个文件)")
            
            # 显示joint目录内容
            if joint_files:
                for joint_file in joint_files:
                    if joint_file.name.endswith('.vcf.gz'):
                        self.logger.info(f"  │   ├── {joint_file.name} (📜 Joint VCF)")
                    elif joint_file.name == 'sample_map.txt':
                        self.logger.info(f"  │   ├── {joint_file.name} (📋 样品映射文件)")
                    else:
                        self.logger.info(f"  │   ├── {joint_file.name}")
        
        self.logger.info(f"  └── tmp/     (🗑️ 临时文件目录)")

def main():
    """主函数 🚀 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 GTX WGS批处理分析脚本 (模块化版本) | GTX WGS Batch Processing Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
💡 示例 | Examples:
  # 基本单样品分析
  %(prog)s -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
  
  # 启用Joint Calling
  %(prog)s -i ./clean_data -o ./gtx_results -r ./reference.fa --enable-joint
  
  # 自定义Joint输出文件名和线程数
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --enable-joint \\
           --joint-output merged_population.vcf.gz --joint-threads 64
  
  # 完整参数示例
  %(prog)s -i ./data -o ./results -r ./ref.fa -t 88 --enable-joint \\
           --min-confidence 25 --min-base-quality 15 --joint-threads 32
        """
    )
    
    # 必需参数 📌 | Required arguments
    parser.add_argument('-i', '--input-dir', required=True, 
                       help='输入目录路径 (包含clean FASTQ文件) 📂 | Input directory path (containing clean FASTQ files)')
    parser.add_argument('-o', '--output-dir', required=True, 
                       help='输出目录路径 📤 | Output directory path')
    parser.add_argument('-r', '--reference', required=True, 
                       help='参考基因组文件路径 🧬 | Reference genome file path')
    parser.add_argument('--force-rebuild-index', action='store_true',
                       help='强制重新构建参考基因组索引 🔨 | Force rebuild reference genome index')
    
    # 可选参数 ✨ | Optional arguments
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='线程数 🧵 | Number of threads')
    parser.add_argument('--gtx-path', default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx', 
                       help='GTX程序路径 💻 | GTX program path')
    parser.add_argument('--tmp-dir', 
                       help='临时目录路径 (默认使用输出目录下的tmp) 🗑️ | Temporary directory path (default: tmp under output directory)')
    
    # Joint calling参数 🤝 | Joint calling parameters
    parser.add_argument('--enable-joint', action='store_true',
                       help='启用Joint Calling 🤝 | Enable Joint Calling')
    parser.add_argument('--joint-output', default='merged_gtx.vcf.gz',
                       help='Joint calling输出文件名 📜 | Joint calling output filename')
    parser.add_argument('--joint-threads', type=int, default=88,
                       help='Joint calling使用的线程数 🧵 | Number of threads for joint calling')
    
    # 质量控制参数 🔬 | Quality control parameters
    parser.add_argument('--min-confidence', type=int, default=30, 
                       help='最小置信度阈值 🎯 | Minimum confidence threshold')
    parser.add_argument('--min-base-quality', type=int, default=20, 
                       help='最小碱基质量阈值 ✨ | Minimum base quality threshold')
    parser.add_argument('--ploidy', type=int, default=2, 
                       help='倍性 🧬 | Ploidy')
    parser.add_argument('--pcr-indel-model', default='CONSERVATIVE', 
                       help='PCR indel模型 🔬 | PCR indel model')
    
    # 文件模式参数 📁 | File pattern parameters
    parser.add_argument('--read1-pattern', default='*_1.fq.gz', 
                       help='R1文件匹配模式 📄 | R1 file pattern')
    parser.add_argument('--read2-pattern', default='*_2.fq.gz', 
                       help='R2文件匹配模式 📄 | R2 file pattern')
    
    args = parser.parse_args()

    # 准备额外的配置参数
    extra_kwargs = {}
    if hasattr(args, 'force_rebuild_index'):
        extra_kwargs['force_rebuild_index'] = args.force_rebuild_index
    
    # 创建分析器并运行 🚀 | Create analyzer and run
    analyzer = GTXAnalyzer(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        reference=args.reference,
        gtx_path=args.gtx_path,
        threads=args.threads,
        tmp_dir=args.tmp_dir,
        min_confidence_threshold=args.min_confidence,
        min_base_quality=args.min_base_quality,
        ploidy=args.ploidy,
        pcr_indel_model=args.pcr_indel_model,
        read1_pattern=args.read1_pattern,
        read2_pattern=args.read2_pattern,
        enable_joint_calling=args.enable_joint,
        joint_output_name=args.joint_output,
        joint_threads=args.joint_threads,
         **extra_kwargs  # 🆕 传递额外参数
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()