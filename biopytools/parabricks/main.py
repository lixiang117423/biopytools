# """
# 🚀 parabricks WGS分析主程序模块 | parabricks WGS Analysis Main Module 🚀
# """

# import argparse
# import sys
# from pathlib import Path
# from .config import parabricksConfig
# from .utils import parabricksLogger, CommandRunner, check_dependencies, FileProcessor
# from .processor import parabricksProcessor
# from .results import SummaryGenerator

# class parabricksAnalyzer:
#     """🧬 parabricks WGS分析主类 | Main parabricks WGS Analyzer Class"""
    
#     # def __init__(self, **kwargs):
#     #     # 初始化配置 ⚙️ | Initialize configuration
#     #     self.config = parabricksConfig(**kwargs)
#     #     self.config.validate()
        
#     #     # 初始化日志 📝 | Initialize logging
#     #     self.logger_manager = parabricksLogger(self.config.output_path)
#     #     self.logger = self.logger_manager.get_logger()
        
#     #     # 初始化命令执行器 ▶️ | Initialize command runner
#     #     self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#     #     # 初始化各个处理器 🛠️ | Initialize processors
#     #     self.file_processor = FileProcessor(self.config, self.logger)
#     #     self.parabricks_processor = parabricksProcessor(self.config, self.logger, self.cmd_runner)
#     #     self.summary_generator = SummaryGenerator(self.config, self.logger)
#     #     # 设置容器环境
#     #     self.cmd_runner.setup_container(Path(self.config.parabricks_path))

#     def __init__(self, **kwargs):
#         # 初始化配置 ⚙️ | Initialize configuration
#         self.config = parabricksConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志 📝 | Initialize logging
#         self.logger_manager = parabricksLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 ▶️ | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#         # 🐳 设置容器环境 | Setup container environment
#         self.cmd_runner.setup_container(Path(self.config.parabricks_path))
        
#         # 初始化各个处理器 🛠️ | Initialize processors
#         self.file_processor = FileProcessor(self.config, self.logger)
#         self.parabricks_processor = parabricksProcessor(self.config, self.logger, self.cmd_runner)
#         self.summary_generator = SummaryGenerator(self.config, self.logger)
    
#     def check_dependencies(self):
#         """检查依赖软件 🧩 | Check dependencies"""
#         return check_dependencies(self.config, self.logger)
    
#     # def run_analysis(self):
#     #     """运行完整的parabricks WGS分析流程 🚀 | Run complete parabricks WGS analysis pipeline"""
#     #     try:
#     #         self.logger.info("🚀 开始parabricks WGS批处理分析流程 | Starting parabricks WGS batch processing analysis pipeline")
#     #         self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
#     #         self.logger.info(f"📤 输出目录 | Output directory: {self.config.output_dir}")
#     #         self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
            
#     #         # 检查依赖 🧩 | Check dependencies
#     #         self.check_dependencies()
            
#     #         # 查找输入文件 🔍 | Find input files
#     #         file_pairs = self.file_processor.find_fastq_files()
            
#     #         if not file_pairs:
#     #             raise RuntimeError("❌ 没有找到需要处理的FASTQ文件对 | No FASTQ file pairs found to process")
            
#     #         # 处理每个样品 🔬 | Process each sample
#     #         self.logger.info("=" * 60)
#     #         self.logger.info("🔬 开始批量处理样品 | Starting batch processing of samples")
            
#     #         total_samples = len(file_pairs)
#     #         successful_samples = 0
#     #         failed_samples = []
            
#     #         for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
#     #             self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing sample: {sample_name}")
                
#     #             success = self.parabricks_processor.process_sample(sample_name, r1_file, r2_file)
                
#     #             if success:
#     #                 successful_samples += 1
#     #                 self.logger.info(f"  ✅ 样品 {sample_name} 分析完成 ({current}/{total_samples}) | Sample {sample_name} analysis completed ({current}/{total_samples})")
#     #             else:
#     #                 failed_samples.append(sample_name)
                
#     #             self.logger.info("-" * 40)
            
#     #         # 输出处理结果统计 📊 | Output processing results statistics
#     #         self.logger.info("=" * 60)
#     #         self.logger.info("🎉 批量处理完成 | Batch processing completed")
#     #         self.logger.info(f"✅ 成功处理样品数 | Successfully processed samples: {successful_samples}")
#     #         self.logger.info(f"❌ 失败样品数 | Failed samples: {len(failed_samples)}")
            
#     #         if failed_samples:
#     #             self.logger.warning(f"⚠️ 失败样品列表 | Failed samples list: {', '.join(failed_samples)}")
            
#     #         # 生成总结报告 📄 | Generate summary report
#     #         self.summary_generator.generate_summary_report()
            
#     #         # 显示最终目录结构 🌳 | Display final directory structure
#     #         self.logger.info("")
#     #         self.logger.info("🌳 最终输出目录结构 | Final output directory structure:")
#     #         self.logger.info(f"  {self.config.output_dir}/")
#     #         self.logger.info(f"  ├── bam/     (💾 BAM文件: {len(list(self.config.bam_output_dir.glob('*.sorted.bam')))} 个)")
#     #         self.logger.info(f"  ├── vcf/     (📜 VCF文件: {len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))} 个)")
#     #         self.logger.info(f"  └── tmp/     (🗑️ 临时文件目录)")
            
#     #         self.logger.info("=" * 60)
#     #         self.logger.info("🎉 parabricks WGS批处理分析全部完成 | parabricks WGS batch processing analysis completed")
#     #         self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
#     #     except Exception as e:
#     #         self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
#     #         sys.exit(1)

#     def run_analysis(self):
#         """运行完整的parabricks WGS分析流程 🚀 | Run complete parabricks WGS analysis pipeline"""
#         try:
#             self.logger.info("🚀 开始parabricks WGS批处理分析流程 | Starting parabricks WGS batch processing analysis pipeline")
#             self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
#             self.logger.info(f"📤 输出目录 | Output directory: {self.config.output_dir}")
#             self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
            
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
#             processed_sample_names = []  # 新增：记录成功处理的样本名
            
#             for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
#                 self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing sample: {sample_name}")
                
#                 success = self.parabricks_processor.process_sample(sample_name, r1_file, r2_file)
                
#                 if success:
#                     successful_samples += 1
#                     processed_sample_names.append(sample_name)  # 新增：记录成功的样本
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
            
#             # 新增：执行Joint Calling（如果启用）
#             if self.config.joint_calling and processed_sample_names:
#                 self.parabricks_processor.joint_calling(processed_sample_names)
            
#             # 生成总结报告 📄 | Generate summary report
#             self.summary_generator.generate_summary_report()
            
#             # 显示最终目录结构 🌳 | Display final directory structure
#             self.logger.info("")
#             self.logger.info("🌳 最终输出目录结构 | Final output directory structure:")
#             self.logger.info(f"  {self.config.output_dir}/")
#             self.logger.info(f"  ├── bam/     (💾 BAM文件: {len(list(self.config.bam_output_dir.glob('*.sorted.bam')))} 个)")
            
#             # 根据是否有joint calling显示不同的VCF统计
#             if self.config.gvcf:
#                 gvcf_count = len(list(self.config.vcf_output_dir.glob('*.g.vcf.gz')))
#                 if self.config.joint_calling:
#                     individual_gvcf = gvcf_count - (1 if (self.config.vcf_output_dir / self.config.combined_output_name).exists() else 0)
#                     self.logger.info(f"  ├── vcf/     (📜 个体GVCF: {individual_gvcf} 个, 🔗 Combined GVCF: {1 if (self.config.vcf_output_dir / self.config.combined_output_name).exists() else 0} 个)")
#                 else:
#                     self.logger.info(f"  ├── vcf/     (📜 GVCF文件: {gvcf_count} 个)")
#             else:
#                 vcf_count = len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))
#                 self.logger.info(f"  ├── vcf/     (📜 VCF文件: {vcf_count} 个)")
            
#             self.logger.info(f"  └── tmp/     (🗑️ 临时文件目录)")
            
#             self.logger.info("=" * 60)
#             self.logger.info("🎉 parabricks WGS批处理分析全部完成 | parabricks WGS batch processing analysis completed")
#             self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
#         except Exception as e:
#             self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
#             sys.exit(1)

# def main():
#     """主函数 🚀 | Main function"""
#     parser = argparse.ArgumentParser(
#         description='🧬 parabricks WGS批处理分析脚本 (模块化版本) | parabricks WGS Batch Processing Analysis Script (Modular Version)',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         epilog="""
# 💡 示例 | Examples:
#   %(prog)s -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
#   %(prog)s -i ./clean_data -o ./parabricks_results -r ./reference.fa -t 64
#   %(prog)s -i /data/clean -o /results -r /genome/ref.fa --parabricks-path /custom/parabricks/path
#   %(prog)s -i ./data -o ./results -r ./ref.fa --min-confidence 25 --min-base-quality 15
#         """
#     )
    
#     # 必需参数 📌 | Required arguments
#     parser.add_argument('-i', '--input-dir', required=True, 
#                        help='📂 输入目录路径 (包含clean FASTQ文件) | Input directory path (containing clean FASTQ files)')
#     parser.add_argument('-o', '--output-dir', required=True, 
#                        help='📤 输出目录路径 | Output directory path')
#     parser.add_argument('-r', '--reference', required=True, 
#                        help='🧬 参考基因组文件路径 | Reference genome file path')
#     ## 新增：GVCF输出控制参数 🧬 | GVCF output control
#     parser.add_argument('--no-gvcf', action='store_false', dest='gvcf',
#                     help='📄 输出VCF而非GVCF | Output VCF instead of GVCF (default: output GVCF)')

#     # 修改：Joint Calling参数（默认开启）🔗 | Joint Calling parameters (enabled by default)
#     parser.add_argument('--no-joint-calling', action='store_false', dest='joint_calling', default=True,
#                     help='🔗 禁用Joint Calling | Disable Joint Calling (default: enabled)')
#     parser.add_argument('--combined-output', default='combined.g.vcf',
#                     help='📜 Joint Calling输出文件名 | Joint Calling output filename (default: combined.g.vcf)')
    
#     # 可选参数 ✨ | Optional arguments
#     parser.add_argument('-t', '--threads', type=int, default=88, 
#                        help='🧵 线程数 | Number of threads')
#     parser.add_argument('--parabricks-path', default='/share/apps/containers/parabricks.sif', 
#                        help='💻 parabricks程序路径 | parabricks program path')
#     parser.add_argument('--tmp-dir', 
#                        help='🗑️ 临时目录路径 (默认使用输出目录下的tmp) | Temporary directory path (default: tmp under output directory)')
    
#     # 质量控制参数 🔬 | Quality control parameters
#     parser.add_argument('--min-confidence', type=int, default=30, 
#                        help='🎯 最小置信度阈值 | Minimum confidence threshold')
#     parser.add_argument('--min-base-quality', type=int, default=20, 
#                        help='✨ 最小碱基质量阈值 | Minimum base quality threshold')
#     parser.add_argument('--ploidy', type=int, default=2, 
#                        help='🧬 倍性 | Ploidy')
#     parser.add_argument('--pcr-indel-model', default='CONSERVATIVE', 
#                        help='🔬 PCR indel模型 | PCR indel model')
    
#     # 文件模式参数 📁 | File pattern parameters
#     parser.add_argument('--read1-pattern', default='*_1.clean.fq.gz', 
#                        help='📄 R1文件匹配模式 | R1 file pattern')
#     parser.add_argument('--read2-pattern', default='*_2.clean.fq.gz', 
#                        help='📄 R2文件匹配模式 | R2 file pattern')
    
    
#     args = parser.parse_args()
    
#     # 创建分析器并运行 🚀 | Create analyzer and run
#     analyzer = parabricksAnalyzer(
#         input_dir=args.input_dir,
#         output_dir=args.output_dir,
#         reference=args.reference,
#         parabricks_path=args.parabricks_path,
#         threads=args.threads,
#         tmp_dir=args.tmp_dir,
#         min_confidence_threshold=args.min_confidence,
#         min_base_quality=args.min_base_quality,
#         ploidy=args.ploidy,
#         pcr_indel_model=args.pcr_indel_model,
#         read1_pattern=args.read1_pattern,
#         read2_pattern=args.read2_pattern,
#         gvcf=args.gvcf,  # 添加
#         joint_calling=args.joint_calling,  # 新增
#         combined_output_name=args.combined_output  # 新增
#     )
    
#     analyzer.run_analysis()

# if __name__ == "__main__":
#     main()


"""
🚀 parabricks WGS分析主程序模块 | parabricks WGS Analysis Main Module 🚀
"""

import argparse
import sys
from pathlib import Path
from .config import parabricksConfig
from .utils import parabricksLogger, CommandRunner, check_dependencies, FileProcessor
from .processor import parabricksProcessor
from .results import SummaryGenerator

class parabricksAnalyzer:
    """🧬 parabricks WGS分析主类 | Main parabricks WGS Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 ⚙️ | Initialize configuration
        self.config = parabricksConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 📝 | Initialize logging
        self.logger_manager = parabricksLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 ▶️ | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 🐳 设置容器环境 | Setup container environment
        self.cmd_runner.setup_container(Path(self.config.parabricks_path))
        
        # 初始化各个处理器 🛠️ | Initialize processors
        self.file_processor = FileProcessor(self.config, self.logger)
        self.parabricks_processor = parabricksProcessor(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 🧩 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的parabricks WGS分析流程 🚀 | Run complete parabricks WGS analysis pipeline"""
        try:
            self.logger.info("🚀 开始parabricks WGS分析 | Starting parabricks WGS analysis")
            self.logger.info(f"📂 输入目录 | Input: {self.config.input_dir}")
            self.logger.info(f"📤 输出目录 | Output: {self.config.output_dir}")
            self.logger.info(f"🧬 参考基因组 | Reference: {self.config.reference}")
            self.logger.info(f"🔀 工作流程 | Workflow: {self.config.workflow}")
            
            # 检查依赖 🧩 | Check dependencies
            self.check_dependencies()
            
            # 特殊处理：如果只运行genotypegvcf，不需要查找FASTQ文件
            if self.config.workflow == "genotypegvcf":
                self.logger.info("=" * 60)
                self.logger.info("🔗 仅执行genotypegvcf步骤 | Running genotypegvcf only")
                
                # 自动检测所有已有的GVCF文件
                gvcf_files = list(self.config.vcf_output_dir.glob("*.g.vcf.gz"))
                sample_names = [f.stem.replace(".g.vcf", "") for f in gvcf_files]
                
                if not sample_names:
                    raise RuntimeError("❌ 未找到GVCF文件 | No GVCF files found")
                
                self.logger.info(f"📊 找到 {len(sample_names)} 个GVCF文件 | Found {len(sample_names)} GVCF files")
                
                success = self.parabricks_processor.run_genotypegvcf(sample_names)
                
                if success:
                    self.logger.info("✅ genotypegvcf完成 | genotypegvcf completed")
                else:
                    self.logger.error("❌ genotypegvcf失败 | genotypegvcf failed")
                    sys.exit(1)
                
                # 生成总结报告
                self.summary_generator.generate_summary_report()
                self._display_final_structure()
                return
            
            # 查找输入文件 🔍 | Find input files
            file_pairs = self.file_processor.find_fastq_files()
            
            if not file_pairs:
                raise RuntimeError("❌ 没有找到FASTQ文件对 | No FASTQ file pairs found")
            
            # 处理每个样品 🔬 | Process each sample
            self.logger.info("=" * 60)
            self.logger.info(f"🔬 开始批量处理 ({self.config.workflow}) | Starting batch processing")
            
            total_samples = len(file_pairs)
            successful_samples = 0
            failed_samples = []
            processed_sample_names = []
            
            for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
                self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing: {sample_name}")
                
                success = self.parabricks_processor.process_sample(sample_name, r1_file, r2_file)
                
                if success:
                    successful_samples += 1
                    processed_sample_names.append(sample_name)
                    self.logger.info(f"  ✅ 完成 ({current}/{total_samples}) | Completed")
                else:
                    failed_samples.append(sample_name)
                
                self.logger.info("-" * 40)
            
            # 输出处理结果统计 📊 | Output processing statistics
            self.logger.info("=" * 60)
            self.logger.info("🎉 批量处理完成 | Batch processing completed")
            self.logger.info(f"✅ 成功: {successful_samples} | Success: {successful_samples}")
            self.logger.info(f"❌ 失败: {len(failed_samples)} | Failed: {len(failed_samples)}")
            
            if failed_samples:
                self.logger.warning(f"⚠️ 失败样品 | Failed: {', '.join(failed_samples)}")
            
            # 执行Joint Calling（如果是all或genotypegvcf workflow且启用了joint_calling）
            if self.config.workflow in ["all", "genotypegvcf"] and self.config.joint_calling:
                if processed_sample_names or self.config.workflow == "genotypegvcf":
                    # 如果是all workflow，使用处理过的样品名
                    # 如果是genotypegvcf，自动检测所有GVCF
                    if self.config.workflow == "genotypegvcf":
                        gvcf_files = list(self.config.vcf_output_dir.glob("*.g.vcf.gz"))
                        sample_names_for_joint = [f.stem.replace(".g.vcf", "") for f in gvcf_files]
                    else:
                        sample_names_for_joint = processed_sample_names
                    
                    self.parabricks_processor.run_genotypegvcf(sample_names_for_joint)
            
            # 生成总结报告 📄 | Generate summary report
            self.summary_generator.generate_summary_report()
            
            # 显示最终目录结构 🌳 | Display final structure
            self._display_final_structure()
            
            self.logger.info("=" * 60)
            self.logger.info("🎉 parabricks WGS分析完成 | Analysis completed")
            self.logger.info(f"📁 结果保存在 | Results in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"💥 分析流程异常终止 | Pipeline terminated: {e}")
            sys.exit(1)
    
    def _display_final_structure(self):
        """显示最终目录结构 🌳 | Display final directory structure"""
        self.logger.info("")
        self.logger.info("🌳 输出目录结构 | Output structure:")
        self.logger.info(f"  {self.config.output_dir}/")
        
        bam_count = len(list(self.config.bam_output_dir.glob('*.sorted.bam')))
        self.logger.info(f"  ├── bam/     (💾 BAM: {bam_count} 个)")
        
        if self.config.gvcf:
            gvcf_count = len(list(self.config.vcf_output_dir.glob('*.g.vcf.gz')))
            if self.config.joint_calling:
                individual_gvcf = gvcf_count - (1 if (self.config.vcf_output_dir / 
                    Path(self.config.combined_output_name).with_suffix('').with_suffix('.vcf.gz').name).exists() else 0)
                combined_exists = (self.config.vcf_output_dir / 
                    Path(self.config.combined_output_name).with_suffix('').with_suffix('.vcf.gz').name).exists()
                self.logger.info(f"  ├── vcf/     (📜 个体GVCF: {individual_gvcf} 个, 🔗 Combined: {'是' if combined_exists else '否'})")
            else:
                self.logger.info(f"  ├── vcf/     (📜 GVCF: {gvcf_count} 个)")
        else:
            vcf_count = len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))
            self.logger.info(f"  ├── vcf/     (📜 VCF: {vcf_count} 个)")
        
        self.logger.info(f"  └── tmp/     (🗑️ 临时文件)")

def main():
    """主函数 🚀 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 parabricks WGS批量分析工具 | parabricks WGS Batch Analysis Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
💡 示例 | Examples:
  # 运行完整流程 (fq2bam + haplotypecaller + genotypegvcf)
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow all
  
  # 仅运行 fq2bam (FASTQ转BAM)
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow fq2bam
  
  # 仅运行 haplotypecaller (BAM转VCF/GVCF)
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow haplotypecaller
  
  # 仅运行 genotypegvcf (Joint Calling)
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow genotypegvcf
  
  # 先运行 fq2bam，稍后再运行 haplotypecaller
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow fq2bam
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --workflow haplotypecaller
        """
    )
    
    # 必需参数 📌 | Required arguments
    parser.add_argument('-i', '--input-dir', required=True,
                       help='📂 输入目录 (FASTQ文件) | Input directory (FASTQ files)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='📤 输出目录 | Output directory')
    parser.add_argument('-r', '--reference', required=True,
                       help='🧬 参考基因组 | Reference genome')
    
    # 工作流程参数 🔀 | Workflow parameters
    parser.add_argument('-w', '--workflow', 
                       choices=['fq2bam', 'haplotypecaller', 'genotypegvcf', 'all'],
                       default='all',
                       help='🔀 工作流程 | Workflow: fq2bam(比对), haplotypecaller(变异检测), genotypegvcf(合并), all(完整流程)')
    
    # 可选参数 ✨ | Optional arguments
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Threads')
    parser.add_argument('--parabricks-path', default='/share/apps/containers/parabricks.sif',
                       help='💻 parabricks路径 | parabricks path')
    parser.add_argument('--tmp-dir',
                       help='🗑️ 临时目录 | Temporary directory')
    
    # 质量控制参数 🔬 | Quality control parameters
    parser.add_argument('--min-confidence', type=int, default=30,
                       help='🎯 最小置信度 | Min confidence')
    parser.add_argument('--min-base-quality', type=int, default=20,
                       help='✨ 最小碱基质量 | Min base quality')
    parser.add_argument('--ploidy', type=int, default=2,
                       help='🧬 倍性 | Ploidy')
    parser.add_argument('--pcr-indel-model', default='CONSERVATIVE',
                       help='🔬 PCR indel模型 | PCR indel model')
    
    # 文件模式参数 📁 | File pattern parameters
    parser.add_argument('--read1-pattern', default='*_1.clean.fq.gz',
                       help='📄 R1文件模式 | R1 pattern')
    parser.add_argument('--read2-pattern', default='*_2.clean.fq.gz',
                       help='📄 R2文件模式 | R2 pattern')
    
    # GVCF和Joint Calling参数 🧬 | GVCF and Joint Calling parameters
    parser.add_argument('--no-gvcf', action='store_false', dest='gvcf',
                       help='📄 输出VCF而非GVCF | Output VCF instead of GVCF')
    parser.add_argument('--no-joint-calling', action='store_false', dest='joint_calling', default=True,
                       help='🔗 禁用Joint Calling | Disable Joint Calling')
    parser.add_argument('--combined-output', default='combined.g.vcf',
                       help='📜 Joint Calling输出文件名 | Joint Calling output filename')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 🚀 | Create analyzer and run
    analyzer = parabricksAnalyzer(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        reference=args.reference,
        workflow=args.workflow,  # 新增
        parabricks_path=args.parabricks_path,
        threads=args.threads,
        tmp_dir=args.tmp_dir,
        min_confidence_threshold=args.min_confidence,
        min_base_quality=args.min_base_quality,
        ploidy=args.ploidy,
        pcr_indel_model=args.pcr_indel_model,
        read1_pattern=args.read1_pattern,
        read2_pattern=args.read2_pattern,
        gvcf=args.gvcf,
        joint_calling=args.joint_calling,
        combined_output_name=args.combined_output
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()