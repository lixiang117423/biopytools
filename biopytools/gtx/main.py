"""
🚀 GTX WGS分析主程序模块 | GTX WGS Analysis Main Module 🚀
"""

import argparse
import sys
from .config import GTXConfig
from .utils import GTXLogger, CommandRunner, check_dependencies, FileProcessor
from .processor import GTXProcessor
from .results import SummaryGenerator

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
    
    def check_dependencies(self):
        """检查依赖软件 🧩 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的GTX WGS分析流程 🚀 | Run complete GTX WGS analysis pipeline"""
        try:
            self.logger.info("🚀 开始GTX WGS批处理分析流程 | Starting GTX WGS batch processing analysis pipeline")
            self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
            self.logger.info(f"📤 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
            
            # 检查依赖 🧩 | Check dependencies
            self.check_dependencies()
            
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
            
            for current, (sample_name, r1_file, r2_file) in enumerate(file_pairs, 1):
                self.logger.info(f"[{current}/{total_samples}] 🧪 处理样品 | Processing sample: {sample_name}")
                
                success = self.gtx_processor.process_sample(sample_name, r1_file, r2_file)
                
                if success:
                    successful_samples += 1
                    self.logger.info(f"  ✅ 样品 {sample_name} 分析完成 ({current}/{total_samples}) | Sample {sample_name} analysis completed ({current}/{total_samples})")
                else:
                    failed_samples.append(sample_name)
                
                self.logger.info("-" * 40)
            
            # 输出处理结果统计 📊 | Output processing results statistics
            self.logger.info("=" * 60)
            self.logger.info("🎉 批量处理完成 | Batch processing completed")
            self.logger.info(f"✅ 成功处理样品数 | Successfully processed samples: {successful_samples}")
            self.logger.info(f"❌ 失败样品数 | Failed samples: {len(failed_samples)}")
            
            if failed_samples:
                self.logger.warning(f"⚠️ 失败样品列表 | Failed samples list: {', '.join(failed_samples)}")
            
            # 生成总结报告 📄 | Generate summary report
            self.summary_generator.generate_summary_report()
            
            # 显示最终目录结构 🌳 | Display final directory structure
            self.logger.info("")
            self.logger.info("🌳 最终输出目录结构 | Final output directory structure:")
            self.logger.info(f"  {self.config.output_dir}/")
            self.logger.info(f"  ├── bam/     (💾 BAM文件: {len(list(self.config.bam_output_dir.glob('*.sorted.bam')))} 个)")
            self.logger.info(f"  ├── vcf/     (📜 VCF文件: {len(list(self.config.vcf_output_dir.glob('*.vcf.gz')))} 个)")
            self.logger.info(f"  └── tmp/     (🗑️ 临时文件目录)")
            
            self.logger.info("=" * 60)
            self.logger.info("🎉 GTX WGS批处理分析全部完成 | GTX WGS batch processing analysis completed")
            self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 🚀 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 GTX WGS批处理分析脚本 (模块化版本) | GTX WGS Batch Processing Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
💡 示例 | Examples:
  %(prog)s -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
  %(prog)s -i ./clean_data -o ./gtx_results -r ./reference.fa -t 64
  %(prog)s -i /data/clean -o /results -r /genome/ref.fa --gtx-path /custom/gtx/path
  %(prog)s -i ./data -o ./results -r ./ref.fa --min-confidence 25 --min-base-quality 15
        """
    )
    
    # 必需参数 📌 | Required arguments
    parser.add_argument('-i', '--input-dir', required=True, 
                       help='输入目录路径 (包含clean FASTQ文件) 📂 | Input directory path (containing clean FASTQ files)')
    parser.add_argument('-o', '--output-dir', required=True, 
                       help='输出目录路径 📤 | Output directory path')
    parser.add_argument('-r', '--reference', required=True, 
                       help='参考基因组文件路径 🧬 | Reference genome file path')
    
    # 可选参数 ✨ | Optional arguments
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='线程数 🧵 | Number of threads')
    parser.add_argument('--gtx-path', default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx', 
                       help='GTX程序路径 💻 | GTX program path')
    parser.add_argument('--tmp-dir', 
                       help='临时目录路径 (默认使用输出目录下的tmp) 🗑️ | Temporary directory path (default: tmp under output directory)')
    
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
    parser.add_argument('--read1-pattern', default='*_1.clean.fq.gz', 
                       help='R1文件匹配模式 📄 | R1 file pattern')
    parser.add_argument('--read2-pattern', default='*_2.clean.fq.gz', 
                       help='R2文件匹配模式 📄 | R2 file pattern')
    
    args = parser.parse_args()
    
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
        read2_pattern=args.read2_pattern
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
