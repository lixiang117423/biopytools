"""
基因组组装主程序模块 🎯 | Genome Assembly Main Module
"""

import argparse
import sys
import os
from pathlib import Path
from datetime import datetime

# 添加当前目录到sys.path，解决相对导入问题
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

try:
    from config import AssemblyConfig
    from logger import AssemblyLogger
    from utils import check_dependencies
    from assembler import HifiasmAssembler
    from report import ReportGenerator
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    sys.path.append(str(current_dir.parent))
    from biopytools.genome_assembler.config import AssemblyConfig
    from biopytools.genome_assembler.logger import AssemblyLogger
    from biopytools.genome_assembler.utils import check_dependencies
    from biopytools.genome_assembler.assembler import HifiasmAssembler
    from biopytools.genome_assembler.report import ReportGenerator

class GenomeAssembler:
    """基因组组装主类 | Main Genome Assembler Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = AssemblyConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        logger_path = Path(self.config.log_dir)
        self.logger_manager = AssemblyLogger(logger_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化组件 | Initialize components
        self.assembler = HifiasmAssembler(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)
    
    def run_assembly(self):
        """运行完整的组装流程 | Run complete assembly pipeline"""
        start_time = datetime.now()
        
        try:
            self.logger.info("=" * 80)
            self.logger.info("🧬 开始基因组组装流程 | Starting Genome Assembly Pipeline")
            self.logger.info("=" * 80)
            
            # 1. 检查依赖 | Check dependencies
            check_dependencies(self.logger)
            
            # 2. 显示配置信息 | Display configuration
            self._print_config()
            
            # 3. 执行组装 | Execute assembly
            if not self.assembler.run_assembly():
                sys.exit(1)
            
            # 4. 转换格式 | Convert formats
            fasta_results = self.assembler.convert_gfa_to_fasta()
            
            end_time = datetime.now()
            
            # 5. 生成报告 | Generate report
            self.report_generator.generate_summary(fasta_results, start_time, end_time)
            
            # 6. 输出结果 | Output results
            duration = int((end_time - start_time).total_seconds())
            self._print_summary(duration, fasta_results)
            
        except Exception as e:
            self.logger.error(f"❌ 组装流程在执行过程中意外终止 | Assembly pipeline terminated unexpectedly: {e}")
            sys.exit(1)
    
    def _print_config(self):
        """打印配置信息 | Print configuration"""
        self.logger.info("组装配置 | Assembly Configuration:")
        self.logger.info(f"  样本名称 | Sample Name: {self.config.prefix}")
        self.logger.info(f"  基因组大小 | Genome Size: {self.config.genome_size}")
        self.logger.info(f"  倍性 | Ploidy: {self.config.n_hap}")
        self.logger.info(f"  线程数 | Threads: {self.config.threads}")
        self.logger.info(f"  HiFi数据 | HiFi Data: {self.config.hifi_data}")
        self.logger.info(f"  Hi-C R1 | Hi-C R1: {self.config.hic_r1}")
        self.logger.info(f"  Hi-C R2 | Hi-C R2: {self.config.hic_r2}")
        self.logger.info(f"  工作目录 | Work Directory: {self.config.work_dir}")
    
    def _print_summary(self, duration: int, fasta_results: dict):
        """打印摘要 | Print summary"""
        self.logger.info("=" * 80)
        self.logger.info("🎉 组装流程完成 | Assembly Pipeline Completed")
        self.logger.info("=" * 80)
        self.logger.info(f"⏱️  总耗时 | Total Duration: {duration} 秒")
        self.logger.info(f"✅ 生成FASTA文件 | Generated FASTA files: {len(fasta_results)}")
        
        # 打印文件结构 | Print file structure
        self.report_generator.print_file_tree()

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 HiFi + Hi-C 基因组组装工具 (hifiasm) | HiFi + Hi-C Genome Assembly Tool (hifiasm)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  %(prog)s --hifi hifi.fq --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz -p sample1
  
  # 指定参数 | Specify parameters
  %(prog)s -i hifi.fq -1 hic_R1.fq.gz -2 hic_R2.fq.gz -p sample1 -t 88 -g 2.5g
  
  # 自定义输出目录 | Custom output directory
  %(prog)s --hifi hifi.fq --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz -p sample1 -o /data/assembly
        """
    )
    
    # 必需参数 | Required arguments
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--hifi', '-i', 
                       help='🧬 HiFi数据文件路径 | Path to HiFi data file')
    
    parser.add_argument('--hic-r1', '-1', required=True,
                       help='🧬 Hi-C R1文件路径 | Path to Hi-C R1 file')
    parser.add_argument('--hic-r2', '-2', required=True,
                       help='🧬 Hi-C R2文件路径 | Path to Hi-C R2 file')
    
    # 样本参数 | Sample parameters
    parser.add_argument('--prefix', '-p', default="genome_sample",
                       help='📝 样本前缀 | Sample prefix')
    
    # 组装参数 | Assembly parameters
    parser.add_argument('--threads', '-t', type=int, default=88,
                       help='🧵 线程数 | Number of threads')
    
    parser.add_argument('--genome-size', '-g', default="1.45g",
                       help='🧬 预估基因组大小 | Estimated genome size (e.g., 1.45g, 250m)')
    
    parser.add_argument('--n-hap', type=int, default=2,
                       help='🔢 倍性 | Ploidy (haploid count)')
    
    # 路径参数 | Path parameters
    parser.add_argument('--output', '-o', default="./assembly_output",
                       help='📁 输出目录 | Output directory')
    
    args = parser.parse_args()
    
    # 创建组装器并运行 | Create assembler and run
    assembler = GenomeAssembler(
        hifi_data=args.hifi,
        hic_r1=args.hic_r1,
        hic_r2=args.hic_r2,
        prefix=args.prefix,
        threads=args.threads,
        genome_size=args.genome_size,
        n_hap=args.n_hap,
        base_dir=args.output
    )
    
    assembler.run_assembly()

if __name__ == "__main__":
    main()
