"""
VCF过滤主程序模块 | VCF Filtering Main Module
"""

import argparse
import sys
from .config import FilterConfig
from .utils import FilterLogger, CommandRunner, check_dependencies
from .separator import VCFSeparator
from .filter import VCFFilter
from .statistics import VCFStatistics

class VCFFilterAnalyzer:
    """VCF过滤分析主类 | Main VCF Filter Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = FilterConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = FilterLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.separator = VCFSeparator(self.config, self.logger, self.cmd_runner)
        self.filter = VCFFilter(self.config, self.logger, self.cmd_runner)
        self.statistics = VCFStatistics(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_filtering(self):
        """运行完整的过滤流程 | Run complete filtering pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("🧬 开始VCF SNP/INDEL过滤流程 | Starting VCF SNP/INDEL Filtering Pipeline")
            self.logger.info("=" * 60)
            
            # 步骤0: 检查依赖 | Step 0: Check dependencies
            self.logger.info("📦 检查依赖软件 | Checking dependencies")
            if not self.check_dependencies():
                raise RuntimeError("❌ 依赖检查失败 | Dependency check failed")
            
            # 步骤1: 分离SNP和INDEL | Step 1: Separate SNPs and INDELs
            if not self.separator.separate_variants():
                raise RuntimeError("❌ SNP/INDEL分离失败 | SNP/INDEL separation failed")
            
            # 步骤2a: 过滤SNP | Step 2a: Filter SNPs
            if not self.filter.filter_snps(self.separator.raw_snp_file):
                raise RuntimeError("❌ SNP过滤失败 | SNP filtering failed")
            
            # 步骤2b: 过滤INDEL | Step 2b: Filter INDELs
            if not self.filter.filter_indels(self.separator.raw_indel_file):
                raise RuntimeError("❌ INDEL过滤失败 | INDEL filtering failed")
            
            # 步骤3: 合并过滤后的变异 | Step 3: Merge filtered variants
            if not self.filter.merge_filtered_variants():
                raise RuntimeError("❌ 变异合并失败 | Variant merging failed")
            
            # 步骤4: 生成统计报告 | Step 4: Generate statistics report
            if not self.statistics.generate_statistics_report(self.separator, self.filter):
                raise RuntimeError("❌ 统计报告生成失败 | Statistics report generation failed")
            
            # 完成 | Complete
            self.logger.info("=" * 60)
            self.logger.info("🎉 VCF过滤流程成功完成 | VCF Filtering Pipeline Completed Successfully!")
            self.logger.info("=" * 60)
            self.logger.info(f"📂 结果保存在 | Results saved in: {self.config.output_dir}")
            self.logger.info("=" * 60 + "\n")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 VCF SNP/INDEL过滤脚本 (模块化版本) | VCF SNP/INDEL Filtering Script (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  %(prog)s -i variants.vcf -o filtered_output
  
  # 自定义SNP过滤参数 | Custom SNP filtering
  %(prog)s -i variants.vcf -o output --snp-qual 40 --snp-dp 15 --snp-mq 50
  
  # 自定义INDEL过滤参数 | Custom INDEL filtering
  %(prog)s -i variants.vcf -o output --indel-qual 35 --indel-fs 250
  
  # 指定线程数 | Specify threads
  %(prog)s -i variants.vcf -o output -t 64

输出文件 | Output Files:
  1. filtered.raw.snp.vcf.gz          - 原始SNP文件 | Raw SNP file
  2. filtered.filtered.snp.vcf.gz     - 过滤后SNP文件 | Filtered SNP file
  3. filtered.raw.indel.vcf.gz        - 原始INDEL文件 | Raw INDEL file
  4. filtered.filtered.indel.vcf.gz   - 过滤后INDEL文件 | Filtered INDEL file
  5. filtered.filtered.merged.vcf.gz  - 合并后文件 | Merged file
  6. filtering_statistics.txt         - 统计报告 | Statistics report
  7. vcf_filtering.log                - 运行日志 | Run log
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, dest='vcf_file',
                       help='输入VCF文件路径 (支持压缩和未压缩) | Input VCF file path (supports compressed and uncompressed)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./filtered_vcf', dest='output_dir',
                       help='输出目录 | Output directory (default: ./filtered_vcf)')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='线程数 | Number of threads (default: 88)')
    
    # SNP过滤参数 | SNP filtering parameters
    snp_group = parser.add_argument_group('🔹 SNP过滤参数 | SNP Filtering Parameters')
    snp_group.add_argument('--snp-qual', type=float, default=30.0,
                          help='SNP最小质量值 | SNP minimum QUAL (default: 30.0)')
    snp_group.add_argument('--snp-dp', type=int, default=10,
                          help='SNP最小测序深度 | SNP minimum DP (default: 10)')
    snp_group.add_argument('--snp-mq', type=float, default=40.0,
                          help='SNP最小比对质量 | SNP minimum MQ (default: 40.0)')
    snp_group.add_argument('--snp-qd', type=float, default=2.0,
                          help='SNP最小质量/深度比 | SNP minimum QD (default: 2.0)')
    snp_group.add_argument('--snp-fs', type=float, default=60.0,
                          help='SNP最大FisherStrand值 | SNP maximum FS (default: 60.0)')
    snp_group.add_argument('--snp-sor', type=float, default=3.0,
                          help='SNP最大StrandOddsRatio | SNP maximum SOR (default: 3.0)')
    snp_group.add_argument('--snp-mqrs', type=float, default=-12.5,
                          help='SNP最小MappingQualityRankSum | SNP minimum MQRankSum (default: -12.5)')
    snp_group.add_argument('--snp-rprs', type=float, default=-8.0,
                          help='SNP最小ReadPosRankSum | SNP minimum ReadPosRankSum (default: -8.0)')
    # [新增] MAF参数
    snp_group.add_argument('--snp-maf', type=float, default=0.05,
                          help='SNP最小次等位基因频率 | SNP minimum MAF (default: 0.05)')
    
    # INDEL过滤参数 | INDEL filtering parameters
    indel_group = parser.add_argument_group('🔸 INDEL过滤参数 | INDEL Filtering Parameters')
    indel_group.add_argument('--indel-qual', type=float, default=30.0,
                            help='INDEL最小质量值 | INDEL minimum QUAL (default: 30.0)')
    indel_group.add_argument('--indel-dp', type=int, default=10,
                            help='INDEL最小测序深度 | INDEL minimum DP (default: 10)')
    indel_group.add_argument('--indel-mq', type=float, default=40.0,
                            help='INDEL最小比对质量 | INDEL minimum MQ (default: 40.0)')
    indel_group.add_argument('--indel-qd', type=float, default=2.0,
                            help='INDEL最小质量/深度比 | INDEL minimum QD (default: 2.0)')
    indel_group.add_argument('--indel-fs', type=float, default=200.0,
                            help='INDEL最大FisherStrand值 | INDEL maximum FS (default: 200.0)')
    indel_group.add_argument('--indel-sor', type=float, default=10.0,
                            help='INDEL最大StrandOddsRatio | INDEL maximum SOR (default: 10.0)')
    indel_group.add_argument('--indel-rprs', type=float, default=-20.0,
                            help='INDEL最小ReadPosRankSum | INDEL minimum ReadPosRankSum (default: -20.0)')
    
    # 工具路径 | Tool paths
    parser.add_argument('--bcftools-path', default='bcftools',
                       help='BCFtools软件路径 | BCFtools software path (default: bcftools)')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = VCFFilterAnalyzer(
        vcf_file=args.vcf_file,
        output_dir=args.output_dir,
        threads=args.threads,
        bcftools_path=args.bcftools_path,
        # SNP parameters
        snp_qual=args.snp_qual,
        snp_dp=args.snp_dp,
        snp_mq=args.snp_mq,
        snp_qd=args.snp_qd,
        snp_fs=args.snp_fs,
        snp_sor=args.snp_sor,
        snp_mqrs=args.snp_mqrs,
        snp_rprs=args.snp_rprs,
        snp_maf=args.snp_maf,  # [新增] 传入MAF参数
        # INDEL parameters
        indel_qual=args.indel_qual,
        indel_dp=args.indel_dp,
        indel_mq=args.indel_mq,
        indel_qd=args.indel_qd,
        indel_fs=args.indel_fs,
        indel_sor=args.indel_sor,
        indel_rprs=args.indel_rprs
    )
    
    analyzer.run_filtering()

if __name__ == "__main__":
    main()
