"""
VCF联合分型主程序模块 | VCF Joint Genotyping Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import VCFJointGenotypingConfig
from .utils import setup_logging, CommandRunner, DependencyChecker
from .processor import VCFJointGenotypingProcessor
from .reporter import VCFJointGenotypingReporter

class VCFJointGenotypingAnalyzer:
    """VCF联合分型分析器 | VCF Joint Genotyping Analyzer"""
    
    def __init__(self, config: VCFJointGenotypingConfig, logger):
        self.config = config
        self.logger = logger
        
        # 初始化组件 | Initialize components
        self.cmd_runner = CommandRunner(logger, config.output_dir)
        self.dependency_checker = DependencyChecker(logger)
        self.processor = VCFJointGenotypingProcessor(config, logger, self.cmd_runner)
        self.reporter = VCFJointGenotypingReporter(config, logger)
    
    def run_analysis(self) -> bool:
        """🚀 运行完整分析流程 | Run complete analysis pipeline"""
        self.logger.info("="*60)
        self.logger.info("🚀 开始VCF联合分型和过滤处理 | Starting VCF joint genotyping and filtering")
        self.logger.info("="*60)
        
        try:
            # 步骤1: 验证配置 | Step 1: Validate configuration
            self.logger.info("📋 步骤1: 验证配置参数 | Step 1: Validating configuration parameters")
            validation_errors = self.config.validate()
            if validation_errors:
                for error in validation_errors:
                    self.logger.error(f"❌ {error}")
                return False
            
            # 步骤2: 检查依赖 | Step 2: Check dependencies
            self.logger.info("🔧 步骤2: 检查依赖软件 | Step 2: Checking dependencies")
            success, missing_deps = self.dependency_checker.check_all_dependencies(self.config)
            if not success:
                self.logger.error(f"❌ 缺少必要的依赖软件，请安装后重试 | Missing required dependencies, please install and retry")
                return False
            
            # 步骤3: 检查并创建VCF索引 | Step 3: Check and create VCF indexes
            self.logger.info("🔍 步骤3: 检查并创建VCF索引 | Step 3: Checking and creating VCF indexes")
            if not self.processor.check_and_index_vcf_files():
                return False
            
            # 步骤4: 创建样本映射文件 | Step 4: Create sample mapping file
            self.logger.info("📋 步骤4: 创建样本映射文件 | Step 4: Creating sample mapping file")
            if not self.processor.create_sample_map():
                return False
            
            # 步骤5: 联合分型 | Step 5: Joint genotyping
            self.logger.info("🔗 步骤5: GTX联合分型 | Step 5: GTX joint genotyping")
            if not self.processor.joint_genotyping():
                return False
            
            # 步骤6: 提取变异 | Step 6: Extract variants
            self.logger.info("🧬 步骤6: 提取SNP和INDEL变异 | Step 6: Extracting SNP and INDEL variants")
            if not self.processor.extract_variants():
                return False
            
            # 步骤7: 过滤变异 | Step 7: Filter variants
            self.logger.info("🔍 步骤7: 过滤变异质量 | Step 7: Filtering variant quality")
            if not self.processor.filter_variants():
                return False
            
            # 步骤8: 压缩和索引 | Step 8: Compress and index
            self.logger.info("📦 步骤8: 压缩和索引结果文件 | Step 8: Compressing and indexing result files")
            if not self.processor.compress_and_index():
                return False
            
            # 步骤9: 生成统计信息 | Step 9: Generate statistics
            self.logger.info("📊 步骤9: 生成统计信息 | Step 9: Generating statistics")
            stats = self.processor.generate_statistics()
            
            # 步骤10: 生成报告 | Step 10: Generate report
            self.logger.info("📋 步骤10: 生成处理报告 | Step 10: Generating processing report")
            if not self.reporter.generate_report(stats):
                return False
            
            # 完成 | Completion
            self.logger.info("="*60)
            self.logger.info("🎉 VCF联合分型和过滤处理完成! | VCF joint genotyping and filtering completed!")
            self.logger.info("="*60)
            
            # 显示最终结果 | Show final results
            self._show_final_results(stats)
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 分析过程中发生错误 | Error occurred during analysis: {e}")
            return False
    
    def _show_final_results(self, stats):
        """📊 显示最终结果 | Show final results"""
        self.logger.info("📊 最终结果摘要 | Final Results Summary:")
        self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
        
        output_paths = self.config.get_output_paths()
        
        if output_paths['final_snp_vcf'].exists():
            self.logger.info(f"🎯 最终SNP结果 | Final SNP results: {output_paths['final_snp_vcf'].name}")
            if 'snp_count' in stats:
                self.logger.info(f"   📊 SNP数量 | SNP count: {stats['snp_count']}")
        
        if output_paths['final_indel_vcf'].exists():
            self.logger.info(f"🎯 最终INDEL结果 | Final INDEL results: {output_paths['final_indel_vcf'].name}")
            if 'indel_count' in stats:
                self.logger.info(f"   📊 INDEL数量 | INDEL count: {stats['indel_count']}")
        
        if 'total_variants' in stats:
            self.logger.info(f"📈 总变异数 | Total variants: {stats['total_variants']}")
        
        self.logger.info(f"📋 详细报告 | Detailed report: {output_paths['report_file'].name}")

def parse_arguments():
    """解析命令行参数 | Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="VCF联合分型和过滤处理工具 | VCF Joint Genotyping and Filtering Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  python run_vcf_joint_genotyping.py --vcf-dir vcf_files --output results --reference genome.fa
  
  # 自定义过滤参数 | Custom filtering parameters
  python run_vcf_joint_genotyping.py --vcf-dir data/vcf --output analysis_results \\
      --reference ref.fa --snp-maf 0.01 --indel-maf 0.02 --snp-max-missing 0.3
  
  # 指定软件路径 | Specify software paths
  python run_vcf_joint_genotyping.py --vcf-dir vcf_input --output vcf_output \\
      --reference genome.fa --gtx-path /usr/local/bin/gtx --gatk-path gatk \\
      --bcftools-path /usr/bin/bcftools
  
  # 跳过某些步骤 | Skip certain steps
  python run_vcf_joint_genotyping.py --vcf-dir vcf_files --output results \\
      --reference genome.fa --skip-joint --skip-compress
  
  # 遇到损坏文件时终止程序 | Fail on corrupted files
  python run_vcf_joint_genotyping.py --vcf-dir vcf_files --output results \\
      --reference genome.fa --fail-on-corrupted
"""
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required arguments')
    required.add_argument(
        '--vcf-dir', '-d',
        type=str,
        required=True,
        help='VCF文件输入目录 | VCF files input directory'
    )
    required.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='输出目录 | Output directory'
    )
    required.add_argument(
        '--reference', '-r',
        type=str,
        required=True,
        help='参考基因组文件 | Reference genome file'
    )
    
    # 软件路径参数 | Software path arguments
    software = parser.add_argument_group('软件路径参数 | Software path arguments')
    software.add_argument(
        '--gtx-path',
        type=str,
        default='gtx',
        help='GTX软件路径 | GTX software path (default: gtx)'
    )
    software.add_argument(
        '--gatk-path',
        type=str,
        default='gatk',
        help='GATK软件路径 | GATK software path (default: gatk)'
    )
    software.add_argument(
        '--vcftools-path',
        type=str,
        default='vcftools',
        help='VCFtools软件路径 | VCFtools software path (default: vcftools)'
    )
    software.add_argument(
        '--bgzip-path',
        type=str,
        default='bgzip',
        help='bgzip软件路径 | bgzip software path (default: bgzip)'
    )
    software.add_argument(
        '--tabix-path',
        type=str,
        default='tabix',
        help='tabix软件路径 | tabix software path (default: tabix)'
    )
    software.add_argument(
        '--bcftools-path',
        type=str,
        default='bcftools',
        help='BCFtools软件路径 | BCFtools software path (default: bcftools)'
    )
    
    # SNP过滤参数 | SNP filtering arguments
    snp_filter = parser.add_argument_group('SNP过滤参数 | SNP filtering arguments')
    snp_filter.add_argument(
        '--snp-maf',
        type=float,
        default=0.05,
        help='SNP最小等位基因频率 | SNP minimum allele frequency (default: 0.05)'
    )
    snp_filter.add_argument(
        '--snp-max-missing',
        type=float,
        default=0.5,
        help='SNP最大缺失率 | SNP maximum missing rate (default: 0.5)'
    )
    snp_filter.add_argument(
        '--snp-hwe',
        type=float,
        default=1e-6,
        help='SNP Hardy-Weinberg平衡p值阈值 | SNP Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)'
    )
    snp_filter.add_argument(
        '--snp-min-dp',
        type=int,
        default=5,
        help='SNP最小平均深度 | SNP minimum mean depth (default: 5)'
    )
    snp_filter.add_argument(
        '--snp-max-dp',
        type=int,
        default=50,
        help='SNP最大平均深度 | SNP maximum mean depth (default: 50)'
    )
    
    # INDEL过滤参数 | INDEL filtering arguments
    indel_filter = parser.add_argument_group('INDEL过滤参数 | INDEL filtering arguments')
    indel_filter.add_argument(
        '--indel-maf',
        type=float,
        default=0.05,
        help='INDEL最小等位基因频率 | INDEL minimum allele frequency (default: 0.05)'
    )
    indel_filter.add_argument(
        '--indel-max-missing',
        type=float,
        default=1.0,
        help='INDEL最大缺失率 | INDEL maximum missing rate (default: 1.0)'
    )
    indel_filter.add_argument(
        '--indel-hwe',
        type=float,
        default=1e-6,
        help='INDEL Hardy-Weinberg平衡p值阈值 | INDEL Hardy-Weinberg equilibrium p-value threshold (default: 1e-6)'
    )
    indel_filter.add_argument(
        '--indel-min-dp',
        type=int,
        default=5,
        help='INDEL最小平均深度 | INDEL minimum mean depth (default: 5)'
    )
    indel_filter.add_argument(
        '--indel-max-dp',
        type=int,
        default=50,
        help='INDEL最大平均深度 | INDEL maximum mean depth (default: 50)'
    )
    
    # 控制参数 | Control arguments
    control = parser.add_argument_group('控制参数 | Control arguments')
    control.add_argument(
        '--skip-joint',
        action='store_true',
        help='跳过联合分型步骤 | Skip joint genotyping step'
    )
    control.add_argument(
        '--skip-extract',
        action='store_true',
        help='跳过变异提取步骤 | Skip variant extraction step'
    )
    control.add_argument(
        '--skip-filter',
        action='store_true',
        help='跳过变异过滤步骤 | Skip variant filtering step'
    )
    control.add_argument(
        '--skip-compress',
        action='store_true',
        help='跳过压缩和索引步骤 | Skip compression and indexing step'
    )
    control.add_argument(
        '--fail-on-corrupted',
        action='store_true',
        help='遇到损坏文件时终止程序（默认跳过损坏文件）| Fail program on corrupted files (default: skip corrupted files)'
    )
    
    # 输出文件前缀参数 | Output prefix arguments
    prefix = parser.add_argument_group('输出文件前缀参数 | Output prefix arguments')
    prefix.add_argument(
        '--merged-prefix',
        type=str,
        default='all.merged',
        help='合并VCF文件前缀 | Merged VCF file prefix (default: all.merged)'
    )
    prefix.add_argument(
        '--snp-prefix',
        type=str,
        default='all.merged.snp',
        help='SNP VCF文件前缀 | SNP VCF file prefix (default: all.merged.snp)'
    )
    prefix.add_argument(
        '--indel-prefix',
        type=str,
        default='all.merged.indel',
        help='INDEL VCF文件前缀 | INDEL VCF file prefix (default: all.merged.indel)'
    )
    prefix.add_argument(
        '--filtered-snp-prefix',
        type=str,
        default='final.filtered.snp',
        help='过滤后SNP文件前缀 | Filtered SNP file prefix (default: final.filtered.snp)'
    )
    prefix.add_argument(
        '--filtered-indel-prefix',
        type=str,
        default='final.filtered.indel',
        help='过滤后INDEL文件前缀 | Filtered INDEL file prefix (default: final.filtered.indel)'
    )
    
    # 系统参数 | System arguments
    system = parser.add_argument_group('系统参数 | System arguments')
    system.add_argument(
        '--memory',
        type=str,
        default='128g',
        help='GATK使用的内存量 | Memory for GATK (default: 128g)'
    )
    system.add_argument(
        '--threads',
        type=int,
        default=1,
        help='线程数 | Number of threads (default: 1)'
    )
    
    # 其他参数 | Other arguments
    other = parser.add_argument_group('其他参数 | Other arguments')
    other.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='详细输出 | Verbose output'
    )
    other.add_argument(
        '--log-file',
        type=str,
        help='日志文件路径（可选，同时输出到控制台和文件）| Log file path (optional, output to both console and file)'
    )
    
    return parser.parse_args()

def main():
    """主函数 | Main function"""
    try:
        # 解析参数 | Parse arguments
        args = parse_arguments()
        
        # 设置日志 | Setup logging
        logger = setup_logging(args.verbose, args.log_file)
        
        # 创建配置对象 | Create configuration object
        config = VCFJointGenotypingConfig(
            # 基本参数 | Basic parameters
            vcf_input_dir=args.vcf_dir,
            output_dir=args.output,
            reference_genome=args.reference,
            
            # 软件路径 | Software paths
            gtx_path=args.gtx_path,
            gatk_path=args.gatk_path,
            vcftools_path=args.vcftools_path,
            bgzip_path=args.bgzip_path,
            tabix_path=args.tabix_path,
            bcftools_path=args.bcftools_path,
            
            # SNP过滤参数 | SNP filtering parameters
            snp_maf=args.snp_maf,
            snp_max_missing=args.snp_max_missing,
            snp_hwe_pvalue=args.snp_hwe,
            snp_min_mean_dp=args.snp_min_dp,
            snp_max_mean_dp=args.snp_max_dp,
            
            # INDEL过滤参数 | INDEL filtering parameters
            indel_maf=args.indel_maf,
            indel_max_missing=args.indel_max_missing,
            indel_hwe_pvalue=args.indel_hwe,
            indel_min_mean_dp=args.indel_min_dp,
            indel_max_mean_dp=args.indel_max_dp,
            
            # 控制参数 | Control parameters
            skip_joint=args.skip_joint,
            skip_extract=args.skip_extract,
            skip_filter=args.skip_filter,
            skip_compress=args.skip_compress,
            skip_corrupted_files=not args.fail_on_corrupted,  # 反转逻辑
            
            # 输出文件前缀 | Output prefixes
            merged_prefix=args.merged_prefix,
            snp_prefix=args.snp_prefix,
            indel_prefix=args.indel_prefix,
            filtered_snp_prefix=args.filtered_snp_prefix,
            filtered_indel_prefix=args.filtered_indel_prefix,
            
            # 系统参数 | System parameters
            memory=args.memory,
            threads=args.threads
        )
        
        # 创建分析器并运行 | Create analyzer and run
        analyzer = VCFJointGenotypingAnalyzer(config, logger)
        success = analyzer.run_analysis()
        
        # 退出 | Exit
        sys.exit(0 if success else 1)
        
    except KeyboardInterrupt:
        print("\n用户中断程序 | User interrupted the program")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错 | Program execution error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
