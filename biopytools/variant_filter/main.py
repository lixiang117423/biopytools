"""
🧬 变异筛选主程序模块 | Variant Filtering Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import VariantFilterConfig
from .utils import FilterLogger, CommandRunner, check_dependencies, get_vcf_stats
from .variant_selection import VariantSelector
from .variant_filtering import VariantFilter
from .chromosome_filter import ChromosomeFilter
from .compression import CompressionHandler

class VariantFilterPipeline:
    """🧬 变异筛选主类 | Main Variant Filter Class"""
    
    def __init__(self, **kwargs):
        # ⚙️ 初始化配置 | Initialize configuration
        self.config = VariantFilterConfig(**kwargs)
        self.config.validate()
        
        # 📝 初始化日志 | Initialize logging
        self.logger_manager = FilterLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 💻 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir, self.config.temp_dir)
        
        # 🧩 初始化各个处理器 | Initialize processors
        self.variant_selector = VariantSelector(self.config, self.logger, self.cmd_runner)
        self.variant_filter = VariantFilter(self.config, self.logger, self.cmd_runner)
        self.chromosome_filter = ChromosomeFilter(self.config, self.logger, self.cmd_runner)
        self.compression_handler = CompressionHandler(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """🔗 检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_pipeline(self):
        """🚀 运行完整的变异筛选流程 | Run complete variant filtering pipeline"""
        try:
            self.logger.info(f"🚀 开始变异筛选流程 - 变异类型: {self.config.variant_type} | Starting variant filtering pipeline - variant type: {self.config.variant_type}")
            
            # 📊 获取输入统计信息 | Get input statistics
            input_stats = get_vcf_stats(self.config.input_vcf)
            self.logger.info(f"📊 输入VCF统计 | Input VCF statistics: {input_stats['samples']} 样本, {input_stats['variants']} 变异")
            
            # 🔗 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 🎯 步骤1: 变异类型选择 | Step 1: Variant type selection
            if not self.variant_selector.select_variants():
                raise RuntimeError("❌ 变异选择失败 | Variant selection failed")
            
            # 📉 步骤2: VCFtools过滤 | Step 2: VCFtools filtering
            if not self.variant_filter.filter_variants():
                raise RuntimeError("❌ VCFtools变异过滤失败 | VCFtools variant filtering failed")
            
            # 📦 步骤3: 压缩和索引 | Step 3: Compression and indexing
            if not self.compression_handler.compress_and_index():
                raise RuntimeError("❌ 压缩和索引失败 | Compression and indexing failed")
            
            # 🧩 步骤4: 染色体过滤（可选）| Step 4: Chromosome filtering (optional)
            if not self.chromosome_filter.filter_chromosomes():
                raise RuntimeError("❌ 染色体过滤失败 | Chromosome filtering failed")
            
            # 📋 生成最终报告 | Generate final report
            self.generate_summary()
            
            self.logger.info("🎉 变异筛选流程完成 | Variant filtering pipeline completed successfully")
            
        except Exception as e:
            self.logger.error(f"💥 流程执行失败 | Pipeline execution failed: {str(e)}")
            raise
    
    def generate_summary(self):
        """📋 生成总结报告 | Generate summary report"""
        summary_file = self.config.output_dir / "filtering_summary.txt"
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("🧬 变异筛选总结报告 | Variant Filtering Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            # 📁 输入信息 | Input information
            f.write("📁 输入信息 | Input Information:\n")
            f.write(f"  - 📄 输入VCF文件 | Input VCF file: {self.config.input_vcf}\n")
            f.write(f"  - 🧬 变异类型 | Variant type: {self.config.variant_type}\n")
            
            # 🔧 变异选择方法 | Variant selection method
            f.write("\n🔧 变异选择方法 | Variant Selection Method:\n")
            if self.config.variant_type == "BOTH":
                f.write("  - 🔄 方法 | Method: 保留所有变异类型 | Keep all variant types\n")
            elif self.config.use_gatk_select:
                f.write("  - 🧬 方法 | Method: GATK SelectVariants\n")
            else:
                f.write("  - 🔧 方法 | Method: BCFtools view (推荐 | Recommended)\n")
            
            # ⚙️ 处理参数 | Processing parameters
            f.write("\n⚙️ 处理参数 | Processing Parameters:\n")
            f.write(f"  - 📊 MAF阈值 | MAF threshold: {self.config.maf}\n")
            f.write(f"  - 📊 最大缺失率 | Max missing rate: {self.config.max_missing}\n")
            f.write(f"  - 📊 HWE p值阈值 | HWE p-value threshold: {self.config.hwe}\n")
            f.write(f"  - 📊 最小平均深度 | Min mean depth: {self.config.min_meanDP}\n")
            f.write(f"  - 📊 最大平均深度 | Max mean depth: {self.config.max_meanDP}\n")
            f.write(f"  - 📊 最小质量值 | Min quality: {self.config.minQ}\n")
            f.write(f"  - 📊 最小基因型质量 | Min genotype quality: {self.config.minGQ}\n")
            
            # 🧩 染色体过滤 | Chromosome filtering
            if self.config.chromosome_filter:
                f.write(f"  - 🧩 染色体过滤模式 | Chromosome filter pattern: {self.config.chromosome_pattern}\n")
            
            # 💾 输出文件 | Output files
            f.write("\n💾 输出文件 | Output Files:\n")
            if self.config.variant_type != "BOTH":
                f.write(f"  - 🎯 变异选择结果 | Variant selection result: {self.config.selected_vcf}\n")
            f.write(f"  - ✅ 最终过滤结果 | Final filtered result: {self.config.final_vcf}\n")
            if self.config.chromosome_filter:
                f.write(f"  - 🧩 染色体过滤结果 | Chromosome filtered result: {self.config.chr_filtered_vcf}\n")
            if self.config.create_index:
                f.write(f"  - 🏷️ 索引文件 | Index files: *.tbi\n")
            
            f.write(f"\n📁 输出目录 | Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"📋 总结报告已生成 | Summary report generated: {summary_file}")

def main():
    """🚀 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🧬 变异筛选工具 (支持SNP/INDEL) | Variant Filtering Tool (Support SNP/INDEL)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📚 示例用法 | Example Usage:
    # 🧬 SNP筛选 | SNP filtering
    python -m variant_filter -i input.vcf.gz -o output_dir -t SNP
    
    # 🧬 INDEL筛选 | INDEL filtering  
    python -m variant_filter -i input.vcf.gz -o output_dir -t INDEL \\
        --max-missing 1.0 --chromosome-filter --chromosome-pattern "^OV"
    
    # 🔧 使用BCFtools进行变异选择（推荐）| Use BCFtools for variant selection (recommended)
    python -m variant_filter -i input.vcf.gz -o output_dir -t INDEL \\
        --no-gatk --maf 0.05 --max-missing 1.0
    
    # 🧬 同时处理SNP和INDEL | Process both SNP and INDEL
    python -m variant_filter -i input.vcf.gz -o output_dir -t BOTH \\
        --maf 0.01 --max-missing 0.95
    
    # 🎯 完整INDEL处理流程 | Complete INDEL processing pipeline
    python -m variant_filter -i input.vcf.gz -o output_dir -t INDEL \\
        --no-gatk --maf 0.05 --max-missing 1.0 --hwe 0.001 \\
        --min-meanDP 5 --max-meanDP 50 --chromosome-filter \\
        --chromosome-pattern "^OV" --verbose
        """
    )
    
    # 📋 必需参数 | Required parameters
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    required.add_argument('-i', '--input-vcf', required=True,
                         help='📄 输入VCF文件路径 | Input VCF file path')
    
    # ⚙️ 基本参数 | Basic parameters
    basic_group = parser.add_argument_group('⚙️ 基本参数 | Basic parameters')
    basic_group.add_argument('-o', '--output-dir', default='./variant_filter_output',
                            help='📁 输出目录 | Output directory (default: ./variant_filter_output)')
    basic_group.add_argument('-t', '--variant-type', choices=['SNP', 'INDEL', 'BOTH'], default='SNP',
                            help='🧬 变异类型 | Variant type (default: SNP)')
    basic_group.add_argument('--output-prefix', default='filtered_variants',
                            help='📝 输出文件前缀 | Output file prefix (default: filtered_variants)')
    
    # 🔧 变异选择参数 | Variant selection parameters
    selection_group = parser.add_argument_group('🔧 变异选择参数 | Variant selection parameters')
    selection_group.add_argument('--use-gatk', action='store_true', dest='use_gatk_select',
                                help='🧬 使用GATK SelectVariants进行变异选择 | Use GATK SelectVariants for variant selection')
    selection_group.add_argument('--no-gatk', action='store_false', dest='use_gatk_select',
                                help='❌ 不使用GATK，优先使用BCFtools (推荐) | Do not use GATK, prefer BCFtools (recommended)')
    selection_group.add_argument('--gatk-path', default='gatk',
                                help='🧬 GATK程序路径 | GATK program path (default: gatk)')
    selection_group.add_argument('--bcftools-path', default='bcftools',
                                help='🔧 BCFtools程序路径 | BCFtools program path (default: bcftools)')
    selection_group.add_argument('--java-options', default='-Xmx128g',
                                help='☕ Java选项 | Java options (default: -Xmx128g)')
    
    # 📉 过滤参数 | Filtering parameters
    filter_group = parser.add_argument_group('📉 过滤参数 | Filtering parameters')
    filter_group.add_argument('--vcftools-path', default='vcftools',
                             help='🔧 VCFtools程序路径 | VCFtools program path (default: vcftools)')
    filter_group.add_argument('--maf', type=float,
                             help='📊 最小等位基因频率 | Minimum allele frequency (default: 0.05 for SNP/INDEL)')
    filter_group.add_argument('--max-missing', type=float,
                             help='📊 最大缺失率 | Maximum missing rate (default: 0.9 for SNP, 1.0 for INDEL)')
    filter_group.add_argument('--hwe', type=float, default=0.001,
                             help='📊 HWE检验p值阈值 | HWE test p-value threshold (default: 0.001)')
    filter_group.add_argument('--min-meanDP', type=int, default=5,
                             help='📊 最小平均深度 | Minimum mean depth (default: 5)')
    filter_group.add_argument('--max-meanDP', type=int, default=50,
                             help='📊 最大平均深度 | Maximum mean depth (default: 50)')
    filter_group.add_argument('--minQ', type=int, default=30,
                             help='📊 最小质量值 | Minimum quality score (default: 30)')
    filter_group.add_argument('--minGQ', type=int, default=20,
                             help='📊 最小基因型质量 | Minimum genotype quality (default: 20)')
    filter_group.add_argument('--max-alleles', type=int, default=2,
                             help='📊 最大等位基因数 | Maximum number of alleles (default: 2)')
    filter_group.add_argument('--min-alleles', type=int, default=2,
                             help='📊 最小等位基因数 | Minimum number of alleles (default: 2)')
    filter_group.add_argument('--thin', type=int,
                             help='📊 稀疏化变异，保留每N个碱基的变异 | Thin variants, keep every Nth base')
    
    # 🧩 染色体过滤参数 | Chromosome filtering parameters
    chr_group = parser.add_argument_group('🧩 染色体过滤参数 | Chromosome filtering parameters')
    chr_group.add_argument('--chromosome-filter', action='store_true',
                          help='🧩 启用染色体过滤 | Enable chromosome filtering')
    chr_group.add_argument('--chromosome-pattern', default='^OV',
                          help='🧩 染色体过滤模式 | Chromosome filter pattern (default: ^OV)')
    
    # 📦 压缩参数 | Compression parameters
    compress_group = parser.add_argument_group('📦 压缩参数 | Compression parameters')
    compress_group.add_argument('--bgzip-path', default='bgzip',
                               help='📦 bgzip程序路径 | bgzip program path (default: bgzip)')
    compress_group.add_argument('--tabix-path', default='tabix',
                               help='🏷️ tabix程序路径 | tabix program path (default: tabix)')
    compress_group.add_argument('--no-compress', action='store_true',
                               help='❌ 不压缩输出文件 | Do not compress output files')
    compress_group.add_argument('--no-index', action='store_true',
                               help='❌ 不创建索引文件 | Do not create index files')
    
    # ⚙️ 其他参数 | Other parameters
    other_group = parser.add_argument_group('⚙️ 其他参数 | Other parameters')
    other_group.add_argument('-j', '--threads', type=int, default=1,
                            help='🧮 线程数 | Number of threads (default: 1)')
    other_group.add_argument('--memory', default='128g',
                            help='💾 内存限制 | Memory limit (default: 128g)')
    other_group.add_argument('--temp-dir',
                            help='📁 临时目录 | Temporary directory')
    other_group.add_argument('--keep-intermediate', action='store_true',
                            help='💾 保留中间文件 | Keep intermediate files')
    other_group.add_argument('-v', '--verbose', action='store_true',
                            help='📝 详细输出 | Verbose output')
    
    args = parser.parse_args()
    
    # ⚙️ 设置默认的use_gatk_select值 | Set default use_gatk_select value
    if not hasattr(args, 'use_gatk_select'):
        args.use_gatk_select = False
    
    # 🚀 创建分析器实例 | Create analyzer instance
    analyzer = VariantFilterPipeline(
        input_vcf=args.input_vcf,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        variant_type=args.variant_type,
        gatk_path=args.gatk_path,
        vcftools_path=args.vcftools_path,
        bcftools_path=args.bcftools_path,
        bgzip_path=args.bgzip_path,
        tabix_path=args.tabix_path,
        java_options=args.java_options,
        use_gatk_select=args.use_gatk_select,
        maf=args.maf,
        max_missing=args.max_missing,
        hwe=args.hwe,
        min_meanDP=args.min_meanDP,
        max_meanDP=args.max_meanDP,
        minQ=args.minQ,
        minGQ=args.minGQ,
        max_alleles=args.max_alleles,
        min_alleles=args.min_alleles,
        thin=args.thin,
        chromosome_filter=args.chromosome_filter,
        chromosome_pattern=args.chromosome_pattern,
        compress_output=not args.no_compress,
        create_index=not args.no_index,
        threads=args.threads,
        memory=args.memory,
        temp_dir=args.temp_dir,
        keep_intermediate=args.keep_intermediate,
        verbose=args.verbose
    )
    
    # 🚀 运行分析 | Run analysis
    analyzer.run_pipeline()

if __name__ == "__main__":
    main()
