"""
VCF筛选主程序模块 | VCF Filtering Main Module
"""

import argparse
import sys
import os
import time
from typing import Union, List, Optional
from .config import VCFFilterConfig
from .core import VCFFilter
from .utils import FilterLogger, check_dependencies

class VCFFilterMain:
    """VCF筛选主类 | Main VCF Filter Class"""
    
    def __init__(self, **kwargs):
        # 处理参数并创建配置 | Process parameters and create config
        self.config = VCFFilterConfig(**kwargs)
        
        # 使用标准日志器
        self.logger_manager = FilterLogger(self.config.output_dir, verbose=self.config.verbose)
        self.logger = self.logger_manager.get_logger()
    
    def run_analysis(self):
        """运行完整的VCF筛选流程 | Run complete VCF filtering pipeline"""
        try:
            self.logger.info("开始VCF文件筛选流程 | Starting VCF file filtering pipeline")
            
            # 检查依赖 | Check dependencies
            check_dependencies(self.logger)
            
            # 创建筛选器 | Create filter
            vcf_filter = VCFFilter(self.config.vcf_file)
            
            # 执行筛选 | Execute filtering
            output_file = vcf_filter.filter_vcf(
                chr_name=self.config.chr_name,
                start=self.config.start,
                end=self.config.end,
                output_file=self.config.output_file,
                convert_format=self.config.convert_format,
                plink_path=self.config.plink_path,
                allow_extra_chr=self.config.allow_extra_chr,
                min_maf=self.config.min_maf,
                max_missing=self.config.max_missing,
                quality_threshold=self.config.quality_threshold,
                min_depth=self.config.min_depth,
                max_depth=self.config.max_depth,
                keep_samples=self.config.keep_samples,
                remove_samples=self.config.remove_samples,
                biallelic_only=self.config.biallelic_only,
                remove_indels=self.config.remove_indels,
                skip_validation=self.config.skip_validation,
                verbose=self.config.verbose
            )
            
            # 完成信息
            self.logger.info("\\n" + "=" * 60)
            self.logger.info("VCF文件筛选完成！| VCF file filtering completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"输出文件: {output_file} | Output file: {output_file}")
            
            return output_file
            
        except Exception as e:
            self.logger.error(f"筛选流程在执行过程中意外终止 | "
                            f"Filtering pipeline terminated unexpectedly: {e}")
            raise

def filter_vcf_file(vcf_file: str, 
                   chr_name: Union[str, List[str]], 
                   start: Optional[int] = None, 
                   end: Optional[int] = None,
                   output_file: Optional[str] = None,
                   convert_format: bool = False,
                   **kwargs) -> str:
    """
    便捷函数：筛选VCF文件（保持原始接口）
    
    Args:
        vcf_file: 输入VCF文件路径
        chr_name: 染色体名称
        start: 起始位置
        end: 结束位置  
        output_file: 输出文件路径
        convert_format: 是否使用plink转换格式
        **kwargs: 其他参数
    
    Returns:
        str: 输出文件路径
    
    Example:
        # 筛选1号染色体1000-2000bp区域（高性能模式）
        output = filter_vcf_file("input.vcf", "1", 1000, 2000, skip_validation=True)
        
        # 使用plink转换格式
        output = filter_vcf_file("input.vcf", ["1", "2"], convert_format=True)
    """
    
    vcf_filter = VCFFilter(vcf_file)
    return vcf_filter.filter_vcf(
        chr_name=chr_name,
        start=start,
        end=end,
        output_file=output_file,
        convert_format=convert_format,
        **kwargs
    )

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="VCF文件筛选工具 (模块化高性能版本) | VCF File Filtering Tool (Modular High-Performance Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
使用示例 | Examples:
  %(prog)s -i input.vcf -c chr1 -s 1000 -e 2000
  %(prog)s -i input.vcf -c chr1 --convert-format --maf 0.05
  %(prog)s -i input.vcf -c "chr1,chr2,chr3" --biallelic-only --verbose
  %(prog)s -i input.vcf -c chr1 --keep-samples "sample1,sample2,sample3"
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入VCF文件路径 | Input VCF file path')
    
    # 输出参数 | Output parameters
    parser.add_argument('-o', '--output',
                       help='输出VCF文件路径 | Output VCF file path')
    
    # 位置筛选参数 | Position filtering parameters
    parser.add_argument('-c', '--chr', '--chromosome', required=True,
                       help='染色体名称 (支持逗号分隔的多个染色体) | '
                            'Chromosome name(s) (comma-separated for multiple)')
    parser.add_argument('-s', '--start', type=int,
                       help='起始位置 | Start position')
    parser.add_argument('-e', '--end', type=int,
                       help='结束位置 | End position')
    
    # 格式转换参数 | Format conversion parameters
    parser.add_argument('--convert-format', action='store_true',
                       help='使用PLINK进行格式转换 | Use PLINK for format conversion')
    parser.add_argument('--plink-path', default='plink',
                       help='PLINK可执行文件路径 | PLINK executable path')
    parser.add_argument('--allow-extra-chr', action='store_true', default=True,
                       help='允许额外染色体 | Allow extra chromosomes')
    
    # 质量控制参数 | Quality control parameters
    parser.add_argument('--maf', type=float,
                       help='最小等位基因频率 | Minimum allele frequency')
    parser.add_argument('--max-missing', type=float,
                       help='最大缺失率 | Maximum missing rate')
    parser.add_argument('--quality-threshold', type=float,
                       help='质量阈值 | Quality threshold')
    parser.add_argument('--min-depth', type=int,
                       help='最小深度 | Minimum depth')
    parser.add_argument('--max-depth', type=int,
                       help='最大深度 | Maximum depth')
    
    # 样本筛选参数 | Sample filtering parameters
    parser.add_argument('--keep-samples',
                       help='保留的样本名称 (逗号分隔) | '
                            'Sample names to keep (comma-separated)')
    parser.add_argument('--remove-samples',
                       help='移除的样本名称 (逗号分隔) | '
                            'Sample names to remove (comma-separated)')
    
    # 变异位点筛选参数 | Variant filtering parameters
    parser.add_argument('--keep-ids',
                       help='保留的变异位点ID (逗号分隔) | '
                            'Variant IDs to keep (comma-separated)')
    parser.add_argument('--remove-ids',
                       help='移除的变异位点ID (逗号分隔) | '
                            'Variant IDs to remove (comma-separated)')
    parser.add_argument('--biallelic-only', action='store_true',
                       help='只保留双等位基因位点 | Keep only biallelic sites')
    parser.add_argument('--remove-indels', action='store_true',
                       help='移除插入缺失变异 | Remove indel variants')
    
    # 性能优化参数 | Performance optimization parameters
    parser.add_argument('--skip-validation', action='store_true', default=True,
                       help='跳过输入验证以提高速度（默认开启）| Skip input validation for speed (default enabled)')
    parser.add_argument('--force-validation', action='store_true',
                       help='强制执行输入验证 | Force input validation')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='显示详细信息 | Show verbose information')
    
    args = parser.parse_args()
    
    # 检查输入文件 | Check input file
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在 | Error: Input file does not exist - {args.input}")
        sys.exit(1)
    
    # 处理染色体参数 | Process chromosome parameter
    chr_name = [c.strip() for c in args.chr.split(',')]
    if len(chr_name) == 1:
        chr_name = chr_name[0]
    
    # 处理样本参数 | Process sample parameters
    keep_samples = None
    if args.keep_samples:
        keep_samples = [s.strip() for s in args.keep_samples.split(',')]
    
    remove_samples = None
    if args.remove_samples:
        remove_samples = [s.strip() for s in args.remove_samples.split(',')]
    
    # 处理ID参数 | Process ID parameters
    keep_ids = None
    if args.keep_ids:
        keep_ids = [i.strip() for i in args.keep_ids.split(',')]
    
    remove_ids = None
    if args.remove_ids:
        remove_ids = [i.strip() for i in args.remove_ids.split(',')]
    
    # 处理验证参数 | Process validation parameters
    skip_validation = args.skip_validation and not args.force_validation
    
    try:
        start_time = time.time()
        
        # 创建筛选器并运行 | Create filter and run
        analyzer = VCFFilterMain(
            vcf_file=args.input,
            output_file=args.output,
            chr_name=chr_name,
            start=args.start,
            end=args.end,
            convert_format=args.convert_format,
            plink_path=args.plink_path,
            allow_extra_chr=args.allow_extra_chr,
            min_maf=args.maf,
            max_missing=args.max_missing,
            quality_threshold=args.quality_threshold,
            min_depth=args.min_depth,
            max_depth=args.max_depth,
            keep_samples=keep_samples,
            remove_samples=remove_samples,
            keep_ids=keep_ids,
            remove_ids=remove_ids,
            biallelic_only=args.biallelic_only,
            remove_indels=args.remove_indels,
            skip_validation=skip_validation,
            verbose=args.verbose
        )
        
        output_file = analyzer.run_analysis()
        
        elapsed_time = time.time() - start_time
        print(f"\n筛选完成！用时 {elapsed_time:.1f} 秒")
        print(f"输出文件: {output_file}")
        
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
