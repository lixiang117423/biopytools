"""
VCF单体型提取主程序模块 | VCF Haplotype Extraction Main Module
"""

import argparse
import sys
from typing import Optional

from .config import HaplotypeConfig
from .utils import HaplotypeLogger, CommandRunner, check_dependencies, TempFileManager
from .data_processing import VCFProcessor

class HaplotypeExtractor:
    """VCF单体型提取器 | VCF Haplotype Extractor"""
    
    def __init__(self, vcf_file: str, output_file: str, 
                 position_file: Optional[str] = None,
                 chromosome: Optional[str] = None, 
                 position: Optional[int] = None,
                 bcftools_path: str = 'bcftools', 
                 verbose: bool = True):
        """
        初始化单体型提取器 | Initialize haplotype extractor
        
        Args:
            vcf_file: VCF文件路径 | VCF file path
            output_file: 输出文件路径 | Output file path
            position_file: 位置文件路径 | Position file path
            chromosome: 染色体标识符 | Chromosome identifier
            position: SNP位置 | SNP position
            bcftools_path: bcftools路径 | bcftools path
            verbose: 详细输出 | Verbose output
        """
        self.config = HaplotypeConfig(
            vcf_file=vcf_file,
            output_file=output_file,
            position_file=position_file,
            chromosome=chromosome,
            position=position,
            bcftools_path=bcftools_path,
            verbose=verbose
        )
        
        # 验证配置 | Validate configuration
        self.config.validate()
        
        # 初始化组件 | Initialize components
        self.logger_manager = HaplotypeLogger(output_file)
        self.logger = self.logger_manager.get_logger()
        self.command_runner = CommandRunner(self.logger)
        self.temp_manager = TempFileManager(self.logger)
        self.vcf_processor = VCFProcessor(self.logger, self.command_runner, self.temp_manager)
    
    def write_matrix_output(self, samples: list, variant_data: dict, target_positions: list):
        """写入矩阵格式输出文件 | Write matrix format output file"""
        self.logger.info(f"写入矩阵格式结果到 | Writing matrix format results to: {self.config.output_file}")
        
        with open(self.config.output_file, 'w') as f:
            # 写入标题行 | Write header
            header = ['CHROM', 'POS', 'REF', 'ALT'] + samples
            f.write('\t'.join(header) + '\n')
            
            # 按照target_positions的顺序写入数据
            for chrom, pos in target_positions:
                variant_key = (chrom, pos)
                if variant_key in variant_data:
                    variant_info = variant_data[variant_key]
                    row = [
                        chrom,
                        str(pos),
                        variant_info['ref'],
                        variant_info['alt']
                    ] + variant_info['genotypes']
                    f.write('\t'.join(row) + '\n')
                else:
                    # 如果该位点在VCF中不存在，写入缺失信息
                    row = [chrom, str(pos), 'N/A', 'N/A'] + ['--'] * len(samples)
                    f.write('\t'.join(row) + '\n')
                    self.logger.warning(f"位点 {chrom}:{pos} 在VCF文件中未找到 | Position {chrom}:{pos} not found in VCF file")
        
        self.logger.info(f"输出文件已保存 | Output file saved: {self.config.output_file}")
    
    def print_summary(self, samples: list, variant_data: dict, target_positions: list):
        """打印分析摘要 | Print analysis summary"""
        if not self.config.verbose:
            return
        
        found_variants = len(variant_data)
        total_requested = len(target_positions)
        
        self.logger.info("="*60)
        self.logger.info("分析摘要 | Analysis Summary")
        self.logger.info("="*60)
        self.logger.info(f"请求的位点数 | Requested positions: {total_requested}")
        self.logger.info(f"找到的变异数 | Found variants: {found_variants}")
        self.logger.info(f"样本数量 | Number of samples: {len(samples)}")
        
        if found_variants < total_requested:
            missing_count = total_requested - found_variants
            self.logger.warning(f"缺失的位点数 | Missing positions: {missing_count}")
        
        # 统计基因型分布
        genotype_stats = {}
        for variant_info in variant_data.values():
            for gt in variant_info['genotypes']:
                genotype_stats[gt] = genotype_stats.get(gt, 0) + 1
        
        self.logger.info("基因型分布 | Genotype distribution:")
        for gt, count in sorted(genotype_stats.items()):
            percentage = (count / sum(genotype_stats.values())) * 100
            self.logger.info(f"  {gt}: {count} ({percentage:.1f}%)")
    
    def run_analysis(self):
        """运行单体型提取分析 | Run haplotype extraction analysis"""
        try:
            self.logger.info("="*60)
            self.logger.info("VCF单体型提取分析开始 | VCF Haplotype Extraction Analysis Started")
            self.logger.info("="*60)
            
            # 检查依赖 | Check dependencies
            check_dependencies(self.config, self.logger)
            
            # 获取目标位点列表 | Get target positions
            target_positions = self.config.get_positions()
            self.logger.info(f"目标位点数量 | Number of target positions: {len(target_positions)}")
            
            # 提取VCF变异 | Extract VCF variants
            extracted_vcf = self.vcf_processor.extract_variants(self.config, target_positions)
            
            # 解析VCF为矩阵格式 | Parse VCF to matrix format
            samples, variant_data = self.vcf_processor.parse_vcf_to_matrix(extracted_vcf, target_positions)
            
            if not variant_data:
                self.logger.warning("未找到任何匹配的变异位点 | No matching variants found")
                return
            
            # 写入输出 | Write output
            self.write_matrix_output(samples, variant_data, target_positions)
            
            # 打印摘要 | Print summary
            self.print_summary(samples, variant_data, target_positions)
            
            self.logger.info("="*60)
            self.logger.info("分析完成 | Analysis completed successfully")
            self.logger.info(f"结果文件 | Results file: {self.config.output_file}")
            self.logger.info("="*60)
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error occurred during analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            sys.exit(1)
        finally:
            # 清理临时文件 | Cleanup temporary files
            self.temp_manager.cleanup()

def create_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description="VCF单体型提取工具 v2.0 | VCF Haplotype Extractor v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 使用位置文件批量提取 | Batch extraction using position file
  %(prog)s -v input.vcf.gz -p positions.txt -o output.txt
  
  # 提取单个SNP | Extract single SNP
  %(prog)s -v input.vcf.gz --single -c OV12 -s 93635286 -o output.txt
  
  # 指定bcftools路径 | Specify bcftools path
  %(prog)s -v input.vcf.gz -p positions.txt -o output.txt --bcftools-path /path/to/bcftools

位置文件格式 | Position File Format:
  # 可以有表头（推荐）| Can have header (recommended)
  CHR	POS
  OV12	93635286
  OV12	93440535
  OV12	93234780
  
  # 也可以没有表头 | Can also without header
  OV12	93635286
  OV12	93440535
  OV12	93234780

输出格式 | Output Format:
  CHROM	POS	REF	ALT	Sample1	Sample2	Sample3
  OV12	93635286	C	A	CA	AA	AA
  OV12	93440535	C	A	CA	AA	CA
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("-v", "--vcf", required=True,
                       help="输入VCF文件路径 | Input VCF file path")
    parser.add_argument("-o", "--output", required=True,
                       help="输出文件路径 | Output file path")
    
    # 位置参数 (二选一) | Position parameters (choose one)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--positions", 
                       help="位置文件路径 (CHR\\tPOS格式，可选表头) | Position file path (CHR\\tPOS format, optional header)")
    group.add_argument("--single", action='store_true',
                       help="使用单个位点模式 | Use single position mode")
    
    # 单个位点参数 | Single position parameters
    parser.add_argument("-c", "--chr", 
                       help="染色体标识符 (单个位点模式) | Chromosome identifier (single position mode)")
    parser.add_argument("-s", "--site", type=int,
                       help="SNP位置 (单个位点模式) | SNP position (single position mode)")
    
    # 可选参数 | Optional arguments
    parser.add_argument("--bcftools-path", default="bcftools",
                       help="bcftools程序路径 | bcftools program path")
    parser.add_argument("--quiet", action="store_true",
                       help="静默模式，减少输出 | Quiet mode, reduce output")
    
    return parser

def main():
    """主函数 | Main function"""
    parser = create_parser()
    args = parser.parse_args()
    
    # 验证参数组合
    if args.single:
        if not args.chr or args.site is None:
            parser.error("单个位点模式需要同时指定 -c/--chr 和 -s/--site 参数")
        position_file = None
        chromosome = args.chr
        position = args.site
    else:
        position_file = args.positions
        chromosome = None
        position = None
    
    try:
        # 创建并运行分析器 | Create and run analyzer
        extractor = HaplotypeExtractor(
            vcf_file=args.vcf,
            output_file=args.output,
            position_file=position_file,
            chromosome=chromosome,
            position=position,
            bcftools_path=args.bcftools_path,
            verbose=not args.quiet
        )
        
        extractor.run_analysis()
        
    except Exception as e:
        print(f"错误 | Error: {e}", file=sys.stderr)
        sys.exit(1)
