"""
序列提取主程序模块 | Sequence Extraction Main Module
"""

import argparse
import sys
import os
from .config import SequenceConfig
from .utils import SequenceLogger, check_dependencies
from .data_processing import GenomeProcessor, VCFProcessor, SequenceBuilder
from .results import SequenceExporter, StatisticsGenerator, SummaryGenerator

class SequenceExtractor:
    """序列提取主类 | Main Sequence Extractor Class"""
    
    def __init__(self, **kwargs):
        # 检查依赖 | Check dependencies
        check_dependencies()
        
        # 初始化配置 | Initialize configuration
        self.config = SequenceConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = SequenceLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化处理器 | Initialize processors
        self.genome_processor = GenomeProcessor(self.config, self.logger)
        self.vcf_processor = VCFProcessor(self.config, self.logger)
        self.sequence_builder = SequenceBuilder(self.config, self.logger)
        self.sequence_exporter = SequenceExporter(self.config, self.logger)
        self.statistics_generator = StatisticsGenerator(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_extraction(self):
        """运行序列提取流程 | Run sequence extraction pipeline"""
        try:
            self.logger.info("开始序列提取流程 | Starting sequence extraction pipeline")
            self.logger.info(f"提取区间 | Extraction region: {self.config.chrom}:{self.config.start}-{self.config.end}")
            
            # 步骤1: 打开文件 | Step 1: Open files
            self.logger.info("步骤1: 打开输入文件 | Step 1: Opening input files")
            if not self.genome_processor.open_genome():
                self.logger.error("无法打开基因组文件 | Cannot open genome file")
                return False
            
            if not self.vcf_processor.open_vcf():
                self.logger.error("无法打开VCF文件 | Cannot open VCF file")
                return False
            
            # 步骤2: 获取参考序列 | Step 2: Get reference sequence
            self.logger.info("步骤2: 获取参考序列 | Step 2: Getting reference sequence")
            reference_seq = self.genome_processor.get_reference_sequence()
            if not reference_seq:
                self.logger.error("无法获取参考序列 | Cannot get reference sequence")
                return False
            
            sequence_length = len(reference_seq)
            self.logger.info(f"参考序列长度 | Reference sequence length: {sequence_length} bp")
            
            # 步骤3: 获取样品信息 | Step 3: Get sample information
            self.logger.info("步骤3: 获取样品信息 | Step 3: Getting sample information")
            sample_names = self.vcf_processor.get_sample_names()
            if not sample_names:
                self.logger.error("未找到样品信息 | No sample information found")
                return False
            
            self.logger.info(f"找到 {len(sample_names)} 个样品 | Found {len(sample_names)} samples")
            
            # 步骤4: 获取变异信息 | Step 4: Get variant information
            self.logger.info("步骤4: 获取变异信息 | Step 4: Getting variant information")
            variants = self.vcf_processor.get_variants_in_region()
            variants_count = len(variants)
            self.logger.info(f"找到 {variants_count} 个变异 | Found {variants_count} variants")
            
            # 步骤5: 构建样品序列 | Step 5: Build sample sequences
            self.logger.info("步骤5: 构建样品序列 | Step 5: Building sample sequences")
            sample_sequences = self.sequence_builder.build_sample_sequences(reference_seq, variants)
            if not sample_sequences:
                self.logger.error("无法构建样品序列 | Cannot build sample sequences")
                return False
            
            # 步骤6: 导出序列 | Step 6: Export sequences
            self.logger.info("步骤6: 导出序列 | Step 6: Exporting sequences")
            if not self.sequence_exporter.export_sequences(sample_sequences, reference_seq):
                self.logger.error("序列导出失败 | Sequence export failed")
                return False
            
            # 步骤7: 生成统计报告 | Step 7: Generate statistics report
            self.logger.info("步骤7: 生成统计报告 | Step 7: Generating statistics report")
            self.statistics_generator.generate_statistics_report(
                sample_sequences, reference_seq, variants_count
            )
            
            # 步骤8: 生成总结报告 | Step 8: Generate summary report
            self.logger.info("步骤8: 生成总结报告 | Step 8: Generating summary report")
            output_files = [
                f"{self.config.chrom}_{self.config.start}_{self.config.end}_sequences.{self.config.export_format}",
                f"{self.config.chrom}_{self.config.start}_{self.config.end}_statistics.txt",
                "extraction_summary.txt",
                "sequence_extraction.log"
            ]
            
            self.summary_generator.generate_summary_report(
                len(sample_sequences), variants_count, sequence_length, output_files
            )
            
            self.logger.info("序列提取完成！| Sequence extraction completed!")
            self.logger.info(f"输出文件位于 | Output files in: {self.config.output_dir}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"序列提取过程中发生错误 | Error during sequence extraction: {e}")
            return False
        
        finally:
            # 清理资源 | Cleanup resources
            self.genome_processor.close_genome()
            self.vcf_processor.close_vcf()

def main():
    """主函数 | Main function"""
    # 防止被误认为是setup.py，清理参数
    if len(sys.argv) > 1 and sys.argv[1] in ['build', 'install', 'sdist', 'bdist_wheel', '--help-commands']:
        print("Error: This is a sequence extraction script, not a setup script.")
        print("Usage: python -m sequence_toolkit.main -h")
        return 1
    
    # 设置程序名称
    prog_name = os.path.basename(sys.argv[0])
    if prog_name in ['setup.py', '__main__.py']:
        prog_name = 'sequence_extractor'
    
    parser = argparse.ArgumentParser(
        prog=prog_name,
        description='从VCF文件和基因组文件中提取特定区间的序列变异信息 | Extract sequence variation information from VCF and genome files for specific regions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-v', '--vcf', required=True,
                       help='VCF文件路径 | VCF file path')
    parser.add_argument('-g', '--genome', required=True,
                       help='基因组FASTA文件路径 | Genome FASTA file path')
    parser.add_argument('-c', '--chrom', required=True,
                       help='染色体名称 | Chromosome name (e.g., chr1, 1)')
    parser.add_argument('-s', '--start', type=int, required=True,
                       help='起始位置 (1-based) | Start position (1-based)')
    parser.add_argument('-e', '--end', type=int, required=True,
                       help='结束位置 (1-based, inclusive) | End position (1-based, inclusive)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output-dir', default='./sequence_output',
                       help='输出目录 | Output directory')
    parser.add_argument('--format', choices=['tab', 'fasta', 'csv'], default='tab',
                       help='输出格式 | Output format')
    parser.add_argument('--second-allele', action='store_true',
                       help='使用第二个等位基因而不是第一个 | Use second allele instead of first')
    parser.add_argument('--no-reference', action='store_true',
                       help='不包含参考序列 | Do not include reference sequence')
    parser.add_argument('--min-qual', type=int,
                       help='最小质量值过滤 | Minimum quality filter')
    parser.add_argument('--samples',
                       help='指定样品列表文件或逗号分隔的样品名称 | Sample list file or comma-separated sample names')
    parser.add_argument('--exclude-samples',
                       help='排除样品列表文件或逗号分隔的样品名称 | Exclude sample list file or comma-separated sample names')
    
    # 解析参数时捕获帮助请求
    try:
        args = parser.parse_args()
    except SystemExit as e:
        # 如果是帮助请求 (-h, --help)，正常退出
        return e.code
    
    # 处理样品列表 | Process sample lists
    sample_list = None
    exclude_samples = None
    
    if args.samples:
        if ',' in args.samples:
            sample_list = [s.strip() for s in args.samples.split(',')]
        else:
            try:
                with open(args.samples, 'r') as f:
                    sample_list = [line.strip() for line in f if line.strip()]
            except:
                sample_list = [args.samples]
    
    if args.exclude_samples:
        if ',' in args.exclude_samples:
            exclude_samples = [s.strip() for s in args.exclude_samples.split(',')]
        else:
            try:
                with open(args.exclude_samples, 'r') as f:
                    exclude_samples = [line.strip() for line in f if line.strip()]
            except:
                exclude_samples = [args.exclude_samples]
    
    # 创建提取器并运行 | Create extractor and run
    try:
        extractor = SequenceExtractor(
            vcf_file=args.vcf,
            genome_file=args.genome,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
            output_dir=args.output_dir,
            export_format=args.format,
            use_first_allele=not args.second_allele,
            include_reference=not args.no_reference,
            min_qual=args.min_qual,
            sample_list=sample_list,
            exclude_samples=exclude_samples
        )
        
        success = extractor.run_extraction()
        return 0 if success else 1
        
    except Exception as e:
        print(f"程序执行错误 | Program execution error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
