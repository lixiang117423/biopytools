"""
GFF3工具主程序模块 | GFF3 Tools Main Module
"""

import argparse
import sys
from .config import GFFConfig
from .utils import GFFLogger
from .gff_extractor import GFFExtractor
from .results import SummaryGenerator

class GFFAnalyzer:
    """GFF3分析主类 | Main GFF3 Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = GFFConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = GFFLogger(self.config.output_file)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化处理器 | Initialize processors
        self.gff_extractor = GFFExtractor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_extraction(self):
        """
        运行完整的提取流程 | Run complete extraction pipeline
        """
        try:
            self.logger.info("开始GFF3基因转录本提取分析 | Starting GFF3 gene transcript extraction analysis")
            self.logger.info(f"输入文件 | Input file: {self.config.gff3_file}")
            self.logger.info(f"输出文件 | Output file: {self.config.output_file}")
            
            # 提取基因转录本信息 | Extract gene transcript information
            transcript_data = self.gff_extractor.extract_gene_transcript_info()
            
            # 写入结果 | Write results
            self.gff_extractor.write_results(transcript_data)
            
            # 生成总结报告 | Generate summary report
            self.summary_generator.generate_summary_report(transcript_data)
            
            # 打印摘要 | Print summary
            self.gff_extractor.print_summary(transcript_data)
            
            self.logger.info("提取完成 | Extraction completed successfully")
            
        except Exception as e:
            self.logger.error(f"提取失败 | Extraction failed: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="从GFF3文件中为每个转录本提取整合的基因和转录本信息 | Extract integrated gene and transcript information for each transcript from GFF3 files",
        epilog="示例 | Example: parse_gene_info -g input.gff3 -o gene_transcript_info.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument(
        '--gff3', '-g',
        required=True,
        help="输入的GFF3文件路径 | Input GFF3 file path"
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help="输出的TSV文件路径 | Output TSV file path"
    )
    
    # 可选参数 | Optional arguments
    parser.add_argument(
        '--gene-type',
        default='gene',
        help="基因特征类型 | Gene feature type"
    )
    
    parser.add_argument(
        '--transcript-types',
        nargs='+',
        default=['mRNA', 'transcript'],
        help="转录本特征类型列表 | Transcript feature types list"
    )
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = GFFAnalyzer(
        gff3_file=args.gff3,
        output_file=args.output,
        gene_type=args.gene_type,
        transcript_types=set(args.transcript_types)
    )
    
    analyzer.run_extraction()

if __name__ == "__main__":
    main()
