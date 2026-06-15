"""
Protein Stats分析主程序模块|Protein Stats Analysis Main Module
"""

import argparse
import os
from .config import ProteinStatsConfig
from .utils import ProteinStatsLogger
from .analyzer import ProteinStatsAnalyzer


class ProteinStatsPipeline:
    """Protein Stats分析管道类|Protein Stats Analysis Pipeline Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ProteinStatsConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = ProteinStatsLogger()
        self.logger = self.logger_manager.get_logger()

        # 初始化分析器|Initialize analyzer
        self.analyzer = ProteinStatsAnalyzer(self.logger, self.config)

    def run_analysis(self):
        """运行分析|Run analysis"""
        self.logger.info("=" * 60)
        self.logger.info("开始Protein Stats分析|Starting Protein Stats Analysis")
        self.logger.info("=" * 60)

        # 1. 分析蛋白质序列|Analyze protein sequences
        self.logger.info("步骤1: 分析蛋白质序列|Step 1: Analyzing protein sequences")
        df = self.analyzer.analyze()

        # 2. 保存结果|Save results
        self.logger.info("步骤2: 保存结果|Step 2: Saving results")
        self.analyzer.save_results(df)

        # 3. 打印摘要|Print summary
        self.logger.info("步骤3: 生成统计摘要|Step 3: Generating summary")
        self.analyzer.print_summary(df)

        self.logger.info("=" * 60)
        self.logger.info("分析完成|Analysis completed")
        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Protein Stats分析工具：计算蛋白质序列理化性质|Protein Stats Analysis Tool: Calculate protein sequence properties',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--protein-fasta', required=True,
                       help='蛋白序列FASTA文件|Protein sequence FASTA file')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-file', default='protein_stats.tsv',
                       help='输出文件路径|Output file path')
    parser.add_argument('--output-format', choices=['tsv', 'csv', 'excel'],
                       default='tsv',
                       help='输出文件格式|Output file format')

    # 计算选项|Calculation options
    parser.add_argument('--no-length', action='store_false', dest='calculate_length',
                       help='不计算序列长度|Do not calculate sequence length')
    parser.add_argument('--no-mw', action='store_false', dest='calculate_mw',
                       help='不计算分子量|Do not calculate molecular weight')
    parser.add_argument('--no-pi', action='store_false', dest='calculate_pi',
                       help='不计算等电点|Do not calculate isoelectric point')

    # 高级选项|Advanced options
    parser.add_argument('--aa-composition', action='store_true',
                       help='计算氨基酸组成|Calculate amino acid composition')
    parser.add_argument('--instability-index', action='store_true',
                       help='计算不稳定指数|Calculate instability index')
    parser.add_argument('--gravy', action='store_true',
                       help='计算脂肪指数(疏水性)|Calculate gravy (hydropathy)')
    parser.add_argument('--aromaticity', action='store_true',
                       help='计算芳香性|Calculate aromaticity')

    args = parser.parse_args()

    # 创建分析器并运行|Create pipeline and run
    pipeline = ProteinStatsPipeline(
        protein_fasta=args.protein_fasta,
        output_file=args.output_file,
        output_format=args.output_format,
        calculate_length=args.calculate_length,
        calculate_mw=args.calculate_mw,
        calculate_pi=args.calculate_pi,
        calculate_aa_composition=args.aa_composition,
        calculate_instability_index=args.instability_index,
        calculate_gravy=args.gravy,
        calculate_aromacity=args.aromaticity
    )

    pipeline.run_analysis()


if __name__ == "__main__":
    main()
