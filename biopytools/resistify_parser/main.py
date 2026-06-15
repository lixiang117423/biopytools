"""
Resistify Parser主程序|Resistify Parser Main Module
"""

import pandas as pd
import argparse
from pathlib import Path
from .config import ResistifyParserConfig
from .parser import ResistifyParser
from .sequence_extractor import SequenceExtractor
from .utils import ResistifyParserLogger


class ResistifyParserPipeline:
    """Resistify Parser主流程|Resistify Parser Pipeline"""

    def __init__(self, config: ResistifyParserConfig, logger):
        self.config = config
        self.logger = logger
        self.parser = ResistifyParser(logger, config)
        self.extractor = SequenceExtractor(logger, config)

    def run(self):
        """运行完整的解析流程|Run complete parsing pipeline"""
        self.logger.info("="*60)
        self.logger.info("开始Resistify结果解析|Starting Resistify results parsing")
        self.logger.info("="*60)

        # 1. 解析results.tsv|Parse results.tsv
        results_df = self.parser.parse_results()

        # 2. 解析domains.tsv|Parse domains.tsv
        domains_df = self.parser.parse_domains()

        # 3. 合并数据|Merge data
        merged_df = self.parser.merge_data(results_df, domains_df)

        # 4. 过滤数据|Filter data
        filtered_df = self.parser.filter_data(merged_df)

        # 5. 准备汇总表格|Prepare summary table
        summary_df = self.parser.prepare_summary_table(filtered_df)

        # 6. 提取序列|Extract sequences
        self.extractor.extract_and_write(filtered_df)

        # 7. 保存结果|Save results
        self._save_results(summary_df)

        self.logger.info("="*60)
        self.logger.info("Resistify结果解析完成|Resistify results parsing completed")
        self.logger.info(f"共处理|Total processed: {len(results_df)}条NLR基因")
        self.logger.info(f"筛选后|After filtering: {len(filtered_df)}条NLR基因")
        self.logger.info("="*60)

        return summary_df

    def _save_results(self, df: pd.DataFrame):
        """
        保存结果到文件|Save results to files

        Args:
            df: 结果DataFrame|Result DataFrame
        """
        self.logger.info("保存结果|Saving results")

        # 保存TSV|Save TSV
        if self.config.output_tsv:
            tsv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.tsv"
            df.to_csv(tsv_file, index=False, sep='\t', encoding='utf-8')
            self.logger.info(f"TSV已保存|TSV saved: {tsv_file}")

        # 保存CSV|Save CSV
        if self.config.output_csv:
            csv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.csv"
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            self.logger.info(f"CSV已保存|CSV saved: {csv_file}")

        # 保存Excel|Save Excel
        if self.config.output_excel:
            try:
                excel_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.xlsx"
                df.to_excel(excel_file, index=False, engine='openpyxl')
                self.logger.info(f"Excel已保存|Excel saved: {excel_file}")
            except ImportError:
                self.logger.warning("openpyxl未安装，跳过Excel输出|openpyxl not installed, skipping Excel output")
            except Exception as e:
                self.logger.error(f"Excel保存失败|Excel save failed: {e}")


def run_resistify_parser(config: ResistifyParserConfig, logger):
    """
    运行Resistify Parser的便捷函数|Convenience function to run Resistify Parser

    Args:
        config: 配置对象|Configuration object
        logger: 日志对象|Logger object

    Returns:
        pd.DataFrame: 结果DataFrame|Result DataFrame
    """
    pipeline = ResistifyParserPipeline(config, logger)
    return pipeline.run()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Resistify Parser工具：解析Resistify软件输出结果|Resistify Parser Tool: Parse Resistify software output results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input-dir', required=True,
                       help='Resistify输出目录|Resistify output directory')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-prefix', default='resistify_results',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--output-dir', default='.',
                       help='输出目录|Output directory')
    parser.add_argument('--no-tsv', action='store_false', dest='output_tsv',
                       help='不输出TSV文件|Do not output TSV file')
    parser.add_argument('--no-csv', action='store_false', dest='output_csv',
                       help='不输出CSV文件|Do not output CSV file')
    parser.add_argument('--no-excel', action='store_false', dest='output_excel',
                       help='不输出Excel文件|Do not output Excel file')

    # 序列提取选项|Sequence extraction options
    parser.add_argument('--extract-nlr', action='store_true', dest='extract_nlr_sequences',
                       help='提取NLR序列|Extract NLR sequences')
    parser.add_argument('--extract-nbarc', action='store_true', dest='extract_nbarc_sequences',
                       help='提取NB-ARC序列|Extract NB-ARC sequences')

    # 筛选选项|Filtering options
    parser.add_argument('--filter-classification', dest='filter_classification',
                       help='按分类筛选(如TN, CNL, NL等)|Filter by classification (e.g., TN, CNL, NL)')
    parser.add_argument('--min-length', type=int, dest='min_length',
                       help='最小序列长度|Minimum sequence length')
    parser.add_argument('--max-length', type=int, dest='max_length',
                       help='最大序列长度|Maximum sequence length')
    parser.add_argument('--min-lrr-length', type=int, dest='min_lrr_length',
                       help='最小LRR长度|Minimum LRR length')

    # 其他选项|Other options
    parser.add_argument('--include-motifs', action='store_true',
                       help='包含motifs详情|Include motifs details')

    args = parser.parse_args()

    # 初始化日志|Initialize logging
    logger_manager = ResistifyParserLogger()
    logger = logger_manager.get_logger()

    # 创建配置|Create configuration
    config = ResistifyParserConfig(
        input_dir=args.input_dir,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        output_tsv=args.output_tsv,
        output_csv=args.output_csv,
        output_excel=args.output_excel,
        extract_nlr_sequences=args.extract_nlr_sequences,
        extract_nbarc_sequences=args.extract_nbarc_sequences,
        filter_classification=args.filter_classification,
        min_length=args.min_length,
        max_length=args.max_length,
        min_lrr_length=args.min_lrr_length,
        include_motifs=args.include_motifs
    )

    # 验证配置|Validate configuration
    try:
        config.validate()
    except ValueError as e:
        logger.error(f"配置错误|Configuration error:\n{e}")
        return 1

    # 运行分析|Run analysis
    try:
        run_resistify_parser(config, logger)
        return 0
    except Exception as e:
        logger.error(f"运行失败|Run failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
