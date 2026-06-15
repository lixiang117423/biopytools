"""
Resistify主程序模块|Resistify Main Module
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from .config import ResistifyConfig
from .parser import ResistifyParser
from .sequence_extractor import SequenceExtractor
from .utils import ResistifyLogger, CommandRunner, build_conda_command


class ResistifyPipeline:
    """Resistify主流程|Resistify Pipeline"""

    def __init__(self, config: ResistifyConfig, logger):
        self.config = config
        self.logger = logger
        self.parser = ResistifyParser(logger, config)
        self.extractor = SequenceExtractor(logger, config)
        self.cmd_runner = CommandRunner(logger)

    def run(self):
        """运行完整的分析流程|Run complete analysis pipeline"""
        if self.config.is_directory:
            return self._run_batch()
        return self._run_single()

    def _run_single(self):
        """单文件分析流程|Single file analysis pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("开始Resistify分析|Starting Resistify analysis")
        self.logger.info("=" * 60)

        # 1. 运行Resistify（如果未跳过）|Run Resistify (if not skipped)
        if not self.config.skip_resistify:
            if not self._run_resistify(self.config.input_file, self.config.output_dir):
                self.logger.error("Resistify运行失败|Resistify run failed")
                return None

        # 2. 解析结果|Parse results
        input_dir = self.config.resistify_output_dir
        results_df = self.parser.parse_results(input_dir)
        domains_df = self.parser.parse_domains(input_dir)

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

        self.logger.info("=" * 60)
        self.logger.info("Resistify分析完成|Resistify analysis completed")
        self.logger.info(f"共处理|Total processed: {len(results_df)}条NLR基因")
        self.logger.info(f"筛选后|After filtering: {len(filtered_df)}条NLR基因")
        self.logger.info("=" * 60)

        return summary_df

    def _run_batch(self):
        """批量分析流程|Batch analysis pipeline"""
        tasks = self.config.get_batch_tasks()

        self.logger.info("=" * 60)
        self.logger.info("开始Resistify批量分析|Starting Resistify batch analysis")
        self.logger.info(f"输入目录|Input directory: {self.config.input_file}")
        self.logger.info(f"共发现|Found {len(tasks)}个样本")
        self.logger.info("=" * 60)

        all_summaries = []
        total_results = 0
        total_filtered = 0

        for sample_name, input_path, resistify_dir in tasks:
            self.logger.info("-" * 60)
            self.logger.info(f"处理样本|Processing sample: {sample_name}")
            self.logger.info(f"输入路径|Input path: {input_path}")

            # 运行Resistify|Run Resistify
            if not self.config.skip_resistify:
                if not self._run_resistify(input_path, resistify_dir):
                    self.logger.error(f"样本{sample_name}的Resistify运行失败|Resistify run failed for sample {sample_name}")
                    continue

            # 解析结果|Parse results
            try:
                results_df = self.parser.parse_results(resistify_dir)
                domains_df = self.parser.parse_domains(resistify_dir)
                merged_df = self.parser.merge_data(results_df, domains_df)
                filtered_df = self.parser.filter_data(merged_df)

                # 添加样本列|Add sample column
                filtered_df['Sample'] = sample_name

                summary_df = self.parser.prepare_summary_table(filtered_df)
                all_summaries.append(summary_df)

                # 提取序列|Extract sequences
                self.extractor.extract_and_write(
                    filtered_df,
                    resistify_output_dir=resistify_dir,
                    output_dir=resistify_dir
                )

                total_results += len(results_df)
                total_filtered += len(filtered_df)

                self.logger.info(
                    f"样本{sample_name}: 处理{len(results_df)}条，筛选后{len(filtered_df)}条|"
                    f"Sample {sample_name}: processed {len(results_df)}, filtered {len(filtered_df)}"
                )
            except Exception as e:
                self.logger.error(f"样本{sample_name}处理失败|Sample {sample_name} processing failed: {e}")

        # 合并并保存汇总结果|Merge and save summary results
        if all_summaries:
            merged_summary = pd.concat(all_summaries, ignore_index=True)
            self._save_results(merged_summary)
        else:
            merged_summary = None

        self.logger.info("=" * 60)
        self.logger.info("Resistify批量分析完成|Resistify batch analysis completed")
        self.logger.info(f"共处理|Total processed: {total_results}条NLR基因")
        self.logger.info(f"筛选后|After filtering: {total_filtered}条NLR基因")
        self.logger.info("=" * 60)

        return merged_summary

    def _run_resistify(self, input_file, output_dir):
        """运行Resistify工具|Run Resistify tool"""
        self.logger.info("运行Resistify|Running Resistify")

        cmd_list = build_conda_command(
            self.config.resistify_path,
            ['nlr', input_file, '-o', output_dir]
        )
        cmd_str = ' '.join(cmd_list)

        return self.cmd_runner.run(cmd_str, "Resistify NLR分析|Resistify NLR analysis")

    def _save_results(self, df):
        """
        保存结果到文件|Save results to files

        Args:
            df: 结果DataFrame|Result DataFrame
        """
        self.logger.info("保存结果|Saving results")

        if self.config.output_tsv:
            tsv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.tsv"
            df.to_csv(tsv_file, index=False, sep='\t', encoding='utf-8')
            self.logger.info(f"TSV已保存|TSV saved: {tsv_file}")

        if self.config.output_csv:
            csv_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.csv"
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            self.logger.info(f"CSV已保存|CSV saved: {csv_file}")

        if self.config.output_excel:
            try:
                excel_file = Path(self.config.output_dir) / f"{self.config.output_prefix}.xlsx"
                df.to_excel(excel_file, index=False, engine='openpyxl')
                self.logger.info(f"Excel已保存|Excel saved: {excel_file}")
            except ImportError:
                self.logger.warning("openpyxl未安装，跳过Excel输出|openpyxl not installed, skipping Excel output")
            except Exception as e:
                self.logger.error(f"Excel保存失败|Excel save failed: {e}")


def run_resistify(config: ResistifyConfig, logger):
    """
    运行Resistify的便捷函数|Convenience function to run Resistify

    Args:
        config: 配置对象|Configuration object
        logger: 日志对象|Logger object

    Returns:
        DataFrame: 结果DataFrame|Result DataFrame
    """
    pipeline = ResistifyPipeline(config, logger)
    return pipeline.run()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Resistify NLR分析工具：运行Resistify并解析结果|Resistify NLR Analysis Tool: Run Resistify and parse results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='输入蛋白质FASTA文件或目录|Input protein FASTA file or directory')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir', default='./resistify_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='并行线程数|Number of parallel threads')

    # Resistify工具参数|Resistify tool parameters
    parser.add_argument('--resistify-path',
                        help='Resistify可执行文件路径|Resistify executable path')
    parser.add_argument('--skip-resistify', action='store_true',
                        help='跳过Resistify运行，仅解析已有结果|Skip Resistify run, only parse existing results')

    # 解析结果输出配置|Parse results output configuration
    parser.add_argument('--output-prefix', default='resistify_results',
                        help='解析结果文件前缀|Parsed results file prefix')
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
    logger_manager = ResistifyLogger()
    logger = logger_manager.get_logger()

    # 创建配置|Create configuration
    config_kwargs = {
        'input_file': args.input,
        'output_dir': args.output_dir,
        'threads': args.threads,
        'output_prefix': args.output_prefix,
        'output_tsv': args.output_tsv,
        'output_csv': args.output_csv,
        'output_excel': args.output_excel,
        'extract_nlr_sequences': args.extract_nlr_sequences,
        'extract_nbarc_sequences': args.extract_nbarc_sequences,
        'filter_classification': args.filter_classification,
        'min_length': args.min_length,
        'max_length': args.max_length,
        'min_lrr_length': args.min_lrr_length,
        'include_motifs': args.include_motifs,
        'skip_resistify': args.skip_resistify,
    }

    # resistify_path仅在用户指定非默认值时传递|Only pass resistify_path if user specified non-default
    if args.resistify_path:
        config_kwargs['resistify_path'] = args.resistify_path

    config = ResistifyConfig(**config_kwargs)

    # 验证配置|Validate configuration
    try:
        config.validate()
    except ValueError as e:
        logger.error(f"配置错误|Configuration error:\n{e}")
        return 1

    # 运行分析|Run analysis
    try:
        run_resistify(config, logger)
        return 0
    except Exception as e:
        logger.error(f"运行失败|Run failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
