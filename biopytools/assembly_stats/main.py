"""
基因组装配统计主程序模块|Genome Assembly Statistics Main Module
"""

import argparse
import sys
import pandas as pd
from pathlib import Path

from .config import AssemblyStatsConfig
from .utils import AssemblyStatsLogger
from .stats_analyzer import AssemblyStatsAnalyzer


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="基因组装配序列长度统计工具|Genome Assembly Sequence Length Statistics Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genome.fa
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='输入文件或文件夹|Input file or directory path')

    parser.add_argument('-l', '--min-length',
                        type=int,
                        default=1,
                        help='最小序列长度过滤|Minimum sequence length cutoff')

    parser.add_argument('-s',
                        action='store_true',
                        help='Grep友好输出格式|Print grep-friendly output')

    parser.add_argument('-t',
                        action='store_true',
                        help='Tab分隔输出|Print tab-delimited output')

    parser.add_argument('-u',
                        action='store_true',
                        help='Tab分隔输出且无header|Print tab-delimited output without header')

    parser.add_argument('-o', '--output-dir',
                        default='./assembly_stats_output',
                        help='输出目录|Output directory')

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 1.0.0')

    return parser.parse_args()


class AssemblyStatsRunner:
    """基因组装配统计运行器|Genome Assembly Statistics Runner"""

    def __init__(self, config: AssemblyStatsConfig, logger):
        self.config = config
        self.logger = logger
        self.analyzer = AssemblyStatsAnalyzer(config, logger)

    def run(self):
        """运行分析|Run analysis"""
        self.logger.info("开始基因组装配统计分析|Starting genome assembly statistics analysis")
        self.logger.info(f"输入文件数量|Number of input files: {len(self.config.input_files)}")

        # 分析所有文件|Analyze all files
        all_stats = []
        for file_path in self.config.input_files:
            stats = self.analyzer.analyze_file(file_path)
            if stats:
                all_stats.append(stats)

                # 打印到标准输出|Print to stdout
                output = self.analyzer.format_output(stats)
                print(output)
                print()

        # 生成报告文件|Generate report files
        if all_stats:
            self._generate_reports(all_stats)

        self.logger.info("分析完成|Analysis completed")

    def _generate_reports(self, all_stats: list):
        """
       生成报告文件|Generate report files

        Args:
            all_stats: 所有统计结果|All statistics results
        """
        # 辅助函数：格式化大数字|Helper function: format large numbers
        def format_number(num: int) -> str:
            """格式化数字，以M为单位|Format number in M unit"""
            if num >= 1_000_000:
                return f"{num / 1_000_000:.2f}M"
            return str(num)

        # 竖着格式：每个文件一行数据，转置为键值对|Vertical format: transpose to key-value pairs
        if self.config.vertical_format:
            # 为每个文件生成独立的竖着格式报告|Generate separate vertical report for each file
            for stats in all_stats:
                file_name = stats['file_name']

                # 创建竖着格式的数据|Create vertical format data
                vertical_data = [
                    ['File|文件', file_name],
                    ['Sum|总和', format_number(stats['sum'])],
                    ['N|序列数', stats['n']],
                    ['Average|平均长度', format_number(int(stats['ave']))],
                    ['Largest|最大长度', format_number(stats['largest'])],
                ]

                # 添加Nx统计|Add Nx statistics
                for percentage in [50, 60, 70, 80, 90, 100]:
                    nx_key = f'N{percentage}'
                    nx_n_key = f'N{percentage}_n'
                    vertical_data.append([nx_key, format_number(stats[nx_key])])
                    vertical_data.append([f'{nx_key}_count|{nx_key}_序列数', stats[nx_n_key]])

                # Gap和N统计|Gap and N statistics
                vertical_data.append(['N_count|N碱基数', format_number(stats['n_count'])])
                vertical_data.append(['Gaps|Gap数', stats['gaps']])

                # 创建DataFrame|Create DataFrame
                df = pd.DataFrame(vertical_data, columns=['指标|Metric', '值|Value'])

                # 输出CSV文件（每个文件单独一个）|Output CSV file (separate for each file)
                if self.config.generate_csv:
                    base_name = Path(file_name).stem
                    csv_file = self.config.output_path / f"{base_name}_stats.csv"
                    df.to_csv(csv_file, index=False, encoding='utf-8')
                    self.logger.info(f"CSV报告已保存|CSV report saved: {csv_file}")

                # 输出Excel文件|Output Excel file
                if self.config.generate_xlsx:
                    base_name = Path(file_name).stem
                    xlsx_file = self.config.output_path / f"{base_name}_stats.xlsx"
                    df.to_excel(xlsx_file, index=False, engine='openpyxl')
                    self.logger.info(f"Excel报告已保存|Excel report saved: {xlsx_file}")

            # 同时生成汇总报告（横着格式，包含所有文件）|Also generate summary report (horizontal, all files)
            summary_data = []
            for stats in all_stats:
                row = {
                    'File|文件': stats['file_name'],
                    'Sum|总和': stats['sum'],
                    'N|序列数': stats['n'],
                    'Average|平均长度': f"{stats['ave']:.2f}",
                    'Largest|最大长度': stats['largest'],
                }

                # 添加Nx统计|Add Nx statistics
                for percentage in [50, 60, 70, 80, 90, 100]:
                    row[f'N{percentage}'] = stats[f'N{percentage}']
                    row[f'N{percentage}_count|N{percentage}_序列数'] = stats[f'N{percentage}_n']

                row['N_count|N碱基数'] = stats['n_count']
                row['Gaps|Gap数'] = stats['gaps']

                summary_data.append(row)

            summary_df = pd.DataFrame(summary_data)

            # 输出汇总CSV文件|Output summary CSV file
            if self.config.generate_csv:
                summary_csv = self.config.output_path / f"{self.config.output_prefix}_summary.csv"
                summary_df.to_csv(summary_csv, index=False, encoding='utf-8')
                self.logger.info(f"汇总CSV报告已保存|Summary CSV report saved: {summary_csv}")

            # 输出汇总Excel文件|Output summary Excel file
            if self.config.generate_xlsx:
                summary_xlsx = self.config.output_path / f"{self.config.output_prefix}_summary.xlsx"
                summary_df.to_excel(summary_xlsx, index=False, engine='openpyxl')
                self.logger.info(f"汇总Excel报告已保存|Summary Excel report saved: {summary_xlsx}")

        else:
            # 原有横着格式|Original horizontal format
            df_data = []
            for stats in all_stats:
                row = {
                    'File|文件': stats['file_name'],
                    'Sum|总和': stats['sum'],
                    'N|序列数': stats['n'],
                    'Average|平均长度': f"{stats['ave']:.2f}",
                    'Largest|最大长度': stats['largest'],
                }

                # 添加Nx统计|Add Nx statistics
                for percentage in [50, 60, 70, 80, 90, 100]:
                    row[f'N{percentage}'] = stats[f'N{percentage}']
                    row[f'N{percentage}_count|N{percentage}_序列数'] = stats[f'N{percentage}_n']

                row['N_count|N碱基数'] = stats['n_count']
                row['Gaps|Gap数'] = stats['gaps']

                df_data.append(row)

            df = pd.DataFrame(df_data)

            # 输出CSV文件|Output CSV file
            if self.config.generate_csv:
                csv_file = self.config.output_path / f"{self.config.output_prefix}.csv"
                df.to_csv(csv_file, index=False, encoding='utf-8')
                self.logger.info(f"CSV报告已保存|CSV report saved: {csv_file}")

            # 输出Excel文件|Output Excel file
            if self.config.generate_xlsx:
                xlsx_file = self.config.output_path / f"{self.config.output_prefix}.xlsx"
                df.to_excel(xlsx_file, index=False, engine='openpyxl')
                self.logger.info(f"Excel报告已保存|Excel report saved: {xlsx_file}")


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = AssemblyStatsConfig(
            input_path=args.input,
            output_dir=args.output_dir,
            min_length=args.min_length,
            grep_friendly=args.s,
            tab_delimited=args.t,
            no_header=args.u
        )
        config.validate()

        # 创建日志|Create logger
        log_file = Path(config.output_dir) / "assembly_stats.log"
        logger_manager = AssemblyStatsLogger(str(log_file))
        logger = logger_manager.get_logger()

        # 运行分析|Run analysis
        runner = AssemblyStatsRunner(config, logger)
        runner.run()

        sys.exit(0)

    except KeyboardInterrupt:
        print("用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
