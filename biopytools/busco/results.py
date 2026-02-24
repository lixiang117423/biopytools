"""
BUSCO结果处理模块|BUSCO Results Processing Module
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict, Any

class ResultsProcessor:
    """结果处理器|Results Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def compile_results(self, results_list: List[Dict[str, Any]]) -> pd.DataFrame:
        """编译结果为表格|Compile results to table"""
        self.logger.info("编译BUSCO分析结果|Compiling BUSCO analysis results")

        if not results_list:
            self.logger.warning("没有成功的分析结果|No successful analysis results")
            return pd.DataFrame()

        # 创建DataFrame|Create DataFrame
        df = pd.DataFrame(results_list)

        # 重新排列列顺序|Reorder columns
        column_order = [
            'sample_name',
            'lineage_name',
            'complete_percentage',
            'complete_buscos',
            'single_copy_percentage',
            'single_copy_buscos',
            'multi_copy_percentage',
            'multi_copy_buscos',
            'fragmented_percentage',
            'fragmented_buscos',
            'missing_percentage',
            'missing_buscos',
            'n_markers',
            'domain',
            'status'
        ]

        # 确保所有列都存在|Ensure all columns exist
        for col in column_order:
            if col not in df.columns:
                df[col] = 'N/A'

        df = df[column_order]

        # 设置列名（双语）|Set column names (bilingual)
        df.columns = [
            '样本名称|Sample_Name',
            '数据库|Database',
            '完整度|Complete_%',
            '完整基因数|Complete_BUSCOs',
            '单拷贝|Single_Copy_%',
            '单拷贝基因数|Single_Copy_BUSCOs',
            '多拷贝|Multi_Copy_%',
            '多拷贝基因数|Multi_Copy_BUSCOs',
            '片段化|Fragmented_%',
            '片段化基因数|Fragmented_BUSCOs',
            '缺失|Missing_%',
            '缺失基因数|Missing_BUSCOs',
            '标记基因总数|Total_Markers',
            '域|Domain',
            '状态|Status'
        ]

        self.logger.info(f"成功编译 {len(df)} 个样本的结果|Successfully compiled results for {len(df)} samples")

        return df
    
    def add_failed_samples(self, df: pd.DataFrame, failed_samples: List[str]) -> pd.DataFrame:
        """添加失败样本到结果表格|Add failed samples to results table"""
        if not failed_samples:
            return df

        self.logger.info(f"添加 {len(failed_samples)} 个失败样本|Adding {len(failed_samples)} failed samples")

        # 为失败样本创建行|Create rows for failed samples
        failed_rows = []
        for sample_name in failed_samples:
            failed_row = {
                '样本名称|Sample_Name': sample_name,
                '数据库|Database': self.config.lineage,
                '完整度|Complete_%': 'N/A',
                '完整基因数|Complete_BUSCOs': 'N/A',
                '单拷贝|Single_Copy_%': 'N/A',
                '单拷贝基因数|Single_Copy_BUSCOs': 'N/A',
                '多拷贝|Multi_Copy_%': 'N/A',
                '多拷贝基因数|Multi_Copy_BUSCOs': 'N/A',
                '片段化|Fragmented_%': 'N/A',
                '片段化基因数|Fragmented_BUSCOs': 'N/A',
                '缺失|Missing_%': 'N/A',
                '缺失基因数|Missing_BUSCOs': 'N/A',
                '标记基因总数|Total_Markers': 'N/A',
                '域|Domain': 'N/A',
                '状态|Status': '失败|Failure'
            }
            failed_rows.append(failed_row)

        # 添加到DataFrame|Add to DataFrame
        failed_df = pd.DataFrame(failed_rows)
        result_df = pd.concat([df, failed_df], ignore_index=True)

        return result_df
    
    def save_results(self, df: pd.DataFrame, output_filename: str = None):
        """保存结果到文件|Save results to file"""
        if df.empty:
            self.logger.warning("结果表格为空，跳过保存|Results table is empty, skipping save")
            return

        if not output_filename:
            output_filename = f"busco_summary_results.{self.config.output_format}"

        output_file = self.config.output_path / output_filename

        try:
            if self.config.output_format.lower() == 'csv':
                df.to_csv(output_file, index=False, encoding='utf-8')
            elif self.config.output_format.lower() == 'xlsx':
                df.to_excel(output_file, index=False, engine='openpyxl')
            else:  # 默认txt格式|Default txt format
                df.to_csv(output_file, index=False, sep='\t', encoding='utf-8')

            self.logger.info(f"结果已保存|Results saved: {output_file}")

            # 打印结果摘要|Print results summary
            self.print_results_summary(df)

        except Exception as e:
            self.logger.error(f"保存结果失败|Failed to save results: {e}")

    def print_results_summary(self, df: pd.DataFrame):
        """打印结果摘要|Print results summary"""
        self.logger.info("="*80)
        self.logger.info("BUSCO分析结果摘要|BUSCO Analysis Results Summary")
        self.logger.info("="*80)

        total_samples = len(df)
        successful_samples = len(df[df['状态|Status'].str.contains('成功|Success', na=False)])
        failed_samples = total_samples - successful_samples

        self.logger.info(f"总样本数|Total samples: {total_samples}")
        self.logger.info(f"成功样本数|Successful samples: {successful_samples}")
        self.logger.info(f"失败样本数|Failed samples: {failed_samples}")

        if successful_samples > 0:
            # 计算成功样本的统计信息|Calculate statistics for successful samples
            success_df = df[df['状态|Status'].str.contains('成功|Success', na=False)]

            # 转换百分比列为数值|Convert percentage columns to numeric
            numeric_cols = [
                '完整度|Complete_%',
                '单拷贝|Single_Copy_%',
                '多拷贝|Multi_Copy_%',
                '片段化|Fragmented_%',
                '缺失|Missing_%'
            ]

            for col in numeric_cols:
                success_df[col] = pd.to_numeric(success_df[col], errors='coerce')

            self.logger.info(f"平均完整度|Average completeness: {success_df['完整度|Complete_%'].mean():.2f}%")
            self.logger.info(f"平均单拷贝比例|Average single-copy ratio: {success_df['单拷贝|Single_Copy_%'].mean():.2f}%")
            self.logger.info(f"平均缺失比例|Average missing ratio: {success_df['缺失|Missing_%'].mean():.2f}%")

        self.logger.info("="*80)

class SummaryGenerator:
    """总结生成器|Summary Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary_report(self, df: pd.DataFrame):
        """生成总结报告|Generate summary report"""
        report_file = self.config.output_path / "busco_analysis_summary.txt"

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("BUSCO质量评估分析总结报告|BUSCO Quality Assessment Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")

            # 输入信息|Input information
            f.write("输入信息|Input Information:\n")
            f.write(f"  - 输入路径|Input path: {self.config.input_path}\n")
            f.write(f"  - 数据库|Database: {self.config.lineage}\n")
            f.write(f"  - 分析模式|Analysis mode: {self.config.mode}\n")
            f.write(f"  - 样本后缀模式|Sample suffix pattern: {self.config.sample_suffix}\n")
            f.write("\n")

            # 分析参数|Analysis parameters
            f.write("分析参数|Analysis Parameters:\n")
            f.write(f"  - CPU线程数|CPU threads: {self.config.threads}\n")
            f.write(f"  - E值阈值|E-value threshold: {self.config.evalue}\n")
            f.write(f"  - 候选区域限制|Candidate region limit: {self.config.limit}\n")
            f.write(f"  - 数据集版本|Dataset version: {self.config.datasets_version}\n")

            if self.config.augustus:
                f.write(f"  - 使用Augustus|Use Augustus: 是|Yes\n")
            if self.config.metaeuk:
                f.write(f"  - 使用Metaeuk|Use Metaeuk: 是|Yes\n")
            if self.config.miniprot:
                f.write(f"  - 使用Miniprot|Use Miniprot: 是|Yes\n")

            f.write("\n")

            # 结果统计|Results statistics
            if not df.empty:
                total_samples = len(df)
                successful_samples = len(df[df['状态|Status'].str.contains('成功|Success', na=False)])
                failed_samples = total_samples - successful_samples

                f.write("结果统计|Results Statistics:\n")
                f.write(f"  - 总样本数|Total samples: {total_samples}\n")
                f.write(f"  - 成功分析|Successful analyses: {successful_samples}\n")
                f.write(f"  - 失败分析|Failed analyses: {failed_samples}\n")
                f.write(f"  - 成功率|Success rate: {(successful_samples/total_samples*100):.2f}%\n")

            f.write(f"\n输出目录|Output directory: {self.config.output_dir}\n")

        self.logger.info(f"总结报告已生成|Summary report generated: {report_file}")
