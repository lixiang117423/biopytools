"""
GEC阈值计算器|GEC Threshold Calculator
功能: 解析GEC输出结果并计算显著性阈值|
Features: Parse GEC output and calculate significance thresholds
"""

import gzip
import os
from pathlib import Path
from typing import Dict, List, Tuple


class GECThresholdCalculator:
    """GEC阈值计算器类|GEC Threshold Calculator Class"""

    def __init__(self, logger):
        self.logger = logger

    def parse_effective_size_file(self, result_file: str) -> List[Dict]:
        """
        解析effective size文件|Parse effective size file

        Args:
            result_file: GEC输出文件路径 (*.effective.size.txt.gz)|GEC output file path

        Returns:
            包含每个LD块信息的字典列表|List of dictionaries containing LD block information
        """
        self.logger.info(f"解析有效检验数文件|Parsing effective test number file: {result_file}")

        blocks = []

        try:
            # 判断是否为gzip文件|Check if file is gzipped
            open_func = gzip.open if result_file.endswith('.gz') else open

            with open_func(result_file, 'rt') as f:
                # 读取并解析表头|Read and parse header
                header = f.readline().strip()
                self.logger.debug(f"文件表头|File header: {header}")

                # 解析列索引|Parse column indices
                columns = header.split('\t')
                col_indices = {col: i for i, col in enumerate(columns)}

                # 读取数据行|Read data rows
                for line_num, line in enumerate(f, 2):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    fields = line.split('\t')

                    block_info = {
                        'chrom': fields[col_indices.get('Chrom', 0)],
                        'start_pos': int(fields[col_indices.get('StartPos', 1)]),
                        'end_pos': int(fields[col_indices.get('EndPos', 2)]),
                        'num': int(fields[col_indices.get('Num', 3)]),
                        'effective_num': float(fields[col_indices.get('EffectiveNum', 4)])
                    }

                    blocks.append(block_info)

            self.logger.info(f"成功解析 {len(blocks)} 个LD块|Successfully parsed {len(blocks)} LD blocks")
            return blocks

        except FileNotFoundError:
            self.logger.error(f"文件不存在|File not found: {result_file}")
            raise
        except Exception as e:
            self.logger.error(f"解析文件失败|Failed to parse file: {str(e)}")
            raise

    def calculate_total_effective_tests(self, blocks: List[Dict]) -> float:
        """
        计算总有效检验数|Calculate total effective test number

        Args:
            blocks: LD块信息列表|List of LD block information

        Returns:
            总有效检验数|Total effective test number
        """
        total_effective = sum(block['effective_num'] for block in blocks)

        self.logger.info(f"总有效检验数|Total effective test number: {total_effective:.2f}")

        return total_effective

    def calculate_significance_threshold(self, total_effective: float, alpha: float = 0.05) -> float:
        """
        计算显著性阈值|Calculate significance threshold

        Args:
            total_effective: 总有效检验数|Total effective test number
            alpha: 显著性水平|Significance level

        Returns:
            校正后的P值阈值|Adjusted P-value threshold
        """
        threshold = alpha / total_effective

        self.logger.info(f"显著性阈值计算|Significance threshold calculation:")
        self.logger.info(f"  Alpha水平|Alpha level: {alpha}")
        self.logger.info(f"  总有效检验数|Total effective tests: {total_effective:.2f}")
        self.logger.info(f"  校正后阈值|Adjusted threshold: {threshold:.2e}")

        return threshold

    def generate_summary_report(self, blocks: List[Dict], total_effective: float,
                                threshold: float, output_file: str):
        """
        生成汇总报告|Generate summary report

        Args:
            blocks: LD块信息列表|List of LD block information
            total_effective: 总有效检验数|Total effective test number
            threshold: 校正后的阈值|Adjusted threshold
            output_file: 输出文件路径|Output file path
        """
        self.logger.info(f"生成汇总报告|Generating summary report: {output_file}")

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("GEC基因组范围多重检验校正报告|GEC Genome-wide Error Correction Report\n")
            f.write("=" * 80 + "\n\n")

            # 基本信息|Basic information
            f.write("一、校正结果汇总|I. Correction Summary\n")
            f.write("-" * 80 + "\n")
            f.write(f"总LD块数|Total LD blocks: {len(blocks)}\n")
            f.write(f"总有效检验数|Total effective tests: {total_effective:.2f}\n")
            f.write(f"显著性水平|Significance level (alpha): 0.05\n")
            f.write(f"校正后阈值|Adjusted threshold: {threshold:.2e}\n")
            f.write("\n")

            # 各染色体统计|Chromosome-wise statistics
            f.write("二、各染色体统计|II. Chromosome-wise Statistics\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'染色体|Chrom':<10} {'LD块数|Blocks':<15} {'有效检验数|EffTests':<20}\n")
            f.write("-" * 80 + "\n")

            chrom_stats = {}
            for block in blocks:
                chrom = block['chrom']
                if chrom not in chrom_stats:
                    chrom_stats[chrom] = {'count': 0, 'effective': 0}
                chrom_stats[chrom]['count'] += 1
                chrom_stats[chrom]['effective'] += block['effective_num']

            for chrom in sorted(chrom_stats.keys(), key=lambda x: (len(x), x)):
                stats = chrom_stats[chrom]
                f.write(f"{chrom:<10} {stats['count']:<15} {stats['effective']:<20.2f}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("说明|Notes:\n")
            f.write("1. 有效检验数考虑了LD结构，显著小于原始SNP数量|"
                   "Effective tests account for LD structure, significantly less than raw SNP count\n")
            f.write("2. 校正阈值 = Alpha / 总有效检验数|Adjusted threshold = Alpha / Total effective tests\n")
            f.write("3. 建议使用此阈值作为GWAS显著性判断标准|"
                   "Recommend using this threshold for GWAS significance determination\n")
            f.write("=" * 80 + "\n")

        self.logger.info(f"汇总报告已保存|Summary report saved: {output_file}")

    def run_analysis(self, result_file: str, alpha: float = 0.05,
                     output_prefix: str = None) -> Dict:
        """
        运行完整的分析流程|Run complete analysis pipeline

        Args:
            result_file: GEC输出文件路径|GEC output file path
            alpha: 显著性水平|Significance level
            output_prefix: 输出文件前缀|Output file prefix

        Returns:
            包含分析结果的字典|Dictionary containing analysis results
        """
        # 解析有效检验数文件|Parse effective test number file
        blocks = self.parse_effective_size_file(result_file)

        # 计算总有效检验数|Calculate total effective test number
        total_effective = self.calculate_total_effective_tests(blocks)

        # 计算显著性阈值|Calculate significance threshold
        threshold = self.calculate_significance_threshold(total_effective, alpha)

        # 生成汇总报告|Generate summary report
        if output_prefix:
            summary_file = f"{output_prefix}_summary.txt"
            self.generate_summary_report(blocks, total_effective, threshold, summary_file)

        return {
            'total_ld_blocks': len(blocks),
            'total_effective_tests': total_effective,
            'significance_threshold': threshold,
            'alpha': alpha,
            'blocks': blocks
        }
