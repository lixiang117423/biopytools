"""
HMMsearch分析主程序模块|HMMsearch Analysis Main Module
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path
from .config import HMMsearchConfig
from .utils import HMMsearchLogger
from .parser import DomtbloutParser
from .sequence_extractor import SequenceExtractor
from .hmmsearch_runner import HMMsearchRunner


class HMMsearchAnalyzer:
    """HMMsearch分析主类|HMMsearch Analyzer Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = HMMsearchConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = HMMsearchLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化解析器|Initialize parser
        self.parser = DomtbloutParser(self.logger)

        # 确定工作模式|Determine working mode
        self.run_hmmsearch_mode = self.config.hmm_file is not None

        # 如果是运行hmmsearch模式，初始化runner
        if self.run_hmmsearch_mode:
            self.hmmsearch_runner = HMMsearchRunner(self.logger, self.config)

        # 初始化序列提取器|Initialize sequence extractor
        if self.config.extract_protein_sequences or self.config.extract_domain_sequences:
            self.extractor = SequenceExtractor(self.logger, self.config.protein_fastas)
        else:
            self.extractor = None

    def run_analysis(self):
        """运行分析|Run analysis"""
        self.logger.info("=" * 60)
        if self.run_hmmsearch_mode:
            self.logger.info("开始HMMsearch完整分析流程|Starting complete HMMsearch pipeline")
        else:
            self.logger.info("开始HMMsearch结果分析|Starting HMMsearch result analysis")
        self.logger.info("=" * 60)

        # 0. 如果需要，运行hmmsearch|Run hmmsearch if needed
        if self.run_hmmsearch_mode:
            self.logger.info("步骤0: 运行HMMsearch|Step 0: Running HMMsearch")
            domtblout_file = self.hmmsearch_runner.run()
            self.config.domtblout_file = domtblout_file
        else:
            domtblout_file = self.config.domtblout_file

        # 1. 解析domtblout文件|Parse domtblout file
        self.logger.info("步骤1: 解析domtblout文件|Step 1: Parsing domtblout file")
        hits = self.parser.parse(self.config.domtblout_file)

        if not hits:
            self.logger.warning("未解析到任何domain命中|No domain hits found")
            return

        # 2. 过滤结果|Filter results
        if self.config.evalue_threshold is not None or self.config.score_threshold is not None:
            self.logger.info("步骤2: 过滤结果|Step 2: Filtering results")
            hits = self.parser.filter_hits(
                hits,
                evalue_threshold=self.config.evalue_threshold,
                score_threshold=self.config.score_threshold
            )

        # 3. 生成统计摘要|Generate summary
        self.logger.info("步骤3: 生成统计摘要|Step 3: Generating summary")
        summary = self.parser.generate_summary(hits)
        self._print_summary(summary)

        # 3.5 多物种模式下生成汇总表格|Generate species summary table in directory mode
        if len(self.config.protein_fastas) > 1:
            self._generate_species_summary()

        # 4. 转换为DataFrame并输出|Convert to DataFrame and output
        self.logger.info("步骤4: 生成表格|Step 4: Generating tables")
        df = self.parser.hits_to_dataframe(hits)

        # 多文件模式下添加来源文件名列|Add source file column in directory mode
        if self.extractor and len(self.config.protein_fastas) > 1:
            df.insert(0, 'source_file', df['target_name'].map(self.extractor.sequence_sources))

        # 输出CSV|Output CSV
        if self.config.output_csv:
            csv_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}.csv")
            df.to_csv(csv_file, index=False)
            self.logger.info(f"CSV文件已保存|CSV file saved: {csv_file}")

        # 输出Excel|Output Excel
        if self.config.output_excel:
            excel_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}.xlsx")
            df.to_excel(excel_file, index=False, engine='openpyxl')
            self.logger.info(f"Excel文件已保存|Excel file saved: {excel_file}")

        # 5. 提取序列|Extract sequences
        if self.extractor:
            if self.config.extract_protein_sequences:
                self.logger.info("步骤5: 提取蛋白序列|Step 5: Extracting protein sequences")
                protein_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_proteins.fa")
                self.extractor.extract_protein_sequences(hits, protein_file)

            if self.config.extract_domain_sequences:
                self.logger.info("步骤6: 提取domain序列|Step 6: Extracting domain sequences")
                domain_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_domains.fa")
                self.extractor.extract_domain_sequences(hits, domain_file)

        self.logger.info("=" * 60)
        self.logger.info("分析完成|Analysis completed")
        self.logger.info("=" * 60)

    def _generate_species_summary(self):
        """生成物种汇总表格（仅目录模式）|Generate species summary table (directory mode only)"""
        self.logger.info("生成物种汇总表格|Generating species summary table")

        species_stats = []
        total_hits = 0
        total_targets = 0

        for fasta_path in self.config.protein_fastas:
            species_name = Path(fasta_path).stem
            per_file_domtblout = os.path.join(
                self.config.output_dir,
                f"{self.config.output_prefix}_{species_name}.domtblout"
            )

            if not os.path.exists(per_file_domtblout):
                species_stats.append({
                    '物种|Species': species_name,
                    '总命中数|Total Hits': 0,
                    '唯一基因数|Unique Targets': 0,
                })
                continue

            hits = self.parser.parse(per_file_domtblout)
            unique_targets = len(set(h.target_name for h in hits))
            species_stats.append({
                '物种|Species': species_name,
                '总命中数|Total Hits': len(hits),
                '唯一基因数|Unique Targets': unique_targets,
            })
            total_hits += len(hits)
            total_targets += unique_targets

        # 添加合计行|Add totals row
        species_stats.append({
            '物种|Species': '合计|Total',
            '总命中数|Total Hits': total_hits,
            '唯一基因数|Unique Targets': total_targets,
        })

        df = pd.DataFrame(species_stats)

        summary_csv = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_species_summary.csv")
        summary_excel = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_species_summary.xlsx")

        if self.config.output_csv:
            df.to_csv(summary_csv, index=False)
            self.logger.info(f"物种汇总CSV已保存|Species summary CSV saved: {summary_csv}")

        if self.config.output_excel:
            df.to_excel(summary_excel, index=False, engine='openpyxl')
            self.logger.info(f"物种汇总Excel已保存|Species summary Excel saved: {summary_excel}")

        self.logger.info(f"共{len(species_stats)-1}个物种，总命中{total_hits}，唯一基因{total_targets}|{len(species_stats)-1} species, {total_hits} hits, {total_targets} unique targets")

    def _print_summary(self, summary: dict):
        """打印统计摘要|Print summary"""
        self.logger.info("-" * 60)
        self.logger.info("统计摘要|Statistical Summary")
        self.logger.info("-" * 60)
        self.logger.info(f"总命中数|Total hits: {summary['total_hits']}")
        self.logger.info(f"唯一基因数|Unique targets: {summary['unique_targets']}")

        self.logger.info("-" * 60)
        self.logger.info("Domain数量分布|Domain Count Distribution")
        self.logger.info("-" * 60)
        for count, genes in sorted(summary['domain_count_distribution'].items()):
            self.logger.info(f"{count} domain(s): {genes} gene(s)")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='HMMsearch分析工具：运行hmmsearch并处理结果|HMMsearch Analysis Tool: Run hmmsearch and process results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 主输入参数（自动判断模式）|Main input parameter (auto-detect mode)
    parser.add_argument('-i', '--input',
                       help='输入文件：domtblout文件（模式1）或HMM文件（模式2，需同时指定-p）|Input file: domtblout (mode 1) or HMM file (mode 2, requires -p)')

    parser.add_argument('-p', '--protein-fasta',
                       help='蛋白序列FASTA文件或目录（模式2必需，模式1提取序列时需要）|Protein FASTA file or directory (required for mode 2, needed for mode 1 if extracting sequences)')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir', default='./hmmsearch_output',
                       help='输出目录|Output directory')
    parser.add_argument('--output-prefix', default='hmmsearch_results',
                       help='输出文件前缀|Output file prefix')

    # HMMsearch软件配置|HMMsearch software configuration
    parser.add_argument('--hmmsearch-path',
                       default='~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch',
                       help='hmmsearch程序路径|hmmsearch program path')
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')

    # HMMsearch运行参数|HMMsearch run parameters
    parser.add_argument('--evalue-cutoff', type=float,
                       help='hmmsearch报告E-value阈值|hmmsearch reporting E-value threshold')
    parser.add_argument('--score-cutoff', type=float,
                       help='hmmsearch报告分数阈值|hmmsearch reporting score threshold')
    parser.add_argument('--cut-tc', action='store_true',
                       help='使用模型的TC trusted cutoff|Use model TC trusted cutoff')
    parser.add_argument('--cut-ga', action='store_true',
                       help='使用模型的GA gathering cutoff|Use model GA gathering cutoff')
    parser.add_argument('--cut-nc', action='store_true',
                       help='使用模型的NC noise cutoff|Use model NC noise cutoff')

    # 结果过滤参数|Result filtering parameters (post-processing)
    parser.add_argument('-e', '--evalue-threshold', type=float,
                       help='Domain E-value阈值(保留小于等于该值的)|Domain E-value threshold (keep <= this value)')
    parser.add_argument('-s', '--score-threshold', type=float,
                       help='Domain分数阈值(保留大于等于该值的)|Domain score threshold (keep >= this value)')

    # 序列提取选项|Sequence extraction options
    parser.add_argument('--extract-proteins', action='store_true', default=True,
                       help='提取匹配的蛋白序列|Extract matched protein sequences (default: True)')
    parser.add_argument('--no-extract-proteins', action='store_false', dest='extract_proteins',
                       help='不提取蛋白序列|Do not extract protein sequences')
    parser.add_argument('--extract-domains', action='store_true', default=True,
                       help='提取domain序列|Extract domain sequences (default: True)')
    parser.add_argument('--no-extract-domains', action='store_false', dest='extract_domains',
                       help='不提取domain序列|Do not extract domain sequences')

    # 输出格式选项|Output format options
    parser.add_argument('--no-csv', action='store_true',
                       help='不输出CSV文件|Do not output CSV file')
    parser.add_argument('--no-excel', action='store_true',
                       help='不输出Excel文件|Do not output Excel file')

    args = parser.parse_args()

    # 智能判断工作模式|Intelligently determine working mode
    domtblout = None
    hmm_file = None

    if args.protein_fasta:
        # 如果指定了-p，说明是模式2（运行hmmsearch）|If -p specified, mode 2 (run hmmsearch)
        hmm_file = args.input
    else:
        # 否则是模式1（处理domtblout）|Otherwise mode 1 (process domtblout)
        domtblout = args.input

    # 创建分析器并运行|Create analyzer and run
    analyzer = HMMsearchAnalyzer(
        domtblout_file=domtblout,
        hmm_file=hmm_file,
        protein_fasta=args.protein_fasta,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        hmmsearch_path=args.hmmsearch_path,
        threads=args.threads,
        evalue_cutoff=args.evalue_cutoff,
        score_cutoff=args.score_cutoff,
        use_cut_tc=args.cut_tc,
        use_cut_ga=args.cut_ga,
        use_cut_nc=args.cut_nc,
        evalue_threshold=args.evalue_threshold,
        score_threshold=args.score_threshold,
        extract_protein_sequences=args.extract_proteins,
        extract_domain_sequences=args.extract_domains,
        output_csv=not args.no_csv,
        output_excel=not args.no_excel
    )

    analyzer.run_analysis()


if __name__ == "__main__":
    main()
