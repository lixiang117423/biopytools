"""
甲基化分析主程序模块 | Methylation Analysis Main Module
"""

import argparse
import os
import sys

from .config import MethylationConfig
from .core_analysis import BasicPipeline, EnhancedAnalysis
from .reports import ReportGenerator
from .utils import (
    CommandRunner,
    MethylationLogger,
    check_completed_steps,
    check_dependencies,
)


class MethylationAnalyzer:
    """甲基化分析主类 | Main Methylation Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = MethylationConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = MethylationLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

        # 初始化处理器 | Initialize processors
        self.basic_pipeline = BasicPipeline(self.config, self.logger, self.cmd_runner)

        # 增强分析处理器 | Enhanced analysis processors
        if self.config.enhanced_mode:
            self.enhanced_analysis = EnhancedAnalysis(
                self.config, self.logger, self.cmd_runner
            )

        # 报告生成器 | Report generator
        self.report_generator = ReportGenerator(self.config, self.logger)

    def check_dependencies(self):
        """🔍 检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def run_analysis(self):
        """🚀 运行完整的甲基化分析流程 | Run complete methylation analysis pipeline"""
        try:
            self.logger.info(
                "🚀 开始甲基化分析流程 | Starting methylation analysis pipeline"
            )
            self.logger.info(
                f"🔧 分析模式 | Analysis mode: {'Enhanced' if self.config.enhanced_mode else 'Basic'}"
            )

            # 检查依赖
            self.check_dependencies()

            # 样品分组信息
            self.logger.info("🧬 样品分组信息 | Sample grouping information:")
            all_groups = self.config.all_groups
            for group in all_groups:
                group_samples = self.config.get_samples_by_group(group)
                if group == "CK":
                    emoji = "✅"
                elif group == "Treatment":
                    emoji = "🦠"
                else:
                    emoji = "🔬"
                self.logger.info(
                    f"  {emoji} {group}组: {group_samples} ({len(group_samples)} 样品)"
                )

            # 显示分组比较计划
            if self.config.enhanced_mode:
                comparisons = self.config.get_pairwise_comparisons()
                if comparisons:
                    self.logger.info(f"🔬 计划进行 {len(comparisons)} 个分组比较:")
                    for i, (group1, group2) in enumerate(comparisons, 1):
                        self.logger.info(f"  {i}. {group1} vs {group2}")
                else:
                    self.logger.warning(
                        "⚠️  没有足够的分组进行比较分析（需要至少2个分组）"
                    )

            # 智能检测已完成的步骤
            status = check_completed_steps(self.config, self.logger)

            # 执行基础流程
            self._run_basic_pipeline(status)

            # 执行增强分析
            if self.config.enhanced_mode:
                self._run_enhanced_analysis()

            # 生成综合报告
            if not self.report_generator.generate_comprehensive_report():
                self.logger.warning("⚠️  报告生成失败，但分析已完成")

            self.logger.info("🎉 === 甲基化分析流程完成 ===")
            self._print_completion_summary()

        except Exception as e:
            self.logger.error(f"❌ 分析流程失败 | Analysis pipeline failed: {e}")
            raise

    def _run_basic_pipeline(self, status):
        """🏗️ 运行基础分析流程 | Run basic analysis pipeline"""
        self.logger.info("🏗️ 执行基础甲基化分析流程...")

        # 步骤1: 文件重命名和清理 (智能跳过)
        if status["need_rename"]:
            if not self.basic_pipeline.rename_files():
                raise RuntimeError("❌ 文件重命名失败")
        else:
            self.logger.info("📝 步骤1: 🚀 跳过文件重命名 - 已完成或无需执行")

        # 步骤2: fastp质控 (智能跳过)
        if status["fastp_completed"]:
            self.logger.info("🧹 步骤2: 🚀 跳过fastp质控 - 已完成")
        else:
            if not self.basic_pipeline.run_fastp():
                raise RuntimeError("❌ fastp质控失败")

        # 步骤3: 构建Bismark索引 (智能跳过)
        if status["index_exists"]:
            self.logger.info("🏗️ 步骤3: 🚀 跳过Bismark索引构建 - 已存在")
        else:
            if not self.basic_pipeline.build_bismark_index():
                raise RuntimeError("❌ Bismark索引构建失败")

        # 步骤4: 甲基化比对分析 (智能续传)
        if not self.basic_pipeline.run_methylation_mapping():
            raise RuntimeError("❌ 甲基化比对分析失败")

    def _run_enhanced_analysis(self):
        """🔬 运行增强分析流程 | Run enhanced analysis pipeline"""
        self.logger.info("🔬 ✅ 基础甲基化分析已完成，开始增强分析...")

        # 步骤A: methylKit差异甲基化分析
        if not self.enhanced_analysis.run_methylkit_analysis():
            self.logger.warning("⚠️  methylKit分析失败，但继续执行其他分析")

        # 步骤C: 基因组bins分析
        if not self.enhanced_analysis.run_genomic_bins_analysis():
            self.logger.warning("⚠️  基因组bins分析失败，但继续执行其他分析")

        # 步骤F: BLAST定位和目标序列分析
        if not self.enhanced_analysis.run_blast_analysis():
            self.logger.warning("⚠️  BLAST分析失败，但继续执行其他分析")

    def _print_completion_summary(self):
        """🎯 打印完成总结 | Print completion summary"""
        self.logger.info("")
        self.logger.info(
            "🎉 分析完成！主要结果文件 | Analysis completed! Main result files:"
        )

        # 基础结果
        self.logger.info(f"  📊 Clean数据 | Clean data: {self.config.clean_dir}")
        self.logger.info(f"  📈 比对结果 | Mapping results: {self.config.mapping_dir}")

        # 增强分析结果
        if self.config.enhanced_mode and hasattr(self.config, "analysis_dir"):
            self.logger.info(
                f"  🧬 差异甲基化 | Differential methylation: {os.path.join(self.config.analysis_dir, 'differential_methylation')}"
            )

            # 显示分组比较结果统计
            comparisons = self.config.get_pairwise_comparisons()
            if comparisons:
                successful_comparisons = 0
                for group1, group2 in comparisons:
                    comparison_name = f"{group1}_vs_{group2}"
                    diff_sites_file = os.path.join(
                        self.config.analysis_dir,
                        "differential_methylation",
                        f"{comparison_name}_differential_sites.txt",
                    )
                    if os.path.exists(diff_sites_file):
                        successful_comparisons += 1

                self.logger.info(
                    f"  📊 分组比较结果: {successful_comparisons}/{len(comparisons)} 个比较成功完成"
                )
                if successful_comparisons > 0:
                    self.logger.info(
                        f"     查看具体结果: ls -la {os.path.join(self.config.analysis_dir, 'differential_methylation')}/*_vs_*"
                    )

            self.logger.info(
                f"  🎯 目标序列分析 | Target sequence analysis: {os.path.join(self.config.analysis_dir, 'target_sequence_analysis')}"
            )
            self.logger.info(
                f"  🔍 BLAST定位 | BLAST localization: {os.path.join(self.config.analysis_dir, 'blast_results')}"
            )

        # 报告
        reports_dir = os.path.join(self.config.output_dir, "reports")
        self.logger.info(
            f"  📋 综合报告 | Comprehensive report: {os.path.join(reports_dir, 'comprehensive_methylation_report.txt')}"
        )

        self.logger.info("")
        self.logger.info("🔍 快速检查命令 | Quick check commands:")
        bismark_results_dir = os.path.join(self.config.mapping_dir, "bismark_results")
        self.logger.info(
            f"  查看比对统计 | Check mapping stats: ls -la {bismark_results_dir}/*/CX_report.txt"
        )

        if self.config.enhanced_mode and hasattr(self.config, "analysis_dir"):
            # 显示批量比较结果的检查命令
            comparisons = self.config.get_pairwise_comparisons()
            if comparisons and len(comparisons) > 0:
                self.logger.info(
                    f"  查看分组比较结果 | Check group comparison results:"
                )
                for group1, group2 in comparisons[:3]:  # 只显示前3个作为示例
                    comparison_name = f"{group1}_vs_{group2}"
                    self.logger.info(
                        f"    {comparison_name}: head -10 {os.path.join(self.config.analysis_dir, 'differential_methylation', f'{comparison_name}_differential_sites.txt')}"
                    )
                if len(comparisons) > 3:
                    self.logger.info(f"    ... 还有 {len(comparisons) - 3} 个比较结果")

                # 显示综合报告检查命令
                self.logger.info(
                    f"  查看综合分析报告 | Check comprehensive analysis report: cat {os.path.join(self.config.analysis_dir, 'reports', 'methylkit_comprehensive_summary.txt')}"
                )
            else:
                self.logger.info(
                    f"  查看差异甲基化 | Check differential methylation: head -10 {os.path.join(self.config.analysis_dir, 'differential_methylation', 'all_differential_sites.txt')}"
                )

            if self.config.target_promoter_fa:
                self.logger.info(
                    f"  查看目标序列定位 | Check target localization: cat {os.path.join(self.config.analysis_dir, 'blast_results', 'target_positions.bed')}"
                )


def main():
    """🎯 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🧬 甲基化分析工具包 | Methylation Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基础分析 | Basic analysis
  python methylation_pipeline.py -r /path/to/raw -g /path/to/genome.fa -o /path/to/output
  
  # 完整增强分析 | Complete enhanced analysis
  python methylation_pipeline.py -r /path/to/raw -g /path/to/genome.fa -o /path/to/output \\
      -p /path/to/target_promoter.fa -a /path/to/genome.gff3 -e
      
  # 自定义参数 | Custom parameters
  python methylation_pipeline.py -r ./raw -g ./genome.fa -o ./results \\
      -j 16 -c 10 -d 0.25 -v 0.001 -e
        """,
    )

    # 必需参数 | Required arguments
    parser.add_argument(
        "-r", "--raw-dir", required=True, help="📁 原始数据目录 | Raw data directory"
    )
    parser.add_argument(
        "-g", "--genome", required=True, help="🧬 基因组FASTA文件 | Genome FASTA file"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="./methylation_output",
        help="📂 输出目录 | Output directory (default: ./methylation_output)",
    )

    # 可选输入文件 | Optional input files
    parser.add_argument(
        "-t", "--target-fa", help="🎯 目标序列文件 | Target sequence file"
    )
    parser.add_argument(
        "-p",
        "--target-promoter-fa",
        help="🎯 目标启动子序列文件 | Target promoter sequence file",
    )
    parser.add_argument(
        "-a", "--annotation-gff", help="📊 基因注释GFF文件 | Gene annotation GFF file"
    )

    # 分析模式 | Analysis mode
    parser.add_argument(
        "-e",
        "--enhanced-mode",
        action="store_true",
        help="🔬 启用增强分析模式 | Enable enhanced analysis mode (default: False)",
    )

    # 基础分析参数 | Basic analysis parameters
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=88,
        help="🔧 线程数 | Number of threads (default: 88)",
    )
    parser.add_argument(
        "-c",
        "--min-coverage",
        type=int,
        default=5,
        help="📈 最小覆盖度 | Minimum coverage (default: 5)",
    )
    parser.add_argument(
        "-n",
        "--min-cytosines",
        type=int,
        default=5,
        help="🧬 最小胞嘧啶数 | Minimum cytosines (default: 5)",
    )

    # 差异分析参数 | Differential analysis parameters
    parser.add_argument(
        "-d",
        "--methylation-diff-threshold",
        type=float,
        default=0.3,
        help="📊 甲基化差异阈值 | Methylation difference threshold (default: 0.3)",
    )
    parser.add_argument(
        "-v",
        "--pvalue-threshold",
        type=float,
        default=0.01,
        help="📈 p值阈值 | P-value threshold (default: 0.01)",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=100,
        help="🪟 窗口大小(bp) | Window size in bp (default: 100)",
    )
    parser.add_argument(
        "-s",
        "--step-size",
        type=int,
        default=50,
        help="👣 步长(bp) | Step size in bp (default: 50)",
    )
    parser.add_argument(
        "-m",
        "--merge-distance",
        type=int,
        default=200,
        help="🔗 差异区域合并距离(bp) | Merge distance for differential regions in bp (default: 200)",
    )

    # bins分析参数 | Bins analysis parameters
    parser.add_argument(
        "--gene-bins",
        type=int,
        default=60,
        help="📦 基因体bins数 | Number of gene body bins (default: 60)",
    )
    parser.add_argument(
        "--flank-bins",
        type=int,
        default=20,
        help="📦 侧翼区域bins数 | Number of flanking region bins (default: 20)",
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=2000,
        help="📏 侧翼区域大小(bp) | Flanking region size in bp (default: 2000)",
    )

    # 工具路径 | Tool paths (基础工具)
    parser.add_argument(
        "--fastp-path",
        default="fastp",
        help="🧹 fastp可执行文件路径 | fastp executable path (default: fastp)",
    )
    parser.add_argument(
        "--bismark-path",
        default="bismark",
        help="🎯 bismark可执行文件路径 | bismark executable path (default: bismark)",
    )
    parser.add_argument(
        "--bismark-genome-preparation-path",
        default="bismark_genome_preparation",
        help="🏗️ bismark_genome_preparation路径 | bismark_genome_preparation path (default: bismark_genome_preparation)",
    )
    parser.add_argument(
        "--bowtie2-path",
        default="bowtie2",
        help="🎯 bowtie2可执行文件路径 | bowtie2 executable path (default: bowtie2)",
    )
    parser.add_argument(
        "--deduplicate-bismark-path",
        default="deduplicate_bismark",
        help="🔄 deduplicate_bismark路径 | deduplicate_bismark path (default: deduplicate_bismark)",
    )
    parser.add_argument(
        "--bismark-methylation-extractor-path",
        default="bismark_methylation_extractor",
        help="🧪 bismark_methylation_extractor路径 | bismark_methylation_extractor path (default: bismark_methylation_extractor)",
    )
    parser.add_argument(
        "--bismark2report-path",
        default="bismark2report",
        help="📊 bismark2report路径 | bismark2report path (default: bismark2report)",
    )
    parser.add_argument(
        "--bismark2summary-path",
        default="bismark2summary",
        help="📋 bismark2summary路径 | bismark2summary path (default: bismark2summary)",
    )

    # 工具路径 | Tool paths (可选工具)
    parser.add_argument(
        "--makeblastdb-path",
        default="makeblastdb",
        help="🔍 makeblastdb可执行文件路径 | makeblastdb executable path (default: makeblastdb)",
    )
    parser.add_argument(
        "--blastn-path",
        default="blastn",
        help="🔍 blastn可执行文件路径 | blastn executable path (default: blastn)",
    )
    parser.add_argument(
        "--multiqc-path",
        default="multiqc",
        help="📊 multiqc可执行文件路径 | multiqc executable path (default: multiqc)",
    )
    parser.add_argument(
        "--bedtools-path",
        default="bedtools",
        help="🛠️ bedtools可执行文件路径 | bedtools executable path (default: bedtools)",
    )
    parser.add_argument(
        "--samtools-path",
        default="samtools",
        help="🛠️ samtools可执行文件路径 | samtools executable path (default: samtools)",
    )

    # R环境路径 | R environment path
    parser.add_argument(
        "--r-executable",
        default="/share/org/YZWL/yzwl_lixg/miniforge3/envs/methylkit/bin/R",
        help="📊 R可执行文件路径 | R executable path (default: /share/org/YZWL/yzwl_lixg/miniforge3/envs/methylkit/bin/R)",
    )

    # 样品分组参数 | Sample grouping parameters (高级用户)
    parser.add_argument(
        "--sample-groups-file",
        help="🧬 样品分组配置文件(制表符分隔的TXT文件或JSON格式) | Sample grouping config file (tab-separated TXT file or JSON format)\n"
        + "📋 功能特性 | Features:\n"
        + "  - 支持多个分组的两两比较分析 | Support pairwise comparison analysis for multiple groups\n"
        + "  - 自动避免重复比较 (A vs B = B vs A) | Automatically avoid duplicate comparisons\n"
        + "  - 批量生成差异甲基化结果 | Batch generation of differential methylation results\n"
        + "📁 文件格式 | File format:\n"
        + "  TXT格式 (推荐): 第一列样品名，第二列分组名，制表符分隔\n"
        + "  TXT format (recommended): column1=sample_name, column2=group_name, tab-separated\n"
        + "  🔬 示例 (3个分组) | Example (3 groups):\n"
        + "    FZY4201    Control\n"
        + "    FZY4202    Treatment_A\n"
        + "    FZY4203    Treatment_B\n"
        + "    FZY4204    Control\n"
        + "  ➡️  将自动进行: Control vs Treatment_A, Control vs Treatment_B, Treatment_A vs Treatment_B\n"
        + '  JSON格式: {"sample1": "group1", "sample2": "group2"}\n',
    )

    # 帮助和调试参数 | Help and debugging parameters
    parser.add_argument(
        "--help-all",
        action="store_true",
        help="❓ 显示所有参数的详细说明 | Show detailed help for all parameters",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="📝 详细输出模式 | Verbose output mode (default: False)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="🔍 试运行模式(仅检查配置不执行分析) | Dry run mode (check config only, default: False)",
    )

    args = parser.parse_args()

    # 处理详细帮助 | Handle detailed help
    if args.help_all:
        print("""
🧬 甲基化分析工具包 - 详细参数说明
===========================================

📁 输入文件参数 | Input File Parameters:
  -r, --raw-dir                原始FASTQ数据目录，应包含*_1.fq.gz和*_2.fq.gz文件
  -g, --genome                 参考基因组FASTA文件，用于构建bismark索引
  -t, --target-fa              目标序列FASTA文件（可选）
  -p, --target-promoter-fa     未知位置的目标启动子序列文件，需要BLAST定位
  -a, --annotation-gff         基因注释GFF3文件，用于基因组区域注释

🔧 分析模式 | Analysis Mode:
  -e, --enhanced-mode          启用增强分析模式，包含methylKit差异分析、BLAST定位等

⚙️  基础参数 | Basic Parameters:
  -o, --output-dir             输出目录，所有结果文件都保存在此目录下
  -j, --threads                并行线程数，建议设置为CPU核心数的一半
  -c, --min-coverage           最小覆盖度，低于此值的甲基化位点将被过滤
  -n, --min-cytosines          每个窗口最小胞嘧啶数

📊 差异分析参数 | Differential Analysis Parameters:
  -d, --methylation-diff-threshold  甲基化差异阈值(0-1)，如0.3表示30%差异
  -v, --pvalue-threshold       统计显著性p值阈值
  -w, --window-size            滑动窗口大小(bp)，用于窗口分析
  -s, --step-size              滑动窗口步长(bp)
  -m, --merge-distance         差异区域合并距离(bp)，距离小于此值的区域将合并

🧬 Bins分析参数 | Bins Analysis Parameters:
  --gene-bins                  基因体划分的bins数量，用于基因体甲基化分析
  --flank-bins                 侧翼区域bins数量
  --flank-size                 基因上下游侧翼区域大小(bp)

🛠️  工具路径 | Tool Paths:
  --fastp-path                 fastp质控工具路径
  --bismark-path               bismark比对工具路径
  --bismark-genome-preparation-path  bismark索引构建工具路径
  --bowtie2-path               bowtie2比对工具路径
  --deduplicate-bismark-path   bismark去重工具路径
  --bismark-methylation-extractor-path  甲基化提取工具路径
  --bismark2report-path        bismark报告生成工具路径
  --bismark2summary-path       bismark汇总工具路径
  --makeblastdb-path           BLAST数据库构建工具路径
  --blastn-path                BLAST搜索工具路径
  --multiqc-path               MultiQC质控报告工具路径
  --bedtools-path              Bedtools基因组区间工具路径
  --samtools-path              Samtools序列处理工具路径
  --r-executable               R语言可执行文件路径，用于methylKit分析

🔍 调试参数 | Debug Parameters:
  --sample-groups-file         自定义样品分组配置文件
                               🔬 批量分组比较功能:
                               - 自动识别所有分组并进行两两比较
                               - 避免重复比较 (A vs B = B vs A)
                               - 为每个比较生成独立的结果文件
                               支持两种格式：
                               1. TXT格式 (推荐): 制表符分隔，第一列样品名，第二列分组名
                                  示例多分组文件:
                                  FZY4201    Control
                                  FZY4202    Treatment_A  
                                  FZY4203    Treatment_B
                                  FZY4204    Control
                                  ➡️ 将进行3个比较: Control vs Treatment_A, Control vs Treatment_B, Treatment_A vs Treatment_B
                               2. JSON格式: {"sample1": "group1", "sample2": "group2"}
  --verbose                    详细输出模式，显示更多调试信息
  --dry-run                    试运行模式，仅验证配置不执行实际分析
  --help-all                   显示此详细帮助信息

📋 使用建议 | Usage Recommendations:
  1. 基础分析：只需指定-r, -g, -o三个必需参数
  2. 增强分析：添加-e参数并提供-p和-a文件以获得完整功能
  3. 高覆盖度数据：建议增加-c参数值（如-c 10）
  4. 严格分析：降低-d和-v参数值以获得更严格的差异标准
  5. 性能优化：根据机器配置调整-j参数

💡 示例命令 | Example Commands:
  # 快速开始
  python run_methylation_pipeline.py -r ./raw -g ./genome.fa -o ./results
  
  # 完整分析
  python run_methylation_pipeline.py -r ./raw -g ./genome.fa -o ./results \\
      -p ./target.fa -a ./genes.gff3 -e -j 16 -c 8 -d 0.25 -v 0.005
""")
        return

    # 处理试运行模式 | Handle dry run mode
    if args.dry_run:
        print(
            "🔍 试运行模式 - 仅验证配置 | Dry run mode - validating configuration only"
        )
        print("=" * 60)

    # 处理样品分组文件 | Handle sample groups file
    sample_groups = None
    if args.sample_groups_file:
        try:
            # 检查文件扩展名来决定解析方式
            file_ext = os.path.splitext(args.sample_groups_file)[1].lower()

            if file_ext == ".json":
                # JSON格式 | JSON format
                import json

                with open(args.sample_groups_file, "r", encoding="utf-8") as f:
                    sample_groups = json.load(f)
                print(f"✅ 已加载JSON格式样品分组配置: {args.sample_groups_file}")

            else:
                # 文本格式 (TXT, TSV等) | Text format (TXT, TSV, etc.)
                sample_groups = {}
                with open(args.sample_groups_file, "r", encoding="utf-8") as f:
                    for line_num, line in enumerate(f, 1):
                        line = line.strip()

                        # 跳过空行和注释行
                        if not line or line.startswith("#"):
                            continue

                        # 尝试制表符分隔，如果没有则尝试空格分隔
                        if "\t" in line:
                            parts = line.split("\t")
                        else:
                            parts = line.split()

                        if len(parts) >= 2:
                            sample_name = parts[0].strip()
                            group_name = parts[1].strip()
                            sample_groups[sample_name] = group_name
                        else:
                            print(f"⚠️  警告: 第{line_num}行格式不正确，已跳过: {line}")

                if sample_groups:
                    print(f"✅ 已加载文本格式样品分组配置: {args.sample_groups_file}")
                    print(f"   共读取 {len(sample_groups)} 个样品分组:")
                    for sample, group in sorted(sample_groups.items()):
                        print(f"     {sample} -> {group}")
                else:
                    raise ValueError("文件中没有找到有效的样品分组信息")

        except Exception as e:
            print(f"❌ 样品分组配置文件加载失败: {e}")
            print("📋 文件格式要求:")
            print("   方式1 (推荐) - 制表符分隔的文本文件:")
            print("     FZY4201    CK")
            print("     FZY4202    Treatment")
            print("     FZY4203    CK")
            print("   方式2 - JSON格式文件:")
            print('     {"FZY4201": "CK", "FZY4202": "Treatment"}')
            sys.exit(1)

    # 创建分析器 | Create analyzer
    analyzer_kwargs = {
        "raw_dir": args.raw_dir,
        "genome_fa": args.genome,
        "output_dir": args.output_dir,
        "target_fa": args.target_fa,
        "target_promoter_fa": args.target_promoter_fa,
        "annotation_gff": args.annotation_gff,
        "enhanced_mode": args.enhanced_mode,
        "threads": args.threads,
        "min_coverage": args.min_coverage,
        "min_cytosines": args.min_cytosines,
        "methylation_diff_threshold": args.methylation_diff_threshold,
        "pvalue_threshold": args.pvalue_threshold,
        "window_size": args.window_size,
        "step_size": args.step_size,
        "merge_distance": args.merge_distance,
        "gene_bins": args.gene_bins,
        "flank_bins": args.flank_bins,
        "flank_size": args.flank_size,
        # 工具路径 | Tool paths
        "fastp_path": args.fastp_path,
        "bismark_path": args.bismark_path,
        "bismark_genome_preparation_path": args.bismark_genome_preparation_path,
        "bowtie2_path": args.bowtie2_path,
        "deduplicate_bismark_path": args.deduplicate_bismark_path,
        "bismark_methylation_extractor_path": args.bismark_methylation_extractor_path,
        "bismark2report_path": args.bismark2report_path,
        "bismark2summary_path": args.bismark2summary_path,
        "makeblastdb_path": args.makeblastdb_path,
        "blastn_path": args.blastn_path,
        "multiqc_path": args.multiqc_path,
        "bedtools_path": args.bedtools_path,
        "samtools_path": args.samtools_path,
        "r_executable": args.r_executable,
    }

    # 添加自定义样品分组 | Add custom sample groups
    if sample_groups:
        analyzer_kwargs["sample_groups"] = sample_groups

    try:
        analyzer = MethylationAnalyzer(**analyzer_kwargs)

        # 试运行模式只验证配置 | Dry run mode only validates configuration
        if args.dry_run:
            print("✅ 配置验证通过 | Configuration validation passed")
            print("📋 分析配置摘要 | Analysis Configuration Summary:")
            print(f"  - 原始数据目录: {analyzer.config.raw_dir}")
            print(f"  - 基因组文件: {analyzer.config.genome_fa}")
            print(f"  - 输出目录: {analyzer.config.output_dir}")
            print(
                f"  - 分析模式: {'增强模式' if analyzer.config.enhanced_mode else '基础模式'}"
            )
            print(f"  - 线程数: {analyzer.config.threads}")
            print(f"  - CK组样品: {analyzer.config.ck_samples}")
            print(f"  - Treatment组样品: {analyzer.config.treatment_samples}")

            if analyzer.config.enhanced_mode:
                print("  - 增强分析参数:")
                print(
                    f"    * 甲基化差异阈值: {analyzer.config.methylation_diff_threshold}"
                )
                print(f"    * p值阈值: {analyzer.config.pvalue_threshold}")
                print(f"    * 窗口大小: {analyzer.config.window_size} bp")
                print(f"    * 基因体bins: {analyzer.config.gene_bins}")
                print(f"    * 侧翼bins: {analyzer.config.flank_bins}")

            print("\n🎯 准备就绪，可以移除 --dry-run 参数开始正式分析")
            return

        # 正常运行分析 | Normal analysis run
        analyzer.run_analysis()

    except Exception as e:
        print(f"❌ 分析失败 | Analysis failed: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
