"""
报告生成模块 | Reports Generation Module
包含各种分析报告和统计的生成功能
"""

import glob
import os
from datetime import datetime


class ReportGenerator:
    """报告生成器 | Report Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_comprehensive_report(self) -> bool:
        """📋 生成综合分析报告 | Generate comprehensive analysis report"""
        self.logger.info("📋 步骤G: 生成综合分析报告...")

        reports_dir = os.path.join(self.config.output_dir, "reports")
        os.makedirs(reports_dir, exist_ok=True)

        report_file = os.path.join(reports_dir, "comprehensive_methylation_report.txt")

        try:
            with open(report_file, "w", encoding="utf-8") as f:
                self._write_report_header(f)
                self._write_sample_info(f)
                self._write_analysis_results(f)
                self._write_output_files(f)
                self._write_recommendations(f)

            self.logger.info(f"✅ 综合分析报告已生成: {report_file}")
            return True

        except Exception as e:
            self.logger.error(f"❌ 报告生成失败: {e}")
            return False

    def _write_report_header(self, f):
        """📝 写入报告头部 | Write report header"""
        f.write("================================\n")
        f.write("🧬 甲基化分析综合报告\n")
        f.write("Comprehensive Methylation Analysis Report\n")
        f.write("================================\n")
        f.write(
            f"📅 生成时间 | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        )
        f.write(
            f"🔧 分析模式 | Analysis Mode: {'Enhanced' if self.config.enhanced_mode else 'Basic'}\n\n"
        )

        f.write("⚙️ 分析参数 | Analysis Parameters:\n")
        f.write(f"  - 线程数 | Threads: {self.config.threads}\n")
        f.write(f"  - 最小覆盖度 | Min Coverage: {self.config.min_coverage}\n")
        f.write(f"  - 最小胞嘧啶数 | Min Cytosines: {self.config.min_cytosines}\n")
        if self.config.enhanced_mode:
            f.write(
                f"  - 甲基化差异阈值 | Methylation Diff Threshold: {self.config.methylation_diff_threshold} ({self.config.methylation_diff_threshold * 100}%)\n"
            )
            f.write(
                f"  - p值阈值 | P-value Threshold: {self.config.pvalue_threshold}\n"
            )
            f.write(f"  - 窗口大小 | Window Size: {self.config.window_size} bp\n")
            f.write(f"  - 步长 | Step Size: {self.config.step_size} bp\n")
            f.write(f"  - 合并距离 | Merge Distance: {self.config.merge_distance} bp\n")
        f.write("\n")

    def _write_sample_info(self, f):
        """🧬 写入样品信息 | Write sample information"""
        f.write("=== 🧬 样品信息 | Sample Information ===\n")

        all_groups = self.config.all_groups
        for group in all_groups:
            group_samples = self.config.get_samples_by_group(group)
            f.write(
                f"{group}组: {', '.join(group_samples)} ({len(group_samples)} 样品)\n"
            )

        # 统计实际处理的样品
        if os.path.exists(self.config.clean_dir):
            clean_samples = []
            for filename in os.listdir(self.config.clean_dir):
                if filename.endswith("_R1_clean.fq.gz"):
                    sample = filename.replace("_R1_clean.fq.gz", "")
                    clean_samples.append(sample)
            f.write(f"实际处理样品数 | Actually Processed: {len(clean_samples)}\n")
            if clean_samples:
                f.write(
                    f"处理的样品列表 | Processed samples: {', '.join(clean_samples)}\n"
                )
        f.write("\n")

    def _write_analysis_results(self, f):
        """📊 写入分析结果统计 | Write analysis results statistics"""
        f.write("=== 📊 分析结果统计 | Analysis Results Statistics ===\n")

        # 基础分析统计
        bismark_results_dir = os.path.join(self.config.mapping_dir, "bismark_results")
        total_samples = len(
            glob.glob(os.path.join(bismark_results_dir, "*", "*CX_report.txt"))
        )
        f.write(f"成功处理样品数 | Successfully Processed Samples: {total_samples}\n")

        # 增强分析统计
        if (
            self.config.enhanced_mode
            and hasattr(self.config, "analysis_dir")
            and os.path.exists(self.config.analysis_dir)
        ):
            # methylKit分析结果
            methylkit_summary = os.path.join(
                self.config.analysis_dir, "reports", "methylkit_summary.txt"
            )
            if os.path.exists(methylkit_summary):
                f.write("\n📊 methylKit分析结果 | methylKit Analysis Results:\n")
                with open(methylkit_summary, "r") as summary_f:
                    f.write(summary_f.read())

            # BLAST定位结果
            positions_file = os.path.join(
                self.config.analysis_dir, "blast_results", "target_positions.bed"
            )
            if os.path.exists(positions_file):
                try:
                    with open(positions_file, "r") as pos_f:
                        blast_matches = sum(
                            1 for line in pos_f if not line.startswith("#")
                        )
                    f.write(
                        f"🎯 BLAST定位结果 | BLAST Localization Results: {blast_matches} 个匹配区域\n"
                    )
                except:
                    pass

        f.write("\n")

    def _write_output_files(self, f):
        """📁 写入输出文件位置 | Write output file locations"""
        f.write("=== 📁 输出文件位置 | Output File Locations ===\n")

        f.write("1. 📊 基础分析结果 | Basic Analysis Results:\n")
        f.write(f"   - Clean数据 | Clean Data: {self.config.clean_dir}\n")
        f.write(f"   - 比对结果 | Alignment Results: {self.config.mapping_dir}\n")

        if self.config.enhanced_mode and hasattr(self.config, "analysis_dir"):
            f.write("\n2. 🔬 增强分析结果 | Enhanced Analysis Results:\n")
            f.write(
                f"   - 差异甲基化分析: {os.path.join(self.config.analysis_dir, 'differential_methylation')}\n"
            )
            f.write(
                f"   - 基因组bins分析: {os.path.join(self.config.analysis_dir, 'genomic_bins')}\n"
            )
            f.write(
                f"   - 目标序列分析: {os.path.join(self.config.analysis_dir, 'target_sequence_analysis')}\n"
            )
            f.write(
                f"   - BLAST定位结果: {os.path.join(self.config.analysis_dir, 'blast_results')}\n"
            )

        f.write(f"\n3. 📈 可视化结果 | Visualization Results:\n")
        if self.config.enhanced_mode and hasattr(self.config, "analysis_dir"):
            f.write(
                f"   - methylKit图表: {os.path.join(self.config.analysis_dir, 'methylkit_plots')}\n"
            )

        multiqc_dir = os.path.join(self.config.mapping_dir, "multiqc_report")
        if os.path.exists(multiqc_dir):
            f.write(f"   - 质量控制图表: {multiqc_dir}\n")

        f.write("\n")

    def _write_recommendations(self, f):
        """💡 写入分析建议 | Write analysis recommendations"""
        f.write("=== 💡 分析建议 | Analysis Recommendations ===\n")
        f.write("1. 查看差异甲基化位点的基因组分布和功能注释\n")
        f.write(
            "   Review genomic distribution and functional annotation of differential methylation sites\n\n"
        )

        f.write("2. 分析基因体及其上下游区域的甲基化模式\n")
        f.write(
            "   Analyze methylation patterns in gene bodies and flanking regions\n\n"
        )

        f.write("3. 结合目标基因的甲基化变化分析生物学意义\n")
        f.write(
            "   Combine target gene methylation changes to analyze biological significance\n\n"
        )

        f.write("4. 考虑进行基因本体(GO)或KEGG通路富集分析\n")
        f.write(
            "   Consider Gene Ontology (GO) or KEGG pathway enrichment analysis\n\n"
        )

        f.write("=== 🚀 下一步分析 | Next Steps ===\n")
        f.write("建议使用R/Python进一步分析:\n")
        f.write("Recommended further analysis using R/Python:\n")
        f.write("1. 📊 可视化甲基化模式 | Visualize methylation patterns\n")
        f.write("2. 🔍 功能富集分析 | Functional enrichment analysis\n")
        f.write("3. 🧬 与转录组数据整合分析 | Integration with transcriptome data\n")
        f.write(
            "4. ✅ 验证关键差异甲基化区域 | Validate key differential methylation regions\n"
        )
