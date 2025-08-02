"""
核心分析模块 | Core Analysis Module
直接执行所有分析步骤，不生成脚本文件
"""

import glob
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import CommandRunner

try:
    import numpy as np
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


class BasicPipeline:
    """基础流程执行器 | Basic Pipeline Executor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def rename_files(self) -> bool:
        """📝 重命名和清理原始文件 | Rename and clean raw files"""
        self.logger.info("📝 步骤1: 开始重命名和清理原始文件...")

        try:
            # 重命名包含URL参数的文件
            for filename in os.listdir(self.config.raw_dir):
                if "?" in filename:
                    old_path = os.path.join(self.config.raw_dir, filename)
                    new_name = filename.split("?")[0]
                    new_path = os.path.join(self.config.raw_dir, new_name)

                    self.logger.info(f"重命名: {filename} -> {new_name}")
                    shutil.move(old_path, new_path)

            # 清理文件名前缀
            for filename in os.listdir(self.config.raw_dir):
                if filename.startswith("FZYM412_"):
                    old_path = os.path.join(self.config.raw_dir, filename)
                    new_name = filename.replace("FZYM412_", "")
                    new_path = os.path.join(self.config.raw_dir, new_name)

                    self.logger.info(f"去前缀: {filename} -> {new_name}")
                    shutil.move(old_path, new_path)

            self.logger.info("✅ 步骤1: 文件重命名和清理完成")
            return True

        except Exception as e:
            self.logger.error(f"❌ 文件重命名失败: {e}")
            return False

    def run_fastp(self) -> bool:
        """🧹 使用fastp进行质控 | Run fastp for quality control"""
        self.logger.info("🧹 步骤2: 开始使用fastp进行质控...")

        # 创建输出目录
        os.makedirs(self.config.clean_dir, exist_ok=True)

        # 获取所有样品名称 - 改进的样品名提取逻辑
        samples = []
        sample_files = {}

        for filename in os.listdir(self.config.raw_dir):
            if filename.endswith(".fq.gz"):
                if "_1.fq.gz" in filename:
                    # 提取基础样品名（去掉_1.fq.gz）
                    base_name = filename.replace("_1.fq.gz", "")
                    if base_name not in sample_files:
                        sample_files[base_name] = {}
                    sample_files[base_name]["R1"] = filename
                elif "_2.fq.gz" in filename:
                    # 提取基础样品名（去掉_2.fq.gz）
                    base_name = filename.replace("_2.fq.gz", "")
                    if base_name not in sample_files:
                        sample_files[base_name] = {}
                    sample_files[base_name]["R2"] = filename

        # 只保留有完整配对文件的样品
        for base_name, files in sample_files.items():
            if "R1" in files and "R2" in files:
                samples.append(base_name)

        if not samples:
            self.logger.error(
                f"❌ 在 {self.config.raw_dir} 中没有找到配对的 *.fq.gz 文件"
            )
            return False

        self.logger.info(f"🎯 发现样品: {samples}")

        # 处理每个样品
        for sample in samples:
            # 清理样品名称（去掉可能的前缀）
            clean_sample = sample.replace("FZYM412_", "")
            group = self.config.sample_groups.get(clean_sample, "Unknown")
            self.logger.info(
                f"🔧 处理样品: {sample} (清理后: {clean_sample}, 分组: {group})"
            )

            # 构建输入文件路径
            r1 = os.path.join(self.config.raw_dir, sample_files[sample]["R1"])
            r2 = os.path.join(self.config.raw_dir, sample_files[sample]["R2"])

            # 构建输出文件路径 - 使用清理后的样品名，但明确标注R1/R2
            clean_r1 = os.path.join(
                self.config.clean_dir, f"{clean_sample}_R1_clean.fq.gz"
            )
            clean_r2 = os.path.join(
                self.config.clean_dir, f"{clean_sample}_R2_clean.fq.gz"
            )

            # 检查是否已经存在clean文件
            if os.path.exists(clean_r1) and os.path.exists(clean_r2):
                self.logger.info(f"  ✅ 跳过 {clean_sample} - clean文件已存在")
                continue

            if os.path.exists(r1) and os.path.exists(r2):
                cmd = (
                    f"{self.config.fastp_path} "
                    f"-i {r1} -I {r2} "
                    f"-o {clean_r1} -O {clean_r2} "
                    f"-h {os.path.join(self.config.clean_dir, f'{clean_sample}_fastp.html')} "
                    f"-j {os.path.join(self.config.clean_dir, f'{clean_sample}_fastp.json')} "
                    f"--thread {self.config.threads} "
                    f"--detect_adapter_for_pe --correction --cut_front --cut_tail "
                    f"--cut_mean_quality 20 --qualified_quality_phred 20 "
                    f"--unqualified_percent_limit 40 --n_base_limit 10 --length_required 36"
                )

                if not self.cmd_runner.run(cmd, f"fastp处理: {clean_sample}"):
                    return False

                self.logger.info(f"  ✅ fastp完成: {clean_sample}")
            else:
                self.logger.warning(f"  ⚠️ 警告: 找不到配对文件 - {sample}")
                self.logger.warning(f"       期望的文件: {r1}, {r2}")

        self.logger.info("✅ 步骤2: fastp质控完成")
        return True

    def build_bismark_index(self) -> bool:
        """🏗️ 构建Bismark索引 | Build Bismark index"""
        self.logger.info("🏗️ 步骤3: 检查并构建Bismark索引...")

        bismark_index_dir = os.path.join(self.config.mapping_dir, "bismark_index")
        os.makedirs(bismark_index_dir, exist_ok=True)

        # 复制基因组文件到索引目录
        genome_dest = os.path.join(
            bismark_index_dir, os.path.basename(self.config.genome_fa)
        )
        if not os.path.exists(genome_dest):
            shutil.copy2(self.config.genome_fa, genome_dest)

        # 获取bowtie2路径
        bowtie2_path = shutil.which(self.config.bowtie2_path)
        if not bowtie2_path:
            self.logger.error("❌ 找不到bowtie2，请检查安装")
            return False

        bowtie2_dir = os.path.dirname(bowtie2_path)
        self.logger.info(f"🔧 使用bowtie2路径: {bowtie2_dir}")

        # 构建索引
        cmd = (
            f"{self.config.bismark_genome_preparation_path} "
            f"--path_to_aligner {bowtie2_dir} "
            f"--verbose {bismark_index_dir}"
        )

        if self.cmd_runner.run(cmd, "构建Bismark索引"):
            self.logger.info("✅ Bismark索引构建完成")
            return True
        else:
            # 尝试让bismark自动查找
            self.logger.info("🔄 尝试让bismark自动查找bowtie2...")
            cmd_auto = (
                f"{self.config.bismark_genome_preparation_path} "
                f"--verbose {bismark_index_dir}"
            )

            if self.cmd_runner.run(cmd_auto, "构建Bismark索引（自动查找）"):
                self.logger.info("✅ Bismark索引构建完成")
                return True
            else:
                self.logger.error("❌ Bismark索引构建失败")
                return False

    def run_methylation_mapping(self) -> bool:
        """🧬 运行甲基化比对分析 | Run methylation mapping analysis"""
        self.logger.info("🧬 步骤4: 开始甲基化比对分析...")

        # 获取clean样品
        clean_samples = []
        if os.path.exists(self.config.clean_dir):
            for filename in os.listdir(self.config.clean_dir):
                if filename.endswith("_R1_clean.fq.gz"):
                    sample = filename.replace("_R1_clean.fq.gz", "")
                    clean_samples.append(sample)

        if not clean_samples:
            self.logger.error("❌ 没有找到clean数据文件")
            return False

        self.logger.info(f"🎯 发现clean样品: {clean_samples}")

        bismark_index_dir = os.path.join(self.config.mapping_dir, "bismark_index")
        bismark_results_dir = os.path.join(self.config.mapping_dir, "bismark_results")

        for sample in clean_samples:
            group = self.config.sample_groups.get(sample, "Unknown")
            sample_output_dir = os.path.join(bismark_results_dir, sample)

            # 检查该样品是否已经完成甲基化分析
            cx_report_files = glob.glob(
                os.path.join(sample_output_dir, "*CX_report.txt")
            )
            if cx_report_files:
                self.logger.info(
                    f"🚀 跳过甲基化分析: {sample} (分组: {group}) - 已完成"
                )
                continue

            self.logger.info(f"🔬 甲基化分析: {sample} (分组: {group})")

            clean_r1 = os.path.join(self.config.clean_dir, f"{sample}_R1_clean.fq.gz")
            clean_r2 = os.path.join(self.config.clean_dir, f"{sample}_R2_clean.fq.gz")

            os.makedirs(sample_output_dir, exist_ok=True)

            if os.path.exists(clean_r1) and os.path.exists(clean_r2):
                # Step 4.1: Bismark比对
                if not self._run_bismark_alignment(
                    sample, clean_r1, clean_r2, bismark_index_dir, sample_output_dir
                ):
                    continue

                # Step 4.2: 去重复
                if not self._run_deduplication(sample, sample_output_dir):
                    continue

                # Step 4.3: 甲基化提取
                if not self._run_methylation_extraction(
                    sample, sample_output_dir, bismark_index_dir
                ):
                    continue

                # Step 4.4: 生成报告
                self._generate_reports(sample, sample_output_dir)

                self.logger.info(f"✅ 样品分析完成: {sample}")
            else:
                self.logger.warning(f"⚠️ 找不到clean文件 - {sample}")
                self.logger.warning(f"    期望的文件: {clean_r1}, {clean_r2}")

        self.logger.info("✅ 步骤4: 甲基化比对分析完成")
        return True

    def _run_bismark_alignment(
        self, sample, clean_r1, clean_r2, bismark_index_dir, output_dir
    ) -> bool:
        """🎯 运行Bismark比对 | Run Bismark alignment"""
        self.logger.info(f"  🎯 Bismark比对: {sample}")

        # 设置环境变量确保bismark能找到bowtie2
        bowtie2_path = shutil.which(self.config.bowtie2_path)
        if bowtie2_path:
            bowtie2_dir = os.path.dirname(bowtie2_path)
            os.environ["PATH"] = f"{bowtie2_dir}:{os.environ.get('PATH', '')}"

        cmd = (
            f"{self.config.bismark_path} "
            f"--genome {bismark_index_dir} "
            f"--multicore {self.config.threads // 4} "
            f"-1 {clean_r1} -2 {clean_r2} "
            f"--output_dir {output_dir} "
            f"--temp_dir {output_dir} "
            f"--bowtie2 --non_directional"
        )

        return self.cmd_runner.run(cmd, f"Bismark比对: {sample}")

    def _run_deduplication(self, sample, output_dir) -> bool:
        """🔄 运行去重复 | Run deduplication"""
        self.logger.info(f"  🔄 去重复: {sample}")

        bam_files = glob.glob(os.path.join(output_dir, "*_pe.bam"))
        if not bam_files:
            self.logger.error(f"❌ 找不到BAM文件 - {sample}")
            return False

        bam_file = bam_files[0]
        cmd = (
            f"{self.config.deduplicate_bismark_path} "
            f"--paired --output_dir {output_dir} {bam_file}"
        )

        return self.cmd_runner.run(cmd, f"去重复: {sample}")

    def _run_methylation_extraction(
        self, sample, output_dir, bismark_index_dir
    ) -> bool:
        """🧪 运行甲基化提取 | Run methylation extraction"""
        self.logger.info(f"  🧪 甲基化提取: {sample}")

        dedup_bam_files = glob.glob(os.path.join(output_dir, "*deduplicated.bam"))
        if not dedup_bam_files:
            self.logger.error(f"❌ 找不到去重复BAM文件 - {sample}")
            return False

        dedup_bam = dedup_bam_files[0]
        cmd = (
            f"{self.config.bismark_methylation_extractor_path} "
            f"--paired-end --comprehensive --merge_non_CpG --report "
            f"--cytosine_report --genome_folder {bismark_index_dir} "
            f"--output {output_dir} --multicore {self.config.threads} {dedup_bam}"
        )

        return self.cmd_runner.run(cmd, f"甲基化提取: {sample}")

    def _generate_reports(self, sample, output_dir):
        """📊 生成HTML报告 | Generate HTML reports"""
        self.logger.info(f"  📊 生成HTML报告: {sample}")

        # 生成HTML报告
        try:
            alignment_reports = glob.glob(os.path.join(output_dir, "*_PE_report.txt"))
            dedup_reports = glob.glob(
                os.path.join(output_dir, "*_deduplication_report.txt")
            )
            splitting_reports = glob.glob(
                os.path.join(output_dir, "*_splitting_report.txt")
            )
            mbias_reports = glob.glob(os.path.join(output_dir, "*_M-bias.txt"))

            if alignment_reports:
                cmd = (
                    f"{self.config.bismark2report_path} "
                    f"--alignment_report {alignment_reports[0]}"
                )
                if dedup_reports:
                    cmd += f" --dedup_report {dedup_reports[0]}"
                if splitting_reports:
                    cmd += f" --splitting_report {splitting_reports[0]}"
                if mbias_reports:
                    cmd += f" --mbias_report {mbias_reports[0]}"

                self.cmd_runner.run(cmd, f"生成HTML报告: {sample}")
        except Exception as e:
            self.logger.warning(f"⚠️ HTML报告生成可能有问题: {e}")

        # 生成汇总报告 - 需要直接指定BAM文件
        try:
            # 查找所有BAM文件
            bam_files = glob.glob(os.path.join(output_dir, "*.bam"))
            if bam_files:
                # 直接指定BAM文件作为参数
                bam_files_str = " ".join(bam_files)
                cmd = f"{self.config.bismark2summary_path} {bam_files_str}"
                self.cmd_runner.run(cmd, f"生成汇总报告: {sample}")
                self.logger.info(
                    f"  📊 汇总报告已生成，使用BAM文件: {len(bam_files)} 个"
                )
            else:
                self.logger.warning(
                    f"⚠️ 样品 {sample} 没有找到BAM文件，跳过汇总报告生成"
                )

        except Exception as e:
            self.logger.warning(f"⚠️ 汇总报告生成失败: {e}")
            # 如果直接指定BAM文件失败，尝试在样品目录中运行
            try:
                self.logger.info(f"  🔄 尝试在样品目录中运行汇总报告生成...")
                sample_cmd_runner = CommandRunner(self.logger, output_dir)
                cmd = f"{self.config.bismark2summary_path}"
                sample_cmd_runner.run(cmd, f"生成汇总报告(备用方法): {sample}")
            except Exception as e2:
                self.logger.warning(f"⚠️ 备用方法也失败，跳过汇总报告: {e2}")


class EnhancedAnalysis:
    """增强分析执行器 | Enhanced Analysis Executor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_methylkit_analysis(self) -> bool:
        """📊 运行methylKit差异甲基化分析 | Run methylKit differential methylation analysis"""
        self.logger.info("📊 步骤A: 开始methylKit差异甲基化分析...")

        # 创建输出目录
        os.makedirs(
            os.path.join(self.config.analysis_dir, "methylkit_plots"), exist_ok=True
        )
        os.makedirs(
            os.path.join(self.config.analysis_dir, "differential_methylation"),
            exist_ok=True,
        )
        os.makedirs(os.path.join(self.config.analysis_dir, "reports"), exist_ok=True)

        # 获取所有可能的分组比较
        comparisons = self.config.get_pairwise_comparisons()

        if not comparisons:
            self.logger.warning("⚠️  没有足够的分组进行比较分析（至少需要2个分组）")
            return False

        self.logger.info(f"🔬 检测到 {len(comparisons)} 个分组比较组合:")
        for i, (group1, group2) in enumerate(comparisons, 1):
            group1_samples = self.config.get_samples_by_group(group1)
            group2_samples = self.config.get_samples_by_group(group2)
            self.logger.info(
                f"  {i}. {group1} ({len(group1_samples)} 样品) vs {group2} ({len(group2_samples)} 样品)"
            )

        # 执行每个分组比较
        successful_comparisons = 0

        for group1, group2 in comparisons:
            comparison_name = f"{group1}_vs_{group2}"
            self.logger.info(f"🔬 开始分组比较: {comparison_name}")

            if self._run_single_methylkit_comparison(group1, group2, comparison_name):
                successful_comparisons += 1
                self.logger.info(f"✅ 分组比较完成: {comparison_name}")
            else:
                self.logger.warning(f"⚠️  分组比较失败: {comparison_name}")

        # 生成综合分析报告
        if successful_comparisons > 0:
            self._generate_methylkit_summary_report(comparisons, successful_comparisons)
            self.logger.info(
                f"🎉 methylKit分析完成! 成功完成 {successful_comparisons}/{len(comparisons)} 个分组比较"
            )
            return True
        else:
            self.logger.error("❌ 所有分组比较都失败了")
            return False

    def _run_single_methylkit_comparison(
        self, group1: str, group2: str, comparison_name: str
    ) -> bool:
        """🔬 运行单个分组比较的methylKit分析 | Run single group comparison methylKit analysis"""

        # 检查是否已有结果
        diff_sites_file = os.path.join(
            self.config.analysis_dir,
            "differential_methylation",
            f"{comparison_name}_differential_sites.txt",
        )
        if os.path.exists(diff_sites_file):
            self.logger.info(f"🚀 跳过 {comparison_name} 分析 - 结果已存在")
            return True

        # 获取分组样品
        group1_samples = self.config.get_samples_by_group(group1)
        group2_samples = self.config.get_samples_by_group(group2)

        if len(group1_samples) == 0 or len(group2_samples) == 0:
            self.logger.warning(f"⚠️  分组 {group1} 或 {group2} 没有样品，跳过比较")
            return False

        # 构建样品列表和处理组信息
        all_samples = group1_samples + group2_samples
        # group1 设为对照组(0)，group2 设为处理组(1)
        treatments = [0] * len(group1_samples) + [1] * len(group2_samples)

        # 构建R命令字符串
        r_script = f'''
suppressPackageStartupMessages({{
    .libPaths("/share/org/YZWL/yzwl_lixg/miniforge3/envs/methylkit/lib/R/library")
    library(methylKit)
    library(GenomicRanges)
    library(ggplot2)
    library(dplyr)
}})

mapping_dir <- "{os.path.join(self.config.mapping_dir, "bismark_results")}"
analysis_dir <- "{self.config.analysis_dir}"
comparison_name <- "{comparison_name}"
min_coverage <- {self.config.min_coverage}
min_cytosines <- {self.config.min_cytosines}
methylation_diff <- {self.config.methylation_diff_threshold}
pvalue_threshold <- {self.config.pvalue_threshold}

cat("=== methylKit分析开始: {comparison_name} ===\\n")
cat(sprintf("比对结果目录: %s\\n", mapping_dir))
cat(sprintf("分析目录: %s\\n", analysis_dir))
cat(sprintf("比较: {group1} (对照组) vs {group2} (处理组)\\n"))

# 定义样品信息
sample_ids <- c({", ".join([f'"{s}"' for s in all_samples])})
treatments <- c({", ".join(map(str, treatments))})
group_names <- c({", ".join([f'"{group1}"'] * len(group1_samples) + [f'"{group2}"'] * len(group2_samples))})

cat(sprintf("对照组 {group1}: %s\\n", paste({", ".join([f'"{s}"' for s in group1_samples])}, collapse=", ")))
cat(sprintf("处理组 {group2}: %s\\n", paste({", ".join([f'"{s}"' for s in group2_samples])}, collapse=", ")))

# 搜索CX报告文件
file_list <- c()
existing_samples <- c()
existing_treatments <- c()
existing_groups <- c()

for (i in seq_along(sample_ids)) {{
    sample_dir <- file.path(mapping_dir, sample_ids[i])
    if (dir.exists(sample_dir)) {{
        cx_files <- list.files(sample_dir, pattern = "*CX_report.txt$", full.names = TRUE)
        if (length(cx_files) > 0) {{
            file_list <- c(file_list, cx_files[1])
            existing_samples <- c(existing_samples, sample_ids[i])
            existing_treatments <- c(existing_treatments, treatments[i])
            existing_groups <- c(existing_groups, group_names[i])
        }}
    }}
}}

if (length(file_list) == 0) {{
    cat("错误: 没有找到任何CX报告文件\\n")
    quit(status = 1)
}}

cat(sprintf("找到 %d 个样品的数据文件\\n", length(file_list)))

# 检查是否两个组都有样品
group1_count <- sum(existing_groups == "{group1}")
group2_count <- sum(existing_groups == "{group2}")

if (group1_count == 0 || group2_count == 0) {{
    cat(sprintf("错误: 分组不完整 - %s组: %d个样品, %s组: %d个样品\\n", 
                "{group1}", group1_count, "{group2}", group2_count))
    quit(status = 1)
}}

cat(sprintf("有效样品分布 - %s组: %d个, %s组: %d个\\n", 
            "{group1}", group1_count, "{group2}", group2_count))

# 读取甲基化数据
cat("读取甲基化数据...\\n")
myobj <- methRead(file_list, sample.id = existing_samples, assembly = "custom", 
                 treatment = existing_treatments, context = "CpG", mincov = min_coverage)

# 生成描述性统计
comparison_plots_dir <- file.path(analysis_dir, "methylkit_plots", comparison_name)
dir.create(comparison_plots_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(comparison_plots_dir, "methylation_stats.pdf"), width = 12, height = 8)
if (length(myobj) > 0) {{
    getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
    getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
}}
dev.off()

# 过滤数据
cat("过滤低质量数据...\\n")
filtered.myobj <- filterByCoverage(myobj, lo.count = min_coverage, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 合并样品数据
cat("合并样品数据...\\n")
meth <- unite(filtered.myobj, destrand = FALSE)
meth.filtered <- filterByCoverage(meth, lo.count = min_coverage, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 标准化
cat("标准化数据...\\n")
meth.norm <- normalizeCoverage(meth.filtered)

# 聚类分析
pdf(file.path(comparison_plots_dir, "clustering_analysis.pdf"), width = 10, height = 8)
getCorrelation(meth.norm, plot = TRUE)
PCASamples(meth.norm)
clusterSamples(meth.norm, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

# 差异甲基化分析
cat("进行差异甲基化分析...\\n")
myDiff <- calculateDiffMeth(meth.norm)

# 获取差异甲基化位点
diff.meth <- getMethylDiff(myDiff, difference = methylation_diff*100, qvalue = pvalue_threshold, type = "all")
hyper.meth <- getMethylDiff(myDiff, difference = methylation_diff*100, qvalue = pvalue_threshold, type = "hyper")
hypo.meth <- getMethylDiff(myDiff, difference = methylation_diff*100, qvalue = pvalue_threshold, type = "hypo")

# 输出统计信息
cat(sprintf("差异甲基化位点统计 ({comparison_name}):​\\n"))
cat(sprintf("  总差异位点数: %d\\n", nrow(diff.meth)))
cat(sprintf("  高甲基化位点 ({group2} > {group1}): %d\\n", nrow(hyper.meth)))
cat(sprintf("  低甲基化位点 ({group2} < {group1}): %d\\n", nrow(hypo.meth)))

# 保存结果
write.table(diff.meth, file = file.path(analysis_dir, "differential_methylation", paste0(comparison_name, "_differential_sites.txt")), sep = "\\t", quote = FALSE, row.names = FALSE)
write.table(hyper.meth, file = file.path(analysis_dir, "differential_methylation", paste0(comparison_name, "_hypermethylated_sites.txt")), sep = "\\t", quote = FALSE, row.names = FALSE)
write.table(hypo.meth, file = file.path(analysis_dir, "differential_methylation", paste0(comparison_name, "_hypomethylated_sites.txt")), sep = "\\t", quote = FALSE, row.names = FALSE)

# 可视化差异甲基化
pdf(file.path(comparison_plots_dir, "differential_methylation_plots.pdf"), width = 12, height = 8)
plot(myDiff, plot.type = "volcano")
plot(myDiff, plot.type = "scatterplot")
dev.off()

# 生成本次比较的汇总报告
summary_stats <- data.frame(
    Comparison = comparison_name,
    Group1 = "{group1}",
    Group2 = "{group2}",
    Group1_Samples = group1_count,
    Group2_Samples = group2_count,
    Total_CpG_Sites = nrow(meth.norm),
    Differential_Sites_All = nrow(diff.meth),
    Hypermethylated_Sites = nrow(hyper.meth),
    Hypomethylated_Sites = nrow(hypo.meth),
    Methylation_Threshold = methylation_diff,
    Pvalue_Threshold = pvalue_threshold
)

write.table(summary_stats, file = file.path(analysis_dir, "reports", paste0(comparison_name, "_summary.txt")), sep = "\\t", quote = FALSE, row.names = FALSE)

cat(sprintf("methylKit分析完成: %s\\n", comparison_name))
'''

        # 执行R命令
        try:
            result = subprocess.run(
                [
                    self.config.r_executable,
                    "--slave",
                    "--no-restore",
                    "--no-save",
                    "-e",
                    r_script,
                ],
                capture_output=True,
                text=True,
                timeout=1800,
            )  # 30分钟超时

            if result.returncode == 0:
                self.logger.info(f"✅ {comparison_name} methylKit分析完成")
                if result.stdout:
                    for line in result.stdout.split("\n"):
                        if line.strip() and (
                            "methylKit分析" in line
                            or "差异甲基化位点统计" in line
                            or "有效样品分布" in line
                        ):
                            self.logger.info(f"  📊 {line}")
                return True
            else:
                self.logger.error(f"❌ {comparison_name} methylKit分析失败")
                self.logger.error(f"错误输出: {result.stderr}")
                return False

        except subprocess.TimeoutExpired:
            self.logger.error(f"⏰ {comparison_name} methylKit分析超时")
            return False
        except Exception as e:
            self.logger.error(f"❌ {comparison_name} methylKit分析异常: {e}")
            return False

    def _generate_methylkit_summary_report(
        self, comparisons: List[Tuple[str, str]], successful_count: int
    ):
        """📋 生成methylKit综合分析报告 | Generate methylKit comprehensive analysis report"""
        try:
            summary_report_file = os.path.join(
                self.config.analysis_dir,
                "reports",
                "methylkit_comprehensive_summary.txt",
            )

            with open(summary_report_file, "w", encoding="utf-8") as f:
                f.write("🔬 methylKit 差异甲基化分析综合报告\n")
                f.write("=" * 50 + "\n")
                f.write(
                    f"生成时间: {subprocess.run(['date'], capture_output=True, text=True).stdout.strip()}\n\n"
                )

                f.write("📊 分析概览:\n")
                f.write(f"  总分组比较数: {len(comparisons)}\n")
                f.write(f"  成功完成数: {successful_count}\n")
                f.write(f"  失败数: {len(comparisons) - successful_count}\n\n")

                f.write("🔬 分组比较详情:\n")
                for i, (group1, group2) in enumerate(comparisons, 1):
                    comparison_name = f"{group1}_vs_{group2}"
                    group1_samples = self.config.get_samples_by_group(group1)
                    group2_samples = self.config.get_samples_by_group(group2)

                    # 检查结果文件是否存在
                    diff_sites_file = os.path.join(
                        self.config.analysis_dir,
                        "differential_methylation",
                        f"{comparison_name}_differential_sites.txt",
                    )
                    status = "✅ 成功" if os.path.exists(diff_sites_file) else "❌ 失败"

                    f.write(f"  {i}. {comparison_name}: {status}\n")
                    f.write(
                        f"     {group1} 组: {', '.join(group1_samples)} ({len(group1_samples)} 样品)\n"
                    )
                    f.write(
                        f"     {group2} 组: {', '.join(group2_samples)} ({len(group2_samples)} 样品)\n"
                    )

                    # 如果有结果文件，读取统计信息
                    if os.path.exists(diff_sites_file):
                        try:
                            with open(diff_sites_file, "r") as diff_f:
                                lines = diff_f.readlines()
                                if len(lines) > 1:  # 有数据行
                                    f.write(f"     差异位点数: {len(lines) - 1}\n")
                        except:
                            pass
                    f.write("\n")

                f.write("📁 结果文件位置:\n")
                f.write(
                    f"  差异甲基化结果: {os.path.join(self.config.analysis_dir, 'differential_methylation')}\n"
                )
                f.write(
                    f"  可视化图表: {os.path.join(self.config.analysis_dir, 'methylkit_plots')}\n"
                )
                f.write(
                    f"  分析报告: {os.path.join(self.config.analysis_dir, 'reports')}\n\n"
                )

                f.write("📋 文件命名规则:\n")
                f.write(
                    "  - {Group1}_vs_{Group2}_differential_sites.txt: 所有差异位点\n"
                )
                f.write(
                    "  - {Group1}_vs_{Group2}_hypermethylated_sites.txt: 高甲基化位点\n"
                )
                f.write(
                    "  - {Group1}_vs_{Group2}_hypomethylated_sites.txt: 低甲基化位点\n"
                )
                f.write("  - {Group1}_vs_{Group2}_summary.txt: 比较统计摘要\n\n")

                f.write("💡 分析参数:\n")
                f.write(f"  最小覆盖度: {self.config.min_coverage}\n")
                f.write(
                    f"  甲基化差异阈值: {self.config.methylation_diff_threshold} ({self.config.methylation_diff_threshold * 100}%)\n"
                )
                f.write(f"  p值阈值: {self.config.pvalue_threshold}\n")

            self.logger.info(f"📋 methylKit综合分析报告已生成: {summary_report_file}")

        except Exception as e:
            self.logger.warning(f"⚠️  综合分析报告生成失败: {e}")

    def run_blast_analysis(self) -> bool:
        """🎯 运行BLAST定位分析 | Run BLAST localization analysis"""
        if not self.config.target_promoter_fa:
            return True

        self.logger.info("🎯 步骤F: 开始目标序列BLAST定位和甲基化分析...")

        blast_output_dir = os.path.join(self.config.analysis_dir, "blast_results")
        os.makedirs(blast_output_dir, exist_ok=True)

        # 检查BLAST工具
        if not (
            shutil.which(self.config.makeblastdb_path)
            and shutil.which(self.config.blastn_path)
        ):
            self.logger.warning("⚠️  BLAST工具未安装，跳过目标序列定位")
            return True

        # 构建BLAST数据库
        blast_db = os.path.join(blast_output_dir, "genome_blastdb")
        cmd_makedb = f"{self.config.makeblastdb_path} -in {self.config.genome_fa} -dbtype nucl -out {blast_db}"

        if not self.cmd_runner.run(cmd_makedb, "构建BLAST数据库"):
            return False

        # 执行BLAST搜索
        blast_result = os.path.join(blast_output_dir, "target_blast_result.txt")
        cmd_blast = (
            f"{self.config.blastn_path} -query {self.config.target_promoter_fa} "
            f"-db {blast_db} -out {blast_result} "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
            f"-evalue 1e-5 -max_target_seqs 10 -word_size 11"
        )

        if not self.cmd_runner.run(cmd_blast, "执行BLAST搜索"):
            return False

        # 解析BLAST结果并提取甲基化信息
        positions_file = os.path.join(blast_output_dir, "target_positions.bed")
        if self._parse_blast_results(blast_result, positions_file):
            self.logger.info("✅ BLAST定位完成")

            # 提取目标区域甲基化信息
            target_analysis_dir = os.path.join(
                self.config.analysis_dir, "target_sequence_analysis"
            )
            os.makedirs(target_analysis_dir, exist_ok=True)

            if self._extract_target_methylation(positions_file, target_analysis_dir):
                self.logger.info("✅ 目标区域甲基化提取完成")
                self._generate_target_methylation_report(target_analysis_dir)

            return True
        else:
            self.logger.warning("⚠️  BLAST结果解析失败")
            return False

    def _parse_blast_results(self, blast_result, positions_file) -> bool:
        """📄 解析BLAST结果 | Parse BLAST results"""
        try:
            with open(blast_result, "r") as f, open(positions_file, "w") as bed:
                bed.write(
                    "# chrom\tstart\tend\tname\tscore\tstrand\tidentity\tlength\tevalue\n"
                )

                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 14:
                        qseqid = fields[0]
                        sseqid = fields[1]
                        pident = float(fields[2])
                        length = int(fields[3])
                        sstart = int(fields[8])
                        send = int(fields[9])
                        evalue = float(fields[10])
                        bitscore = float(fields[11])

                        # 确定方向和坐标
                        if sstart <= send:
                            start = sstart - 1
                            end = send
                            strand = "+"
                        else:
                            start = send - 1
                            end = sstart
                            strand = "-"

                        bed.write(
                            f"{sseqid}\t{start}\t{end}\t{qseqid}\t{bitscore}\t{strand}\t{pident}\t{length}\t{evalue}\n"
                        )

            return os.path.getsize(positions_file) > 100  # 检查文件是否有内容

        except Exception as e:
            self.logger.error(f"❌ BLAST结果解析失败: {e}")
            return False

    def _extract_target_methylation(self, positions_file, target_analysis_dir) -> bool:
        """🧪 提取目标区域甲基化信息 | Extract target region methylation"""
        try:
            # 获取所有样品的CX报告文件
            samples = []
            bismark_results_dir = os.path.join(
                self.config.mapping_dir, "bismark_results"
            )
            for sample in self.config.sample_groups.keys():
                cx_files = glob.glob(
                    os.path.join(bismark_results_dir, sample, "*CX_report.txt")
                )
                if cx_files:
                    samples.append((sample, cx_files[0]))

            if not samples:
                self.logger.warning("⚠️ 未找到甲基化数据文件")
                return False

            # 读取位置信息
            positions = []
            with open(positions_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 4:
                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])
                        name = fields[3]
                        positions.append((chrom, start, end, name))

            # 为每个样品提取目标区域甲基化
            for sample, cx_report in samples:
                output_file = os.path.join(
                    target_analysis_dir, f"{sample}_target_region_methylation.txt"
                )

                with open(cx_report, "r") as f, open(output_file, "w") as out:
                    out.write(
                        "chrom\tposition\tstrand\tmeth_count\tunmeth_count\tcontext\ttotal_count\tmeth_level\ttarget_region\tsample\n"
                    )

                    for line in f:
                        fields = line.strip().split("\t")
                        if len(fields) >= 6:
                            chrom = fields[0]
                            pos = int(fields[1])
                            strand = fields[2]
                            meth_count = int(fields[3])
                            unmeth_count = int(fields[4])
                            context = fields[5]

                            # 检查是否在目标区域内
                            for (
                                target_chrom,
                                target_start,
                                target_end,
                                target_name,
                            ) in positions:
                                if (
                                    chrom == target_chrom
                                    and target_start <= pos <= target_end
                                ):
                                    total_count = meth_count + unmeth_count
                                    if total_count > 0:
                                        meth_level = (meth_count / total_count) * 100
                                        out.write(
                                            f"{chrom}\t{pos}\t{strand}\t{meth_count}\t{unmeth_count}\t{context}\t{total_count}\t{meth_level}\t{target_name}\t{sample}\n"
                                        )
                                    break

            return True

        except Exception as e:
            self.logger.error(f"❌ 目标区域甲基化信息提取失败: {e}")
            return False

    def _generate_target_methylation_report(self, target_analysis_dir):
        """📋 生成目标序列甲基化报告 | Generate target methylation report"""
        if not HAS_PANDAS:
            self.logger.warning("⚠️  pandas未安装，跳过目标序列甲基化报告生成")
            self.logger.warning("    安装命令: pip install pandas numpy")
            return

        try:
            report_file = os.path.join(
                target_analysis_dir, "target_methylation_report.txt"
            )

            # 收集所有甲基化数据文件
            meth_files = glob.glob(
                os.path.join(target_analysis_dir, "*_target_region_methylation.txt")
            )

            if not meth_files:
                self.logger.warning("⚠️ 未找到目标区域甲基化数据文件")
                return

            # 读取和合并数据
            all_data = []
            for file in meth_files:
                try:
                    df = pd.read_csv(file, sep="\t")
                    if not df.empty:
                        all_data.append(df)
                except Exception as e:
                    self.logger.error(f"❌ 读取文件失败 {file}: {e}")

            if not all_data:
                self.logger.warning("⚠️ 没有有效的甲基化数据")
                return

            # 合并所有数据
            combined_data = pd.concat(all_data, ignore_index=True)
            cpg_data = combined_data[combined_data["context"] == "CpG"].copy()

            if cpg_data.empty:
                self.logger.warning("⚠️ 没有CpG甲基化数据")
                return

            # 生成报告
            with open(report_file, "w", encoding="utf-8") as report:
                report.write("🎯 目标序列甲基化分析报告\n")
                report.write("=" * 50 + "\n")
                report.write(f"生成时间: {pd.Timestamp.now()}\n\n")

                # 基本统计
                report.write("=== 📊 基本统计 ===\n")
                report.write(f"总CpG位点数: {len(cpg_data)}\n")
                report.write(f"样品数: {cpg_data['sample'].nunique()}\n")
                report.write(f"目标区域数: {cpg_data['target_region'].nunique()}\n")
                report.write(f"染色体数: {cpg_data['chrom'].nunique()}\n\n")

                # 样品统计
                report.write("=== 🧬 样品统计 ===\n")
                for sample in cpg_data["sample"].unique():
                    sample_data = cpg_data[cpg_data["sample"] == sample]
                    report.write(f"样品 {sample}:\n")
                    report.write(f"  CpG位点数: {len(sample_data)}\n")
                    report.write(
                        f"  平均甲基化水平: {sample_data['meth_level'].mean():.2f}%\n"
                    )
                    report.write(
                        f"  甲基化水平标准差: {sample_data['meth_level'].std():.2f}%\n"
                    )
                    report.write(
                        f"  平均覆盖度: {sample_data['total_count'].mean():.1f}\n\n"
                    )

                # 分组比较
                cpg_data["group"] = cpg_data["sample"].map(self.config.sample_groups)

                if cpg_data["group"].notna().any():
                    report.write("=== 📈 分组比较 ===\n")
                    for group in cpg_data["group"].dropna().unique():
                        group_data = cpg_data[cpg_data["group"] == group]
                        report.write(f"{group}组:\n")
                        report.write(f"  CpG位点数: {len(group_data)}\n")
                        report.write(
                            f"  平均甲基化水平: {group_data['meth_level'].mean():.2f}%\n"
                        )
                        report.write(
                            f"  甲基化水平标准差: {group_data['meth_level'].std():.2f}%\n\n"
                        )

            self.logger.info(f"✅ 目标序列甲基化报告已生成: {report_file}")

        except Exception as e:
            self.logger.error(f"❌ 目标序列甲基化报告生成失败: {e}")

    def run_genomic_bins_analysis(self) -> bool:
        """📦 运行基因组bins甲基化分析 | Run genomic bins methylation analysis"""
        if not self.config.annotation_gff:
            self.logger.warning("⚠️  基因注释文件不存在，跳过bins分析")
            return True

        if not HAS_PANDAS:
            self.logger.warning("⚠️  pandas未安装，跳过基因组bins分析")
            self.logger.warning("    安装命令: pip install pandas numpy")
            return True

        self.logger.info("📦 步骤C: 开始基因组bins甲基化分析...")

        bins_output_dir = os.path.join(self.config.analysis_dir, "genomic_bins")
        os.makedirs(bins_output_dir, exist_ok=True)

        # 检查是否已有结果
        bins_result_file = os.path.join(bins_output_dir, "genomic_bins_methylation.txt")
        if os.path.exists(bins_result_file):
            self.logger.info("🚀 跳过基因组bins分析 - 结果已存在")
            return True

        try:
            # 直接在Python中执行bins分析
            return self._execute_bins_analysis(bins_output_dir)

        except Exception as e:
            self.logger.error(f"❌ 基因组bins分析失败: {e}")
            return False

    def _execute_bins_analysis(self, output_dir) -> bool:
        """🔧 执行bins分析的具体实现 | Execute bins analysis implementation"""
        self.logger.info("  🔧 执行基因组bins分析...")

        # 读取基因注释
        genes = self._read_gff_genes()
        if not genes:
            self.logger.warning("⚠️ 没有找到基因注释，跳过bins分析")
            return False

        self.logger.info(f"🧬 找到 {len(genes)} 个基因")

        # 找到所有CX报告文件
        cx_files = []
        bismark_results_dir = os.path.join(self.config.mapping_dir, "bismark_results")
        for sample in self.config.sample_groups.keys():
            sample_dir = os.path.join(bismark_results_dir, sample)
            if os.path.exists(sample_dir):
                for file in os.listdir(sample_dir):
                    if file.endswith("CX_report.txt"):
                        cx_files.append((sample, os.path.join(sample_dir, file)))
                        break

        if not cx_files:
            self.logger.warning("⚠️ 没有找到CX报告文件")
            return False

        self.logger.info(f"📊 找到 {len(cx_files)} 个CX报告文件")

        # 处理每个样品
        all_results = []

        for sample, cx_file in cx_files:
            self.logger.info(f"  🔬 处理样品: {sample}")
            self.logger.info(print("=" * 80))

            # 读取甲基化数据
            methylation_data = self._read_cx_report(cx_file)
            if methylation_data is None or (HAS_PANDAS and methylation_data.empty):
                self.logger.warning(f"⚠️ 跳过样品 {sample} (无有效数据)")
                continue

            self.logger.info(f"    📊 样品 {sample}: {len(methylation_data)} 个CpG位点")

            # 为每个基因计算bins甲基化
            for gene in genes:
                bins = self._create_gene_bins(gene)
                bin_methylation = self._calculate_bin_methylation(
                    bins, methylation_data
                )

                # 添加样品信息
                for result in bin_methylation:
                    result["sample"] = sample
                    result["gene_name"] = gene["gene_name"]
                    all_results.append(result)

        # 保存结果
        if all_results:
            results_df = pd.DataFrame(all_results)
            output_file = os.path.join(output_dir, "genomic_bins_methylation.txt")
            results_df.to_csv(output_file, sep="\t", index=False)
            self.logger.info(f"✅ 结果已保存到: {output_file}")

            # 生成统计摘要
            self._generate_bins_summary(results_df, output_dir)
            return True
        else:
            self.logger.warning("⚠️ 没有生成任何结果")
            return False

    def _read_gff_genes(self) -> List[Dict]:
        """📄 读取GFF文件中的基因信息 | Read gene information from GFF file"""
        genes = []
        try:
            with open(self.config.annotation_gff, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 9 and parts[2] == "gene":
                        chrom = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6]

                        # 解析属性
                        attrs = {}
                        for attr in parts[8].split(";"):
                            if "=" in attr:
                                key, value = attr.split("=", 1)
                                attrs[key] = value

                        gene_id = attrs.get("ID", f"gene_{len(genes)}")
                        gene_name = attrs.get("Name", gene_id)

                        genes.append(
                            {
                                "chrom": chrom,
                                "start": start,
                                "end": end,
                                "strand": strand,
                                "gene_id": gene_id,
                                "gene_name": gene_name,
                            }
                        )
        except Exception as e:
            self.logger.error(f"❌ 读取GFF文件失败: {e}")
            return []

        return genes

    def _read_cx_report(self, cx_file):
        """📊 读取CX report文件 | Read CX report file"""
        if not HAS_PANDAS:
            self.logger.error("❌ pandas未安装，无法读取CX报告文件")
            return None

        try:
            df = pd.read_csv(
                cx_file,
                sep="\t",
                header=None,
                names=[
                    "chr",
                    "position",
                    "strand",
                    "count_methylated",
                    "count_unmethylated",
                    "context",
                    "sequence",
                ],
            )

            # 只保留CpG上下文
            df = df[df["context"] == "CpG"].copy()

            # 计算甲基化水平
            df["total_count"] = df["count_methylated"] + df["count_unmethylated"]
            df["methylation_level"] = df["count_methylated"] / df["total_count"] * 100

            # 过滤低覆盖度位点
            df = df[df["total_count"] >= self.config.min_coverage].copy()

            return df
        except Exception as e:
            self.logger.error(f"❌ 读取CX报告文件失败 {cx_file}: {e}")
            return None if not HAS_PANDAS else pd.DataFrame()

    def _create_gene_bins(self, gene) -> List[Dict]:
        """🔢 为基因创建bins | Create bins for gene"""
        bins = []

        # 上游区域bins
        upstream_start = gene["start"] - self.config.flank_size
        upstream_end = gene["start"]
        upstream_bin_size = self.config.flank_size // self.config.flank_bins

        for i in range(self.config.flank_bins):
            bin_start = upstream_start + i * upstream_bin_size
            bin_end = upstream_start + (i + 1) * upstream_bin_size
            bins.append(
                {
                    "chrom": gene["chrom"],
                    "start": max(1, bin_start),  # 避免负坐标
                    "end": bin_end,
                    "region": "upstream",
                    "bin_id": i,
                    "gene_id": gene["gene_id"],
                }
            )

        # 基因体bins
        gene_length = gene["end"] - gene["start"] + 1
        gene_bin_size = gene_length / self.config.gene_bins

        for i in range(self.config.gene_bins):
            bin_start = gene["start"] + int(i * gene_bin_size)
            bin_end = gene["start"] + int((i + 1) * gene_bin_size)
            bins.append(
                {
                    "chrom": gene["chrom"],
                    "start": bin_start,
                    "end": min(bin_end, gene["end"]),
                    "region": "gene_body",
                    "bin_id": i,
                    "gene_id": gene["gene_id"],
                }
            )

        # 下游区域bins
        downstream_start = gene["end"]
        downstream_end = gene["end"] + self.config.flank_size
        downstream_bin_size = self.config.flank_size // self.config.flank_bins

        for i in range(self.config.flank_bins):
            bin_start = downstream_start + i * downstream_bin_size
            bin_end = downstream_start + (i + 1) * downstream_bin_size
            bins.append(
                {
                    "chrom": gene["chrom"],
                    "start": bin_start,
                    "end": bin_end,
                    "region": "downstream",
                    "bin_id": i,
                    "gene_id": gene["gene_id"],
                }
            )

        return bins

    def _calculate_bin_methylation(self, bins, methylation_data) -> List[Dict]:
        """🧮 计算每个bin的甲基化水平 | Calculate methylation level for each bin"""
        bin_methylation = []

        for bin_info in bins:
            chrom = bin_info["chrom"]
            start = bin_info["start"]
            end = bin_info["end"]

            # 找到该bin内的甲基化位点
            bin_sites = methylation_data[
                (methylation_data["chr"] == chrom)
                & (methylation_data["position"] >= start)
                & (methylation_data["position"] <= end)
            ].copy()

            if len(bin_sites) > 0:
                # 计算加权平均甲基化水平
                total_methylated = bin_sites["count_methylated"].sum()
                total_unmethylated = bin_sites["count_unmethylated"].sum()
                total_count = total_methylated + total_unmethylated

                if total_count > 0:
                    avg_methylation = (total_methylated / total_count) * 100
                else:
                    if HAS_PANDAS:
                        avg_methylation = np.nan
                    else:
                        avg_methylation = None

                num_sites = len(bin_sites)
            else:
                if HAS_PANDAS:
                    avg_methylation = np.nan
                else:
                    avg_methylation = None
                num_sites = 0
                total_count = 0

            bin_methylation.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "region": bin_info["region"],
                    "bin_id": bin_info["bin_id"],
                    "gene_id": bin_info["gene_id"],
                    "avg_methylation": avg_methylation,
                    "num_sites": num_sites,
                    "total_coverage": total_count,
                }
            )

        return bin_methylation

    def _generate_bins_summary(self, results_df, output_dir):
        """📈 生成bins分析统计摘要 | Generate bins analysis summary"""
        if not HAS_PANDAS:
            self.logger.warning("⚠️  pandas未安装，跳过bins统计摘要生成")
            return

        try:
            summary_stats = []
            for sample in results_df["sample"].unique():
                sample_data = results_df[results_df["sample"] == sample]

                for region in ["upstream", "gene_body", "downstream"]:
                    region_data = sample_data[sample_data["region"] == region]
                    valid_data = region_data.dropna(subset=["avg_methylation"])

                    if len(valid_data) > 0:
                        summary_stats.append(
                            {
                                "sample": sample,
                                "region": region,
                                "num_bins": len(region_data),
                                "bins_with_data": len(valid_data),
                                "avg_methylation": valid_data["avg_methylation"].mean(),
                                "median_methylation": valid_data[
                                    "avg_methylation"
                                ].median(),
                            }
                        )

            summary_df = pd.DataFrame(summary_stats)
            summary_file = os.path.join(output_dir, "bins_methylation_summary.txt")
            summary_df.to_csv(summary_file, sep="\t", index=False)
            self.logger.info(f"📊 统计摘要已保存到: {summary_file}")

        except Exception as e:
            self.logger.error(f"❌ 生成bins统计摘要失败: {e}")
