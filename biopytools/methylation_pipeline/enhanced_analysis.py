"""
增强分析模块 | Enhanced Analysis Module
包含methylKit差异分析、基因组注释、BLAST定位等高级功能
"""

import os
import subprocess
from pathlib import Path

from .utils import CommandRunner


class MethylKitAnalyzer:
    """methylKit分析器 | methylKit Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_methylkit_analysis(self) -> bool:
        """运行methylKit差异甲基化分析 | Run methylKit differential methylation analysis"""
        if not self.config.enhanced_mode:
            return True

        self.logger.info("步骤A: 开始methylKit差异甲基化分析...")

        # 创建methylKit分析目录
        methylkit_dir = os.path.join(self.config.analysis_dir, "methylkit_analysis")
        os.makedirs(methylkit_dir, exist_ok=True)

        # 检查是否已有结果
        diff_sites_file = os.path.join(
            self.config.analysis_dir,
            "differential_methylation",
            "all_differential_sites.txt",
        )
        if os.path.exists(diff_sites_file):
            self.logger.info("🚀 跳过methylKit分析 - 结果已存在")
            return True

        # 创建R脚本
        r_script_path = os.path.join(methylkit_dir, "methylkit_analysis.R")
        self._create_methylkit_script(r_script_path)

        # 运行R脚本
        cmd = (
            f"{self.config.r_executable} --slave --no-restore --no-save "
            f"--file={r_script_path} --args "
            f"{self.config.mapping_dir} {self.config.analysis_dir} "
            f"{self.config.min_coverage} {self.config.min_cytosines} "
            f"{self.config.methylation_diff_threshold} {self.config.pvalue_threshold}"
        )

        if self.cmd_runner.run(cmd, "methylKit差异甲基化分析"):
            self.logger.info("✅ methylKit分析完成")
            return True
        else:
            self.logger.error("❌ methylKit分析失败")
            return False

    def _create_methylkit_script(self, script_path):
        """创建methylKit R脚本 | Create methylKit R script"""
        script_content = """#!/usr/bin/env Rscript

# methylKit分析脚本
suppressPackageStartupMessages({
    library(methylKit)
    library(GenomicRanges)
    library(ggplot2)
    library(dplyr)
})

# 从命令行获取参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    cat("Usage: Rscript methylkit_analysis.R <mapping_dir> <analysis_dir> <min_coverage> <min_cytosines> <methylation_diff> <pvalue_threshold>\\n")
    quit(status = 1)
}

mapping_dir <- args[1]
analysis_dir <- args[2] 
min_coverage <- as.numeric(args[3])
min_cytosines <- as.numeric(args[4])
methylation_diff <- as.numeric(args[5])
pvalue_threshold <- as.numeric(args[6])

cat("=== methylKit分析开始 ===\\n")
cat("参数设置:\\n")
cat(sprintf("  映射目录: %s\\n", mapping_dir))
cat(sprintf("  分析目录: %s\\n", analysis_dir))
cat(sprintf("  最小覆盖度: %d\\n", min_coverage))
cat(sprintf("  最小胞嘧啶数: %d\\n", min_cytosines))
cat(sprintf("  甲基化差异阈值: %.2f\\n", methylation_diff))
cat(sprintf("  p值阈值: %.3f\\n", pvalue_threshold))

# 确保输出目录存在
output_dirs <- c(
    file.path(analysis_dir, "methylkit_analysis"),
    file.path(analysis_dir, "differential_methylation"),
    file.path(analysis_dir, "reports")
)

for (dir_path in output_dirs) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
}

# 定义样品信息
sample_ids <- c("FZY4201", "FZY4203", "FZY4205", "FZY4207", 
                "FZY4202", "FZY4204", "FZY4206", "FZY4208")
treatments <- c(0, 0, 0, 0, 1, 1, 1, 1)  # 0=CK, 1=Treatment

# 构建文件路径列表
cat("搜索CX报告文件...\\n")
file_list <- c()
existing_samples <- c()
existing_treatments <- c()

for (i in seq_along(sample_ids)) {
    sample_dir <- file.path(mapping_dir, sample_ids[i])
    
    if (dir.exists(sample_dir)) {
        cx_files <- list.files(sample_dir, pattern = "*CX_report.txt$", full.names = TRUE)
        
        if (length(cx_files) > 0) {
            file_list <- c(file_list, cx_files[1])
            existing_samples <- c(existing_samples, sample_ids[i])
            existing_treatments <- c(existing_treatments, treatments[i])
        }
    }
}

if (length(file_list) == 0) {
    cat("错误: 没有找到任何CX报告文件\\n")
    quit(status = 1)
}

cat(sprintf("找到 %d 个样品的数据文件\\n", length(file_list)))

# 读取甲基化数据
cat("读取甲基化数据...\\n")
myobj <- methRead(file_list,
                 sample.id = existing_samples,
                 assembly = "custom",
                 treatment = existing_treatments,
                 context = "CpG",
                 mincov = min_coverage)

# 生成描述性统计
pdf(file.path(analysis_dir, "methylkit_analysis", "methylation_stats.pdf"), width = 12, height = 8)
getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
dev.off()

# 过滤数据
cat("过滤低质量数据...\\n")
filtered.myobj <- filterByCoverage(myobj, lo.count = min_coverage, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 合并样品数据
cat("合并样品数据...\\n")
meth <- unite(filtered.myobj, destrand = FALSE)

# 过滤低覆盖位点
meth.filtered <- filterByCoverage(meth, lo.count = min_coverage, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 标准化
cat("标准化数据...\\n")
meth.norm <- normalizeCoverage(meth.filtered)

# 聚类分析
pdf(file.path(analysis_dir, "methylkit_analysis", "clustering_analysis.pdf"), width = 10, height = 8)
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
cat(sprintf("差异甲基化位点统计:\\n"))
cat(sprintf("  总差异位点数: %d\\n", nrow(diff.meth)))
cat(sprintf("  高甲基化位点: %d\\n", nrow(hyper.meth)))
cat(sprintf("  低甲基化位点: %d\\n", nrow(hypo.meth)))

# 保存结果
write.table(diff.meth, file = file.path(analysis_dir, "differential_methylation", "all_differential_sites.txt"), sep = "\\t", quote = FALSE, row.names = FALSE)
write.table(hyper.meth, file = file.path(analysis_dir, "differential_methylation", "hypermethylated_sites.txt"), sep = "\\t", quote = FALSE, row.names = FALSE)
write.table(hypo.meth, file = file.path(analysis_dir, "differential_methylation", "hypomethylated_sites.txt"), sep = "\\t", quote = FALSE, row.names = FALSE)

# 可视化差异甲基化
pdf(file.path(analysis_dir, "methylkit_analysis", "differential_methylation_plots.pdf"), width = 12, height = 8)
plot(myDiff, plot.type = "volcano")
plot(myDiff, plot.type = "scatterplot")
dev.off()

# 生成汇总报告
summary_stats <- data.frame(
    Metric = c("Total CpG sites analyzed", "Differential sites (all)", "Hypermethylated sites", "Hypomethylated sites", "Samples analyzed"),
    Count = c(nrow(meth.norm), nrow(diff.meth), nrow(hyper.meth), nrow(hypo.meth), length(existing_samples))
)

write.table(summary_stats, file = file.path(analysis_dir, "reports", "methylkit_summary.txt"), sep = "\\t", quote = FALSE, row.names = FALSE)

cat("methylKit分析完成！\\n")
"""

        with open(script_path, "w", encoding="utf-8") as f:
            f.write(script_content)

        os.chmod(script_path, 0o755)


class BLASTAnalyzer:
    """BLAST分析器 | BLAST Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_blast_analysis(self) -> bool:
        """运行BLAST定位分析 | Run BLAST localization analysis"""
        if not self.config.target_promoter_fa:
            return True

        self.logger.info("步骤F: 开始目标序列BLAST定位和甲基化分析...")

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

        # 解析BLAST结果
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

            return True
        else:
            self.logger.warning("⚠️  BLAST结果解析失败")
            return False

    def _parse_blast_results(self, blast_result, positions_file) -> bool:
        """解析BLAST结果 | Parse BLAST results"""
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
            self.logger.error(f"BLAST结果解析失败: {e}")
            return False

    def _extract_target_methylation(self, positions_file, target_analysis_dir) -> bool:
        """提取目标区域甲基化信息 | Extract target region methylation"""
        try:
            # 获取所有样品的CX报告文件
            samples = []
            for sample in self.config.sample_groups.keys():
                cx_files = glob.glob(
                    os.path.join(self.config.mapping_dir, sample, "*CX_report.txt")
                )
                if cx_files:
                    samples.append((sample, cx_files[0]))

            if not samples:
                self.logger.warning("未找到甲基化数据文件")
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
            self.logger.error(f"目标区域甲基化信息提取失败: {e}")
            return False
