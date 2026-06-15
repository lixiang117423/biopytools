"""可视化处理器|Visualization Processor"""

import os
import subprocess
from ..utils import CommandRunner


class VisualizationProcessor:
    """可视化处理器|Visualization Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行可视化|Run visualization"""
        self.logger.info("开始生成可视化结果|Starting visualization generation")

        # 1. 曼哈顿图|Manhattan plot
        self._generate_manhattan_plot()

        # 2. QQ图|QQ plot
        self._generate_qq_plot()

        self.logger.info("可视化生成完成|Visualization generation completed")

        return True

    def _generate_manhattan_plot(self):
        """生成曼哈顿图|Generate Manhattan plot"""
        asso_dir = self.config.dirs['association']
        vis_dir = self.config.dirs['visualization']

        if not os.path.exists(asso_dir):
            self.logger.warning("关联结果目录不存在，跳过曼哈顿图|Association directory not found, skip Manhattan plot")
            return

        # 查找关联结果文件|Find association result files
        # 支持多种格式：*.assoc*, *.result*, *.ps (k-mer GWAS)
        # Support multiple formats: *.assoc*, *.result*, *.ps (k-mer GWAS)
        from pathlib import Path
        asso_files = list(Path(asso_dir).glob('*assoc*'))
        asso_files += list(Path(asso_dir).glob('*result*'))
        asso_files += list(Path(asso_dir).glob('*.ps'))

        if not asso_files:
            self.logger.warning("未找到关联结果文件|No association result files found")
            return

        # 使用R脚本生成曼哈顿图（同时生成PDF和PNG）|Use R script to generate Manhattan plot (both PDF and PNG)
        # 支持标准GWAS格式和k-mer GWAS格式（.ps文件）
        # Support standard GWAS format and k-mer GWAS format (.ps files)
        r_script = '''
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript manhattan.r <input_file> <output_prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# 读取数据（尝试检测格式）
# 首先尝试标准GWAS格式（有表头）
data <- tryCatch({
  read.table(input_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
}, error = function(e) {
  # 如果失败，尝试k-mer GWAS格式（无表头：kmer, effect, pval）
  read.table(input_file, header=FALSE, sep="\\t", stringsAsFactors=FALSE)
})

# 检测文件格式并处理
is_kmer_format <- FALSE
if (!("CHR" %in% names(data)) && ncol(data) >= 3) {
  # k-mer GWAS格式：无列名，三列为 [kmer, effect, pval]
  is_kmer_format <- TRUE
  names(data) <- c("kmer", "effect", "pval")
  # 添加索引作为位置
  data$index <- 1:nrow(data)
  data$neglog10p <- -log10(data$pval)
} else if ("CHR" %in% names(data) && "BP" %in% names(data)) {
  # 标准GWAS格式
  if ("P" %in% names(data)) {
    data$neglog10p <- -log10(data$P)
  } else if ("pval" %in% names(data)) {
    names(data)[names(data) == "pval"] <- "P"
    data$neglog10p <- -log10(data$P)
  }
} else {
  stop("Unsupported file format")
}

# 创建图形
if (is_kmer_format) {
  # k-mer GWAS: 按索引的p-value分布图
  p <- ggplot(data, aes(x=index, y=neglog10p)) +
    geom_point(alpha=0.3, size=0.5) +
    theme_bw() +
    labs(x="K-mer Index", y="-log10(P-value)", title="K-mer Association P-value Distribution")
} else {
  # 标准GWAS: 曼哈顿图
  p <- ggplot(data, aes(x=BP, y=neglog10p, color=as.factor(CHR))) +
    geom_point(alpha=0.5, size=1) +
    theme_bw() +
    theme(legend.position="none") +
    labs(x="Chromosomal Position", y="-log10(P-value)", title="Manhattan Plot")
}

# 保存为PDF（矢量图，适合出版）
pdf_file <- paste0(output_prefix, "_manhattan.pdf")
pdf(pdf_file, width=12, height=6)
print(p)
dev.off()
cat("PDF saved:", pdf_file, "\\n")

# 保存为PNG（位图，适合快速查看）
png_file <- paste0(output_prefix, "_manhattan.png")
png(png_file, width=20, height=10, units="in", res=150)
print(p)
dev.off()
cat("PNG saved:", png_file, "\\n")

# 保存高分辨率PNG（适合大数据集）
png_file_hr <- paste0(output_prefix, "_manhattan_hires.png")
png(png_file_hr, width=24, height=12, units="in", res=300)
print(p)
dev.off()
cat("High-res PNG saved:", png_file_hr, "\\n")
'''

        r_script_file = os.path.join(vis_dir, 'manhattan.r')
        with open(r_script_file, 'w') as f:
            f.write(r_script)

        # 为每个结果文件生成图|Generate plot for each result file
        for asso_file in asso_files[:3]:  # 限制处理数量|Limit number of files
            output_prefix = os.path.join(vis_dir, asso_file.stem)

            cmd = ['Rscript', r_script_file, str(asso_file), output_prefix]

            success, _ = self.cmd_runner.run_command(
                cmd,
                description=f"生成曼哈顿图|Generating Manhattan plot: {asso_file.name}",
                check=False
            )

            if success:
                self.logger.info(f"曼哈顿图已保存|Manhattan plots saved: {output_prefix}_manhattan.pdf/png")

    def _generate_qq_plot(self):
        """生成QQ图|Generate QQ plot"""
        asso_dir = self.config.dirs['association']
        vis_dir = self.config.dirs['visualization']

        if not os.path.exists(asso_dir):
            self.logger.warning("关联结果目录不存在，跳过QQ图|Association directory not found, skip QQ plot")
            return

        # 查找关联结果文件|Find association result files
        # 支持多种格式：*.assoc*, *.result*, *.ps (k-mer GWAS)
        # Support multiple formats: *.assoc*, *.result*, *.ps (k-mer GWAS)
        from pathlib import Path
        asso_files = list(Path(asso_dir).glob('*assoc*'))
        asso_files += list(Path(asso_dir).glob('*result*'))
        asso_files += list(Path(asso_dir).glob('*.ps'))

        if not asso_files:
            self.logger.warning("未找到关联结果文件|No association result files found")
            return

        # 使用R脚本生成QQ图（同时生成PDF和PNG）|Use R script to generate QQ plot (both PDF and PNG)
        # 支持标准GWAS格式和k-mer GWAS格式（.ps文件）
        # Support standard GWAS format and k-mer GWAS format (.ps files)
        r_script = '''
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript qqplot.r <input_file> <output_prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# 读取数据（尝试检测格式）
# 首先尝试标准GWAS格式（有表头）
data <- tryCatch({
  read.table(input_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
}, error = function(e) {
  # 如果失败，尝试k-mer GWAS格式（无表头：kmer, effect, pval）
  read.table(input_file, header=FALSE, sep="\\t", stringsAsFactors=FALSE)
})

# 检测文件格式并提取p值列
if (!("P" %in% names(data)) && ncol(data) >= 3) {
  # k-mer GWAS格式：无列名，第三列是pval
  data$pval <- data[, 3]
} else if ("pval" %in% names(data)) {
  data$P <- data$pval
}

# 移除NA值
data <- data[!is.na(data$P), ]

# 计算期望值和观测值
n <- nrow(data)
data$expected <- -log10((rank(data$P, na.last="keep") - 0.5) / n)
data$observed <- -log10(data$P)

# 创建QQ图对象
p <- ggplot(data, aes(x=expected, y=observed)) +
  geom_point(alpha=0.5, size=1) +
  geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
  theme_bw() +
  labs(x="Expected -log10(P-value)", y="Observed -log10(P-value)", title="QQ Plot")

# 保存为PDF（矢量图，适合出版）
pdf_file <- paste0(output_prefix, "_qqplot.pdf")
pdf(pdf_file, width=6, height=6)
print(p)
dev.off()
cat("PDF saved:", pdf_file, "\\n")

# 保存为PNG（位图，适合快速查看）
png_file <- paste0(output_prefix, "_qqplot.png")
png(png_file, width=8, height=8, units="in", res=150)
print(p)
dev.off()
cat("PNG saved:", png_file, "\\n")

# 保存高分辨率PNG（适合大数据集）
png_file_hr <- paste0(output_prefix, "_qqplot_hires.png")
png(png_file_hr, width=10, height=10, units="in", res=300)
print(p)
dev.off()
cat("High-res PNG saved:", png_file_hr, "\\n")
'''

        r_script_file = os.path.join(vis_dir, 'qqplot.r')
        with open(r_script_file, 'w') as f:
            f.write(r_script)

        # 为每个结果文件生成图|Generate plot for each result file
        for asso_file in asso_files[:3]:  # 限制处理数量|Limit number of files
            output_prefix = os.path.join(vis_dir, asso_file.stem)

            cmd = ['Rscript', r_script_file, str(asso_file), output_prefix]

            success, _ = self.cmd_runner.run_command(
                cmd,
                description=f"生成QQ图|Generating QQ plot: {asso_file.name}",
                check=False
            )

            if success:
                self.logger.info(f"QQ图已保存|QQ plots saved: {output_prefix}_qqplot.pdf/png")
