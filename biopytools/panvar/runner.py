"""泛基因组变异分析运行器|Pan-genome Variant Analysis Runner"""

import os
import subprocess
from pathlib import Path
from collections import Counter

import pysam

from .config import PanvarConfig
from .utils import PanvarLogger

VERSION = "1.0.0"


class PanvarRunner:
    """泛基因组变异分析运行器|Pan-genome Variant Analysis Runner"""

    def __init__(self, config: PanvarConfig):
        self.config = config
        log_file = os.path.join(self.config.logs_dir, "panvar.log")
        self.logger_manager = PanvarLogger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

    def run(self) -> bool:
        """运行完整流程|Run full pipeline"""
        self.config.validate()

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(f"Panvar v{VERSION}")
        self.logger.info(f"输入类型|Input type: {self.config.input_type}")
        self.logger.info("=" * 60)

        # Step 1: vg deconstruct (仅GBZ输入)|vg deconstruct (GBZ input only)
        vcf_file = self.config.input_file
        if self.config.input_type == 'gbz':
            vcf_file = self._run_deconstruct()
            if vcf_file is None:
                return False

        # Step 2: 变异分类统计|Variant classification summary
        tsv_output = os.path.join(self.config.summarize_dir, "pangenome_variants_summary.tsv")
        self._summarize_variants(vcf_file, tsv_output)

        # Step 3: 泛基因组增长曲线 + gamma值 (R/ggplot2)|Growth curve + gamma (R/ggplot2)
        png_output = os.path.join(self.config.summarize_dir, "pangenome_growth_curve.png")
        self._plot_growth_curve(tsv_output, png_output)

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("分析完成|Analysis completed")
        self.logger.info("=" * 60)
        return True

    # ============ Step 1: vg deconstruct ============

    def _run_deconstruct(self) -> str:
        """从GBZ图提取变异VCF|Extract variant VCF from GBZ graph"""
        output_vcf = os.path.join(
            self.config.deconstruct_dir,
            f"{Path(self.config.input_file).stem}.vcf"
        )

        if os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0:
            self.logger.info(f"VCF已存在，跳过deconstruct|VCF exists, skipping: {output_vcf}")
            return output_vcf

        self.logger.info(f"执行vg deconstruct|Running vg deconstruct")
        cmd = (
            f"conda run -n {self.config.vg_env} --no-capture-output "
            f"vg deconstruct -t {self.config.threads} -a "
            f"{self.config.ref_path} {self.config.input_file}"
            f" > {output_vcf}"
        )
        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=None)
            if result.returncode != 0:
                self.logger.error(f"vg deconstruct失败|vg deconstruct failed")
                self.logger.error(f"错误信息|Error: {result.stderr[:500]}")
                return None
            size_mb = os.path.getsize(output_vcf) / (1024 * 1024)
            self.logger.info(f"输出VCF|Output VCF: {output_vcf} ({size_mb:.2f} MB)")
            return output_vcf
        except Exception as e:
            self.logger.error(f"vg deconstruct异常|Exception: {e}")
            return None

    # ============ Step 2: 变异分类统计 ============

    def _summarize_variants(self, vcf_path: str, output_path: str):
        """VCF变异分类统计|Classify and summarize VCF variants"""
        self.logger.info(f"变异分类统计|Variant classification: {vcf_path}")

        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)

        if not samples:
            self.logger.warning("VCF中没有样本|No samples in VCF")
            return

        self.logger.info(f"样本数|Sample count: {len(samples)}")

        type_counts = Counter()
        total_records = 0

        with open(output_path, "w") as f:
            header = ["CHROM", "POS", "ID", "TYPE", "LENGTH", "REF",
                      "ALT_COUNT", "ALT_ALLELES", "FREQ_DETAIL"] + samples
            f.write("\t".join(header) + "\n")

            for record in vcf:
                allele_counts = {i: 0 for i in range(len(record.alleles))}
                sample_alleles = []

                for s in samples:
                    gt = record.samples[s]['GT']
                    if gt and gt[0] is not None:
                        allele_counts[gt[0]] += 1
                        sample_alleles.append(record.alleles[gt[0]])
                    else:
                        sample_alleles.append(".")

                alts = record.alleles[1:]
                if not alts:
                    continue

                v_type, v_len = _get_variant_info(record.ref, alts[0])
                type_counts[v_type] += 1
                total_records += 1

                freq_list = []
                for i in range(1, len(alts) + 1):
                    count = allele_counts[i]
                    perc = (count / len(samples)) * 100
                    freq_list.append(f"ALT{i}:{count}({perc:.1f}%)")

                row = [
                    record.chrom, str(record.pos), str(record.id or "."),
                    v_type, str(v_len), record.ref, str(len(alts)),
                    ",".join(alts), "; ".join(freq_list)
                ] + sample_alleles
                f.write("\t".join(row) + "\n")

        self.logger.info(f"变异总数|Total variants: {total_records}")
        self.logger.info(f"SNP: {type_counts['SNP']} | InDel: {type_counts['InDel']} | SV: {type_counts['SV']}")
        self.logger.info(f"结果|Output: {output_path}")

    # ============ Step 3: 泛基因组增长曲线 (R/ggplot2) ============

    def _plot_growth_curve(self, tsv_path: str, png_output: str):
        """用R/ggplot2绘制泛基因组增长曲线并计算gamma值|Plot growth curve with R/ggplot2 and calculate gamma"""
        self.logger.info(f"计算增长曲线与gamma值|Calculating growth curve & gamma: {self.config.n_permutations}次置换")

        r_script = _build_r_script(tsv_path, png_output,
                                   self.config.ref_size_mb,
                                   self.config.n_permutations)

        script_path = os.path.join(self.config.summarize_dir, "plot_growth.R")
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(r_script)

        cmd = f"{self.config.r_path} {script_path}"
        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=None)
            # 解析R输出中的gamma值|Parse gamma value from R output
            for line in result.stdout.strip().split('\n'):
                line = line.strip()
                if line:
                    self.logger.info(f"R输出|R output: {line}")
            if result.returncode != 0:
                self.logger.error(f"R脚本执行失败|R script failed: {result.stderr[:500]}")
        except Exception as e:
            self.logger.error(f"R脚本异常|R script exception: {e}")


# ============ 纯函数 ============

def _get_variant_info(ref: str, alt: str) -> tuple:
    """判断变异类型及长度|Classify variant type and length"""
    ref_len = len(ref)
    alt_len = len(alt)
    abs_ilen = abs(alt_len - ref_len)

    if ref_len == 1 and alt_len == 1:
        return "SNP", 0
    elif abs_ilen <= 50:
        return "InDel", alt_len - ref_len
    else:
        return "SV", alt_len - ref_len


def _build_r_script(tsv_path: str, png_output: str,
                    ref_size_mb: float, n_permutations: int) -> str:
    """构建R绘图脚本|Build R plotting script"""
    # R代码中用{{ }}表示f-string的转义|{{ }} escapes f-string braces
    ref_size_arg = str(ref_size_mb) if ref_size_mb > 0 else "NULL"

    return r'''#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

tsv_path <- "''' + tsv_path + r'''"
output_png <- "''' + png_output + r'''"
ref_size <- ''' + ref_size_arg + r'''
n_perm <- ''' + str(n_permutations) + r'''

# 读取数据|Read data
data <- read.delim(tsv_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
sample_matrix <- data[, 10:ncol(data)]
variant_lengths <- abs(data$LENGTH)
n_samples <- ncol(sample_matrix)

# 自动推断参考基因组大小|Auto-infer ref size
if (is.null(ref_size) || ref_size == 0) {
  chrom_max <- aggregate(POS ~ CHROM, data = data, FUN = max)
  ref_size <- sum(as.numeric(chrom_max$POS)) / 1e6
  if (ref_size < 1.0) ref_size <- 87.0
}
message(sprintf("参考基因组大小|Reference size: %.2f Mb", ref_size))

# 构建0/1存在/缺失矩阵|Build 0/1 presence/absence matrix
is_present <- sample_matrix != "." & sample_matrix != data$REF

# 随机置换计算增长曲线|Random permutation growth curve
growth_results <- matrix(0, nrow = n_perm, ncol = n_samples)

message(sprintf("正在进行%d次随机置换|%d random permutations...", n_perm, n_perm))
for (i in 1:n_perm) {
  shuffled_idx <- sample(1:n_samples)
  current_union <- rep(FALSE, nrow(is_present))
  for (j in 1:n_samples) {
    current_union <- current_union | is_present[, shuffled_idx[j]]
    growth_results[i, j] <- ref_size + (sum(variant_lengths[current_union]) / 1e6)
  }
}

# 统计汇总|Summary statistics
mean_size <- colMeans(growth_results)
sd_size <- apply(growth_results, 2, sd)
summary_df <- data.frame(
  Samples = 1:n_samples,
  Mean_Size = mean_size,
  SD_Size = sd_size
)

# 拟合Heaps' Law: y = a * x^gamma
# log(y) = log(a) + gamma * log(x)
fit <- lm(log(Mean_Size) ~ log(Samples), data = summary_df)
gamma_val <- coef(fit)["log(Samples)"]
a_val <- exp(coef(fit)["(Intercept)"])
r_sq <- summary(fit)$r.squared

message(sprintf("gamma = %.4f", gamma_val))
message(sprintf("R-squared = %.4f", r_sq))
message(sprintf("Heaps' Law: y = %.4f * x^%.4f", a_val, gamma_val))

# 绘图|Plotting
p <- ggplot(summary_df, aes(x = Samples, y = Mean_Size)) +
  geom_ribbon(aes(ymin = Mean_Size - SD_Size, ymax = Mean_Size + SD_Size),
              fill = "skyblue", alpha = 0.3) +
  geom_line(color = "#2c3e50", linewidth = 1) +
  geom_point(color = "#e74c3c", size = 2) +
  scale_x_continuous(breaks = seq(0, n_samples + 10, by = 10), expand = c(0.02, 0)) +
  labs(
    title = "Pan-genome Growth Curve (Heaps' Law)",
    subtitle = sprintf("Reference: %.1f Mb | gamma = %.4f | R-squared = %.4f",
                       ref_size, gamma_val, r_sq),
    x = "Number of Genomes",
    y = "Total Pan-genome Size (Mb)"
  ) +
  stat_function(
    fun = function(x) a_val * x^gamma_val,
    color = "red", linetype = "dashed", linewidth = 1
  ) +
  theme_bw()

ggsave(output_png, p, width = 8, height = 6, dpi = 300)
message(sprintf("图片已保存|Plot saved: %s", output_png))
'''
