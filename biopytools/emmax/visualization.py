"""
📊 GWAS可视化模块 | GWAS Visualization Module
"""

import os
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")

# 设置中文字体 | Set Chinese font
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "SimHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False


class GWASVisualizer:
    """GWAS可视化器 | GWAS Visualizer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def create_manhattan_plot(self, results_file: Path) -> bool:
        """创建曼哈顿图 | Create Manhattan plot"""
        self.logger.info("📊 创建曼哈顿图 | Creating Manhattan plot")

        try:
            # 读取结果 | Read results
            df = pd.read_csv(results_file, sep="\t")

            # 准备数据 | Prepare data
            df["chrom"] = df["chrom"].astype(str)
            df["bp"] = pd.to_numeric(df["bp"], errors="coerce")
            df["neg_log10_p"] = pd.to_numeric(df["neg_log10_p"], errors="coerce")

            # 去除无效数据 | Remove invalid data
            df = df.dropna(subset=["chrom", "bp", "neg_log10_p"])

            # 计算累积位置 | Calculate cumulative position
            df["cumulative_pos"] = 0
            chrom_lengths = {}
            cumulative_pos = 0

            for chrom in sorted(df["chrom"].unique()):
                chrom_data = df[df["chrom"] == chrom]
                chrom_lengths[chrom] = chrom_data["bp"].max()
                df.loc[df["chrom"] == chrom, "cumulative_pos"] = (
                    chrom_data["bp"] + cumulative_pos
                )
                cumulative_pos += chrom_data["bp"].max() + 1e6  # 添加1Mb间隔

            # 创建图形 | Create figure
            plt.figure(figsize=(16, 8))

            # 绘制每个染色体的点 | Plot points for each chromosome
            colors = [
                "#1f77b4",
                "#ff7f0e",
                "#2ca02c",
                "#d62728",
                "#9467bd",
                "#8c564b",
                "#e377c2",
                "#7f7f7f",
                "#bcbd22",
                "#17becf",
            ]

            for i, chrom in enumerate(sorted(df["chrom"].unique())):
                chrom_data = df[df["chrom"] == chrom]
                color = colors[i % len(colors)]
                plt.scatter(
                    chrom_data["cumulative_pos"],
                    chrom_data["neg_log10_p"],
                    c=color,
                    alpha=0.7,
                    s=20,
                )

            # 添加显著性阈值线 | Add significance threshold line
            threshold = -np.log10(self.config.p_value_threshold)
            plt.axhline(
                y=threshold,
                color="red",
                linestyle="--",
                alpha=0.8,
                label=f"Significance threshold (-log10(p) = {threshold:.1f})",
            )

            # 设置图形属性 | Set figure properties
            plt.xlabel("Chromosome Position")
            plt.ylabel("-log10(p-value)")
            plt.title("Manhattan Plot")
            plt.legend()
            plt.grid(True, alpha=0.3)

            # 设置x轴标签 | Set x-axis labels
            tick_positions = []
            tick_labels = []
            cumulative_pos = 0

            for chrom in sorted(df["chrom"].unique()):
                chrom_data = df[df["chrom"] == chrom]
                center_pos = cumulative_pos + chrom_data["bp"].max() / 2
                tick_positions.append(center_pos)
                tick_labels.append(chrom)
                cumulative_pos += chrom_data["bp"].max() + 1e6

            plt.xticks(tick_positions, tick_labels)

            # 保存图形 | Save figure
            plt.tight_layout()
            plt.savefig(self.config.manhattan_plot, dpi=300, bbox_inches="tight")
            plt.close()

            self.logger.info(
                f"✅ 曼哈顿图已保存 | Manhattan plot saved: {self.config.manhattan_plot}"
            )
            return True

        except Exception as e:
            self.logger.error(
                f"❌ 曼哈顿图创建失败 | Manhattan plot creation failed: {e}"
            )
            return False

    def create_qq_plot(self, results_file: Path) -> bool:
        """创建QQ图 | Create QQ plot"""
        self.logger.info("📊 创建QQ图 | Creating QQ plot")

        try:
            # 读取结果 | Read results
            df = pd.read_csv(results_file, sep="\t")

            # 准备数据 | Prepare data
            p_values = df["p_value"].dropna()
            p_values = p_values[p_values > 0]  # 去除0值

            # 计算期望的p值 | Calculate expected p-values
            n = len(p_values)
            expected_p = np.arange(1, n + 1) / (n + 1)
            observed_p = np.sort(p_values)

            # 转换为-log10 | Convert to -log10
            expected_log = -np.log10(expected_p)
            observed_log = -np.log10(observed_p)

            # 创建图形 | Create figure
            plt.figure(figsize=(8, 8))

            # 绘制QQ图 | Plot QQ plot
            plt.scatter(expected_log, observed_log, alpha=0.6, s=20)

            # 添加对角线 | Add diagonal line
            max_val = max(expected_log.max(), observed_log.max())
            plt.plot(
                [0, max_val],
                [0, max_val],
                "r--",
                alpha=0.8,
                label="Expected under null",
            )

            # 设置图形属性 | Set figure properties
            plt.xlabel("Expected -log10(p-value)")
            plt.ylabel("Observed -log10(p-value)")
            plt.title("QQ Plot")
            plt.legend()
            plt.grid(True, alpha=0.3)

            # 保存图形 | Save figure
            plt.tight_layout()
            plt.savefig(self.config.qq_plot, dpi=300, bbox_inches="tight")
            plt.close()

            self.logger.info(f"✅ QQ图已保存 | QQ plot saved: {self.config.qq_plot}")
            return True

        except Exception as e:
            self.logger.error(f"❌ QQ图创建失败 | QQ plot creation failed: {e}")
            return False

    def extract_significant_snps(self, results_file: Path) -> bool:
        """提取显著SNP | Extract significant SNPs"""
        self.logger.info("🎯 提取显著SNP | Extracting significant SNPs")

        try:
            # 读取结果 | Read results
            df = pd.read_csv(results_file, sep="\t")

            # 筛选显著SNP | Filter significant SNPs
            significant_snps = df[df["p_value"] < self.config.p_value_threshold].copy()
            significant_snps = significant_snps.sort_values("p_value")

            # 添加排名 | Add ranking
            significant_snps["rank"] = range(1, len(significant_snps) + 1)

            # 保存结果 | Save results
            significant_snps.to_csv(self.config.significant_snps, sep="\t", index=False)

            self.logger.info(
                f"✅ 提取到 {len(significant_snps)} 个显著SNP | Extracted {len(significant_snps)} significant SNPs"
            )
            self.logger.info(
                f"📁 显著SNP已保存 | Significant SNPs saved: {self.config.significant_snps}"
            )

            return True

        except Exception as e:
            self.logger.error(
                f"❌ 显著SNP提取失败 | Significant SNP extraction failed: {e}"
            )
            return False
