"""
单倍型可视化绘图模块|Haplotype Visualization Plotter Module

使用matplotlib绘制单倍型图|Plot haplotype figures using matplotlib
"""

from pathlib import Path
from typing import List, Dict, Optional, Tuple
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端|Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle, FancyBboxPatch

from .parser import Gene, MRNA, CDS, Exon, UTR
from .builder import Haplotype, SNPInfo
from .config import HapTypeConfig


class HaplotypePlotter:
    """单倍型绘图器|Haplotype Plotter"""

    def __init__(self, config: HapTypeConfig, logger):
        """初始化绘图器|Initialize plotter

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 绘图参数|Plotting parameters
        self.fig_width = config.width
        self.fig_height = config.height
        self.dpi = config.dpi

        # 颜色参数|Color parameters
        self.ref_color = config.ref_color
        self.alt_color = config.alt_color
        self.missing_color = config.missing_color

        # 显示参数|Display parameters
        self.show_gene_name = config.show_gene_name
        self.show_position = config.show_position
        self.show_sample_count = config.show_sample_count

    def plot(
        self,
        genes: List[Gene],
        haplotypes: List[Haplotype],
        snp_info_list: List[SNPInfo],
        output_file: str
    ) -> bool:
        """绘制单倍型图|Plot haplotype figure

        Args:
            genes: 基因列表|Gene list
            haplotypes: 单倍型列表|Haplotype list
            snp_info_list: SNP信息列表|SNP information list
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 是否成功|Success
        """
        try:
            self.logger.info("开始绘制单倍型图|Start plotting haplotype figure")

            # 创建图形|Create figure
            fig = plt.figure(
                figsize=(self.fig_width, self.fig_height),
                dpi=self.dpi
            )

            # 计算子图布局|Calculate subplot layout
            n_haplotypes = len(haplotypes)
            gene_height_ratio = 0.3  # 基因模型占30%高度|Gene model takes 30% height
            hap_height_ratio = 0.7  # 单倍型表占70%高度|Haplotype table takes 70% height

            # 创建子图|Create subplots
            gs = fig.add_gridspec(
                2, 1,
                height_ratios=[gene_height_ratio, hap_height_ratio],
                hspace=0.05
            )

            ax_gene = fig.add_subplot(gs[0])
            ax_hap = fig.add_subplot(gs[1])

            # 绘制基因模型|Draw gene model
            self._draw_gene_model(ax_gene, genes)

            # 绘制单倍型表|Draw haplotype table
            self._draw_haplotype_table(ax_hap, haplotypes, snp_info_list)

            # 保存图形|Save figure
            plt.savefig(
                output_file,
                dpi=self.dpi,
                bbox_inches='tight',
                pad_inches=0.1
            )
            plt.close(fig)

            self.logger.info(f"单倍型图已保存到|Haplotype figure saved to: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"绘制单倍型图失败|Failed to plot haplotype figure: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _draw_gene_model(self, ax, genes: List[Gene]):
        """绘制基因模型|Draw gene model

        Args:
            ax: matplotlib轴对象|matplotlib axis object
            genes: 基因列表|Gene list
        """
        ax.set_xlim(self.config.start, self.config.end)
        ax.set_ylim(-1, len(genes))

        # 隐藏y轴刻度|Hide y-axis ticks
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels([g.gene_name for g in genes] if self.show_gene_name else [])

        # 设置x轴|Set x-axis
        ax.set_xlabel(f'Position ({self.config.chrom})', fontsize=10)
        ax.tick_params(axis='x', labelsize=8)

        # 隐藏顶部和右侧边框|Hide top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # 绘制每个基因|Draw each gene
        for i, gene in enumerate(genes):
            self._draw_single_gene(ax, gene, i)

    def _draw_single_gene(self, ax, gene: Gene, y_pos: int):
        """绘制单个基因|Draw single gene

        Args:
            ax: matplotlib轴对象|matplotlib axis object
            gene: 基因对象|Gene object
            y_pos: y轴位置|y-axis position
        """
        # 选择最长的转录本|Select longest transcript
        if not gene.mrnas:
            # 如果没有mRNA，只画基因线|If no mRNA, draw gene line only
            ax.plot(
                [gene.start, gene.end],
                [y_pos, y_pos],
                'k-',
                linewidth=1
            )
            return

        longest_mrna = max(gene.mrnas, key=lambda m: m.end - m.start)

        # 绘制基因线|Draw gene line
        ax.plot(
            [longest_mrna.start, longest_mrna.end],
            [y_pos, y_pos],
            'k-',
            linewidth=1
        )

        # 绘制UTR|Draw UTRs
        if longest_mrna.five_prime_utr:
            utr = longest_mrna.five_prime_utr
            self._draw_utr(ax, utr, y_pos, is_five_prime=True)

        if longest_mrna.three_prime_utr:
            utr = longest_mrna.three_prime_utr
            self._draw_utr(ax, utr, y_pos, is_five_prime=False)

        # 绘制CDS|Draw CDS
        for cds in longest_mrna.cds_list:
            self._draw_cds(ax, cds, y_pos)

        # 绘制内含子线|Draw intron lines
        if len(longest_mrna.cds_list) > 1:
            for i in range(len(longest_mrna.cds_list) - 1):
                current_cds = longest_mrna.cds_list[i]
                next_cds = longest_mrna.cds_list[i + 1]

                # 画内含子线|Draw intron line
                ax.plot(
                    [current_cds.end, next_cds.start],
                    [y_pos, y_pos],
                    'k-',
                    linewidth=0.5
                )

    def _draw_cds(self, ax, cds: CDS, y_pos: int):
        """绘制CDS|Draw CDS

        Args:
            ax: matplotlib轴对象|matplotlib axis object
            cds: CDS对象|CDS object
            y_pos: y轴位置|y-axis position
        """
        height = 0.6
        rect = Rectangle(
            (cds.start, y_pos - height / 2),
            cds.end - cds.start,
            height,
            facecolor='#2E8B57',  # SeaGreen
            edgecolor='black',
            linewidth=0.5
        )
        ax.add_patch(rect)

    def _draw_utr(self, ax, utr: UTR, y_pos: int, is_five_prime: bool):
        """绘制UTR|Draw UTR

        Args:
            ax: matplotlib轴对象|matplotlib axis object
            utr: UTR对象|UTR object
            y_pos: y轴位置|y-axis position
            is_five_prime: 是否为5'UTR|Whether it's 5' UTR
        """
        height = 0.4
        # 绘制UTR为空心矩形|Draw UTR as hollow rectangle
        rect = Rectangle(
            (utr.start, y_pos - height / 2),
            utr.end - utr.start,
            height,
            facecolor='none',
            edgecolor='black',
            linewidth=0.5,
            linestyle='--'
        )
        ax.add_patch(rect)

    def _draw_haplotype_table(
        self,
        ax,
        haplotypes: List[Haplotype],
        snp_info_list: List[SNPInfo]
    ):
        """绘制单倍型表（显示实际碱基）|Draw haplotype table (show actual nucleotides)

        Args:
            ax: matplotlib轴对象|matplotlib axis object
            haplotypes: 单倍型列表|Haplotype list
            snp_info_list: SNP信息列表|SNP information list
        """
        n_haplotypes = len(haplotypes)
        n_snps = len(snp_info_list)

        # 设置坐标轴范围|Set axis limits
        ax.set_xlim(self.config.start - 0.5, self.config.end + 0.5)
        ax.set_ylim(-0.5, n_haplotypes - 0.5)

        # 隐藏坐标轴|Hide axes
        ax.set_xticks([])
        ax.set_yticks([])

        # 隐藏边框|Hide spines
        for spine in ax.spines.values():
            spine.set_visible(False)

        # 绘制SNP位置标记|Draw SNP position markers
        if self.show_position:
            snp_positions = [snp.pos for snp in snp_info_list]
            ax2 = ax.twiny()
            ax2.set_xlim(self.config.start, self.config.end)
            ax2.set_xticks(snp_positions)
            ax2.set_xticklabels([str(pos) for pos in snp_positions], rotation=90, fontsize=6)
            ax2.tick_params(axis='x', which='major', pad=2)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.spines['left'].set_visible(False)

        # 计算每个SNP的显示宽度|Calculate display width for each SNP
        snp_width = (self.config.end - self.config.start) / max(n_snps, 1) * 0.8

        # 绘制单倍型格子（带碱基文字）|Draw haplotype cells (with nucleotide text)
        for i, haplotype in enumerate(haplotypes):
            for j, allele in enumerate(haplotype.allele_pattern):
                if j >= len(snp_info_list):
                    break

                snp = snp_info_list[j]

                # 确定碱基和颜色|Determine nucleotide and color
                if allele == -1 or allele == '.':
                    nucleotide = 'N'
                    color = self.missing_color
                    text_color = 'black'
                elif allele == 0:
                    nucleotide = snp.ref
                    color = self.ref_color
                    text_color = 'black'
                elif allele == 1:
                    nucleotide = snp.alt
                    color = self.alt_color
                    text_color = 'black'
                else:
                    nucleotide = 'N'
                    color = self.missing_color
                    text_color = 'black'

                # 绘制格子背景|Draw cell background
                rect = Rectangle(
                    (snp.pos - snp_width / 2, i - 0.4),
                    snp_width,
                    0.8,
                    facecolor=color,
                    edgecolor='white',
                    linewidth=0.5
                )
                ax.add_patch(rect)

                # 添加碱基文字|Add nucleotide text
                ax.text(
                    snp.pos,
                    i,
                    nucleotide,
                    ha='center',
                    va='center',
                    fontsize=8,
                    color=text_color,
                    weight='bold'
                )

        # 添加标签|Add labels
        for i, haplotype in enumerate(haplotypes):
            # 单倍型ID和样本数|Haplotype ID and sample count
            if self.show_sample_count:
                label = f"{haplotype.haplotype_id} (n={haplotype.sample_count})"
            else:
                label = haplotype.haplotype_id

            ax.text(
                self.config.start - (self.config.end - self.config.start) * 0.02,
                i,
                label,
                ha='right',
                va='center',
                fontsize=8
            )

        # 添加图例|Add legend
        self._add_legend(ax)

    def _add_legend(self, ax):
        """添加图例|Add legend

        Args:
            ax: matplotlib轴对象|matplotlib axis object
        """
        # 创建图例元素|Create legend elements
        ref_patch = Rectangle((0, 0), 1, 1, facecolor=self.ref_color, edgecolor='black')
        alt_patch = Rectangle((0, 0), 1, 1, facecolor=self.alt_color, edgecolor='black')
        missing_patch = Rectangle((0, 0), 1, 1, facecolor=self.missing_color, edgecolor='black')

        # 添加图例|Add legend
        ax.legend(
            [ref_patch, alt_patch, missing_patch],
            ['Ref', 'Alt', 'Missing'],
            loc='upper right',
            fontsize=8,
            frameon=True,
            fancybox=True,
            shadow=True
        )
