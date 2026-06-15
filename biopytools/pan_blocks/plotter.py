"""泛基因组Block构建 - 可视化模块|Pan-Blocks Construction - Visualization Module"""

import os
from pathlib import Path
from typing import List, Dict, Tuple
import logging

from .config import PanBlocksConfig
from .utils import parse_coords_file


class PanBlocksPlotter:
    """泛基因组Block可视化|Pan-Genome Block Plotter"""

    def __init__(self, config: PanBlocksConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.coords_dir = Path(config.coords_dir)
        self.blocks_dir = Path(config.blocks_dir)
        self.plots_dir = Path(config.plots_dir)

    def plot_all_chromosomes(self):
        """绘制所有染色体Pan-Blocks图|Plot pan-blocks for all chromosomes"""
        import matplotlib
        matplotlib.use('Agg')

        chromosomes = self.config.get_target_chromosomes()
        self.logger.info(f"开始绘制|Starting plotting: {len(chromosomes)} 条染色体|chromosomes")

        success = 0
        for chrom in chromosomes:
            try:
                self.plot_chromosome(chrom)
                success += 1
            except Exception as e:
                self.logger.error(f"绘图失败|Plot failed for {chrom}: {e}")

        self.logger.info(f"绘图完成|Plotting completed: {success}/{len(chromosomes)}")
        return success == len(chromosomes)

    def plot_chromosome(self, chrom: str):
        """绘制单条染色体Pan-Blocks图|Plot pan-blocks for a single chromosome"""
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import matplotlib.lines as mlines

        genome_order = self.config.genome_order_list
        chrlen_data = self.config.chrlen_data

        # 读取 Pan-Blocks 数据（含 contributor 信息）
        genome_blocks = {}
        pan_blocks_file = self.blocks_dir / f"{chrom}.pan_blocks.bed"
        if not pan_blocks_file.exists():
            self.logger.warning(f"Pan-Blocks文件不存在|Pan-blocks file not found: {chrom}, 跳过|skipping")
            return

        with open(pan_blocks_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('Chr'):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                genome_name = parts[3]
                if genome_name not in genome_blocks:
                    genome_blocks[genome_name] = []
                genome_blocks[genome_name].append((int(parts[1]), int(parts[2]), parts[3]))

        # 颜色映射
        colors = self._generate_colors(len(genome_order))
        genome_color = {name: colors[i] for i, name in enumerate(genome_order)}

        # 计算最长染色体（用于确定X轴范围）
        max_chr_len = max(
            chrlen_data.get(g, {}).get(chrom, 0) for g in genome_order
        )

        # 图形参数
        num_genomes = len(genome_order)
        bar_height = 0.6
        spacing = max(0.8, 300.0 / num_genomes) if num_genomes > 20 else 1.0
        fig_width = self.config.plot_width
        fig_height = max(self.config.plot_height, num_genomes * spacing * 0.08)

        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

        for i, genome in enumerate(genome_order):
            y_center = i * spacing
            chr_len = chrlen_data.get(genome, {}).get(chrom, max_chr_len)

            # 背景条
            ax.barh(y_center, chr_len, height=bar_height,
                    color='#f6f5ec', edgecolor='black', linewidth=0.5)

            # Pan-Blocks 彩色块
            blocks = genome_blocks.get(genome, [])
            for start, end, contributor in blocks:
                color = genome_color.get(contributor, '#cccccc')
                rect = patches.Rectangle(
                    (start, y_center - bar_height / 2), end - start, bar_height,
                    facecolor=color, edgecolor='none', alpha=0.85
                )
                ax.add_patch(rect)

            # 基因组标签
            ax.text(-max_chr_len * 0.01, y_center, genome,
                    ha='right', va='center', fontsize=max(5, 8 - num_genomes * 0.1))

            # 长度标签
            ax.text(-max_chr_len * 0.005, y_center - bar_height / 2 - spacing * 0.15,
                    f"{chr_len / 1e6:.1f}Mb",
                    ha='right', va='center', fontsize=max(4, 6 - num_genomes * 0.05))

        # 绘制共线性连线（相邻基因组之间）
        for i in range(len(genome_order) - 1):
            genome_a = genome_order[i]
            genome_b = genome_order[i + 1]
            self._draw_synteny_links(
                ax, genome_a, genome_b, i, spacing, bar_height,
                chrom, max_chr_len, genome_color
            )

        # 格式化
        ax.set_xlim(0, max_chr_len)
        ax.set_ylim(-spacing * 0.5, num_genomes * spacing - spacing * 0.5)
        ax.set_xlabel(f'{chrom} position (Mb)', fontsize=10)
        ax.set_title(f'Pan-genome blocks: {chrom}', fontsize=12)

        # X轴 Mb 标签
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x / 1e6:.0f}'))
        ax.set_yticks([])

        # 图例
        legend_handles = [
            mlines.Line2D([], [], color=colors[i], label=genome_order[i])
            for i in range(len(genome_order))
        ]
        ncol = max(1, len(genome_order) // 8 + 1)
        ax.legend(handles=legend_handles, loc='upper right',
                  fontsize=max(5, 7 - num_genomes * 0.1), ncol=ncol)

        plt.tight_layout()

        # 保存
        fmt = self.config.plot_format
        output_file = self.plots_dir / f"{chrom}.{fmt}"
        plt.savefig(str(output_file), format=fmt, dpi=300)
        plt.close()
        self.logger.info(f"图表已保存|Plot saved: {output_file}")

    def _draw_synteny_links(self, ax, genome_a: str, genome_b: str, idx_a: int,
                            spacing: float, bar_height: float, chrom: str,
                            max_chr_len: int, genome_color: dict):
        """绘制两基因组间的共线性连线|Draw synteny links between two genomes"""
        coords_file = self.coords_dir / f"{genome_a}.vs.{genome_b}.filtered.coords"
        if not coords_file.exists():
            return

        alignments = parse_coords_file(str(coords_file))

        y_bottom = idx_a * spacing + bar_height / 2
        y_top = (idx_a + 1) * spacing - bar_height / 2

        for aln in alignments:
            if aln['ref_chr'] != chrom or aln['qry_chr'] != chrom:
                continue

            ref_mid = (aln['ref_start'] + aln['ref_end']) / 2
            qry_mid = (aln['qry_start'] + aln['qry_end']) / 2

            ref_forward = aln['ref_start'] < aln['ref_end']
            qry_forward = aln['qry_start'] < aln['qry_end']

            if ref_forward == qry_forward:
                color, alpha = '#d3d7d4', 0.3
            else:
                color, alpha = '#00EE00', 0.4

            ax.plot([ref_mid, qry_mid], [y_bottom, y_top],
                    color=color, alpha=alpha, linewidth=0.3)

    def _generate_colors(self, n: int) -> List[str]:
        """生成n个可区分的颜色|Generate n distinguishable colors"""
        default_colors = [
            "#F8766D", "#00BFC4", "#B79F00", "#619CFF", "#00BA38",
            "#F564E3", "#694d9f", "#f36c21", "#007d65", "#fdb933",
            "#ADFF2F", "#BC8F8F", "#A0522D", "#DEB887", "#BDB76B",
            "#F08080", "#5F9EA0", "#EE82EE", "#9370DB", "#DB7093",
        ]
        if n <= len(default_colors):
            return default_colors[:n]
        # 超过预设颜色数时重复
        return default_colors * ((n // len(default_colors)) + 1)
