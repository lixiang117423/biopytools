"""
泛基因组可视化模块|Pangenome Visualization Module
"""

import pandas as pd
from pathlib import Path
from typing import Dict


class PangenomeVisualizer:
    """泛基因组可视化器|Pangenome Visualizer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def create_rarefaction_plot(self, rarefaction_results: Dict, output_dir: Path):
        """创建稀释曲线图|Create rarefaction curve plot"""
        if not self.config.generate_plots:
            return

        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
        except ImportError:
            self.logger.warning("matplotlib未安装，跳过可视化|matplotlib not installed, skipping visualization")
            return

        self.logger.info("创建稀释曲线图|Creating rarefaction curve plot")

        viz_dir = Path(output_dir) / "05_visualization"
        viz_dir.mkdir(parents=True, exist_ok=True)

        summary_file = Path(output_dir) / "03_rarefaction_analysis" / "rarefaction_curve_summary.tsv"
        if not summary_file.exists():
            self.logger.warning("稀释分析摘要文件不存在，跳过稀释曲线图|Rarefaction summary file not found, skipping plot")
            return

        try:
            df = pd.read_csv(summary_file, sep='\t')

            fig, ax = plt.subplots(1, 1, figsize=(10, 8))

            sample_sizes = df['Sample_Size']
            pan_means = df['Pan_Mean']
            pan_stds = df['Pan_Std']
            core_means = df['Core_Mean']
            core_stds = df['Core_Std']

            ax.plot(sample_sizes, pan_means, 'b-', linewidth=2, marker='o', markersize=4, label='Pan')
            ax.fill_between(sample_sizes, pan_means - pan_stds, pan_means + pan_stds, alpha=0.2, color='blue')

            ax.plot(sample_sizes, core_means, 'r-', linewidth=2, marker='s', markersize=4, label='Core')
            ax.fill_between(sample_sizes, core_means - core_stds, core_means + core_stds, alpha=0.2, color='red')

            ax.set_xlabel('Sample number', fontsize=12)
            ax.set_ylabel('Family number', fontsize=12)
            ax.legend(fontsize=11)
            ax.grid(True, alpha=0.3)

            ax.set_xlim(0, max(sample_sizes) + 1)
            ax.set_ylim(0, max(pan_means) * 1.05)

            plt.tight_layout()

            plot_file = viz_dir / f"pangenome_rarefaction_curve.{self.config.plot_format}"
            plt.savefig(plot_file, dpi=self.config.figure_dpi, bbox_inches='tight')
            plt.close()

            self.logger.info(f"稀释曲线图已保存|Rarefaction curve plot saved: {plot_file}")

        except Exception as e:
            self.logger.error(f"生成稀释曲线图失败|Failed to generate rarefaction curve plot: {e}")

    def create_classification_pie_chart(self, classification_results: Dict, output_dir: Path):
        """创建分类饼图|Create classification pie chart"""
        if not self.config.generate_plots:
            return

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            self.logger.warning("matplotlib未安装，跳过可视化|matplotlib not installed, skipping visualization")
            return

        self.logger.info("创建分类饼图|Creating classification pie chart")

        viz_dir = Path(output_dir) / "05_visualization"
        viz_dir.mkdir(parents=True, exist_ok=True)

        categories = ['core', 'softcore', 'dispensable', 'private']
        counts = [len(classification_results[cat]) for cat in categories]
        labels = ['Core', 'Softcore', 'Dispensable', 'Private']
        colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4']

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        wedges, texts, autotexts = ax.pie(counts, labels=labels, colors=colors, autopct='%1.1f%%',
                                         startangle=90, textprops={'fontsize': 11})

        ax.legend(wedges, [f'{label}: {count:,}' for label, count in zip(labels, counts)],
                 title="Gene Family Categories",
                 loc="center left",
                 bbox_to_anchor=(1, 0, 0.5, 1))

        ax.set_title('Pangenome Gene Family Distribution', fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        plot_file = viz_dir / f"pangenome_classification_pie.{self.config.plot_format}"
        plt.savefig(plot_file, dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

        self.logger.info(f"分类饼图已保存|Classification pie chart saved: {plot_file}")

    def create_frequency_distribution_plot(self, frequency_distributions: Dict, output_dir: Path):
        """创建频率分布图|Create frequency distribution plot"""
        if not self.config.generate_plots:
            return

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            self.logger.warning("matplotlib未安装，跳过可视化|matplotlib not installed, skipping visualization")
            return

        self.logger.info("创建频率分布图|Creating frequency distribution plot")

        viz_dir = Path(output_dir) / "05_visualization"
        viz_dir.mkdir(parents=True, exist_ok=True)

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()

        categories = ['core', 'softcore', 'dispensable', 'private']
        titles = ['Core Gene Families', 'Softcore Gene Families', 'Dispensable Gene Families', 'Private Gene Families']
        colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4']

        for i, (category, title, color) in enumerate(zip(categories, titles, colors)):
            ax = axes[i]
            freq_dist = frequency_distributions.get(category, {})

            if freq_dist:
                frequencies = list(freq_dist.keys())
                counts = list(freq_dist.values())

                ax.bar(frequencies, counts, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
                ax.set_xlabel('Frequency (genomes)', fontsize=10)
                ax.set_ylabel('Number of orthogroups', fontsize=10)
                ax.set_title(title, fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3)

                for freq, count in zip(frequencies, counts):
                    ax.text(freq, count + max(counts) * 0.01, str(count),
                           ha='center', va='bottom', fontsize=8)
            else:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(title, fontsize=11, fontweight='bold')

        plt.tight_layout()

        plot_file = viz_dir / f"frequency_distribution.{self.config.plot_format}"
        plt.savefig(plot_file, dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

        self.logger.info(f"频率分布图已保存|Frequency distribution plot saved: {plot_file}")
