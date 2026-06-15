"""
DeepBSA绘图模块|DeepBSA Plot Generator Module
使用Python绘制BSA结果图|Plot BSA results using Python
"""

import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as mpl
from matplotlib.gridspec import GridSpec


# 设置中文字体支持|Set Chinese font support
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False


class DeepBSAPlotter:
    """DeepBSA绘图类|DeepBSA Plotter Class"""

    # 默认颜色|Default colors
    COLOR_THRESHOLD = "#1F77B4FF"  # 蓝色|Blue
    COLOR_SMOOTH = "#FF7F0EFF"     # 橙色|Orange
    COLOR_POINT = "#CCCCCC"        # 灰色|Grey (grey80 = #CCCCCC)

    def __init__(self, data_file: Path, output_file: Path, logger=None):
        """初始化绘图器|Initialize plotter

        Args:
            data_file: 绘图数据CSV文件|Plot data CSV file
            output_file: 输出图片文件|Output image file
            logger: 日志器|Logger
        """
        self.data_file = data_file
        self.output_file = output_file
        self.logger = logger

        # 颜色配置|Color configuration
        self.color_threshold = self.COLOR_THRESHOLD
        self.color_smooth = self.COLOR_SMOOTH
        self.color_point = self.COLOR_POINT

        # 图片尺寸|Image size
        self.width = 12
        self.height = 8
        self.dpi = 300

    def set_colors(self, threshold: str = None, smooth: str = None, point: str = None):
        """设置颜色|Set colors

        Args:
            threshold: 阈值线颜色|Threshold line color
            smooth: 平滑曲线颜色|Smooth curve color
            point: 散点颜色|Point color
        """
        if threshold:
            self.color_threshold = threshold
        if smooth:
            self.color_smooth = smooth
        if point:
            self.color_point = point

    def set_size(self, width: float = 12, height: float = 8, dpi: int = 300):
        """设置图片尺寸|Set image size

        Args:
            width: 宽度（英寸）|Width (inches)
            height: 高度（英寸）|Height (inches)
            dpi: 分辨率|Resolution
        """
        self.width = width
        self.height = height
        self.dpi = dpi

    def load_data(self) -> pd.DataFrame:
        """加载数据|Load data

        Returns:
            DataFrame: 绘图数据|Plot data
        """
        if self.logger:
            self.logger.info(f"读取绘图数据|Reading plot data: {self.data_file}")

        df = pd.read_csv(self.data_file)

        # 转换位置为Mb|Convert position to Mb
        df['Position_Mb'] = df['Position'] / 1_000_000

        if self.logger:
            self.logger.info(f"  加载|Loaded {len(df)} 个数据点|data points")
            self.logger.info(f"  方法|Methods: {', '.join(df['Method'].unique())}")
            self.logger.info(f"  染色体|Chromosomes: {', '.join(df['Chromosome'].unique())}")

        return df

    def plot(self, df: pd.DataFrame):
        """绘制BSA结果图|Plot BSA results

        Args:
            df: 绘图数据|Plot data
        """
        if self.logger:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("开始绘制BSA结果图|Plotting BSA Results")
            self.logger.info("=" * 60)

        # 获取唯一的methods和chromosomes
        # Get unique methods and chromosomes
        methods = sorted(df['Method'].unique())
        chromosomes = sorted(df['Chromosome'].unique(), key=lambda x: int(x.replace('Chr', '')) if x.replace('Chr', '').isdigit() else x)

        n_methods = len(methods)
        n_chromosomes = len(chromosomes)

        if self.logger:
            self.logger.info(f"网格布局|Grid layout: {n_methods} 行|rows × {n_chromosomes} 列|cols")

        # 创建图|Create figure
        fig, axes = plt.subplots(
            n_methods,
            n_chromosomes,
            figsize=(self.width, self.height),
            sharex='col',
            sharey='row',
            squeeze=False
        )

        # 为每个method和chromosome绘制子图
        # Plot subplot for each method and chromosome
        for i, method in enumerate(methods):
            for j, chromosome in enumerate(chromosomes):
                ax = axes[i, j]

                # 获取当前method和chromosome的数据
                # Get data for current method and chromosome
                subset = df[(df['Method'] == method) & (df['Chromosome'] == chromosome)]

                if len(subset) == 0:
                    # 无数据时显示空白|Show blank when no data
                    ax.text(0.5, 0.5, 'No Data',
                           transform=ax.transAxes,
                           ha='center', va='center',
                           fontsize=12)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue

                # 获取阈值|Get threshold
                threshold = subset['Threshold'].iloc[0]

                # 绘制散点（原始值）|Plot points (raw values)
                if self.logger:
                    self.logger.debug(f"color_point值|color_point value: {self.color_point!r}")
                ax.scatter(
                    subset['Position_Mb'],
                    subset['Raw_Value'],
                    c=self.color_point,
                    s=1,
                    alpha=0.5,
                    label='Raw Value'
                )

                # 绘制阈值线|Plot threshold line
                ax.axhline(
                    y=threshold,
                    c=self.color_threshold,
                    linestyle='--',
                    linewidth=1,
                    alpha=0.8,
                    label='Threshold'
                )

                # 绘制平滑曲线|Plot smooth curve
                # 按Position排序以避免线条混乱|Sort by Position to avoid messy lines
                subset_sorted = subset.sort_values('Position_Mb')
                ax.plot(
                    subset_sorted['Position_Mb'],
                    subset_sorted['Smooth_Value'],
                    c=self.color_smooth,
                    linewidth=1.5,
                    alpha=0.9,
                    label='Smooth Value'
                )

                # 设置标题和标签|Set title and labels
                if i == 0:
                    ax.set_title(chromosome, fontsize=10, fontweight='bold')

                if j == 0:
                    ax.set_ylabel(method, fontsize=10, fontweight='bold')

                if i == n_methods - 1:
                    ax.set_xlabel('Position (Mb)', fontsize=9)

                # 设置刻度标签字体大小|Set tick label font size
                ax.tick_params(axis='both', which='major', labelsize=8)

                # 添加网格|Add grid
                ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

                # 只在第一行最后一列添加图例|Add legend only in first row last column
                if i == 0 and j == n_chromosomes - 1:
                    ax.legend(
                        loc='upper right',
                        fontsize=7,
                        framealpha=0.9,
                        edgecolor='gray'
                    )

        # 调整布局|Adjust layout
        plt.tight_layout()

        # 保存图片|Save figure
        if self.logger:
            self.logger.info(f"保存图片到|Saving plot to: {self.output_file}")

        # 确保输出目录存在|Ensure output directory exists
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        plt.savefig(
            self.output_file,
            dpi=self.dpi,
            bbox_inches='tight',
            facecolor='white'
        )

        if self.logger:
            self.logger.info(f"图片尺寸|Image size: {self.width}×{self.height} inches, {self.dpi} DPI")
            self.logger.info("绘图完成|Plotting completed")

        plt.close(fig)

    def run(self):
        """运行绘图流程|Run plotting process"""
        try:
            # 加载数据|Load data
            df = self.load_data()

            # 绘图|Plot
            self.plot(df)

            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"绘图失败|Plotting failed: {e}")
            import traceback
            traceback.print_exc()
            return False


def plot_deepbsa_results(
    data_file: Path,
    output_file: Path,
    width: float = 12,
    height: float = 8,
    dpi: int = 300,
    color_threshold: str = None,
    color_smooth: str = None,
    color_point: str = None,
    logger=None
) -> bool:
    """绘制DeepBSA结果的便捷函数|Convenience function to plot DeepBSA results

    Args:
        data_file: 绘图数据CSV文件|Plot data CSV file
        output_file: 输出图片文件|Output image file
        width: 宽度（英寸）|Width (inches)
        height: 高度（英寸）|Height (inches)
        dpi: 分辨率|Resolution
        color_threshold: 阈值线颜色|Threshold line color
        color_smooth: 平滑曲线颜色|Smooth curve color
        color_point: 散点颜色|Point color
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Success or not
    """
    plotter = DeepBSAPlotter(data_file, output_file, logger)

    # 设置尺寸|Set size
    plotter.set_size(width=width, height=height, dpi=dpi)

    # 设置颜色|Set colors
    plotter.set_colors(
        threshold=color_threshold,
        smooth=color_smooth,
        point=color_point
    )

    # 运行绘图|Run plotting
    return plotter.run()
