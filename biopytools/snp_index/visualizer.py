"""
SNP Index可视化模块 | SNP Index Visualization Module
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Optional, Tuple
from .config import SNPIndexConfig
from .sliding_window import SlidingWindowAnalyzer


class SNPIndexVisualizer:
    """SNP Index可视化器 | SNP Index Visualizer"""

    def __init__(self, data: List[Dict], sample_names: List[str],
                 config: Optional[SNPIndexConfig] = None):
        """
        初始化可视化器 | Initialize visualizer

        Args:
            data: 数据列表 | Data list
            sample_names: 样本名称列表 | Sample names list
            config: 配置对象 | Configuration object
        """
        self.data = data
        self.sample_names = sample_names
        self.config = config or SNPIndexConfig()
        self.logger = logging.getLogger(self.__class__.__name__)

        # 设置matplotlib字体 | Set matplotlib font
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
        plt.rcParams['axes.unicode_minus'] = False

        # 提取数据数组 | Extract data arrays
        self._extract_data_arrays()

    def _extract_data_arrays(self) -> None:
        """提取数据数组用于可视化 | Extract data arrays for visualization"""
        if not self.data:
            self.delta_snp_indices = []
            self.snp_indices1 = []
            self.snp_indices2 = []
            self.chromosomes = []
            self.positions = []
            return

        self.delta_snp_indices = [d['delta_snp_index'] for d in self.data]
        self.snp_indices1 = [d['sample1_snp_index'] for d in self.data]
        self.snp_indices2 = [d['sample2_snp_index'] for d in self.data]
        self.chromosomes = [d['chromosome'] for d in self.data]
        self.positions = [d['position'] for d in self.data]

    def create_comprehensive_plot(self, output_file: str, figsize: Tuple[int, int] = (15, 12)) -> bool:
        """
        创建综合分析图 | Create comprehensive analysis plot

        Args:
            output_file: 输出文件路径 | Output file path
            figsize: 图像大小 | Figure size

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            self.logger.info("创建综合分析图 | Creating comprehensive analysis plot")

            fig, axes = plt.subplots(2, 2, figsize=figsize)
            fig.suptitle('SNP Index Analysis Results', fontsize=16, fontweight='bold')

            # 1. ΔSNP index分布 | ΔSNP index distribution
            self._plot_delta_distribution(axes[0, 0])

            # 2. SNP index分布对比 | SNP index distribution comparison
            self._plot_snp_comparison(axes[0, 1])

            # 3. 散点图：Sample1 vs Sample2 | Scatter plot: Sample1 vs Sample2
            self._plot_scatter_comparison(axes[1, 0])

            # 4. 按染色体分布 | Chromosome distribution
            self._plot_chromosome_distribution(axes[1, 1])

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"综合分析图已保存 | Comprehensive analysis plot saved: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"创建综合分析图时出错 | Error creating comprehensive analysis plot: {str(e)}")
            return False

    def create_manhattan_plot(self, output_file: str, figsize: Tuple[int, int] = (12, 6)) -> bool:
        """
        创建曼哈顿图 | Create Manhattan plot

        Args:
            output_file: 输出文件路径 | Output file path
            figsize: 图像大小 | Figure size

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            self.logger.info("创建曼哈顿图 | Creating Manhattan plot")

            plt.figure(figsize=figsize)

            # 按染色体分组 | Group by chromosome
            chrom_data = {}
            for i, chrom in enumerate(self.chromosomes):
                if chrom not in chrom_data:
                    chrom_data[chrom] = {'positions': [], 'deltas': []}
                chrom_data[chrom]['positions'].append(self.positions[i])
                chrom_data[chrom]['deltas'].append(self.delta_snp_indices[i])

            # 排序染色体 | Sort chromosomes
            sorted_chroms = self._sort_chromosomes(list(chrom_data.keys()))

            # 绘制每个染色体 | Plot each chromosome
            colors = plt.cm.tab20(np.linspace(0, 1, len(sorted_chroms)))
            x_tick_positions = []
            x_tick_labels = []
            current_x = 0

            for i, chrom in enumerate(sorted_chroms):
                positions = chrom_data[chrom]['positions']
                deltas = chrom_data[chrom]['deltas']

                # 创建x坐标 | Create x coordinates
                x_coords = [current_x + pos for pos in positions]

                # 绘制散点 | Draw scatter
                plt.scatter(x_coords, deltas, c=[colors[i]], alpha=0.6, s=1, label=chrom)

                # 记录染色体中间位置用于标签 | Record chromosome middle position for label
                x_tick_positions.append(current_x + (max(positions) + min(positions)) / 2)
                x_tick_labels.append(chrom)

                # 更新x偏移 | Update x offset
                if positions:
                    current_x += max(positions) - min(positions) + 1000000  # 1MB gap

            # 添加阈值线 | Add threshold lines
            plt.axhline(y=self.config.extreme_threshold, color='red', linestyle='--', alpha=0.7, linewidth=1)
            plt.axhline(y=-self.config.extreme_threshold, color='red', linestyle='--', alpha=0.7, linewidth=1)
            plt.axhline(y=0, color='gray', linestyle='-', alpha=0.5, linewidth=0.5)

            plt.xlabel('Genomic Position')
            plt.ylabel('Delta SNP Index')
            plt.title('Delta SNP Index Manhattan Plot')

            # 设置x轴标签 | Set x-axis labels
            plt.xticks(x_tick_positions, x_tick_labels, rotation=45, ha='right')

            # 添加网格 | Add grid
            plt.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"曼哈顿图已保存 | Manhattan plot saved: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"创建曼哈顿图时出错 | Error creating Manhattan plot: {str(e)}")
            return False

    def _plot_delta_distribution(self, ax) -> None:
        """绘制ΔSNP index分布图 | Plot ΔSNP index distribution"""
        ax.hist(self.delta_snp_indices, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.5, linewidth=2)
        ax.axvline(x=np.mean(self.delta_snp_indices), color='green', linestyle='--', alpha=0.5, linewidth=2, label='Mean')
        ax.set_xlabel('Delta SNP Index')
        ax.set_ylabel('Frequency')
        ax.set_title('Delta SNP Index Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_snp_comparison(self, ax) -> None:
        """绘制SNP index对比图 | Plot SNP index comparison"""
        if len(self.sample_names) >= 2:
            ax.hist(self.snp_indices1, bins=50, alpha=0.5, label=self.sample_names[0], color='green')
            ax.hist(self.snp_indices2, bins=50, alpha=0.5, label=self.sample_names[1], color='orange')
        else:
            ax.hist(self.snp_indices1, bins=50, alpha=0.5, label='Sample1', color='green')
            ax.hist(self.snp_indices2, bins=50, alpha=0.5, label='Sample2', color='orange')

        ax.set_xlabel('SNP Index')
        ax.set_ylabel('Frequency')
        ax.set_title('SNP Index Distribution Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)

    def _plot_scatter_comparison(self, ax) -> None:
        """绘制散点对比图 | Plot scatter comparison"""
        if len(self.sample_names) >= 2:
            ax.scatter(self.snp_indices2, self.snp_indices1, alpha=0.1, s=1)
            ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=2)  # 对角线 | diagonal
            ax.set_xlabel(f'{self.sample_names[1]} SNP Index')
            ax.set_ylabel(f'{self.sample_names[0]} SNP Index')
            ax.set_title(f'SNP Index: {self.sample_names[0]} vs {self.sample_names[1]}')
        else:
            ax.scatter(self.snp_indices2, self.snp_indices1, alpha=0.1, s=1)
            ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=2)
            ax.set_xlabel('Sample2 SNP Index')
            ax.set_ylabel('Sample1 SNP Index')
            ax.set_title('SNP Index: Sample1 vs Sample2')

        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def _plot_chromosome_distribution(self, ax) -> None:
        """绘制按染色体分布图 | Plot chromosome distribution"""
        # 计算每个染色体的平均ΔSNP index | Calculate mean ΔSNP index for each chromosome
        from collections import defaultdict
        chrom_deltas = defaultdict(list)

        for i, chrom in enumerate(self.chromosomes):
            chrom_deltas[chrom].append(self.delta_snp_indices[i])

        chrom_means = {}
        for chrom, deltas in chrom_deltas.items():
            chrom_means[chrom] = np.mean(deltas)

        # 排序染色体 | Sort chromosomes
        sorted_chroms = self._sort_chromosomes(list(chrom_means.keys()))
        sorted_means = [chrom_means[chrom] for chrom in sorted_chroms]

        # 简化染色体名称 | Simplify chromosome names
        simplified_names = [chrom.split('_')[0] for chrom in sorted_chroms]

        bars = ax.bar(simplified_names, sorted_means, alpha=0.7)

        # 根据值的正负设置颜色 | Set colors based on positive/negative values
        for bar, mean_val in zip(bars, sorted_means):
            if mean_val > 0:
                bar.set_color('red')
            else:
                bar.set_color('blue')

        ax.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
        ax.set_xlabel('Chromosome')
        ax.set_ylabel('Mean Delta SNP Index')
        ax.set_title('Mean Delta SNP Index by Chromosome')

        # 旋转x轴标签 | Rotate x-axis labels
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        ax.grid(True, alpha=0.3)

    def _sort_chromosomes(self, chromosomes: List[str]) -> List[str]:
        """
        排序染色体列表 | Sort chromosome list

        Args:
            chromosomes: 染色体列表 | Chromosome list

        Returns:
            list: 排序后的染色体列表 | Sorted chromosome list
        """
        def chrom_key(chrom):
            # 尝试提取数字部分 | Try to extract numeric part
            parts = chrom.split('_')
            for part in parts:
                try:
                    return int(part)
                except ValueError:
                    continue
            return chrom  # 如果无法提取数字，返回原始字符串 | If cannot extract number, return original string

        try:
            return sorted(chromosomes, key=chrom_key)
        except:
            return sorted(chromosomes)

    def create_distribution_plots(self, output_prefix: str) -> bool:
        """
        创建分布图集合 | Create distribution plot collection

        Args:
            output_prefix: 输出文件前缀 | Output file prefix

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            # ΔSNP index 分布 | ΔSNP index distribution
            plt.figure(figsize=(10, 6))
            plt.hist(self.delta_snp_indices, bins=100, alpha=0.7, color='blue', edgecolor='black')
            plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, linewidth=2, label='Zero line')
            plt.axvline(x=np.mean(self.delta_snp_indices), color='green', linestyle='--', alpha=0.7, linewidth=2, label='Mean')
            plt.axvline(x=np.median(self.delta_snp_indices), color='orange', linestyle='--', alpha=0.7, linewidth=2, label='Median')

            # 添加阈值线 | Add threshold lines
            plt.axvline(x=self.config.extreme_threshold, color='red', linestyle=':', alpha=0.7, linewidth=1, label=f'Threshold (±{self.config.extreme_threshold})')
            plt.axvline(x=-self.config.extreme_threshold, color='red', linestyle=':', alpha=0.7, linewidth=1)

            plt.xlabel('Delta SNP Index')
            plt.ylabel('Frequency')
            plt.title('Delta SNP Index Distribution')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_delta_distribution.png", dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            # SNP index 对比分布 | SNP index comparison distribution
            plt.figure(figsize=(12, 6))
            plt.hist(self.snp_indices1, bins=100, alpha=0.6, label=self.sample_names[0] if len(self.sample_names) >= 2 else 'Sample1', color='green', edgecolor='black')
            plt.hist(self.snp_indices2, bins=100, alpha=0.6, label=self.sample_names[1] if len(self.sample_names) >= 2 else 'Sample2', color='orange', edgecolor='black')
            plt.xlabel('SNP Index')
            plt.ylabel('Frequency')
            plt.title('SNP Index Distribution Comparison')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_snp_comparison.png", dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"分布图已保存 | Distribution plots saved with prefix: {output_prefix}")
            return True

        except Exception as e:
            self.logger.error(f"创建分布图时出错 | Error creating distribution plots: {str(e)}")
            return False

    def create_correlation_plot(self, output_file: str) -> bool:
        """
        创建相关性散点图 | Create correlation scatter plot

        Args:
            output_file: 输出文件路径 | Output file path

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            plt.figure(figsize=(8, 8))

            # 创建密度散点图 | Create density scatter plot
            from matplotlib.colors import LogNorm

            if len(self.sample_names) >= 2:
                plt.hist2d(self.snp_indices2, self.snp_indices1, bins=100, cmap='YlOrRd', norm=LogNorm())
                plt.colorbar(label='Point Density')
                plt.plot([0, 1], [0, 1], 'b--', alpha=0.7, linewidth=2, label='y = x')
                plt.xlabel(f'{self.sample_names[1]} SNP Index')
                plt.ylabel(f'{self.sample_names[0]} SNP Index')
                plt.title(f'SNP Index Correlation: {self.sample_names[0]} vs {self.sample_names[1]}')

                # 计算相关系数 | Calculate correlation coefficient
                correlation = np.corrcoef(self.snp_indices2, self.snp_indices1)[0, 1]
                plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}',
                        transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            else:
                plt.hist2d(self.snp_indices2, self.snp_indices1, bins=100, cmap='YlOrRd', norm=LogNorm())
                plt.colorbar(label='Point Density')
                plt.plot([0, 1], [0, 1], 'b--', alpha=0.7, linewidth=2, label='y = x')
                plt.xlabel('Sample2 SNP Index')
                plt.ylabel('Sample1 SNP Index')
                plt.title('SNP Index Correlation: Sample1 vs Sample2')

                correlation = np.corrcoef(self.snp_indices2, self.snp_indices1)[0, 1]
                plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}',
                        transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"相关性图已保存 | Correlation plot saved: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"创建相关性图时出错 | Error creating correlation plot: {str(e)}")
            return False

    def create_sliding_window_plot(self, output_file: str,
                                  window_size: int = 1000000,
                                  step_size: int = 100000,
                                  confidence_level: float = 0.95,
                                  show_confidence: bool = True,
                                  show_threshold: bool = True) -> bool:
        """
        创建滑动窗口折线图 | Create sliding window line plot

        Args:
            output_file: 输出文件路径 | Output file path
            window_size: 窗口大小(bp) | Window size in bp
            step_size: 步长(bp) | Step size in bp
            confidence_level: 置信水平 | Confidence level
            show_confidence: 是否显示置信区间 | Whether to show confidence intervals
            show_threshold: 是否显示阈值线 | Whether to show threshold lines

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            self.logger.info("创建滑动窗口折线图 | Creating sliding window line plot")

            # 设置滑动窗口参数 | Set sliding window parameters
            if hasattr(self.config, 'window_size'):
                window_size = self.config.window_size
            if hasattr(self.config, 'step_size'):
                step_size = self.config.step_size

            # 创建滑动窗口分析器 | Create sliding window analyzer
            self.config.window_size = window_size
            self.config.step_size = step_size
            analyzer = SlidingWindowAnalyzer(self.data, self.config)

            # 创建窗口 | Create windows
            windows = analyzer.create_sliding_windows()

            if not windows:
                self.logger.warning("没有足够的窗口数据 | Insufficient window data")
                return False

            # 计算置信区间 | Calculate confidence intervals
            confidence_intervals = None
            if show_confidence:
                try:
                    from scipy import stats
                    confidence_intervals = analyzer.calculate_confidence_intervals(confidence_level)
                except ImportError:
                    self.logger.warning("scipy未安装，跳过置信区间计算 | scipy not installed, skipping confidence intervals")
                    show_confidence = False

            # 识别候选区域 | Identify candidate regions
            candidate_regions = analyzer.identify_candidate_regions()

            # 创建图形 | Create figure
            fig, ax = plt.subplots(figsize=(14, 8))

            # 按染色体分组绘制 | Draw by chromosome groups
            self._plot_chromosome_lines(ax, windows, confidence_intervals, show_confidence)

            # 添加阈值线 | Add threshold lines
            if show_threshold:
                threshold = self.config.region_threshold
                ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.7, linewidth=2,
                          label=f'Threshold (±{threshold})')
                ax.axhline(y=-threshold, color='red', linestyle='--', alpha=0.7, linewidth=2)

            # 添加置信区间 | Add confidence intervals
            if show_confidence and confidence_intervals:
                ci_lower = confidence_intervals.get('ci_lower_percentile', 0)
                ci_upper = confidence_intervals.get('ci_upper_percentile', 0)
                ax.axhline(y=ci_lower, color='blue', linestyle=':', alpha=0.7, linewidth=1,
                          label=f'{confidence_level*100:.0f}% CI')
                ax.axhline(y=ci_upper, color='blue', linestyle=':', alpha=0.7, linewidth=1)

            # 标记候选区域 | Mark candidate regions
            self._highlight_candidate_regions(ax, candidate_regions)

            # 设置图形属性 | Set figure properties
            ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
            ax.set_ylabel('Mean Delta SNP index', fontsize=12)
            ax.set_title('Sliding Window Delta SNP index Analysis', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right')

            # 设置y轴范围 | Set y-axis limits
            y_min = min([w['mean_delta'] for w in windows])
            y_max = max([w['mean_delta'] for w in windows])
            y_range = y_max - y_min
            ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"滑动窗口折线图已保存 | Sliding window line plot saved: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"创建滑动窗口折线图时出错 | Error creating sliding window line plot: {str(e)}")
            return False

    def _plot_chromosome_lines(self, ax, windows: List[Dict],
                              confidence_intervals: Optional[Dict],
                              show_confidence: bool) -> None:
        """
        按染色体绘制折线 | Draw lines by chromosome

        Args:
            ax: matplotlib轴对象 | matplotlib axis object
            windows: 窗口数据 | Window data
            confidence_intervals: 置信区间 | Confidence intervals
            show_confidence: 是否显示置信区间 | Whether to show confidence intervals
        """
        # 按染色体分组 | Group by chromosome
        chrom_windows = {}
        for window in windows:
            chrom = window['chromosome']
            if chrom not in chrom_windows:
                chrom_windows[chrom] = []
            chrom_windows[chrom].append(window)

        # 获取所有染色体并排序 | Get and sort all chromosomes
        chromosomes = self._sort_chromosomes(list(chrom_windows.keys()))

        # 颜色映射 | Color mapping
        colors = plt.cm.tab10(np.linspace(0, 1, len(chromosomes)))

        # 计算基因组累积位置 | Calculate cumulative genomic positions
        cumulative_positions = {}
        cumulative_length = 0

        for i, chrom in enumerate(chromosomes):
            chrom_windows[chrom].sort(key=lambda x: x['start'])

            # 为该染色体的窗口计算累积位置 | Calculate cumulative positions for this chromosome's windows
            for window in chrom_windows[chrom]:
                # 转换为Mb | Convert to Mb
                window['cumulative_pos_mb'] = (cumulative_length + window['center']) / 1000000

            # 更新累积长度 | Update cumulative length
            if chrom_windows[chrom]:
                max_pos = max(w['end'] for w in chrom_windows[chrom])
                cumulative_length += max_pos + 1000000  # 添加染色体间间隔 | Add inter-chromosome gap

        # 绘制每个染色体的折线 | Draw line for each chromosome
        for i, chrom in enumerate(chromosomes):
            windows_list = chrom_windows[chrom]
            if not windows_list:
                continue

            positions = [w['cumulative_pos_mb'] for w in windows_list]
            deltas = [w['mean_delta'] for w in windows_list]

            # 绘制主折线 | Draw main line
            ax.plot(positions, deltas, color=colors[i], linewidth=2,
                   alpha=0.8, label=chrom, marker='o', markersize=3, markevery=max(1, len(positions)//50))

            # 绘制置信区间阴影 | Draw confidence interval shadow
            if show_confidence and confidence_intervals:
                std_multiplier = 1.96 if confidence_intervals.get('confidence_level') == 0.95 else 2.576
                std_err = confidence_intervals.get('std', 0) / np.sqrt(len(windows_list))

                upper_bound = [d + std_multiplier * std_err for d in deltas]
                lower_bound = [d - std_multiplier * std_err for d in deltas]

                ax.fill_between(positions, lower_bound, upper_bound,
                              color=colors[i], alpha=0.2)

        # 添加染色体分隔线 | Add chromosome separators
        current_pos = 0
        for i, chrom in enumerate(chromosomes[:-1]):
            if chrom in chrom_windows and chrom_windows[chrom]:
                max_pos = max(w['end'] for w in chrom_windows[chrom])
                current_pos += (max_pos + 1000000) / 1000000  # 转换为Mb | Convert to Mb
                ax.axvline(x=current_pos, color='gray', linestyle=':', alpha=0.5)

    def _highlight_candidate_regions(self, ax, candidate_regions: List[Dict]) -> None:
        """
        高亮显示候选区域 | Highlight candidate regions

        Args:
            ax: matplotlib轴对象 | matplotlib axis object
            candidate_regions: 候选区域列表 | Candidate region list
        """
        for region in candidate_regions:
            # 计算区域的累积位置 | Calculate cumulative position of region
            # 这里需要根据实际的染色体顺序来计算 | This needs to be calculated based on actual chromosome order
            center_pos = (region['start'] + region['end']) / 2 / 1000000  # 转换为Mb | Convert to Mb

            # 简化处理：直接使用center position | Simplified: use center position directly
            ax.scatter(center_pos, region['mean_delta'],
                      color='red', s=100, alpha=0.7, marker='*',
                      edgecolors='darkred', linewidth=1,
                      label='Candidate Region' if region == candidate_regions[0] else '')

    def create_multi_chromosome_sliding_plot(self, output_file: str,
                                            chromosomes: Optional[List[str]] = None,
                                            window_size: int = 1000000,
                                            step_size: int = 100000) -> bool:
        """
        创建多染色体分离的滑动窗口图 | Create multi-chromosome separated sliding window plot

        Args:
            output_file: 输出文件路径 | Output file path
            chromosomes: 要显示的染色体列表 | Chromosome list to display
            window_size: 窗口大小 | Window size
            step_size: 步长 | Step size

        Returns:
            bool: 创建是否成功 | Whether creation succeeded
        """
        try:
            self.logger.info("创建多染色体分离的滑动窗口图 | Creating multi-chromosome separated sliding window plot")

            # 设置滑动窗口参数 | Set sliding window parameters
            self.config.window_size = window_size
            self.config.step_size = step_size
            analyzer = SlidingWindowAnalyzer(self.data, self.config)
            windows = analyzer.create_sliding_windows()

            # 按染色体分组 | Group by chromosome
            chrom_windows = {}
            for window in windows:
                chrom = window['chromosome']
                if chrom not in chrom_windows:
                    chrom_windows[chrom] = []
                chrom_windows[chrom].append(window)

            # 如果没有指定染色体，显示前几个 | If no chromosomes specified, show first few
            if chromosomes is None:
                chromosomes = self._sort_chromosomes(list(chrom_windows.keys()))[:12]  # 最多12个染色体

            # 创建子图 | Create subplots
            n_chroms = len(chromosomes)
            cols = 3
            rows = (n_chroms + cols - 1) // cols

            fig, axes = plt.subplots(rows, cols, figsize=(15, 3 * rows))
            if n_chroms == 1:
                axes = [axes]
            elif rows == 1:
                axes = axes.reshape(1, -1)
            elif cols == 1:
                axes = axes.reshape(-1, 1)

            for i, chrom in enumerate(chromosomes):
                row = i // cols
                col = i % cols

                if rows == 1:
                    ax = axes[col] if cols > 1 else axes[0]
                elif cols == 1:
                    ax = axes[row]
                else:
                    ax = axes[row, col]

                if chrom in chrom_windows:
                    chrom_data = chrom_windows[chrom]
                    chrom_data.sort(key=lambda x: x['start'])

                    positions_mb = [w['center'] / 1000000 for w in chrom_data]
                    deltas = [w['mean_delta'] for w in chrom_data]

                    ax.plot(positions_mb, deltas, linewidth=2, color='steelblue', marker='o', markersize=2)
                    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.5)

                    # 添加阈值线 | Add threshold lines
                    threshold = self.config.region_threshold
                    ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.7)
                    ax.axhline(y=-threshold, color='red', linestyle='--', alpha=0.7)

                    ax.set_title(f'Chromosome {chrom}', fontsize=10)
                    ax.set_xlabel('Position (Mb)', fontsize=8)
                    ax.set_ylabel('Delta SNP index', fontsize=8)
                    ax.grid(True, alpha=0.3)

                # 隐藏多余的子图 | Hide extra subplots
                if i >= n_chroms:
                    ax.set_visible(False)

            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()

            self.logger.info(f"多染色体滑动窗口图已保存 | Multi-chromosome sliding window plot saved: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"创建多染色体滑动窗口图时出错 | Error creating multi-chromosome sliding window plot: {str(e)}")
            return False