"""
图表生成模块 | Plot Generation Module
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

class PlotGenerator:
    """图表生成器 | Plot Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def create_ld_heatmap(self, ld_matrix: np.ndarray, positions: np.ndarray):
        """创建LD热图 | Create LD heatmap"""
        if self.config.verbose:
            self.logger.info("正在生成热图 | Generating heatmap...")
        
        # 应用LD阈值 | Apply LD threshold
        if self.config.ld_threshold > 0:
            ld_matrix[ld_matrix < self.config.ld_threshold] = 0
        
        # 创建子图布局 | Create subplot layout
        fig = plt.figure(figsize=self.config.figsize, dpi=self.config.dpi)
        
        # 创建两个子图：上半部分线条图，下半部分三角形热图
        # Create two subplots: upper line plot, lower triangular heatmap
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 2], hspace=0.05)
        
        # 上半部分：绘制竖直线条图 | Upper part: Draw vertical line plot
        ax_lines = fig.add_subplot(gs[0])
        self._draw_line_plot(ax_lines, ld_matrix, positions)
        
        # 下半部分：绘制三角形热图 | Lower part: Draw triangular heatmap
        ax_heatmap = fig.add_subplot(gs[1])
        self._draw_triangular_heatmap(ax_heatmap, ld_matrix, positions)
        
        # 设置总标题 | Set overall title
        if self.config.title:
            fig.suptitle(self.config.title, fontsize=16, y=0.95)
        else:
            fig.suptitle('Linkage Disequilibrium Analysis', fontsize=16, y=0.95)
        
        # 使用手动调整布局 | Use manual layout adjustment
        plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1, hspace=0.05)
        
        # 保存图片 | Save figure
        plt.savefig(self.config.output_file, dpi=self.config.dpi, bbox_inches='tight')
        if self.config.verbose:
            self.logger.info(f"热图已保存到 | Heatmap saved to: {self.config.output_file}")
        
        plt.close()
    
    def _draw_line_plot(self, ax, ld_matrix: np.ndarray, positions: np.ndarray):
        """绘制上半部分的竖直线条图 | Draw upper vertical line plot"""
        n_variants = len(positions)
        
        # 设置颜色映射 | Set colormap
        cmap = plt.cm.get_cmap(self.config.colormap)
        
        # 为每个SNP位置绘制竖直线条 | Draw vertical lines for each SNP position
        for i in range(n_variants):
            # 计算该位置的平均LD值作为线条高度 | Calculate average LD as line height
            avg_ld = np.nanmean([ld_matrix[i, j] for j in range(n_variants) if j != i])
            
            # 绘制竖直线条 | Draw vertical line
            color = cmap(avg_ld)
            ax.plot([i, i], [0, avg_ld], color=color, linewidth=2, alpha=0.8)
            
            # 在顶部添加小圆点 | Add dot at top
            ax.scatter(i, avg_ld, color=color, s=20, zorder=10)
        
        # 设置坐标轴 | Set axes
        ax.set_xlim(-0.5, n_variants-0.5)
        ax.set_ylim(0, 1)
        
        # 隐藏x轴标签 | Hide x-axis labels
        ax.set_xticks([])
        ax.set_ylabel('Average LD', fontsize=10)
        
        # 美化边框 | Beautify spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        # 添加网格 | Add grid
        ax.grid(True, axis='y', alpha=0.3)
    
    def _draw_triangular_heatmap(self, ax, ld_matrix: np.ndarray, positions: np.ndarray):
        """绘制下半部分的三角形热图 | Draw lower triangular heatmap"""
        n_variants = len(positions)
        
        # 设置颜色映射 | Set colormap
        cmap = plt.cm.get_cmap(self.config.colormap)
        
        # 绘制三角形热图 | Draw triangular heatmap
        for i in range(n_variants):
            for j in range(i, n_variants):
                r2_value = ld_matrix[i, j]
                
                # 跳过低于阈值的值 | Skip values below threshold
                if r2_value < self.config.ld_threshold:
                    continue
                
                # 计算三角形位置 | Calculate triangle position
                # 将矩形坐标转换为三角形坐标
                x_center = (i + j) / 2
                y_center = (j - i) / 2
                
                # 绘制小正方形 | Draw small square
                color = cmap(r2_value)
                rect = plt.Rectangle((x_center - 0.4, y_center - 0.4), 0.8, 0.8,
                                   facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(rect)
        
        # 设置坐标轴 | Set axes
        ax.set_xlim(-1, n_variants)
        ax.set_ylim(-1, n_variants/2 + 1)
        ax.set_aspect('equal')
        
        # 隐藏坐标轴 | Hide axes
        ax.set_xticks([])
        ax.set_yticks([])
        
        # 移除边框 | Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # 添加颜色条 | Add colorbar
        # 创建一个假的mappable对象用于颜色条
        import matplotlib.cm as cm
        mappable = cm.ScalarMappable(cmap=cmap)
        mappable.set_array([0, 1])
        cbar = plt.colorbar(mappable, ax=ax, shrink=0.6, aspect=20, pad=0.02)
        cbar.set_label('r²', rotation=0, labelpad=15)
        
        # 添加位置标签 | Add position labels
        if n_variants <= 20:
            step = max(1, n_variants // 10)
            tick_positions = range(0, n_variants, step)
            tick_labels = [f"{positions[i]:,}" for i in tick_positions]
        else:
            n_ticks = min(8, n_variants)
            tick_indices = np.linspace(0, n_variants-1, n_ticks, dtype=int)
            tick_positions = tick_indices
            tick_labels = [f"{positions[i]:,}" for i in tick_indices]
        
        # 在底部添加标签 | Add labels at bottom
        for i, (pos, label) in enumerate(zip(tick_positions, tick_labels)):
            ax.text(pos, -0.8, label, rotation=45, ha='right', va='top', fontsize=8)
    
    def _draw_arc_plot(self, ax, ld_matrix: np.ndarray, positions: np.ndarray):
        """绘制上半部分的弧线连接图 | Draw upper arc connection plot"""
        n_variants = len(positions)
        
        # 设置颜色映射 | Set colormap
        cmap = plt.cm.get_cmap(self.config.colormap)
        
        # 计算最大弧线高度 | Calculate maximum arc height
        max_height = n_variants / 6  # 调整比例使弧线更合适
        
        # 绘制弧线 | Draw arcs
        for i in range(n_variants):
            for j in range(i+1, n_variants):
                r2_value = ld_matrix[i, j]
                
                # 跳过低于阈值的值 | Skip values below threshold
                if r2_value < self.config.ld_threshold:
                    continue
                
                # 计算弧线参数 | Calculate arc parameters
                x1, x2 = i, j
                xc = (x1 + x2) / 2  # 弧线中心x坐标 | Arc center x coordinate
                width = x2 - x1     # 弧线宽度 | Arc width
                height = min(width / 2, max_height)  # 弧线高度，限制最大高度
                
                # 根据r²值设置颜色和透明度 | Set color and alpha based on r² value
                color = cmap(r2_value)
                alpha = 0.3 + 0.7 * r2_value  # 透明度与r²值成正比 | Alpha proportional to r²
                linewidth = 0.5 + r2_value * 2  # 线宽与r²值成正比
                
                # 绘制弧线 | Draw arc
                arc = patches.Arc((xc, 0), width, height, 
                                angle=0, theta1=0, theta2=180,
                                color=color, alpha=alpha, linewidth=linewidth)
                ax.add_patch(arc)
        
        # 在底部绘制位置标记 | Draw position markers at bottom
        ax.scatter(range(n_variants), [0]*n_variants, 
                  c='black', s=15, zorder=10, alpha=0.7)
        
        # 设置坐标轴 | Set axes
        ax.set_xlim(-0.5, n_variants-0.5)
        ax.set_ylim(0, max_height)
        
        # 移除 set_aspect('equal') 以避免布局问题
        # Remove set_aspect('equal') to avoid layout issues
        
        # 添加位置标签 | Add position labels
        if n_variants <= 20:  # 只在SNP数量较少时显示所有标签 | Only show all labels when SNP count is small
            tick_positions = range(n_variants)
            tick_labels = [f"{pos:,}" for pos in positions]
            rotation = 45 if n_variants > 10 else 0
        else:
            # 选择代表性位置显示标签 | Select representative positions for labels
            n_ticks = min(8, n_variants)  # 减少标签数量
            tick_indices = np.linspace(0, n_variants-1, n_ticks, dtype=int)
            tick_positions = tick_indices
            tick_labels = [f"{positions[i]:,}" for i in tick_indices]
            rotation = 45
        
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=rotation, ha='right', fontsize=9)
        
        # 隐藏y轴和边框 | Hide y-axis and spines
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_position(('data', 0))
        
        # 添加标题 | Add title
        ax.set_title('LD Connections (r² values)', fontsize=11, pad=10)
    
    def _draw_heatmap(self, ax, ld_matrix: np.ndarray, positions: np.ndarray):
        """绘制下半部分的热图 | Draw lower heatmap"""
        # 如果只显示三角形，屏蔽下三角 | Mask lower triangle if triangle only
        if self.config.triangle_only:
            mask = np.tril(np.ones_like(ld_matrix, dtype=bool), k=-1)
        else:
            mask = None
        
        # 创建热图 | Create heatmap
        im = ax.imshow(ld_matrix, cmap=self.config.colormap, 
                      aspect='equal', interpolation='nearest',
                      vmin=0, vmax=1)
        
        # 如果有mask，应用它 | Apply mask if exists
        if mask is not None:
            masked_data = np.ma.array(ld_matrix, mask=mask)
            im = ax.imshow(masked_data, cmap=self.config.colormap, 
                          aspect='equal', interpolation='nearest',
                          vmin=0, vmax=1)
        
        # 添加颜色条 | Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8, aspect=20)
        cbar.set_label('r²', rotation=0, labelpad=15)
        
        # 设置坐标轴标签 | Set axis labels
        n_variants = len(positions)
        
        # 选择显示的刻度 | Select ticks to display
        if n_variants <= 20:
            tick_positions = range(n_variants)
            tick_labels = [f"{pos:,}" for pos in positions]
        else:
            n_ticks = min(10, n_variants)
            tick_indices = np.linspace(0, n_variants-1, n_ticks, dtype=int)
            tick_positions = tick_indices
            tick_labels = [f"{positions[i]:,}" for i in tick_indices]
        
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(tick_labels, fontsize=8)
        
        # 设置标签 | Set labels
        ax.set_xlabel('Genomic Position', fontsize=10)
        ax.set_ylabel('Genomic Position', fontsize=10)
        ax.set_title('LD Heatmap (r² values)', fontsize=12, pad=10)
