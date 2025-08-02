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
        
        # 设置图形 | Setup figure
        plt.figure(figsize=self.config.figsize, dpi=self.config.dpi)
        
        # 如果只显示三角形，屏蔽下三角 | Mask lower triangle if triangle only
        if self.config.triangle_only:
            mask = np.tril(np.ones_like(ld_matrix, dtype=bool), k=-1)
        else:
            mask = None
        
        # 创建热图 | Create heatmap
        sns.heatmap(ld_matrix, 
                    mask=mask,
                    cmap=self.config.colormap,
                    square=True,
                    cbar_kws={"shrink": 0.8, "label": "r²"},
                    xticklabels=False,
                    yticklabels=False)
        
        # 设置标题 | Set title
        if self.config.title:
            plt.title(self.config.title, fontsize=14, pad=20)
        else:
            plt.title('Linkage Disequilibrium Heatmap', fontsize=14, pad=20)
        
        # 添加位置信息 | Add position information
        n_ticks = min(10, len(positions))
        tick_indices = np.linspace(0, len(positions)-1, n_ticks, dtype=int)
        tick_positions = positions[tick_indices]
        tick_labels = [f"{pos:,}" for pos in tick_positions]
        
        plt.xticks(tick_indices, tick_labels, rotation=45, ha='right')
        plt.yticks(tick_indices, tick_labels)
        plt.xlabel('Genomic Position')
        plt.ylabel('Genomic Position')
        
        # 调整布局 | Adjust layout
        plt.tight_layout()
        
        # 保存图片 | Save figure
        plt.savefig(self.config.output_file, dpi=self.config.dpi, bbox_inches='tight')
        if self.config.verbose:
            self.logger.info(f"热图已保存到 | Heatmap saved to: {self.config.output_file}")
        
        plt.close()
