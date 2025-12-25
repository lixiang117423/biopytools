"""
PCA可视化模块 | PCA Visualization Module
"""

import numpy as np
import pandas as pd
from pathlib import Path

class PCAVisualizer:
    """PCA可视化器 | PCA Visualizer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def create_plots(self):
        """创建PCA图表 | Create PCA plots"""
        if not self.config.plot:
            return True
        
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            self.logger.warning("matplotlib/seaborn未安装，跳过可视化 | matplotlib/seaborn not installed, skipping visualization")
            return True
        
        self.logger.info("创建PCA可视化图表 | Creating PCA visualization plots")
        
        # 读取PCA数据 | Read PCA data
        if hasattr(self.config, 'merged_pca_data'):
            pca_data = self.config.merged_pca_data
        else:
            eigenvec_file = self.config.output_path / "pca_eigenvectors_formatted.txt"
            if eigenvec_file.exists():
                pca_data = pd.read_csv(eigenvec_file, sep='\t')
            else:
                self.logger.warning("未找到PCA数据文件 | PCA data file not found")
                return False
        
        # 创建图表 | Create plots
        self.plot_scree_plot()
        self.plot_pca_scatter(pca_data)
        self.plot_pca_pairs(pca_data)
        
        return True
    
    def plot_scree_plot(self):
        """绘制碎石图 | Plot scree plot"""
        try:
            import matplotlib.pyplot as plt
            
            eigenval_file = Path(f"{self.config.pca_prefix}.eigenval")
            if not eigenval_file.exists():
                return
            
            eigenvalues = np.loadtxt(eigenval_file)
            explained_variance_ratio = eigenvalues / np.sum(eigenvalues)
            
            plt.figure(figsize=(12, 5))
            
            # 子图1: 解释方差比 | Subplot 1: Explained variance ratio
            plt.subplot(1, 2, 1)
            plt.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, 'bo-', linewidth=2, markersize=6)
            plt.xlabel('主成分 | Principal Component', fontsize=12)
            plt.ylabel('解释方差比 | Explained Variance Ratio', fontsize=12)
            plt.title('碎石图 | Scree Plot', fontsize=14, fontweight='bold')
            plt.grid(True, alpha=0.3)
            
            # 子图2: 累积解释方差 | Subplot 2: Cumulative explained variance
            plt.subplot(1, 2, 2)
            cumulative_variance = np.cumsum(explained_variance_ratio)
            plt.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, 'ro-', linewidth=2, markersize=6)
            plt.xlabel('主成分 | Principal Component', fontsize=12)
            plt.ylabel('累积解释方差比 | Cumulative Explained Variance Ratio', fontsize=12)
            plt.title('累积解释方差 | Cumulative Explained Variance', fontsize=14, fontweight='bold')
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            scree_file = self.config.output_path / "pca_scree_plot.png"
            plt.savefig(scree_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"碎石图已保存 | Scree plot saved: {scree_file}")
            
        except Exception as e:
            self.logger.error(f"绘制碎石图失败 | Failed to plot scree plot: {e}")
    
    def plot_pca_scatter(self, pca_data):
        """绘制PCA散点图 | Plot PCA scatter plot"""
        try:
            import matplotlib.pyplot as plt
            
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            
            # PC1 vs PC2
            self.scatter_subplot(axes[0, 0], pca_data, 'PC1', 'PC2')
            
            # PC1 vs PC3
            self.scatter_subplot(axes[0, 1], pca_data, 'PC1', 'PC3')
            
            # PC2 vs PC3
            self.scatter_subplot(axes[1, 0], pca_data, 'PC2', 'PC3')
            
            # PC1 vs PC4 (如果存在)
            if 'PC4' in pca_data.columns:
                self.scatter_subplot(axes[1, 1], pca_data, 'PC1', 'PC4')
            else:
                axes[1, 1].axis('off')
            
            plt.tight_layout()
            scatter_file = self.config.output_path / "pca_scatter_plots.png"
            plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"PCA散点图已保存 | PCA scatter plots saved: {scatter_file}")
            
        except Exception as e:
            self.logger.error(f"绘制PCA散点图失败 | Failed to plot PCA scatter: {e}")
    
    def scatter_subplot(self, ax, data, x_col, y_col):
        """绘制单个散点图子图 | Plot single scatter subplot"""
        if self.config.group_column and self.config.group_column in data.columns:
            # 按组分色 | Color by group
            groups = data[self.config.group_column].unique()
            colors = plt.cm.Set1(np.linspace(0, 1, len(groups)))
            
            for group, color in zip(groups, colors):
                group_data = data[data[self.config.group_column] == group]
                ax.scatter(group_data[x_col], group_data[y_col], 
                          c=[color], label=group, alpha=0.7, s=50)
            
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            # 单色散点图 | Single color scatter
            ax.scatter(data[x_col], data[y_col], alpha=0.7, s=50)
        
        ax.set_xlabel(x_col, fontsize=11)
        ax.set_ylabel(y_col, fontsize=11)
        ax.set_title(f'{x_col} vs {y_col}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    def plot_pca_pairs(self, pca_data):
        """绘制PCA成对图 | Plot PCA pair plots"""
        try:
            import seaborn as sns
            import matplotlib.pyplot as plt
            
            # 选择前4个主成分 | Select first 4 PCs
            pc_cols = [col for col in pca_data.columns if col.startswith('PC')][:4]
            
            if len(pc_cols) < 2:
                return
            
            # 创建成对图 | Create pair plot
            if self.config.group_column and self.config.group_column in pca_data.columns:
                g = sns.pairplot(pca_data[pc_cols + [self.config.group_column]], 
                               hue=self.config.group_column, 
                               plot_kws={'alpha': 0.7, 's': 40})
            else:
                g = sns.pairplot(pca_data[pc_cols], 
                               plot_kws={'alpha': 0.7, 's': 40})
            
            pairs_file = self.config.output_path / "pca_pairs_plot.png"
            g.savefig(pairs_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"PCA成对图已保存 | PCA pairs plot saved: {pairs_file}")
            
        except Exception as e:
            self.logger.error(f"绘制PCA成对图失败 | Failed to plot PCA pairs: {e}")
