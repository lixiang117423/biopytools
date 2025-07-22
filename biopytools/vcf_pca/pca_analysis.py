"""
PCA分析核心模块 | PCA Analysis Core Module
"""

import numpy as np
import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class PCAAnalyzer:
    """PCA分析器 | PCA Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_pca_analysis(self):
        """运行PCA分析 | Run PCA analysis"""
        input_prefix = self.config.plink_prefix_qc
        output_prefix = self.config.output_path / "pca_results"
        
        # 确保输出目录存在 | Ensure output directory exists
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix_abs = output_prefix.resolve()
        
        self.logger.info(f"开始PCA分析 | Starting PCA analysis")
        self.logger.info(f"主成分数量 | Number of components: {self.config.components}")
        
        # 运行PLINK PCA | Run PLINK PCA
        cmd = (
            f"{self.config.plink_path} --bfile {input_prefix} "
            f"--pca {self.config.components} --out {output_prefix_abs} --allow-extra-chr"
        )
        
        success = self.cmd_runner.run(cmd, f"PCA分析 | PCA analysis ({self.config.components} components)")
        
        if success:
            self.config.pca_prefix = str(output_prefix_abs)
            self.analyze_pca_results()
        
        return success
    
    def analyze_pca_results(self):
        """分析PCA结果 | Analyze PCA results"""
        self.logger.info("分析PCA结果 | Analyzing PCA results")
        
        # 读取特征值 | Read eigenvalues
        eigenval_file = Path(f"{self.config.pca_prefix}.eigenval")
        eigenvec_file = Path(f"{self.config.pca_prefix}.eigenvec")
        
        if eigenval_file.exists() and eigenvec_file.exists():
            # 计算解释方差比 | Calculate explained variance ratio
            eigenvalues = np.loadtxt(eigenval_file)
            total_variance = np.sum(eigenvalues)
            explained_variance_ratio = eigenvalues / total_variance
            cumulative_variance = np.cumsum(explained_variance_ratio)
            
            # 记录解释方差 | Log explained variance
            self.logger.info("主成分解释方差 | Principal component explained variance:")
            for i, (var_ratio, cum_var) in enumerate(zip(explained_variance_ratio, cumulative_variance)):
                self.logger.info(f"  PC{i+1}: {var_ratio:.4f} ({var_ratio*100:.2f}%), 累积 | cumulative: {cum_var:.4f} ({cum_var*100:.2f}%)")
            
            # 保存详细结果 | Save detailed results
            self.save_detailed_results(eigenvalues, explained_variance_ratio, cumulative_variance)
        
        else:
            self.logger.warning("PCA结果文件不存在 | PCA result files not found")
    
    def save_detailed_results(self, eigenvalues, explained_variance_ratio, cumulative_variance):
        """保存详细结果 | Save detailed results"""
        # 保存解释方差表 | Save explained variance table
        variance_df = pd.DataFrame({
            'PC': [f'PC{i+1}' for i in range(len(eigenvalues))],
            'Eigenvalue': eigenvalues,
            'Explained_Variance_Ratio': explained_variance_ratio,
            'Cumulative_Variance_Ratio': cumulative_variance,
            'Explained_Variance_Percent': explained_variance_ratio * 100,
            'Cumulative_Variance_Percent': cumulative_variance * 100
        })
        
        variance_file = self.config.output_path / "pca_explained_variance.txt"
        variance_df.to_csv(variance_file, sep='\t', index=False)
        self.logger.info(f"解释方差表已保存 | Explained variance table saved: {variance_file}")
        
        # 读取并处理特征向量 | Read and process eigenvectors
        eigenvec_df = pd.read_csv(f"{self.config.pca_prefix}.eigenvec", sep=r'\s+', header=None)
        
        # 重命名列 | Rename columns
        columns = ['FID', 'IID'] + [f'PC{i+1}' for i in range(len(eigenvec_df.columns) - 2)]
        eigenvec_df.columns = columns
        
        # 保存格式化的特征向量 | Save formatted eigenvectors
        eigenvec_formatted_file = self.config.output_path / "pca_eigenvectors_formatted.txt"
        eigenvec_df.to_csv(eigenvec_formatted_file, sep='\t', index=False)
        self.logger.info(f"格式化特征向量已保存 | Formatted eigenvectors saved: {eigenvec_formatted_file}")
        
        # 如果有样本信息，进行合并 | Merge with sample info if available
        if self.config.sample_info_file:
            self.merge_sample_info(eigenvec_df)
    
    def merge_sample_info(self, eigenvec_df):
        """合并样本信息 | Merge sample information"""
        try:
            self.logger.info("合并样本信息 | Merging sample information")
            
            # 读取样本信息 | Read sample information
            sample_info = pd.read_csv(self.config.sample_info_file, sep='\t')
            
            # 确定合并键 | Determine merge key
            if 'IID' in sample_info.columns:
                merge_key = 'IID'
            elif 'sample' in sample_info.columns:
                sample_info = sample_info.rename(columns={'sample': 'IID'})
                merge_key = 'IID'
            else:
                # 假设第一列是样本ID | Assume first column is sample ID
                first_col = sample_info.columns[0]
                sample_info = sample_info.rename(columns={first_col: 'IID'})
                merge_key = 'IID'
            
            # 合并数据 | Merge data
            merged_df = pd.merge(eigenvec_df, sample_info, on=merge_key, how='left')
            
            # 保存合并结果 | Save merged results
            merged_file = self.config.output_path / "pca_with_sample_info.txt"
            merged_df.to_csv(merged_file, sep='\t', index=False)
            self.logger.info(f"PCA结果与样本信息已合并 | PCA results merged with sample info: {merged_file}")
            
            # 保存到配置用于可视化 | Save to config for visualization
            self.config.merged_pca_data = merged_df
            
        except Exception as e:
            self.logger.error(f"合并样本信息失败 | Failed to merge sample information: {e}")
