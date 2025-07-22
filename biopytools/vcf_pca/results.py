"""
PCA结果处理模块 | PCA Results Processing Module
"""

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        report_file = self.config.output_path / "pca_analysis_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("VCF PCA分析总结报告 | VCF PCA Analysis Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            # 输入文件信息 | Input file information
            f.write("输入文件 | Input Files:\n")
            f.write(f"  - VCF文件 | VCF file: {self.config.vcf_file}\n")
            if self.config.sample_info_file:
                f.write(f"  - 样本信息文件 | Sample info file: {self.config.sample_info_file}\n")
            f.write("\n")
            
            # 分析参数 | Analysis parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write(f"  - 主成分数量 | Number of components: {self.config.components}\n")
            f.write(f"  - 跳过质控 | Skip QC: {'是 | Yes' if self.config.skip_qc else '否 | No'}\n")
            
            if not self.config.skip_qc:
                f.write(f"  - MAF阈值 | MAF threshold: {self.config.maf}\n")
                f.write(f"  - 缺失率阈值 | Missing rate threshold: {self.config.missing_rate}\n")
                f.write(f"  - HWE p值阈值 | HWE p-value threshold: {self.config.hwe_pvalue}\n")
            
            f.write(f"  - 生成可视化 | Generate visualization: {'是 | Yes' if self.config.plot else '否 | No'}\n")
            if self.config.group_column:
                f.write(f"  - 分组列 | Group column: {self.config.group_column}\n")
            f.write("\n")
            
            # 输出文件 | Output files
            f.write("输出文件 | Output Files:\n")
            f.write(f"  - pca_results.eigenval: 特征值文件 | Eigenvalues file\n")
            f.write(f"  - pca_results.eigenvec: 特征向量文件 | Eigenvectors file\n")
            f.write(f"  - pca_explained_variance.txt: 解释方差表 | Explained variance table\n")
            f.write(f"  - pca_eigenvectors_formatted.txt: 格式化的特征向量 | Formatted eigenvectors\n")
            
            if self.config.sample_info_file:
                f.write(f"  - pca_with_sample_info.txt: PCA结果与样本信息合并 | PCA results merged with sample info\n")
            
            if self.config.plot:
                f.write(f"  - pca_scree_plot.png: 碎石图 | Scree plot\n")
                f.write(f"  - pca_scatter_plots.png: PCA散点图 | PCA scatter plots\n")
                f.write(f"  - pca_pairs_plot.png: PCA成对图 | PCA pairs plot\n")
            
            f.write(f"\n输出目录 | Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
