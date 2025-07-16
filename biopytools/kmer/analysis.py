"""
K-mer分析模块 | K-mer Analysis Module
"""

import pandas as pd
from typing import Dict, List, Tuple

try:
    from sklearn.mixture import BayesianGaussianMixture
    import matplotlib.pyplot as plt
    import seaborn as sns
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

class KmerMatrixProcessor:
    """K-mer矩阵处理器 | K-mer Matrix Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_query_results(self, gene_kmers: Dict[str, List[Tuple[str, int]]]) -> pd.DataFrame:
        """处理查询结果生成最终矩阵 | Process query results to generate final matrix"""
        self.logger.info("处理查询结果生成k-mer矩阵 | Processing query results to generate k-mer matrix...")
        
        # 读取查询结果
        query_results = {}
        with open(self.config.query_result_file, 'r') as f:
            header = f.readline().strip().split()
            sample_names = header[1:]  # 跳过第一列(k-mer)
            
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    kmer = parts[0]
                    values = [int(x) for x in parts[1:]]
                    query_results[kmer] = values
        
        # 生成最终矩阵
        matrix_data = []
        kmer_id = 1
        
        for gene_name, kmer_list in gene_kmers.items():
            for kmer, position in kmer_list:
                row = {
                    'gene_id': gene_name,
                    'kmer_id': kmer_id,
                    'kmer_sequence': kmer
                }
                
                # 添加样本数据
                if kmer in query_results:
                    for i, sample_name in enumerate(sample_names):
                        if i < len(query_results[kmer]):
                            row[sample_name] = query_results[kmer][i]
                        else:
                            row[sample_name] = 0
                else:
                    # 如果k-mer未找到，全部设为0
                    for sample_name in sample_names:
                        row[sample_name] = 0
                
                matrix_data.append(row)
                kmer_id += 1
        
        # 创建DataFrame
        df = pd.DataFrame(matrix_data)
        
        # 重新排列列的顺序
        sample_columns = [col for col in df.columns if col not in ['gene_id', 'kmer_id', 'kmer_sequence']]
        final_columns = ['gene_id', 'kmer_id', 'kmer_sequence'] + sorted(sample_columns)
        df = df[final_columns]
        
        # 保存矩阵
        df.to_csv(self.config.kmer_matrix_final, sep='\t', index=False)
        self.logger.info(f"✓ 最终k-mer矩阵已保存 | Final k-mer matrix saved: {self.config.kmer_matrix_final}")
        self.logger.info(f"矩阵维度 | Matrix dimensions: {df.shape[0]} k-mer × {len(sample_columns)} 样本 | samples")
        
        return df

class HaplotypeAnalyzer:
    """单倍型分析器 | Haplotype Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        if not VISUALIZATION_AVAILABLE:
            self.logger.warning("可视化功能不可用，请安装: pip install scikit-learn matplotlib seaborn | Visualization unavailable, install: pip install scikit-learn matplotlib seaborn")
    
    def run_haplotype_analysis(self, kmer_matrix: pd.DataFrame) -> bool:
        """运行单倍型聚类分析 | Run haplotype clustering analysis"""
        if not self.config.run_haplotype:
            return True
        
        if not VISUALIZATION_AVAILABLE:
            self.logger.error("单倍型分析需要scikit-learn等包，请安装 | Haplotype analysis requires scikit-learn packages, please install")
            return False
            
        self.logger.info("运行单倍型聚类分析 | Running haplotype clustering analysis...")
        
        try:
            # 准备数据
            sample_columns = [col for col in kmer_matrix.columns if col not in ['gene_id', 'kmer_id', 'kmer_sequence']]
            genotype_matrix = kmer_matrix[sample_columns].T.values
            
            # BGMM聚类
            bgmm = BayesianGaussianMixture(
                n_components=self.config.n_components,
                random_state=42,
                weight_concentration_prior=1e-2
            )
            sample_groups = bgmm.fit_predict(genotype_matrix)
            
            # 保存结果
            results_df = pd.DataFrame({
                'Sample_ID': sample_columns,
                'Haplotype_Group': sample_groups
            })
            
            results_df.to_csv(self.config.haplotype_results, sep='\t', index=False)
            self.logger.info(f"✓ 单倍型分析结果已保存 | Haplotype analysis results saved: {self.config.haplotype_results}")
            
            # 生成热图 (可选)
            if len(sample_columns) <= 1000:  # 避免过大的热图
                self.generate_haplotype_heatmap(kmer_matrix, results_df)
            
            return True
            
        except Exception as e:
            self.logger.error(f"单倍型分析失败 | Haplotype analysis failed: {e}")
            return False
    
    def generate_haplotype_heatmap(self, kmer_matrix: pd.DataFrame, haplotype_results: pd.DataFrame):
        """生成单倍型热图 | Generate haplotype heatmap"""
        try:
            self.logger.info("生成单倍型热图 | Generating haplotype heatmap...")
            
            # 准备数据
            sample_columns = [col for col in kmer_matrix.columns if col not in ['gene_id', 'kmer_id', 'kmer_sequence']]
            plot_data = kmer_matrix[sample_columns].T
            
            # 按单倍型分组排序
            sample_order = haplotype_results.sort_values('Haplotype_Group')['Sample_ID'].tolist()
            plot_data = plot_data.loc[sample_order]
            
            # 生成热图
            plt.figure(figsize=(12, 8))
            sns.heatmap(plot_data.iloc[:100, :50], cmap='viridis', cbar=True)  # 只显示前100个样本和50个k-mer
            plt.title('K-mer Haplotype Heatmap')
            plt.xlabel('K-mer Position')
            plt.ylabel('Sample ID')
            plt.tight_layout()
            plt.savefig(self.config.heatmap_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"✓ 单倍型热图已保存 | Haplotype heatmap saved: {self.config.heatmap_file}")
            
        except Exception as e:
            self.logger.warning(f"热图生成失败 | Heatmap generation failed: {e}")
