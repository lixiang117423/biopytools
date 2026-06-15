"""
VCF2PCACluster后端分析器|VCF2PCACluster Backend Analyzer
"""

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional


class V2PAnalyzer:
    """VCF2PCACluster分析器|VCF2PCACluster Analyzer"""

    def __init__(self, config, logger):
        """
        初始化V2P分析器|Initialize V2P analyzer

        Args:
            config: VCF2PCAConfig配置对象
            logger: 日志器
        """
        self.config = config
        self.logger = logger

    def run_analysis(self) -> bool:
        """
        运行完整的VCF2PCACluster分析流程|Run complete VCF2PCACluster analysis

        Returns:
            bool: 分析是否成功
        """
        self.logger.info("=" * 80)
        self.logger.info("开始VCF2PCACluster分析|Starting VCF2PCACluster analysis")
        self.logger.info("=" * 80)

        try:
            # 构建命令|Build command
            cmd = self._build_command()

            # 记录命令|Log command
            self.logger.info(f"执行命令|Executing command: {' '.join(cmd)}")

            # 执行命令|Execute command
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=self.config.output_dir
            )

            # 记录输出|Log output
            if result.stdout:
                for line in result.stdout.splitlines():
                    self.logger.info(f"VCF2PCACluster| {line}")

            # 检查输出文件|Check output files
            self._check_outputs()

            # 处理聚类结果|Process clustering results
            if self.config.cluster:
                self._process_clustering_results()

            self.logger.info("=" * 80)
            self.logger.info("VCF2PCACluster分析完成|VCF2PCACluster analysis completed")
            self.logger.info("=" * 80)

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"VCF2PCACluster执行失败|VCF2PCACluster execution failed: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误输出|Error output: {e.stderr}")
            return False

        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {e}")
            return False

    def _build_command(self) -> list:
        """
        构建VCF2PCACluster命令|Build VCF2PCACluster command

        Returns:
            list: 命令列表
        """
        cmd = [
            self.config.vcf2pca_path,
            '-InVCF', self.config.vcf_file,
            '-OutPut', os.path.join(self.config.output_dir, 'vcf2pca'),
            '-Threads', str(self.config.threads),
            '-PCnum', str(self.config.components)
        ]

        # 聚类参数|Clustering parameters
        if self.config.cluster:
            cmd.extend(['-ClusterMethod', self._get_cluster_method()])

            if self.config.cluster_method == 'kmeans':
                if self.config.cluster_k:
                    cmd.extend(['-BestKManually', str(self.config.cluster_k)])

        return cmd

    def _get_cluster_method(self) -> str:
        """
        获取聚类方法名称|Get clustering method name

        Returns:
            str: VCF2PCACluster格式的聚类方法
        """
        method_map = {
            'kmeans': 'Kmean',
            'dbscan': 'DBSCAN',
            'em': 'EM'
        }
        return method_map.get(self.config.cluster_method, 'Kmean')

    def _check_outputs(self):
        """检查并记录输出文件|Check and log output files"""
        output_prefix = os.path.join(self.config.output_dir, 'vcf2pca')

        # PCA结果文件（VCF2PCACluster实际输出的文件名）|PCA result files (actual VCF2PCACluster output names)
        pca_eig = f"{output_prefix}.eigenval"
        pca_eigvec = f"{output_prefix}.eigenvec"
        pca_eigvec_renamed = os.path.join(self.config.output_dir, 'pca_sample_point.txt')  # 重命名后的特征向量文件|Renamed eigenvector file
        pca_evaluation = f"{output_prefix}.evaluation"  # 聚类评估文件|Clustering evaluation file

        # Kinship文件|Kinship file
        kinship_file = f"{output_prefix}.Normalized_IBS.Kinship"

        # 检查文件存在性|Check file existence
        if os.path.exists(pca_eig):
            self._analyze_pca_eigenvectors(pca_eig)
        else:
            self.logger.warning(f"PCA特征值文件未找到|PCA eigenvalue file not found: {pca_eig}")

        if os.path.exists(pca_eigvec):
            self.logger.info(f"PCA特征向量文件|PCA eigenvector file: {pca_eigvec}")
            # 重命名特征向量文件|Rename eigenvector file
            try:
                os.rename(pca_eigvec, pca_eigvec_renamed)
                self.logger.info(f"特征向量文件已重命名|Eigenvector file renamed: {pca_eigvec} -> {pca_eigvec_renamed}")
            except Exception as e:
                self.logger.error(f"重命名特征向量文件失败|Failed to rename eigenvector file: {e}")
        else:
            self.logger.warning(f"PCA特征向量文件未找到|PCA eigenvector file not found: {pca_eigvec}")

        if os.path.exists(pca_evaluation):
            self.logger.info(f"聚类评估文件|Clustering evaluation file: {pca_evaluation}")

        if os.path.exists(kinship_file):
            self.logger.info(f"Kinship矩阵文件|Kinship matrix file: {kinship_file}")
        else:
            self.logger.warning(f"Kinship矩阵文件未找到|Kinship matrix file not found: {kinship_file}")

    def _analyze_pca_eigenvectors(self, eig_file: str):
        """
        分析PCA特征值|Analyze PCA eigenvalues

        Args:
            eig_file: 特征值文件路径
        """
        try:
            # 读取特征值|Read eigenvalues
            eigenvalues = []
            with open(eig_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # 跳过空行和注释行|Skip empty lines and comment lines
                    if not line or line.startswith('#'):
                        continue
                    # 分割并取第一列（特征值）|Split and take first column (eigenvalue)
                    parts = line.split()
                    if parts:
                        try:
                            eigenvalues.append(float(parts[0]))
                        except ValueError:
                            self.logger.warning(f"无法解析特征值行|Failed to parse eigenvalue line: {line}")
                            continue

            if not eigenvalues:
                self.logger.warning("未能读取特征值|Failed to read eigenvalues")
                return

            # 计算统计信息|Calculate statistics
            total_variance = sum(eigenvalues)
            explained_variance_ratio = [ev / total_variance for ev in eigenvalues]
            cumulative_variance = np.cumsum(explained_variance_ratio)

            # 记录结果|Log results
            self.logger.info("PCA分析结果|PCA Analysis Results:")
            self.logger.info(f"总方差|Total variance: {total_variance:.6f}")
            self.logger.info("主成分解释方差|Principal component explained variance:")

            for i, (ev, var_ratio, cum_var) in enumerate(zip(
                eigenvalues[:self.config.components],
                explained_variance_ratio[:self.config.components],
                cumulative_variance[:self.config.components]
            )):
                self.logger.info(
                    f"  PC{i+1}: 特征值|Eigenvalue = {ev:.6f}, "
                    f"解释方差比|Explained variance ratio = {var_ratio:.6f} ({var_ratio*100:.2f}%), "
                    f"累积|Cumulative = {cum_var:.6f} ({cum_var*100:.2f}%)"
                )

            # 保存详细结果|Save detailed results
            self._save_variance_table(explained_variance_ratio, cumulative_variance)

        except Exception as e:
            self.logger.error(f"分析PCA特征值失败|Failed to analyze PCA eigenvalues: {e}")

    def _save_variance_table(self, explained_variance, cumulative_variance):
        """
        保存解释方差表|Save explained variance table

        Args:
            explained_variance: 解释方差比列表
            cumulative_variance: 累积方差比列表
        """
        try:
            df = pd.DataFrame({
                'PC': [f'PC{i+1}' for i in range(len(explained_variance))],
                'Explained_Variance_Ratio': explained_variance,
                'Cumulative_Variance_Ratio': cumulative_variance,
                'Explained_Variance_Percent': [v * 100 for v in explained_variance],
                'Cumulative_Variance_Percent': [v * 100 for v in cumulative_variance]
            })

            output_file = os.path.join(self.config.output_dir, 'pca_explained_variance.txt')
            df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
            self.logger.info(f"解释方差表已保存|Explained variance table saved: {output_file}")

        except Exception as e:
            self.logger.error(f"保存解释方差表失败|Failed to save explained variance table: {e}")

    def _process_clustering_results(self):
        """处理聚类结果|Process clustering results"""
        cluster_file = os.path.join(self.config.output_dir, 'vcf2pca.cluster')

        if os.path.exists(cluster_file):
            self.logger.info(f"聚类结果文件|Clustering result file: {cluster_file}")

            try:
                # 读取聚类结果|Read clustering results
                df = pd.read_csv(cluster_file, sep=r'\s+', header=0)
                self.logger.info(f"聚类结果包含|Clustering result contains: {len(df)} 样本|samples")

                # 统计每个聚类的样本数|Count samples per cluster
                if 'Cluster' in df.columns:
                    cluster_counts = df['Cluster'].value_counts().sort_index()
                    self.logger.info("各聚类样本数|Samples per cluster:")
                    for cluster, count in cluster_counts.items():
                        self.logger.info(f"  聚类|Cluster {cluster}: {count} 样本|samples")

            except Exception as e:
                self.logger.error(f"读取聚类结果失败|Failed to read clustering results: {e}")
        else:
            self.logger.warning(f"聚类结果文件未找到|Clustering result file not found: {cluster_file}")

    def get_pca_results(self) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
        """
        获取PCA结果|Get PCA results

        Returns:
            Tuple[DataFrame, DataFrame]: (特征向量, 解释方差) 或 None
        """
        try:
            output_prefix = os.path.join(self.config.output_dir, 'vcf2pca')

            # 读取特征向量（优先读取重命名后的文件）|Read eigenvectors (prefer renamed file)
            eigvec_file = os.path.join(self.config.output_dir, 'pca_sample_point.txt')
            if not os.path.exists(eigvec_file):
                # 如果重命名文件不存在，尝试原始文件名|If renamed file doesn't exist, try original name
                eigvec_file = f"{output_prefix}.eigenvec"
            if not os.path.exists(eigvec_file):
                return None

            eigvec_df = pd.read_csv(eigvec_file, sep=r'\s+', header=None)
            eigvec_df.columns = ['Sample'] + [f'PC{i+1}' for i in range(eigvec_df.shape[1] - 1)]

            # 读取特征值|Read eigenvalues
            eig_file = f"{output_prefix}.eigenval"
            if not os.path.exists(eig_file):
                return None

            eigenvalues = []
            with open(eig_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # 跳过空行和注释行|Skip empty lines and comment lines
                    if not line or line.startswith('#'):
                        continue
                    # 分割并取第一列（特征值）|Split and take first column (eigenvalue)
                    parts = line.split()
                    if parts:
                        try:
                            eigenvalues.append(float(parts[0]))
                        except ValueError:
                            continue

            total_variance = sum(eigenvalues)
            explained_variance = pd.DataFrame({
                'PC': [f'PC{i+1}' for i in range(len(eigenvalues))],
                'Eigenvalue': eigenvalues,
                'Explained_Variance_Ratio': [ev / total_variance for ev in eigenvalues]
            })

            return eigvec_df, explained_variance

        except Exception as e:
            self.logger.error(f"获取PCA结果失败|Failed to get PCA results: {e}")
            return None
