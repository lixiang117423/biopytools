"""
最大化系统发育多样性算法|Maximize Phylogenetic Diversity Algorithm
基于AI建议实现：使用Cophenetic Distance和贪婪算法
Implement based on AI suggestion: Use Cophenetic Distance and greedy algorithm
"""

import numpy as np
from typing import Dict, List, Tuple, Set
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


class PhylogeneticDiversitySelector:
    """系统发育多样性选择器|Phylogenetic Diversity Selector"""
    
    def __init__(self, tree_analyzer, logger=None):
        """初始化|Initialize
        
        Args:
            tree_analyzer: 树结构分析器|Tree structure analyzer
            logger: 日志记录器|Logger
        """
        self.tree_analyzer = tree_analyzer
        self.logger = logger
        
    def calculate_cophenetic_distance_matrix(self, samples: List[str]) -> np.ndarray:
        """计算Cophenetic距离矩阵（树上距离矩阵）|Calculate Cophenetic distance matrix
        
        Args:
            samples: 样品列表|Sample list
            
        Returns:
            距离矩阵|Distance matrix
        """
        n = len(samples)
        dist_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = self.tree_analyzer.get_tree_distance(samples[i], samples[j])
                if dist is not None:
                    dist_matrix[i][j] = dist
                    dist_matrix[j][i] = dist
        
        return dist_matrix
    
    def deduplicate_by_pca_2d(self, samples: List[Dict],
                              sample_to_pca: Dict[str, List[float]],
                              threshold: float = 0.01) -> List[Dict]:
        """在2D PCA空间（PC1+PC2）中去重|Deduplicate in 2D PCA space (PC1+PC2)

        如果两个样品在PCA平面上的距离小于阈值，保留枝长较长的|If two samples are too close in PCA space, keep longer BL

        Args:
            samples: 样品列表|Sample list
            sample_to_pca: PCA坐标|PCA coordinates
            threshold: PCA距离阈值|PCA distance threshold

        Returns:
            去重后的样品|Deduplicated samples
        """
        # 提取PC1和PC2|Extract PC1 and PC2
        sample_pc = {}
        for s in samples:
            if s['name'] in sample_to_pca:
                pca = sample_to_pca[s['name']]
                if pca and len(pca) >= 2:
                    sample_pc[s['name']] = (pca[0], pca[1])

        if len(sample_pc) < 2:
            return samples

        # 按枝长排序（优先保留长枝长）|Sort by branch length (keep longer BL)
        sorted_samples = sorted(samples, key=lambda x: x['branch_length'], reverse=True)

        # 贪心去重|Greedy deduplication
        kept = []
        removed = []

        from scipy.spatial.distance import euclidean
        for candidate in sorted_samples:
            if candidate['name'] not in sample_pc:
                kept.append(candidate)
                continue

            # 检查与已保留样品的距离|Check distance to kept samples
            is_too_close = False
            for kept_sample in kept:
                if kept_sample['name'] in sample_pc:
                    dist = euclidean(sample_pc[candidate['name']], sample_pc[kept_sample['name']])
                    if dist < threshold:
                        is_too_close = True
                        removed.append(candidate['name'])
                        break

            if not is_too_close:
                kept.append(candidate)

        if self.logger:
            self.logger.info(f"  2D PCA去重（阈值|threshold={threshold}）:")
            self.logger.info(f"    移除|Removed: {len(removed)} 个样品|samples")
            self.logger.info(f"    保留|Kept: {len(kept)} 个样品|samples")

        return kept

    def adaptive_deduplicate_by_pca_2d(self, samples: List[Dict],
                                       sample_to_pca: Dict[str, List[float]],
                                       max_close_pairs_ratio: float = 0.1) -> List[Dict]:
        """自适应2D PCA去重|Adaptive 2D PCA deduplication

        自动分析PCA距离分布，根据分位数动态调整去重阈值
        Automatically analyze PCA distance distribution and adjust threshold based on quantiles

        Args:
            samples: 样品列表|Sample list
            sample_to_pca: PCA坐标|PCA coordinates
            max_close_pairs_ratio: 允许的最大close pairs比例|Max allowed ratio of close pairs

        Returns:
            去重后的样品|Deduplicated samples
        """
        # 提取PC1和PC2|Extract PC1 and PC2
        sample_pc = {}
        for s in samples:
            if s['name'] in sample_to_pca:
                pca = sample_to_pca[s['name']]
                if pca and len(pca) >= 2:
                    sample_pc[s['name']] = (pca[0], pca[1])

        if len(sample_pc) < 2:
            return samples

        # 计算所有PCA距离|Calculate all PCA distances
        from scipy.spatial.distance import euclidean
        sample_names = list(sample_pc.keys())
        all_distances = []

        for i in range(len(sample_names)):
            for j in range(i+1, len(sample_names)):
                dist = euclidean(sample_pc[sample_names[i]], sample_pc[sample_names[j]])
                all_distances.append(dist)

        all_distances = np.array(all_distances)

        # 分析距离分布|Analyze distance distribution
        min_dist = all_distances.min()
        p5_dist = np.percentile(all_distances, 5)
        p10_dist = np.percentile(all_distances, 10)
        median_dist = np.percentile(all_distances, 50)

        if self.logger:
            self.logger.info(f"  PCA距离分布分析|PCA distance distribution:")
            self.logger.info(f"    最小值|Min: {min_dist:.4f}")
            self.logger.info(f"    5%分位数|5th percentile: {p5_dist:.4f}")
            self.logger.info(f"    10%分位数|10th percentile: {p10_dist:.4f}")
            self.logger.info(f"    中位数|Median: {median_dist:.4f}")

        # 统计close pairs（距离<0.01）|Count close pairs (dist < 0.01)
        n_total_pairs = len(all_distances)
        n_close_pairs = np.sum(all_distances < 0.01)
        close_ratio = n_close_pairs / n_total_pairs

        if self.logger:
            self.logger.info(f"  Close pairs统计|Close pairs statistics:")
            self.logger.info(f"    距离<0.01的对数|Pairs with dist<0.01: {n_close_pairs}/{n_total_pairs} ({close_ratio*100:.1f}%)")

        # 如果close pairs比例过高，进行去重|If too many close pairs, deduplicate
        if close_ratio > max_close_pairs_ratio:
            # 使用5%分位数作为阈值（保留至少95%的样品对）|Use 5th percentile as threshold (keep 95% of pairs)
            threshold = max(p5_dist, 0.01)  # 最小0.01|Minimum 0.01

            if self.logger:
                self.logger.info(f"  Close pairs过多，进行自适应去重|Too many close pairs, performing adaptive dedup:")
                self.logger.info(f"    使用阈值|Using threshold: {threshold:.4f}")

            # 按枝长排序（优先保留长枝长）|Sort by branch length (keep longer BL)
            sorted_samples = sorted(samples, key=lambda x: x['branch_length'], reverse=True)

            # 贪心去重|Greedy deduplication
            kept = []
            removed = []

            for candidate in sorted_samples:
                if candidate['name'] not in sample_pc:
                    kept.append(candidate)
                    continue

                # 检查与已保留样品的距离|Check distance to kept samples
                is_too_close = False
                for kept_sample in kept:
                    if kept_sample['name'] in sample_pc:
                        dist = euclidean(sample_pc[candidate['name']], sample_pc[kept_sample['name']])
                        if dist < threshold:
                            is_too_close = True
                            removed.append(candidate['name'])
                            break

                if not is_too_close:
                    kept.append(candidate)

            if self.logger:
                self.logger.info(f"    移除|Removed: {len(removed)} 个样品|samples")
                self.logger.info(f"    保留|Kept: {len(kept)} 个样品|samples")

            return kept
        else:
            if self.logger:
                self.logger.info(f"  Close pairs比例可接受，无需去重|Close pairs ratio acceptable, no dedup needed")
            return samples

    def maximize_phylogenetic_diversity_greedy(self, samples: List[Dict],
                                              sample_to_pca: Dict[str, List[float]],
                                              target_count: int,
                                              sample_to_group: Dict[str, str] = None,
                                              group_samples: List[Dict] = None) -> List[Dict]:
        """贪婪算法最大化系统发育多样性（带迭代优化）|Greedy algorithm with iterative optimization

        算法流程|Algorithm flow:
        1. 如果有分组预选样品，先锁定它们|If group samples exist, lock them first
        2. PCA空间分层采样（每个象限选择2倍名额）|PCA stratified sampling (2x quota per quadrant)
        3. 迭代优化（类似梯度下降）|Iterative optimization (like gradient descent):
           - 评估当前选择的PCA距离分布|Evaluate PCA distance distribution
           - 如果close pairs过多，增加PCA权重重新选择|If too many close pairs, increase PCA weight and reselect
           - 重复直到满足目标或达到最大迭代次数|Repeat until target met or max iterations
        4. 计算Cophenetic距离矩阵|Calculate Cophenetic distance matrix
        5. 贪心迭代：每次选择对已有集合多样性贡献最大的样品|Greedy iteration: select sample with max contribution

        Args:
            samples: 候选样品列表|Candidate sample list
            sample_to_pca: PCA坐标|PCA coordinates
            target_count: 目标数量|Target count
            sample_to_group: 分组映射|Group mapping (optional)
            group_samples: 预选的分组样品|Pre-selected group samples (optional)

        Returns:
            选择的样品|Selected samples
        """
        if self.logger:
            self.logger.info(f"最大化系统发育多样性（迭代优化）|Maximizing PD with Iterative Optimization")
            self.logger.info(f"  候选样品|Candidates: {len(samples)}")
            self.logger.info(f"  目标数量|Target: {target_count}")
            self.logger.info(f"  使用PCA空间分层采样|Using PCA-stratified sampling for balance")

        # 按PCA空间分层，确保均衡分布|Stratify by PCA space to ensure balanced distribution
        if sample_to_pca:
            samples = self._stratify_by_pca_space(samples, sample_to_pca, target_count)
            if self.logger:
                self.logger.info(f"  PCA分层后样品数|After PCA stratification: {len(samples)}")
        
        sample_names = [s['name'] for s in samples]
        n_samples = len(sample_names)
        
        # 计算Cophenetic距离矩阵|Calculate Cophenetic distance matrix
        if self.logger:
            self.logger.info(f"  计算Cophenetic距离矩阵|Calculating Cophenetic distance matrix...")
        
        tree_dist_matrix = self.calculate_cophenetic_distance_matrix(sample_names)
        
        # 计算PCA距离矩阵（用于辅助筛选）|Calculate PCA distance matrix (for auxiliary filtering)
        pca_dist_matrix = None
        if sample_to_pca:
            pca_dist_matrix = np.zeros((n_samples, n_samples))
            pca_values = [sample_to_pca.get(name) for name in sample_names]
            valid_pca = [i for i, v in enumerate(pca_values) if v is not None]
            
            for i in valid_pca:
                for j in valid_pca:
                    if i < j:
                        from scipy.spatial.distance import euclidean
                        dist = euclidean(pca_values[i], pca_values[j])
                        pca_dist_matrix[i][j] = dist
                        pca_dist_matrix[j][i] = dist
        
        # 初始化已选集合|Initialize selected set
        selected = set()
        selected_list = []
        
        # 如果有分组预选样品，先加入|If group samples exist, add them first
        if group_samples:
            for gs in group_samples:
                if gs['name'] in sample_names:
                    selected.add(gs['name'])
                    selected_list.append(gs)
            if self.logger:
                self.logger.info(f"  锁定分组样品|Locked group samples: {len(selected)}")
        
        # 贪心迭代|Greedy iteration
        if self.logger:
            self.logger.info(f"  开始贪心选择|Starting greedy selection...")
        
        iteration = 0
        while len(selected) < target_count and len(selected) < n_samples:
            iteration += 1
            
            # 找到贡献最大的样品|Find sample with maximum contribution
            best_sample = None
            best_contribution = -float('inf')
            
            for i, name in enumerate(sample_names):
                if name in selected:
                    continue
                
                # 计算这个样品的贡献：到所有已选样品的最小树距离
                # Contribution: min tree distance to all selected samples
                if len(selected) == 0:
                    # 第一个样品：选择枝长最大的|First sample: select max branch length
                    contribution = samples[i]['branch_length']
                else:
                    # 计算到所有已选样品的最小树距离|Min tree distance to all selected
                    min_tree_dist = float('inf')
                    min_pca_dist = float('inf')

                    for sel_name in selected:
                        if sel_name in sample_names:
                            j = sample_names.index(sel_name)
                            tree_dist = tree_dist_matrix[i][j]
                            if tree_dist < min_tree_dist:
                                min_tree_dist = tree_dist

                            if pca_dist_matrix is not None:
                                pca_dist = pca_dist_matrix[i][j]
                                if pca_dist < min_pca_dist:
                                    min_pca_dist = pca_dist

                    # 综合贡献：树距离和PCA距离同等重要|Tree distance and PCA distance equally important
                    contribution = min_tree_dist
                    if pca_dist_matrix is not None:
                        # PCA距离权重加倍|PCA distance with 2x weight
                        contribution += min_pca_dist * 2.0  # 200%的权重给PCA多样性|200% weight for PCA diversity
                
                if contribution > best_contribution:
                    best_contribution = contribution
                    best_sample = name
            
            if best_sample:
                selected.add(best_sample)
                idx = sample_names.index(best_sample)
                selected_list.append(samples[idx])
                
                if iteration % 10 == 0 or len(selected) >= target_count:
                    if self.logger:
                        self.logger.info(f"    迭代|Iteration {iteration}: 已选|Selected {len(selected)}/{target_count}, 最新添加|Latest: {best_sample}")
            else:
                break
        
        if self.logger:
            self.logger.info(f"  贪心选择完成|Greedy selection completed: {len(selected_list)} 个样品|samples")
        
        return selected_list
    
    def tree_based_clustering_selection(self, samples: List[Dict],
                                       sample_to_pca: Dict[str, List[float]],
                                       target_count: int,
                                       n_clusters: int = None,
                                       sample_to_group: Dict[str, str] = None,
                                       group_samples: List[Dict] = None) -> List[Dict]:
        """基于系统发育树聚类的选择|Tree-based clustering selection
        
        算法流程|Algorithm flow:
        1. 使用Cophenetic距离进行层次聚类|Hierarchical clustering using Cophenetic distance
        2. 将树切割成n_clusters个簇|Cut tree into n_clusters clusters
        3. 每个簇内选择代表性样品|Select representative from each cluster
        4. 结合PCA避免选择离群点|Combine PCA to avoid outliers
        
        Args:
            samples: 候选样品|Candidate samples
            sample_to_pca: PCA坐标|PCA coordinates
            target_count: 目标数量|Target count
            n_clusters: 聚类数量|Number of clusters (auto-calculate if None)
            sample_to_group: 分组映射|Group mapping
            group_samples: 预选分组样品|Pre-selected group samples
            
        Returns:
            选择的样品|Selected samples
        """
        if self.logger:
            self.logger.info(f"基于树的聚类选择|Tree-based Clustering Selection")
            self.logger.info(f"  候选样品|Candidates: {len(samples)}")
            self.logger.info(f"  目标数量|Target: {target_count}")
        
        sample_names = [s['name'] for s in samples]
        
        # 计算Cophenetic距离矩阵|Calculate Cophenetic distance matrix
        if self.logger:
            self.logger.info(f"  计算Cophenetic距离矩阵|Calculating Cophenetic distance matrix...")
        
        tree_dist_matrix = self.calculate_cophenetic_distance_matrix(sample_names)
        
        # 层次聚类|Hierarchical clustering
        if n_clusters is None:
            # 自动计算聚类数：目标数量的1/2到2/3之间|Auto-calculate: 1/2 to 2/3 of target
            n_clusters = max(10, target_count // 2)
        
        if self.logger:
            self.logger.info(f"  层次聚类|Hierarchical clustering into {n_clusters} clusters...")
        
        # 将距离矩阵转换为condensed形式|Convert to condensed form
        condensed_dist = squareform(tree_dist_matrix)
        
        # 使用平均链接聚类|Use average linkage clustering
        Z = linkage(condensed_dist, method='average')
        
        # 切割树得到聚类|Cut tree to get clusters
        cluster_labels = fcluster(Z, t=n_clusters, criterion='maxclust')
        
        # 按聚类组织样品|Organize samples by cluster
        clusters = defaultdict(list)
        for idx, label in enumerate(cluster_labels):
            clusters[label].append(samples[idx])
        
        if self.logger:
            self.logger.info(f"  聚类结果|Clustering result:")
            for label, cluster_samples in sorted(clusters.items()):
                self.logger.info(f"    聚类|Cluster {label}: {len(cluster_samples)} 个样品|samples")
        
        # 为每个聚类分配名额|Allocate quota for each cluster
        total_available = len(samples)
        locked_count = len(group_samples) if group_samples else 0
        remaining_slots = target_count - locked_count
        
        cluster_quotas = {}
        for label, cluster_samples in clusters.items():
            quota = max(1, int(len(cluster_samples) / total_available * remaining_slots))
            cluster_quotas[label] = quota
        
        # 调整名额总和|Adjust total quota
        total_quota = sum(cluster_quotas.values())
        if total_quota < remaining_slots:
            # 增加大聚类的名额|Increase quota for large clusters
            while total_quota < remaining_slots:
                max_cluster = max(clusters.items(), key=lambda x: len(x[1]))[0]
                cluster_quotas[max_cluster] += 1
                total_quota += 1
        elif total_quota > remaining_slots:
            # 减少小聚类的名额|Decrease quota for small clusters
            while total_quota > remaining_slots:
                min_cluster = min(cluster_quotas.items(), key=lambda x: clusters[x[0]][0])[0]
                if cluster_quotas[min_cluster] > 1:
                    cluster_quotas[min_cluster] -= 1
                    total_quota -= 1
                else:
                    break
        
        # 每个聚类内选择代表|Select representative from each cluster
        selected = []
        if group_samples:
            selected.extend(group_samples)
        
        if self.logger:
            self.logger.info(f"  每个聚类内选择代表|Selecting representatives from each cluster...")
        
        for label, cluster_samples in sorted(clusters.items()):
            quota = cluster_quotas[label]
            if len(cluster_samples) == 0:
                continue
            
            if self.logger:
                self.logger.info(f"    聚类|Cluster {label}: 选择|selecting {quota}/{len(cluster_samples)}")
            
            # 在聚类内选择：优先选择靠近PCA中心的样品|Select within cluster: prioritize PCA center
            cluster_selected = self._select_from_cluster(
                cluster_samples, sample_to_pca, quota
            )
            selected.extend(cluster_selected)
        
        if self.logger:
            self.logger.info(f"  聚类选择完成|Cluster selection completed: {len(selected)} 个样品|samples")
        
        return selected[:target_count]  # 确保不超过目标|Ensure not exceeding target
    
    def _stratify_by_pca_space(self, samples: List[Dict],
                               sample_to_pca: Dict[str, List[float]],
                               target_count: int) -> List[Dict]:
        """按PCA空间分层采样（5x5网格：箱线图风格）|Stratified sampling by PCA space (5x5 grid: boxplot style)

        在PC1和PC2的2D空间中按箱线图风格划分5x5网格，更好地处理极端值
        Divide PC1+PC2 2D space into 5x5 grid using boxplot-style bins to better handle extreme values

        5个区间|5 bins:
        1. 下离群区间|Lower outliers (< Q1 - 1.5*IQR)
        2. Q1边界到Q1|Q1 boundary to Q1
        3. Q1到Q3（中间50%，包含中位数）|Q1 to Q3 (central 50%, contains median)
        4. Q3到Q3边界|Q3 to Q3 boundary
        5. 上离群区间|Upper outliers (> Q3 + 1.5*IQR)

        Args:
            samples: 候选样品|Candidate samples
            sample_to_pca: PCA坐标|PCA coordinates
            target_count: 目标数量|Target count

        Returns:
            分层后的样品列表|Stratified sample list
        """
        # 提取PC1和PC2值|Extract PC1 and PC2 values
        sample_pc = {}
        for s in samples:
            if s['name'] in sample_to_pca:
                pca = sample_to_pca[s['name']]
                if pca and len(pca) >= 2:
                    sample_pc[s['name']] = (pca[0], pca[1])  # (PC1, PC2)

        if len(sample_pc) == 0:
            return samples

        # 提取PC1和PC2的值|Extract PC1 and PC2 values
        pc1_values = [v[0] for v in sample_pc.values()]
        pc2_values = [v[1] for v in sample_pc.values()]

        # 计算箱线图参数|Calculate boxplot parameters
        pc1_q25 = np.percentile(pc1_values, 25)  # Q1
        pc1_q75 = np.percentile(pc1_values, 75)  # Q3
        pc1_iqr = pc1_q75 - pc1_q25  # IQR
        pc1_lower_outlier = pc1_q25 - 1.5 * pc1_iqr  # 下离群阈值|Lower outlier threshold
        pc1_upper_outlier = pc1_q75 + 1.5 * pc1_iqr  # 上离群阈值|Upper outlier threshold

        pc2_q25 = np.percentile(pc2_values, 25)  # Q1
        pc2_q75 = np.percentile(pc2_values, 75)  # Q3
        pc2_iqr = pc2_q75 - pc2_q25  # IQR
        pc2_lower_outlier = pc2_q25 - 1.5 * pc2_iqr  # 下离群阈值|Lower outlier threshold
        pc2_upper_outlier = pc2_q75 + 1.5 * pc2_iqr  # 上离群阈值|Upper outlier threshold

        if self.logger:
            self.logger.info(f"  PC1箱线图参数|PC1 boxplot: Q25={pc1_q25:.4f}, Q75={pc1_q75:.4f}, "
                           f"IQR={pc1_iqr:.4f}, 下离群|low_out={pc1_lower_outlier:.4f}, "
                           f"上离群|high_out={pc1_upper_outlier:.4f}")
            self.logger.info(f"  PC2箱线图参数|PC2 boxplot: Q25={pc2_q25:.4f}, Q75={pc2_q75:.4f}, "
                           f"IQR={pc2_iqr:.4f}, 下离群|low_out={pc2_lower_outlier:.4f}, "
                           f"上离群|high_out={pc2_upper_outlier:.4f}")

        # 创建5x5网格|Create 5x5 grid
        # 区间1: <下离群阈值, 区间2: 下离群阈值~Q1, 区间3: Q1~Q3, 区间4: Q3~上离群阈值, 区间5: >上离群阈值
        # Bin 1: <lower_outlier, Bin 2: lower_outlier~Q1, Bin 3: Q1~Q3, Bin 4: Q3~upper_outlier, Bin 5: >upper_outlier
        strata = {}
        for i in range(5):
            for j in range(5):
                strata[f'B{i+1}B{j+1}'] = []

        for s in samples:
            if s['name'] in sample_pc:
                pc1, pc2 = sample_pc[s['name']]

                # 确定PC1所在的行|Determine PC1 row
                if pc1 < pc1_lower_outlier:
                    row = 0  # 下离群区间|Lower outliers
                elif pc1 < pc1_q25:
                    row = 1  # Q1边界到Q1|Q1 boundary to Q1
                elif pc1 < pc1_q75:
                    row = 2  # Q1到Q3（中间50%）|Q1 to Q3 (central 50%)
                elif pc1 < pc1_upper_outlier:
                    row = 3  # Q3到Q3边界|Q3 to Q3 boundary
                else:
                    row = 4  # 上离群区间|Upper outliers

                # 确定PC2所在的列|Determine PC2 column
                if pc2 < pc2_lower_outlier:
                    col = 0  # 下离群区间|Lower outliers
                elif pc2 < pc2_q25:
                    col = 1  # Q1边界到Q1|Q1 boundary to Q1
                elif pc2 < pc2_q75:
                    col = 2  # Q1到Q3（中间50%）|Q1 to Q3 (central 50%)
                elif pc2 < pc2_upper_outlier:
                    col = 3  # Q3到Q3边界|Q3 to Q3 boundary
                else:
                    col = 4  # 上离群区间|Upper outliers

                strata[f'B{row+1}B{col+1}'].append(s)
            else:
                # 没有PCA数据的放入中间|No PCA data -> center
                strata['B3B3'].append(s)

        if self.logger:
            self.logger.info(f"  PCA空间分层结果（5x5网格：箱线图风格）|PCA space strata (5x5 grid: boxplot style):")
            for stratum_name in [f'B{i}B{j}' for i in range(1, 6) for j in range(1, 6)]:
                stratum_samples = strata.get(stratum_name, [])
                if stratum_samples:
                    pc1s = [sample_pc.get(s['name'], (0, 0))[0] for s in stratum_samples]
                    pc2s = [sample_pc.get(s['name'], (0, 0))[1] for s in stratum_samples]
                    self.logger.info(f"    {stratum_name}: {len(stratum_samples):3d} 个|samples "
                                   f"(PC1: {np.min(pc1s):7.4f}~{np.max(pc1s):7.4f}, "
                                   f"PC2: {np.min(pc2s):7.4f}~{np.max(pc2s):7.4f})")

        # 计算每个分层应该选择的数量|Calculate quota for each stratum
        total_with_pca = sum(len(s) for s in strata.values())
        strata_quotas = {}
        for stratum_name, stratum_samples in strata.items():
            if total_with_pca > 0 and len(stratum_samples) > 0:
                quota = max(1, int(len(stratum_samples) / total_with_pca * target_count))
                strata_quotas[stratum_name] = quota

        # 调整名额总和|Adjust total quota
        total_quota = sum(strata_quotas.values())
        if total_quota < target_count:
            while total_quota < target_count:
                max_stratum = max(strata.items(), key=lambda x: len(x[1]))[0]
                strata_quotas[max_stratum] += 1
                total_quota += 1
        elif total_quota > target_count:
            while total_quota > target_count:
                candidates = [(name, quota) for name, quota in strata_quotas.items() if quota > 1]
                if not candidates:
                    break
                min_stratum = min(candidates, key=lambda x: x[1])[0]
                strata_quotas[min_stratum] -= 1
                total_quota -= 1

        if self.logger:
            self.logger.info(f"  网格名额|Grid quotas: {strata_quotas}")

        # 从每个分层中采样|Sample from each stratum
        # 每个网格选择2倍名额|Select 2x quota from each grid cell
        stratified_samples = []
        for stratum_name in [f'B{i}B{j}' for i in range(1, 6) for j in range(1, 6)]:
            stratum_samples = strata.get(stratum_name, [])
            quota = strata_quotas.get(stratum_name, 0) * 2  # 2倍名额|2x quota

            if len(stratum_samples) == 0 or quota == 0:
                continue

            # 按枝长排序，选择前quota个|Sort by BL, select top quota
            sorted_samples = sorted(stratum_samples, key=lambda x: x['branch_length'], reverse=True)
            if len(sorted_samples) > quota:
                sorted_samples = sorted_samples[:quota]

            stratified_samples.extend(sorted_samples)

        if self.logger:
            self.logger.info(f"  分层采样完成|Stratified sampling completed: {len(stratified_samples)} 个样品|samples")

        return stratified_samples

    def _deduplicate_quadrant_pca(self, quadrant_samples: List[Dict],
                                  sample_pc: Dict[str, Tuple[float, float]],
                                  threshold: float = 0.01) -> List[Dict]:
        """在单个象限内进行PCA去重|PCA deduplication within a single quadrant

        Args:
            quadrant_samples: 象限内的样品|Samples in quadrant
            sample_pc: PCA坐标|PCA coordinates
            threshold: 距离阈值|Distance threshold

        Returns:
            去重后的样品|Deduplicated samples
        """
        # 按枝长排序（优先保留长枝长）|Sort by branch length (keep longer BL)
        sorted_samples = sorted(quadrant_samples, key=lambda x: x['branch_length'], reverse=True)

        # 贪心去重|Greedy deduplication
        kept = []

        from scipy.spatial.distance import euclidean
        for candidate in sorted_samples:
            if candidate['name'] not in sample_pc:
                kept.append(candidate)
                continue

            # 检查与已保留样品的距离|Check distance to kept samples
            is_too_close = False
            for kept_sample in kept:
                if kept_sample['name'] in sample_pc:
                    dist = euclidean(sample_pc[candidate['name']], sample_pc[kept_sample['name']])
                    if dist < threshold:
                        is_too_close = True
                        break

            if not is_too_close:
                kept.append(candidate)

        return kept

    def _greedy_select_max_pca_distance(self, stratum_samples: List[Dict],
                                       sample_pc: Dict[str, Tuple[float, float]],
                                       n_select: int) -> List[Dict]:
        """在分层内使用贪心算法最大化PCA距离|Greedy selection to maximize PCA distance within stratum

        Args:
            stratum_samples: 分层内的样品|Samples in stratum
            sample_pc: 样品PCA坐标|Sample PCA coordinates
            n_select: 选择数量|Number to select

        Returns:
            选择的样品|Selected samples
        """
        if len(stratum_samples) <= n_select:
            return stratum_samples

        # 贪心选择：每次选择离已选样品最远的|Greedy: select farthest from selected
        selected = []
        remaining = list(stratum_samples)

        # 第一个：选择枝长最大的|First: select max branch length
        remaining.sort(key=lambda x: x['branch_length'], reverse=True)
        selected.append(remaining.pop(0))

        # 迭代选择|Iterative selection
        from scipy.spatial.distance import euclidean
        while len(selected) < n_select and remaining:
            best_sample = None
            max_min_dist = -1

            for candidate in remaining:
                if candidate['name'] not in sample_pc:
                    continue

                # 计算到所有已选样品的最小距离|Min distance to all selected
                min_dist = float('inf')
                for sel in selected:
                    if sel['name'] in sample_pc:
                        dist = euclidean(sample_pc[candidate['name']], sample_pc[sel['name']])
                        if dist < min_dist:
                            min_dist = dist

                if min_dist > max_min_dist:
                    max_min_dist = min_dist
                    best_sample = candidate

            if best_sample:
                selected.append(best_sample)
                remaining.remove(best_sample)
            else:
                # 如果找不到合适的，按枝长选择|If not found, select by BL
                selected.append(remaining.pop(0))

        return selected

    def _select_from_cluster(self, cluster_samples: List[Dict],
                            sample_to_pca: Dict[str, List[float]],
                            n_select: int) -> List[Dict]:
        """在聚类内选择样品|Select samples within a cluster

        策略|Strategy:
        1. 如果有PCA数据，优先选择靠近PCA中心的样品|If PCA available, prioritize samples near PCA center
        2. 其次选择枝长较大的样品|Then select samples with larger branch length

        Args:
            cluster_samples: 聚类内样品|Samples in cluster
            sample_to_pca: PCA坐标|PCA coordinates
            n_select: 选择数量|Number to select

        Returns:
            选择的样品|Selected samples
        """
        if len(cluster_samples) <= n_select:
            return cluster_samples

        # 计算每个样品到PCA中心的距离|Calculate distance to PCA center
        if sample_to_pca:
            pca_values = []
            for s in cluster_samples:
                if s['name'] in sample_to_pca:
                    pca_values.append(sample_to_pca[s['name']])

            if len(pca_values) > 0:
                # 计算PCA中心|Calculate PCA center
                pca_center = np.mean(pca_values, axis=0)

                # 计算每个样品到中心的距离|Calculate distance to center for each sample
                scored_samples = []
                for s in cluster_samples:
                    if s['name'] in sample_to_pca:
                        from scipy.spatial.distance import euclidean
                        dist_to_center = euclidean(sample_to_pca[s['name']], pca_center)
                        # 综合得分：靠近中心 + 枝长大|Combined score: near center + large BL
                        score = -dist_to_center + s['branch_length'] * 10  # 枝长权重|BL weight
                        scored_samples.append((s, score))
                    else:
                        scored_samples.append((s, s['branch_length'] * 10))

                # 按得分排序，选择最高的|Sort by score, select highest
                scored_samples.sort(key=lambda x: x[1], reverse=True)
                return [s for s, score in scored_samples[:n_select]]

        # 如果没有PCA数据，按枝长选择|If no PCA, select by branch length
        sorted_samples = sorted(cluster_samples, key=lambda x: x['branch_length'], reverse=True)
        return sorted_samples[:n_select]

    # ========== 混合距离算法（拓扑距离 + Cophenetic距离）|Hybrid Distance Algorithm ==========

    def maximize_pd_with_topology_filter(self, samples: List[Dict],
                                         sample_to_pca: Dict[str, List[float]],
                                         target_count: int,
                                         sample_to_group: Dict[str, str] = None,
                                         group_samples: List[Dict] = None,
                                         topology_multiplier: float = 2.0) -> List[Dict]:
        """使用拓扑距离预筛选 + Cophenetic距离精确选择|Use topology distance for pre-filtering + Cophenetic for precise selection

        方案B实现|Scheme B Implementation:
        1. PCA空间分层（5x5箱线图网格）|PCA space stratification (5x5 boxplot grid)
        2. 每个网格内：
           - 用拓扑距离快速预筛选（2倍名额）|Use topology distance for fast pre-filtering (2x quota)
           - 用Cophenetic距离精确贪心选择|Use Cophenetic distance for precise greedy selection
        3. 自适应PCA去重|Adaptive PCA deduplication

        Args:
            samples: 候选样品列表|Candidate sample list
            sample_to_pca: PCA坐标|PCA coordinates
            target_count: 目标数量|Target count
            sample_to_group: 分组映射|Group mapping (optional)
            group_samples: 预选的分组样品|Pre-selected group samples (optional)
            topology_multiplier: 拓扑距离预筛选倍数|Topology distance pre-filter multiplier

        Returns:
            选择的样品|Selected samples
        """
        if self.logger:
            self.logger.info(f"混合距离算法|Hybrid Distance Algorithm (Topology + Cophenetic)")
            self.logger.info(f"  候选样品|Candidates: {len(samples)}")
            self.logger.info(f"  目标数量|Target: {target_count}")
            self.logger.info(f"  拓扑距离预筛选倍数|Topology pre-filter multiplier: {topology_multiplier}x")

        # ========== 第1步：PCA空间分层|Step 1: PCA Space Stratification ==========
        if sample_to_pca:
            stratified_samples = self._stratify_by_pca_space(samples, sample_to_pca, target_count)
            if self.logger:
                self.logger.info(f"  PCA分层后样品数|After PCA stratification: {len(stratified_samples)}")
        else:
            stratified_samples = samples

        # ========== 第2步：每个PCA网格内使用拓扑距离预筛选|Step 2: Topology Pre-filter within each PCA grid ==========
        # 获取层级关系字典|Get hierarchy dictionary
        hierarchy = self.tree_analyzer.extract_hierarchy_dict()

        if not hierarchy:
            if self.logger:
                self.logger.warning(f"  无法提取层级关系，回退到纯Cophenetic距离算法|Cannot extract hierarchy, fallback to pure Cophenetic algorithm")
            return self.maximize_phylogenetic_diversity_greedy(
                samples, sample_to_pca, target_count, sample_to_group, group_samples
            )

        # 按PCA网格分组|Group by PCA grid
        sample_pc = {}
        for s in stratified_samples:
            if s['name'] in sample_to_pca:
                pca = sample_to_pca[s['name']]
                if pca and len(pca) >= 2:
                    sample_pc[s['name']] = (pca[0], pca[1])

        # 计算箱线图参数|Calculate boxplot parameters
        pc1_values = [v[0] for v in sample_pc.values()]
        pc2_values = [v[1] for v in sample_pc.values()]
        pc1_q25, pc1_q75 = np.percentile(pc1_values, [25, 75])
        pc1_iqr = pc1_q75 - pc1_q25
        pc1_lower_outlier = pc1_q25 - 1.5 * pc1_iqr
        pc1_upper_outlier = pc1_q75 + 1.5 * pc1_iqr
        pc2_q25, pc2_q75 = np.percentile(pc2_values, [25, 75])
        pc2_iqr = pc2_q75 - pc2_q25
        pc2_lower_outlier = pc2_q25 - 1.5 * pc2_iqr
        pc2_upper_outlier = pc2_q75 + 1.5 * pc2_iqr

        # 分配样品到5x5网格|Assign samples to 5x5 grid
        grid_samples = {}
        for s in stratified_samples:
            if s['name'] not in sample_pc:
                continue
            pc1, pc2 = sample_pc[s['name']]

            # 确定行|Determine row
            if pc1 < pc1_lower_outlier:
                row = 0
            elif pc1 < pc1_q25:
                row = 1
            elif pc1 < pc1_q75:
                row = 2
            elif pc1 < pc1_upper_outlier:
                row = 3
            else:
                row = 4

            # 确定列|Determine column
            if pc2 < pc2_lower_outlier:
                col = 0
            elif pc2 < pc2_q25:
                col = 1
            elif pc2 < pc2_q75:
                col = 2
            elif pc2 < pc2_upper_outlier:
                col = 3
            else:
                col = 4

            grid_key = f'B{row+1}B{col+1}'
            if grid_key not in grid_samples:
                grid_samples[grid_key] = []
            grid_samples[grid_key].append(s)

        if self.logger:
            self.logger.info(f"  PCA网格分布|PCA grid distribution: {len(grid_samples)} 个网格|grids")

        # 计算每个网格的目标名额|Calculate quota for each grid
        total_stratified = len(stratified_samples)
        grid_quotas = {}
        for grid_key, grid_samps in grid_samples.items():
            quota = max(1, int(len(grid_samps) / total_stratified * target_count))
            grid_quotas[grid_key] = quota

        # 调整名额总和|Adjust total quota
        total_quota = sum(grid_quotas.values())
        if total_quota < target_count:
            while total_quota < target_count:
                max_grid = max(grid_samples.items(), key=lambda x: len(x[1]))[0]
                grid_quotas[max_grid] += 1
                total_quota += 1

        # ========== 第3步：每个网格内使用拓扑距离预筛选，再用Cophenetic距离精确选择|Step 3: Topology pre-filter + Cophenetic precise selection ==========
        final_selected = []
        if group_samples:
            final_selected.extend(group_samples)

        for grid_key in sorted(grid_samples.keys()):
            grid_samps = grid_samples[grid_key]
            quota = grid_quotas.get(grid_key, 0)

            if len(grid_samps) == 0 or quota == 0:
                continue

            if self.logger:
                self.logger.info(f"    网格|Grid {grid_key}: {len(grid_samps)} 个候选|candidates, 目标|target {quota}")

            # 使用拓扑距离预筛选|Use topology distance for pre-filtering
            pre_filter_count = min(int(quota * topology_multiplier), len(grid_samps))
            pre_filtered = self._topology_based_selection(
                grid_samps, hierarchy, pre_filter_count
            )

            if self.logger:
                self.logger.info(f"      拓扑预筛选|Topology pre-filter: {len(pre_filtered)} 个|samples")

            # 使用Cophenetic距离精确选择|Use Cophenetic distance for precise selection
            if len(pre_filtered) > quota:
                pre_filtered_names = [s['name'] for s in pre_filtered]

                # 计算Cophenetic距离矩阵|Calculate Cophenetic distance matrix
                tree_dist_matrix = self.calculate_cophenetic_distance_matrix(pre_filtered_names)

                # 贪心选择|Greedy selection
                selected = set()
                selected_list = []

                for _ in range(quota):
                    best_sample = None
                    best_contribution = -float('inf')

                    for i, s in enumerate(pre_filtered):
                        if s['name'] in selected:
                            continue

                        if len(selected) == 0:
                            contribution = s['branch_length']
                        else:
                            # 计算到已选样品的最小Cophenetic距离|Min Cophenetic distance to selected
                            min_dist = float('inf')
                            for selected_s in selected_list:
                                j = pre_filtered_names.index(selected_s['name'])
                                dist = tree_dist_matrix[i][j]
                                if dist < min_dist:
                                    min_dist = dist
                            contribution = min_dist

                        if contribution > best_contribution:
                            best_contribution = contribution
                            best_sample = s

                    if best_sample:
                        selected.add(best_sample['name'])
                        selected_list.append(best_sample)

                final_selected.extend(selected_list)
            else:
                final_selected.extend(pre_filtered)

        if self.logger:
            self.logger.info(f"  混合距离选择完成|Hybrid distance selection completed: {len(final_selected)} 个样品|samples")

        # ========== 第4步：自适应PCA去重|Step 4: Adaptive PCA Deduplication ==========
        if sample_to_pca:
            final_selected = self.adaptive_deduplicate_by_pca_2d(final_selected, sample_to_pca)

        return final_selected[:target_count]

    def _topology_based_selection(self, samples: List[Dict],
                                   hierarchy: Dict[str, List[int]],
                                   n_select: int) -> List[Dict]:
        """基于拓扑距离的选择|Topology distance based selection

        使用贪心算法最大化拓扑距离（LCA距离）|Use greedy algorithm to maximize topology distance (LCA distance)

        Args:
            samples: 候选样品|Candidate samples
            hierarchy: 层级关系字典|Hierarchy dictionary
            n_select: 选择数量|Number to select

        Returns:
            选择的样品|Selected samples
        """
        if len(samples) <= n_select:
            return samples

        # 辅助函数：计算拓扑距离|Helper: calculate topology distance
        def calc_topo_dist(name1, name2):
            if name1 not in hierarchy or name2 not in hierarchy:
                return float('inf')
            path1 = hierarchy[name1]
            path2 = hierarchy[name2]

            # 找LCA|Find LCA
            lca_depth = 0
            for n1, n2 in zip(path1, path2):
                if n1 == n2:
                    lca_depth += 1
                else:
                    break

            return len(path1) + len(path2) - 2 * lca_depth

        # 贪心选择|Greedy selection
        selected = []
        remaining = list(samples)

        # 第一个：选择枝长最大的|First: select max branch length
        remaining.sort(key=lambda x: x['branch_length'], reverse=True)
        selected.append(remaining.pop(0))

        # 迭代选择|Iterative selection
        while len(selected) < n_select and remaining:
            best_sample = None
            best_min_dist = -1

            for candidate in remaining:
                # 计算到已选样品的最小拓扑距离|Min topology distance to selected
                min_dist = min(calc_topo_dist(candidate['name'], s['name']) for s in selected)

                if min_dist > best_min_dist:
                    best_min_dist = min_dist
                    best_sample = candidate

            if best_sample:
                selected.append(best_sample)
                remaining.remove(best_sample)
            else:
                break

        return selected
