"""
系统发育树样品选择计算模块|Phylogenetic Tree Sample Selector Calculator Module

基于均匀间隔和PCA去重的智能选择算法|Intelligent selection algorithm based on uniform interval and PCA deduplication
"""

import numpy as np
from typing import List, Dict, Optional
from scipy.spatial.distance import euclidean, pdist


class PhyloSelectorCalculator:
    """系统发育树样品选择计算器|Phylogenetic Tree Sample Selector Calculator"""

    def __init__(self, config, logger):
        """
        初始化计算器|Initialize calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

    def hierarchy_based_selection(self,
                                  all_samples: List[Dict],
                                  hierarchy_df,
                                  sample_to_pca: Dict[str, List[float]],
                                  sample_to_group: Optional[Dict[str, str]] = None) -> List[Dict]:
        """基于均匀间隔的智能选择|Uniform interval-based intelligent selection

        算法流程|Algorithm Flow:
        1. 按固定间隔均匀选择样品
        2. 检查顺序相邻性，避免选择树结构上挨着的样品
        3. PCA去重检查并替换

        Args:
            all_samples: 所有候选样品|All candidate samples
            hierarchy_df: 层级关系DataFrame|Hierarchy DataFrame (未使用|unused)
            sample_to_pca: PCA坐标映射|PCA coordinates mapping
            sample_to_group: 分组映射（可选，暂未使用）|Group mapping (optional, unused for now)

        Returns:
            List[Dict]: 选择的样品|Selected samples
        """
        target_count = min(self.config.n_samples, len(all_samples))

        self.logger.info("="*70)
        self.logger.info("基于均匀间隔的智能选择算法|Uniform Interval-Based Selection")
        self.logger.info("="*70)
        self.logger.info(f"总样品数|Total samples: {len(all_samples)}")
        self.logger.info(f"目标数量|Target count: {target_count}")
        self.logger.info(f"选择比例|Selection ratio: {target_count/len(all_samples)*100:.2f}%")

        # ========== 步骤1: 均匀间隔选择|Step 1: Uniform interval selection ==========
        self.logger.info(f"\n[步骤1|Step 1] 均匀间隔选择|Uniform interval selection...")

        # 计算间隔
        interval = len(all_samples) / target_count
        self.logger.info(f"  间隔|Interval: 每 {interval:.1f} 个样品选1个|Select 1 per {interval:.1f} samples")

        selected = []
        selected_indices = []

        for i in range(target_count):
            base_idx = int(i * interval)
            # 寻找合适的索引，避免与已选样品在顺序上相邻
            idx = self._find_non_adjacent_index(base_idx, selected_indices, len(all_samples))
            selected_indices.append(idx)
            selected.append(all_samples[idx])

        self.logger.info(f"  已选|Selected: {len(selected)} 个样品|samples")

        # ========== 步骤2: PCA去重检查|Step 2: PCA deduplication ==========
        self.logger.info(f"\n[步骤2|Step 2] PCA去重检查|PCA deduplication check...")
        selected = self._pca_dedup_with_replacement(
            selected,
            all_samples,
            sample_to_pca,
            target_count
        )

        # ========== 步骤3: 最终统计|Step 3: Final statistics ==========
        self.logger.info(f"\n[步骤3|Step 3] 最终统计|Final statistics...")

        # 验证最终结果中无相邻样品
        selected = self._ensure_no_adjacent_samples(selected, all_samples)

        # 计算PCA距离
        if sample_to_pca:
            final_pca = [sample_to_pca.get(s['name']) for s in selected]
            valid_pca = [p for p in final_pca if p is not None]
            if len(valid_pca) > 1:
                pca_dists = pdist(valid_pca)
                self.logger.info(f"  最小PCA距离|Min PCA distance: {np.min(pca_dists):.6f}")
                self.logger.info(f"  平均PCA距离|Avg PCA distance: {np.mean(pca_dists):.6f}")

                # 检查距离为0的对
                zero_dist_count = np.sum(np.array(pca_dists) < self.config.pca_dedup_threshold)
                if zero_dist_count > 0:
                    self.logger.warning(f"  仍有{zero_dist_count}对样品PCA距离<{self.config.pca_dedup_threshold}|Warning: {zero_dist_count} pairs with PCA distance < {self.config.pca_dedup_threshold}")
                else:
                    self.logger.info(f"  无重复样品（PCA距离都>{self.config.pca_dedup_threshold}）|No duplicate samples (all PCA distances > {self.config.pca_dedup_threshold})")

        self.logger.info(f"  最终样品数|Final sample count: {len(selected)}")

        return selected[:target_count]

    def _ensure_no_adjacent_samples(self, selected_samples: List[Dict],
                                    all_samples: List[Dict],
                                    min_gap: int = 2) -> List[Dict]:
        """确保最终结果中没有相邻的样品|Ensure no adjacent samples in final result

        Args:
            selected_samples: 已选样品|Selected samples
            all_samples: 所有候选样品|All candidate samples
            min_gap: 最小间隔|Minimum gap

        Returns:
            List[Dict]: 调整后的样品列表|Adjusted sample list
        """
        # 构建样品名到索引的映射
        sample_to_idx = {s['name']: i for i, s in enumerate(all_samples)}

        # 找到选中样品的索引
        selected_indices = [sample_to_idx.get(s['name'], -1) for s in selected_samples]
        selected_indices = [idx for idx in selected_indices if idx >= 0]

        # 按索引排序
        sorted_pairs = sorted(zip(selected_samples, selected_indices), key=lambda x: x[1])

        # 检查相邻性并调整
        adjusted_samples = []
        removed_indices = set()

        for i, (sample, idx) in enumerate(sorted_pairs):
            # 检查是否与已保留的样品相邻
            is_adjacent = False
            for kept_idx in removed_indices:
                if abs(idx - kept_idx) < min_gap:
                    is_adjacent = True
                    self.logger.info(f"  跳过相邻样品|Skip adjacent sample: {sample['name']} (索引|Index {idx}) 与 索引|Index {kept_idx} 相邻")
                    break

            if not is_adjacent:
                adjusted_samples.append(sample)
                removed_indices.add(idx)

        if len(adjusted_samples) < len(selected_samples):
            self.logger.info(f"  调整|Adjusted: {len(selected_samples)} -> {len(adjusted_samples)} (移除|removed {len(selected_samples) - len(adjusted_samples)} 个相邻样品|adjacent samples)")

        return adjusted_samples

    def _find_non_adjacent_index(self, base_idx: int, selected_indices: List[int],
                                 total_samples: int, min_gap: int = 2) -> int:
        """寻找不与已选样品相邻的索引|Find index not adjacent to selected samples

        Args:
            base_idx: 基础索引|Base index
            selected_indices: 已选样品的索引列表|List of selected sample indices
            total_samples: 总样品数|Total number of samples
            min_gap: 最小间隔|Minimum gap (default: 2)

        Returns:
            int: 合适的索引|Suitable index
        """
        # 检查base_idx是否与任何已选索引相邻
        for sel_idx in selected_indices:
            if abs(base_idx - sel_idx) < min_gap:
                # 相邻了，尝试在附近找替代索引
                # 向后搜索
                for offset in range(1, 10):
                    new_idx = min(base_idx + offset, total_samples - 1)
                    # 检查新索引是否与所有已选索引都相隔足够远
                    if all(abs(new_idx - s) >= min_gap for s in selected_indices):
                        return new_idx
                    # 也向前搜索
                    new_idx = max(base_idx - offset, 0)
                    if all(abs(new_idx - s) >= min_gap for s in selected_indices):
                        return new_idx
                # 如果找不到合适的，返回base_idx
                return base_idx
        return base_idx

    def _pca_dedup_with_replacement(self,
                                    selected_samples: List[Dict],
                                    all_samples: List[Dict],
                                    sample_to_pca: Dict[str, List[float]],
                                    target_count: int,
                                    max_iterations: int = 10) -> List[Dict]:
        """PCA去重并替换|PCA deduplication with replacement

        检查已选样品的PCA距离，如果太近则替换

        Args:
            selected_samples: 已选样品|Selected samples
            all_samples: 所有候选样品（用于替换）|All candidate samples (for replacement)
            sample_to_pca: PCA坐标映射|PCA coordinates mapping
            target_count: 目标数量|Target count
            max_iterations: 最大迭代次数|Max iterations

        Returns:
            List[Dict]: 去重后的样品|Deduplicated samples
        """
        iteration = 0
        selected = list(selected_samples)
        selected_names = set(s['name'] for s in selected)
        available_for_replacement = [s for s in all_samples if s['name'] not in selected_names]

        while iteration < max_iterations:
            iteration += 1

            # 计算已选样品的PCA距离矩阵
            close_pairs = []

            for i in range(len(selected)):
                for j in range(i+1, len(selected)):
                    name1 = selected[i]['name']
                    name2 = selected[j]['name']

                    if name1 in sample_to_pca and name2 in sample_to_pca:
                        dist = euclidean(sample_to_pca[name1], sample_to_pca[name2])

                        if dist < self.config.pca_dedup_threshold:
                            close_pairs.append((i, j, dist, name1, name2))

            if not close_pairs:
                self.logger.info(f"  迭代|Iteration {iteration}: 无PCA距离<{self.config.pca_dedup_threshold}的样品对|No pairs with PCA distance < {self.config.pca_dedup_threshold}")
                break

            self.logger.info(f"  迭代|Iteration {iteration}: 发现|Found {len(close_pairs)} 对PCA距离<{self.config.pca_dedup_threshold}的样品|pairs with PCA distance < {self.config.pca_dedup_threshold}")

            # 替换距离太近的样品中的一个
            to_remove = set()
            replacements = []

            for i, j, dist, name1, name2 in close_pairs:
                if name1 in to_remove or name2 in to_remove:
                    continue

                # 找到对应的样品对象
                sample1 = next(s for s in selected if s['name'] == name1)
                sample2 = next(s for s in selected if s['name'] == name2)

                # 保留索引较小的（即在层级文件中更靠前的）
                if i < j:
                    to_remove.add(name2)
                    removed_sample = sample2
                else:
                    to_remove.add(name1)
                    removed_sample = sample1

                self.logger.info(f"    移除|Remove {removed_sample['name']} (与|with {selected[i if name1 == removed_sample['name'] else j]['name']} 距离|distance={dist:.6f})")

                # 从可用池中找一个替代品
                replacement = self._find_replacement_sample(
                    removed_sample,
                    selected,
                    available_for_replacement,
                    sample_to_pca,
                    all_samples
                )

                if replacement:
                    replacements.append((removed_sample, replacement))

            # 执行替换
            for removed, repl in replacements:
                selected = [s for s in selected if s['name'] != removed['name']]
                selected.append(repl)
                selected_names.add(repl['name'])
                available_for_replacement = [s for s in available_for_replacement if s['name'] != repl['name']]
                self.logger.info(f"    替换为|Replace with {repl['name']}")

        # 最终检查
        if sample_to_pca:
            final_pca = [sample_to_pca.get(s['name']) for s in selected]
            valid_pca = [p for p in final_pca if p is not None]
            if len(valid_pca) > 1:
                pca_dists = pdist(valid_pca)
                zero_count = np.sum(np.array(pca_dists) < self.config.pca_dedup_threshold)
                if zero_count > 0:
                    self.logger.warning(f"  最终仍有{zero_count}对样品PCA距离<{self.config.pca_dedup_threshold}|Warning: {zero_count} pairs still with PCA distance < {self.config.pca_dedup_threshold}")
                else:
                    self.logger.info(f"  PCA去重完成|PCA dedup completed")

        return selected

    def _find_replacement_sample(self,
                                 removed_sample: Dict,
                                 current_selected: List[Dict],
                                 available_pool: List[Dict],
                                 sample_to_pca: Dict[str, List[float]],
                                 all_samples: List[Dict] = None,
                                 min_gap: int = 2) -> Dict:
        """寻找替换样品|Find replacement sample

        从可用池中找一个与所有已选样品PCA距离都>=阈值且不与已选样品相邻的样品

        Args:
            removed_sample: 被移除的样品|Removed sample
            current_selected: 当前已选样品|Current selected samples
            available_pool: 可用池|Available pool
            sample_to_pca: PCA坐标映射|PCA coordinates mapping
            all_samples: 所有样品列表（用于检查相邻性）|All samples list (for adjacency check)
            min_gap: 最小索引间隔|Minimum index gap

        Returns:
            Dict: 替换样品|Replacement sample (or None)
        """
        # 构建样品名到索引的映射（如果需要检查相邻性）
        sample_to_idx = {}
        if all_samples:
            sample_to_idx = {s['name']: i for i, s in enumerate(all_samples)}

        # 获取已选样品的索引
        selected_indices = set()
        if all_samples:
            for sel in current_selected:
                if sel['name'] in sample_to_idx:
                    selected_indices.add(sample_to_idx[sel['name']])

        # 尝试找一个合适的替换品
        for candidate in available_pool:
            # 检查相邻性
            if all_samples and candidate['name'] in sample_to_idx:
                candidate_idx = sample_to_idx[candidate['name']]
                is_adjacent = False
                for sel_idx in selected_indices:
                    if abs(candidate_idx - sel_idx) < min_gap:
                        is_adjacent = True
                        break
                if is_adjacent:
                    continue  # 跳过相邻的样品

            # 检查PCA距离
            if candidate['name'] not in sample_to_pca:
                return candidate  # 没有PCA数据，直接使用

            # 检查与所有已选样品的PCA距离
            min_dist = float('inf')
            for sel in current_selected:
                if sel['name'] in sample_to_pca:
                    dist = euclidean(sample_to_pca[candidate['name']], sample_to_pca[sel['name']])
                    if dist < min_dist:
                        min_dist = dist

            if min_dist >= self.config.pca_dedup_threshold:
                return candidate

        # 如果找不到完全满足的，返回第一个不违反相邻性的
        if available_pool and all_samples:
            for candidate in available_pool:
                if candidate['name'] in sample_to_idx:
                    candidate_idx = sample_to_idx[candidate['name']]
                    is_adjacent = any(abs(candidate_idx - sel_idx) < min_gap for sel_idx in selected_indices)
                    if not is_adjacent:
                        return candidate

        # 最后的备选：返回第一个
        if available_pool:
            return available_pool[0]

        return None
