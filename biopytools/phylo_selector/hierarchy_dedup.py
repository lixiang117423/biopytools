#!/usr/bin/env python3
"""
层级约束去重|Hierarchy-based Deduplication

确保在进化树的深层分支上不会选择过多样品
Ensure that not too many samples are selected from deep branches of the phylogenetic tree
"""

import pandas as pd
from typing import List, Dict, Set, Tuple


def hierarchy_based_dedup(samples: List[Dict],
                          hierarchy_file: str,
                          start_level: int = 25,
                          max_samples_per_node: int = 1,
                          logger=None) -> List[Dict]:
    """基于层级关系的去重|Deduplication based on hierarchy

    从指定层级开始，检查每个节点，如果超过max_samples_per_node个样品被选中，
    则只保留枝长最大的样品

    Args:
        samples: 候选样品列表（已按枝长降序排序）|Candidate samples (sorted by BL descending)
        hierarchy_file: 层级关系Excel文件|Hierarchy Excel file
        start_level: 开始检查的层级（从深层开始）|Level to start checking (from deep levels)
        max_samples_per_node: 每个节点最多保留的样品数|Max samples to keep per node
        logger: 日志对象|Logger object

    Returns:
        去重后的样品列表|Deduplicated sample list
    """
    if logger:
        logger.info(f"[层级约束去重|Hierarchy-based Deduplication]")
        logger.info(f"  起始层级|Start level: parent_{start_level}")
        logger.info(f"  每节点最大样品数|Max samples per node: {max_samples_per_node}")

    # 读取层级关系|Read hierarchy
    hierarchy = pd.read_excel(hierarchy_file)
    sample_names = {s['name'] for s in samples}

    # 过滤出在候选列表中的样品|Filter samples in candidate list
    hierarchy = hierarchy[hierarchy['label'].isin(sample_names)].copy()

    # 标记样品是否已保留|Mark if sample is kept
    sample_kept = {s['name']: False for s in samples}
    sample_to_bl = {s['name']: s['branch_length'] for s in samples}

    # 按层级从深到浅检查|Check from deep to shallow levels
    n_removed = 0
    n_kept = 0

    for level in range(start_level, 0, -1):
        col = f'parent_{level}'
        if col not in hierarchy.columns:
            continue

        # 按该层级的节点分组|Group by nodes at this level
        grouped = hierarchy[hierarchy[col].notna()].groupby(col)

        for node_id, group in grouped:
            # 该节点下候选样品中已经保留的|Already kept from this node
            kept_in_node = [s for s in group['label'] if sample_kept[s]]

            if len(kept_in_node) >= max_samples_per_node:
                # 已经达到上限，跳过|Already reached limit, skip
                continue

            # 该节点下所有候选样品（包括未保留的）|All candidates from this node
            candidates = list(group['label'])

            # 按枝长排序（从大到小）|Sort by branch length (descending)
            candidates.sort(key=lambda x: sample_to_bl[x], reverse=True)

            # 选择前max_samples_per_node个|Select top max_samples_per_node
            to_keep = candidates[:max_samples_per_node]

            # 标记保留|Mark as kept
            for sample in to_keep:
                if not sample_kept[sample]:
                    sample_kept[sample] = True
                    n_kept += 1

            # 移除其他的|Remove others
            to_remove = [s for s in candidates if s not in to_keep and not sample_kept[s]]
            for sample in to_remove:
                sample_kept[sample] = False  # 这行其实不用，因为默认就是False
                # 只是在第一次标记时才计数
                if sample in sample_to_bl:
                    pass  # 已经被移除了
            n_removed += len(to_remove)

        if logger and level % 5 == 0:
            logger.info(f"  处理层级|Processing level parent_{level}: 已保留|kept {n_kept}, 移除|removed {n_removed}")

    # 构建结果|Build result
    result = [s for s in samples if sample_kept[s['name']]]

    if logger:
        logger.info(f"  层级去重完成|Hierarchy deduplication completed:")
        logger.info(f"    原始|Original: {len(samples)} 个样品|samples")
        logger.info(f"    保留|Kept: {len(result)} 个样品|samples")
        logger.info(f"    移除|Removed: {len(samples) - len(result)} 个样品|samples")

    return result


def adaptive_hierarchy_dedup(samples: List[Dict],
                             hierarchy_file: str,
                             target_count: int,
                             min_level: int = 15,
                             max_level: int = 38,
                             logger=None) -> List[Dict]:
    """自适应层级去重|Adaptive hierarchy deduplication

    从深层到浅层逐步去重，确保每个深层节点最多选择1个样品

    Args:
        samples: 候选样品列表（按枝长降序排序）|Candidate samples (sorted by BL descending)
        hierarchy_file: 层级关系Excel文件|Hierarchy Excel file
        target_count: 目标数量（仅供参考，去重不考虑此限制）|Target count (for reference only, dedup ignores this)
        min_level: 最浅检查层级|Minimum level to check
        max_level: 最深检查层级|Maximum level to check
        logger: 日志对象|Logger object

    Returns:
        去重后的样品列表|Deduplicated sample list
    """
    if logger:
        logger.info(f"[自适应层级去重|Adaptive Hierarchy Deduplication]")
        logger.info(f"  检查范围|Check range: parent_{min_level} ~ parent_{max_level}")
        logger.info(f"  原始数量|Original: {len(samples)}")

    # 读取层级关系|Read hierarchy
    hierarchy = pd.read_excel(hierarchy_file)
    sample_names = {s['name'] for s in samples}
    hierarchy = hierarchy[hierarchy['label'].isin(sample_names)].copy()

    sample_kept = {s['name']: True for s in samples}
    sample_to_bl = {s['name']: s['branch_length'] for s in samples}

    # 从最深层开始去重|Start deduplication from deepest level
    n_removed = 0

    for level in range(max_level, min_level - 1, -1):
        col = f'parent_{level}'
        if col not in hierarchy.columns:
            continue

        # 找出该层级有多个样品的节点|Find nodes with multiple samples at this level
        grouped = hierarchy[hierarchy[col].notna()].groupby(col)

        for node_id, group in grouped:
            # 该节点下仍然保留的样品|Samples still kept in this node
            kept_in_node = [s for s in group['label'] if sample_kept[s]]

            if len(kept_in_node) <= 1:
                continue  # 只有1个或没有，不需要去重|Only 1 or none, no need to dedup

            # 按枝长排序|Sort by branch length
            kept_in_node.sort(key=lambda x: sample_to_bl[x], reverse=True)

            # 只保留第一个|Keep only the first one
            to_remove = kept_in_node[1:]
            for sample in to_remove:
                if sample_kept[sample]:  # 只在之前保留的情况下才计数
                    sample_kept[sample] = False
                    n_removed += 1

        # 定期报告进度|Periodically report progress
        if logger and level % 5 == 0:
            current_count = sum(1 for s in samples if sample_kept[s['name']])
            logger.info(f"  层级|Level parent_{level}: 当前保留|current kept {current_count}, 已移除|removed {n_removed}")

    # 构建结果|Build result
    result = [s for s in samples if sample_kept[s['name']]]

    if logger:
        logger.info(f"  层级去重完成|Hierarchy deduplication completed:")
        logger.info(f"    最终保留|Finally kept: {len(result)} 个样品|samples")
        logger.info(f"    移除|Removed: {n_removed} 个样品|samples")

    return result


def find_clustering_level(hierarchy_file: str,
                         selected_samples: List[Dict],
                         max_cluster_threshold: int = 2,
                         min_check_level: int = 10,
                         max_check_level: int = None,
                         logger=None) -> int:
    """动态检测聚集层级|Dynamically detect clustering level

    从深层到浅层扫描，找到第一个存在聚集的层级（节点选中数>阈值）
    Scan from deep to shallow levels, find the first level with clustering (node selected count > threshold)

    Args:
        hierarchy_file: 层级关系Excel文件|Hierarchy Excel file
        selected_samples: 已选样品列表|Selected samples list
        max_cluster_threshold: 最大聚集阈值（超过此值认为存在聚集）|Max cluster threshold
        min_check_level: 最浅检查层级|Minimum level to check
        max_check_level: 最深检查层级（None=自动检测）|Maximum level to check (None=auto-detect)
        logger: 日志对象|Logger object

    Returns:
        建议的起始去重层级|Recommended deduplication start level
    """
    import pandas as pd

    hierarchy = pd.read_excel(hierarchy_file)
    selected_names = {s['name'] for s in selected_samples}
    hierarchy_selected = hierarchy[hierarchy['label'].isin(selected_names)].copy()

    # 自动检测最大深度|Auto-detect maximum depth
    if max_check_level is None:
        max_check_level = 0
        for col in hierarchy.columns:
            if col.startswith('parent_'):
                level = int(col.split('_')[1])
                if level > max_check_level:
                    max_check_level = level
        if logger:
            logger.info(f"  自动检测最大深度|Auto-detected max depth: parent_{max_check_level}")

    if logger:
        logger.info(f"[动态检测聚集层级|Dynamic Clustering Level Detection]")
        logger.info(f"  检查范围|Check range: parent_{min_check_level} ~ parent_{max_check_level}")
        logger.info(f"  聚集阈值|Cluster threshold: >{max_cluster_threshold} 个选中|selected")

    # 从深层到浅层扫描|Scan from deep to shallow
    for level in range(max_check_level, min_check_level - 1, -1):
        col = f'parent_{level}'
        if col not in hierarchy_selected.columns:
            continue

        # 统计该层级每个节点的选中数|Count selected samples per node at this level
        grouped = hierarchy_selected[hierarchy_selected[col].notna()].groupby(col)
        max_selected = 0
        max_node = None

        for node_id, group in grouped:
            n_selected = len(group)
            if n_selected > max_selected:
                max_selected = n_selected
                max_node = node_id

        if max_selected > max_cluster_threshold:
            if logger:
                logger.info(f"  发现聚集|Found clustering at parent_{level}:")
                logger.info(f"    最大聚集|Max cluster: {max_selected} 个选中|selected (超过阈值|exceeds threshold {max_cluster_threshold})")
            return level

    # 如果没有发现聚集，返回min_check_level
    if logger:
        logger.info(f"  未发现聚集|No clustering found，使用默认层级|using default level: parent_{min_check_level}")

    return min_check_level


def adaptive_hierarchy_dedup_with_replacement(selected_samples: List[Dict],
                                               candidate_pool: List[Dict],
                                               hierarchy_file: str,
                                               target_count: int,
                                               min_level: int = 20,
                                               max_level: int = 38,
                                               sample_to_pca: Dict[str, List[float]] = None,
                                               tree_analyzer = None,
                                               logger=None) -> Tuple[List[Dict], List[str]]:
    """层级去重 + 替换|Hierarchy dedup with replacement

    1. 从selected_samples中移除在深层节点聚集的样品
    2. 从candidate_pool中选择新样品来补充，确保最终数量=target_count

    Args:
        selected_samples: 已选样品列表|Selected samples
        candidate_pool: 候选池（未选中的样品）|Candidate pool (unselected samples)
        hierarchy_file: 层级关系Excel文件|Hierarchy Excel file
        target_count: 目标数量|Target count
        min_level: 最浅检查层级|Minimum level to check
        max_level: 最深检查层级|Maximum level to check
        sample_to_pca: PCA坐标|PCA coordinates
        tree_analyzer: 树分析器|Tree analyzer
        logger: 日志对象|Logger object

    Returns:
        (去重并替换后的样品列表, 被替换的样品名列表)|(Deduplicated and replaced samples, list of replaced sample names)
    """
    if logger:
        logger.info(f"[层级去重+替换|Hierarchy Dedup + Replacement]")
        logger.info(f"  检查范围|Check range: parent_{min_level} ~ parent_{max_level}")
        logger.info(f"  已选数量|Selected: {len(selected_samples)}")
        logger.info(f"  候选池|Candidate pool: {len(candidate_pool)}")
        logger.info(f"  目标数量|Target: {target_count}")

    # 读取层级关系|Read hierarchy
    import pandas as pd
    hierarchy = pd.read_excel(hierarchy_file)

    selected_names = {s['name'] for s in selected_samples}
    candidate_names = {s['name'] for s in candidate_pool}

    # 第1步：层级去重|Step 1: Hierarchy dedup
    sample_kept = {s['name']: True for s in selected_samples}
    sample_to_bl = {s['name']: s['branch_length'] for s in selected_samples}
    removed_samples = []

    # 只在selected_samples中筛选
    hierarchy_selected = hierarchy[hierarchy['label'].isin(selected_names)].copy()

    for level in range(max_level, min_level - 1, -1):
        col = f'parent_{level}'
        if col not in hierarchy_selected.columns:
            continue

        grouped = hierarchy_selected[hierarchy_selected[col].notna()].groupby(col)

        for node_id, group in grouped:
            kept_in_node = [s for s in group['label'] if sample_kept.get(s, False)]

            if len(kept_in_node) <= 1:
                continue

            # 按枝长排序，保留最大的|Sort by BL, keep the largest
            kept_in_node.sort(key=lambda x: sample_to_bl[x], reverse=True)

            # 移除其他的|Remove others
            to_remove = kept_in_node[1:]
            for sample in to_remove:
                if sample_kept.get(sample, False):
                    sample_kept[sample] = False
                    removed_samples.append(sample)

    # 构建去重后的结果|Build deduplicated result
    deduplicated = [s for s in selected_samples if sample_kept[s['name']]]

    if logger:
        logger.info(f"  层级去重|Hierarchy dedup:")
        logger.info(f"    移除|Removed: {len(removed_samples)} 个聚集样品|clustered samples")

    # 第2步：计算需要补充的数量|Step 2: Calculate how many to add
    n_to_add = target_count - len(deduplicated)

    if n_to_add <= 0:
        if logger:
            logger.info(f"  已满足目标，无需补充|Target already met, no replacement needed")
        return deduplicated, []

    if logger:
        logger.info(f"  需要补充|Need to add: {n_to_add} 个样品|samples")

    # 第3步：从候选池中选择新样品|Step 3: Select new samples from candidate pool
    # 排除已选中的|Exclude already selected
    available_candidates = [s for s in candidate_pool if s['name'] not in selected_names]

    if len(available_candidates) == 0:
        if logger:
            logger.warning(f"  候选池为空，无法补充|Candidate pool empty, cannot replace")
        return deduplicated, removed_samples

    # 策略：贪心选择，每次选择对已选集合PD贡献最大的候选
    # Strategy: Greedy selection, select candidate with max PD contribution
    final_selected = list(deduplicated)
    current_names = {s['name'] for s in final_selected}

    # 构建候选到层级的映射|Build candidate to hierarchy mapping
    hierarchy_candidates = hierarchy[hierarchy['label'].isin(candidate_names)].copy()

    for iteration in range(n_to_add):
        if len(available_candidates) == 0:
            break

        best_candidate = None
        best_contribution = -float('inf')

        for candidate in available_candidates:
            # 检查是否会在深层节点上与已选样品聚集|Check if will cluster with selected at deep level
            candidate_row = hierarchy_candidates[hierarchy_candidates['label'] == candidate['name']]

            if len(candidate_row) == 0:
                continue

            will_cluster = False
            for level in range(max_level, min_level - 1, -1):
                col = f'parent_{level}'
                if col not in candidate_row.columns:
                    continue

                node_id = candidate_row[col].values[0]
                if pd.isna(node_id):
                    continue

                # 检查该节点是否已有选中样品|Check if node already has selected sample
                node_samples = hierarchy_selected[hierarchy_selected[col] == node_id]['label'].values
                if any(s in current_names for s in node_samples):
                    will_cluster = True
                    break

            if will_cluster:
                continue  # 跳过会聚集的候选|Skip candidates that will cluster

            # 计算PD贡献|Calculate PD contribution
            if len(final_selected) == 0:
                contribution = candidate['branch_length']
            else:
                # 计算到所有已选样品的最小Cophenetic距离|Min Cophenetic distance to selected
                min_dist = float('inf')
                for selected_s in final_selected:
                    if tree_analyzer:
                        dist = tree_analyzer.get_tree_distance(candidate['name'], selected_s['name'])
                        if dist is not None and dist < min_dist:
                            min_dist = dist
                contribution = min_dist

            if contribution > best_contribution:
                best_contribution = contribution
                best_candidate = candidate

        if best_candidate:
            final_selected.append(best_candidate)
            current_names.add(best_candidate['name'])
            available_candidates.remove(best_candidate)

            if logger and (iteration + 1) % 10 == 0:
                logger.info(f"    进度|Progress: {iteration + 1}/{n_to_add} 已补充|added")

    # 第4步：如果还是不够，放宽约束|Step 4: If still not enough, relax constraints
    if len(final_selected) < target_count and len(available_candidates) > 0:
        n_still_needed = target_count - len(final_selected)
        if logger:
            logger.info(f"  仍需要|Still need: {n_still_needed} 个，放宽层级约束|samples, relaxing hierarchy constraints")

        # 放宽到parent_15|Relax to parent_15
        relaxed_min_level = 15

        for _ in range(n_still_needed):
            if len(available_candidates) == 0:
                break

            best_candidate = None
            best_contribution = -float('inf')

            for candidate in available_candidates:
                # 检查是否会在更深层节点上聚集|Check clustering at deeper levels
                candidate_row = hierarchy_candidates[hierarchy_candidates['label'] == candidate['name']]

                if len(candidate_row) == 0:
                    continue

                will_cluster = False
                for level in range(max_level, relaxed_min_level - 1, -1):
                    col = f'parent_{level}'
                    if col not in candidate_row.columns:
                        continue

                    node_id = candidate_row[col].values[0]
                    if pd.isna(node_id):
                        continue

                    node_samples = hierarchy_selected[hierarchy_selected[col] == node_id]['label'].values
                    if any(s in current_names for s in node_samples):
                        will_cluster = True
                        break

                if will_cluster:
                    continue

                # 计算PD贡献|Calculate PD contribution
                min_dist = float('inf')
                for selected_s in final_selected:
                    if tree_analyzer:
                        dist = tree_analyzer.get_tree_distance(candidate['name'], selected_s['name'])
                        if dist is not None and dist < min_dist:
                            min_dist = dist
                contribution = min_dist

                if contribution > best_contribution:
                    best_contribution = contribution
                    best_candidate = candidate

            if best_candidate:
                final_selected.append(best_candidate)
                current_names.add(best_candidate['name'])
                available_candidates.remove(best_candidate)

    replaced = removed_samples[:len(final_selected) - len(deduplicated)]

    if logger:
        logger.info(f"  替换完成|Replacement completed:")
        logger.info(f"    最终保留|Finally kept: {len(final_selected)} 个样品|samples")
        logger.info(f"    移除|Removed: {len(removed_samples)} 个聚集样品|clustered samples")
        logger.info(f"    新增|Added: {len(final_selected) - len(deduplicated)} 个新样品|new samples")
        logger.info(f"    被替换|Replaced: {len(replaced)} 个样品|samples")

    return final_selected, replaced

