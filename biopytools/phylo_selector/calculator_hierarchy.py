"""
基于层级关系的智能选择算法|Hierarchy-based Intelligent Selection Algorithm
解决按行顺序选择导致的"树紧挨着"问题
"""

import pandas as pd
from collections import defaultdict
from typing import List, Dict


def hierarchy_based_selection(samples: List[Dict],
                              hierarchy_df: pd.DataFrame,
                              target_count: int,
                              hierarchy_level: int = 10,
                              logger=None) -> List[Dict]:
    """基于层级关系的选择|Hierarchy-based selection

    按照指定的层级深度（如parent_10）进行分组，从每个分组中按比例选择样品
    这样可以确保选出的样品在进化树的不同分支上均匀分布

    Args:
        samples: 样品列表|Sample list
        hierarchy_df: 层级关系DataFrame|Hierarchy DataFrame
        target_count: 目标数量|Target count
        hierarchy_level: 使用的层级深度|Hierarchy level to use (default: 10)
        logger: 日志器|Logger

    Returns:
        List[Dict]: 选择的样品|Selected samples
    """
    # 创建样品名到样品对象的映射
    sample_map = {s['name']: s for s in samples}

    # 找到样品在层级文件中的行
    label_to_idx = {label: idx for idx, label in enumerate(hierarchy_df['label'])}

    # 只处理在层级文件中的样品
    valid_samples = [s for s in samples if s['name'] in label_to_idx]

    if len(valid_samples) == 0:
        if logger:
            logger.warning(f"  没有样品在层级文件中|No samples in hierarchy file")
        return samples[:target_count]

    # 按照指定层级分组
    col_name = f'parent_{hierarchy_level}'
    group_to_samples = defaultdict(list)

    for sample in valid_samples:
        idx = label_to_idx[sample['name']]
        row = hierarchy_df.iloc[idx]
        parent_id = row[col_name]

        if pd.isna(parent_id):
            # 如果该层级为空，尝试更浅的层级
            for level in range(hierarchy_level - 1, 0, -1):
                col = f'parent_{level}'
                val = row[col]
                if not pd.isna(val):
                    parent_id = val
                    break

        if not pd.isna(parent_id):
            group_to_samples[parent_id].append(sample)

    # 计算每个分组应该选择多少个样品（按比例分配）
    total_valid = len(valid_samples)
    group_quotas = {}
    remaining = target_count

    # 先给每个分组至少分配1个
    for group_id, group_samples in group_to_samples.items():
        group_quotas[group_id] = 1
        remaining -= 1

    # 按比例分配剩余名额
    if remaining > 0:
        for group_id, group_samples in group_to_samples.items():
            additional_quota = int(len(group_samples) / total_valid * remaining)
            group_quotas[group_id] += additional_quota

    # 调整到目标数量
    current_total = sum(group_quotas.values())
    if current_total < target_count:
        # 如果还不够，分配给样品最多的组
        sorted_groups = sorted(group_to_samples.items(),
                             key=lambda x: len(x[1]),
                             reverse=True)
        for i in range(target_count - current_total):
            if i < len(sorted_groups):
                group_id = sorted_groups[i][0]
                group_quotas[group_id] += 1
    elif current_total > target_count:
        # 如果太多了，从样品最多的组减少
        sorted_groups = sorted(group_to_samples.items(),
                             key=lambda x: len(x[1]),
                             reverse=True)
        for i in range(current_total - target_count):
            if i < len(sorted_groups):
                group_id = sorted_groups[i][0]
                if group_quotas[group_id] > 1:
                    group_quotas[group_id] -= 1

    if logger:
        logger.info(f"  层级|Level: parent_{hierarchy_level}")
        logger.info(f"  分组数|Groups: {len(group_to_samples)}")
        logger.info(f"  分组大小|Group sizes: {sorted([len(g) for g in group_to_samples.values()])[:10]}...")
        logger.info(f"  分组配额|Group quotas: {sorted(list(group_quotas.values()))[:10]}...")

    # 从每个分组中选择样品（均匀选择）
    selected = []
    for group_id, quota in group_quotas.items():
        group_samples = group_to_samples[group_id]
        group_sample_count = len(group_samples)

        if group_sample_count <= quota:
            # 如果分组样品数<=配额，全选
            selected.extend(group_samples)
        else:
            # 均匀选择
            step = group_sample_count / quota
            for i in range(quota):
                idx = int(i * step)
                selected.append(group_samples[idx])

    if logger:
        logger.info(f"  层级选择完成|Hierarchy selection completed: {len(selected)} 个样品|samples")

    return selected[:target_count]
