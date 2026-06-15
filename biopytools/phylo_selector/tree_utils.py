"""
进化树结构分析工具|Phylogenetic Tree Structure Analysis Tools
"""

from Bio import Phylo
import numpy as np
from typing import Dict, List, Tuple, Set
from collections import defaultdict


class TreeStructureAnalyzer:
    """进化树结构分析器|Phylogenetic Tree Structure Analyzer"""

    def __init__(self, newick_file: str, logger=None):
        """初始化树分析器|Initialize tree analyzer

        Args:
            newick_file: Newick格式树文件|Newick format tree file
            logger: 日志记录器|Logger
        """
        self.newick_file = newick_file
        self.logger = logger

        # 解析树|Parse tree
        self.tree = Phylo.read(newick_file, "newick")

        # 获取所有叶节点|Get all terminals
        self.terminals = self.tree.get_terminals()
        self.terminal_map = {t.name: t for t in self.terminals if t.name}

        if self.logger:
            self.logger.info(f"  加载进化树|Loaded phylogenetic tree: {len(self.terminals)} 个叶节点|terminals")
            self.logger.info(f"  总节点数|Total nodes: {len(self.tree.get_nonterminals()) + len(self.terminals)}")

    def get_tree_distance(self, sample1: str, sample2: str) -> float:
        """计算两个样品之间的树距离（进化距离）|Calculate tree distance between two samples

        树距离 = 从样品1到样品2沿树的最短路径长度
        Tree distance = shortest path length from sample1 to sample2 along the tree

        Args:
            sample1: 样品1名称|Sample 1 name
            sample2: 样品2名称|Sample 2 name

        Returns:
            树距离|Tree distance
        """
        if sample1 not in self.terminal_map or sample2 not in self.terminal_map:
            return None

        t1 = self.terminal_map[sample1]
        t2 = self.terminal_map[sample2]

        # 计算树距离：distance(t1) + distance(t2) - 2*distance(common_ancestor)
        tree_dist = self.tree.distance(t1) + self.tree.distance(t2)
        common_anc = self.tree.common_ancestor(t1, t2)

        if common_anc:
            tree_dist -= 2 * self.tree.distance(common_anc)

        return tree_dist

    def get_path_to_root(self, sample: str) -> List[str]:
        """获取样品到根节点的路径上的所有节点名称|Get path from sample to root

        Args:
            sample: 样品名称|Sample name

        Returns:
            路径上的节点名称列表（从样品到根）|List of node names on path (sample to root)
        """
        if sample not in self.terminal_map:
            return None

        terminal = self.terminal_map[sample]
        path = self.tree.get_path(terminal)

        # 获取路径上所有节点的标识（用clades或创建ID）
        path_ids = []
        for i, node in enumerate(path):
            if node.name:
                path_ids.append(node.name)
            else:
                # 对于没有名字的内部节点，用其位置创建唯一ID
                path_ids.append(f"internal_{i}_{id(node)}")

        return path_ids

    def find_common_ancestor_depth(self, sample1: str, sample2: str) -> float:
        """找到两个样品的共同祖先的深度（到根节点的距离）|Find depth of common ancestor

        Args:
            sample1: 样品1名称|Sample 1 name
            sample2: 样品2名称|Sample 2 name

        Returns:
            共同祖先到根节点的距离|Distance from common ancestor to root
        """
        if sample1 not in self.terminal_map or sample2 not in self.terminal_map:
            return None

        t1 = self.terminal_map[sample1]
        t2 = self.terminal_map[sample2]
        common_anc = self.tree.common_ancestor(t1, t2)

        if common_anc:
            return self.tree.distance(common_anc)
        return 0.0

    def calculate_branch_length(self, sample: str) -> float:
        """获取样品的枝长（到根节点的距离）|Get branch length of sample

        Args:
            sample: 样品名称|Sample name

        Returns:
            枝长|Branch length
        """
        if sample not in self.terminal_map:
            return None

        return self.tree.distance(self.terminal_map[sample])

    def calculate_all_tree_distances(self, samples: List[str]) -> Dict[Tuple[str, str], float]:
        """计算所有样品对之间的树距离|Calculate tree distances for all sample pairs

        Args:
            samples: 样品列表|Sample list

        Returns:
            字典：{(sample1, sample2): tree_distance}|Dict: {(sample1, sample2): tree_distance}
        """
        distances = {}
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                dist = self.get_tree_distance(samples[i], samples[j])
                if dist is not None:
                    distances[(samples[i], samples[j])] = dist
        return distances

    def find_closest_pairs_by_tree_distance(self, samples: List[str], n_pairs: int = 10) -> List[Tuple[str, str, float]]:
        """找出树距离最近的样品对|Find sample pairs with smallest tree distances

        Args:
            samples: 样品列表|Sample list
            n_pairs: 返回的对数|Number of pairs to return

        Returns:
            [(sample1, sample2, tree_distance), ...]|List of tuples
        """
        all_distances = self.calculate_all_tree_distances(samples)
        sorted_pairs = sorted(all_distances.items(), key=lambda x: x[1])

        return [(k[0], k[1], v) for k, v in sorted_pairs[:n_pairs]]

    def deduplicate_by_tree_distance(self, samples: List[Dict], threshold: float = 0.01) -> Tuple[List[Dict], List[Dict]]:
        """基于树距离去重|Deduplicate samples based on tree distance

        移除树距离小于阈值的样品对中的一个（保留枝长较长的）|Remove one sample from pairs with tree distance < threshold (keep longer BL)

        Args:
            samples: 样品列表，每个样品包含'name'和'branch_length'|Sample list with 'name' and 'branch_length'
            threshold: 树距离阈值|Tree distance threshold

        Returns:
            (保留的样品, 移除的样品)|(kept_samples, removed_samples)
        """
        sample_names = [s['name'] for s in samples]

        # 计算所有样品对的树距离
        all_distances = self.calculate_all_tree_distances(sample_names)

        # 找出距离小于阈值的样品对
        close_pairs = [(s1, s2, dist) for (s1, s2), dist in all_distances.items() if dist < threshold]

        # 按树距离排序，从小到大
        close_pairs.sort(key=lambda x: x[2])

        to_remove = set()
        removed_samples = []

        for s1_name, s2_name, tree_dist in close_pairs:
            if s1_name in to_remove or s2_name in to_remove:
                continue

            # 找到对应的样品对象
            s1 = next(s for s in samples if s['name'] == s1_name)
            s2 = next(s for s in samples if s['name'] == s2_name)

            # 移除枝长较短的一个
            if s1['branch_length'] >= s2['branch_length']:
                to_remove.add(s2_name)
                removed_samples.append(s2)
            else:
                to_remove.add(s1_name)
                removed_samples.append(s1)

        # 保留未被移除的样品
        kept_samples = [s for s in samples if s['name'] not in to_remove]

        if self.logger:
            self.logger.info(f"    基于树距离去重|Deduplication by tree distance:")
            self.logger.info(f"      阈值|Threshold: {threshold}")
            self.logger.info(f"      移除|Removed: {len(removed_samples)} 个样品|samples")
            self.logger.info(f"      保留|Kept: {len(kept_samples)} 个样品|samples")

        return kept_samples, removed_samples

    def get_subtree_samples(self, sample: str, max_distance: float = 0.01) -> Set[str]:
        """获取样品附近子树中的所有样品|Get all samples in the subtree near a sample

        找出与给定样品的树距离小于max_distance的所有样品|Find all samples with tree distance < max_distance

        Args:
            sample: 样品名称|Sample name
            max_distance: 最大树距离|Maximum tree distance

        Returns:
            附近样品的名称集合|Set of nearby sample names
        """
        nearby = set()

        if sample not in self.terminal_map:
            return nearby

        for other_name, other_terminal in self.terminal_map.items():
            if other_name == sample:
                continue

            dist = self.get_tree_distance(sample, other_name)
            if dist is not None and dist < max_distance:
                nearby.add(other_name)

        return nearby

    def calculate_tree_spread(self, samples: List[str]) -> Dict[str, float]:
        """计算一组样品在树上的分布广度|Calculate spread of samples on tree

        Args:
            samples: 样品列表|Sample list

        Returns:
            统计信息字典|Statistics dictionary
        """
        if len(samples) < 2:
            return {}

        # 计算所有样品对的树距离
        all_distances = self.calculate_all_tree_distances(samples)
        distances = list(all_distances.values())

        return {
            'min_tree_distance': min(distances),
            'max_tree_distance': max(distances),
            'mean_tree_distance': np.mean(distances),
            'median_tree_distance': np.median(distances),
            'std_tree_distance': np.std(distances)
        }

    # ========== 层级关系提取|Hierarchy Extraction ==========

    def _build_node_hierarchy(self, node, parent_id=None, depth=0, node_counter=None, node_map=None):
        """递归构建节点层级关系|Recursively build node hierarchy

        Args:
            node: Bio.Phylo Clade对象|Bio.Phylo Clade object
            parent_id: 父节点ID|Parent node ID
            depth: 当前深度|Current depth
            node_counter: 节点计数器字典|Node counter dictionary
            node_map: 节点映射|Node mapping

        Returns:
            当前节点ID|Current node ID
        """
        if node_counter is None:
            node_counter = {'count': 0}
        if node_map is None:
            node_map = {}

        # 为当前节点分配ID|Assign ID to current node
        current_id = node_counter['count']
        node_counter['count'] += 1

        # 存储节点信息|Store node information
        node_map[current_id] = {
            'parent': parent_id,
            'depth': depth,
            'clade': node,
            'name': node.name if hasattr(node, 'name') and node.name else None
        }

        # 递归处理子节点|Recursively process children
        for child in node.clades:
            self._build_node_hierarchy(child, current_id, depth + 1, node_counter, node_map)

        return current_id

    def extract_hierarchy_dict(self, max_depth: int = None) -> Dict[str, List[int]]:
        """提取层级关系字典|Extract hierarchy dictionary

        返回格式|Return format: {sample_name: [node_id_1, node_id_2, ...]}

        Args:
            max_depth: 最大深度（可选）|Maximum depth (optional)

        Returns:
            样品名到节点路径的映射|Sample name to node path mapping
        """
        # 构建节点层级|Build node hierarchy
        node_map = {}
        node_counter = {'count': 0}
        self._build_node_hierarchy(self.tree.root, node_map=node_map, node_counter=node_counter)

        # 为每个叶节点构建路径|Build path for each terminal
        hierarchy = {}
        for terminal in self.terminals:
            if not terminal.name:
                continue

            # 找到叶节点的ID|Find terminal ID
            terminal_id = None
            for node_id, info in node_map.items():
                if info['clade'] == terminal:
                    terminal_id = node_id
                    break

            if terminal_id is None:
                continue

            # 追踪路径到根|Trace path to root
            path = []
            current_id = terminal_id
            while current_id is not None:
                path.append(current_id)
                current_id = node_map[current_id]['parent']

            # 反转路径（从根到叶）|Reverse path (root to leaf)
            path = list(reversed(path))

            # 应用最大深度限制|Apply max depth limit
            if max_depth is not None and len(path) > max_depth:
                path = path[:max_depth]

            hierarchy[terminal.name] = path

        return hierarchy

    def extract_hierarchy_dataframe(self, max_depth: int = None):
        """提取层级关系DataFrame|Extract hierarchy DataFrame

        返回格式类似于R函数的输出|Return format similar to R function output:
        parent_1, parent_2, ..., parent_N, label

        Args:
            max_depth: 最大深度（可选）|Maximum depth (optional)

        Returns:
            pandas DataFrame
        """
        try:
            import pandas as pd
        except ImportError:
            if self.logger:
                self.logger.error("pandas未安装|pandas not installed")
            return None

        hierarchy = self.extract_hierarchy_dict(max_depth=max_depth)

        if not hierarchy:
            return None

        # 确定最大深度|Determine max depth
        if max_depth is None:
            max_depth = max(len(path) for path in hierarchy.values())
        else:
            max_depth = min(max_depth, max(len(path) for path in hierarchy.values()))

        # 构建DataFrame|Build DataFrame
        data = []
        for label, path in hierarchy.items():
            row = {'label': label}
            for i in range(max_depth):
                col_name = f'parent_{i+1}'
                if i < len(path):
                    row[col_name] = path[i]
                else:
                    row[col_name] = None
            data.append(row)

        df = pd.DataFrame(data)

        # 调整列顺序|Adjust column order
        parent_cols = [f'parent_{i+1}' for i in range(max_depth)]
        df = df[parent_cols + ['label']]

        return df

    def calculate_topology_distance(self, sample1: str, sample2: str) -> int:
        """计算两个样品之间的拓扑距离|Calculate topology distance between two samples

        拓扑距离 = depth1 + depth2 - 2 * lca_depth
        Topology distance = depth1 + depth2 - 2 * lca_depth

        这是节点数距离，不包含枝长|This is node count distance, does not include branch length

        Args:
            sample1: 样品1名称|Sample 1 name
            sample2: 样品2名称|Sample 2 name

        Returns:
            拓扑距离（整数）|Topology distance (integer)
        """
        if sample1 not in self.terminal_map or sample2 not in self.terminal_map:
            return float('inf')

        # 获取路径|Get paths
        path1 = self.get_path_to_root(sample1)
        path2 = self.get_path_to_root(sample2)

        if path1 is None or path2 is None:
            return float('inf')

        # 找到LCA（最近共同祖先）|Find LCA (Lowest Common Ancestor)
        lca_depth = 0
        for node1, node2 in zip(path1, path2):
            if node1 == node2:
                lca_depth += 1
            else:
                break

        # 计算距离|Calculate distance
        distance = len(path1) + len(path2) - 2 * lca_depth

        return distance
