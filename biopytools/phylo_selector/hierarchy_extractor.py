#!/usr/bin/env python3
"""
系统发育树层级关系提取器|Phylogenetic Tree Hierarchy Extractor

从Newick格式树文件中提取层级关系，生成父子路径表
Extract parent-child path hierarchy from Newick format tree files

Author: Claude Code
Date: 2026-01-30
"""

import os
import json
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

try:
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Clade
except ImportError:
    Phylo = None
    Clade = None

import numpy as np


class HierarchyExtractor:
    """系统发育树层级关系提取器|Hierarchy Extractor for Phylogenetic Trees"""

    def __init__(self, logger=None):
        """
        初始化|Initialize

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger
        self.tree = None
        self.node_to_parent = {}  # 节点ID -> 父节点ID|node_id -> parent_id
        self.node_to_depth = {}   # 节点ID -> 深度|node_id -> depth
        self.label_to_path = {}   # 样品名 -> 路径|sample_name -> path
        self.node_counter = 0     # 节点计数器|Node counter

    def load_tree(self, newick_file: str) -> bool:
        """
        加载Newick格式树文件|Load Newick format tree file

        Args:
            newick_file: Newick文件路径|Path to Newick file

        Returns:
            是否成功|Success status
        """
        if Phylo is None:
            if self.logger:
                self.logger.error("Bio.Phylo未安装|Bio.Phylo not installed")
            return False

        if not os.path.exists(newick_file):
            if self.logger:
                self.logger.error(f"树文件不存在|Tree file not found: {newick_file}")
            return False

        try:
            self.tree = Phylo.read(newick_file, 'newick')
            if self.logger:
                self.logger.info(f"成功加载树文件|Tree loaded successfully: {newick_file}")
            return True
        except Exception as e:
            if self.logger:
                self.logger.error(f"加载树文件失败|Failed to load tree: {e}")
            return False

    def _assign_node_ids(self, clade: Clade, parent_id: Optional[int] = None, depth: int = 0):
        """
        递归分配节点ID|Recursively assign node IDs

        Args:
            clade: Bio.Phylo Clade对象|Bio.Phylo Clade object
            parent_id: 父节点ID|Parent node ID
            depth: 当前深度|Current depth
        """
        current_id = self.node_counter
        self.node_counter += 1

        # 存储父子关系|Store parent-child relationship
        self.node_to_parent[current_id] = parent_id
        self.node_to_depth[current_id] = depth

        # 如果有标签，存储标签到ID的映射|If has label, store label to ID mapping
        if clade.name:
            self.label_to_node[clade.name] = current_id

        # 递归处理子节点|Recursively process children
        for child in clade.clades:
            self._assign_node_ids(child, current_id, depth + 1)

    def _build_hierarchy_map(self):
        """构建层级映射|Build hierarchy mapping"""
        if self.tree is None:
            return False

        # 初始化|Initialize
        self.node_to_parent = {}
        self.node_to_depth = {}
        self.label_to_node = {}
        self.node_counter = 0
        self.label_to_path = {}

        # 分配节点ID|Assign node IDs
        if self.tree.root is not None:
            self._assign_node_ids(self.tree.root)

        if self.logger:
            self.logger.info(f"分配了 {self.node_counter} 个节点|Assigned {self.node_counter} nodes")

        # 为每个标签构建路径|Build path for each label
        for label, node_id in self.label_to_node.items():
            path = []
            current_id = node_id
            while current_id is not None:
                path.append(current_id)
                current_id = self.node_to_parent.get(current_id)
            # 路径从根到叶|Path from root to leaf
            path = list(reversed(path))
            self.label_to_path[label] = path

        if self.logger:
            self.logger.info(f"构建了 {len(self.label_to_path)} 个标签的路径|Built paths for {len(self.label_to_path)} labels")

        return True

    def extract_hierarchy(self, max_depth: Optional[int] = None) -> Dict[str, List[int]]:
        """
        提取层级关系|Extract hierarchy relationships

        Args:
            max_depth: 最大深度（可选）|Maximum depth (optional)

        Returns:
            标签到路径的映射|Label to path mapping
        """
        if self.tree is None:
            if self.logger:
                self.logger.error("树未加载|Tree not loaded")
            return {}

        # 构建层级映射|Build hierarchy mapping
        if not self._build_hierarchy_map():
            return {}

        # 如果指定了最大深度，截断路径|If max_depth specified, truncate paths
        if max_depth is not None:
            for label in self.label_to_path:
                if len(self.label_to_path[label]) > max_depth:
                    self.label_to_path[label] = self.label_to_path[label][:max_depth]

        return self.label_to_path

    def export_to_dataframe(self, max_depth: Optional[int] = None):
        """
        导出为pandas DataFrame格式|Export to pandas DataFrame format

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

        # 提取层级关系|Extract hierarchy
        hierarchy = self.extract_hierarchy(max_depth)

        if not hierarchy:
            return None

        # 确定最大深度|Determine max depth
        if max_depth is None:
            max_depth = max(len(path) for path in hierarchy.values())
        else:
            max_depth = min(max_depth, max(len(path) for path in hierarchy.values()))

        # 构建DataFrame|Build DataFrame
        data = {}
        for label, path in hierarchy.items():
            row = {}
            for i in range(max_depth):
                col_name = f'parent_{i+1}'
                if i < len(path):
                    row[col_name] = path[i]
                else:
                    row[col_name] = None
            row['label'] = label
            data[label] = row

        # 转换为DataFrame|Convert to DataFrame
        df_data = []
        for label, row in data.items():
            df_data.append(row)

        df = pd.DataFrame(df_data)

        # 调整列顺序|Adjust column order
        parent_cols = [f'parent_{i+1}' for i in range(max_depth)]
        df = df[parent_cols + ['label']]

        return df

    def export_to_excel(self, output_file: str, max_depth: Optional[int] = None) -> bool:
        """
        导出为Excel文件|Export to Excel file

        Args:
            output_file: 输出文件路径|Output file path
            max_depth: 最大深度（可选）|Maximum depth (optional)

        Returns:
            是否成功|Success status
        """
        try:
            import pandas as pd
        except ImportError:
            if self.logger:
                self.logger.error("pandas未安装|pandas not installed")
            return False

        df = self.export_to_dataframe(max_depth)

        if df is None:
            return False

        try:
            df.to_excel(output_file, index=False)
            if self.logger:
                self.logger.info(f"成功导出到|Successfully exported to: {output_file}")
            return True
        except Exception as e:
            if self.logger:
                self.logger.error(f"导出失败|Export failed: {e}")
            return False

    def calculate_topology_distance(self, label1: str, label2: str) -> float:
        """
        计算两个标签之间的拓扑距离|Calculate topology distance between two labels

        拓扑距离 = depth1 + depth2 - 2 * lca_depth
        Topology distance = depth1 + depth2 - 2 * lca_depth

        Args:
            label1: 标签1|Label 1
            label2: 标签2|Label 2

        Returns:
            拓扑距离|Topology distance
        """
        if label1 not in self.label_to_path or label2 not in self.label_to_path:
            return float('inf')

        path1 = self.label_to_path[label1]
        path2 = self.label_to_path[label2]

        # 找到LCA（最近共同祖先）|Find LCA (Lowest Common Ancestor)
        lca_depth = 0
        for i, (node1, node2) in enumerate(zip(path1, path2)):
            if node1 == node2:
                lca_depth = i + 1
            else:
                break

        # 计算距离|Calculate distance
        distance = len(path1) + len(path2) - 2 * lca_depth

        return distance

    def calculate_topology_distance_matrix(self, labels: List[str]) -> Dict[str, Dict[str, float]]:
        """
        计算标签间的拓扑距离矩阵|Calculate topology distance matrix

        Args:
            labels: 标签列表|List of labels

        Returns:
            距离矩阵|Distance matrix
        """
        matrix = {}
        for label1 in labels:
            matrix[label1] = {}
            for label2 in labels:
                if label1 == label2:
                    matrix[label1][label2] = 0.0
                else:
                    matrix[label1][label2] = self.calculate_topology_distance(label1, label2)
        return matrix


def extract_hierarchy_from_newick(newick_file: str,
                                  output_file: Optional[str] = None,
                                  max_depth: Optional[int] = None,
                                  logger=None) -> Optional[Dict]:
    """
    便捷函数：从Newick文件提取层级关系|Convenience function: Extract hierarchy from Newick file

    Args:
        newick_file: Newick文件路径|Path to Newick file
        output_file: 输出Excel文件路径（可选）|Output Excel file path (optional)
        max_depth: 最大深度（可选）|Maximum depth (optional)
        logger: 日志对象|Logger object

    Returns:
        标签到路径的映射|Label to path mapping
    """
    extractor = HierarchyExtractor(logger=logger)

    if not extractor.load_tree(newick_file):
        return None

    hierarchy = extractor.extract_hierarchy(max_depth=max_depth)

    if output_file:
        extractor.export_to_excel(output_file, max_depth=max_depth)

    return hierarchy


# ========== 命令行接口|Command Line Interface ==========
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='从Newick树提取层级关系|Extract hierarchy from Newick tree',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 从Newick文件提取层级并保存为Excel
  python hierarchy_extractor.py -i tree.nwk -o hierarchy.xlsx

  # 限制最大深度
  python hierarchy_extractor.py -i tree.nwk -o hierarchy.xlsx --max-depth 20

  # 输出到JSON
  python hierarchy_extractor.py -i tree.nwk -o hierarchy.json --format json
        '''
    )

    parser.add_argument('-i', '--input', required=True,
                       help='Newick树文件|Newick tree file')
    parser.add_argument('-o', '--output', required=True,
                       help='输出文件路径|Output file path')
    parser.add_argument('--max-depth', type=int, default=None,
                       help='最大深度|Maximum depth (default: auto)')
    parser.add_argument('--format', choices=['excel', 'json'], default='excel',
                       help='输出格式|Output format (default: excel)')

    args = parser.parse_args()

    # 提取层级关系|Extract hierarchy
    extractor = HierarchyExtractor(logger=None)

    if not extractor.load_tree(args.input):
        print(f"错误：无法加载树文件|Error: Failed to load tree file: {args.input}")
        return 1

    if args.format == 'excel':
        success = extractor.export_to_excel(args.output, max_depth=args.max_depth)
    else:  # json
        import json
        hierarchy = extractor.extract_hierarchy(max_depth=args.max_depth)
        if hierarchy:
            with open(args.output, 'w') as f:
                json.dump(hierarchy, f, indent=2)
            success = True
        else:
            success = False

    if success:
        print(f"成功导出到|Successfully exported to: {args.output}")
        return 0
    else:
        print(f"导出失败|Export failed")
        return 1


if __name__ == '__main__':
    main()
