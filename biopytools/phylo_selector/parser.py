"""
系统发育树样品选择解析模块|Phylogenetic Tree Sample Selector Parser Module
负责解析Newick树文件、分组表文件和PCA文件|Responsible for parsing Newick tree files, group table files, and PCA files
"""

import re
from typing import Dict, List, Tuple


class NewickParser:
    """Newick树文件解析器|Newick Tree File Parser"""

    def __init__(self, logger):
        """初始化解析器|Initialize parser

        Args:
            logger: 日志器|Logger object
        """
        self.logger = logger

    def parse_newick(self, newick_str: str) -> List[Dict]:
        """解析Newick格式字符串，提取所有样品及其枝长|Parse Newick format string and extract all samples with branch lengths

        保持Newick文件中的原始顺序（拓扑顺序）|Maintain original order from Newick file (topological order)

        Args:
            newick_str: Newick格式字符串|Newick format string

        Returns:
            List[Dict]: 样品列表，每个样品包含name和branch_length|List of samples, each containing name and branch_length
        """
        # 正则表达式匹配样品名称和枝长|Regex to match sample names and branch lengths
        # 匹配模式: 样品名:枝长|Match pattern: sample_name:branch_length
        pattern = r'([A-Za-z0-9_\-\.\[\]]+):([0-9.eE+-]+)'
        matches = re.findall(pattern, newick_str)

        samples = []
        seen = set()

        # 保持Newick字符串中的原始顺序（即拓扑顺序）
        # Maintain original order from Newick string (topological order)
        for name, length in matches:
            # 过滤掉内部节点（如内部节点编号）|Filter out internal nodes
            if name and name not in seen:
                try:
                    branch_length = float(length)
                    samples.append({
                        'name': name,
                        'branch_length': branch_length
                    })
                    seen.add(name)
                except ValueError:
                    self.logger.warning(f"无效的枝长值|Invalid branch length value: {length}")
                    continue

        self.logger.info(f"从Newick树中提取到|Extracted from Newick tree: {len(samples)} 个样品|samples")
        return samples

    def load_newick_file(self, file_path: str) -> str:
        """加载Newick文件|Load Newick file

        Args:
            file_path: Newick文件路径|Newick file path

        Returns:
            str: Newick格式字符串|Newick format string
        """
        with open(file_path, 'r', encoding='utf-8') as f:
            newick_str = f.read().strip()

        # 移除可能的换行符和多余空格|Remove possible newlines and extra spaces
        newick_str = ' '.join(newick_str.split())

        self.logger.info(f"成功加载Newick文件|Successfully loaded Newick file: {file_path}")
        return newick_str


class GroupTableParser:
    """分组表文件解析器|Group Table File Parser"""

    def __init__(self, logger):
        """初始化解析器|Initialize parser

        Args:
            logger: 日志器|Logger object
        """
        self.logger = logger

    def parse_group_file(self, file_path: str, delimiter: str = None) -> Dict[str, str]:
        """解析分组表文件|Parse group table file

        支持多种格式:
        - CSV格式 (逗号分隔)|CSV format (comma-separated)
        - TSV格式 (制表符分隔)|TSV format (tab-separated)
        - 空格分隔|Space-separated

        Args:
            file_path: 分组表文件路径|Group table file path
            delimiter: 分隔符，None表示自动检测|Delimiter, None means auto-detect

        Returns:
            Dict[str, str]: 样品名到分组的映射|Sample name to group mapping
        """
        # 检测文件格式|Detect file format
        if delimiter is None:
            delimiter = self._detect_delimiter(file_path)

        sample_to_group = {}

        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # 跳过空行和注释行|Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue

                parts = line.split(delimiter)

                if len(parts) < 2:
                    self.logger.warning(f"第{line_num}行格式错误，跳过|Line {line_num} format error, skipped")
                    continue

                sample_name = parts[0].strip()
                group_name = parts[1].strip()

                if sample_name and group_name:
                    sample_to_group[sample_name] = group_name

        self.logger.info(f"从分组文件中加载了|Loaded from group file: {len(sample_to_group)} 个样品的分组信息|sample group assignments")
        self.logger.info(f"检测到|Detected: {len(set(sample_to_group.values()))} 个分组|groups")

        return sample_to_group

    def _detect_delimiter(self, file_path: str) -> str:
        """自动检测文件分隔符|Auto-detect file delimiter

        Args:
            file_path: 文件路径|File path

        Returns:
            str: 检测到的分隔符|Detected delimiter
        """
        delimiters = ['\t', ',', ' ']

        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                for delimiter in delimiters:
                    if delimiter in line and line.count(delimiter) >= 1:
                        return delimiter

                break

        # 默认使用制表符|Default to tab
        return '\t'

    def get_group_statistics(self, sample_to_group: Dict[str, str]) -> Dict[str, int]:
        """获取分组统计信息|Get group statistics

        Args:
            sample_to_group: 样品到分组的映射|Sample to group mapping

        Returns:
            Dict[str, int]: 分组名到样品数的映射|Group name to sample count mapping
        """
        group_counts = {}

        for group in sample_to_group.values():
            group_counts[group] = group_counts.get(group, 0) + 1

        return group_counts


class PCAFileParser:
    """PCA文件解析器|PCA File Parser"""

    def __init__(self, logger):
        """初始化解析器|Initialize parser

        Args:
            logger: 日志器|Logger object
        """
        self.logger = logger

    def parse_pca_file(self, file_path: str, n_components: int = 10) -> Dict[str, List[float]]:
        """解析PCA文件|Parse PCA file

        Args:
            file_path: PCA文件路径|PCA file path
            n_components: 使用的PC组件数|Number of PC components to use

        Returns:
            Dict[str, List[float]]: 样品名到PC坐标的映射|Sample name to PC coordinates mapping
        """
        sample_to_pca = {}

        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # 跳过表头|Skip header
        if lines[0].startswith('FID') or lines[0].startswith('#'):
            data_lines = lines[1:]
        else:
            data_lines = lines

        for line_num, line in enumerate(data_lines, start=2):
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                # 尝试空格分隔|Try space delimiter
                parts = line.split()

            if len(parts) < 3:
                self.logger.warning(f"第{line_num}行格式错误，跳过|Line {line_num} format error, skipped")
                continue

            # IID列通常是第二列|IID column is usually the second column
            sample_name = parts[1].strip()

            # 提取PC值|Extract PC values
            try:
                pc_values = []
                for i in range(min(n_components, len(parts) - 2)):
                    pc_val = float(parts[2 + i])
                    pc_values.append(pc_val)

                sample_to_pca[sample_name] = pc_values
            except (ValueError, IndexError) as e:
                self.logger.warning(f"第{line_num}行PC值解析错误|Line {line_num} PC value parsing error: {e}")
                continue

        self.logger.info(f"从PCA文件中加载了|Loaded from PCA file: {len(sample_to_pca)} 个样品的PCA坐标|sample PCA coordinates")
        self.logger.info(f"使用|Using: {n_components} 个主成分|principal components")

        return sample_to_pca

    def calculate_pca_distance(self, pc1: List[float], pc2: List[float]) -> float:
        """计算两个PCA坐标之间的欧氏距离|Calculate Euclidean distance between two PCA coordinates

        Args:
            pc1: 第一个样品的PC坐标|PC coordinates of first sample
            pc2: 第二个样品的PC坐标|PC coordinates of second sample

        Returns:
            float: 欧氏距离|Euclidean distance
        """
        if len(pc1) != len(pc2):
            self.logger.warning(f"PC坐标维度不匹配|PC coordinate dimension mismatch: {len(pc1)} vs {len(pc2)}")
            min_len = min(len(pc1), len(pc2))
            pc1 = pc1[:min_len]
            pc2 = pc2[:min_len]

        # 计算欧氏距离|Calculate Euclidean distance
        distance = sum((a - b) ** 2 for a, b in zip(pc1, pc2)) ** 0.5

        return distance
