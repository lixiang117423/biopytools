"""
基因组挂载率统计核心计算模块|Genome Mount Rate Core Calculation Module
"""

import os
from typing import List, Dict, Any


class GenomeMountRateCalculator:
    """基因组挂载率计算器|Genome Mount Rate Calculator"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: GenomeMountRateConfig配置对象|GenomeMountRateConfig object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger
        self.lengths = []

    def read_fasta_lengths(self) -> List[int]:
        """读取FASTA文件并返回所有序列长度的列表|Read FASTA file and return list of sequence lengths

        Returns:
            List[int]: 序列长度列表|List of sequence lengths
        """
        lengths = []
        current_len = 0
        in_seq = False

        fasta_file = self.config.fasta_file

        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue  # 跳过空行|Skip empty lines

                    if line.startswith(">"):
                        if in_seq:
                            lengths.append(current_len)
                        current_len = 0
                        in_seq = True
                    else:
                        current_len += len(line)

                # 添加最后一条序列|Add the last sequence
                if in_seq:
                    lengths.append(current_len)

            self.logger.info(f"成功读取|Successfully read {len(lengths)} 条序列|sequences")

        except FileNotFoundError:
            self.logger.error(f"找不到文件|File not found: {fasta_file}")
            raise
        except Exception as e:
            self.logger.error(f"读取文件时出错|Error reading file: {str(e)}")
            raise

        self.lengths = lengths
        return lengths

    def calculate(self) -> Dict[str, Any]:
        """计算挂载率|Calculate mount rate

        Returns:
            Dict[str, Any]: 计算结果字典|Calculation results dictionary
            包含|Contains:
                - total_seqs: 总序列数|Total number of sequences
                - total_bp: 总基因组大小|Total genome size in bp
                - target_n: 目标序列数|Target number of sequences
                - target_bp: 目标序列总长|Total length of target sequences
                - percentage: 占比百分比|Percentage
                - sorted: 是否排序|Whether sorted by length
        """
        # 1. 读取所有长度|Read all lengths
        self.logger.info(f"正在读取文件|Reading file: {self.config.fasta_file}")
        lengths = self.read_fasta_lengths()

        total_seqs = len(lengths)
        total_bp = sum(lengths)

        if total_seqs == 0:
            self.logger.error("文件中没有读取到序列|No sequences found in file")
            raise ValueError("文件中没有序列|No sequences in file")

        # 2. 处理排序逻辑|Handle sorting logic
        if self.config.sort_by_length:
            self.logger.info("正在按序列长度从大到小排序|Sorting sequences by length (descending)")
            lengths.sort(reverse=True)
            sorted_flag = True
        else:
            self.logger.info("保持文件原有顺序|Keeping original file order")
            sorted_flag = False

        # 3. 获取前 N 条|Get top N sequences
        target_n = min(self.config.number, total_seqs)
        top_n_lengths = lengths[:target_n]
        top_n_bp = sum(top_n_lengths)

        # 4. 计算百分比|Calculate percentage
        percentage = (top_n_bp / total_bp) * 100

        # 5. 返回结果|Return results
        results = {
            'total_seqs': total_seqs,
            'total_bp': total_bp,
            'target_n': target_n,
            'target_bp': top_n_bp,
            'percentage': percentage,
            'sorted': sorted_flag
        }

        return results

    def format_number(self, num: int) -> str:
        """格式化数字，添加千位分隔符|Format number with thousands separator

        Args:
            num: 要格式化的数字|Number to format

        Returns:
            str: 格式化后的字符串|Formatted string
        """
        return f"{num:,}"
