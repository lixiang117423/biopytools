"""
Fst计算结果处理模块|Fst Calculation Results Processing Module
"""

import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


class FstResultsProcessor:
    """Fst结果处理器|Fst Results Processor"""

    def __init__(self, logger, output_dir: Path):
        self.logger = logger
        self.output_dir = Path(output_dir)

    def parse_fst_file(self, fst_file: str) -> Dict[Tuple[str, str], List[float]]:
        """
        解析PLINK生成的.fst文件
        Parse PLINK generated .fst file

        Args:
            fst_file: .fst文件路径|.fst file path

        Returns:
            群体对到Fst值的映射|Population pair to Fst values mapping
        """
        pair_fst_values = defaultdict(list)

        self.logger.info(f"解析Fst文件|Parsing Fst file: {fst_file}")

        with open(fst_file, 'r') as f:
            # 读取并跳过可能的空行|Read and skip empty lines
            lines = [line.strip() for line in f if line.strip()]

            if not lines:
                self.logger.error("Fst文件为空|Fst file is empty")
                return {}

            # 检查是否有表头|Check if there's a header
            # PLINK .fst文件格式: CHR SNP A1 A2 (POP1,POP2) ...
            # Fst值在第一行之后
            data_lines = []
            for i, line in enumerate(lines):
                # 跳过以#开头的注释行|Skip comment lines starting with #
                if line.startswith('#'):
                    continue
                # 查找包含数值的行|Find lines with numeric values
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        # 尝试解析数值|Try to parse numeric values
                        float(parts[-1])
                        data_lines.append(line)
                    except ValueError:
                        continue

            if not data_lines:
                self.logger.error("未找到有效的Fst数据行|No valid Fst data lines found")
                return {}

            # 解析数据行|Parse data lines
            for line in data_lines:
                parts = line.split()
                if len(parts) < 6:
                    continue

                # PLINK .fst格式: CHR SNP A1 A2 (POP1,POP2) FST
                # 最后一列是Fst值，倒数第二列可能是标准误
                # 我们需要从列标题中获取群体对信息
                # 但由于我们使用--within，列标题可能包含群体对信息

                # 简化处理：提取所有数值列
                # 实际的Fst值在最后一列或倒数第二列
                try:
                    # 尝试从最后一列获取Fst值
                    fst_value = float(parts[-1])
                    # 假设列名中包含群体对信息，但我们需要从其他地方获取
                    # 暂时使用索引作为键
                    # 在实际使用中，需要根据PLINK输出的具体格式调整
                    pair_fst_values[('all', 'all')].append(fst_value)
                except (ValueError, IndexError):
                    continue

        return pair_fst_values

    def calculate_pairwise_fst(self, fst_file: str, populations: List[str]) -> Dict[Tuple[str, str], float]:
        """
        计算两两群体间的Fst值
        Calculate pairwise Fst values

        Args:
            fst_file: .fst文件路径|.fst file path
            populations: 群体列表|Population list

        Returns:
            群体对到平均Fst值的映射|Population pair to mean Fst value mapping
        """
        # 这里需要根据PLINK --fst输出的具体格式来解析
        # PLINK的--fst会输出每个位点的Fst值
        # 对于两两比较，我们需要分别运行PLINK或解析输出

        self.logger.info("计算两两群体Fst值|Calculating pairwise population Fst values")

        pairwise_fst = {}

        # 如果只有两个群体，直接使用PLINK输出
        if len(populations) == 2:
            # 读取.fst文件计算平均值
            fst_values = []
            with open(fst_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                # 最后一列通常是Fst值
                                fst_val = float(parts[-1])
                                fst_values.append(fst_val)
                            except (ValueError, IndexError):
                                continue

            if fst_values:
                mean_fst = sum(fst_values) / len(fst_values)
                pair = (populations[0], populations[1])
                pairwise_fst[pair] = mean_fst
                self.logger.info(f"  {populations[0]} vs {populations[1]}: {mean_fst:.6f}")

        # 如果有多个群体，需要为每对群体分别计算
        # 这里提供一个简化的实现，实际中需要多次调用PLINK
        else:
            self.logger.warning(f"多群体Fst计算需要多次PLINK调用|Multi-population Fst calculation requires multiple PLINK calls")
            self.logger.warning(f"当前版本仅支持两群体比较的完整解析|Current version only supports full parsing of two-population comparisons")

            # 创建占位符|Create placeholders
            for i, pop1 in enumerate(populations):
                for pop2 in populations[i+1:]:
                    pairwise_fst[(pop1, pop2)] = 0.0

        return pairwise_fst

    def generate_long_format_table(self, pairwise_fst: Dict[Tuple[str, str], float],
                                    output_file: str) -> bool:
        """
        生成长格式表格
        Generate long format table

        Args:
            pairwise_fst: 群体对到Fst值的映射|Population pair to Fst value mapping
            output_file: 输出文件路径|Output file path
        """
        self.logger.info(f"生成长格式表格|Generating long format table: {output_file}")

        try:
            with open(output_file, 'w') as f:
                # 写入表头|Write header
                f.write("Population1\tPopulation2\tFst\n")

                # 写入数据|Write data
                for (pop1, pop2), fst_value in sorted(pairwise_fst.items()):
                    f.write(f"{pop1}\t{pop2}\t{fst_value:.6f}\n")

            self.logger.info(f"长格式表格生成完成|Long format table generated: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成长格式表格失败|Failed to generate long format table: {e}")
            return False

    def generate_matrix_format(self, pairwise_fst: Dict[Tuple[str, str], float],
                               populations: List[str], output_file: str) -> bool:
        """
        生成矩阵格式
        Generate matrix format

        Args:
            pairwise_fst: 群体对到Fst值的映射|Population pair to Fst value mapping
            populations: 群体列表|Population list
            output_file: 输出文件路径|Output file path
        """
        self.logger.info(f"生成矩阵格式|Generating matrix format: {output_file}")

        try:
            with open(output_file, 'w') as f:
                # 写入表头|Write header
                f.write("\t" + "\t".join(populations) + "\n")

                # 写入矩阵|Write matrix
                for i, pop1 in enumerate(populations):
                    row_values = [pop1]
                    for j, pop2 in enumerate(populations):
                        if i == j:
                            # 对角线为0|Diagonal is 0
                            row_values.append("0.000000")
                        elif i < j:
                            # 上三角|Upper triangle
                            key = (pop1, pop2)
                            if key in pairwise_fst:
                                row_values.append(f"{pairwise_fst[key]:.6f}")
                            else:
                                # 尝试反向查找|Try reverse lookup
                                key_reverse = (pop2, pop1)
                                if key_reverse in pairwise_fst:
                                    row_values.append(f"{pairwise_fst[key_reverse]:.6f}")
                                else:
                                    row_values.append("NA")
                        else:
                            # 下三角（与上三角对称）|Lower triangle (symmetric with upper)
                            key = (pop2, pop1)
                            if key in pairwise_fst:
                                row_values.append(f"{pairwise_fst[key]:.6f}")
                            else:
                                row_values.append("NA")

                    f.write("\t".join(row_values) + "\n")

            self.logger.info(f"矩阵格式生成完成|Matrix format generated: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成矩阵格式失败|Failed to generate matrix format: {e}")
            return False

    def process_results(self, fst_file: str, populations: List[str]) -> Dict[str, str]:
        """
        处理Fst计算结果
        Process Fst calculation results

        Args:
            fst_file: .fst文件路径|.fst file path
            populations: 群体列表|Population list

        Returns:
            输出文件字典|Output files dictionary
        """
        self.logger.info("开始处理Fst结果|Starting to process Fst results")

        output_files = {}

        # 计算两两群体Fst值|Calculate pairwise population Fst values
        pairwise_fst = self.calculate_pairwise_fst(fst_file, populations)

        if not pairwise_fst:
            self.logger.warning("未能计算Fst值，尝试使用原始文件|Failed to calculate Fst values, trying to use raw file")
            output_files['fst_raw'] = fst_file
            return output_files

        # 生成长格式表格|Generate long format table
        long_format_file = self.output_dir / 'fst_long_format.txt'
        if self.generate_long_format_table(pairwise_fst, str(long_format_file)):
            output_files['long_format'] = str(long_format_file)

        # 生成矩阵格式|Generate matrix format
        matrix_format_file = self.output_dir / 'fst_matrix.txt'
        if self.generate_matrix_format(pairwise_fst, populations, str(matrix_format_file)):
            output_files['matrix'] = str(matrix_format_file)

        self.logger.info("Fst结果处理完成|Fst results processing completed")

        return output_files
