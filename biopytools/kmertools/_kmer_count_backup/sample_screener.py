"""样本筛选器|Sample Screener

基于连续0数量筛选样本
Screen samples based on consecutive zeros
"""

import csv
import os
from pathlib import Path
from typing import Dict, Set, List


class SampleScreener:
    """样本筛选器|Sample Screener

    筛选符合条件的样本（全1或连续0不超过阈值）
    Screen samples that meet criteria (all 1s or consecutive 0s not exceeding threshold)
    """

    def __init__(self, config, logger):
        """初始化样本筛选器|Initialize sample screener

        Args:
            config: KmerPAVConfig配置对象|KmerPAVConfig object
            logger: 日志器|Logger instance
        """
        self.config = config
        self.logger = logger

    def process_matrix_file(self, matrix_file_path: str, gene_name: str) -> List[str]:
        """处理单个矩阵文件|Process single matrix file

        Args:
            matrix_file_path: 矩阵文件路径|Matrix file path
            gene_name: 基因名称|Gene name

        Returns:
            List[str]: 符合条件的样本名列表|List of sample names that meet criteria
        """
        self.logger.info(f"处理基因|Processing gene: {gene_name}")

        sample_names = []

        try:
            with open(matrix_file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)  # 跳过表头|Skip header

                for row in reader:
                    if not row:
                        continue

                    sample_name = row[0]
                    matrix_data = row[1:]

                    all_ones = True
                    zero_segments_valid = True
                    zero_segment_length = 0
                    has_zeros = False

                    for data_point in matrix_data:
                        if data_point == '0':
                            all_ones = False
                            has_zeros = True
                            zero_segment_length += 1
                        elif data_point == '1':
                            if zero_segment_length > self.config.max_zeros:
                                zero_segments_valid = False
                                break
                            zero_segment_length = 0
                        else:
                            continue  # 忽略其他值|Ignore other values

                    # 检查行末尾的连续0段|Check consecutive 0s at end of row
                    if zero_segment_length > self.config.max_zeros:
                        zero_segments_valid = False

                    if all_ones:
                        sample_names.append(sample_name)
                        self.logger.debug(f"  基因|Gene: {gene_name}, 样本|Sample: {sample_name} - 符合条件|Meets criteria (全1|all 1s)")
                    elif has_zeros and zero_segments_valid:
                        sample_names.append(sample_name)
                        self.logger.debug(f"  基因|Gene: {gene_name}, 样本|Sample: {sample_name} - 符合条件|Meets criteria (连续0不超过|consecutive 0s not exceeding {self.config.max_zeros})")
                    else:
                        self.logger.debug(f"  基因|Gene: {gene_name}, 样本|Sample: {sample_name} - 不符合|Does not meet criteria")

        except FileNotFoundError:
            self.logger.error(f"矩阵文件未找到|Matrix file not found: {matrix_file_path}")

        return sample_names

    def screen_all_genes(self) -> Dict[str, List[str]]:
        """筛选所有基因|Screen all genes

        Returns:
            Dict[str, List[str]]: 基因到样本列表的映射|Gene to sample list mapping
        """
        self.logger.info("开始筛选所有基因|Starting to screen all genes")

        # 读取基因列表|Read gene list
        gene_list = []
        try:
            with open(self.config.gene_list, 'r') as f:
                for line in f:
                    gene_list.append(line.strip())
        except FileNotFoundError:
            self.logger.error(f"基因列表文件未找到|Gene list file not found: {self.config.gene_list}")
            return {}

        self.logger.info(f"待筛选基因数量|Genes to screen: {len(gene_list)}")

        # 处理每个基因|Process each gene
        output_data = {}
        for gene_name in gene_list:
            gene_dir = Path(self.config.matrix_dir) / gene_name
            matrix_filename = f"{gene_name}{self.config.matrix_suffix}"
            matrix_filepath = gene_dir / matrix_filename

            if gene_dir.exists() and matrix_filepath.exists():
                sample_names = self.process_matrix_file(str(matrix_filepath), gene_name)
                output_data[gene_name] = sample_names
            else:
                self.logger.warning(f"基因|Gene: {gene_name} 的目录或矩阵文件不存在|directory or matrix file not found")
                output_data[gene_name] = []

        # 写入筛选结果|Write screening results
        output_csv = self.output_path / "screening_results.csv"
        self._write_screening_results(output_data, output_csv)

        return output_data

    def _write_screening_results(self, output_data: Dict[str, List[str]], output_file: Path) -> bool:
        """写入筛选结果|Write screening results

        Args:
            output_data: 筛选结果数据|Screening result data
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"写入筛选结果|Writing screening results to: {output_file}")

        try:
            with open(output_file, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Gene', 'Sample Names', 'Sample Count'])

                for gene_name, sample_names in output_data.items():
                    sample_names_str = ";".join(sample_names)
                    sample_count = len(sample_names)
                    csv_writer.writerow([gene_name, sample_names_str, sample_count])

            self.logger.info("筛选结果已写入|Screening results written")
            return True

        except Exception as e:
            self.logger.error(f"写入筛选结果失败|Failed to write screening results: {e}")
            return False

    @property
    def output_path(self):
        """获取输出路径|Get output path"""
        return Path(self.config.output_dir)

    def generate_gene_sample_matrix(self) -> bool:
        """生成样本-基因矩阵|Generate sample-gene matrix

        从筛选结果生成样本×基因的矩阵
        Generate sample×gene matrix from screening results

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("生成样本-基因矩阵|Generating sample-gene matrix")

        output_csv = self.output_path / "screening_results.csv"
        matrix_output = self.output_path / "gene_sample_matrix.csv"

        if not output_csv.exists():
            self.logger.error(f"筛选结果文件不存在|Screening results file not found: {output_csv}")
            return False

        # 读取筛选结果|Read screening results
        gene_sample_dict = {}
        all_genes = []
        all_samples = []

        try:
            with open(output_csv, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    gene_name = row['Gene']
                    sample_names_str = row['Sample Names']
                    sample_names = []
                    if sample_names_str:
                        sample_names = sample_names_str.split(';')

                    gene_sample_dict[gene_name] = set(sample_names)
                    all_genes.append(gene_name)

            self.logger.info(f"读取基因数量|Genes read: {len(all_genes)}")

        except Exception as e:
            self.logger.error(f"读取筛选结果失败|Failed to read screening results: {e}")
            return False

        # 读取样本列表|Read sample list
        try:
            with open(self.config.sample_list, 'r') as f:
                for line in f:
                    sample_name = line.strip()
                    all_samples.append(sample_name)

            self.logger.info(f"读取样本数量|Samples read: {len(all_samples)}")

        except FileNotFoundError:
            self.logger.error(f"样本列表文件未找到|Sample list file not found: {self.config.sample_list}")
            return False

        # 生成矩阵|Generate matrix
        try:
            matrix_data = []
            header_row = ['Sample'] + all_genes
            matrix_data.append(header_row)

            for sample_name in all_samples:
                sample_row = [sample_name]
                for gene_name in all_genes:
                    sample_set = gene_sample_dict.get(gene_name, set())
                    if sample_name in sample_set:
                        sample_row.append('1')
                    else:
                        sample_row.append('0')
                matrix_data.append(sample_row)

            # 写入矩阵文件|Write matrix file
            with open(matrix_output, 'w', newline='') as csv_matrix_file:
                writer = csv.writer(csv_matrix_file)
                writer.writerows(matrix_data)

            self.logger.info(f"样本-基因矩阵已生成|Sample-gene matrix generated: {matrix_output}")
            return True

        except Exception as e:
            self.logger.error(f"生成矩阵失败|Failed to generate matrix: {e}")
            return False
