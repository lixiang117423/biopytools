"""
互作转录组结果处理模块|Dual RNA-seq Results Processing Module
"""

import os
import pandas as pd
from typing import List


class DualMatrixMerger:
    """双物种表达矩阵合并器|Dual-species Expression Matrix Merger"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def merge_expression_matrix(self, fpkm_files: List[str], species_name: str) -> bool:
        """合并表达矩阵|Merge expression matrix"""
        output_dir = self.config.output_dir
        matrix_dir = os.path.join(output_dir, "05.expression_matrix")
        os.makedirs(matrix_dir, exist_ok=True)

        self.logger.info(f"合并{species_name}表达矩阵|Merging {species_name} expression matrix")

        all_data = []

        # 读取FPKM文件|Read all FPKM files
        for fpkm_file in fpkm_files:
            if os.path.exists(fpkm_file):
                self.logger.info(f"读取文件|Reading file: {fpkm_file}")
                try:
                    # gene_id, transcript_id, cov, FPKM, TPM, sample
                    # 读取文件，假设六列|Read file with six columns
                    df = pd.read_csv(fpkm_file, sep="\t")

                    # 检查数据类型并转换|Check data types and convert
                    df["cov"] = pd.to_numeric(df["cov"], errors="coerce")
                    df["FPKM"] = pd.to_numeric(df["FPKM"], errors="coerce")
                    df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")

                    # 移除NaN值|Remove NaN values
                    df = df.dropna()

                    all_data.append(df)
                    self.logger.info(f"成功读取|Successfully read {len(df)}行|rows of data")

                except Exception as e:
                    self.logger.error(f"读取文件错误|Error reading file {fpkm_file}: {e}")
                    continue
            else:
                self.logger.warning(f"文件不存在|File does not exist: {fpkm_file}")

        if not all_data:
            self.logger.error(f"未找到有效的FPKM数据|No valid FPKM data found for {species_name}")
            return False

        self.logger.info(f"找到|Found {len(all_data)}个有效数据文件|valid data files")

        # 合并所有数据|Merge all data
        combined_df = pd.concat(all_data, ignore_index=True)
        self.logger.info(f"合并后总计|Total {len(combined_df)}行|rows after merging")

        # 获取样本名列表|Get sample name list
        sample_names = sorted([str(x) for x in combined_df["sample"].unique()])
        self.logger.info(f"样本名称|Sample names: {', '.join(sample_names)}")

        # 生成物种表达矩阵|Generate species expression matrix
        output_file = os.path.join(matrix_dir, f"{species_name}_matrix.txt")
        self._generate_matrix(combined_df, output_file, species_name)

        return True

    def _generate_matrix(self, df: pd.DataFrame, output_file: str, species_name: str):
        """生成表达矩阵文件|Generate expression matrix file"""
        self.logger.info(f"生成{species_name}表达矩阵|Generating {species_name} expression matrix")

        # 重新排列列|Reorder columns
        combined_df_ordered = df[["gene_id", "transcript_id", "cov", "FPKM", "TPM", "sample"]]
        combined_df_ordered.to_csv(output_file, sep="\t", index=False, header=True)

        self.logger.info(f"保存完成|Save completed: {output_file} ({len(combined_df_ordered)}行|rows)")


class SummaryGenerator:
    """总结报告生成器|Summary Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary_report(self) -> bool:
        """生成总结报告|Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "dual_rnaseq_summary.txt")

        try:
            self.logger.info("生成总结报告|Generating summary report")

            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("=" * 60 + "\n")
                f.write("互作转录组分析总结报告|Dual RNA-seq Analysis Summary Report\n")
                f.write("=" * 60 + "\n\n")

                # 输入文件信息|Input file information
                f.write("输入文件|Input Files:\n")
                f.write(f"  物种1|Species 1 ({self.config.species1_name}):\n")
                f.write(f"    - 基因组文件|Genome file: {self.config.species1_genome}\n")
                f.write(f"    - GTF文件|GTF file: {self.config.species1_gtf}\n")
                f.write(f"  物种2|Species 2 ({self.config.species2_name}):\n")
                f.write(f"    - 基因组文件|Genome file: {self.config.species2_genome}\n")
                f.write(f"    - GTF文件|GTF file: {self.config.species2_gtf}\n")
                f.write(f"  输入路径|Input path: {self.config.input_path}\n\n")

                # 配置参数|Configuration parameters
                f.write("配置参数|Configuration Parameters:\n")
                f.write(f"  - 线程数|Thread count: {self.config.threads}\n")
                f.write(f"  - 输出目录|Output directory: {self.config.output_dir}\n")
                f.write(f"  - 最小MAPQ|Minimum MAPQ: {self.config.min_mapq}\n")
                f.write(f"  - 仅保留唯一比对|Unique only: {self.config.unique_only}\n")

                if self.config.fastq_pattern:
                    f.write(f"  - FASTQ模式|FASTQ pattern: {self.config.fastq_pattern}\n")

                # 样本信息|Sample information
                if hasattr(self.config, 'samples') and self.config.samples:
                    f.write(f"\n样本信息|Sample Information:\n")
                    f.write(f"  - 样本数量|Sample count: {len(self.config.samples)}\n")
                    f.write(f"  - 样本列表|Sample list:\n")
                    for sample in self.config.samples:
                        f.write(f"    * {sample['name']}\n")

                f.write(f"\n输出文件|Output Files:\n")
                f.write(f"  - 01.index/: HISAT2索引文件|HISAT2 index files\n")
                f.write(f"  - 02.classification/: Reads物种分类结果|Reads species classification results\n")
                f.write(f"  - 03.alignment_statistics/: 比对统计报告|Alignment statistics report\n")
                f.write(f"    - mapping_statistics.tsv: 比对统计汇总表|Alignment statistics summary\n")
                f.write(f"  - 04.quantification/: 定量结果|Quantification results\n")
                f.write(f"  - 05.expression_matrix/: 双物种表达矩阵|Dual-species expression matrix\n")
                f.write(f"    - {self.config.species1_name}_matrix.txt: 物种1表达矩阵|Species 1 expression matrix\n")
                f.write(f"    - {self.config.species2_name}_matrix.txt: 物种2表达矩阵|Species 2 expression matrix\n")

                f.write("=" * 60 + "\n")

            self.logger.info(f"总结报告已生成|Summary report generated: {report_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成总结报告时出错|Error generating summary report: {e}")
            return False
