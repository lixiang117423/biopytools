# """
# RNA-seq | RNA-seq Results Processing Module
# """

# import os
# import pandas as pd
# from pathlib import Path
# from typing import List

# class ExpressionMatrixMerger:
#     """ | Expression Matrix Merger"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def merge_expression_matrix(self, fpkm_files: List[str]) -> bool:
#         """ | Merge expression matrix from all samples"""
#         output_dir = self.config.output_dir
        
#         self.logger.info(" | Merging expression matrix")
        
#         all_data = []

#         # FPKM | Read all FPKM files
#         for fpkm_file in fpkm_files:
#             if os.path.exists(fpkm_file):
#                 self.logger.info(f" | Reading file: {fpkm_file}")
#                 try:
#                     # gene_id, transcript_id, cov, FPKM, TPM, sample
#                     # Read file, assuming six columns: gene_id, transcript_id, cov, FPKM, TPM, sample
#                     df = pd.read_csv(fpkm_file, sep="\t")

#                     #  | Check data types and convert
#                     df["cov"] = pd.to_numeric(df["cov"], errors="coerce")
#                     df["FPKM"] = pd.to_numeric(df["FPKM"], errors="coerce")
#                     df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")

#                     # NaN | Remove NaN values
#                     df = df.dropna()

#                     all_data.append(df)
#                     self.logger.info(f"   {len(df)}  | Successfully read {len(df)} rows of data")

#                 except Exception as e:
#                     self.logger.error(f"   | Error reading file {fpkm_file}: {e}")
#                     continue
#             else:
#                 self.logger.warning(f" | File does not exist: {fpkm_file}")

#         if not all_data:
#             self.logger.error(" FPKM | No valid FPKM data found")
#             return False

#         self.logger.info(f" {len(all_data)}  | Found {len(all_data)} valid data files")

#         #  | Merge all data
#         combined_df = pd.concat(all_data, ignore_index=True)
#         self.logger.info(f" {len(combined_df)}  | Total {len(combined_df)} rows after merging")

#         #  | Get sample name list
#         sample_names = sorted(combined_df["sample"].unique())
#         self.logger.info(f" | Sample names: {', '.join(sample_names)}")

#         #  | Ensure output directory exists
#         os.makedirs(output_dir, exist_ok=True)

#         #  all.fpkm.tpm.txt () | Generate all.fpkm.tpm.txt (direct merge of all data)
#         self.logger.info(" all.fpkm.tpm.txt | Generating all.fpkm.tpm.txt...")
#         all_output_file = os.path.join(output_dir, "all.fpkm.tpm.txt")
#         # gene_id, transcript_id, cov, FPKM, TPM, sample | Reorder columns
#         combined_df_ordered = combined_df[["gene_id", "transcript_id", "cov", "FPKM", "TPM", "sample"]]
#         combined_df_ordered.to_csv(all_output_file, sep="\t", index=False, header=True)
#         self.logger.info(f"  | Save completed: {all_output_file} ({len(combined_df_ordered)}  | rows)")

#         return True

# class SummaryGenerator:
#     """ | Summary Generator"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_summary_report(self) -> bool:
#         """ | Generate summary report"""
#         report_file = os.path.join(self.config.output_dir, "rnaseq_summary.txt")
        
#         try:
#             with open(report_file, 'w', encoding='utf-8') as f:
#                 f.write("RNA-seq | RNA-seq Analysis Summary Report\n")
#                 f.write("=" * 50 + "\n\n")
                
#                 #  | Input file information
#                 f.write(" | Input Files:\n")
#                 f.write(f"  -  | Genome file: {self.config.genome_file}\n")
#                 f.write(f"  - GTF | GTF file: {self.config.gtf_file}\n")
#                 f.write(f"  -  | Input path: {self.config.input_path}\n\n")
                
#                 #  | Configuration parameters
#                 f.write(" | Configuration Parameters:\n")
#                 f.write(f"  -  | Thread count: {self.config.threads}\n")
#                 f.write(f"  -  | Output directory: {self.config.output_dir}\n")
#                 f.write(f"  - BAM | Remove BAM files: {' | Yes' if self.config.remove_bam.lower() in ['yes', 'y'] else ' | No'}\n")
                
#                 if self.config.fastq_pattern:
#                     f.write(f"  - FASTQ | FASTQ pattern: {self.config.fastq_pattern}\n")
                
#                 #  | Sample information
#                 if hasattr(self.config, 'samples') and self.config.samples:
#                     f.write(f"\n | Sample Information:\n")
#                     f.write(f"  -  | Sample count: {len(self.config.samples)}\n")
#                     f.write(f"  -  | Sample list:\n")
#                     for sample in self.config.samples:
#                         f.write(f"    * {sample['name']}\n")
                
#                 f.write(f"\n | Output Files:\n")
#                 f.write(f"  - all.fpkm.tpm.txt: FPKMTPM | FPKM and TPM data for all samples\n")
#                 f.write(f"  - stringtie_output/: StringTie | StringTie quantification results\n")
#                 f.write(f"  - fpkm_output/: FPKM | FPKM files for each sample\n")
                
#                 if self.config.remove_bam.lower() not in ['yes', 'y']:
#                     f.write(f"  - *.sorted.bam: BAM | Alignment result BAM files\n")
            
#             self.logger.info(f" | Summary report generated: {report_file}")
#             return True
            
#         except Exception as e:
#             self.logger.error(f" | Error generating summary report: {e}")
#             return False

"""
RNA-seq | RNA-seq Results Processing Module
"""

import os
import pandas as pd
from pathlib import Path
from typing import List

class ExpressionMatrixMerger:
    """ | Expression Matrix Merger"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def merge_expression_matrix(self, fpkm_files: List[str]) -> bool:
        """合并表达矩阵|Merge expression matrix from all samples"""
        output_dir = self.config.output_dir

        self.logger.info("合并表达矩阵|Merging expression matrix")

        all_data = []

        # 读取FPKM文件|Read all FPKM files
        for fpkm_file in fpkm_files:
            if os.path.exists(fpkm_file):
                self.logger.info(f"读取文件|Reading file: {fpkm_file}")
                try:
                    # gene_id, transcript_id, cov, FPKM, TPM, sample
                    # 读取文件，假设六列：gene_id, transcript_id, cov, FPKM, TPM, sample
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
            self.logger.error("未找到有效的FPKM数据|No valid FPKM data found")
            return False

        self.logger.info(f"找到|Found {len(all_data)}个有效数据文件|valid data files")

        # 合并所有数据|Merge all data
        combined_df = pd.concat(all_data, ignore_index=True)
        self.logger.info(f"合并后总计|Total {len(combined_df)}行|rows after merging")

        # 获取样本名列表|Get sample name list
        sample_names = sorted([str(x) for x in combined_df["sample"].unique()])
        self.logger.info(f"样本名称|Sample names: {', '.join(sample_names)}")

        # 确保输出目录存在|Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # 生成all.fpkm.tpm.txt（直接合并所有数据）|Generate all.fpkm.tpm.txt (direct merge of all data)
        self.logger.info("生成all.fpkm.tpm.txt文件|Generating all.fpkm.tpm.txt...")
        all_output_file = os.path.join(output_dir, "all.fpkm.tpm.txt")
        # 重新排列列：gene_id, transcript_id, cov, FPKM, TPM, sample|Reorder columns
        combined_df_ordered = combined_df[["gene_id", "transcript_id", "cov", "FPKM", "TPM", "sample"]]
        combined_df_ordered.to_csv(all_output_file, sep="\t", index=False, header=True)
        self.logger.info(f"保存完成|Save completed: {all_output_file} ({len(combined_df_ordered)}行|rows)")

        return True

class SummaryGenerator:
    """总结报告生成器|Summary Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary_report(self) -> bool:
        """生成总结报告|Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "rnaseq_summary.txt")

        try:
            self.logger.info("生成总结报告|Generating summary report")

            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("=" * 60 + "\n")
                f.write("RNA-seq分析总结报告|RNA-seq Analysis Summary Report\n")
                f.write("=" * 60 + "\n\n")

                # 输入文件信息|Input file information
                f.write("输入文件|Input Files:\n")
                f.write(f"  - 基因组文件|Genome file: {self.config.genome_file}\n")
                f.write(f"  - GTF文件|GTF file: {self.config.gtf_file}\n")
                f.write(f"  - 输入路径|Input path: {self.config.input_path}\n\n")

                # 配置参数|Configuration parameters
                f.write("配置参数|Configuration Parameters:\n")
                f.write(f"  - 线程数|Thread count: {self.config.threads}\n")
                f.write(f"  - 输出目录|Output directory: {self.config.output_dir}\n")
                remove_bam = "是|Yes" if self.config.remove_bam.lower() in ['yes', 'y'] else "否|No"
                f.write(f"  - 删除BAM文件|Remove BAM files: {remove_bam}\n")

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
                f.write(f"  - all.fpkm.tpm.txt: 所有样本的FPKM和TPM数据|FPKM and TPM data for all samples\n")
                f.write(f"  - stringtie_output/: StringTie定量结果|StringTie quantification results\n")
                f.write(f"  - fpkm_output/: 各样本的FPKM文件|FPKM files for each sample\n")

                if self.config.remove_bam.lower() not in ['yes', 'y']:
                    f.write(f"  - *.sorted.bam: 比对结果BAM文件|Alignment result BAM files\n")

                f.write("=" * 60 + "\n")

            self.logger.info(f"总结报告已生成|Summary report generated: {report_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成总结报告时出错|Error generating summary report: {e}")
            return False