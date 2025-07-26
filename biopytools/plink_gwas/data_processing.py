# """
# PLINK GWAS数据处理模块 | PLINK GWAS Data Processing Module
# """

# import pandas as pd
# import numpy as np
# from pathlib import Path
# from .utils import CommandRunner, FileManager

# class PhenotypeProcessor:
#     """表型数据处理器 | Phenotype Data Processor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def convert_phenotype(self) -> pd.DataFrame:
#         """转换表型文件格式 | Convert phenotype file format"""
#         self.logger.info("转换表型文件格式 | Converting phenotype file format...")
#         self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
        
#         # 读取原始表型文件 | Read original phenotype file
#         df = pd.read_csv("phenotype.txt", sep=r"\s+")
#         self.logger.info(f"原始表型文件形状 | Original phenotype file shape: {df.shape}")
        
#         # 检查列名 | Check column names
#         if df.shape[1] < 2:
#             raise ValueError("表型文件至少需要2列（样本ID和表型值） | Phenotype file needs at least 2 columns (sample ID and phenotype)")
        
#         # 假设第一列是样本ID，第二列是表型值 | Assume first column is sample ID, second is phenotype
#         sample_col = df.columns[0]
#         phenotype_col = df.columns[1]
        
#         # 统计原始表型分布 | Count original phenotype distribution
#         pheno_counts = df[phenotype_col].value_counts()
#         self.logger.info(f"原始表型分布 | Original phenotype distribution: {dict(pheno_counts)}")
        
#         if self.config.trait_type == "qualitative":
#             # 质量性状处理 | Qualitative trait processing
#             self.logger.info("处理质量性状：将0转换为1（对照），1转换为2（病例） | Processing qualitative trait: converting 0 to 1 (control), 1 to 2 (case)")
            
#             def convert_qualitative_value(val):
#                 """转换质量性状值 | Convert qualitative trait values"""
#                 if pd.isna(val):
#                     return -9
#                 elif val == 0:
#                     return 1  # 抗病 -> 对照 | Resistant -> Control
#                 elif val == 1:
#                     return 2  # 感病 -> 病例 | Susceptible -> Case
#                 else:
#                     return -9  # 其他值设为缺失 | Other values as missing
            
#             converted_phenotype = df[phenotype_col].apply(convert_qualitative_value)
            
#             # 统计转换后的分布 | Count converted distribution
#             converted_counts = converted_phenotype.value_counts()
#             self.logger.info(f"转换后表型分布 | Converted phenotype distribution: {dict(converted_counts)}")
#             self.logger.info(f"对照数 (抗病, 原值0) | Controls (resistant, original 0): {(converted_phenotype == 1).sum()}")
#             self.logger.info(f"病例数 (感病, 原值1) | Cases (susceptible, original 1): {(converted_phenotype == 2).sum()}")
            
#         else:
#             # 数量性状处理 | Quantitative trait processing
#             self.logger.info("处理数量性状：保持原始数值，只处理缺失值 | Processing quantitative trait: keeping original values, only handling missing values")
            
#             def convert_quantitative_value(val):
#                 """转换数量性状值 | Convert quantitative trait values"""
#                 if pd.isna(val):
#                     return -9
#                 else:
#                     return val  # 保持原值 | Keep original value
            
#             converted_phenotype = df[phenotype_col].apply(convert_quantitative_value)
            
#             # 统计数量性状的基本统计 | Basic statistics for quantitative trait
#             valid_values = converted_phenotype[converted_phenotype != -9]
#             if len(valid_values) > 0:
#                 self.logger.info(f"数量性状统计 | Quantitative trait statistics:")
#                 self.logger.info(f"  有效值数量 | Valid values: {len(valid_values)}")
#                 self.logger.info(f"  均值 | Mean: {valid_values.mean():.4f}")
#                 self.logger.info(f"  标准差 | Standard deviation: {valid_values.std():.4f}")
#                 self.logger.info(f"  最小值 | Min: {valid_values.min():.4f}")
#                 self.logger.info(f"  最大值 | Max: {valid_values.max():.4f}")
#                 self.logger.info(f"  缺失值数量 | Missing values: {(converted_phenotype == -9).sum()}")
        
#         # 创建PLINK格式的表型文件 | Create PLINK format phenotype file
#         plink_pheno = pd.DataFrame({
#             'FID': df[sample_col],
#             'IID': df[sample_col],
#             'Phenotype': converted_phenotype
#         })
        
#         # 保存转换后的表型文件 | Save converted phenotype file
#         plink_pheno.to_csv("phenotype_formatted.txt", sep=' ', index=False, header=False)
#         self.logger.info("表型文件转换完成 | Phenotype file conversion completed")
        
#         return plink_pheno

# class VCFProcessor:
#     """VCF文件处理器 | VCF File Processor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def convert_vcf_to_plink(self) -> pd.DataFrame:
#         """转换VCF文件为PLINK格式 | Convert VCF file to PLINK format"""
#         self.logger.info("转换VCF文件为PLINK格式 | Converting VCF file to PLINK format...")
        
#         cmd = [
#             "plink",
#             "--vcf", "input.vcf.gz",
#             "--make-bed",
#             "--out", "raw_data",
#             "--allow-extra-chr",
#             "--set-missing-var-ids", "@:#",
#             "--keep-allele-order"
#         ]
        
#         if self.config.threads > 1:
#             cmd.extend(["--threads", str(self.config.threads)])
        
#         self.cmd_runner.run(cmd, "VCF转换为PLINK格式 | Converting VCF to PLINK format")
        
#         # 检查转换结果 | Check conversion results
#         bed_file = Path("raw_data.bed")
#         if not bed_file.exists():
#             raise RuntimeError("VCF转换失败 | VCF conversion failed")
        
#         # 统计染色体信息 | Count chromosome information
#         bim_df = pd.read_csv("raw_data.bim", sep="\t", header=None,
#                            names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
#         chr_counts = bim_df['CHR'].value_counts()
#         self.logger.info(f"染色体分布 | Chromosome distribution: {dict(chr_counts.head(10))}")
        
#         return bim_df
    
#     def merge_phenotype(self):
#         """合并表型信息 | Merge phenotype information"""
#         self.logger.info("合并表型信息 | Merging phenotype information...")
        
#         cmd = [
#             "plink",
#             "--bfile", "raw_data",
#             "--pheno", "phenotype_formatted.txt",
#             "--make-bed",
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "data_with_pheno"
#         ]
        
#         self.cmd_runner.run(cmd, "合并表型信息 | Merging phenotype information")
        
#         bed_file = Path("data_with_pheno.bed")
#         if not bed_file.exists():
#             raise RuntimeError("表型合并失败 | Phenotype merging failed")
        
#         self.logger.info("表型合并成功 | Phenotype merging successful")

# class QualityController:
#     """质量控制器 | Quality Controller"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def quality_control(self):
#         """数据质量控制 | Data quality control"""
#         self.logger.info("开始数据质量控制 | Starting data quality control...")
        
#         # 计算基本统计 | Calculate basic statistics
#         self.logger.info("计算基本统计 | Calculating basic statistics...")
#         cmd = [
#             "plink",
#             "--bfile", "data_with_pheno",
#             "--freq",
#             "--missing", 
#             "--hardy",
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "qc_stats"
#         ]
#         self.cmd_runner.run(cmd, "计算基本统计 | Calculating basic statistics")
        
#         # 质量控制过滤 | Apply quality control filters
#         self.logger.info("应用质量控制过滤 | Applying quality control filters...")
#         cmd = [
#             "plink",
#             "--bfile", "data_with_pheno",
#             "--mind", str(self.config.mind),
#             "--geno", str(self.config.geno),
#             "--maf", str(self.config.maf),
#             "--hwe", str(self.config.hwe),
#             "--make-bed",
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "data_qc1"
#         ]
        
#         if self.config.threads > 1:
#             cmd.extend(["--threads", str(self.config.threads)])
        
#         self.cmd_runner.run(cmd, "质量控制过滤 | Quality control filtering")
        
#         bed_file = Path("data_qc1.bed")
#         if not bed_file.exists():
#             raise RuntimeError("质量控制失败 | Quality control failed")
        
#         # 统计质控结果 | Count QC results
#         bim_qc = pd.read_csv("data_qc1.bim", sep="\t", header=None)
#         fam_qc = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
        
#         self.logger.info(f"质控后SNP数 | SNPs after QC: {len(bim_qc)}")
#         self.logger.info(f"质控后样本数 | Samples after QC: {len(fam_qc)}")

# 20250726添加显性模型和隐性模型分析
"""
PLINK GWAS数据处理模块 | PLINK GWAS Data Processing Module
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from .utils import CommandRunner

class PhenotypeProcessor:
    """表型处理器 | Phenotype Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_phenotype(self):
        """转换表型文件 | Convert phenotype file"""
        self.logger.info("开始转换表型文件 | Starting phenotype file conversion...")
        self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
        
        # 读取原始表型文件 | Read original phenotype file
        try:
            pheno_df = pd.read_csv("phenotype.txt", sep='\t')
            self.logger.info(f"成功读取表型文件 | Successfully read phenotype file: {pheno_df.shape}")
            self.logger.info(f"表型文件列名 | Phenotype file columns: {list(pheno_df.columns)}")
        except Exception as e:
            self.logger.error(f"读取表型文件失败 | Failed to read phenotype file: {e}")
            raise
        
        # 检查表型文件格式 | Check phenotype file format
        if pheno_df.shape[1] < 2:
            raise ValueError("表型文件至少需要2列（样本ID和表型值） | Phenotype file needs at least 2 columns (sample ID and phenotype value)")
        
        # 确定样本ID列和表型列 | Determine sample ID column and phenotype column
        sample_col = pheno_df.columns[0]  # 第一列作为样本ID | First column as sample ID
        phenotype_col = pheno_df.columns[1]  # 第二列作为表型值 | Second column as phenotype value
        
        self.logger.info(f"样本ID列 | Sample ID column: {sample_col}")
        self.logger.info(f"表型值列 | Phenotype column: {phenotype_col}")
        
        # 统计原始表型分布 | Count original phenotype distribution
        pheno_counts = pheno_df[phenotype_col].value_counts()
        self.logger.info(f"原始表型分布 | Original phenotype distribution: {dict(pheno_counts)}")
        
        # 根据表型类型进行转换 | Convert based on trait type
        if self.config.trait_type == "qualitative":
            converted_df = self._convert_qualitative_phenotype(pheno_df, sample_col, phenotype_col)
        else:
            converted_df = self._convert_quantitative_phenotype(pheno_df, sample_col, phenotype_col)
        
        # 保存转换后的表型文件 | Save converted phenotype file
        converted_df.to_csv("phenotype_converted.txt", sep=' ', index=False, header=False)
        self.logger.info("表型文件转换完成 | Phenotype file conversion completed")
        
        # 输出转换统计 | Output conversion statistics
        self._report_conversion_statistics(pheno_df, converted_df, sample_col, phenotype_col)
        
        return converted_df
    
    def _convert_qualitative_phenotype(self, pheno_df, sample_col, phenotype_col):
        """转换质量性状表型 | Convert qualitative trait phenotype"""
        self.logger.info("处理质量性状：将0转换为1（对照），1转换为2（病例） | Processing qualitative trait: converting 0 to 1 (control), 1 to 2 (case)")
        
        def convert_qualitative_value(val):
            """转换质量性状值 | Convert qualitative trait values"""
            if pd.isna(val):
                return -9
            elif val == 0:
                return 1  # 抗病 -> 对照 | Resistant -> Control
            elif val == 1:
                return 2  # 感病 -> 病例 | Susceptible -> Case
            else:
                return -9  # 其他值设为缺失 | Other values as missing
        
        converted_phenotype = pheno_df[phenotype_col].apply(convert_qualitative_value)
        
        # 创建PLINK格式的表型文件 | Create PLINK format phenotype file
        plink_pheno = pd.DataFrame({
            'FID': pheno_df[sample_col],
            'IID': pheno_df[sample_col],
            'Phenotype': converted_phenotype
        })
        
        # 统计转换后的分布 | Count converted distribution
        converted_counts = converted_phenotype.value_counts()
        self.logger.info(f"转换后表型分布 | Converted phenotype distribution: {dict(converted_counts)}")
        self.logger.info(f"对照数 (抗病, 原值0) | Controls (resistant, original 0): {(converted_phenotype == 1).sum()}")
        self.logger.info(f"病例数 (感病, 原值1) | Cases (susceptible, original 1): {(converted_phenotype == 2).sum()}")
        
        return plink_pheno
    
    def _convert_quantitative_phenotype(self, pheno_df, sample_col, phenotype_col):
        """转换数量性状表型 | Convert quantitative trait phenotype"""
        self.logger.info("处理数量性状：保持原始数值，只处理缺失值 | Processing quantitative trait: keeping original values, only handling missing values")
        
        def convert_quantitative_value(val):
            """转换数量性状值 | Convert quantitative trait values"""
            if pd.isna(val):
                return -9
            else:
                return val  # 保持原值 | Keep original value
        
        converted_phenotype = pheno_df[phenotype_col].apply(convert_quantitative_value)
        
        # 创建PLINK格式的表型文件 | Create PLINK format phenotype file
        plink_pheno = pd.DataFrame({
            'FID': pheno_df[sample_col],
            'IID': pheno_df[sample_col],
            'Phenotype': converted_phenotype
        })
        
        # 统计数量性状的基本统计 | Basic statistics for quantitative trait
        valid_values = converted_phenotype[converted_phenotype != -9]
        if len(valid_values) > 0:
            self.logger.info(f"数量性状统计 | Quantitative trait statistics:")
            self.logger.info(f"  有效值数量 | Valid values: {len(valid_values)}")
            self.logger.info(f"  均值 | Mean: {valid_values.mean():.4f}")
            self.logger.info(f"  标准差 | Standard deviation: {valid_values.std():.4f}")
            self.logger.info(f"  最小值 | Min: {valid_values.min():.4f}")
            self.logger.info(f"  最大值 | Max: {valid_values.max():.4f}")
            self.logger.info(f"  缺失值数量 | Missing values: {(converted_phenotype == -9).sum()}")
        
        return plink_pheno
    
    def _report_conversion_statistics(self, original_df, converted_df, sample_col, phenotype_col):
        """报告转换统计信息 | Report conversion statistics"""
        self.logger.info("=== 表型转换统计 | Phenotype Conversion Statistics ===")
        
        # 原始数据统计 | Original data statistics
        original_counts = original_df[phenotype_col].value_counts(dropna=False).sort_index()
        self.logger.info(f"原始表型分布 | Original phenotype distribution:")
        for value, count in original_counts.items():
            self.logger.info(f"  {value}: {count}")
        
        # 转换后数据统计 | Converted data statistics
        converted_counts = converted_df['Phenotype'].value_counts(dropna=False).sort_index()
        self.logger.info(f"转换后表型分布 | Converted phenotype distribution:")
        for value, count in converted_counts.items():
            self.logger.info(f"  {value}: {count}")
        
        total_samples = len(original_df)
        valid_original = original_df[phenotype_col].notna().sum()
        valid_converted = len(converted_df[converted_df['Phenotype'] != -9])
        
        self.logger.info(f"总样本数 | Total samples: {total_samples}")
        self.logger.info(f"原始有效样本数 | Original valid samples: {valid_original}")
        self.logger.info(f"转换后有效样本数 | Converted valid samples: {valid_converted}")
        
        if self.config.trait_type == "qualitative":
            # 质量性状特定统计 | Qualitative trait specific statistics
            controls = (converted_df['Phenotype'] == 1).sum()
            cases = (converted_df['Phenotype'] == 2).sum()
            missing = (converted_df['Phenotype'] == -9).sum()
            
            self.logger.info(f"转换后统计 | Post-conversion statistics:")
            self.logger.info(f"  对照组(抗病) | Controls(resistant): {controls}")
            self.logger.info(f"  病例组(感病) | Cases(susceptible): {cases}")
            self.logger.info(f"  缺失值 | Missing: {missing}")
            
            if controls > 0 and cases > 0:
                case_control_ratio = cases / controls
                self.logger.info(f"  病例/对照比例 | Case/control ratio: {case_control_ratio:.3f}")
                
                # 检查样本量是否足够 | Check if sample size is sufficient
                min_group_size = min(controls, cases)
                if min_group_size < 10:
                    self.logger.warning(f"样本量较小 | Small sample size: 最小组仅有{min_group_size}个样本 | minimum group has only {min_group_size} samples")
                elif min_group_size < 50:
                    self.logger.info(f"样本量适中 | Moderate sample size: 最小组有{min_group_size}个样本 | minimum group has {min_group_size} samples")
                else:
                    self.logger.info(f"样本量充足 | Sufficient sample size: 最小组有{min_group_size}个样本 | minimum group has {min_group_size} samples")
        
        # 保存转换过程的详细信息用于调试 | Save detailed conversion info for debugging
        debug_info = pd.DataFrame({
            'Original_Sample': original_df[sample_col],
            'Original_Phenotype': original_df[phenotype_col],
            'Converted_FID': converted_df['FID'],
            'Converted_IID': converted_df['IID'],
            'Converted_Phenotype': converted_df['Phenotype']
        })
        debug_info.to_csv("phenotype_conversion_debug.txt", sep='\t', index=False)
        self.logger.info("转换详情已保存到 phenotype_conversion_debug.txt | Conversion details saved to phenotype_conversion_debug.txt")


class VCFProcessor:
    """VCF文件处理器 | VCF File Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_vcf_to_plink(self):
        """转换VCF文件为PLINK格式 | Convert VCF file to PLINK format"""
        self.logger.info("转换VCF文件为PLINK格式 | Converting VCF file to PLINK format...")
        
        # 检查输入文件 | Check input file
        input_vcf = "input.vcf.gz" if os.path.exists("input.vcf.gz") else "input.vcf"
        if not os.path.exists(input_vcf):
            raise FileNotFoundError(f"VCF文件不存在 | VCF file does not exist: {input_vcf}")
        
        cmd = [
            "plink",
            "--vcf", input_vcf,
            "--make-bed",
            "--out", "raw_data",
            "--allow-extra-chr",
            "--set-missing-var-ids", "@:#",
            "--keep-allele-order"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(cmd, "VCF转换为PLINK格式 | Converting VCF to PLINK format")
        
        # 检查转换结果 | Check conversion results
        bed_file = Path("raw_data.bed")
        if not bed_file.exists():
            raise RuntimeError("VCF转换失败 | VCF conversion failed")
        
        # 统计染色体信息 | Count chromosome information
        if os.path.exists("raw_data.bim"):
            bim_df = pd.read_csv("raw_data.bim", sep="\t", header=None,
                               names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
            chr_counts = bim_df['CHR'].value_counts()
            self.logger.info(f"染色体分布 | Chromosome distribution: {dict(chr_counts.head(10))}")
            self.logger.info(f"总SNP数 | Total SNPs: {len(bim_df)}")
        
        self.logger.info("VCF转换完成 | VCF conversion completed")
    
    def merge_phenotype(self):
        """合并表型信息 | Merge phenotype information"""
        self.logger.info("合并表型信息 | Merging phenotype information...")
        
        # 检查必要文件是否存在 | Check if necessary files exist
        if not os.path.exists("raw_data.bed"):
            raise FileNotFoundError("PLINK数据文件不存在 | PLINK data files do not exist: raw_data.*")
        
        if not os.path.exists("phenotype_converted.txt"):
            raise FileNotFoundError("转换后的表型文件不存在 | Converted phenotype file does not exist")
        
        cmd = [
            "plink",
            "--bfile", "raw_data",
            "--pheno", "phenotype_converted.txt",
            "--make-bed",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "data_with_pheno"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(cmd, "合并表型信息 | Merging phenotype information")
        
        # 检查合并结果 | Check merge results
        bed_file = Path("data_with_pheno.bed")
        if not bed_file.exists():
            raise RuntimeError("表型合并失败 | Phenotype merging failed")
        
        # 统计合并后的数据 | Count merged data
        if os.path.exists("data_with_pheno.fam"):
            fam_df = pd.read_csv("data_with_pheno.fam", sep=r"\s+", header=None)
            self.logger.info(f"合并后样本数 | Samples after merge: {len(fam_df)}")
            
            # 统计表型分布 | Count phenotype distribution
            pheno_counts = fam_df.iloc[:, 5].value_counts()
            self.logger.info(f"表型分布 | Phenotype distribution: {dict(pheno_counts)}")
        
        self.logger.info("表型合并成功 | Phenotype merging successful")


class QualityController:
    """质量控制器 | Quality Controller"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def quality_control(self):
        """数据质量控制 | Data quality control"""
        self.logger.info("开始数据质量控制 | Starting data quality control...")
        
        # 检查输入文件 | Check input files
        if not os.path.exists("data_with_pheno.bed"):
            raise FileNotFoundError("质控输入文件不存在 | QC input files do not exist: data_with_pheno.*")
        
        # 计算基本统计 | Calculate basic statistics
        self.logger.info("计算基本统计 | Calculating basic statistics...")
        cmd = [
            "plink",
            "--bfile", "data_with_pheno",
            "--freq",
            "--missing", 
            "--hardy",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "qc_stats"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        # 使用 check=False 来避免因统计计算失败而中断流程
        result = self.cmd_runner.run(cmd, "计算基本统计 | Calculating basic statistics", check=False)
        
        if result.returncode != 0:
            self.logger.warning("基本统计计算失败，继续质控过滤 | Basic statistics calculation failed, proceeding with QC filtering")
        
        # 质量控制过滤 | Apply quality control filters
        self.logger.info("应用质量控制过滤 | Applying quality control filters...")
        self.logger.info(f"质控参数 | QC parameters: mind={self.config.mind}, geno={self.config.geno}, maf={self.config.maf}, hwe={self.config.hwe}")
        
        cmd = [
            "plink",
            "--bfile", "data_with_pheno",
            "--mind", str(self.config.mind),
            "--geno", str(self.config.geno),
            "--maf", str(self.config.maf),
            "--hwe", str(self.config.hwe),
            "--make-bed",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "data_qc1"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(cmd, "质量控制过滤 | Quality control filtering")
        
        # 检查质控结果 | Check QC results
        bed_file = Path("data_qc1.bed")
        if not bed_file.exists():
            raise RuntimeError("质量控制失败 | Quality control failed")
        
        # 统计质控结果 | Count QC results
        self._report_qc_statistics()
        
        self.logger.info("数据质量控制完成 | Data quality control completed")
    
    def _report_qc_statistics(self):
        """报告质控统计信息 | Report QC statistics"""
        try:
            # 质控前统计 | Pre-QC statistics
            if os.path.exists("data_with_pheno.bim") and os.path.exists("data_with_pheno.fam"):
                bim_before = pd.read_csv("data_with_pheno.bim", sep="\t", header=None)
                fam_before = pd.read_csv("data_with_pheno.fam", sep=r"\s+", header=None)
                snps_before = len(bim_before)
                samples_before = len(fam_before)
            else:
                snps_before = "未知"
                samples_before = "未知"
            
            # 质控后统计 | Post-QC statistics
            if os.path.exists("data_qc1.bim") and os.path.exists("data_qc1.fam"):
                bim_after = pd.read_csv("data_qc1.bim", sep="\t", header=None)
                fam_after = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
                snps_after = len(bim_after)
                samples_after = len(fam_after)
                
                # 统计染色体数 | Count chromosomes
                chromosomes = bim_after.iloc[:, 0].nunique()
                
                # 统计表型分布 | Count phenotype distribution
                if self.config.trait_type == "qualitative":
                    cases = (fam_after.iloc[:, 5] == 2).sum()
                    controls = (fam_after.iloc[:, 5] == 1).sum()
                    missing_pheno = (fam_after.iloc[:, 5] == -9).sum()
                    
                    self.logger.info("=== 质控统计结果 | QC Statistics Results ===")
                    self.logger.info(f"质控前SNP数 | SNPs before QC: {snps_before}")
                    self.logger.info(f"质控后SNP数 | SNPs after QC: {snps_after}")
                    self.logger.info(f"质控前样本数 | Samples before QC: {samples_before}")
                    self.logger.info(f"质控后样本数 | Samples after QC: {samples_after}")
                    self.logger.info(f"染色体数 | Number of chromosomes: {chromosomes}")
                    self.logger.info(f"病例数(感病) | Cases(susceptible): {cases}")
                    self.logger.info(f"对照数(抗病) | Controls(resistant): {controls}")
                    self.logger.info(f"表型缺失 | Missing phenotype: {missing_pheno}")
                    
                    if isinstance(snps_before, int) and isinstance(samples_before, int):
                        snp_retention = (snps_after / snps_before) * 100
                        sample_retention = (samples_after / samples_before) * 100
                        self.logger.info(f"SNP保留率 | SNP retention rate: {snp_retention:.2f}%")
                        self.logger.info(f"样本保留率 | Sample retention rate: {sample_retention:.2f}%")
                
            else:
                self.logger.warning("无法读取质控后文件进行统计 | Cannot read post-QC files for statistics")
                
        except Exception as e:
            self.logger.warning(f"质控统计报告生成失败 | QC statistics report generation failed: {e}")