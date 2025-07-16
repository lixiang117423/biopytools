"""
PLINK GWAS数据处理模块 | PLINK GWAS Data Processing Module
"""

import pandas as pd
import numpy as np
from pathlib import Path
from .utils import CommandRunner, FileManager

class PhenotypeProcessor:
    """表型数据处理器 | Phenotype Data Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_phenotype(self) -> pd.DataFrame:
        """转换表型文件格式 | Convert phenotype file format"""
        self.logger.info("转换表型文件格式 | Converting phenotype file format...")
        
        # 读取原始表型文件 | Read original phenotype file
        df = pd.read_csv("phenotype.txt", sep=r"\s+")
        self.logger.info(f"原始表型文件形状 | Original phenotype file shape: {df.shape}")
        
        # 检查列名 | Check column names
        if df.shape[1] < 2:
            raise ValueError("表型文件至少需要2列（样本ID和表型值） | Phenotype file needs at least 2 columns (sample ID and phenotype)")
        
        # 假设第一列是样本ID，第二列是表型值 | Assume first column is sample ID, second is phenotype
        sample_col = df.columns[0]
        phenotype_col = df.columns[1]
        
        # 统计原始表型分布 | Count original phenotype distribution
        pheno_counts = df[phenotype_col].value_counts()
        self.logger.info(f"原始表型分布 | Original phenotype distribution: {dict(pheno_counts)}")
        
        def convert_phenotype_value(val):
            """转换表型值 | Convert phenotype values"""
            if pd.isna(val):
                return -9
            elif val == 0:
                return 1  # 抗病 -> 对照 | Resistant -> Control
            elif val == 1:
                return 2  # 感病 -> 病例 | Susceptible -> Case
            else:
                return -9  # 其他值设为缺失 | Other values as missing
        
        # 创建PLINK格式的表型文件 | Create PLINK format phenotype file
        plink_pheno = pd.DataFrame({
            'FID': df[sample_col],
            'IID': df[sample_col],
            'Phenotype': df[phenotype_col].apply(convert_phenotype_value)
        })
        
        # 保存转换后的表型文件 | Save converted phenotype file
        plink_pheno.to_csv("phenotype_formatted.txt", sep=' ', index=False, header=False)
        
        # 统计转换后的表型分布 | Count converted phenotype distribution
        converted_counts = plink_pheno['Phenotype'].value_counts()
        self.logger.info(f"转换后表型分布 | Converted phenotype distribution: {dict(converted_counts)}")
        self.logger.info(f"对照数 (抗病) | Controls (resistant): {(plink_pheno['Phenotype'] == 1).sum()}")
        self.logger.info(f"病例数 (感病) | Cases (susceptible): {(plink_pheno['Phenotype'] == 2).sum()}")
        
        return plink_pheno

class VCFProcessor:
    """VCF文件处理器 | VCF File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_vcf_to_plink(self) -> pd.DataFrame:
        """转换VCF文件为PLINK格式 | Convert VCF file to PLINK format"""
        self.logger.info("转换VCF文件为PLINK格式 | Converting VCF file to PLINK format...")
        
        cmd = [
            "plink",
            "--vcf", "input.vcf.gz",
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
        bim_df = pd.read_csv("raw_data.bim", sep="\t", header=None,
                           names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        chr_counts = bim_df['CHR'].value_counts()
        self.logger.info(f"染色体分布 | Chromosome distribution: {dict(chr_counts.head(10))}")
        
        return bim_df
    
    def merge_phenotype(self):
        """合并表型信息 | Merge phenotype information"""
        self.logger.info("合并表型信息 | Merging phenotype information...")
        
        cmd = [
            "plink",
            "--bfile", "raw_data",
            "--pheno", "phenotype_formatted.txt",
            "--make-bed",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "data_with_pheno"
        ]
        
        self.cmd_runner.run(cmd, "合并表型信息 | Merging phenotype information")
        
        bed_file = Path("data_with_pheno.bed")
        if not bed_file.exists():
            raise RuntimeError("表型合并失败 | Phenotype merging failed")
        
        self.logger.info("表型合并成功 | Phenotype merging successful")

class QualityController:
    """质量控制器 | Quality Controller"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def quality_control(self):
        """数据质量控制 | Data quality control"""
        self.logger.info("开始数据质量控制 | Starting data quality control...")
        
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
        self.cmd_runner.run(cmd, "计算基本统计 | Calculating basic statistics")
        
        # 质量控制过滤 | Apply quality control filters
        self.logger.info("应用质量控制过滤 | Applying quality control filters...")
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
        
        bed_file = Path("data_qc1.bed")
        if not bed_file.exists():
            raise RuntimeError("质量控制失败 | Quality control failed")
        
        # 统计质控结果 | Count QC results
        bim_qc = pd.read_csv("data_qc1.bim", sep="\t", header=None)
        fam_qc = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
        
        self.logger.info(f"质控后SNP数 | SNPs after QC: {len(bim_qc)}")
        self.logger.info(f"质控后样本数 | Samples after QC: {len(fam_qc)}")
