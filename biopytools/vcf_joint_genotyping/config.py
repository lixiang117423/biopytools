"""
VCF联合分型配置模块 | VCF Joint Genotyping Configuration Module
"""

import os
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

@dataclass
class VCFJointGenotypingConfig:
    """VCF联合分型配置类 | VCF Joint Genotyping Configuration Class"""
    
    # 基本输入输出参数 | Basic I/O parameters
    vcf_input_dir: str
    output_dir: str
    reference_genome: str
    
    # 软件路径 | Software paths
    gtx_path: str = "gtx"
    gatk_path: str = "gatk"
    vcftools_path: str = "vcftools"
    bgzip_path: str = "bgzip"
    tabix_path: str = "tabix"
    bcftools_path: str = "bcftools"
    
    # SNP过滤参数 | SNP filtering parameters
    snp_maf: float = 0.05
    snp_max_missing: float = 0.5
    snp_hwe_pvalue: float = 1e-6
    snp_min_mean_dp: int = 5
    snp_max_mean_dp: int = 50
    
    # INDEL过滤参数 | INDEL filtering parameters
    indel_maf: float = 0.05  # 注意：原脚本中是1，这里改为合理默认值
    indel_max_missing: float = 1.0
    indel_hwe_pvalue: float = 1e-6
    indel_min_mean_dp: int = 5
    indel_max_mean_dp: int = 50
    
    # 控制参数 | Control parameters
    skip_joint: bool = False
    skip_extract: bool = False
    skip_filter: bool = False
    skip_compress: bool = False
    skip_corrupted_files: bool = True  # 默认跳过损坏文件而不是终止程序
    
    # 输出文件前缀 | Output file prefixes
    merged_prefix: str = "all.merged"
    snp_prefix: str = "all.merged.snp"
    indel_prefix: str = "all.merged.indel"
    filtered_snp_prefix: str = "final.filtered.snp"
    filtered_indel_prefix: str = "final.filtered.indel"
    
    # 系统参数 | System parameters
    memory: str = "128g"
    threads: int = 1
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 确保路径使用绝对路径 | Ensure paths use absolute paths
        self.vcf_input_dir = str(Path(self.vcf_input_dir).resolve())
        self.output_dir = str(Path(self.output_dir).resolve())
        self.reference_genome = str(Path(self.reference_genome).resolve())
        
        # 创建输出目录 | Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
    
    def get_output_paths(self):
        """获取输出文件路径 | Get output file paths"""
        output_base = Path(self.output_dir)
        
        paths = {
            'sample_map_file': output_base / "sample_vcf_map.txt",
            'merged_vcf': output_base / f"{self.merged_prefix}.vcf.gz",
            'snp_vcf': output_base / f"{self.snp_prefix}.vcf.gz",
            'indel_vcf': output_base / f"{self.indel_prefix}.vcf.gz",
            'filtered_snp_vcf': output_base / f"{self.filtered_snp_prefix}.recode.vcf",
            'filtered_indel_vcf': output_base / f"{self.filtered_indel_prefix}.recode.vcf",
            'final_snp_vcf': output_base / f"{self.filtered_snp_prefix}.recode.vcf.gz",
            'final_indel_vcf': output_base / f"{self.filtered_indel_prefix}.recode.vcf.gz",
            'report_file': output_base / "processing_report.txt"
        }
        
        return paths
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需的文件和目录 | Check required files and directories
        if not Path(self.vcf_input_dir).exists():
            errors.append(f"VCF输入目录不存在 | VCF input directory does not exist: {self.vcf_input_dir}")
        
        if not Path(self.reference_genome).exists():
            errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {self.reference_genome}")
        
        # 检查参数范围 | Check parameter ranges
        if not (0 <= self.snp_maf <= 1):
            errors.append(f"SNP MAF必须在0-1之间 | SNP MAF must be between 0-1: {self.snp_maf}")
        
        if not (0 <= self.snp_max_missing <= 1):
            errors.append(f"SNP最大缺失率必须在0-1之间 | SNP max missing rate must be between 0-1: {self.snp_max_missing}")
        
        if not (0 <= self.indel_maf <= 1):
            errors.append(f"INDEL MAF必须在0-1之间 | INDEL MAF must be between 0-1: {self.indel_maf}")
        
        if not (0 <= self.indel_max_missing <= 1):
            errors.append(f"INDEL最大缺失率必须在0-1之间 | INDEL max missing rate must be between 0-1: {self.indel_max_missing}")
        
        if self.snp_min_mean_dp >= self.snp_max_mean_dp:
            errors.append(f"SNP最小深度必须小于最大深度 | SNP min depth must be less than max depth")
        
        if self.indel_min_mean_dp >= self.indel_max_mean_dp:
            errors.append(f"INDEL最小深度必须小于最大深度 | INDEL min depth must be less than max depth")
        
        return errors
