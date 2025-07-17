"""
ADMIXTURE数据处理模块 | ADMIXTURE Data Processing Module
"""

import os
import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class VCFProcessor:
    """VCF处理器 | VCF Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def preprocess_vcf(self):
        """预处理VCF文件 | Preprocess VCF file"""
        if self.config.skip_preprocessing:
            self.logger.info("跳过VCF预处理 | Skipping VCF preprocessing")
            return
        
        vcf_file = self.config.vcf_file
        output_dir = self.config.output_dir
        
        self.logger.info(f"预处理VCF文件 | Preprocessing VCF file: {vcf_file}")
        
        # 检查VCF文件格式 | Check VCF file format
        self._check_vcf_format(vcf_file)
        
        # 移除多等位基因位点 | Remove multi-allelic sites
        biallelic_vcf = os.path.join(output_dir, "biallelic.vcf.gz")
        self._remove_multiallelic(vcf_file, biallelic_vcf)
        
        # 检查染色体命名 | Check chromosome naming
        self._check_chromosome_naming(biallelic_vcf)
        
        return biallelic_vcf
    
    def _check_vcf_format(self, vcf_file: str):
        """检查VCF文件格式 | Check VCF file format"""
        if vcf_file.endswith('.gz'):
            cmd = f"zcat {vcf_file} | head -20"
        else:
            cmd = f"head -20 {vcf_file}"
        
        output = self.cmd_runner.run(cmd, "检查VCF文件头信息 | Check VCF header")
        self.logger.info(f"VCF文件头信息 | VCF header info:\n{output}")
    
    def _remove_multiallelic(self, input_vcf: str, output_vcf: str):
        """移除多等位基因位点 | Remove multi-allelic sites"""
        cmd = f"bcftools view -m2 -M2 -v snps {input_vcf} -Oz -o {output_vcf}"
        self.cmd_runner.run(cmd, "移除多等位基因位点 | Remove multi-allelic sites")
    
    def _check_chromosome_naming(self, vcf_file: str):
        """检查染色体命名 | Check chromosome naming"""
        if vcf_file.endswith('.gz'):
            cmd = f"zcat {vcf_file} | grep -v '^#' | cut -f1 | sort | uniq -c"
        else:
            cmd = f"grep -v '^#' {vcf_file} | cut -f1 | sort | uniq -c"
        
        output = self.cmd_runner.run(cmd, "检查染色体命名 | Check chromosome naming")
        self.logger.info(f"染色体信息 | Chromosome info:\n{output}")

class PlinkProcessor:
    """PLINK处理器 | PLINK Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_vcf_to_plink(self, vcf_file: str):
        """VCF转换为PLINK格式 | Convert VCF to PLINK format"""
        output_prefix = os.path.join(self.config.output_dir, "raw_data")
        
        cmd = (
            f"plink --vcf {vcf_file} --make-bed --out {output_prefix} "
            f"--allow-extra-chr --double-id"
        )
        
        self.cmd_runner.run(cmd, "VCF转换为PLINK格式 | Convert VCF to PLINK format")
        return output_prefix
    
    def quality_control(self, input_prefix: str):
        """质量控制 | Quality control"""
        output_prefix = os.path.join(self.config.output_dir, self.config.base_name)
        
        # MAF过滤 | MAF filtering
        maf_prefix = f"{output_prefix}_maf"
        cmd_maf = (
            f"plink --bfile {input_prefix} --maf {self.config.maf} "
            f"--make-bed --out {maf_prefix} --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_maf, f"MAF过滤 | MAF filtering (>={self.config.maf})")
        
        # HWE过滤 | HWE filtering
        hwe_prefix = f"{output_prefix}_hwe"
        cmd_hwe = (
            f"plink --bfile {maf_prefix} --hwe {self.config.hwe_pvalue} "
            f"--make-bed --out {hwe_prefix} --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_hwe, f"HWE过滤 | HWE filtering (p>{self.config.hwe_pvalue})")
        
        # 缺失率过滤 | Missing rate filtering
        cmd_missing = (
            f"plink --bfile {hwe_prefix} "
            f"--geno {self.config.missing_rate} --mind {self.config.missing_rate} "
            f"--make-bed --out {output_prefix} --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_missing, f"缺失率过滤 | Missing rate filtering (<{self.config.missing_rate})")
        
        # 清理中间文件 | Clean intermediate files
        if not self.config.keep_intermediate:
            self._cleanup_intermediate_files([maf_prefix, hwe_prefix])
        
        return output_prefix
    
    def _cleanup_intermediate_files(self, prefixes: list):
        """清理中间文件 | Clean intermediate files"""
        for prefix in prefixes:
            for ext in ['.bed', '.bim', '.fam', '.log']:
                file_path = f"{prefix}{ext}"
                if os.path.exists(file_path):
                    os.remove(file_path)
        
        self.logger.info("清理中间文件完成 | Intermediate files cleaned")
