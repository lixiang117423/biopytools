"""
VCF数据处理模块|VCF Data Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class VCFProcessor:
    """VCF文件处理器|VCF File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def is_compressed_vcf(self, vcf_file: str) -> bool:
        """检查VCF文件是否压缩|Check if VCF file is compressed"""
        return vcf_file.endswith('.gz')
    
    def get_vcf_stats(self, vcf_file: str) -> Dict[str, int]:
        """获取VCF文件统计信息|Get VCF file statistics"""
        self.logger.info("获取VCF文件统计信息|Getting VCF file statistics")
        
        stats = {}
        
        # 获取样本数量|Get sample count
        if self.is_compressed_vcf(vcf_file):
            cmd = f"{self.config.bcftools_path} query -l {vcf_file}|wc -l"
        else:
            cmd = f"grep '^#CHROM' {vcf_file}|cut -f10-|tr '\\t' '\\n'|wc -l"
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            stats['samples'] = int(result.stdout.strip())
        except:
            stats['samples'] = 0
        
        # 获取变异数量|Get variant count
        if self.is_compressed_vcf(vcf_file):
            cmd = f"{self.config.bcftools_path} view -H {vcf_file}|wc -l"
        else:
            cmd = f"grep -v '^#' {vcf_file}|wc -l"
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            stats['variants'] = int(result.stdout.strip())
        except:
            stats['variants'] = 0
        
        self.logger.info(f"VCF统计|VCF statistics: {stats['samples']} 样本|samples, {stats['variants']} 变异|variants")
        
        return stats
    
    def vcf_to_plink(self):
        """VCF转换为PLINK格式|Convert VCF to PLINK format"""
        vcf_path = self.config.vcf_file
        output_prefix = self.config.output_path / self.config.base_name
        
        # 确保输出目录存在并获取绝对路径|Ensure output directory exists and get absolute path
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix_abs = output_prefix.resolve()
        
        # 获取VCF统计信息|Get VCF statistics
        stats = self.get_vcf_stats(vcf_path)
        self.logger.info(f"开始处理包含 {stats['samples']} 个样本和 {stats['variants']} 个变异的VCF文件")
        
        # 检查目录是否可写|Check if directory is writable
        try:
            test_file = output_prefix.parent / "test_write.tmp"
            test_file.touch()
            test_file.unlink()
            self.logger.info(f"输出目录可写|Output directory is writable: {output_prefix.parent}")
        except Exception as e:
            self.logger.error(f"输出目录不可写|Output directory is not writable: {output_prefix.parent}, 错误: {e}")
            return False
        
        cmd = (
            f"{self.config.plink_path} --vcf {vcf_path} "
            f"--make-bed --out {output_prefix_abs} "
            f"--allow-extra-chr --double-id"
        )
        
        success = self.cmd_runner.run(cmd, "VCF转换为PLINK格式|Convert VCF to PLINK format")
        
        if success:
            self.config.plink_prefix = str(output_prefix_abs)
            self.logger.info(f"PLINK文件已生成|PLINK files generated: {output_prefix_abs}")
        
        return success

class QualityController:
    """质量控制器|Quality Controller"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def apply_quality_filters(self):
        """应用质量过滤|Apply quality filters"""
        # 检查是否跳过质控|Check if skipping QC
        if self.config.skip_qc:
            self.logger.info("跳过质量控制过滤 (用户指定)|Skipping quality control filtering (user specified)")
            self.config.plink_prefix_qc = self.config.plink_prefix
            self.log_qc_stats()
            return True
        
        input_prefix = self.config.plink_prefix
        output_prefix = self.config.output_path / f"{self.config.base_name}_qc"
        
        self.logger.info("应用质量控制过滤|Applying quality control filters")
        self.logger.info(f"MAF阈值|MAF threshold: {self.config.maf}")
        self.logger.info(f"缺失率阈值|Missing rate threshold: {self.config.missing_rate}")
        self.logger.info(f"HWE p值阈值|HWE p-value threshold: {self.config.hwe_pvalue}")
        
        # Step 1: MAF过滤|MAF filtering
        maf_prefix = self.config.output_path / f"{self.config.base_name}_maf"
        maf_prefix.parent.mkdir(parents=True, exist_ok=True)
        maf_prefix_abs = maf_prefix.resolve()
        cmd_maf = (
            f"{self.config.plink_path} --bfile {input_prefix} "
            f"--maf {self.config.maf} --make-bed --out {maf_prefix_abs} --allow-extra-chr"
        )
        
        if not self.cmd_runner.run(cmd_maf, f"MAF过滤|MAF filtering (>={self.config.maf})"):
            return False
        
        # Step 2: HWE过滤|HWE filtering
        hwe_prefix = self.config.output_path / f"{self.config.base_name}_hwe"
        hwe_prefix.parent.mkdir(parents=True, exist_ok=True)
        hwe_prefix_abs = hwe_prefix.resolve()
        cmd_hwe = (
            f"{self.config.plink_path} --bfile {maf_prefix_abs} "
            f"--hwe {self.config.hwe_pvalue} --make-bed --out {hwe_prefix_abs} --allow-extra-chr"
        )
        
        if not self.cmd_runner.run(cmd_hwe, f"HWE过滤|HWE filtering (p>{self.config.hwe_pvalue})"):
            return False
        
        # Step 3: 缺失率过滤|Missing rate filtering
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix_abs = output_prefix.resolve()
        cmd_missing = (
            f"{self.config.plink_path} --bfile {hwe_prefix_abs} "
            f"--geno {self.config.missing_rate} --mind {self.config.missing_rate} "
            f"--make-bed --out {output_prefix_abs} --allow-extra-chr"
        )
        
        success = self.cmd_runner.run(cmd_missing, f"缺失率过滤|Missing rate filtering (<{self.config.missing_rate})")
        
        if success:
            self.config.plink_prefix_qc = str(output_prefix_abs)
            self.log_qc_stats()
        
        return success
        cmd_missing = (
            f"{self.config.plink_path} --bfile {hwe_prefix} "
            f"--geno {self.config.missing_rate} --mind {self.config.missing_rate} "
            f"--make-bed --out {output_prefix} --allow-extra-chr"
        )
        
        success = self.cmd_runner.run(cmd_missing, f"缺失率过滤|Missing rate filtering (<{self.config.missing_rate})")
        
        if success:
            self.config.plink_prefix_qc = str(output_prefix)
            self.log_qc_stats()
        
        return success
    
    def log_qc_stats(self):
        """记录质控统计|Log QC statistics"""
        try:
            fam_file = Path(f"{self.config.plink_prefix_qc}.fam")
            bim_file = Path(f"{self.config.plink_prefix_qc}.bim")
            
            if fam_file.exists():
                sample_count = sum(1 for line in open(fam_file))
                if self.config.skip_qc:
                    self.logger.info(f"输入样本数|Input samples: {sample_count}")
                else:
                    self.logger.info(f"质控后样本数|Samples after QC: {sample_count}")
            
            if bim_file.exists():
                snp_count = sum(1 for line in open(bim_file))
                if self.config.skip_qc:
                    self.logger.info(f"输入SNP数|Input SNPs: {snp_count}")
                else:
                    self.logger.info(f"质控后SNP数|SNPs after QC: {snp_count}")
                
        except Exception as e:
            self.logger.warning(f"无法读取统计信息|Cannot read statistics: {e}")
