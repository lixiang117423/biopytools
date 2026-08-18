"""
变异过滤模块|Variant Filtering Module
"""

from pathlib import Path
from .utils import CommandRunner
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command

class VariantFilter:
    """变异过滤器|Variant Filter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def select_snps(self):
        """提取SNP|Extract SNPs"""
        self.logger.info(" 提取SNP变异|Extracting SNP variants")
        
        snp_vcf = self.config.output_path / f"{self.config.base_name}_snps_raw.vcf.gz"
        
        args = [
            "--java-options", f"-Xmx{self.config.memory}",
            "SelectVariants",
            "-R", self.config.reference,
            "-V", self.config.raw_vcf,
            "-select-type", "SNP",
            "-O", str(snp_vcf),
        ]
        cmd = build_conda_command(self.config.gatk_path, args)
        
        success = self.cmd_runner.run(cmd, "提取SNP|Extract SNPs")
        
        if success:
            self.config.snp_raw_vcf = str(snp_vcf)
        
        return success
    
    def select_indels(self):
        """提取INDEL|Extract INDELs"""
        self.logger.info(" 提取INDEL变异|Extracting INDEL variants")
        
        indel_vcf = self.config.output_path / f"{self.config.base_name}_indels_raw.vcf.gz"
        
        args = [
            "--java-options", f"-Xmx{self.config.memory}",
            "SelectVariants",
            "-R", self.config.reference,
            "-V", self.config.raw_vcf,
            "-select-type", "INDEL",
            "-O", str(indel_vcf),
        ]
        cmd = build_conda_command(self.config.gatk_path, args)
        
        success = self.cmd_runner.run(cmd, "提取INDEL|Extract INDELs")
        
        if success:
            self.config.indel_raw_vcf = str(indel_vcf)
        
        return success
    
    def filter_snps(self):
        """过滤SNP|Filter SNPs"""
        self.logger.info(" 过滤SNP变异|Filtering SNP variants")
        
        snp_filtered = self.config.output_path / f"{self.config.base_name}_snps_filtered.vcf.gz"
        
        # 构建过滤表达式|Build filter expression
        filters = [
            (f"QD < {self.config.snp_qd}", "QD_filter"),
            (f"FS > {self.config.snp_fs}", "FS_filter"),
            (f"MQ < {self.config.snp_mq}", "MQ_filter"),
            (f"MQRankSum < {self.config.snp_mqrs}", "MQRS_filter"),
            (f"ReadPosRankSum < {self.config.snp_rprs}", "RPRS_filter"),
            (f"SOR > {self.config.snp_sor}", "SOR_filter"),
        ]

        args = [
            "--java-options", f"-Xmx{self.config.memory}",
            "VariantFiltration",
            "-R", self.config.reference,
            "-V", self.config.snp_raw_vcf,
            "-O", str(snp_filtered),
        ]
        
        for expression, name in filters:
            args += ["--filter-expression", expression, "--filter-name", name]
        
        cmd = build_conda_command(self.config.gatk_path, args)
        
        success = self.cmd_runner.run(cmd, "过滤SNP|Filter SNPs")
        
        if success:
            self.config.snp_filtered_vcf = str(snp_filtered)
        
        return success
    
    def filter_indels(self):
        """过滤INDEL|Filter INDELs"""
        self.logger.info(" 过滤INDEL变异|Filtering INDEL variants")
        
        indel_filtered = self.config.output_path / f"{self.config.base_name}_indels_filtered.vcf.gz"
        
        # 构建过滤表达式|Build filter expression
        filters = [
            (f"QD < {self.config.indel_qd}", "QD_filter"),
            (f"FS > {self.config.indel_fs}", "FS_filter"),
            (f"ReadPosRankSum < {self.config.indel_rprs}", "RPRS_filter"),
            (f"SOR > {self.config.indel_sor}", "SOR_filter"),
        ]
        
        args = [
            "--java-options", f"-Xmx{self.config.memory}",
            "VariantFiltration",
            "-R", self.config.reference,
            "-V", self.config.indel_raw_vcf,
            "-O", str(indel_filtered),
        ]

        for expression, name in filters:
            args += ["--filter-expression", expression, "--filter-name", name]
        
        cmd = build_conda_command(self.config.gatk_path, args)
        
        success = self.cmd_runner.run(cmd, "过滤INDEL|Filter INDELs")
        
        if success:
            self.config.indel_filtered_vcf = str(indel_filtered)
        
        return success
    
    def merge_variants(self):
        """合并SNP和INDEL|Merge SNPs and INDELs"""
        self.logger.info(" 合并SNP和INDEL|Merging SNPs and INDELs")
        
        merged_vcf = self.config.output_path / f"{self.config.base_name}_merged_filtered.vcf.gz"
        
        # 使用bcftools合并|Use bcftools to merge
        args = [
            "concat",
            "-a", "-O", "z",
            "-o", str(merged_vcf),
            self.config.snp_filtered_vcf,
            self.config.indel_filtered_vcf,
        ]
        cmd = build_conda_command(self.config.bcftools_path, args)
        
        success = self.cmd_runner.run(cmd, "合并变异|Merge variants")
        
        if success:
            self.config.merged_vcf = str(merged_vcf)
            # 建立索引|Create index
            self._index_vcf(merged_vcf)
        
        return success
    
    def _index_vcf(self, vcf_file):
        """为VCF文件建立索引|Index VCF file"""
        cmd = build_conda_command(self.config.bcftools_path, ["index", "-t", str(vcf_file)])
        self.cmd_runner.run(cmd, f"索引VCF|Index VCF: {vcf_file}")
