"""
变异过滤模块 | Variant Filtering Module
"""

from pathlib import Path
from .utils import CommandRunner

class VariantFilter:
    """变异过滤器 | Variant Filter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def select_snps(self):
        """提取SNP | Extract SNPs"""
        self.logger.info("🔬 提取SNP变异 | Extracting SNP variants")
        
        snp_vcf = self.config.output_path / f"{self.config.base_name}_snps_raw.vcf.gz"
        
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' SelectVariants"
        cmd += f" -R {self.config.reference}"
        cmd += f" -V {self.config.raw_vcf}"
        cmd += f" -select-type SNP"
        cmd += f" -O {snp_vcf}"
        
        success = self.cmd_runner.run(cmd, "提取SNP | Extract SNPs")
        
        if success:
            self.config.snp_raw_vcf = str(snp_vcf)
        
        return success
    
    def select_indels(self):
        """提取INDEL | Extract INDELs"""
        self.logger.info("🔬 提取INDEL变异 | Extracting INDEL variants")
        
        indel_vcf = self.config.output_path / f"{self.config.base_name}_indels_raw.vcf.gz"
        
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' SelectVariants"
        cmd += f" -R {self.config.reference}"
        cmd += f" -V {self.config.raw_vcf}"
        cmd += f" -select-type INDEL"
        cmd += f" -O {indel_vcf}"
        
        success = self.cmd_runner.run(cmd, "提取INDEL | Extract INDELs")
        
        if success:
            self.config.indel_raw_vcf = str(indel_vcf)
        
        return success
    
    def filter_snps(self):
        """过滤SNP | Filter SNPs"""
        self.logger.info("🧹 过滤SNP变异 | Filtering SNP variants")
        
        snp_filtered = self.config.output_path / f"{self.config.base_name}_snps_filtered.vcf.gz"
        
        # 构建过滤表达式 | Build filter expression
        filters = [
            f'"QD < {self.config.snp_qd}" --filter-name "QD_filter"',
            f'"FS > {self.config.snp_fs}" --filter-name "FS_filter"',
            f'"MQ < {self.config.snp_mq}" --filter-name "MQ_filter"',
            f'"MQRankSum < {self.config.snp_mqrs}" --filter-name "MQRS_filter"',
            f'"ReadPosRankSum < {self.config.snp_rprs}" --filter-name "RPRS_filter"',
            f'"SOR > {self.config.snp_sor}" --filter-name "SOR_filter"'
        ]
        
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' VariantFiltration"
        cmd += f" -R {self.config.reference}"
        cmd += f" -V {self.config.snp_raw_vcf}"
        cmd += f" -O {snp_filtered}"
        
        for filt in filters:
            cmd += f" --filter-expression {filt}"
        
        success = self.cmd_runner.run(cmd, "过滤SNP | Filter SNPs")
        
        if success:
            self.config.snp_filtered_vcf = str(snp_filtered)
        
        return success
    
    def filter_indels(self):
        """过滤INDEL | Filter INDELs"""
        self.logger.info("🧹 过滤INDEL变异 | Filtering INDEL variants")
        
        indel_filtered = self.config.output_path / f"{self.config.base_name}_indels_filtered.vcf.gz"
        
        # 构建过滤表达式 | Build filter expression
        filters = [
            f'"QD < {self.config.indel_qd}" --filter-name "QD_filter"',
            f'"FS > {self.config.indel_fs}" --filter-name "FS_filter"',
            f'"ReadPosRankSum < {self.config.indel_rprs}" --filter-name "RPRS_filter"',
            f'"SOR > {self.config.indel_sor}" --filter-name "SOR_filter"'
        ]
        
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' VariantFiltration"
        cmd += f" -R {self.config.reference}"
        cmd += f" -V {self.config.indel_raw_vcf}"
        cmd += f" -O {indel_filtered}"
        
        for filt in filters:
            cmd += f" --filter-expression {filt}"
        
        success = self.cmd_runner.run(cmd, "过滤INDEL | Filter INDELs")
        
        if success:
            self.config.indel_filtered_vcf = str(indel_filtered)
        
        return success
    
    def merge_variants(self):
        """合并SNP和INDEL | Merge SNPs and INDELs"""
        self.logger.info("🔗 合并SNP和INDEL | Merging SNPs and INDELs")
        
        merged_vcf = self.config.output_path / f"{self.config.base_name}_merged_filtered.vcf.gz"
        
        # 使用bcftools合并 | Use bcftools to merge
        cmd = f"{self.config.bcftools_path} concat"
        cmd += f" -a -O z"
        cmd += f" -o {merged_vcf}"
        cmd += f" {self.config.snp_filtered_vcf} {self.config.indel_filtered_vcf}"
        
        success = self.cmd_runner.run(cmd, "合并变异 | Merge variants")
        
        if success:
            self.config.merged_vcf = str(merged_vcf)
            # 建立索引 | Create index
            self._index_vcf(merged_vcf)
        
        return success
    
    def _index_vcf(self, vcf_file):
        """为VCF文件建立索引 | Index VCF file"""
        cmd = f"{self.config.bcftools_path} index -t {vcf_file}"
        self.cmd_runner.run(cmd, f"索引VCF | Index VCF: {vcf_file}")
