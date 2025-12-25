"""
🎯 变异过滤模块（硬过滤+软过滤）| Variant Filtering Module (Hard + Soft Filtering)
"""

from pathlib import Path
from .utils import CommandRunner, check_file_exists

class VariantFilter:
    """变异过滤器 | Variant Filter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def filter_variants(self, raw_vcf: Path):
        """过滤变异 | Filter variants"""
        self.logger.info("=" * 80)
        self.logger.info("🎯 变异过滤 | Variant filtering")
        self.logger.info("=" * 80)
        
        # 分离SNP和InDel | Separate SNPs and InDels
        snp_vcf = self._select_variants(raw_vcf, "SNP")
        indel_vcf = self._select_variants(raw_vcf, "INDEL")
        
        # 硬过滤 | Hard filtering
        self.logger.info("🔨 执行硬过滤 | Performing hard filtering")
        filtered_snp_hard = self._hard_filter(snp_vcf, "SNP")
        filtered_indel_hard = self._hard_filter(indel_vcf, "INDEL")
        
        # 软过滤（标记而不删除）| Soft filtering (mark but not remove)
        self.logger.info("🏷️  执行软过滤 | Performing soft filtering")
        filtered_snp_soft = self._soft_filter(snp_vcf, "SNP")
        filtered_indel_soft = self._soft_filter(indel_vcf, "INDEL")
        
        # 建立索引 | Build indices
        self._index_vcf(filtered_snp_hard)
        self._index_vcf(filtered_indel_hard)
        self._index_vcf(filtered_snp_soft)
        self._index_vcf(filtered_indel_soft)
        
        return {
            'snp_hard': filtered_snp_hard,
            'indel_hard': filtered_indel_hard,
            'snp_soft': filtered_snp_soft,
            'indel_soft': filtered_indel_soft
        }
    
    def _select_variants(self, vcf_file: Path, variant_type: str) -> Path:
        """选择特定类型的变异 | Select specific type of variants"""
        self.logger.info(f"📌 分离{variant_type} | Separating {variant_type}")
        
        output_vcf = self.config.vcf_dir / f"raw_{variant_type.lower()}.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(output_vcf, self.logger):
            return output_vcf
        
        cmd = (f"{self.config.gatk_path} SelectVariants "
               f"-R {self.config.reference} "
               f"-V {vcf_file} "
               f"-O {output_vcf} "
               f"--select-type-to-include {variant_type}")
        
        self.cmd_runner.run(cmd, f"分离{variant_type} | Select {variant_type}")
        return output_vcf
    
    def _hard_filter(self, vcf_file: Path, variant_type: str) -> Path:
        """硬过滤 | Hard filtering"""
        output_vcf = self.config.vcf_dir / f"filtered.hard.{variant_type}.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(output_vcf, self.logger):
            return output_vcf
        
        filter_expr = self.config.get_filter_expression(variant_type)
        
        cmd = (f"{self.config.gatk_path} VariantFiltration "
               f"-R {self.config.reference} "
               f"-V {vcf_file} "
               f"-O {output_vcf} "
               f"--filter-expression \"{filter_expr}\" "
               f"--filter-name \"GATK_HARD_FILTER\"")
        
        self.cmd_runner.run(cmd, f"硬过滤{variant_type} | Hard filter {variant_type}")
        
        # 移除过滤掉的位点 | Remove filtered sites
        passed_vcf = self.config.vcf_dir / f"all_samples.filtered.{variant_type}.vcf.gz"
        
        cmd_exclude = (f"{self.config.gatk_path} SelectVariants "
                      f"-R {self.config.reference} "
                      f"-V {output_vcf} "
                      f"-O {passed_vcf} "
                      f"--exclude-filtered")
        
        self.cmd_runner.run(cmd_exclude, f"移除过滤位点 | Remove filtered sites")
        return passed_vcf
    
    def _soft_filter(self, vcf_file: Path, variant_type: str) -> Path:
        """软过滤（仅标记）| Soft filtering (mark only)"""
        output_vcf = self.config.vcf_dir / f"filtered.soft.{variant_type}.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(output_vcf, self.logger):
            return output_vcf
        
        filter_expr = self.config.get_filter_expression(variant_type)
        
        cmd = (f"{self.config.gatk_path} VariantFiltration "
               f"-R {self.config.reference} "
               f"-V {vcf_file} "
               f"-O {output_vcf} "
               f"--filter-expression \"{filter_expr}\" "
               f"--filter-name \"GATK_SOFT_FILTER\"")
        
        self.cmd_runner.run(cmd, f"软过滤{variant_type} | Soft filter {variant_type}")
        return output_vcf
    
    def _index_vcf(self, vcf_file: Path):
        """建立VCF索引 | Index VCF"""
        tbi_file = Path(str(vcf_file) + '.tbi')
        
        if not self.config.force_restart and check_file_exists(tbi_file, self.logger):
            return
        
        cmd = f"{self.config.gatk_path} IndexFeatureFile -I {vcf_file}"
        self.cmd_runner.run(cmd, f"建立VCF索引 | Index VCF: {vcf_file.name}")
