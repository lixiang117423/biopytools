"""
🗂️ 变异选择模块 | Variant Selection Module
"""

from .utils import CommandRunner

class VariantSelector:
    """🔍 变异选择器 | Variant Selector"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def select_variants(self) -> bool:
        """🎯 选择指定类型的变异 | Select specified type of variants"""
        if self.config.variant_type == "BOTH":
            # 🔄 如果选择所有变异，直接复制文件 | If selecting all variants, copy file directly
            self.logger.info("🔄 选择所有变异类型，跳过变异选择步骤 | Selecting all variant types, skipping variant selection")
            return self._copy_input_file()
        
        if self.config.use_gatk_select:
            return self._select_with_gatk()
        else:
            return self._select_with_bcftools()
    
    def _copy_input_file(self) -> bool:
        """📋 复制输入文件 | Copy input file"""
        import shutil
        try:
            shutil.copy2(self.config.input_vcf, self.config.selected_vcf)
            self.logger.info(f"✅ 已复制输入文件 | Input file copied: {self.config.selected_vcf}")
            return True
        except Exception as e:
            self.logger.error(f"❌ 复制文件失败 | File copy failed: {e}")
            return False
    
    def _select_with_gatk(self) -> bool:
        """🧬 使用GATK SelectVariants选择变异 | Select variants using GATK SelectVariants"""
        self.logger.info(f"🧬 使用GATK SelectVariants选择{self.config.variant_type}变异 | Using GATK SelectVariants to select {self.config.variant_type} variants")
        
        cmd_parts = [
            self.config.gatk_path,
            f"--java-options \"{self.config.java_options}\"",
            "SelectVariants",
            f"-V {self.config.input_vcf}",
            f"-select-type {self.config.variant_type}",
            f"-O {self.config.selected_vcf}"
        ]
        
        # ⚙️ 添加可选参数 | Add optional parameters
        if self.config.exclude_filtered:
            cmd_parts.append("--exclude-filtered")
            
        if self.config.exclude_non_variants:
            cmd_parts.append("--exclude-non-variants")
            
        if self.config.remove_unused_alternates:
            cmd_parts.append("--remove-unused-alternates")
            
        if self.config.restrict_alleles_to != "ALL":
            cmd_parts.append(f"--restrict-alleles-to {self.config.restrict_alleles_to}")
            
        if self.config.select_random_fraction is not None:
            cmd_parts.append(f"--select-random-fraction {self.config.select_random_fraction}")
            
        if self.config.max_indel_size is not None:
            cmd_parts.append(f"--max-indel-size {self.config.max_indel_size}")
            
        if self.config.min_indel_size is not None:
            cmd_parts.append(f"--min-indel-size {self.config.min_indel_size}")
            
        if self.config.max_filtered_genotypes is not None:
            cmd_parts.append(f"--max-filtered-genotypes {self.config.max_filtered_genotypes}")
            
        if self.config.min_filtered_genotypes is not None:
            cmd_parts.append(f"--min-filtered-genotypes {self.config.min_filtered_genotypes}")
            
        if self.config.max_fraction_filtered_genotypes is not None:
            cmd_parts.append(f"--max-fraction-filtered-genotypes {self.config.max_fraction_filtered_genotypes}")
            
        if self.config.min_fraction_filtered_genotypes is not None:
            cmd_parts.append(f"--min-fraction-filtered-genotypes {self.config.min_fraction_filtered_genotypes}")
            
        if self.config.max_nocall_number is not None:
            cmd_parts.append(f"--max-nocall-number {self.config.max_nocall_number}")
            
        if self.config.max_nocall_fraction is not None:
            cmd_parts.append(f"--max-nocall-fraction {self.config.max_nocall_fraction}")
        
        cmd = " \\\n    ".join(cmd_parts)
        success = self.cmd_runner.run(cmd, f"🧬 GATK {self.config.variant_type}选择 | GATK {self.config.variant_type} selection")
        
        if success:
            self.logger.info(f"✅ {self.config.variant_type}变异选择完成 | {self.config.variant_type} variant selection completed: {self.config.selected_vcf}")
        
        return success
    
    def _select_with_bcftools(self) -> bool:
        """🔧 使用BCFtools选择变异 | Select variants using BCFtools"""
        self.logger.info(f"🔧 使用BCFtools选择{self.config.variant_type}变异 | Using BCFtools to select {self.config.variant_type} variants")
        
        if self.config.variant_type == "SNP":
            # 🧬 选择SNP：排除INDEL | Select SNPs: exclude INDELs
            cmd = (
                f"{self.config.bcftools_path} view "
                f"--exclude-types indels "
                f"-O z -o {self.config.selected_vcf} "
                f"{self.config.input_vcf}"
            )
        elif self.config.variant_type == "INDEL":
            # 🧬 选择INDEL：排除SNP | Select INDELs: exclude SNPs
            cmd = (
                f"{self.config.bcftools_path} view "
                f"--include-types indels "
                f"-O z -o {self.config.selected_vcf} "
                f"{self.config.input_vcf}"
            )
        else:
            self.logger.error(f"❌ 不支持的变异类型 | Unsupported variant type: {self.config.variant_type}")
            return False
        
        success = self.cmd_runner.run(cmd, f"🔧 BCFtools {self.config.variant_type}选择 | BCFtools {self.config.variant_type} selection")
        
        if success:
            self.logger.info(f"✅ {self.config.variant_type}变异选择完成 | {self.config.variant_type} variant selection completed: {self.config.selected_vcf}")
        
        return success
