"""
🔬 变异检测模块（GVCF模式）| Variant Calling Module (GVCF Mode)
"""

from pathlib import Path
from typing import List
from .utils import CommandRunner, check_file_exists

class VariantCaller:
    """变异检测器 | Variant Caller"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def call_variants_gvcf_mode(self, samples: List[dict], bam_files: List[Path]):
        """GVCF模式变异检测 | GVCF mode variant calling"""
        self.logger.info("=" * 80)
        self.logger.info("🔬 GVCF模式变异检测 | GVCF mode variant calling")
        self.logger.info("=" * 80)
        
        # Step 1: 生成GVCF | Generate GVCF for each sample
        gvcf_files = []
        for sample, bam_file in zip(samples, bam_files):
            gvcf_file = self._haplotype_caller(sample['name'], bam_file)
            gvcf_files.append(gvcf_file)
        
        # Step 2: 合并GVCF | Combine GVCFs
        combined_gvcf = self._combine_gvcfs(gvcf_files)
        
        # Step 3: 联合分型 | Joint genotyping
        raw_vcf = self._genotype_gvcfs(combined_gvcf)
        
        return raw_vcf
    
    def _haplotype_caller(self, sample_name: str, bam_file: Path) -> Path:
        """HaplotypeCaller生成GVCF | HaplotypeCaller to generate GVCF"""
        self.logger.info(f"📊 HaplotypeCaller: {sample_name}")
        
        gvcf_file = self.config.gvcf_dir / f"{sample_name}.g.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(gvcf_file, self.logger):
            self.logger.info("⏩ 跳过GVCF生成（文件已存在）| Skipping GVCF generation (file exists)")
            return gvcf_file
        
        # 构建命令 | Build command
        cmd = (f"{self.config.gatk_path} HaplotypeCaller "
               f"-R {self.config.reference} "
               f"-I {bam_file} "
               f"-O {gvcf_file} "
               f"-ERC GVCF "
               f"--sample-ploidy {self.config.ploidy} "
               f"--native-pair-hmm-threads {self.config.threads}")
        
        # 添加区间限定 | Add intervals if specified
        if self.config.intervals:
            cmd += f" -L {self.config.intervals}"
        
        self.cmd_runner.run(cmd, f"生成GVCF | Generate GVCF: {sample_name}")
        return gvcf_file
    
    def _combine_gvcfs(self, gvcf_files: List[Path]) -> Path:
        """合并GVCF文件 | Combine GVCF files"""
        self.logger.info("🔗 合并GVCF文件 | Combining GVCF files")
        
        combined_gvcf = self.config.gvcf_dir / "combined.g.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(combined_gvcf, self.logger):
            self.logger.info("⏩ 跳过GVCF合并（文件已存在）| Skipping GVCF combining (file exists)")
            return combined_gvcf
        
        # 构建输入参数 | Build input arguments
        variant_inputs = " ".join([f"-V {gvcf}" for gvcf in gvcf_files])
        
        cmd = (f"{self.config.gatk_path} CombineGVCFs "
               f"-R {self.config.reference} "
               f"{variant_inputs} "
               f"-O {combined_gvcf}")
        
        self.cmd_runner.run(cmd, "合并GVCF | Combine GVCFs")
        return combined_gvcf
    
    def _genotype_gvcfs(self, combined_gvcf: Path) -> Path:
        """联合分型 | Joint genotyping"""
        self.logger.info("🧬 联合分型 | Joint genotyping")
        
        raw_vcf = self.config.vcf_dir / "raw_variants.vcf.gz"
        
        if not self.config.force_restart and check_file_exists(raw_vcf, self.logger):
            self.logger.info("⏩ 跳过联合分型（文件已存在）| Skipping joint genotyping (file exists)")
            return raw_vcf
        
        cmd = (f"{self.config.gatk_path} GenotypeGVCFs "
               f"-R {self.config.reference} "
               f"-V {combined_gvcf} "
               f"-O {raw_vcf} "
               f"--sample-ploidy {self.config.ploidy}")
        
        # 添加区间限定 | Add intervals if specified
        if self.config.intervals:
            cmd += f" -L {self.config.intervals}"
        
        self.cmd_runner.run(cmd, "联合分型 | Joint genotyping")
        return raw_vcf
