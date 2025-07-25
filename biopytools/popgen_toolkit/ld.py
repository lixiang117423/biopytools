"""
LD计算模块 | LD Calculation Module
"""

from .utils import CommandRunner

class LDCalculator:
    """LD计算器 | LD Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_ld(self, vcf_file: str) -> str:
        """计算连锁不平衡 | Calculate Linkage Disequilibrium"""
        self.logger.info("计算连锁不平衡 (LD) | Calculating Linkage Disequilibrium (LD)")
        
        output_prefix = self.config.output_path / "ld_analysis"
        
        # 使用VCFtools计算LD
        cmd = (
            f"{self.config.vcftools_path} --gzvcf {vcf_file} "
            f"--geno-r2 --ld-window-bp 50000 "
            f"--out {output_prefix}"
        )
        
        success = self.cmd_runner.run(cmd, "计算连锁不平衡 | Calculating linkage disequilibrium")
        
        if success:
            # 也计算基于SNP数量的LD
            cmd_snp = (
                f"{self.config.vcftools_path} --gzvcf {vcf_file} "
                f"--geno-r2 --ld-window 100 "
                f"--out {output_prefix}_snp_based"
            )
            
            self.cmd_runner.run(cmd_snp, "计算基于SNP数量的LD | Calculating SNP-based LD")
        
        return str(output_prefix)
