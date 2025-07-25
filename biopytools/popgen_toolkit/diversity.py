"""
多样性参数计算模块 | Diversity Parameters Calculation Module
"""

from .utils import CommandRunner

class DiversityCalculator:
    """多样性参数计算器 | Diversity Parameters Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_pi_theta_w_tajima_d(self, vcf_file: str, window_size: int) -> str:
        """计算π、θw和Tajima's D | Calculate π, θw and Tajima's D"""
        self.logger.info(f"计算多样性参数 (窗口大小: {window_size} bp) | Calculating diversity parameters (window size: {window_size} bp)")
        
        output_prefix = self.config.output_path / f"diversity_window_{window_size}"
        
        cmd = (
            f"{self.config.vcftools_path} --gzvcf {vcf_file} "
            f"--window-pi {window_size} "
            f"--window-pi-step {int(window_size * (1 - self.config.window_overlap))} "
            f"--out {output_prefix}"
        )
        
        success = self.cmd_runner.run(cmd, f"计算π值 (窗口: {window_size}) | Calculating π values (window: {window_size})")
        
        if success:
            cmd_tajima = (
                f"{self.config.vcftools_path} --gzvcf {vcf_file} "
                f"--TajimaD {window_size} "
                f"--out {output_prefix}"
            )
            
            self.cmd_runner.run(cmd_tajima, f"计算Tajima's D (窗口: {window_size}) | Calculating Tajima's D (window: {window_size})")
        
        return str(output_prefix)
