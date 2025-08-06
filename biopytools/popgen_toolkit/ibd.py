"""
IBD计算模块 | IBD Calculation Module
"""

from .utils import CommandRunner

class IBDCalculator:
    """IBD计算器 | IBD Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_ibd(self, vcf_file: str) -> str:
        """计算IBD | Calculate IBD"""
        self.logger.info("计算同源性分析 (IBD) | Calculating Identity-by-Descent (IBD)")
        
        # 首先转换VCF到PLINK格式 | First convert VCF to PLINK format
        plink_prefix = self.config.output_path / "ibd_analysis"
        
        # VCF转PLINK
        cmd_convert = (
            f"{self.config.plink_path} --vcf {vcf_file} "
            f"--make-bed --out {plink_prefix} --allow-extra-chr"
        )
        
        success = self.cmd_runner.run(cmd_convert, "VCF转PLINK格式用于IBD分析 | Converting VCF to PLINK format for IBD")
        
        if success:
            # 运行IBD分析
            cmd_ibd = (
                f"{self.config.plink_path} --bfile {plink_prefix} "
                f"--genome --out {plink_prefix} --allow-extra-chr"
            )
            
            self.cmd_runner.run(cmd_ibd, "执行IBD分析 | Performing IBD analysis")
        
        return str(plink_prefix)
