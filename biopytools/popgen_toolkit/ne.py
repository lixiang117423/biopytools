"""
有效群体大小计算模块 | Effective Population Size Calculation Module
"""

from typing import Dict, List
from .utils import CommandRunner

class EffectivePopSizeCalculator:
    """有效群体大小计算器 | Effective Population Size Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_ne_smcpp(self, vcf_file: str, group_dict: Dict[str, str]) -> List[str]:
        """使用SMC++计算有效群体大小 | Calculate effective population size using SMC++"""
        if not group_dict:
            self.logger.warning("没有分组信息，跳过SMC++分析 | No group information, skipping SMC++ analysis")
            return []
        
        self.logger.info("使用SMC++计算有效群体大小 | Calculating effective population size using SMC++")
        
        # 获取所有群体 | Get all populations
        groups = list(set(group_dict.values()))
        ne_results = []
        
        for group in groups:
            self.logger.info(f"计算群体 {group} 的有效群体大小 | Calculating Ne for population {group}")
            
            # 获取该群体的样本 | Get samples for this population
            samples = [sample for sample, g in group_dict.items() if g == group]
            
            if len(samples) < 2:
                self.logger.warning(f"群体 {group} 样本数量不足，跳过SMC++分析 | Insufficient samples in population {group}, skipping SMC++")
                continue
            
            # 创建SMC++输入文件 | Create SMC++ input file
            output_prefix = self.config.output_path / f"smc_input_{group}"
            
            # 转换为SMC++格式
            cmd_convert = (
                f"smc++ vcf2smc {vcf_file} {output_prefix}.smc.gz "
                f"{group} {':'.join(samples)}"
            )
            
            success = self.cmd_runner.run(cmd_convert, f"转换VCF为SMC++格式 (群体: {group}) | Converting VCF to SMC++ format (population: {group})")
            
            if success:
                # 运行SMC++估计
                smc_output = self.config.output_path / f"smc_results_{group}"
                cmd_estimate = (
                    f"smc++ estimate 1.25e-8 {output_prefix}.smc.gz "
                    f"--output-dir {smc_output} --cores {self.config.threads}"
                )
                
                success_est = self.cmd_runner.run(cmd_estimate, f"SMC++有效群体大小估计 (群体: {group}) | SMC++ Ne estimation (population: {group})")
                
                if success_est:
                    ne_results.append(str(smc_output))
        
        return ne_results
