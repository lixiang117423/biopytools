# """
# 多样性参数计算模块 | Diversity Parameters Calculation Module
# """

# from .utils import CommandRunner

# class DiversityCalculator:
#     """多样性参数计算器 | Diversity Parameters Calculator"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def calculate_pi_theta_w_tajima_d(self, vcf_file: str, window_size: int) -> str:
#         """计算π、θw和Tajima's D | Calculate π, θw and Tajima's D"""
#         self.logger.info(f"计算多样性参数 (窗口大小: {window_size} bp) | Calculating diversity parameters (window size: {window_size} bp)")
        
#         output_prefix = self.config.output_path / f"diversity_window_{window_size}"
        
#         cmd = (
#             f"{self.config.vcftools_path} --gzvcf {vcf_file} "
#             f"--window-pi {window_size} "
#             f"--window-pi-step {int(window_size * (1 - self.config.window_overlap))} "
#             f"--out {output_prefix}"
#         )
        
#         success = self.cmd_runner.run(cmd, f"计算π值 (窗口: {window_size}) | Calculating π values (window: {window_size})")
        
#         if success:
#             cmd_tajima = (
#                 f"{self.config.vcftools_path} --gzvcf {vcf_file} "
#                 f"--TajimaD {window_size} "
#                 f"--out {output_prefix}"
#             )
            
#             self.cmd_runner.run(cmd_tajima, f"计算Tajima's D (窗口: {window_size}) | Calculating Tajima's D (window: {window_size})")
        
#         return str(output_prefix)

# diversity.py (完整修改版)

"""
多样性参数计算模块 | Diversity Parameters Calculation Module
"""

from typing import List
from .utils import CommandRunner

class DiversityCalculator:
    """多样性参数计算器 | Diversity Parameters Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def _calculate_pi(self, vcf_file: str, window_size: int, output_prefix_basename: str) -> bool:
        """(私有方法) 计算π值 | (Private method) Calculate π"""
        cmd_pi = (
            f"{self.config.vcftools_path} --gzvcf {vcf_file} "
            f"--window-pi {window_size} "
            f"--window-pi-step {int(window_size * (1 - self.config.window_overlap))} "
            f"--out {output_prefix_basename}"
        )
        return self.cmd_runner.run(cmd_pi, f"计算π值 (窗口: {window_size}) | Calculating π values (window: {window_size})")

    def _calculate_tajima_d(self, vcf_file: str, window_size: int, output_prefix_basename: str) -> bool:
        """(私有方法) 计算Tajima's D | (Private method) Calculate Tajima's D"""
        cmd_tajima = (
            f"{self.config.vcftools_path} --gzvcf {vcf_file} "
            f"--TajimaD {window_size} "
            f"--out {output_prefix_basename}"
        )
        return self.cmd_runner.run(cmd_tajima, f"计算Tajima's D (窗口: {window_size}) | Calculating Tajima's D (window: {window_size})")

    def run_diversity_calculations(self, vcf_file: str) -> List[str]:
        """根据配置运行多样性计算 (π, Tajima's D) | Run diversity calculations based on config (π, Tajima's D)"""
        output_prefixes = []
        
        for window_size in self.config.window_sizes:
            self.logger.info(f"处理多样性参数 (窗口大小: {window_size} bp) | Processing diversity parameters (window size: {window_size} bp)")
            
            output_prefix_basename = f"diversity_window_{window_size}"
            full_output_prefix = self.config.output_path / output_prefix_basename
            
            # 标志位，记录当前窗口是否有任何计算被执行
            calculation_performed = False

            # 根据配置决定是否计算 π
            if self.config.calculate_pi:
                self._calculate_pi(vcf_file, window_size, output_prefix_basename)
                calculation_performed = True

            # 根据配置决定是否计算 Tajima's D
            if self.config.calculate_tajima_d:
                self._calculate_tajima_d(vcf_file, window_size, output_prefix_basename)
                calculation_performed = True

            # 如果为当前窗口执行了任何计算，就将输出前缀添加到结果列表
            if calculation_performed:
                output_prefixes.append(str(full_output_prefix))
                
        return output_prefixes