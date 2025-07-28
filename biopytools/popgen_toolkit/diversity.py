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

"""
多样性参数计算模块 | Diversity Parameters Calculation Module
"""

import os
from pathlib import Path
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
        
        # 使用绝对路径 | Use absolute path
        output_prefix = self.config.output_path / f"diversity_window_{window_size}"
        output_prefix_abs = output_prefix.resolve()
        
        # 检查VCF文件是否存在 | Check if VCF file exists
        if not os.path.exists(vcf_file):
            self.logger.error(f"VCF文件不存在 | VCF file does not exist: {vcf_file}")
            return str(output_prefix_abs)
        
        # 确定VCF文件格式并选择正确的参数 | Determine VCF format and choose correct parameter
        vcf_param = "--gzvcf" if vcf_file.endswith('.gz') else "--vcf"
        
        self.logger.info(f"输入VCF文件: {vcf_file} | Input VCF file: {vcf_file}")
        self.logger.info(f"输出前缀: {output_prefix_abs} | Output prefix: {output_prefix_abs}")
        self.logger.info(f"VCF参数: {vcf_param} | VCF parameter: {vcf_param}")
        
        # 计算窗口步长 | Calculate window step
        window_step = int(window_size * (1 - self.config.window_overlap))
        
        # 构建π值计算命令 | Build π calculation command
        cmd = (
            f"{self.config.vcftools_path} {vcf_param} {vcf_file} "
            f"--window-pi {window_size} "
            f"--window-pi-step {window_step} "
            f"--out {output_prefix_abs}"
        )
        
        success = self.cmd_runner.run(cmd, f"计算π值 (窗口: {window_size}) | Calculating π values (window: {window_size})")
        
        if success:
            # 构建Tajima's D计算命令 | Build Tajima's D calculation command
            cmd_tajima = (
                f"{self.config.vcftools_path} {vcf_param} {vcf_file} "
                f"--TajimaD {window_size} "
                f"--out {output_prefix_abs}"
            )
            
            tajima_success = self.cmd_runner.run(cmd_tajima, f"计算Tajima's D (窗口: {window_size}) | Calculating Tajima's D (window: {window_size})")
            
            if tajima_success:
                self.logger.info(f"多样性参数计算完成 | Diversity parameters calculation completed: {output_prefix_abs}")
            else:
                self.logger.warning(f"Tajima's D计算失败 | Tajima's D calculation failed")
        else:
            self.logger.error(f"π值计算失败 | π value calculation failed")
        
        return str(output_prefix_abs)
    
    def check_output_files(self, output_prefix: str) -> dict:
        """检查输出文件 | Check output files"""
        self.logger.info("检查多样性分析输出文件 | Checking diversity analysis output files")
        
        output_files = {
            'pi': f"{output_prefix}.windowed.pi",
            'tajima_d': f"{output_prefix}.Tajima.D"
        }
        
        results = {}
        for analysis_type, file_path in output_files.items():
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path)
                self.logger.info(f"✓ {analysis_type}文件存在: {file_path} ({file_size} bytes)")
                results[analysis_type] = file_path
            else:
                self.logger.warning(f"✗ {analysis_type}文件不存在: {file_path}")
        
        return results