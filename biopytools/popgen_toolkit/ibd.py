# """
# IBD计算模块 | IBD Calculation Module
# """

# from .utils import CommandRunner

# class IBDCalculator:
#     """IBD计算器 | IBD Calculator"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def calculate_ibd(self, vcf_file: str) -> str:
#         """计算IBD | Calculate IBD"""
#         self.logger.info("计算同源性分析 (IBD) | Calculating Identity-by-Descent (IBD)")
        
#         # 首先转换VCF到PLINK格式 | First convert VCF to PLINK format
#         plink_prefix = self.config.output_path / "ibd_analysis"
        
#         # VCF转PLINK
#         cmd_convert = (
#             f"{self.config.plink_path} --vcf {vcf_file} "
#             f"--make-bed --out {plink_prefix} --allow-extra-chr"
#         )
        
#         success = self.cmd_runner.run(cmd_convert, "VCF转PLINK格式用于IBD分析 | Converting VCF to PLINK format for IBD")
        
#         if success:
#             # 运行IBD分析
#             cmd_ibd = (
#                 f"{self.config.plink_path} --bfile {plink_prefix} "
#                 f"--genome --out {plink_prefix} --allow-extra-chr"
#             )
            
#             self.cmd_runner.run(cmd_ibd, "执行IBD分析 | Performing IBD analysis")
        
#         return str(plink_prefix)

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
        
        # --- 核心修改开始 ---

        # 1. 定义不含路径的文件名前缀 (basename)
        plink_prefix_basename = "ibd_analysis"
        # 定义包含完整路径的前缀，用于函数返回和内部逻辑
        full_plink_prefix = self.config.output_path / plink_prefix_basename
        
        # (额外改进) 检查VCF文件是否压缩，以选择正确的plink参数
        # --vcf 用于 .vcf, --bcf 用于 .vcf.gz (在PLINK 1.9中)
        vcf_param = "--bcf" if vcf_file.endswith('.gz') else "--vcf"

        # 2. VCF转PLINK，在命令中只使用basename
        cmd_convert = (
            f"{self.config.plink_path} {vcf_param} {vcf_file} "
            f"--make-bed --out {plink_prefix_basename} --allow-extra-chr" # <-- 关键修改
        )
        
        success = self.cmd_runner.run(cmd_convert, "VCF转PLINK格式用于IBD分析 | Converting VCF to PLINK format for IBD")
        
        if success:
            # 3. 运行IBD分析，在命令中也只使用basename
            cmd_ibd = (
                f"{self.config.plink_path} --bfile {plink_prefix_basename} " # <-- 关键修改
                f"--genome --out {plink_prefix_basename} --allow-extra-chr"   # <-- 关键修改
            )
            
            self.cmd_runner.run(cmd_ibd, "执行IBD分析 | Performing IBD analysis")
        
        # 4. 返回完整的路径前缀
        return str(full_plink_prefix)
        # --- 核心修改结束 ---