# """
# LD计算模块 | LD Calculation Module
# """

# from .utils import CommandRunner

# class LDCalculator:
#     """LD计算器 | LD Calculator"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def calculate_ld(self, vcf_file: str) -> str:
#         """计算连锁不平衡 | Calculate Linkage Disequilibrium"""
#         self.logger.info("计算连锁不平衡 (LD) | Calculating Linkage Disequilibrium (LD)")
        
#         output_prefix = self.config.output_path / "ld_analysis"
        
#         # 使用VCFtools计算LD
#         cmd = (
#             f"{self.config.vcftools_path} --gzvcf {vcf_file} "
#             f"--geno-r2 --ld-window-bp 50000 "
#             f"--out {output_prefix}"
#         )
        
#         success = self.cmd_runner.run(cmd, "计算连锁不平衡 | Calculating linkage disequilibrium")
        
#         if success:
#             # 也计算基于SNP数量的LD
#             cmd_snp = (
#                 f"{self.config.vcftools_path} --gzvcf {vcf_file} "
#                 f"--geno-r2 --ld-window 100 "
#                 f"--out {output_prefix}_snp_based"
#             )
            
#             self.cmd_runner.run(cmd_snp, "计算基于SNP数量的LD | Calculating SNP-based LD")
        
#         return str(output_prefix)

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
        
        # --- 核心修改开始 ---

        # 1. 定义不含路径的文件名前缀 (basename)
        output_prefix_basename_bp = "ld_analysis_bp_based"
        # 定义包含完整路径的前缀，用于函数返回
        full_output_prefix_bp = self.config.output_path / output_prefix_basename_bp
        
        # 2. 构建并执行基于物理距离的LD计算命令
        cmd_bp = (
            f"{self.config.vcftools_path} --gzvcf {vcf_file} "
            f"--geno-r2 --ld-window-bp 50000 "
            f"--out {output_prefix_basename_bp}"  # <-- 关键修改：使用不含路径的前缀
        )
        
        success = self.cmd_runner.run(cmd_bp, "计算基于物理距离的LD | Calculating bp-based linkage disequilibrium")
        
        # 仅在第一次成功后才进行第二次计算
        if success:
            # 3. 为基于SNP数量的计算定义其独立的basename
            output_prefix_basename_snp = "ld_analysis_snp_based"
            
            cmd_snp = (
                f"{self.config.vcftools_path} --gzvcf {vcf_file} "
                f"--geno-r2 --ld-window 100 "
                f"--out {output_prefix_basename_snp}" # <-- 关键修改：使用不含路径的前缀
            )
            
            self.cmd_runner.run(cmd_snp, "计算基于SNP数量的LD | Calculating SNP-based LD")

        # --- 核心修改结束 ---
        
        # 按照原设计，返回第一个分析的完整路径前缀
        # 注意：这个函数生成了两种结果，但只返回一个前缀。
        # 这可能需要根据后续ResultsProcessor的处理逻辑来决定是否需要修改。
        # 但为了解决当前路径问题，我们保持原有的返回逻辑。
        return str(full_output_prefix_bp)