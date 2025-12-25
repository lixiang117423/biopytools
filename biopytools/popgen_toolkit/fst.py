# """
# Fst计算模块 | Fst Calculation Module
# """

# from typing import Dict, List
# from .utils import CommandRunner

# class FstCalculator:
#     """Fst计算器 | Fst Calculator"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def calculate_fst(self, vcf_file: str, group_dict: Dict[str, str]) -> List[str]:
#         """计算群体间Fst | Calculate Fst between populations"""
#         if not group_dict:
#             self.logger.warning("没有分组信息，跳过Fst计算 | No group information, skipping Fst calculation")
#             return []
        
#         self.logger.info("计算群体间Fst | Calculating Fst between populations")
        
#         # 获取所有群体 | Get all populations
#         groups = list(set(group_dict.values()))
#         fst_results = []
        
#         # 为每个群体创建样本列表文件 | Create sample list files for each population
#         group_files = {}
#         for group in groups:
#             samples = [sample for sample, g in group_dict.items() if g == group]
#             group_file = self.config.output_path / f"samples_{group}.txt"
            
#             with open(group_file, 'w') as f:
#                 for sample in samples:
#                     f.write(f"{sample}\n")
            
#             group_files[group] = str(group_file)
#             self.logger.info(f"群体 {group}: {len(samples)} 个样本 | Population {group}: {len(samples)} samples")
        
#         # 计算每对群体间的Fst | Calculate Fst for each pair of populations
#         for i, group1 in enumerate(groups):
#             for j, group2 in enumerate(groups[i+1:], i+1):
#                 self.logger.info(f"计算 {group1} vs {group2} 的Fst | Calculating Fst for {group1} vs {group2}")
                
#                 output_prefix = self.config.output_path / f"fst_{group1}_vs_{group2}"
                
#                 cmd = (
#                     f"{self.config.vcftools_path} --gzvcf {vcf_file} "
#                     f"--weir-fst-pop {group_files[group1]} "
#                     f"--weir-fst-pop {group_files[group2]} "
#                     f"--fst-window-size {self.config.display_window_size} "
#                     f"--fst-window-step {int(self.config.display_window_size * (1 - self.config.window_overlap))} "
#                     f"--out {output_prefix}"
#                 )
                
#                 success = self.cmd_runner.run(cmd, f"Fst计算: {group1} vs {group2}")
                
#                 if success:
#                     fst_results.append(str(output_prefix))
        
#         return fst_results

"""
Fst计算模块 | Fst Calculation Module
"""

from typing import Dict, List
from pathlib import Path
from .utils import CommandRunner
import itertools # 导入itertools以更优雅地处理配对

class FstCalculator:
    """Fst计算器 | Fst Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_fst(self, vcf_file: str, group_dict: Dict[str, str]) -> List[str]:
        """计算群体间Fst | Calculate Fst between populations"""
        if not group_dict or len(set(group_dict.values())) < 2:
            self.logger.warning("需要至少两个群体才能计算Fst，跳过此步骤 | At least two populations are required for Fst calculation, skipping.")
            return []
        
        self.logger.info("计算群体间Fst | Calculating Fst between populations")
        
        # 获取所有群体 | Get all populations
        groups = sorted(list(set(group_dict.values()))) #排序以保证输出顺序一致
        fst_results = []
        
        # 为每个群体创建样本列表文件 | Create sample list files for each population
        # group_files 将存储 {群体名: Path对象}
        group_files = {}
        for group in groups:
            samples = [sample for sample, g in group_dict.items() if g == group]
            # 存储完整的Path对象，方便后续获取文件名和删除文件
            group_file_path = self.config.output_path / f"samples_{group}.txt"
            
            with open(group_file_path, 'w') as f:
                f.write("\n".join(samples))
            
            group_files[group] = group_file_path
            self.logger.info(f"群体 {group}: {len(samples)} 个样本 | Population {group}: {len(samples)} samples")
        
        # 计算每对群体间的Fst | Calculate Fst for each pair of populations
        for group1, group2 in itertools.combinations(groups, 2):
            self.logger.info(f"计算 {group1} vs {group2} 的Fst | Calculating Fst for {group1} vs {group2}")
            
            # --- 核心修改开始 ---
            
            # a. 获取不含路径的文件名 (basename)
            pop1_file_basename = group_files[group1].name
            pop2_file_basename = group_files[group2].name
            output_prefix_basename = f"fst_{group1}_vs_{group2}"
            
            # b. 构建命令时，只使用 basename
            cmd = (
                f"{self.config.vcftools_path} --gzvcf {vcf_file} "
                f"--weir-fst-pop {pop1_file_basename} "     # <-- 修改点
                f"--weir-fst-pop {pop2_file_basename} "     # <-- 修改点
                f"--fst-window-size {self.config.display_window_size} "
                f"--fst-window-step {int(self.config.display_window_size * (1 - self.config.window_overlap))} "
                f"--out {output_prefix_basename}"          # <-- 修改点
            )
            
            success = self.cmd_runner.run(cmd, f"Fst计算: {group1} vs {group2}")
            
            if success:
                # c. 将完整的路径前缀添加到结果列表，供后续模块使用
                full_output_prefix = self.config.output_path / output_prefix_basename
                fst_results.append(str(full_output_prefix))
                
            # --- 核心修改结束 ---

        # (推荐) 清理临时的群体文件
        self.logger.info("清理临时的群体样本文件 | Cleaning up temporary population sample files")
        try:
            for file_path in group_files.values():
                if file_path.exists():
                    file_path.unlink()
        except Exception as e:
            self.logger.warning(f"清理临时文件时出错 | Error during temporary file cleanup: {e}")

        return fst_results