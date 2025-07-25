"""
Fst计算模块 | Fst Calculation Module
"""

from typing import Dict, List
from .utils import CommandRunner

class FstCalculator:
    """Fst计算器 | Fst Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_fst(self, vcf_file: str, group_dict: Dict[str, str]) -> List[str]:
        """计算群体间Fst | Calculate Fst between populations"""
        if not group_dict:
            self.logger.warning("没有分组信息，跳过Fst计算 | No group information, skipping Fst calculation")
            return []
        
        self.logger.info("计算群体间Fst | Calculating Fst between populations")
        
        # 获取所有群体 | Get all populations
        groups = list(set(group_dict.values()))
        fst_results = []
        
        # 为每个群体创建样本列表文件 | Create sample list files for each population
        group_files = {}
        for group in groups:
            samples = [sample for sample, g in group_dict.items() if g == group]
            group_file = self.config.output_path / f"samples_{group}.txt"
            
            with open(group_file, 'w') as f:
                for sample in samples:
                    f.write(f"{sample}\n")
            
            group_files[group] = str(group_file)
            self.logger.info(f"群体 {group}: {len(samples)} 个样本 | Population {group}: {len(samples)} samples")
        
        # 计算每对群体间的Fst | Calculate Fst for each pair of populations
        for i, group1 in enumerate(groups):
            for j, group2 in enumerate(groups[i+1:], i+1):
                self.logger.info(f"计算 {group1} vs {group2} 的Fst | Calculating Fst for {group1} vs {group2}")
                
                output_prefix = self.config.output_path / f"fst_{group1}_vs_{group2}"
                
                cmd = (
                    f"{self.config.vcftools_path} --gzvcf {vcf_file} "
                    f"--weir-fst-pop {group_files[group1]} "
                    f"--weir-fst-pop {group_files[group2]} "
                    f"--fst-window-size {self.config.display_window_size} "
                    f"--fst-window-step {int(self.config.display_window_size * (1 - self.config.window_overlap))} "
                    f"--out {output_prefix}"
                )
                
                success = self.cmd_runner.run(cmd, f"Fst计算: {group1} vs {group2}")
                
                if success:
                    fst_results.append(str(output_prefix))
        
        return fst_results
