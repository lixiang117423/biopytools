"""
PLINK GWAS群体分析模块 | PLINK GWAS Population Analysis Module
"""

import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class PopulationAnalyzer:
    """群体分析器 | Population Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def population_stratification(self) -> bool:
        """群体分层控制 | Population stratification control"""
        self.logger.info("群体分层控制 | Population stratification control...")
        
        # LD剪枝 | LD pruning
        self.logger.info("LD剪枝 | LD pruning...")
        cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--indep-pairwise", str(self.config.ld_window_size), 
            str(self.config.ld_step_size), str(self.config.ld_r2_threshold),
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "pruned_snps"
        ]
        self.cmd_runner.run(cmd, "LD剪枝 | LD pruning")
        
        # 主成分分析 | Principal component analysis
        self.logger.info("主成分分析 | Principal component analysis...")
        cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--extract", "pruned_snps.prune.in",
            "--pca", str(self.config.pca_components),
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "pca_results"
        ]
        self.cmd_runner.run(cmd, "主成分分析 | Principal component analysis")
        
        # 为主成分文件添加header | Add header to PCA file
        pca_file = Path("pca_results.eigenvec")
        if pca_file.exists():
            self.logger.info("为主成分文件添加header | Adding header to PCA file...")
            header = ["FID", "IID"] + [f"PC{i+1}" for i in range(self.config.pca_components)]
            
            # 读取原始文件 | Read original file
            pca_data = pd.read_csv(pca_file, sep=r"\s+", header=None)
            pca_data.columns = header[:pca_data.shape[1]]
            
            # 保存带header的文件 | Save file with header
            pca_data.to_csv("pca_with_header.txt", sep=' ', index=False)
            self.logger.info("主成分分析完成 | Principal component analysis completed")
            return True
        else:
            self.logger.warning("主成分分析失败 | Principal component analysis failed")
            return False
