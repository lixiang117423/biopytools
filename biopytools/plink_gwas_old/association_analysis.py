"""
PLINK GWAS关联分析模块 | PLINK GWAS Association Analysis Module
"""

from pathlib import Path
from .utils import CommandRunner

class AssociationAnalyzer:
    """关联分析器 | Association Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def association_analysis(self, use_pca: bool = True) -> str:
        """关联分析 | Association analysis"""
        self.logger.info("开始关联分析 | Starting association analysis...")
        
        # 基本分析（无协变量） | Basic analysis (no covariates)
        self.logger.info("基本关联分析（无协变量） | Basic association analysis (no covariates)...")
        cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--logistic",
            "--ci", "0.95",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "gwas_basic"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(cmd, "基本关联分析 | Basic association analysis")
        
        # 主成分校正分析 | PCA-adjusted analysis
        if use_pca and Path("pca_with_header.txt").exists():
            self.logger.info("主成分校正关联分析 | PCA-adjusted association analysis...")
            
            # 选择要使用的主成分 | Select PCs to use
            pc_names = [f"PC{i+1}" for i in range(min(self.config.pca_use, self.config.pca_components))]
            
            cmd = [
                "plink",
                "--bfile", "data_qc1", 
                "--logistic",
                "--covar", "pca_with_header.txt",
                "--covar-name", ",".join(pc_names),
                "--ci", "0.95",
                "--allow-extra-chr",
                "--allow-no-sex",
                "--out", "gwas_adjusted"
            ]
            
            if self.config.threads > 1:
                cmd.extend(["--threads", str(self.config.threads)])
            
            result = self.cmd_runner.run(cmd, "主成分校正关联分析 | PCA-adjusted association analysis", check=False)
            
            if Path("gwas_adjusted.assoc.logistic").exists():
                self.logger.info("主成分校正分析完成 | PCA-adjusted analysis completed")
                return "gwas_adjusted"
            else:
                self.logger.warning("主成分校正分析失败，使用基本分析结果 | PCA-adjusted analysis failed, using basic analysis results")
                return "gwas_basic"
        else:
            self.logger.info("跳过主成分校正分析 | Skipping PCA-adjusted analysis")
            return "gwas_basic"
