# """
# PLINK GWAS群体分析模块 | PLINK GWAS Population Analysis Module
# """

# import pandas as pd
# from pathlib import Path
# from .utils import CommandRunner

# class PopulationAnalyzer:
#     """群体分析器 | Population Analyzer"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def population_stratification(self) -> bool:
#         """群体分层控制 | Population stratification control"""
#         self.logger.info("群体分层控制 | Population stratification control...")
        
#         # LD剪枝 | LD pruning
#         self.logger.info("LD剪枝 | LD pruning...")
#         cmd = [
#             "plink",
#             "--bfile", "data_qc1",
#             "--indep-pairwise", str(self.config.ld_window_size), 
#             str(self.config.ld_step_size), str(self.config.ld_r2_threshold),
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "pruned_snps"
#         ]
#         self.cmd_runner.run(cmd, "LD剪枝 | LD pruning")
        
#         # 主成分分析 | Principal component analysis
#         self.logger.info("主成分分析 | Principal component analysis...")
#         cmd = [
#             "plink",
#             "--bfile", "data_qc1",
#             "--extract", "pruned_snps.prune.in",
#             "--pca", str(self.config.pca_components),
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "pca_results"
#         ]
#         self.cmd_runner.run(cmd, "主成分分析 | Principal component analysis")
        
#         # 为主成分文件添加header | Add header to PCA file
#         pca_file = Path("pca_results.eigenvec")
#         if pca_file.exists():
#             self.logger.info("为主成分文件添加header | Adding header to PCA file...")
#             header = ["FID", "IID"] + [f"PC{i+1}" for i in range(self.config.pca_components)]
            
#             # 读取原始文件 | Read original file
#             pca_data = pd.read_csv(pca_file, sep=r"\s+", header=None)
#             pca_data.columns = header[:pca_data.shape[1]]
            
#             # 保存带header的文件 | Save file with header
#             pca_data.to_csv("pca_with_header.txt", sep=' ', index=False)
#             self.logger.info("主成分分析完成 | Principal component analysis completed")
#             return True
#         else:
#             self.logger.warning("主成分分析失败 | Principal component analysis failed")
#             return False

# 20250726添加显性模型和隐性模型的配置选项 | Added options for dominant and recessive models
"""
PLINK GWAS群体分析模块 | PLINK GWAS Population Analysis Module
"""

import pandas as pd
import os
from pathlib import Path
from .utils import CommandRunner

class PopulationAnalyzer:
    """群体分析器 | Population Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def population_analysis(self) -> bool:
        """群体结构分析 | Population structure analysis"""
        self.logger.info("开始群体结构分析 | Starting population structure analysis...")
        
        # 检查输入文件 | Check input files
        if not os.path.exists("data_qc1.bed"):
            raise FileNotFoundError("质控后数据文件不存在 | QC data files do not exist: data_qc1.*")
        
        # LD剪枝 | LD pruning
        self.logger.info("进行LD剪枝 | Performing LD pruning...")
        cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--indep-pairwise", str(self.config.ld_window_size), 
            str(self.config.ld_step_size), str(self.config.ld_r2_threshold),
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "pruned_snps"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(cmd, "LD剪枝 | LD pruning")
        
        # 检查剪枝结果 | Check pruning results
        if not os.path.exists("pruned_snps.prune.in"):
            self.logger.warning("LD剪枝失败，跳过主成分分析 | LD pruning failed, skipping PCA")
            return False
        
        # 统计剪枝结果 | Count pruning results
        try:
            with open("pruned_snps.prune.in", 'r') as f:
                pruned_snps = len(f.readlines())
            self.logger.info(f"LD剪枝后SNP数 | SNPs after LD pruning: {pruned_snps}")
        except:
            self.logger.warning("无法统计剪枝后SNP数 | Cannot count pruned SNPs")
        
        # 主成分分析 | Principal component analysis
        self.logger.info("进行主成分分析 | Performing principal component analysis...")
        cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--extract", "pruned_snps.prune.in",
            "--pca", str(self.config.pca_components),
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "pca_results"
        ]
        
        if self.config.threads > 1:
            cmd.extend(["--threads", str(self.config.threads)])
        
        result = self.cmd_runner.run(cmd, "主成分分析 | Principal component analysis", check=False)
        
        # 处理主成分分析结果 | Process PCA results
        pca_file = Path("pca_results.eigenvec")
        if pca_file.exists() and result.returncode == 0:
            self.logger.info("为主成分文件添加header | Adding header to PCA file...")
            header = ["FID", "IID"] + [f"PC{i+1}" for i in range(self.config.pca_components)]
            
            try:
                # 读取原始文件 | Read original file
                pca_data = pd.read_csv(pca_file, sep=r"\s+", header=None)
                
                # 确保列数匹配 | Ensure column count matches
                actual_columns = min(len(header), pca_data.shape[1])
                pca_data = pca_data.iloc[:, :actual_columns]
                pca_data.columns = header[:actual_columns]
                
                # 保存带header的文件 | Save file with header
                pca_data.to_csv("pca_with_header.txt", sep='\t', index=False)
                
                self.logger.info(f"主成分分析完成 | Principal component analysis completed")
                self.logger.info(f"生成主成分数 | Generated PCs: {actual_columns - 2}")
                self.logger.info(f"样本数 | Number of samples: {len(pca_data)}")
                
                return True
                
            except Exception as e:
                self.logger.error(f"处理主成分分析结果失败 | Failed to process PCA results: {e}")
                return False
        else:
            self.logger.warning("主成分分析失败，将不使用PCA校正 | Principal component analysis failed, will not use PCA correction")
            return False