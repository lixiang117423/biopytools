# """
# PLINK GWAS关联分析模块 | PLINK GWAS Association Analysis Module
# """

# from pathlib import Path
# from .utils import CommandRunner

# class AssociationAnalyzer:
#     """关联分析器 | Association Analyzer"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def association_analysis(self, use_pca: bool = True) -> str:
#         """关联分析 | Association analysis"""
#         self.logger.info("开始关联分析 | Starting association analysis...")
#         self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
        
#         # 根据表型类型选择分析方法 | Choose analysis method based on trait type
#         if self.config.trait_type == "qualitative":
#             return self._qualitative_association_analysis(use_pca)
#         else:
#             return self._quantitative_association_analysis(use_pca)
    
#     def _qualitative_association_analysis(self, use_pca: bool = True) -> str:
#         """质量性状关联分析 | Qualitative trait association analysis"""
#         self.logger.info("质量性状关联分析（逻辑回归） | Qualitative trait association analysis (logistic regression)")
        
#         # 基本分析（无协变量） | Basic analysis (no covariates)
#         self.logger.info("基本关联分析（无协变量） | Basic association analysis (no covariates)...")
#         cmd = [
#             "plink",
#             "--bfile", "data_qc1",
#             "--logistic",
#             "--ci", "0.95",
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "gwas_basic"
#         ]
        
#         if self.config.threads > 1:
#             cmd.extend(["--threads", str(self.config.threads)])
        
#         self.cmd_runner.run(cmd, "基本关联分析 | Basic association analysis")
        
#         # 主成分校正分析 | PCA-adjusted analysis
#         if use_pca and Path("pca_with_header.txt").exists():
#             self.logger.info("主成分校正关联分析 | PCA-adjusted association analysis...")
            
#             # 选择要使用的主成分 | Select PCs to use
#             pc_names = [f"PC{i+1}" for i in range(min(self.config.pca_use, self.config.pca_components))]
            
#             cmd = [
#                 "plink",
#                 "--bfile", "data_qc1", 
#                 "--logistic",
#                 "--covar", "pca_with_header.txt",
#                 "--covar-name", ",".join(pc_names),
#                 "--ci", "0.95",
#                 "--allow-extra-chr",
#                 "--allow-no-sex",
#                 "--out", "gwas_adjusted"
#             ]
            
#             if self.config.threads > 1:
#                 cmd.extend(["--threads", str(self.config.threads)])
            
#             result = self.cmd_runner.run(cmd, "主成分校正关联分析 | PCA-adjusted association analysis", check=False)
            
#             if Path("gwas_adjusted.assoc.logistic").exists():
#                 self.logger.info("主成分校正分析完成 | PCA-adjusted analysis completed")
#                 return "gwas_adjusted"
#             else:
#                 self.logger.warning("主成分校正分析失败，使用基本分析结果 | PCA-adjusted analysis failed, using basic analysis results")
#                 return "gwas_basic"
#         else:
#             self.logger.info("跳过主成分校正分析 | Skipping PCA-adjusted analysis")
#             return "gwas_basic"
    
#     def _quantitative_association_analysis(self, use_pca: bool = True) -> str:
#         """数量性状关联分析 | Quantitative trait association analysis"""
#         self.logger.info("数量性状关联分析（线性回归） | Quantitative trait association analysis (linear regression)")
        
#         # 基本分析（无协变量） | Basic analysis (no covariates)
#         self.logger.info("基本关联分析（无协变量） | Basic association analysis (no covariates)...")
#         cmd = [
#             "plink",
#             "--bfile", "data_qc1",
#             "--linear",
#             "--ci", "0.95",
#             "--allow-extra-chr",
#             "--allow-no-sex",
#             "--out", "gwas_basic"
#         ]
        
#         if self.config.threads > 1:
#             cmd.extend(["--threads", str(self.config.threads)])
        
#         self.cmd_runner.run(cmd, "基本关联分析 | Basic association analysis")
        
#         # 主成分校正分析 | PCA-adjusted analysis
#         if use_pca and Path("pca_with_header.txt").exists():
#             self.logger.info("主成分校正关联分析 | PCA-adjusted association analysis...")
            
#             # 选择要使用的主成分 | Select PCs to use
#             pc_names = [f"PC{i+1}" for i in range(min(self.config.pca_use, self.config.pca_components))]
            
#             cmd = [
#                 "plink",
#                 "--bfile", "data_qc1", 
#                 "--linear",
#                 "--covar", "pca_with_header.txt",
#                 "--covar-name", ",".join(pc_names),
#                 "--ci", "0.95",
#                 "--allow-extra-chr",
#                 "--allow-no-sex",
#                 "--out", "gwas_adjusted"
#             ]
            
#             if self.config.threads > 1:
#                 cmd.extend(["--threads", str(self.config.threads)])
            
#             result = self.cmd_runner.run(cmd, "主成分校正关联分析 | PCA-adjusted association analysis", check=False)
            
#             if Path("gwas_adjusted.assoc.linear").exists():
#                 self.logger.info("主成分校正分析完成 | PCA-adjusted analysis completed")
#                 return "gwas_adjusted"
#             else:
#                 self.logger.warning("主成分校正分析失败，使用基本分析结果 | PCA-adjusted analysis failed, using basic analysis results")
#                 return "gwas_basic"
#         else:
#             self.logger.info("跳过主成分校正分析 | Skipping PCA-adjusted analysis")
#             return "gwas_basic"

# 20250726 添加显性模型和隐性模型的配置选项 | Added options for dominant and recessive models
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
        self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
        self.logger.info(f"遗传模型 | Genetic model: {self.config.genetic_model}")
        
        # 根据表型类型选择分析方法 | Choose analysis method based on trait type
        if self.config.trait_type == "qualitative":
            return self._qualitative_association_analysis(use_pca)
        else:
            return self._quantitative_association_analysis(use_pca)
    
    def _qualitative_association_analysis(self, use_pca: bool = True) -> str:
        """质量性状关联分析 | Qualitative trait association analysis"""
        self.logger.info("质量性状关联分析（逻辑回归） | Qualitative trait association analysis (logistic regression)")
        
        # 根据遗传模型选择分析策略 | Choose analysis strategy based on genetic model
        if self.config.genetic_model == "all":
            return self._run_all_models_analysis(use_pca, analysis_type="logistic")
        else:
            return self._run_single_model_analysis(use_pca, analysis_type="logistic")
    
    def _quantitative_association_analysis(self, use_pca: bool = True) -> str:
        """数量性状关联分析 | Quantitative trait association analysis"""
        self.logger.info("数量性状关联分析（线性回归） | Quantitative trait association analysis (linear regression)")
        
        # 根据遗传模型选择分析策略 | Choose analysis strategy based on genetic model
        if self.config.genetic_model == "all":
            return self._run_all_models_analysis(use_pca, analysis_type="linear")
        else:
            return self._run_single_model_analysis(use_pca, analysis_type="linear")
    
    def _run_all_models_analysis(self, use_pca: bool, analysis_type: str) -> str:
        """运行所有遗传模型分析 | Run all genetic models analysis"""
        self.logger.info("运行所有遗传模型分析 | Running all genetic models analysis...")
        
        # 基本分析（所有模型，无协变量） | Basic analysis (all models, no covariates)
        self.logger.info("基本关联分析（所有模型，无协变量） | Basic association analysis (all models, no covariates)...")
        
        base_cmd = [
            "plink",
            "--bfile", "data_qc1",
            "--model",  # 使用--model参数来测试所有模型
            "--ci", "0.95",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", "gwas_basic_all_models"
        ]
        
        if self.config.threads > 1:
            base_cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(base_cmd, "基本关联分析（所有模型） | Basic association analysis (all models)")
        
        # 主成分校正分析 | PCA-adjusted analysis
        if use_pca and Path("pca_with_header.txt").exists():
            self.logger.info("主成分校正关联分析（所有模型） | PCA-adjusted association analysis (all models)...")
            
            # 选择要使用的主成分 | Select PCs to use
            pc_names = [f"PC{i+1}" for i in range(min(self.config.pca_use, self.config.pca_components))]
            
            pca_cmd = [
                "plink",
                "--bfile", "data_qc1",
                "--model",
                "--covar", "pca_with_header.txt",
                "--covar-name", ",".join(pc_names),
                "--ci", "0.95",
                "--allow-extra-chr",
                "--allow-no-sex",
                "--out", "gwas_adjusted_all_models"
            ]
            
            if self.config.threads > 1:
                pca_cmd.extend(["--threads", str(self.config.threads)])
            
            result = self.cmd_runner.run(pca_cmd, "主成分校正关联分析（所有模型） | PCA-adjusted association analysis (all models)", check=False)
            
            if Path("gwas_adjusted_all_models.model").exists():
                self.logger.info("主成分校正分析（所有模型）完成 | PCA-adjusted analysis (all models) completed")
                return "gwas_adjusted_all_models"
            else:
                self.logger.warning("主成分校正分析失败，使用基本分析结果 | PCA-adjusted analysis failed, using basic analysis results")
                return "gwas_basic_all_models"
        else:
            self.logger.info("跳过主成分校正分析 | Skipping PCA-adjusted analysis")
            return "gwas_basic_all_models"
    
    def _run_single_model_analysis(self, use_pca: bool, analysis_type: str) -> str:
        """运行单一遗传模型分析 | Run single genetic model analysis"""
        self.logger.info(f"运行{self.config.genetic_model}模型分析 | Running {self.config.genetic_model} model analysis...")
        
        # 确定PLINK命令参数 | Determine PLINK command parameters
        if analysis_type == "logistic":
            analysis_flag = "--logistic"
        else:
            analysis_flag = "--linear"
        
        # 添加遗传模型参数 | Add genetic model parameters
        model_params = self._get_model_parameters()
        
        # 基本分析（指定模型，无协变量） | Basic analysis (specified model, no covariates)
        self.logger.info(f"基本关联分析（{self.config.genetic_model}模型，无协变量） | Basic association analysis ({self.config.genetic_model} model, no covariates)...")
        
        base_cmd = [
            "plink",
            "--bfile", "data_qc1",
            analysis_flag,
            "--ci", "0.95",
            "--allow-extra-chr",
            "--allow-no-sex",
            "--out", f"gwas_basic_{self.config.genetic_model}"
        ]
        
        # 添加模型特定参数 | Add model-specific parameters
        base_cmd.extend(model_params)
        
        if self.config.threads > 1:
            base_cmd.extend(["--threads", str(self.config.threads)])
        
        self.cmd_runner.run(base_cmd, f"基本关联分析（{self.config.genetic_model}模型） | Basic association analysis ({self.config.genetic_model} model)")
        
        # 主成分校正分析 | PCA-adjusted analysis
        if use_pca and Path("pca_with_header.txt").exists():
            self.logger.info(f"主成分校正关联分析（{self.config.genetic_model}模型） | PCA-adjusted association analysis ({self.config.genetic_model} model)...")
            
            # 选择要使用的主成分 | Select PCs to use
            pc_names = [f"PC{i+1}" for i in range(min(self.config.pca_use, self.config.pca_components))]
            
            pca_cmd = [
                "plink",
                "--bfile", "data_qc1",
                analysis_flag,
                "--covar", "pca_with_header.txt",
                "--covar-name", ",".join(pc_names),
                "--ci", "0.95",
                "--allow-extra-chr",
                "--allow-no-sex",
                "--out", f"gwas_adjusted_{self.config.genetic_model}"
            ]
            
            # 添加模型特定参数 | Add model-specific parameters
            pca_cmd.extend(model_params)
            
            if self.config.threads > 1:
                pca_cmd.extend(["--threads", str(self.config.threads)])
            
            expected_file = f"gwas_adjusted_{self.config.genetic_model}.assoc.{analysis_type}"
            
            result = self.cmd_runner.run(pca_cmd, f"主成分校正关联分析（{self.config.genetic_model}模型） | PCA-adjusted association analysis ({self.config.genetic_model} model)", check=False)
            
            if Path(expected_file).exists():
                self.logger.info(f"主成分校正分析（{self.config.genetic_model}模型）完成 | PCA-adjusted analysis ({self.config.genetic_model} model) completed")
                return f"gwas_adjusted_{self.config.genetic_model}"
            else:
                self.logger.warning("主成分校正分析失败，使用基本分析结果 | PCA-adjusted analysis failed, using basic analysis results")
                return f"gwas_basic_{self.config.genetic_model}"
        else:
            self.logger.info("跳过主成分校正分析 | Skipping PCA-adjusted analysis")
            return f"gwas_basic_{self.config.genetic_model}"
    
    def _get_model_parameters(self) -> list:
        """获取模型特定参数 | Get model-specific parameters"""
        if self.config.genetic_model == "additive":
            # 加性模型是默认的，不需要额外参数 | Additive model is default, no extra parameters needed
            return []
        elif self.config.genetic_model == "dominant":
            # 显性模型 | Dominant model
            return ["--dominant"]
        elif self.config.genetic_model == "recessive":
            # 隐性模型 | Recessive model  
            return ["--recessive"]
        else:
            return []
    
    def get_result_file_suffix(self) -> str:
        """获取结果文件后缀 | Get result file suffix"""
        if self.config.genetic_model == "all":
            return ".model"
        elif self.config.trait_type == "qualitative":
            return ".assoc.logistic"
        else:
            return ".assoc.linear"