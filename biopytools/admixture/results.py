"""
 ADMIXTURE结果处理和可视化模块|ADMIXTURE Results Processing and Visualization Module
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

class CovariateGenerator:
    """ 协变量生成器|Covariate Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_gwas_covariates(self, best_k: int):
        """ 生成GWAS协变量文件|Generate GWAS covariate file"""
        self.logger.info(f" 生成GWAS协变量文件|Generating GWAS covariate file (K={best_k})")
        
        #  读取ADMIXTURE结果|Read ADMIXTURE results
        q_file = os.path.join(self.config.output_dir, f"{self.config.base_name}.{best_k}.Q")
        q_data = pd.read_csv(q_file, sep=r'\s+', header=None)
        
        #  读取个体信息|Read individual information
        fam_file = os.path.join(self.config.output_dir, f"{self.config.base_name}.fam")
        fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None)
        
        #  创建协变量文件|Create covariate file
        #  使用前K-1个祖先成分作为协变量（避免共线性）| Use first K-1 ancestry components as covariates
        covariates = pd.DataFrame()
        covariates['FID'] = fam_data.iloc[:, 0]
        covariates['IID'] = fam_data.iloc[:, 1]
        
        #  添加祖先成分协变量|Add ancestry component covariates
        for i in range(best_k - 1):
            covariates[f'PC{i+1}'] = q_data.iloc[:, i]
        
        #  保存协变量文件|Save covariate file
        covar_file = os.path.join(self.config.output_dir, "gwas_covariates.txt")
        covariates.to_csv(covar_file, sep='\t', header=True, index=False) # Use header=True for clarity
        
        self.logger.info(f"GWAS协变量文件已保存|GWAS covariate file saved: {covar_file}")
        self.logger.info("在PLINK中使用|Use in PLINK: --covar gwas_covariates.txt --covar-name PC1,PC2,...")
        
        return covar_file

class PlotGenerator:
    """ 绘图生成器|Plot Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_plots(self, q_data: pd.DataFrame, best_k: int):
        """ 生成可视化图表|Generate visualization plots"""
        self.logger.info("正在生成可视化图表|Generating visualization plots...")
        
        #  生成R脚本|Generate R script
        r_script = self._create_r_script(best_k)
        
        #  保存R脚本|Save R script
        r_file = os.path.join(self.config.output_dir, "generate_plots.R")
        with open(r_file, 'w', encoding='utf-8') as f:
            f.write(r_script)
        
        #  尝试运行R脚本|Try to run R script
        try:
            import subprocess
            from shutil import which
            if which("Rscript") is None:
                self.logger.warning("Rscript未安装或不在PATH中，跳过图表生成|Rscript not installed or not in PATH, skipping plot generation")
                return None

            result = subprocess.run(
                ["Rscript", r_file],
                cwd=self.config.output_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            self.logger.info(" R图表生成成功|R plots generated successfully")
            if result.stdout:
                self.logger.info(f"R脚本输出|R script output:\n{result.stdout}")

        except FileNotFoundError:
            self.logger.warning(" Rscript未安装，跳过图表生成|Rscript not installed, skipping plot generation")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"R图表生成失败|R plot generation failed.")
            self.logger.warning(f"R脚本错误信息|R script error message:\n{e.stderr}")
        
        return r_file
    
    def _create_r_script(self, best_k: int):
        """ 创建R脚本|Create R script"""
        return f'''#!/usr/bin/env Rscript
#  ADMIXTURE结果可视化脚本|ADMIXTURE Results Visualization Script

# 检查并安装必要的包|Check and install necessary packages
packages <- c("ggplot2", "reshape2", "dplyr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="https://cran.rstudio.com/")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

# 读取数据|Read data
q_file_path <- paste0("{self.config.base_name}.{best_k}.Q")
fam_file_path <- paste0("{self.config.base_name}.fam")
cv_file_path <- "cv_results.csv"

if (!file.exists(q_file_path) || !file.exists(fam_file_path)) {{
    stop(" Q或FAM文件不存在|Q or FAM file not found.")
}}

q_data <- read.table(q_file_path)
colnames(q_data) <- paste0("Pop", 1:{best_k})

#  添加个体信息|Add individual information
fam_data <- read.table(fam_file_path)
q_data$Individual_ID <- 1:nrow(q_data)
q_data$FID <- fam_data$V1
q_data$IID <- fam_data$V2

# 1.  个体祖先成分条形图|Individual ancestry bar plot
pdf("admixture_barplot.pdf", width = 12, height = 6)
q_long <- melt(q_data, id.vars = c("Individual_ID", "FID", "IID"))

p1 <- ggplot(q_long, aes(x = Individual_ID, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Individual Ancestry Proportions (K={best_k})",
       x = "Individual", y = "Ancestry Proportion",
       fill = "Ancestral Population") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
print(p1)
dev.off()

# 2. 交叉验证误差图|Cross-validation error plot
if (file.exists(cv_file_path)) {{
    cv_data <- read.csv(cv_file_path)

    pdf("cv_error_plot.pdf", width = 8, height = 6)
    p2 <- ggplot(cv_data, aes(x = K, y = CV_error)) +
      geom_line(color = "blue") +
      geom_point(color = "red", size = 3) +
      geom_vline(xintercept = {best_k}, linetype = "dashed", color = "red") +
      annotate("text", x = {best_k}, y = max(cv_data$CV_error), label = paste("Best K = {best_k}"), vjust = -0.5, color="red") +
      labs(title = "Cross-Validation Error vs K",
           x = "K (Number of Ancestral Populations)", y = "Cross-Validation Error") +
      theme_bw() +
      scale_x_continuous(breaks = min(cv_data$K):max(cv_data$K))
    print(p2)
    dev.off()
}}

# 3. 混合程度分布图|Admixture level distribution
max_ancestry <- apply(q_data[,1:{best_k}], 1, max)
admixture_level <- 1 - max_ancestry

pdf("admixture_distribution.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

hist(admixture_level, 
     main = "Admixture Level Distribution", 
     xlab = "Admixture Level (1 - Max Ancestry)", 
     ylab = "Number of Individuals",
     col = "lightblue", 
     breaks = 20)

hist(max_ancestry, 
     main = "Maximum Ancestry Distribution", 
     xlab = "Maximum Ancestry Proportion", 
     ylab = "Number of Individuals",
     col = "lightgreen", 
     breaks = 20)

dev.off()

cat(" 图表已生成|Plots generated:\\n")
cat("- admixture_barplot.pdf\\n")
if(file.exists(cv_file_path)) {{ cat("- cv_error_plot.pdf\\n") }}
cat("- admixture_distribution.pdf\\n")
'''

class SummaryGenerator:
    """ 总结生成器|Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, best_k: int, stats: dict):
        """ 生成分析总结|Generate analysis summary"""
        summary_file = os.path.join(self.config.output_dir, "analysis_summary.txt")
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(" ADMIXTURE分析总结报告|ADMIXTURE Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(" 输入文件|Input Files:\n")
            f.write(f"  -  VCF文件|VCF file: {self.config.vcf_file}\n\n")
            
            f.write(" 分析参数|Analysis Parameters:\n")
            f.write(f"  -  K值范围|K range: {self.config.min_k} - {self.config.max_k}\n")
            f.write(f"  -  交叉验证折数|CV folds: {self.config.cv_folds}\n")
            f.write(f"  -  线程数|Threads: {self.config.threads}\n")
            f.write(f"  -  MAF阈值|MAF threshold: {self.config.maf}\n")
            f.write(f"  -  缺失率阈值|Missing rate threshold: {self.config.missing_rate}\n")
            f.write(f"  -  HWE p值阈值|HWE p-value threshold: {self.config.hwe_pvalue}\n\n")
            
            f.write(" 分析结果|Analysis Results:\n")
            f.write(f"  -  最优K值|Best K value (lowest CV error): {best_k}\n")
            f.write(f"  -  总个体数|Total individuals: {stats['total_individuals']}\n")
            f.write(f"  -  高度混合个体数 (max ancestry < 0.7)|Highly admixed individuals: {stats['highly_admixed']}\n")
            f.write(f"  -  纯合个体数 (max ancestry > 0.9)|Pure individuals: {stats['pure_individuals']}\n")
            f.write(f"  -  平均混合程度|Mean admixture level: {stats['mean_admixture_level']:.4f}\n")
            f.write(f"  -  平均最大祖先成分|Mean max ancestry: {stats['mean_max_ancestry']:.4f}\n\n")
            
            f.write("输出文件|Output Files:\n")
            f.write(f"  -  个体祖先成分|Individual ancestry proportions: admixture_proportions.csv\n")
            f.write(f"  -  GWAS协变量|GWAS covariates: gwas_covariates.txt\n")
            f.write(f"  -  交叉验证结果|Cross-validation results: cv_results.csv\n")
            f.write(f"  -  可视化图表|Visualization plots: *.pdf\n")
            f.write(f"  -  分析总结|Analysis Summary: analysis_summary.txt\n")
            f.write(f"输出目录|Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"分析总结已保存|Analysis summary saved: {summary_file}")
        return summary_file