# """
# PLINK GWAS结果处理模块 | PLINK GWAS Results Processing Module
# """

# import pandas as pd
# import numpy as np
# from pathlib import Path
# from datetime import datetime

# try:
#     from scipy.stats import false_discovery_control
#     HAS_SCIPY = True
# except ImportError:
#     HAS_SCIPY = False

# class ResultsProcessor:
#     """结果处理器 | Results Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def process_results(self, main_result: str) -> pd.DataFrame:
#         """处理分析结果 | Process analysis results"""
#         self.logger.info("处理分析结果 | Processing analysis results...")
        
#         # 根据表型类型选择结果文件 | Select result file based on trait type
#         if self.config.trait_type == "qualitative":
#             result_file = f"{main_result}.assoc.logistic"
#         else:
#             result_file = f"{main_result}.assoc.linear"
        
#         if not Path(result_file).exists():
#             raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
#         # 读取结果文件 | Read results file
#         try:
#             results_df = pd.read_csv(result_file, sep=r"\s+")
#         except Exception as e:
#             self.logger.error(f"读取结果文件失败 | Failed to read results file: {e}")
#             return None
        
#         # 提取ADD模型结果 | Extract ADD model results
#         add_results = results_df[results_df['TEST'] == 'ADD'].copy()
        
#         # 移除无效的P值 | Remove invalid P values
#         add_results = add_results.dropna(subset=['P'])
#         add_results = add_results[add_results['P'] > 0]
        
#         add_results.to_csv("gwas_results_ADD.txt", sep='\t', index=False)
        
#         # 应用显著性校正 | Apply significance correction
#         self._apply_significance_correction(add_results)
        
#         return add_results
    
#     def _apply_significance_correction(self, results_df: pd.DataFrame):
#         """应用显著性校正方法 | Apply significance correction methods"""
#         self.logger.info(f"应用显著性校正方法 | Applying significance correction method: {self.config.correction_method}")
        
#         total_snps = len(results_df)
#         self.logger.info(f"用于校正的SNP总数 | Total SNPs for correction: {total_snps}")
        
#         # 根据用户选择的校正方法 | Based on user selected correction method
#         if self.config.correction_method == "all":
#             # 应用所有三种方法 | Apply all three methods
#             self._bonferroni_correction(results_df, total_snps)
#             self._suggestive_correction(results_df)
#             self._fdr_correction(results_df)
#         elif self.config.correction_method == "bonferroni":
#             self._bonferroni_correction(results_df, total_snps)
#         elif self.config.correction_method == "suggestive":
#             self._suggestive_correction(results_df)
#         elif self.config.correction_method == "fdr":
#             self._fdr_correction(results_df)
    
#     def _bonferroni_correction(self, results_df: pd.DataFrame, total_snps: int):
#         """Bonferroni校正 | Bonferroni correction"""
#         self.logger.info("执行Bonferroni校正 | Performing Bonferroni correction")
        
#         bonferroni_threshold = self.config.bonferroni_alpha / total_snps
#         significant_bonferroni = results_df[results_df['P'] < bonferroni_threshold].copy()
        
#         # 按P值排序 | Sort by P value
#         significant_bonferroni = significant_bonferroni.sort_values('P')
        
#         # 保存结果 | Save results
#         significant_bonferroni.to_csv("significant_bonferroni.txt", sep='\t', index=False)
        
#         self.logger.info(f"Bonferroni校正阈值 | Bonferroni threshold: {bonferroni_threshold:.2e}")
#         self.logger.info(f"Bonferroni显著SNP数 | Bonferroni significant SNPs: {len(significant_bonferroni)}")
        
#         # 保存阈值信息 | Save threshold information
#         with open("bonferroni_info.txt", "w") as f:
#             f.write(f"Bonferroni Correction Information\n")
#             f.write(f"Alpha level: {self.config.bonferroni_alpha}\n")
#             f.write(f"Total SNPs: {total_snps}\n")
#             f.write(f"Corrected threshold: {bonferroni_threshold:.2e}\n")
#             f.write(f"Significant SNPs: {len(significant_bonferroni)}\n")
    
#     def _suggestive_correction(self, results_df: pd.DataFrame):
#         """提示性关联阈值 | Suggestive significance threshold"""
#         self.logger.info("执行提示性关联筛选 | Performing suggestive significance filtering")
        
#         suggestive_threshold = self.config.suggestive_threshold
#         significant_suggestive = results_df[results_df['P'] < suggestive_threshold].copy()
        
#         # 按P值排序 | Sort by P value
#         significant_suggestive = significant_suggestive.sort_values('P')
        
#         # 保存结果 | Save results
#         significant_suggestive.to_csv("significant_suggestive.txt", sep='\t', index=False)
        
#         self.logger.info(f"提示性关联阈值 | Suggestive threshold: {suggestive_threshold:.2e}")
#         self.logger.info(f"提示性关联SNP数 | Suggestive significant SNPs: {len(significant_suggestive)}")
        
#         # 保存阈值信息 | Save threshold information
#         with open("suggestive_info.txt", "w") as f:
#             f.write(f"Suggestive Significance Information\n")
#             f.write(f"Threshold: {suggestive_threshold:.2e}\n")
#             f.write(f"Significant SNPs: {len(significant_suggestive)}\n")
    
#     def _fdr_correction(self, results_df: pd.DataFrame):
#         """FDR校正 | FDR correction"""
#         self.logger.info("执行FDR校正 | Performing FDR correction")
        
#         if not HAS_SCIPY:
#             self.logger.warning("未安装scipy，跳过FDR校正 | scipy not installed, skipping FDR correction")
#             self.logger.warning("请安装scipy: pip install scipy | Please install scipy: pip install scipy")
#             return
        
#         # 获取P值 | Get P values
#         p_values = results_df['P'].values
        
#         try:
#             # 使用scipy进行FDR校正 | Use scipy for FDR correction
#             rejected, p_corrected = false_discovery_control(p_values, alpha=self.config.fdr_alpha, method='bh')
            
#             # 添加校正后的P值和是否显著的标记 | Add corrected P values and significance flag
#             results_df['P_FDR'] = p_corrected
#             results_df['FDR_significant'] = rejected
            
#             # 筛选显著的SNP | Filter significant SNPs
#             significant_fdr = results_df[results_df['FDR_significant']].copy()
            
#             # 按原始P值排序 | Sort by original P value
#             significant_fdr = significant_fdr.sort_values('P')
            
#             # 保存结果 | Save results
#             significant_fdr.to_csv("significant_fdr.txt", sep='\t', index=False)
            
#             self.logger.info(f"FDR校正q值阈值 | FDR q-value threshold: {self.config.fdr_alpha}")
#             self.logger.info(f"FDR显著SNP数 | FDR significant SNPs: {len(significant_fdr)}")
            
#             # 保存阈值信息 | Save threshold information
#             with open("fdr_info.txt", "w") as f:
#                 f.write(f"FDR Correction Information\n")
#                 f.write(f"Method: Benjamini-Hochberg\n")
#                 f.write(f"Alpha (q-value): {self.config.fdr_alpha}\n")
#                 f.write(f"Significant SNPs: {len(significant_fdr)}\n")
#                 f.write(f"Highest corrected P-value: {significant_fdr['P_FDR'].max():.2e}\n" if len(significant_fdr) > 0 else "No significant SNPs\n")
            
#         except Exception as e:
#             self.logger.error(f"FDR校正失败 | FDR correction failed: {e}")

# class ReportGenerator:
#     """报告生成器 | Report Generator"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_report(self, results_df: pd.DataFrame, main_result: str):
#         """生成分析报告 | Generate analysis report"""
#         self.logger.info("生成分析报告 | Generating analysis report...")
        
#         # 统计信息 | Statistics
#         stats = {
#             'analysis_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
#             'main_analysis': main_result,
#             'trait_type': self.config.trait_type,
#             'correction_method': self.config.correction_method,
#             'total_snps': 0,
#             'total_samples': 0,
#             'cases': 0,
#             'controls': 0,
#             'chromosomes': 0
#         }
        
#         # 读取文件统计 | Read file statistics
#         if Path("data_qc1.bim").exists():
#             bim_df = pd.read_csv("data_qc1.bim", sep="\t", header=None)
#             stats['total_snps'] = len(bim_df)
#             stats['chromosomes'] = bim_df.iloc[:, 0].nunique()
        
#         if Path("data_qc1.fam").exists():
#             fam_df = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
#             stats['total_samples'] = len(fam_df)
#             if self.config.trait_type == "qualitative":
#                 stats['cases'] = (fam_df.iloc[:, 5] == 2).sum()
#                 stats['controls'] = (fam_df.iloc[:, 5] == 1).sum()
#             else:
#                 # 数量性状统计 | Quantitative trait statistics
#                 pheno_values = fam_df.iloc[:, 5][fam_df.iloc[:, 5] != -9]
#                 stats['mean_phenotype'] = pheno_values.mean()
#                 stats['std_phenotype'] = pheno_values.std()
        
#         # 统计各种校正方法的结果 | Count results for each correction method
#         correction_results = {}
        
#         if self.config.correction_method in ["bonferroni", "all"] and Path("significant_bonferroni.txt").exists():
#             bonferroni_df = pd.read_csv("significant_bonferroni.txt", sep='\t')
#             correction_results['bonferroni'] = len(bonferroni_df)
        
#         if self.config.correction_method in ["suggestive", "all"] and Path("significant_suggestive.txt").exists():
#             suggestive_df = pd.read_csv("significant_suggestive.txt", sep='\t')
#             correction_results['suggestive'] = len(suggestive_df)
        
#         if self.config.correction_method in ["fdr", "all"] and Path("significant_fdr.txt").exists():
#             fdr_df = pd.read_csv("significant_fdr.txt", sep='\t')
#             correction_results['fdr'] = len(fdr_df)
        
#         # 生成报告 | Generate report
#         report = f"""
# === PLINK GWAS分析完成报告 | PLINK GWAS Analysis Report ===
# 分析时间 | Analysis time: {stats['analysis_time']}
# 主要分析 | Main analysis: {stats['main_analysis']}
# 表型类型 | Trait type: {stats['trait_type']}
# 显著性校正方法 | Significance correction method: {stats['correction_method']}
# 染色体类型 | Chromosome type: 支持非标准编号（如OV开头） | Support non-standard naming (e.g., OV prefix)

# === 样本统计 | Sample Statistics ===
# 质控后总样本数 | Total samples after QC: {stats['total_samples']}"""
        
#         if self.config.trait_type == "qualitative":
#             report += f"""
# 病例数 (感病/原值1) | Cases (susceptible/original 1): {stats['cases']}
# 对照数 (抗病/原值0) | Controls (resistant/original 0): {stats['controls']}"""
#         else:
#             report += f"""
# 表型均值 | Phenotype mean: {stats.get('mean_phenotype', 'N/A'):.4f}
# 表型标准差 | Phenotype std: {stats.get('std_phenotype', 'N/A'):.4f}"""
        
#         report += f"""

# === SNP统计 | SNP Statistics ===
# 质控后总SNP数 | Total SNPs after QC: {stats['total_snps']}
# 染色体数 | Number of chromosomes: {stats['chromosomes']}

# === 显著性校正结果 | Significance Correction Results ==="""
        
#         if 'bonferroni' in correction_results:
#             bonferroni_threshold = self.config.bonferroni_alpha / stats['total_snps']
#             report += f"""
# Bonferroni校正结果 | Bonferroni Correction Results:
#   - 校正阈值 | Corrected threshold: {bonferroni_threshold:.2e}
#   - 显著SNP数 | Significant SNPs: {correction_results['bonferroni']}"""
        
#         if 'suggestive' in correction_results:
#             report += f"""
# 提示性关联结果 | Suggestive Association Results:
#   - 阈值 | Threshold: {self.config.suggestive_threshold:.2e}
#   - 显著SNP数 | Significant SNPs: {correction_results['suggestive']}"""
        
#         if 'fdr' in correction_results:
#             report += f"""
# FDR校正结果 | FDR Correction Results:
#   - q值阈值 | q-value threshold: {self.config.fdr_alpha}
#   - 显著SNP数 | Significant SNPs: {correction_results['fdr']}"""
        
#         report += f"""

# === 分析参数 | Analysis Parameters ===
# 表型类型 | Trait type: {self.config.trait_type}
# 显著性校正方法 | Correction method: {self.config.correction_method}
# 个体缺失率阈值 | Individual missing rate threshold: {self.config.mind}
# SNP缺失率阈值 | SNP missing rate threshold: {self.config.geno}
# 最小等位基因频率 | Minor allele frequency: {self.config.maf}
# Hardy-Weinberg平衡P值 | Hardy-Weinberg equilibrium p-value: {self.config.hwe}
# LD r²阈值 | LD r² threshold: {self.config.ld_r2_threshold}
# 使用主成分数 | Number of PCs used: {self.config.pca_use}

# === 校正方法说明 | Correction Method Description ==="""
        
#         if self.config.correction_method in ["bonferroni", "all"]:
#             report += f"""
# Bonferroni校正: 控制家族错误率(FWER)，阈值 = {self.config.bonferroni_alpha}/总SNP数
# Bonferroni correction: Controls family-wise error rate (FWER), threshold = {self.config.bonferroni_alpha}/total SNPs"""
        
#         if self.config.correction_method in ["suggestive", "all"]:
#             report += f"""
# 提示性关联: 传统GWAS中的提示性显著性阈值，P < {self.config.suggestive_threshold:.0e}
# Suggestive association: Traditional suggestive significance threshold in GWAS, P < {self.config.suggestive_threshold:.0e}"""
        
#         if self.config.correction_method in ["fdr", "all"]:
#             report += f"""
# FDR校正: 控制假发现率，使用Benjamini-Hochberg方法，q值 < {self.config.fdr_alpha}
# FDR correction: Controls false discovery rate using Benjamini-Hochberg method, q-value < {self.config.fdr_alpha}"""
        
#         report += f"""

# === 结果解释 | Results Interpretation ==="""
        
#         if self.config.trait_type == "qualitative":
#             report += """
# - OR > 1: 该等位基因增加感病风险 | This allele increases susceptibility risk
# - OR < 1: 该等位基因具有保护作用（促进抗病） | This allele has protective effect (promotes resistance)"""
#         else:
#             report += """
# - BETA > 0: 该等位基因增加表型值 | This allele increases phenotype value
# - BETA < 0: 该等位基因降低表型值 | This allele decreases phenotype value"""
        
#         report += """
# - 显著位点需要在独立群体中验证 | Significant loci need validation in independent populations

# === 输出文件 | Output Files ===
# - gwas_results_ADD.txt: 主要关联分析结果 | Main association analysis results"""
        
#         if self.config.correction_method in ["bonferroni", "all"]:
#             report += """
# - significant_bonferroni.txt: Bonferroni校正显著位点 | Bonferroni significant loci
# - bonferroni_info.txt: Bonferroni校正信息 | Bonferroni correction information"""
        
#         if self.config.correction_method in ["suggestive", "all"]:
#             report += """
# - significant_suggestive.txt: 提示性关联位点 | Suggestive association loci
# - suggestive_info.txt: 提示性关联信息 | Suggestive association information"""
        
#         if self.config.correction_method in ["fdr", "all"]:
#             report += """
# - significant_fdr.txt: FDR校正显著位点 | FDR significant loci
# - fdr_info.txt: FDR校正信息 | FDR correction information"""
        
#         report += """
# - manhattan_plot.png: Manhattan图 | Manhattan plot
# - qq_plot.png: QQ图 | QQ plot
# """
        
#         # 如果有显著位点，显示每种方法的top hits | Show top hits for each method if significant loci exist
#         if correction_results:
#             report += "\n=== 各校正方法的前5个最显著位点 | Top 5 Most Significant Loci by Each Method ===\n"
            
#             for method, count in correction_results.items():
#                 if count > 0:
#                     method_file = f"significant_{method}.txt"
#                     if Path(method_file).exists():
#                         method_df = pd.read_csv(method_file, sep='\t')
#                         top_hits = method_df.head(5)
                        
#                         report += f"\n{method.upper()}方法 | {method.upper()} Method:\n"
#                         if self.config.trait_type == "qualitative":
#                             report += "CHR\tSNP\tBP\tA1\tOR\tP\n"
#                             for _, row in top_hits.iterrows():
#                                 report += f"{row['CHR']}\t{row['SNP']}\t{row['BP']}\t{row['A1']}\t{row['OR']:.3f}\t{row['P']:.2e}\n"
#                         else:
#                             report += "CHR\tSNP\tBP\tA1\tBETA\tP\n"
#                             for _, row in top_hits.iterrows():
#                                 report += f"{row['CHR']}\t{row['SNP']}\t{row['BP']}\t{row['A1']}\t{row['BETA']:.3f}\t{row['P']:.2e}\n"
        
#         # 保存报告 | Save report
#         with open("analysis_report.txt", "w", encoding='utf-8') as f:
#             f.write(report)
        
#         self.logger.info("分析报告已生成 | Analysis report generated: analysis_report.txt")

# class VisualizationGenerator:
#     """可视化生成器 | Visualization Generator"""
    
#     def __init__(self, config, logger, cmd_runner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def generate_plots(self):
#         """生成可视化图形 | Generate visualization plots"""
#         self.logger.info("生成可视化图形 | Generating visualization plots...")
        
#         try:
#             # 创建R脚本 | Create R script
#             r_script = f'''
# # 安装和载入必要的包 | Install and load required packages
# if (!require("qqman", quietly=TRUE)) {{
#     install.packages("qqman", repos="https://cran.rstudio.com/")
#     library(qqman)
# }}

# # 检查结果文件是否存在 | Check if results file exists
# if (!file.exists("gwas_results_ADD.txt")) {{
#     stop("错误: 找不到 gwas_results_ADD.txt 文件 | Error: Cannot find gwas_results_ADD.txt file")
# }}

# # 读取GWAS结果 | Read GWAS results
# cat("读取GWAS结果... | Reading GWAS results...\\n")
# gwas <- read.table("gwas_results_ADD.txt", header=TRUE, stringsAsFactors=FALSE)

# if (nrow(gwas) == 0) {{
#     stop("错误: GWAS结果文件为空 | Error: GWAS results file is empty")
# }}

# cat("数据行数 | Data rows:", nrow(gwas), "\\n")
# cat("表型类型 | Trait type: {self.config.trait_type}\\n")
# cat("校正方法 | Correction method: {self.config.correction_method}\\n")

# # 处理非标准染色体编号 | Handle non-standard chromosome naming
# gwas$CHR_original <- gwas$CHR
# gwas$CHR_numeric <- as.numeric(gsub("OV0?", "", gwas$CHR))

# # 如果转换失败，使用序号 | If conversion fails, use index
# if (any(is.na(gwas$CHR_numeric))) {{
#     unique_chrs <- unique(gwas$CHR)
#     chr_mapping <- setNames(1:length(unique_chrs), unique_chrs)
#     gwas$CHR_numeric <- chr_mapping[gwas$CHR]
# }}

# # 准备绘图数据 | Prepare plotting data
# gwas_plot <- data.frame(
#     SNP = gwas$SNP,
#     CHR = gwas$CHR_numeric,
#     BP = gwas$BP,
#     P = gwas$P
# )

# # 移除无效数据 | Remove invalid data
# gwas_plot <- gwas_plot[!is.na(gwas_plot$P) & gwas_plot$P > 0 & gwas_plot$P <= 1, ]

# if (nrow(gwas_plot) == 0) {{
#     stop("错误: 没有有效的P值数据 | Error: No valid P-value data")
# }}

# cat("有效数据点数 | Valid data points:", nrow(gwas_plot), "\\n")

# # 设置显著性阈值 | Set significance thresholds
# total_snps <- nrow(gwas_plot)
# bonferroni_threshold <- {self.config.bonferroni_alpha} / total_snps
# suggestive_threshold <- {self.config.suggestive_threshold}

# # 生成Manhattan图 | Generate Manhattan plot
# cat("生成Manhattan图... | Generating Manhattan plot...\\n")
# trait_label <- if ("{self.config.trait_type}" == "qualitative") "Disease Resistance (Qualitative)" else "Quantitative Trait"
# correction_label <- "{self.config.correction_method}"

# png("manhattan_plot.png", width=1400, height=700, res=100)
# tryCatch({{
#     manhattan(gwas_plot, 
#               main=paste("PLINK GWAS Manhattan Plot -", trait_label, "- Correction:", correction_label), 
#               genomewideline = -log10(bonferroni_threshold),
#               suggestiveline = -log10(suggestive_threshold),
#               col = c("blue4", "orange3"))
    
#     # 添加阈值标签 | Add threshold labels
#     abline(h = -log10(bonferroni_threshold), col = "red", lty = 2, lwd = 1)
#     abline(h = -log10(suggestive_threshold), col = "blue", lty = 2, lwd = 1)
    
#     # 添加图例 | Add legend
#     legend("topright", 
#            legend = c(paste("Bonferroni:", format(bonferroni_threshold, scientific=TRUE, digits=2)),
#                      paste("Suggestive:", format(suggestive_threshold, scientific=TRUE, digits=2))),
#            col = c("red", "blue"), lty = 2, cex = 0.8)
# }}, error = function(e) {{
#     plot(1, type="n", main="Manhattan Plot Error", xlab="Position", ylab="-log10(P)")
#     text(1, 1, paste("Error:", e$message), cex=0.8)
# }})
# dev.off()

# # 生成QQ图 | Generate QQ plot
# cat("生成QQ图... | Generating QQ plot...\\n")
# png("qq_plot.png", width=600, height=600, res=100)
# tryCatch({{
#     qq(gwas_plot$P, main=paste("QQ Plot -", trait_label, "- Correction:", correction_label))
# }}, error = function(e) {{
#     plot(1, type="n", main="QQ Plot Error", xlab="Expected", ylab="Observed")
#     text(1, 1, paste("Error:", e$message), cex=0.8)
# }})
# dev.off()

# # 保存染色体映射 | Save chromosome mapping
# write.table(data.frame(
#     Original_CHR = unique(gwas$CHR_original),
#     Numeric_CHR = unique(gwas$CHR_numeric)
# ), "chromosome_mapping.txt", row.names=FALSE, quote=FALSE, sep="\\t")

# # 保存阈值信息 | Save threshold information
# cat("\\n=== 显著性阈值信息 | Significance Threshold Information ===\\n")
# cat("Bonferroni阈值 | Bonferroni threshold:", format(bonferroni_threshold, scientific=TRUE, digits=3), "\\n")
# cat("提示性阈值 | Suggestive threshold:", format(suggestive_threshold, scientific=TRUE, digits=3), "\\n")

# # 统计各阈值下的显著SNP数 | Count significant SNPs under each threshold
# bonferroni_count <- sum(gwas_plot$P < bonferroni_threshold)
# suggestive_count <- sum(gwas_plot$P < suggestive_threshold)

# cat("Bonferroni显著SNP数 | Bonferroni significant SNPs:", bonferroni_count, "\\n")
# cat("提示性显著SNP数 | Suggestive significant SNPs:", suggestive_count, "\\n")

# cat("\\n=== 可视化完成 | Visualization completed ===\\n")
# cat("生成的文件 | Generated files:\\n")
# cat("- manhattan_plot.png\\n") 
# cat("- qq_plot.png\\n")
# cat("- chromosome_mapping.txt\\n")

# cat("\\n=== 数据统计 | Data statistics ===\\n")
# cat("总SNP数 | Total SNPs:", nrow(gwas_plot), "\\n")
# cat("染色体数 | Number of chromosomes:", length(unique(gwas_plot$CHR)), "\\n")
# cat("最小P值 | Minimum P-value:", min(gwas_plot$P), "\\n")
# '''
            
#             # 保存R脚本 | Save R script
#             with open("plot_results.R", "w") as f:
#                 f.write(r_script)
            
#             # 运行R脚本 | Run R script
#             result = self.cmd_runner.run(["Rscript", "plot_results.R"], 
#                                        "生成可视化图形 | Generating visualization plots", 
#                                        check=False)
            
#             if result.returncode == 0:
#                 self.logger.info("可视化图形生成成功 | Visualization plots generated successfully")
#             else:
#                 self.logger.warning("可视化图形生成失败，请手动运行 | Visualization generation failed, please run manually: Rscript plot_results.R")
        
#         except Exception as e:
#             self.logger.error(f"生成可视化脚本失败 | Failed to generate visualization script: {e}")

# 20250726添加显性模型和隐性模型的结果处理 | Added processing for dominant and recessive models
"""
PLINK GWAS结果处理模块 | PLINK GWAS Results Processing Module
"""

# import pandas as pd
# import numpy as np
# import os
# from pathlib import Path
# from datetime import datetime

# try:
#     from scipy.stats import false_discovery_control
#     HAS_SCIPY = True
# except ImportError:
#     HAS_SCIPY = False

# class ResultsProcessor:
#     """结果处理器 | Results Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def process_results(self, main_result: str) -> dict:
#         """处理分析结果 | Process analysis results"""
#         self.logger.info("🔄 处理分析结果 | Processing analysis results...")
#         self.logger.info(f"📁 主结果文件前缀 | Main result file prefix: {main_result}")
        
#         if self.config.genetic_model == "all":
#             return self._process_all_models_results(main_result)
#         else:
#             return self._process_single_model_results(main_result)
    
#     def _process_all_models_results(self, main_result: str) -> dict:
#         """处理所有模型的结果 | Process all models results"""
#         self.logger.info("📊 处理所有遗传模型结果 | Processing all genetic models results...")
        
#         result_file = f"{main_result}.model"
#         if not os.path.exists(result_file):
#             self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
#             raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
#         # 检查文件大小 | Check file size
#         file_size_mb = os.path.getsize(result_file) / (1024 * 1024)
#         self.logger.info(f"📏 结果文件大小 | Result file size: {file_size_mb:.1f} MB")
        
#         # 读取结果文件 | Read results file
#         self.logger.info("📖 读取.model文件，这可能需要一些时间... | Reading .model file, this may take some time...")
#         try:
#             results_df = pd.read_csv(result_file, sep=r"\s+")
#             self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
#         except Exception as e:
#             self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
#             return {}
        
#         # 检查数据结构 | Check data structure
#         self.logger.info(f"📋 结果文件列名 | Result file columns: {list(results_df.columns)}")
#         unique_tests = results_df['TEST'].unique()
#         self.logger.info(f"🧪 包含的检验类型 | Test types included: {unique_tests}")
        
#         # 分别处理各种模型的结果 | Process results for different models separately
#         model_results = {}
        
#         # 处理加性模型结果 | Process additive model results
#         if 'ADD' in results_df['TEST'].values:
#             add_results = results_df[results_df['TEST'] == 'ADD'].copy()
#             add_results = self._clean_results(add_results)
#             if len(add_results) > 0:
#                 add_results.to_csv("gwas_results_ADD.txt", sep='\t', index=False)
#                 self._apply_significance_correction(add_results, "ADD")
#                 model_results['additive'] = add_results
#                 self.logger.info(f"✅ 加性模型结果 | Additive model results: {len(add_results)} SNPs")
#             else:
#                 self.logger.warning("⚠️ 加性模型无有效结果 | No valid additive model results")
        
#         # 处理显性模型结果 | Process dominant model results
#         if 'DOMDEV' in results_df['TEST'].values:
#             dom_results = results_df[results_df['TEST'] == 'DOMDEV'].copy()
#             dom_results = self._clean_results(dom_results)
#             if len(dom_results) > 0:
#                 dom_results.to_csv("gwas_results_DOM.txt", sep='\t', index=False)
#                 self._apply_significance_correction(dom_results, "DOM")
#                 model_results['dominant'] = dom_results
#                 self.logger.info(f"✅ 显性模型结果 | Dominant model results: {len(dom_results)} SNPs")
#             else:
#                 self.logger.warning("⚠️ 显性模型无有效结果 | No valid dominant model results")
        
#         # 处理基因型模型结果 | Process genotypic model results
#         if 'GENO_2DF' in results_df['TEST'].values:
#             geno_results = results_df[results_df['TEST'] == 'GENO_2DF'].copy()
#             geno_results = self._clean_results(geno_results)
#             if len(geno_results) > 0:
#                 geno_results.to_csv("gwas_results_GENO.txt", sep='\t', index=False)
#                 self._apply_significance_correction(geno_results, "GENO")
#                 model_results['genotypic'] = geno_results
#                 self.logger.info(f"✅ 基因型模型结果 | Genotypic model results: {len(geno_results)} SNPs")
#             else:
#                 self.logger.warning("⚠️ 基因型模型无有效结果 | No valid genotypic model results")
        
#         # 处理隐性模型结果（如果存在）| Process recessive model results (if exists)
#         if 'REC' in results_df['TEST'].values:
#             rec_results = results_df[results_df['TEST'] == 'REC'].copy()
#             rec_results = self._clean_results(rec_results)
#             if len(rec_results) > 0:
#                 rec_results.to_csv("gwas_results_REC.txt", sep='\t', index=False)
#                 self._apply_significance_correction(rec_results, "REC")
#                 model_results['recessive'] = rec_results
#                 self.logger.info(f"✅ 隐性模型结果 | Recessive model results: {len(rec_results)} SNPs")
        
#         # 生成综合比较报告 | Generate comprehensive comparison report
#         if model_results:
#             self._generate_model_comparison_report(model_results)
#         else:
#             self.logger.warning("⚠️ 未找到有效的模型结果 | No valid model results found")
        
#         return model_results
    
#     def _process_single_model_results(self, main_result: str) -> dict:
#         """处理单一模型的结果 | Process single model results"""
#         self.logger.info(f"📊 处理{self.config.genetic_model}模型结果 | Processing {self.config.genetic_model} model results...")
        
#         # 根据表型类型确定结果文件 | Determine result file based on trait type
#         if self.config.trait_type == "qualitative":
#             result_file = f"{main_result}.assoc.logistic"
#         else:
#             result_file = f"{main_result}.assoc.linear"
        
#         if not os.path.exists(result_file):
#             self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
#             raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
#         # 读取结果文件 | Read results file
#         try:
#             results_df = pd.read_csv(result_file, sep=r"\s+")
#             self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
#         except Exception as e:
#             self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
#             return {}
        
#         # 根据遗传模型提取相应结果 | Extract results based on genetic model
#         if self.config.genetic_model == "additive":
#             model_results = results_df[results_df['TEST'] == 'ADD'].copy()
#             file_suffix = "ADD"
#         elif self.config.genetic_model == "dominant":
#             model_results = results_df[results_df['TEST'] == 'DOMDEV'].copy()
#             file_suffix = "DOM"
#         elif self.config.genetic_model == "recessive":
#             if 'REC' in results_df['TEST'].values:
#                 model_results = results_df[results_df['TEST'] == 'REC'].copy()
#                 file_suffix = "REC"
#             else:
#                 model_results = results_df[results_df['TEST'] == 'GENO_2DF'].copy() if 'GENO_2DF' in results_df['TEST'].values else results_df[results_df['TEST'] == 'ADD'].copy()
#                 file_suffix = "GENO"
#         else:
#             model_results = results_df.copy()
#             file_suffix = "ALL"
        
#         # 清理和保存结果 | Clean and save results
#         model_results = self._clean_results(model_results)
#         if len(model_results) > 0:
#             model_results.to_csv(f"gwas_results_{file_suffix}.txt", sep='\t', index=False)
#             self._apply_significance_correction(model_results, file_suffix)
#             return {self.config.genetic_model: model_results}
#         else:
#             self.logger.warning(f"⚠️ {self.config.genetic_model}模型无有效结果 | No valid {self.config.genetic_model} model results")
#             return {}
    
#     def _clean_results(self, results_df: pd.DataFrame) -> pd.DataFrame:
#         """清理结果数据 | Clean results data"""
#         if len(results_df) == 0:
#             return results_df
        
#         original_count = len(results_df)
        
#         # 移除无效的P值 | Remove invalid P values
#         results_df = results_df.dropna(subset=['P'])
#         results_df = results_df[results_df['P'] > 0]
#         results_df = results_df[results_df['P'] <= 1]
        
#         # 移除缺失的OR值（对于logistic回归）| Remove missing OR values (for logistic regression)
#         if 'OR' in results_df.columns:
#             results_df = results_df.dropna(subset=['OR'])
#             results_df = results_df[results_df['OR'] > 0]
        
#         cleaned_count = len(results_df)
#         removed_count = original_count - cleaned_count
        
#         if removed_count > 0:
#             self.logger.info(f"🧹 数据清理：移除了{removed_count}个无效记录 | Data cleaning: removed {removed_count} invalid records")
        
#         return results_df
    
#     def _apply_significance_correction(self, results_df: pd.DataFrame, model_name: str):
#         """应用显著性校正方法 | Apply significance correction methods"""
#         self.logger.info(f"📈 对{model_name}模型应用显著性校正 | Applying significance correction for {model_name} model: {self.config.correction_method}")
        
#         total_snps = len(results_df)
#         self.logger.info(f"🔢 用于校正的SNP总数 | Total SNPs for correction: {total_snps}")
        
#         if total_snps == 0:
#             self.logger.warning(f"⚠️ {model_name}模型无SNP数据可用于校正 | No SNP data available for correction in {model_name} model")
#             return
        
#         # 根据用户选择的校正方法 | Based on user selected correction method
#         if self.config.correction_method == "all":
#             self._bonferroni_correction(results_df, total_snps, model_name)
#             self._suggestive_correction(results_df, model_name)
#             self._fdr_correction(results_df, model_name)
#         elif self.config.correction_method == "bonferroni":
#             self._bonferroni_correction(results_df, total_snps, model_name)
#         elif self.config.correction_method == "suggestive":
#             self._suggestive_correction(results_df, model_name)
#         elif self.config.correction_method == "fdr":
#             self._fdr_correction(results_df, model_name)
    
#     def _bonferroni_correction(self, results_df: pd.DataFrame, total_snps: int, model_name: str):
#         """Bonferroni校正 | Bonferroni correction"""
#         bonferroni_threshold = self.config.bonferroni_alpha / total_snps
#         significant = results_df[results_df['P'] < bonferroni_threshold]
        
#         self.logger.info(f"🎯 {model_name}模型Bonferroni校正阈值 | {model_name} model Bonferroni threshold: {bonferroni_threshold:.2e}")
#         self.logger.info(f"⭐ {model_name}模型Bonferroni显著位点数 | {model_name} model Bonferroni significant loci: {len(significant)}")
        
#         significant.to_csv(f"significant_hits_bonferroni_{model_name}.txt", sep='\t', index=False)
        
#         # 如果有显著位点，显示前几个 | If there are significant loci, show top ones
#         if len(significant) > 0:
#             top_hits = significant.nsmallest(min(5, len(significant)), 'P')
#             self.logger.info(f"🏆 {model_name}模型前{len(top_hits)}个Bonferroni显著位点 | Top {len(top_hits)} Bonferroni significant loci for {model_name} model:")
#             for idx, row in top_hits.iterrows():
#                 # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
#                 bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
#                 self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}{bp_info}")
    
#     def _suggestive_correction(self, results_df: pd.DataFrame, model_name: str):
#         """提示性关联阈值 | Suggestive association threshold"""
#         suggestive = results_df[results_df['P'] < self.config.suggestive_threshold]
        
#         self.logger.info(f"🎯 {model_name}模型提示性关联阈值 | {model_name} model suggestive threshold: {self.config.suggestive_threshold:.2e}")
#         self.logger.info(f"🔍 {model_name}模型提示性关联位点数 | {model_name} model suggestive association loci: {len(suggestive)}")
        
#         suggestive.to_csv(f"suggestive_hits_{model_name}.txt", sep='\t', index=False)
        
#         # 如果有提示性位点，显示前几个 | If there are suggestive loci, show top ones
#         if len(suggestive) > 0:
#             top_hits = suggestive.nsmallest(min(10, len(suggestive)), 'P')
#             self.logger.info(f"🔍 {model_name}模型前{len(top_hits)}个提示性关联位点 | Top {len(top_hits)} suggestive loci for {model_name} model:")
#             for idx, row in top_hits.iterrows():
#                 # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
#                 bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
#                 self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}{bp_info}")
    
#     def _fdr_correction(self, results_df: pd.DataFrame, model_name: str):
#         """FDR校正 | FDR correction"""
#         if not HAS_SCIPY:
#             self.logger.warning("⚠️ scipy未安装，跳过FDR校正 | scipy not installed, skipping FDR correction")
#             return
        
#         try:
#             p_values = results_df['P'].values
#             fdr_corrected = false_discovery_control(p_values, alpha=self.config.fdr_alpha)
            
#             results_df_fdr = results_df.copy()
#             results_df_fdr['P_FDR'] = fdr_corrected
            
#             fdr_significant = results_df_fdr[results_df_fdr['P_FDR'] < self.config.fdr_alpha]
            
#             self.logger.info(f"🎯 {model_name}模型FDR校正阈值 | {model_name} model FDR threshold: {self.config.fdr_alpha}")
#             self.logger.info(f"📊 {model_name}模型FDR显著位点数 | {model_name} model FDR significant loci: {len(fdr_significant)}")
            
#             results_df_fdr.to_csv(f"gwas_results_fdr_{model_name}.txt", sep='\t', index=False)
#             fdr_significant.to_csv(f"fdr_significant_hits_{model_name}.txt", sep='\t', index=False)
            
#             if len(fdr_significant) > 0:
#                 top_hits = fdr_significant.nsmallest(min(5, len(fdr_significant)), 'P')
#                 self.logger.info(f"📊 {model_name}模型前{len(top_hits)}个FDR显著位点 | Top {len(top_hits)} FDR significant loci for {model_name} model:")
#                 for idx, row in top_hits.iterrows():
#                     # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
#                     bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
#                     self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, FDR_P={row['P_FDR']:.2e}, CHR={row['CHR']}{bp_info}")
#         except Exception as e:
#             self.logger.error(f"❌ FDR校正失败 | FDR correction failed: {e}")
    
#     def _generate_model_comparison_report(self, model_results: dict):
#         """生成模型比较报告 | Generate model comparison report"""
#         self.logger.info("📋 生成模型比较报告 | Generating model comparison report...")
        
#         comparison_report = []
#         comparison_report.append("=" * 80)
#         comparison_report.append("PLINK GWAS 多模型比较报告 | Multi-Model Comparison Report")
#         comparison_report.append("=" * 80)
#         comparison_report.append(f"分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
#         comparison_report.append(f"表型类型 | Trait type: {self.config.trait_type}")
#         comparison_report.append(f"显著性校正方法 | Correction method: {self.config.correction_method}")
#         comparison_report.append("")
        
#         # 统计各模型结果 | Statistics for each model
#         for model_name, results_df in model_results.items():
#             comparison_report.append(f"=== {model_name.upper()}模型结果 | {model_name.upper()} Model Results ===")
#             comparison_report.append(f"总SNP数 | Total SNPs: {len(results_df)}")
            
#             if len(results_df) > 0:
#                 min_p = results_df['P'].min()
#                 comparison_report.append(f"最小P值 | Minimum P-value: {min_p:.2e}")
                
#                 if self.config.correction_method in ["bonferroni", "all"]:
#                     bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
#                     bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
#                     comparison_report.append(f"Bonferroni显著位点 | Bonferroni significant loci (P<{bonferroni_threshold:.2e}): {bonferroni_hits}")
                
#                 if self.config.correction_method in ["suggestive", "all"]:
#                     suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
#                     comparison_report.append(f"提示性关联位点 | Suggestive association loci (P<{self.config.suggestive_threshold:.2e}): {suggestive_hits}")
                
#                 top_snps = results_df.nsmallest(5, 'P')
#                 if len(top_snps) > 0:
#                     comparison_report.append("前5个最显著SNP | Top 5 most significant SNPs:")
#                     for idx, row in top_snps.iterrows():
#                         or_info = f", OR={row['OR']:.3f}" if 'OR' in row and pd.notna(row['OR']) else ""
#                         comparison_report.append(f"  {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}, BP={row['BP']}{or_info}")
#             else:
#                 comparison_report.append("❌ 无有效结果 | No valid results")
            
#             comparison_report.append("")
        
#         # 模型解释 | Model interpretation
#         comparison_report.append("=== 模型解释 | Model Interpretation ===")
#         comparison_report.append("加性模型 (ADD) | Additive Model:")
#         comparison_report.append("  - 假设杂合子效应是两个纯合子效应的中间值")
#         comparison_report.append("  - 编码: AA=0, Aa=1, aa=2")
#         comparison_report.append("  - 适用于大多数复杂性状")
#         comparison_report.append("")
#         comparison_report.append("显性模型 (DOMDEV) | Dominant Model:")
#         comparison_report.append("  - 假设一个风险等位基因拷贝就足以产生效应")
#         comparison_report.append("  - 编码: AA=0, Aa=1, aa=1")
#         comparison_report.append("  - 适用于显性遗传病")
#         comparison_report.append("")
#         comparison_report.append("基因型模型 (GENO_2DF) | Genotypic Model:")
#         comparison_report.append("  - 2自由度检验，不假设特定的遗传模式")
#         comparison_report.append("  - 可以检测非加性效应")
#         comparison_report.append("  - 适用于未知遗传模式的性状")
#         comparison_report.append("")
        
#         # 建议 | Recommendations
#         comparison_report.append("=== 分析建议 | Analysis Recommendations ===")
#         if len(model_results) > 1:
#             best_model = self._identify_best_model(model_results)
#             comparison_report.append(f"推荐模型 | Recommended model: {best_model}")
#             comparison_report.append("基于最小P值和生物学合理性进行选择")
#         comparison_report.append("建议进行独立群体验证 | Recommend validation in independent populations")
#         comparison_report.append("考虑功能验证和精细定位 | Consider functional validation and fine mapping")
        
#         # 保存报告 | Save report
#         report_content = "\n".join(comparison_report)
#         with open("model_comparison_report.txt", 'w', encoding='utf-8') as f:
#             f.write(report_content)
        
#         self.logger.info("📄 模型比较报告已保存 | Model comparison report saved: model_comparison_report.txt")
    
#     def _identify_best_model(self, model_results: dict) -> str:
#         """识别最佳模型 | Identify best model"""
#         min_p_values = {}
        
#         for model_name, results_df in model_results.items():
#             if len(results_df) > 0:
#                 min_p_values[model_name] = results_df['P'].min()
#             else:
#                 min_p_values[model_name] = 1.0
        
#         if min_p_values:
#             best_model = min(min_p_values.items(), key=lambda x: x[1])[0]
#             return best_model
#         else:
#             return "无有效模型 | No valid model"


# class ReportGenerator:
#     """报告生成器 | Report Generator"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_summary_report(self, model_results: dict, stats: dict):
#         """生成总结报告 | Generate summary report"""
#         self.logger.info("📋 生成总结报告 | Generating summary report...")
        
#         report_lines = []
#         report_lines.append("=" * 80)
#         report_lines.append("PLINK GWAS分析总结报告 | PLINK GWAS Analysis Summary Report")
#         report_lines.append("=" * 80)
#         report_lines.append(f"分析完成时间 | Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
#         report_lines.append(f"遗传模型 | Genetic model: {self.config.genetic_model}")
#         report_lines.append(f"表型类型 | Trait type: {self.config.trait_type}")
#         report_lines.append(f"显著性校正方法 | Correction method: {self.config.correction_method}")
#         report_lines.append("")
        
#         # 基本统计 | Basic statistics
#         report_lines.append("=== 基本统计 | Basic Statistics ===")
#         report_lines.append(f"质控后样本数 | Samples after QC: {stats.get('total_samples', 'N/A')}")
#         report_lines.append(f"质控后SNP数 | SNPs after QC: {stats.get('total_snps', 'N/A')}")
        
#         if self.config.trait_type == "qualitative":
#             report_lines.append(f"病例数 | Cases: {stats.get('cases', 'N/A')}")
#             report_lines.append(f"对照数 | Controls: {stats.get('controls', 'N/A')}")
        
#         report_lines.append("")
        
#         # 模型特定结果 | Model-specific results
#         if self.config.genetic_model == "all":
#             report_lines.append("=== 多模型分析结果 | Multi-Model Analysis Results ===")
#             for model_name, results_df in model_results.items():
#                 report_lines.append(f"{model_name.upper()}模型:")
#                 report_lines.append(f"  分析SNP数: {len(results_df)}")
#                 if len(results_df) > 0:
#                     min_p = results_df['P'].min()
#                     report_lines.append(f"  最小P值: {min_p:.2e}")
#                     top_snp = results_df.loc[results_df['P'].idxmin()]
#                     report_lines.append(f"  最显著SNP: {top_snp['SNP']} (CHR{top_snp['CHR']}:{top_snp['BP']})")
                    
#                     if self.config.correction_method in ["bonferroni", "all"]:
#                         bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
#                         bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
#                         report_lines.append(f"  Bonferroni显著位点: {bonferroni_hits}")
                    
#                     if self.config.correction_method in ["suggestive", "all"]:
#                         suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
#                         report_lines.append(f"  提示性关联位点: {suggestive_hits}")
#                 else:
#                     report_lines.append(f"  ❌ 无有效结果")
#                 report_lines.append("")
#         else:
#             report_lines.append(f"=== {self.config.genetic_model.upper()}模型分析结果 | {self.config.genetic_model.upper()} Model Analysis Results ===")
#             if model_results and len(model_results) > 0:
#                 results_df = list(model_results.values())[0]
#                 report_lines.append(f"分析SNP数 | Analyzed SNPs: {len(results_df)}")
#                 if len(results_df) > 0:
#                     min_p = results_df['P'].min()
#                     report_lines.append(f"最小P值 | Minimum P-value: {min_p:.2e}")
                    
#                     if self.config.correction_method in ["bonferroni", "all"]:
#                         bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
#                         bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
#                         report_lines.append(f"Bonferroni显著位点 | Bonferroni significant loci (P<{bonferroni_threshold:.2e}): {bonferroni_hits}")
                    
#                     if self.config.correction_method in ["suggestive", "all"]:
#                         suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
#                         report_lines.append(f"提示性关联位点 | Suggestive association loci (P<{self.config.suggestive_threshold:.2e}): {suggestive_hits}")
#             else:
#                 report_lines.append("❌ 无有效结果 | No valid results")
        
#         report_lines.append("")
        
#         # 质量控制参数 | Quality control parameters
#         report_lines.append("=== 质量控制参数 | Quality Control Parameters ===")
#         report_lines.append(f"个体缺失率阈值 | Individual missing rate: {self.config.mind}")
#         report_lines.append(f"SNP缺失率阈值 | SNP missing rate: {self.config.geno}")
#         report_lines.append(f"最小等位基因频率 | Minor allele frequency: {self.config.maf}")
#         report_lines.append(f"Hardy-Weinberg平衡P值 | HWE p-value: {self.config.hwe}")
#         report_lines.append("")
        
#         # 输出文件列表 | Output file list
#         report_lines.append("=== 输出文件 | Output Files ===")
#         if self.config.genetic_model == "all":
#             report_lines.append("- gwas_results_ADD.txt: 加性模型结果")
#             report_lines.append("- gwas_results_DOM.txt: 显性模型结果")
#             report_lines.append("- gwas_results_GENO.txt: 基因型模型结果")
#             report_lines.append("- model_comparison_report.txt: 模型比较报告")
#         else:
#             model_suffix = {"additive": "ADD", "dominant": "DOM", "recessive": "REC"}.get(self.config.genetic_model, "ALL")
#             report_lines.append(f"- gwas_results_{model_suffix}.txt: {self.config.genetic_model}模型结果")
        
#         if self.config.correction_method in ["bonferroni", "all"]:
#             report_lines.append("- significant_hits_bonferroni_*.txt: Bonferroni校正显著位点")
#         if self.config.correction_method in ["suggestive", "all"]:
#             report_lines.append("- suggestive_hits_*.txt: 提示性关联位点")
#         if self.config.correction_method in ["fdr", "all"]:
#             report_lines.append("- fdr_significant_hits_*.txt: FDR校正显著位点")
        
#         report_lines.append("")
        
#         # 建议和注意事项 | Recommendations and notes
#         report_lines.append("=== 建议和注意事项 | Recommendations and Notes ===")
#         report_lines.append("1. 显著位点需要在独立群体中验证")
#         report_lines.append("2. 考虑进行功能注释和通路分析")
#         report_lines.append("3. 注意群体分层和家系结构的影响")
#         if self.config.genetic_model == "all":
#             report_lines.append("4. 比较不同模型结果，选择生物学上最合理的模型")
#         report_lines.append("5. 建议结合其他组学数据进行综合分析")
        
#         # 保存报告 | Save report
#         report_content = "\n".join(report_lines)
#         with open("gwas_summary_report.txt", 'w', encoding='utf-8') as f:
#             f.write(report_content)
        
#         self.logger.info("📄 GWAS总结报告已保存 | GWAS summary report saved: gwas_summary_report.txt")


# class VisualizationGenerator:
#     """可视化生成器 | Visualization Generator"""
    
#     def __init__(self, config, logger, cmd_runner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def generate_plots(self, model_results: dict):
#         """生成可视化图表 | Generate visualization plots"""
#         self.logger.info("📊 生成可视化图表 | Generating visualization plots...")
        
#         if not model_results:
#             self.logger.warning("⚠️ 无结果数据，跳过可视化 | No result data, skipping visualization")
#             return
        
#         if self.config.genetic_model == "all":
#             # 为每个模型生成图表 | Generate plots for each model
#             for model_name, results_df in model_results.items():
#                 if len(results_df) > 0:
#                     self._generate_model_plots(results_df, model_name)
            
#             # 生成模型比较图表 | Generate model comparison plots
#             self._generate_comparison_plots(model_results)
#         else:
#             # 为单一模型生成图表 | Generate plots for single model
#             if model_results and len(model_results) > 0:
#                 results_df = list(model_results.values())[0]
#                 if len(results_df) > 0:
#                     self._generate_model_plots(results_df, self.config.genetic_model)
    
#     def _generate_model_plots(self, results_df: pd.DataFrame, model_name: str):
#         """为特定模型生成图表 | Generate plots for specific model"""
#         try:
#             self._create_manhattan_plot(results_df, model_name)
#             self._create_qq_plot(results_df, model_name)
#             self.logger.info(f"📈 {model_name}模型图表生成完成 | {model_name} model plots generated")
#         except Exception as e:
#             self.logger.warning(f"⚠️ {model_name}模型图表生成失败 | Plot generation failed for {model_name}: {e}")
    
#     def _create_manhattan_plot(self, results_df: pd.DataFrame, model_name: str):
#         """创建Manhattan图 | Create Manhattan plot"""
#         model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
        
#         r_script = f"""
# library(ggplot2)
# library(dplyr)

# # 读取数据
# tryCatch({{
#     data <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
    
#     # 检查数据
#     if(nrow(data) == 0) {{
#         cat("No data available for Manhattan plot\\n")
#         quit()
#     }}
    
#     # 计算-log10(P)
#     data$log_p <- -log10(data$P)
    
#     # 处理染色体信息
#     data$CHR <- as.factor(data$CHR)
    
#     # 创建Manhattan图
#     png("manhattan_plot_{model_name}.png", width=1200, height=800, res=150)
#     p <- ggplot(data, aes(x=BP, y=log_p, color=CHR)) +
#         geom_point(alpha=0.6, size=0.8) +
#         theme_minimal() +
#         labs(title="Manhattan Plot - {model_name.upper()} Model",
#              x="Genomic Position", 
#              y="-log10(P-value)",
#              color="Chromosome") +
#         theme(legend.position="bottom",
#               axis.text.x = element_text(angle=45, hjust=1))
    
#     # 添加显著性阈值线
#     if(max(data$log_p, na.rm=TRUE) > 5) {{
#         p <- p + geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed", alpha=0.7)
#     }}
#     if(max(data$log_p, na.rm=TRUE) > 3) {{
#         p <- p + geom_hline(yintercept=-log10(1e-5), color="blue", linetype="dashed", alpha=0.7)
#     }}
    
#     print(p)
#     dev.off()
    
#     cat("Manhattan plot created successfully\\n")
# }}, error = function(e) {{
#     cat("Error creating Manhattan plot:", conditionMessage(e), "\\n")
# }})
# """
        
#         with open(f"manhattan_{model_name}.R", 'w') as f:
#             f.write(r_script)
    
#     def _create_qq_plot(self, results_df: pd.DataFrame, model_name: str):
#         """创建QQ图 | Create QQ plot"""
#         model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
        
#         r_script = f"""
# library(ggplot2)

# # 读取数据
# tryCatch({{
#     data <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
    
#     # 检查数据
#     if(nrow(data) == 0) {{
#         cat("No data available for QQ plot\\n")
#         quit()
#     }}
    
#     # 移除无效P值
#     data <- data[!is.na(data$P) & data$P > 0 & data$P <= 1, ]
    
#     if(nrow(data) == 0) {{
#         cat("No valid P-values for QQ plot\\n")
#         quit()
#     }}
    
#     # 计算QQ图数据
#     observed <- sort(-log10(data$P))
#     expected <- sort(-log10(ppoints(length(observed))))
    
#     qq_data <- data.frame(expected=expected, observed=observed)
    
#     # 创建QQ图
#     png("qq_plot_{model_name}.png", width=800, height=800, res=150)
#     p <- ggplot(qq_data, aes(x=expected, y=observed)) +
#         geom_point(alpha=0.6) +
#         geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
#         theme_minimal() +
#         labs(title="Q-Q Plot - {model_name.upper()} Model",
#              x="Expected -log10(P)",
#              y="Observed -log10(P)")
    
#     print(p)
#     dev.off()
    
#     cat("QQ plot created successfully\\n")
# }}, error = function(e) {{
#     cat("Error creating QQ plot:", conditionMessage(e), "\\n")
# }})
# """
        
#         with open(f"qq_{model_name}.R", 'w') as f:
#             f.write(r_script)
    
#     def _generate_comparison_plots(self, model_results: dict):
#         """生成模型比较图表 | Generate model comparison plots"""
#         self.logger.info("📊 生成模型比较图表 | Generating model comparison plots...")
        
#         # 生成模型比较的R脚本
#         r_script = """
# library(ggplot2)
# library(dplyr)
# library(reshape2)

# # 比较不同模型的P值分布
# models <- c()
# p_values <- c()

# """
        
#         for model_name in model_results.keys():
#             model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
#             r_script += f"""
# if(file.exists("gwas_results_{model_suffix}.txt")) {{
#     tryCatch({{
#         data_{model_name} <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
#         if(nrow(data_{model_name}) > 0) {{
#             models <- c(models, rep("{model_name}", nrow(data_{model_name})))
#             p_values <- c(p_values, data_{model_name}$P)
#         }}
#     }}, error = function(e) {{
#         cat("Error reading {model_suffix} data:", conditionMessage(e), "\\n")
#     }})
# }}
# """
        
#         r_script += """
# # 检查是否有数据
# if(length(models) == 0 || length(p_values) == 0) {
#     cat("No data available for comparison plots\\n")
#     quit()
# }

# # 创建比较数据框
# comparison_data <- data.frame(Model=models, P_value=p_values)
# comparison_data <- comparison_data[!is.na(comparison_data$P_value) & 
#                                   comparison_data$P_value > 0 & 
#                                   comparison_data$P_value <= 1, ]

# if(nrow(comparison_data) == 0) {
#     cat("No valid data for comparison plots\\n")
#     quit()
# }

# comparison_data$log_p <- -log10(comparison_data$P_value)

# # 生成P值分布比较图
# tryCatch({
#     png("model_comparison_pvalue_distribution.png", width=1200, height=800, res=150)
#     p <- ggplot(comparison_data, aes(x=log_p, fill=Model)) +
#         geom_histogram(alpha=0.7, bins=50, position="identity") +
#         facet_wrap(~Model, scales="free_y") +
#         theme_minimal() +
#         labs(title="P-value Distribution Comparison Across Models",
#              x="-log10(P-value)",
#              y="Frequency") +
#         theme(legend.position="bottom")
#     print(p)
#     dev.off()
    
#     cat("Model comparison plot created successfully\\n")
# }, error = function(e) {
#     cat("Error creating comparison plot:", conditionMessage(e), "\\n")
# })
# """
        
#         with open("model_comparison.R", 'w') as f:
#             f.write(r_script)
        
#         self.logger.info("📄 可视化R脚本已生成 | Visualization R scripts generated")

# Google AI version
# results.py - Part 1/3

import pandas as pd
import numpy as np
import os
from pathlib import Path
from datetime import datetime

# 使用 statsmodels 进行 FDR 校正，这是一个更标准、功能更强大的库
try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False

class ResultsProcessor:
    """结果处理器 | Results Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def process_results(self, main_result: str) -> dict:
        """
        处理分析结果的主入口。
        Main entry point for processing analysis results.
        """
        self.logger.info("🔄 处理分析结果 | Processing analysis results...")
        self.logger.info(f"📁 主结果文件前缀 | Main result file prefix: {main_result}")

        if self.config.genetic_model == "all":
            return self._process_all_models_results(main_result)
        else:
            return self._process_single_model_results(main_result)

    def _process_all_models_results(self, main_result: str) -> dict:
        """处理所有遗传模型（--model）的结果 | Process results for all genetic models (--model)"""
        self.logger.info("📊 处理所有遗传模型结果 | Processing all genetic models results...")

        result_file = f"{main_result}.model"
        if not os.path.exists(result_file):
            self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
            raise FileNotFoundError(f"结果文件不存在: {result_file}")

        try:
            results_df = pd.read_csv(result_file, sep=r"\s+")
            self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
        except Exception as e:
            self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
            return {}

        # 健壮地创建 BP 列
        results_df = self._ensure_bp_column(results_df)

        model_results = {}
        model_map = {
            'ADD': 'additive',
            'DOMDEV': 'dominant',
            'REC': 'recessive',
            'GENO_2DF': 'genotypic'
        }

        for test_name, model_name in model_map.items():
            if test_name in results_df['TEST'].values:
                model_df = results_df[results_df['TEST'] == test_name].copy()
                model_df = self._clean_results(model_df)
                if not model_df.empty:
                    file_suffix = model_name[:3].upper()
                    model_df.to_csv(f"gwas_results_{file_suffix}.txt", sep='\t', index=False)
                    self._apply_significance_correction(model_df, model_name)
                    model_results[model_name] = model_df
                    self.logger.info(f"✅ {model_name.capitalize()}模型结果: {len(model_df)} SNPs")

        if model_results:
            self._generate_model_comparison_report(model_results)

        return model_results

    def _process_single_model_results(self, main_result: str) -> dict:
        """处理单一遗传模型（--logistic/--linear）的结果 | Process results for a single genetic model"""
        self.logger.info(f"📊 处理{self.config.genetic_model}模型结果 | Processing {self.config.genetic_model} model results")

        suffix = "assoc.logistic" if self.config.trait_type == "qualitative" else "assoc.linear"
        result_file = f"{main_result}.{suffix}"
        
        if not os.path.exists(result_file):
            self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
            raise FileNotFoundError(f"结果文件不存在: {result_file}")

        try:
            results_df = pd.read_csv(result_file, sep=r"\s+")
            self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
        except Exception as e:
            self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
            return {}
        
        results_df = self._ensure_bp_column(results_df)
        model_results_df = self._clean_results(results_df)

        if not model_results_df.empty:
            file_suffix = self.config.genetic_model[:3].upper()
            model_results_df.to_csv(f"gwas_results_{file_suffix}.txt", sep='\t', index=False)
            self._apply_significance_correction(model_results_df, self.config.genetic_model)
            return {self.config.genetic_model: model_results_df}
        else:
            self.logger.warning(f"⚠️ {self.config.genetic_model}模型无有效结果 | No valid results for {self.config.genetic_model} model")
            return {}

    def _ensure_bp_column(self, df: pd.DataFrame) -> pd.DataFrame:
        """确保DataFrame中存在BP列，如果不存在则尝试从SNP列创建"""
        if 'BP' in df.columns:
            return df

        self.logger.warning("⚠️ 结果文件中缺少BP列，尝试从SNP名称提取位置信息 (格式: CHR:BP)")
        if not df.empty and df['SNP'].astype(str).str.contains(':').any():
            try:
                # n=1 确保只按第一个冒号分割，以处理如 "chr1:12345:A:T" 的情况
                bp_series = df['SNP'].str.split(':', n=1).str[1]
                # 再次分割以处理 "12345:A:T" 格式
                bp_series = bp_series.str.split(':').str[0]
                df['BP'] = pd.to_numeric(bp_series, errors='coerce')
                df['BP'] = df['BP'].fillna(0).astype(int)
                self.logger.info("✅ 成功从SNP名称提取BP位置信息")
            except Exception as e:
                self.logger.warning(f"⚠️ 无法从SNP名称提取BP信息: {e}。BP列将填充为0。")
                df['BP'] = 0
        else:
            self.logger.warning("⚠️ SNP列不包含':'分隔符，无法提取位置。BP列将填充为0。")
            df['BP'] = 0
        return df

    def _clean_results(self, df: pd.DataFrame) -> pd.DataFrame:
        """清理结果数据，移除无效的P值等"""
        if df.empty:
            return df
        
        original_count = len(df)
        # 确保P值为数字类型并移除NA
        df['P'] = pd.to_numeric(df['P'], errors='coerce')
        df.dropna(subset=['P'], inplace=True)
        # 筛选有效的P值范围
        df = df[(df['P'] > 0) & (df['P'] <= 1)]
        
        cleaned_count = len(df)
        removed_count = original_count - cleaned_count
        if removed_count > 0:
            self.logger.info(f"🧹 数据清理：移除了{removed_count}个无效记录 | Data cleaning: removed {removed_count} invalid records")
        
        return df

    def _apply_significance_correction(self, df: pd.DataFrame, model_name: str):
        """应用所有选择的显著性校正方法"""
        self.logger.info(f"📈 对{model_name}模型应用显著性校正: {self.config.correction_method} | Applying significance corrections for {model_name} model: {self.config.correction_method}")
        if df.empty:
            return

        methods = [self.config.correction_method] if self.config.correction_method != "all" else ["bonferroni", "suggestive", "fdr"]
        
        if "bonferroni" in methods:
            self._bonferroni_correction(df, model_name)
        if "suggestive" in methods:
            self._suggestive_correction(df, model_name)
        if "fdr" in methods:
            self._fdr_correction(df, model_name)

    def _bonferroni_correction(self, df: pd.DataFrame, model_name: str):
        total_snps = len(df)
        if total_snps == 0: return
        
        bonferroni_threshold = self.config.bonferroni_alpha / total_snps
        significant = df[df['P'] < bonferroni_threshold].copy()
        
        self.logger.info(f"🎯 {model_name}模型Bonferroni校正阈值: {bonferroni_threshold:.2e} | Bonferroni threshold: {bonferroni_threshold:.2e}")
        self.logger.info(f"⭐ {model_name}模型Bonferroni显著位点数: {len(significant)} | Bonferroni significant loci: {len(significant)}")
        significant.to_csv(f"significant_hits_bonferroni_{model_name}.txt", sep='\t', index=False)

    def _suggestive_correction(self, df: pd.DataFrame, model_name: str):
        suggestive = df[df['P'] < self.config.suggestive_threshold].copy()
        
        self.logger.info(f"🎯 {model_name}模型提示性关联阈值: {self.config.suggestive_threshold:.2e} | Suggestive threshold: {self.config.suggestive_threshold:.2e}")
        self.logger.info(f"🔍 {model_name}模型提示性关联位点数: {len(suggestive)} | Suggestive loci: {len(suggestive)}")
        suggestive.to_csv(f"suggestive_hits_{model_name}.txt", sep='\t', index=False)

    def _fdr_correction(self, df: pd.DataFrame, model_name: str):
        """使用 Benjamini/Hochberg 方法进行FDR校正"""
        if not HAS_STATSMODELS:
            self.logger.warning("⚠️ 未安装statsmodels库，跳过FDR校正 | statsmodels library not installed, skipping FDR correction.")
            return

        try:
            p_values = df['P'].dropna().values
            if len(p_values) == 0:
                self.logger.warning(f"⚠️ {model_name}模型没有有效的P值进行FDR校正 | No valid p-values for FDR correction in {model_name} model")
                return

            is_significant, p_corrected, _, _ = multipletests(
                p_values, 
                alpha=self.config.fdr_alpha, 
                method='fdr_bh'
            )
            
            df_fdr = df.dropna(subset=['P']).copy()
            df_fdr['is_significant_fdr'] = is_significant
            df_fdr['q_value'] = p_corrected
            
            fdr_significant = df_fdr[df_fdr['is_significant_fdr']].copy()
            
            self.logger.info(f"🎯 {model_name}模型FDR校正q值阈值: {self.config.fdr_alpha} | FDR q-value threshold: {self.config.fdr_alpha}")
            self.logger.info(f"📊 {model_name}模型FDR显著位点数: {len(fdr_significant)} | FDR significant loci: {len(fdr_significant)}")
            
            df_fdr.to_csv(f"gwas_results_fdr_{model_name}.txt", sep='\t', index=False)
            fdr_significant.to_csv(f"fdr_significant_hits_{model_name}.txt", sep='\t', index=False)
            
            if not fdr_significant.empty:
                top_hits = fdr_significant.nsmallest(min(5, len(fdr_significant)), 'P')
                self.logger.info(f"📊 {model_name}模型前{len(top_hits)}个FDR显著位点 | Top {len(top_hits)} FDR significant loci:")
                for _, row in top_hits.iterrows():
                    bp_info = f", BP={row.get('BP', 'N/A')}"
                    self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, q-value={row['q_value']:.2e}, CHR={row['CHR']}{bp_info}")
        except Exception as e:
            self.logger.error(f"❌ FDR校正失败 | FDR correction failed: {e}", exc_info=True)

    def _generate_model_comparison_report(self, model_results: dict):
        """生成多模型比较报告"""
        self.logger.info("📋 生成模型比较报告 | Generating model comparison report...")
        
        report = []
        report.append("=" * 80)
        report.append("PLINK GWAS 多模型比较报告 | Multi-Model Comparison Report")
        report.append("=" * 80)
        report.append(f"分析时间 | Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"表型类型 | Trait Type: {self.config.trait_type}")
        
        for model_name, df in model_results.items():
            report.append(f"\n--- {model_name.upper()} 模型结果 ---")
            if df.empty:
                report.append("  ❌ 无有效结果 | No valid results")
                continue
            
            report.append(f"  分析SNP数 | SNPs Analyzed: {len(df)}")
            report.append(f"  最小P值 | Minimum P-value: {df['P'].min():.2e}")
            
            top_snps = df.nsmallest(5, 'P')
            if not top_snps.empty:
                report.append("  前5个最显著SNP | Top 5 Most Significant SNPs:")
                for _, row in top_snps.iterrows():
                    or_val = row.get('OR')
                    or_info = f", OR={or_val:.3f}" if pd.notna(or_val) else ""
                    bp_info = f", BP={row.get('BP', 'N/A')}"
                    report.append(f"    - {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}{bp_info}{or_info}")
        
        report_content = "\n".join(report)
        with open("model_comparison_report.txt", 'w', encoding='utf-8') as f:
            f.write(report_content)
        self.logger.info("📄 模型比较报告已保存 | Model comparison report saved: model_comparison_report.txt")

# --- End of Part 1/3 ---

# results.py - Part 2/3

class ReportGenerator:
    """报告生成器 | Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary_report(self, model_results: dict, stats: dict):
        """生成最终的总结报告"""
        self.logger.info("📋 生成总结报告 | Generating summary report...")
        
        report = []
        report.append("=" * 80)
        report.append("PLINK GWAS分析总结报告 | PLINK GWAS Analysis Summary Report")
        report.append("=" * 80)
        report.append(f"完成时间 | Completion Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"输出目录 | Output Directory: {Path.cwd().resolve()}")
        report.append("\n--- 分析配置 | Analysis Configuration ---")
        report.append(f"遗传模型 | Genetic Model: {self.config.genetic_model}")
        report.append(f"表型类型 | Trait Type: {self.config.trait_type}")
        report.append(f"显著性校正 | Significance Correction: {self.config.correction_method}")
        
        report.append("\n--- 数据统计 | Data Statistics ---")
        report.append(f"质控后样本数 | Samples after QC: {stats.get('total_samples', 'N/A')}")
        report.append(f"质控后SNP数 | SNPs after QC: {stats.get('total_snps', 'N/A')}")
        if self.config.trait_type == "qualitative":
            report.append(f"病例数 (Cases): {stats.get('cases', 'N/A')}")
            report.append(f"对照数 (Controls): {stats.get('controls', 'N/A')}")
            
        report.append("\n--- 关联分析结果摘要 | Association Results Summary ---")
        if not model_results:
            report.append("  ❌ 未生成任何有效结果 | No valid results were generated.")
        else:
            for model_name, df in model_results.items():
                report.append(f"  模型 | Model: {model_name.upper()}")
                if df.empty:
                    report.append("    - 无有效结果 | No valid results")
                    continue
                
                report.append(f"    - 分析SNP数 | SNPs Analyzed: {len(df)}")
                report.append(f"    - 最小P值 | Min P-value: {df['P'].min():.2e}")
                
                top_snp = df.loc[df['P'].idxmin()]
                bp_info = f":{top_snp.get('BP', 'N/A')}"
                report.append(f"    - 最显著SNP | Top SNP: {top_snp['SNP']} (CHR{top_snp['CHR']}{bp_info})")

        report_content = "\n".join(report)
        with open("gwas_summary_report.txt", 'w', encoding='utf-8') as f:
            f.write(report_content)
        self.logger.info("📄 GWAS总结报告已保存 | GWAS summary report saved: gwas_summary_report.txt")

# --- End of Part 2/3 ---

# results.py - Part 3/3

# class VisualizationGenerator:
#     """可视化生成器 | Visualization Generator"""

#     def __init__(self, config, logger, cmd_runner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner

#     def generate_plots(self, model_results: dict):
#         """为所有有效结果生成图表"""
#         self.logger.info("📊 生成可视化图表 | Generating visualization plots...")
#         if not model_results:
#             self.logger.warning("⚠️ 无结果数据，跳过可视化 | No results data, skipping visualization.")
#             return

#         for model_name, df in model_results.items():
#             if not df.empty:
#                 self.logger.info(f"--- 正在为 {model_name} 模型生成图表 | Generating plots for {model_name} model ---")
#                 self._generate_model_plots(df, model_name)
        
#         if len(model_results) > 1:
#              self._generate_comparison_plots(model_results)

#     def _generate_model_plots(self, df: pd.DataFrame, model_name: str):
#         """为特定模型生成Manhattan图和QQ图"""
#         try:
#             # 确保数据已保存，以便R脚本可以读取
#             file_suffix = model_name[:3].upper()
#             results_file = f"gwas_results_{file_suffix}.txt"
#             if not Path(results_file).exists():
#                 df.to_csv(results_file, sep='\t', index=False)

#             self._create_manhattan_plot(results_file, model_name)
#             self._create_qq_plot(results_file, model_name)
#         except Exception as e:
#             self.logger.warning(f"⚠️ {model_name}模型图表生成失败 | Plot generation failed for {model_name} model: {e}")

#     def _run_r_script(self, script_name: str, script_content: str):
#         """执行R脚本"""
#         try:
#             with open(script_name, 'w', encoding='utf-8') as f:
#                 f.write(script_content)
            
#             cmd = ["Rscript", script_name]
#             result = self.cmd_runner.run(cmd, f"执行R脚本 | Executing R script: {script_name}", check=False)
#             if result.returncode != 0:
#                 self.logger.error(f"❌ R脚本 {script_name} 执行失败。请检查R环境和脚本内容。")
#                 self.logger.error(f"R脚本错误输出:\n{result.stderr}")
#         except FileNotFoundError:
#             self.logger.error("❌ Rscript 命令未找到。请确保R已安装并配置在系统PATH中。 | Rscript command not found. Please ensure R is installed and in the system's PATH.")
#         except Exception as e:
#             self.logger.error(f"❌ 执行R脚本时发生意外错误: {e}")

#     def _create_manhattan_plot(self, results_file: str, model_name: str):
#         r_script = f"""
#         # 加载必要的库，如果不存在则自动安装
#         if (!require(ggplot2)) install.packages("ggplot2", repos="http://cran.us.r-project.org")
#         if (!require(dplyr)) install.packages("dplyr", repos="http://cran.us.r-project.org")
#         library(ggplot2)
#         library(dplyr)
        
#         tryCatch({{
#             data <- read.table("{results_file}", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
#             if(nrow(data) == 0) quit()

#             data <- data %>% filter(!is.na(P) & P > 0 & P <= 1)
#             if(nrow(data) == 0) quit()

#             data$log_p <- -log10(data$P)
            
#             # 兼容数字和字符型染色体名称
#             # 提取染色体中的数字部分用于排序
#             data$CHR_num <- as.numeric(gsub("[^0-9]+", "", data$CHR))
#             data <- data[order(data$CHR_num, data$BP), ]
#             # 将原始CHR列作为因子，以保持其原始顺序（例如 'chr1', 'chr2'）
#             data$CHR <- factor(data$CHR, levels=unique(data$CHR))
            
#             # 计算累积BP位置
#             data_cum <- data %>% 
#               group_by(CHR) %>% 
#               summarise(max_bp = max(as.numeric(BP), na.rm=TRUE)) %>% 
#               mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default=0)) %>% 
#               select(CHR, bp_add)
            
#             data <- data %>% inner_join(data_cum, by="CHR") %>%
#               mutate(bp_cum = as.numeric(BP) + bp_add)
              
#             axis_set <- data %>% group_by(CHR) %>% summarize(center = mean(bp_cum))
            
#             p <- ggplot(data, aes(x=bp_cum, y=log_p, color=CHR)) +
#               geom_point(alpha=0.7, size=1.3) +
#               scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
#               scale_color_manual(values = rep(c("#276FBF", "#183059"), length.out=nlevels(data$CHR))) +
#               geom_hline(yintercept=-log10(1e-5), color="blue", linetype="dashed", size=0.8) +
#               geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed", size=0.8) +
#               labs(x="Chromosome", y="-log10(P-value)", title=paste("Manhattan Plot -", toupper("{model_name}"), "Model")) +
#               theme_minimal(base_size=16) +
#               theme(
#                 legend.position="none", 
#                 panel.grid.major.x=element_blank(), 
#                 panel.grid.minor.x=element_blank(),
#                 axis.text.x = element_text(angle=60, vjust=0.5, size=10)
#               )
            
#             ggsave("manhattan_plot_{model_name}.png", plot=p, width=14, height=7, dpi=300)
            
#         }}, error=function(e) {{
#             cat("Error in Manhattan plot generation for {model_name}:", conditionMessage(e), "\\n")
#             # 创建一个空白的错误图片
#             png("manhattan_plot_{model_name}.png", width=1400, height=700)
#             plot(1, type="n", xlab="", ylab="", main="Manhattan Plot Generation Failed")
#             text(1, 1, conditionMessage(e), cex=1.2)
#             dev.off()
#         }})
#         """
#         self._run_r_script(f"manhattan_{model_name}.R", r_script)

#     def _create_qq_plot(self, results_file: str, model_name: str):
#         r_script = f"""
#         # 加载必要的库
#         if (!require(ggplot2)) install.packages("ggplot2", repos="http://cran.us.r-project.org")
#         library(ggplot2)
        
#         tryCatch({{
#             data <- read.table("{results_file}", header=TRUE, sep="\\t", comment.char="", check.names=FALSE)
#             if(nrow(data) == 0) quit()

#             data <- data[!is.na(data$P) & data$P > 0 & data$P <= 1, ]
#             if(nrow(data) == 0) quit()

#             observed <- sort(-log10(data$P))
#             expected <- sort(-log10(ppoints(length(observed))))
            
#             # 计算 lambda (基因组膨胀因子)
#             chisq <- qchisq(1 - data$P, 1)
#             lambda = median(chisq, na.rm=TRUE) / qchisq(0.5, 1)
            
#             qq_df <- data.frame(observed=observed, expected=expected)
            
#             p <- ggplot(qq_df, aes(x=expected, y=observed)) +
#               geom_point(alpha=0.5, color="black") +
#               geom_abline(intercept=0, slope=1, color="red", linetype="dashed", size=1) +
#               labs(
#                 title=paste("Q-Q Plot -", toupper("{model_name}"), "Model"),
#                 subtitle=paste("Lambda (λ) =", format(lambda, digits=4)),
#                 x="Expected -log10(P)", 
#                 y="Observed -log10(P)"
#               ) +
#               theme_minimal(base_size=16) +
#               coord_fixed()
            
#             ggsave("qq_plot_{model_name}.png", plot=p, width=7, height=7, dpi=300)
        
#         }}, error=function(e) {{
#             cat("Error in QQ plot generation for {model_name}:", conditionMessage(e), "\\n")
#             # 创建一个空白的错误图片
#             png("qq_plot_{model_name}.png", width=700, height=700)
#             plot(1, type="n", xlab="", ylab="", main="QQ Plot Generation Failed")
#             text(1, 1, conditionMessage(e), cex=1.2)
#             dev.off()
#         }})
#         """
#         self._run_r_script(f"qq_{model_name}.R", r_script)

#     def _generate_comparison_plots(self, model_results: dict):
#         """此功能当前未激活，但保留框架以便未来扩展"""
#         self.logger.info("模型比较图表功能暂未实现 | Model comparison plot feature is not implemented yet.")
#         pass

# --- End of Part 3/3 - End of file ---

# python版本的可视化
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

class VisualizationGenerator:
    """可视化生成器 | Visualization Generator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 设置matplotlib中文支持和样式 | Setup matplotlib Chinese support and styles
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
        plt.rcParams['axes.unicode_minus'] = False
        sns.set_style("whitegrid")

    def generate_plots(self, model_results: dict):
        """为所有有效结果生成图表 | Generate plots for all valid results"""
        self.logger.info("📊 生成可视化图表 | Generating visualization plots...")
        if not model_results:
            self.logger.warning("⚠️ 无结果数据，跳过可视化 | No results data, skipping visualization.")
            return

        for model_name, df in model_results.items():
            if not df.empty:
                self.logger.info(f"--- 正在为 {model_name} 模型生成图表 | Generating plots for {model_name} model ---")
                self._generate_model_plots(df, model_name)
        
        if len(model_results) > 1:
            self._generate_comparison_plots(model_results)

    def _generate_model_plots(self, df: pd.DataFrame, model_name: str):
        """为特定模型生成Manhattan图和QQ图 | Generate Manhattan and QQ plots for specific model"""
        try:
            # 确保数据已保存，以便后续处理 | Ensure data is saved for subsequent processing
            file_suffix = model_name[:3].upper()
            results_file = f"gwas_results_{file_suffix}.txt"
            if not Path(results_file).exists():
                df.to_csv(results_file, sep='\t', index=False)

            self._create_manhattan_plot(df, model_name)
            self._create_qq_plot(df, model_name)
        except Exception as e:
            self.logger.warning(f"⚠️ {model_name}模型图表生成失败 | Plot generation failed for {model_name} model: {e}")

    def _create_manhattan_plot(self, df: pd.DataFrame, model_name: str):
        """创建Manhattan图 | Create Manhattan plot"""
        try:
            # 数据预处理 | Data preprocessing
            data = df.copy()
            
            # 过滤无效的P值 | Filter invalid P-values
            data = data.dropna(subset=['P'])
            data = data[(data['P'] > 0) & (data['P'] <= 1)]
            
            if len(data) == 0:
                self.logger.warning(f"⚠️ {model_name}模型无有效P值数据 | No valid P-value data for {model_name} model")
                return
            
            # 计算-log10(P值) | Calculate -log10(P-values)
            data['log_p'] = -np.log10(data['P'])
            
            # 处理染色体信息 | Process chromosome information
            data['CHR_str'] = data['CHR'].astype(str)
            # 提取数字部分用于排序 | Extract numeric part for sorting
            data['CHR_num'] = data['CHR_str'].str.extract('(\d+)').astype(float)
            data = data.dropna(subset=['CHR_num'])
            data = data.sort_values(['CHR_num', 'BP'])
            
            # 计算累积位置 | Calculate cumulative positions
            chr_lengths = data.groupby('CHR_str')['BP'].max().reset_index()
            chr_lengths['chr_start'] = chr_lengths['BP'].cumsum().shift(1).fillna(0)
            
            # 合并累积位置信息 | Merge cumulative position information
            data = data.merge(chr_lengths[['CHR_str', 'chr_start']], 
                            left_on='CHR_str', right_on='CHR_str')
            data['bp_cum'] = data['BP'] + data['chr_start']
            
            # 计算染色体中心位置用于x轴标签 | Calculate chromosome center positions for x-axis labels
            chr_centers = data.groupby('CHR_str')['bp_cum'].mean().reset_index()
            
            # 创建图形 | Create figure
            fig, ax = plt.subplots(figsize=(14, 7))
            
            # 为每个染色体分配颜色 | Assign colors for each chromosome
            unique_chrs = data['CHR_str'].unique()
            colors = ['#1F77B4FF', '#FF7F0EFF'] * (len(unique_chrs) // 2 + 1)
            
            # 绘制散点图 | Plot scatter points
            for i, chr_name in enumerate(unique_chrs):
                chr_data = data[data['CHR_str'] == chr_name]
                ax.scatter(chr_data['bp_cum'], chr_data['log_p'], 
                          c=colors[i], alpha=0.7, s=8)
            
            # 计算Y轴最大值 | Calculate Y-axis maximum value
            max_log_p = data['log_p'].max()
            
            # 添加显著性水平线 | Add significance threshold lines
            # 从5开始，隔1个数字画一条横线，直到最大值下面那个值
            # Draw lines starting from 5, every 1 unit, up to the value below maximum
            significance_colors = ['blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive']
            line_start = 5
            max_line_value = int(np.floor(max_log_p))
            
            for i, line_value in enumerate(range(line_start, max_line_value + 1)):
                color = significance_colors[i % len(significance_colors)]
                ax.axhline(y=line_value, color=color, linestyle='--', 
                          linewidth=1, alpha=0.8)
            
            # 设置x轴标签 | Set x-axis labels
            ax.set_xticks(chr_centers['bp_cum'])
            ax.set_xticklabels(chr_centers['CHR_str'], rotation=60)
            
            # 设置Y轴 | Set Y-axis
            ax.set_ylim(0, max_log_p * 1.05)  # Y轴从0开始 | Y-axis starts from 0
            
            # 设置Y轴刻度 | Set Y-axis ticks
            if max_log_p <= 10:
                # 如果最大值不超过10，从0到最大值，间隔1显示 | If max <= 10, show from 0 to max with interval 1
                y_ticks = np.arange(0, int(np.ceil(max_log_p)) + 1, 1)
                ax.set_yticks(y_ticks)
            else:
                # 如果超过10，灵活处理 | If > 10, flexible handling
                # 可以根据最大值调整间隔 | Adjust interval based on maximum value
                if max_log_p <= 20:
                    interval = 2
                elif max_log_p <= 50:
                    interval = 5
                else:
                    interval = 10
                y_ticks = np.arange(0, int(np.ceil(max_log_p)) + 1, interval)
                ax.set_yticks(y_ticks)
            
            # 设置标签和标题 | Set labels and title
            ax.set_xlabel('Chromosome', fontsize=12)
            ax.set_ylabel('-log10(P-value)', fontsize=12)
            ax.set_title(f'Manhattan Plot - {model_name.upper()} Model', fontsize=14, fontweight='bold')
            
            # 优化布局 | Optimize layout
            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            
            # 保存图片 | Save figure
            output_file = f"manhattan_plot_{model_name}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"✅ Manhattan图已保存 | Manhattan plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ Manhattan图生成失败 | Manhattan plot generation failed: {e}")
            # 创建错误占位图 | Create error placeholder image
            fig, ax = plt.subplots(figsize=(14, 7))
            ax.text(0.5, 0.5, f'Manhattan Plot Generation Failed\n{str(e)}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title(f'Manhattan Plot - {model_name.upper()} Model (Error)', fontsize=14)
            plt.savefig(f"manhattan_plot_{model_name}.png", dpi=300, bbox_inches='tight')
            plt.close()

    def _create_qq_plot(self, df: pd.DataFrame, model_name: str):
        """创建QQ图 | Create QQ plot"""
        try:
            # 数据预处理 | Data preprocessing
            data = df.copy()
            
            # 过滤无效的P值 | Filter invalid P-values
            data = data.dropna(subset=['P'])
            data = data[(data['P'] > 0) & (data['P'] <= 1)]
            
            if len(data) == 0:
                self.logger.warning(f"⚠️ {model_name}模型无有效P值数据 | No valid P-value data for {model_name} model")
                return
            
            # 计算观察值和期望值 | Calculate observed and expected values
            observed = np.sort(-np.log10(data['P']))
            n = len(observed)
            expected = np.sort(-np.log10(np.arange(1, n + 1) / (n + 1)))
            
            # 计算基因组膨胀因子 lambda | Calculate genomic inflation factor lambda
            chisq = stats.chi2.ppf(1 - data['P'], df=1)
            lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
            
            # 创建图形 | Create figure
            fig, ax = plt.subplots(figsize=(7, 7))
            
            # 绘制QQ图 | Plot QQ plot
            ax.scatter(expected, observed, alpha=0.6, s=20, color='black', edgecolors='none')
            
            # 添加对角线 | Add diagonal line
            max_val = max(max(expected), max(observed))
            ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, alpha=0.8)
            
            # 设置标签和标题 | Set labels and title
            ax.set_xlabel('Expected -log10(P)', fontsize=12)
            ax.set_ylabel('Observed -log10(P)', fontsize=12)
            ax.set_title(f'Q-Q Plot - {model_name.upper()} Model', fontsize=14, fontweight='bold')
            
            # 添加lambda值 | Add lambda value
            ax.text(0.05, 0.95, f'λ = {lambda_gc:.4f}', 
                   transform=ax.transAxes, fontsize=11, 
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))
            
            # 设置相等的纵横比 | Set equal aspect ratio
            ax.set_aspect('equal', adjustable='box')
            ax.grid(True, alpha=0.3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            plt.tight_layout()
            
            # 保存图片 | Save figure
            output_file = f"qq_plot_{model_name}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"✅ QQ图已保存 | QQ plot saved: {output_file} (λ = {lambda_gc:.4f})")
            
        except Exception as e:
            self.logger.error(f"❌ QQ图生成失败 | QQ plot generation failed: {e}")
            # 创建错误占位图 | Create error placeholder image
            fig, ax = plt.subplots(figsize=(7, 7))
            ax.text(0.5, 0.5, f'QQ Plot Generation Failed\n{str(e)}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title(f'Q-Q Plot - {model_name.upper()} Model (Error)', fontsize=14)
            plt.savefig(f"qq_plot_{model_name}.png", dpi=300, bbox_inches='tight')
            plt.close()

    def _generate_comparison_plots(self, model_results: dict):
        """生成模型比较图表 | Generate model comparison plots"""
        try:
            self.logger.info("📈 生成模型比较图表 | Generating model comparison plots...")
            
            # P值分布比较 | P-value distribution comparison
            self._create_p_value_distribution_comparison(model_results)
            
            # 效应大小比较（如果有BETA列） | Effect size comparison (if BETA column exists)
            if all('BETA' in df.columns for df in model_results.values()):
                self._create_effect_size_comparison(model_results)
                
        except Exception as e:
            self.logger.error(f"❌ 模型比较图表生成失败 | Model comparison plot generation failed: {e}")

    def _create_p_value_distribution_comparison(self, model_results: dict):
        """创建P值分布比较图 | Create P-value distribution comparison plot"""
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # P值分布直方图 | P-value distribution histogram
            for model_name, df in model_results.items():
                if not df.empty and 'P' in df.columns:
                    valid_p = df['P'].dropna()
                    valid_p = valid_p[(valid_p > 0) & (valid_p <= 1)]
                    if len(valid_p) > 0:
                        ax1.hist(-np.log10(valid_p), bins=50, alpha=0.6, 
                                label=f'{model_name} (n={len(valid_p)})', density=True)
            
            ax1.set_xlabel('-log10(P-value)')
            ax1.set_ylabel('Density')
            ax1.set_title('P-value Distribution Comparison')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # 显著性SNP数量比较 | Significant SNPs count comparison
            significance_thresholds = [1e-5, 1e-6, 1e-7, 5e-8]
            model_names = list(model_results.keys())
            sig_counts = {threshold: [] for threshold in significance_thresholds}
            
            for model_name, df in model_results.items():
                if not df.empty and 'P' in df.columns:
                    valid_p = df['P'].dropna()
                    valid_p = valid_p[(valid_p > 0) & (valid_p <= 1)]
                    for threshold in significance_thresholds:
                        count = np.sum(valid_p < threshold)
                        sig_counts[threshold].append(count)
                else:
                    for threshold in significance_thresholds:
                        sig_counts[threshold].append(0)
            
            x = np.arange(len(model_names))
            width = 0.2
            
            for i, threshold in enumerate(significance_thresholds):
                ax2.bar(x + i * width, sig_counts[threshold], width, 
                       label=f'P < {threshold}', alpha=0.8)
            
            ax2.set_xlabel('Models')
            ax2.set_ylabel('Number of Significant SNPs')
            ax2.set_title('Significant SNPs Comparison')
            ax2.set_xticks(x + width * 1.5)
            ax2.set_xticklabels(model_names)
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig("model_comparison_p_values.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info("✅ P值比较图已保存 | P-value comparison plot saved: model_comparison_p_values.png")
            
        except Exception as e:
            self.logger.error(f"❌ P值比较图生成失败 | P-value comparison plot generation failed: {e}")

    def _create_effect_size_comparison(self, model_results: dict):
        """创建效应大小比较图 | Create effect size comparison plot"""
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            effect_sizes = []
            model_labels = []
            
            for model_name, df in model_results.items():
                if not df.empty and 'BETA' in df.columns:
                    valid_beta = df['BETA'].dropna()
                    if len(valid_beta) > 0:
                        effect_sizes.append(valid_beta)
                        model_labels.append(f'{model_name}\n(n={len(valid_beta)})')
            
            if effect_sizes:
                # 创建箱线图 | Create box plot
                bp = ax.boxplot(effect_sizes, labels=model_labels, patch_artist=True)
                
                # 设置颜色 | Set colors
                colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow']
                for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                ax.set_ylabel('Effect Size (BETA)')
                ax.set_title('Effect Size Distribution Comparison')
                ax.grid(True, alpha=0.3)
                ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
                
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig("model_comparison_effect_sizes.png", dpi=300, bbox_inches='tight')
                plt.close()
                
                self.logger.info("✅ 效应大小比较图已保存 | Effect size comparison plot saved: model_comparison_effect_sizes.png")
            else:
                self.logger.warning("⚠️ 无有效的效应大小数据用于比较 | No valid effect size data for comparison")
                
        except Exception as e:
            self.logger.error(f"❌ 效应大小比较图生成失败 | Effect size comparison plot generation failed: {e}")