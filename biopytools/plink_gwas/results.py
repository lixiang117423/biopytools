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

import pandas as pd
import numpy as np
import os
from pathlib import Path
from datetime import datetime

try:
    from scipy.stats import false_discovery_control
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_results(self, main_result: str) -> dict:
        """处理分析结果 | Process analysis results"""
        self.logger.info("🔄 处理分析结果 | Processing analysis results...")
        self.logger.info(f"📁 主结果文件前缀 | Main result file prefix: {main_result}")
        
        if self.config.genetic_model == "all":
            return self._process_all_models_results(main_result)
        else:
            return self._process_single_model_results(main_result)
    
    def _process_all_models_results(self, main_result: str) -> dict:
        """处理所有模型的结果 | Process all models results"""
        self.logger.info("📊 处理所有遗传模型结果 | Processing all genetic models results...")
        
        result_file = f"{main_result}.model"
        if not os.path.exists(result_file):
            self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
            raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
        # 检查文件大小 | Check file size
        file_size_mb = os.path.getsize(result_file) / (1024 * 1024)
        self.logger.info(f"📏 结果文件大小 | Result file size: {file_size_mb:.1f} MB")
        
        # 读取结果文件 | Read results file
        self.logger.info("📖 读取.model文件，这可能需要一些时间... | Reading .model file, this may take some time...")
        try:
            results_df = pd.read_csv(result_file, sep=r"\s+")
            self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
        except Exception as e:
            self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
            return {}
        
        # 检查数据结构 | Check data structure
        self.logger.info(f"📋 结果文件列名 | Result file columns: {list(results_df.columns)}")
        unique_tests = results_df['TEST'].unique()
        self.logger.info(f"🧪 包含的检验类型 | Test types included: {unique_tests}")
        
        # 分别处理各种模型的结果 | Process results for different models separately
        model_results = {}
        
        # 处理加性模型结果 | Process additive model results
        if 'ADD' in results_df['TEST'].values:
            add_results = results_df[results_df['TEST'] == 'ADD'].copy()
            add_results = self._clean_results(add_results)
            if len(add_results) > 0:
                add_results.to_csv("gwas_results_ADD.txt", sep='\t', index=False)
                self._apply_significance_correction(add_results, "ADD")
                model_results['additive'] = add_results
                self.logger.info(f"✅ 加性模型结果 | Additive model results: {len(add_results)} SNPs")
            else:
                self.logger.warning("⚠️ 加性模型无有效结果 | No valid additive model results")
        
        # 处理显性模型结果 | Process dominant model results
        if 'DOMDEV' in results_df['TEST'].values:
            dom_results = results_df[results_df['TEST'] == 'DOMDEV'].copy()
            dom_results = self._clean_results(dom_results)
            if len(dom_results) > 0:
                dom_results.to_csv("gwas_results_DOM.txt", sep='\t', index=False)
                self._apply_significance_correction(dom_results, "DOM")
                model_results['dominant'] = dom_results
                self.logger.info(f"✅ 显性模型结果 | Dominant model results: {len(dom_results)} SNPs")
            else:
                self.logger.warning("⚠️ 显性模型无有效结果 | No valid dominant model results")
        
        # 处理基因型模型结果 | Process genotypic model results
        if 'GENO_2DF' in results_df['TEST'].values:
            geno_results = results_df[results_df['TEST'] == 'GENO_2DF'].copy()
            geno_results = self._clean_results(geno_results)
            if len(geno_results) > 0:
                geno_results.to_csv("gwas_results_GENO.txt", sep='\t', index=False)
                self._apply_significance_correction(geno_results, "GENO")
                model_results['genotypic'] = geno_results
                self.logger.info(f"✅ 基因型模型结果 | Genotypic model results: {len(geno_results)} SNPs")
            else:
                self.logger.warning("⚠️ 基因型模型无有效结果 | No valid genotypic model results")
        
        # 处理隐性模型结果（如果存在）| Process recessive model results (if exists)
        if 'REC' in results_df['TEST'].values:
            rec_results = results_df[results_df['TEST'] == 'REC'].copy()
            rec_results = self._clean_results(rec_results)
            if len(rec_results) > 0:
                rec_results.to_csv("gwas_results_REC.txt", sep='\t', index=False)
                self._apply_significance_correction(rec_results, "REC")
                model_results['recessive'] = rec_results
                self.logger.info(f"✅ 隐性模型结果 | Recessive model results: {len(rec_results)} SNPs")
        
        # 生成综合比较报告 | Generate comprehensive comparison report
        if model_results:
            self._generate_model_comparison_report(model_results)
        else:
            self.logger.warning("⚠️ 未找到有效的模型结果 | No valid model results found")
        
        return model_results
    
    def _process_single_model_results(self, main_result: str) -> dict:
        """处理单一模型的结果 | Process single model results"""
        self.logger.info(f"📊 处理{self.config.genetic_model}模型结果 | Processing {self.config.genetic_model} model results...")
        
        # 根据表型类型确定结果文件 | Determine result file based on trait type
        if self.config.trait_type == "qualitative":
            result_file = f"{main_result}.assoc.logistic"
        else:
            result_file = f"{main_result}.assoc.linear"
        
        if not os.path.exists(result_file):
            self.logger.error(f"❌ 结果文件不存在 | Result file does not exist: {result_file}")
            raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
        # 读取结果文件 | Read results file
        try:
            results_df = pd.read_csv(result_file, sep=r"\s+")
            self.logger.info(f"✅ 成功读取结果文件 | Successfully read results file: {results_df.shape}")
        except Exception as e:
            self.logger.error(f"❌ 读取结果文件失败 | Failed to read results file: {e}")
            return {}
        
        # 根据遗传模型提取相应结果 | Extract results based on genetic model
        if self.config.genetic_model == "additive":
            model_results = results_df[results_df['TEST'] == 'ADD'].copy()
            file_suffix = "ADD"
        elif self.config.genetic_model == "dominant":
            model_results = results_df[results_df['TEST'] == 'DOMDEV'].copy()
            file_suffix = "DOM"
        elif self.config.genetic_model == "recessive":
            if 'REC' in results_df['TEST'].values:
                model_results = results_df[results_df['TEST'] == 'REC'].copy()
                file_suffix = "REC"
            else:
                model_results = results_df[results_df['TEST'] == 'GENO_2DF'].copy() if 'GENO_2DF' in results_df['TEST'].values else results_df[results_df['TEST'] == 'ADD'].copy()
                file_suffix = "GENO"
        else:
            model_results = results_df.copy()
            file_suffix = "ALL"
        
        # 清理和保存结果 | Clean and save results
        model_results = self._clean_results(model_results)
        if len(model_results) > 0:
            model_results.to_csv(f"gwas_results_{file_suffix}.txt", sep='\t', index=False)
            self._apply_significance_correction(model_results, file_suffix)
            return {self.config.genetic_model: model_results}
        else:
            self.logger.warning(f"⚠️ {self.config.genetic_model}模型无有效结果 | No valid {self.config.genetic_model} model results")
            return {}
    
    def _clean_results(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """清理结果数据 | Clean results data"""
        if len(results_df) == 0:
            return results_df
        
        original_count = len(results_df)
        
        # 移除无效的P值 | Remove invalid P values
        results_df = results_df.dropna(subset=['P'])
        results_df = results_df[results_df['P'] > 0]
        results_df = results_df[results_df['P'] <= 1]
        
        # 移除缺失的OR值（对于logistic回归）| Remove missing OR values (for logistic regression)
        if 'OR' in results_df.columns:
            results_df = results_df.dropna(subset=['OR'])
            results_df = results_df[results_df['OR'] > 0]
        
        cleaned_count = len(results_df)
        removed_count = original_count - cleaned_count
        
        if removed_count > 0:
            self.logger.info(f"🧹 数据清理：移除了{removed_count}个无效记录 | Data cleaning: removed {removed_count} invalid records")
        
        return results_df
    
    def _apply_significance_correction(self, results_df: pd.DataFrame, model_name: str):
        """应用显著性校正方法 | Apply significance correction methods"""
        self.logger.info(f"📈 对{model_name}模型应用显著性校正 | Applying significance correction for {model_name} model: {self.config.correction_method}")
        
        total_snps = len(results_df)
        self.logger.info(f"🔢 用于校正的SNP总数 | Total SNPs for correction: {total_snps}")
        
        if total_snps == 0:
            self.logger.warning(f"⚠️ {model_name}模型无SNP数据可用于校正 | No SNP data available for correction in {model_name} model")
            return
        
        # 根据用户选择的校正方法 | Based on user selected correction method
        if self.config.correction_method == "all":
            self._bonferroni_correction(results_df, total_snps, model_name)
            self._suggestive_correction(results_df, model_name)
            self._fdr_correction(results_df, model_name)
        elif self.config.correction_method == "bonferroni":
            self._bonferroni_correction(results_df, total_snps, model_name)
        elif self.config.correction_method == "suggestive":
            self._suggestive_correction(results_df, model_name)
        elif self.config.correction_method == "fdr":
            self._fdr_correction(results_df, model_name)
    
    def _bonferroni_correction(self, results_df: pd.DataFrame, total_snps: int, model_name: str):
        """Bonferroni校正 | Bonferroni correction"""
        bonferroni_threshold = self.config.bonferroni_alpha / total_snps
        significant = results_df[results_df['P'] < bonferroni_threshold]
        
        self.logger.info(f"🎯 {model_name}模型Bonferroni校正阈值 | {model_name} model Bonferroni threshold: {bonferroni_threshold:.2e}")
        self.logger.info(f"⭐ {model_name}模型Bonferroni显著位点数 | {model_name} model Bonferroni significant loci: {len(significant)}")
        
        significant.to_csv(f"significant_hits_bonferroni_{model_name}.txt", sep='\t', index=False)
        
        # 如果有显著位点，显示前几个 | If there are significant loci, show top ones
        if len(significant) > 0:
            top_hits = significant.nsmallest(min(5, len(significant)), 'P')
            self.logger.info(f"🏆 {model_name}模型前{len(top_hits)}个Bonferroni显著位点 | Top {len(top_hits)} Bonferroni significant loci for {model_name} model:")
            for idx, row in top_hits.iterrows():
                # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
                bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
                self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}{bp_info}")
    
    def _suggestive_correction(self, results_df: pd.DataFrame, model_name: str):
        """提示性关联阈值 | Suggestive association threshold"""
        suggestive = results_df[results_df['P'] < self.config.suggestive_threshold]
        
        self.logger.info(f"🎯 {model_name}模型提示性关联阈值 | {model_name} model suggestive threshold: {self.config.suggestive_threshold:.2e}")
        self.logger.info(f"🔍 {model_name}模型提示性关联位点数 | {model_name} model suggestive association loci: {len(suggestive)}")
        
        suggestive.to_csv(f"suggestive_hits_{model_name}.txt", sep='\t', index=False)
        
        # 如果有提示性位点，显示前几个 | If there are suggestive loci, show top ones
        if len(suggestive) > 0:
            top_hits = suggestive.nsmallest(min(10, len(suggestive)), 'P')
            self.logger.info(f"🔍 {model_name}模型前{len(top_hits)}个提示性关联位点 | Top {len(top_hits)} suggestive loci for {model_name} model:")
            for idx, row in top_hits.iterrows():
                # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
                bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
                self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}{bp_info}")
    
    def _fdr_correction(self, results_df: pd.DataFrame, model_name: str):
        """FDR校正 | FDR correction"""
        if not HAS_SCIPY:
            self.logger.warning("⚠️ scipy未安装，跳过FDR校正 | scipy not installed, skipping FDR correction")
            return
        
        try:
            p_values = results_df['P'].values
            fdr_corrected = false_discovery_control(p_values, alpha=self.config.fdr_alpha)
            
            results_df_fdr = results_df.copy()
            results_df_fdr['P_FDR'] = fdr_corrected
            
            fdr_significant = results_df_fdr[results_df_fdr['P_FDR'] < self.config.fdr_alpha]
            
            self.logger.info(f"🎯 {model_name}模型FDR校正阈值 | {model_name} model FDR threshold: {self.config.fdr_alpha}")
            self.logger.info(f"📊 {model_name}模型FDR显著位点数 | {model_name} model FDR significant loci: {len(fdr_significant)}")
            
            results_df_fdr.to_csv(f"gwas_results_fdr_{model_name}.txt", sep='\t', index=False)
            fdr_significant.to_csv(f"fdr_significant_hits_{model_name}.txt", sep='\t', index=False)
            
            if len(fdr_significant) > 0:
                top_hits = fdr_significant.nsmallest(min(5, len(fdr_significant)), 'P')
                self.logger.info(f"📊 {model_name}模型前{len(top_hits)}个FDR显著位点 | Top {len(top_hits)} FDR significant loci for {model_name} model:")
                for idx, row in top_hits.iterrows():
                    # 安全地访问列，如果BP列不存在就不显示 | Safely access columns, don't show BP if it doesn't exist
                    bp_info = f", BP={row['BP']}" if 'BP' in row and pd.notna(row['BP']) else ""
                    self.logger.info(f"   📍 {row['SNP']}: P={row['P']:.2e}, FDR_P={row['P_FDR']:.2e}, CHR={row['CHR']}{bp_info}")
        except Exception as e:
            self.logger.error(f"❌ FDR校正失败 | FDR correction failed: {e}")
    
    def _generate_model_comparison_report(self, model_results: dict):
        """生成模型比较报告 | Generate model comparison report"""
        self.logger.info("📋 生成模型比较报告 | Generating model comparison report...")
        
        comparison_report = []
        comparison_report.append("=" * 80)
        comparison_report.append("PLINK GWAS 多模型比较报告 | Multi-Model Comparison Report")
        comparison_report.append("=" * 80)
        comparison_report.append(f"分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        comparison_report.append(f"表型类型 | Trait type: {self.config.trait_type}")
        comparison_report.append(f"显著性校正方法 | Correction method: {self.config.correction_method}")
        comparison_report.append("")
        
        # 统计各模型结果 | Statistics for each model
        for model_name, results_df in model_results.items():
            comparison_report.append(f"=== {model_name.upper()}模型结果 | {model_name.upper()} Model Results ===")
            comparison_report.append(f"总SNP数 | Total SNPs: {len(results_df)}")
            
            if len(results_df) > 0:
                min_p = results_df['P'].min()
                comparison_report.append(f"最小P值 | Minimum P-value: {min_p:.2e}")
                
                if self.config.correction_method in ["bonferroni", "all"]:
                    bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
                    bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
                    comparison_report.append(f"Bonferroni显著位点 | Bonferroni significant loci (P<{bonferroni_threshold:.2e}): {bonferroni_hits}")
                
                if self.config.correction_method in ["suggestive", "all"]:
                    suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
                    comparison_report.append(f"提示性关联位点 | Suggestive association loci (P<{self.config.suggestive_threshold:.2e}): {suggestive_hits}")
                
                top_snps = results_df.nsmallest(5, 'P')
                if len(top_snps) > 0:
                    comparison_report.append("前5个最显著SNP | Top 5 most significant SNPs:")
                    for idx, row in top_snps.iterrows():
                        or_info = f", OR={row['OR']:.3f}" if 'OR' in row and pd.notna(row['OR']) else ""
                        comparison_report.append(f"  {row['SNP']}: P={row['P']:.2e}, CHR={row['CHR']}, BP={row['BP']}{or_info}")
            else:
                comparison_report.append("❌ 无有效结果 | No valid results")
            
            comparison_report.append("")
        
        # 模型解释 | Model interpretation
        comparison_report.append("=== 模型解释 | Model Interpretation ===")
        comparison_report.append("加性模型 (ADD) | Additive Model:")
        comparison_report.append("  - 假设杂合子效应是两个纯合子效应的中间值")
        comparison_report.append("  - 编码: AA=0, Aa=1, aa=2")
        comparison_report.append("  - 适用于大多数复杂性状")
        comparison_report.append("")
        comparison_report.append("显性模型 (DOMDEV) | Dominant Model:")
        comparison_report.append("  - 假设一个风险等位基因拷贝就足以产生效应")
        comparison_report.append("  - 编码: AA=0, Aa=1, aa=1")
        comparison_report.append("  - 适用于显性遗传病")
        comparison_report.append("")
        comparison_report.append("基因型模型 (GENO_2DF) | Genotypic Model:")
        comparison_report.append("  - 2自由度检验，不假设特定的遗传模式")
        comparison_report.append("  - 可以检测非加性效应")
        comparison_report.append("  - 适用于未知遗传模式的性状")
        comparison_report.append("")
        
        # 建议 | Recommendations
        comparison_report.append("=== 分析建议 | Analysis Recommendations ===")
        if len(model_results) > 1:
            best_model = self._identify_best_model(model_results)
            comparison_report.append(f"推荐模型 | Recommended model: {best_model}")
            comparison_report.append("基于最小P值和生物学合理性进行选择")
        comparison_report.append("建议进行独立群体验证 | Recommend validation in independent populations")
        comparison_report.append("考虑功能验证和精细定位 | Consider functional validation and fine mapping")
        
        # 保存报告 | Save report
        report_content = "\n".join(comparison_report)
        with open("model_comparison_report.txt", 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.logger.info("📄 模型比较报告已保存 | Model comparison report saved: model_comparison_report.txt")
    
    def _identify_best_model(self, model_results: dict) -> str:
        """识别最佳模型 | Identify best model"""
        min_p_values = {}
        
        for model_name, results_df in model_results.items():
            if len(results_df) > 0:
                min_p_values[model_name] = results_df['P'].min()
            else:
                min_p_values[model_name] = 1.0
        
        if min_p_values:
            best_model = min(min_p_values.items(), key=lambda x: x[1])[0]
            return best_model
        else:
            return "无有效模型 | No valid model"


class ReportGenerator:
    """报告生成器 | Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, model_results: dict, stats: dict):
        """生成总结报告 | Generate summary report"""
        self.logger.info("📋 生成总结报告 | Generating summary report...")
        
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("PLINK GWAS分析总结报告 | PLINK GWAS Analysis Summary Report")
        report_lines.append("=" * 80)
        report_lines.append(f"分析完成时间 | Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"遗传模型 | Genetic model: {self.config.genetic_model}")
        report_lines.append(f"表型类型 | Trait type: {self.config.trait_type}")
        report_lines.append(f"显著性校正方法 | Correction method: {self.config.correction_method}")
        report_lines.append("")
        
        # 基本统计 | Basic statistics
        report_lines.append("=== 基本统计 | Basic Statistics ===")
        report_lines.append(f"质控后样本数 | Samples after QC: {stats.get('total_samples', 'N/A')}")
        report_lines.append(f"质控后SNP数 | SNPs after QC: {stats.get('total_snps', 'N/A')}")
        
        if self.config.trait_type == "qualitative":
            report_lines.append(f"病例数 | Cases: {stats.get('cases', 'N/A')}")
            report_lines.append(f"对照数 | Controls: {stats.get('controls', 'N/A')}")
        
        report_lines.append("")
        
        # 模型特定结果 | Model-specific results
        if self.config.genetic_model == "all":
            report_lines.append("=== 多模型分析结果 | Multi-Model Analysis Results ===")
            for model_name, results_df in model_results.items():
                report_lines.append(f"{model_name.upper()}模型:")
                report_lines.append(f"  分析SNP数: {len(results_df)}")
                if len(results_df) > 0:
                    min_p = results_df['P'].min()
                    report_lines.append(f"  最小P值: {min_p:.2e}")
                    top_snp = results_df.loc[results_df['P'].idxmin()]
                    report_lines.append(f"  最显著SNP: {top_snp['SNP']} (CHR{top_snp['CHR']}:{top_snp['BP']})")
                    
                    if self.config.correction_method in ["bonferroni", "all"]:
                        bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
                        bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
                        report_lines.append(f"  Bonferroni显著位点: {bonferroni_hits}")
                    
                    if self.config.correction_method in ["suggestive", "all"]:
                        suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
                        report_lines.append(f"  提示性关联位点: {suggestive_hits}")
                else:
                    report_lines.append(f"  ❌ 无有效结果")
                report_lines.append("")
        else:
            report_lines.append(f"=== {self.config.genetic_model.upper()}模型分析结果 | {self.config.genetic_model.upper()} Model Analysis Results ===")
            if model_results and len(model_results) > 0:
                results_df = list(model_results.values())[0]
                report_lines.append(f"分析SNP数 | Analyzed SNPs: {len(results_df)}")
                if len(results_df) > 0:
                    min_p = results_df['P'].min()
                    report_lines.append(f"最小P值 | Minimum P-value: {min_p:.2e}")
                    
                    if self.config.correction_method in ["bonferroni", "all"]:
                        bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
                        bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
                        report_lines.append(f"Bonferroni显著位点 | Bonferroni significant loci (P<{bonferroni_threshold:.2e}): {bonferroni_hits}")
                    
                    if self.config.correction_method in ["suggestive", "all"]:
                        suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
                        report_lines.append(f"提示性关联位点 | Suggestive association loci (P<{self.config.suggestive_threshold:.2e}): {suggestive_hits}")
            else:
                report_lines.append("❌ 无有效结果 | No valid results")
        
        report_lines.append("")
        
        # 质量控制参数 | Quality control parameters
        report_lines.append("=== 质量控制参数 | Quality Control Parameters ===")
        report_lines.append(f"个体缺失率阈值 | Individual missing rate: {self.config.mind}")
        report_lines.append(f"SNP缺失率阈值 | SNP missing rate: {self.config.geno}")
        report_lines.append(f"最小等位基因频率 | Minor allele frequency: {self.config.maf}")
        report_lines.append(f"Hardy-Weinberg平衡P值 | HWE p-value: {self.config.hwe}")
        report_lines.append("")
        
        # 输出文件列表 | Output file list
        report_lines.append("=== 输出文件 | Output Files ===")
        if self.config.genetic_model == "all":
            report_lines.append("- gwas_results_ADD.txt: 加性模型结果")
            report_lines.append("- gwas_results_DOM.txt: 显性模型结果")
            report_lines.append("- gwas_results_GENO.txt: 基因型模型结果")
            report_lines.append("- model_comparison_report.txt: 模型比较报告")
        else:
            model_suffix = {"additive": "ADD", "dominant": "DOM", "recessive": "REC"}.get(self.config.genetic_model, "ALL")
            report_lines.append(f"- gwas_results_{model_suffix}.txt: {self.config.genetic_model}模型结果")
        
        if self.config.correction_method in ["bonferroni", "all"]:
            report_lines.append("- significant_hits_bonferroni_*.txt: Bonferroni校正显著位点")
        if self.config.correction_method in ["suggestive", "all"]:
            report_lines.append("- suggestive_hits_*.txt: 提示性关联位点")
        if self.config.correction_method in ["fdr", "all"]:
            report_lines.append("- fdr_significant_hits_*.txt: FDR校正显著位点")
        
        report_lines.append("")
        
        # 建议和注意事项 | Recommendations and notes
        report_lines.append("=== 建议和注意事项 | Recommendations and Notes ===")
        report_lines.append("1. 显著位点需要在独立群体中验证")
        report_lines.append("2. 考虑进行功能注释和通路分析")
        report_lines.append("3. 注意群体分层和家系结构的影响")
        if self.config.genetic_model == "all":
            report_lines.append("4. 比较不同模型结果，选择生物学上最合理的模型")
        report_lines.append("5. 建议结合其他组学数据进行综合分析")
        
        # 保存报告 | Save report
        report_content = "\n".join(report_lines)
        with open("gwas_summary_report.txt", 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.logger.info("📄 GWAS总结报告已保存 | GWAS summary report saved: gwas_summary_report.txt")


class VisualizationGenerator:
    """可视化生成器 | Visualization Generator"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def generate_plots(self, model_results: dict):
        """生成可视化图表 | Generate visualization plots"""
        self.logger.info("📊 生成可视化图表 | Generating visualization plots...")
        
        if not model_results:
            self.logger.warning("⚠️ 无结果数据，跳过可视化 | No result data, skipping visualization")
            return
        
        if self.config.genetic_model == "all":
            # 为每个模型生成图表 | Generate plots for each model
            for model_name, results_df in model_results.items():
                if len(results_df) > 0:
                    self._generate_model_plots(results_df, model_name)
            
            # 生成模型比较图表 | Generate model comparison plots
            self._generate_comparison_plots(model_results)
        else:
            # 为单一模型生成图表 | Generate plots for single model
            if model_results and len(model_results) > 0:
                results_df = list(model_results.values())[0]
                if len(results_df) > 0:
                    self._generate_model_plots(results_df, self.config.genetic_model)
    
    def _generate_model_plots(self, results_df: pd.DataFrame, model_name: str):
        """为特定模型生成图表 | Generate plots for specific model"""
        try:
            self._create_manhattan_plot(results_df, model_name)
            self._create_qq_plot(results_df, model_name)
            self.logger.info(f"📈 {model_name}模型图表生成完成 | {model_name} model plots generated")
        except Exception as e:
            self.logger.warning(f"⚠️ {model_name}模型图表生成失败 | Plot generation failed for {model_name}: {e}")
    
    def _create_manhattan_plot(self, results_df: pd.DataFrame, model_name: str):
        """创建Manhattan图 | Create Manhattan plot"""
        model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
        
        r_script = f"""
library(ggplot2)
library(dplyr)

# 读取数据
tryCatch({{
    data <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
    
    # 检查数据
    if(nrow(data) == 0) {{
        cat("No data available for Manhattan plot\\n")
        quit()
    }}
    
    # 计算-log10(P)
    data$log_p <- -log10(data$P)
    
    # 处理染色体信息
    data$CHR <- as.factor(data$CHR)
    
    # 创建Manhattan图
    png("manhattan_plot_{model_name}.png", width=1200, height=800, res=150)
    p <- ggplot(data, aes(x=BP, y=log_p, color=CHR)) +
        geom_point(alpha=0.6, size=0.8) +
        theme_minimal() +
        labs(title="Manhattan Plot - {model_name.upper()} Model",
             x="Genomic Position", 
             y="-log10(P-value)",
             color="Chromosome") +
        theme(legend.position="bottom",
              axis.text.x = element_text(angle=45, hjust=1))
    
    # 添加显著性阈值线
    if(max(data$log_p, na.rm=TRUE) > 5) {{
        p <- p + geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed", alpha=0.7)
    }}
    if(max(data$log_p, na.rm=TRUE) > 3) {{
        p <- p + geom_hline(yintercept=-log10(1e-5), color="blue", linetype="dashed", alpha=0.7)
    }}
    
    print(p)
    dev.off()
    
    cat("Manhattan plot created successfully\\n")
}}, error = function(e) {{
    cat("Error creating Manhattan plot:", conditionMessage(e), "\\n")
}})
"""
        
        with open(f"manhattan_{model_name}.R", 'w') as f:
            f.write(r_script)
    
    def _create_qq_plot(self, results_df: pd.DataFrame, model_name: str):
        """创建QQ图 | Create QQ plot"""
        model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
        
        r_script = f"""
library(ggplot2)

# 读取数据
tryCatch({{
    data <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
    
    # 检查数据
    if(nrow(data) == 0) {{
        cat("No data available for QQ plot\\n")
        quit()
    }}
    
    # 移除无效P值
    data <- data[!is.na(data$P) & data$P > 0 & data$P <= 1, ]
    
    if(nrow(data) == 0) {{
        cat("No valid P-values for QQ plot\\n")
        quit()
    }}
    
    # 计算QQ图数据
    observed <- sort(-log10(data$P))
    expected <- sort(-log10(ppoints(length(observed))))
    
    qq_data <- data.frame(expected=expected, observed=observed)
    
    # 创建QQ图
    png("qq_plot_{model_name}.png", width=800, height=800, res=150)
    p <- ggplot(qq_data, aes(x=expected, y=observed)) +
        geom_point(alpha=0.6) +
        geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
        theme_minimal() +
        labs(title="Q-Q Plot - {model_name.upper()} Model",
             x="Expected -log10(P)",
             y="Observed -log10(P)")
    
    print(p)
    dev.off()
    
    cat("QQ plot created successfully\\n")
}}, error = function(e) {{
    cat("Error creating QQ plot:", conditionMessage(e), "\\n")
}})
"""
        
        with open(f"qq_{model_name}.R", 'w') as f:
            f.write(r_script)
    
    def _generate_comparison_plots(self, model_results: dict):
        """生成模型比较图表 | Generate model comparison plots"""
        self.logger.info("📊 生成模型比较图表 | Generating model comparison plots...")
        
        # 生成模型比较的R脚本
        r_script = """
library(ggplot2)
library(dplyr)
library(reshape2)

# 比较不同模型的P值分布
models <- c()
p_values <- c()

"""
        
        for model_name in model_results.keys():
            model_suffix = {"additive": "ADD", "dominant": "DOM", "genotypic": "GENO", "recessive": "REC"}.get(model_name, model_name.upper())
            r_script += f"""
if(file.exists("gwas_results_{model_suffix}.txt")) {{
    tryCatch({{
        data_{model_name} <- read.table("gwas_results_{model_suffix}.txt", header=TRUE, sep="\\t")
        if(nrow(data_{model_name}) > 0) {{
            models <- c(models, rep("{model_name}", nrow(data_{model_name})))
            p_values <- c(p_values, data_{model_name}$P)
        }}
    }}, error = function(e) {{
        cat("Error reading {model_suffix} data:", conditionMessage(e), "\\n")
    }})
}}
"""
        
        r_script += """
# 检查是否有数据
if(length(models) == 0 || length(p_values) == 0) {
    cat("No data available for comparison plots\\n")
    quit()
}

# 创建比较数据框
comparison_data <- data.frame(Model=models, P_value=p_values)
comparison_data <- comparison_data[!is.na(comparison_data$P_value) & 
                                  comparison_data$P_value > 0 & 
                                  comparison_data$P_value <= 1, ]

if(nrow(comparison_data) == 0) {
    cat("No valid data for comparison plots\\n")
    quit()
}

comparison_data$log_p <- -log10(comparison_data$P_value)

# 生成P值分布比较图
tryCatch({
    png("model_comparison_pvalue_distribution.png", width=1200, height=800, res=150)
    p <- ggplot(comparison_data, aes(x=log_p, fill=Model)) +
        geom_histogram(alpha=0.7, bins=50, position="identity") +
        facet_wrap(~Model, scales="free_y") +
        theme_minimal() +
        labs(title="P-value Distribution Comparison Across Models",
             x="-log10(P-value)",
             y="Frequency") +
        theme(legend.position="bottom")
    print(p)
    dev.off()
    
    cat("Model comparison plot created successfully\\n")
}, error = function(e) {
    cat("Error creating comparison plot:", conditionMessage(e), "\\n")
})
"""
        
        with open("model_comparison.R", 'w') as f:
            f.write(r_script)
        
        self.logger.info("📄 可视化R脚本已生成 | Visualization R scripts generated")