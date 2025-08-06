"""
PLINK GWAS结果处理模块 | PLINK GWAS Results Processing Module
"""

import pandas as pd
import numpy as np
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
    
    def process_results(self, main_result: str) -> pd.DataFrame:
        """处理分析结果 | Process analysis results"""
        self.logger.info("处理分析结果 | Processing analysis results...")
        
        # 根据表型类型选择结果文件 | Select result file based on trait type
        if self.config.trait_type == "qualitative":
            result_file = f"{main_result}.assoc.logistic"
        else:
            result_file = f"{main_result}.assoc.linear"
        
        if not Path(result_file).exists():
            raise FileNotFoundError(f"结果文件不存在 | Result file does not exist: {result_file}")
        
        # 读取结果文件 | Read results file
        try:
            results_df = pd.read_csv(result_file, sep=r"\s+")
        except Exception as e:
            self.logger.error(f"读取结果文件失败 | Failed to read results file: {e}")
            return None
        
        # 提取ADD模型结果 | Extract ADD model results
        add_results = results_df[results_df['TEST'] == 'ADD'].copy()
        
        # 移除无效的P值 | Remove invalid P values
        add_results = add_results.dropna(subset=['P'])
        add_results = add_results[add_results['P'] > 0]
        
        add_results.to_csv("gwas_results_ADD.txt", sep='\t', index=False)
        
        # 应用显著性校正 | Apply significance correction
        self._apply_significance_correction(add_results)
        
        return add_results
    
    def _apply_significance_correction(self, results_df: pd.DataFrame):
        """应用显著性校正方法 | Apply significance correction methods"""
        self.logger.info(f"应用显著性校正方法 | Applying significance correction method: {self.config.correction_method}")
        
        total_snps = len(results_df)
        self.logger.info(f"用于校正的SNP总数 | Total SNPs for correction: {total_snps}")
        
        # 根据用户选择的校正方法 | Based on user selected correction method
        if self.config.correction_method == "all":
            # 应用所有三种方法 | Apply all three methods
            self._bonferroni_correction(results_df, total_snps)
            self._suggestive_correction(results_df)
            self._fdr_correction(results_df)
        elif self.config.correction_method == "bonferroni":
            self._bonferroni_correction(results_df, total_snps)
        elif self.config.correction_method == "suggestive":
            self._suggestive_correction(results_df)
        elif self.config.correction_method == "fdr":
            self._fdr_correction(results_df)
    
    def _bonferroni_correction(self, results_df: pd.DataFrame, total_snps: int):
        """Bonferroni校正 | Bonferroni correction"""
        self.logger.info("执行Bonferroni校正 | Performing Bonferroni correction")
        
        bonferroni_threshold = self.config.bonferroni_alpha / total_snps
        significant_bonferroni = results_df[results_df['P'] < bonferroni_threshold].copy()
        
        # 按P值排序 | Sort by P value
        significant_bonferroni = significant_bonferroni.sort_values('P')
        
        # 保存结果 | Save results
        significant_bonferroni.to_csv("significant_bonferroni.txt", sep='\t', index=False)
        
        self.logger.info(f"Bonferroni校正阈值 | Bonferroni threshold: {bonferroni_threshold:.2e}")
        self.logger.info(f"Bonferroni显著SNP数 | Bonferroni significant SNPs: {len(significant_bonferroni)}")
        
        # 保存阈值信息 | Save threshold information
        with open("bonferroni_info.txt", "w") as f:
            f.write(f"Bonferroni Correction Information\n")
            f.write(f"Alpha level: {self.config.bonferroni_alpha}\n")
            f.write(f"Total SNPs: {total_snps}\n")
            f.write(f"Corrected threshold: {bonferroni_threshold:.2e}\n")
            f.write(f"Significant SNPs: {len(significant_bonferroni)}\n")
    
    def _suggestive_correction(self, results_df: pd.DataFrame):
        """提示性关联阈值 | Suggestive significance threshold"""
        self.logger.info("执行提示性关联筛选 | Performing suggestive significance filtering")
        
        suggestive_threshold = self.config.suggestive_threshold
        significant_suggestive = results_df[results_df['P'] < suggestive_threshold].copy()
        
        # 按P值排序 | Sort by P value
        significant_suggestive = significant_suggestive.sort_values('P')
        
        # 保存结果 | Save results
        significant_suggestive.to_csv("significant_suggestive.txt", sep='\t', index=False)
        
        self.logger.info(f"提示性关联阈值 | Suggestive threshold: {suggestive_threshold:.2e}")
        self.logger.info(f"提示性关联SNP数 | Suggestive significant SNPs: {len(significant_suggestive)}")
        
        # 保存阈值信息 | Save threshold information
        with open("suggestive_info.txt", "w") as f:
            f.write(f"Suggestive Significance Information\n")
            f.write(f"Threshold: {suggestive_threshold:.2e}\n")
            f.write(f"Significant SNPs: {len(significant_suggestive)}\n")
    
    def _fdr_correction(self, results_df: pd.DataFrame):
        """FDR校正 | FDR correction"""
        self.logger.info("执行FDR校正 | Performing FDR correction")
        
        if not HAS_SCIPY:
            self.logger.warning("未安装scipy，跳过FDR校正 | scipy not installed, skipping FDR correction")
            self.logger.warning("请安装scipy: pip install scipy | Please install scipy: pip install scipy")
            return
        
        # 获取P值 | Get P values
        p_values = results_df['P'].values
        
        try:
            # 使用scipy进行FDR校正 | Use scipy for FDR correction
            rejected, p_corrected = false_discovery_control(p_values, alpha=self.config.fdr_alpha, method='bh')
            
            # 添加校正后的P值和是否显著的标记 | Add corrected P values and significance flag
            results_df['P_FDR'] = p_corrected
            results_df['FDR_significant'] = rejected
            
            # 筛选显著的SNP | Filter significant SNPs
            significant_fdr = results_df[results_df['FDR_significant']].copy()
            
            # 按原始P值排序 | Sort by original P value
            significant_fdr = significant_fdr.sort_values('P')
            
            # 保存结果 | Save results
            significant_fdr.to_csv("significant_fdr.txt", sep='\t', index=False)
            
            self.logger.info(f"FDR校正q值阈值 | FDR q-value threshold: {self.config.fdr_alpha}")
            self.logger.info(f"FDR显著SNP数 | FDR significant SNPs: {len(significant_fdr)}")
            
            # 保存阈值信息 | Save threshold information
            with open("fdr_info.txt", "w") as f:
                f.write(f"FDR Correction Information\n")
                f.write(f"Method: Benjamini-Hochberg\n")
                f.write(f"Alpha (q-value): {self.config.fdr_alpha}\n")
                f.write(f"Significant SNPs: {len(significant_fdr)}\n")
                f.write(f"Highest corrected P-value: {significant_fdr['P_FDR'].max():.2e}\n" if len(significant_fdr) > 0 else "No significant SNPs\n")
            
        except Exception as e:
            self.logger.error(f"FDR校正失败 | FDR correction failed: {e}")

class ReportGenerator:
    """报告生成器 | Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_report(self, results_df: pd.DataFrame, main_result: str):
        """生成分析报告 | Generate analysis report"""
        self.logger.info("生成分析报告 | Generating analysis report...")
        
        # 统计信息 | Statistics
        stats = {
            'analysis_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'main_analysis': main_result,
            'trait_type': self.config.trait_type,
            'correction_method': self.config.correction_method,
            'total_snps': 0,
            'total_samples': 0,
            'cases': 0,
            'controls': 0,
            'chromosomes': 0
        }
        
        # 读取文件统计 | Read file statistics
        if Path("data_qc1.bim").exists():
            bim_df = pd.read_csv("data_qc1.bim", sep="\t", header=None)
            stats['total_snps'] = len(bim_df)
            stats['chromosomes'] = bim_df.iloc[:, 0].nunique()
        
        if Path("data_qc1.fam").exists():
            fam_df = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
            stats['total_samples'] = len(fam_df)
            if self.config.trait_type == "qualitative":
                stats['cases'] = (fam_df.iloc[:, 5] == 2).sum()
                stats['controls'] = (fam_df.iloc[:, 5] == 1).sum()
            else:
                # 数量性状统计 | Quantitative trait statistics
                pheno_values = fam_df.iloc[:, 5][fam_df.iloc[:, 5] != -9]
                stats['mean_phenotype'] = pheno_values.mean()
                stats['std_phenotype'] = pheno_values.std()
        
        # 统计各种校正方法的结果 | Count results for each correction method
        correction_results = {}
        
        if self.config.correction_method in ["bonferroni", "all"] and Path("significant_bonferroni.txt").exists():
            bonferroni_df = pd.read_csv("significant_bonferroni.txt", sep='\t')
            correction_results['bonferroni'] = len(bonferroni_df)
        
        if self.config.correction_method in ["suggestive", "all"] and Path("significant_suggestive.txt").exists():
            suggestive_df = pd.read_csv("significant_suggestive.txt", sep='\t')
            correction_results['suggestive'] = len(suggestive_df)
        
        if self.config.correction_method in ["fdr", "all"] and Path("significant_fdr.txt").exists():
            fdr_df = pd.read_csv("significant_fdr.txt", sep='\t')
            correction_results['fdr'] = len(fdr_df)
        
        # 生成报告 | Generate report
        report = f"""
=== PLINK GWAS分析完成报告 | PLINK GWAS Analysis Report ===
分析时间 | Analysis time: {stats['analysis_time']}
主要分析 | Main analysis: {stats['main_analysis']}
表型类型 | Trait type: {stats['trait_type']}
显著性校正方法 | Significance correction method: {stats['correction_method']}
染色体类型 | Chromosome type: 支持非标准编号（如OV开头） | Support non-standard naming (e.g., OV prefix)

=== 样本统计 | Sample Statistics ===
质控后总样本数 | Total samples after QC: {stats['total_samples']}"""
        
        if self.config.trait_type == "qualitative":
            report += f"""
病例数 (感病/原值1) | Cases (susceptible/original 1): {stats['cases']}
对照数 (抗病/原值0) | Controls (resistant/original 0): {stats['controls']}"""
        else:
            report += f"""
表型均值 | Phenotype mean: {stats.get('mean_phenotype', 'N/A'):.4f}
表型标准差 | Phenotype std: {stats.get('std_phenotype', 'N/A'):.4f}"""
        
        report += f"""

=== SNP统计 | SNP Statistics ===
质控后总SNP数 | Total SNPs after QC: {stats['total_snps']}
染色体数 | Number of chromosomes: {stats['chromosomes']}

=== 显著性校正结果 | Significance Correction Results ==="""
        
        if 'bonferroni' in correction_results:
            bonferroni_threshold = self.config.bonferroni_alpha / stats['total_snps']
            report += f"""
Bonferroni校正结果 | Bonferroni Correction Results:
  - 校正阈值 | Corrected threshold: {bonferroni_threshold:.2e}
  - 显著SNP数 | Significant SNPs: {correction_results['bonferroni']}"""
        
        if 'suggestive' in correction_results:
            report += f"""
提示性关联结果 | Suggestive Association Results:
  - 阈值 | Threshold: {self.config.suggestive_threshold:.2e}
  - 显著SNP数 | Significant SNPs: {correction_results['suggestive']}"""
        
        if 'fdr' in correction_results:
            report += f"""
FDR校正结果 | FDR Correction Results:
  - q值阈值 | q-value threshold: {self.config.fdr_alpha}
  - 显著SNP数 | Significant SNPs: {correction_results['fdr']}"""
        
        report += f"""

=== 分析参数 | Analysis Parameters ===
表型类型 | Trait type: {self.config.trait_type}
显著性校正方法 | Correction method: {self.config.correction_method}
个体缺失率阈值 | Individual missing rate threshold: {self.config.mind}
SNP缺失率阈值 | SNP missing rate threshold: {self.config.geno}
最小等位基因频率 | Minor allele frequency: {self.config.maf}
Hardy-Weinberg平衡P值 | Hardy-Weinberg equilibrium p-value: {self.config.hwe}
LD r²阈值 | LD r² threshold: {self.config.ld_r2_threshold}
使用主成分数 | Number of PCs used: {self.config.pca_use}

=== 校正方法说明 | Correction Method Description ==="""
        
        if self.config.correction_method in ["bonferroni", "all"]:
            report += f"""
Bonferroni校正: 控制家族错误率(FWER)，阈值 = {self.config.bonferroni_alpha}/总SNP数
Bonferroni correction: Controls family-wise error rate (FWER), threshold = {self.config.bonferroni_alpha}/total SNPs"""
        
        if self.config.correction_method in ["suggestive", "all"]:
            report += f"""
提示性关联: 传统GWAS中的提示性显著性阈值，P < {self.config.suggestive_threshold:.0e}
Suggestive association: Traditional suggestive significance threshold in GWAS, P < {self.config.suggestive_threshold:.0e}"""
        
        if self.config.correction_method in ["fdr", "all"]:
            report += f"""
FDR校正: 控制假发现率，使用Benjamini-Hochberg方法，q值 < {self.config.fdr_alpha}
FDR correction: Controls false discovery rate using Benjamini-Hochberg method, q-value < {self.config.fdr_alpha}"""
        
        report += f"""

=== 结果解释 | Results Interpretation ==="""
        
        if self.config.trait_type == "qualitative":
            report += """
- OR > 1: 该等位基因增加感病风险 | This allele increases susceptibility risk
- OR < 1: 该等位基因具有保护作用（促进抗病） | This allele has protective effect (promotes resistance)"""
        else:
            report += """
- BETA > 0: 该等位基因增加表型值 | This allele increases phenotype value
- BETA < 0: 该等位基因降低表型值 | This allele decreases phenotype value"""
        
        report += """
- 显著位点需要在独立群体中验证 | Significant loci need validation in independent populations

=== 输出文件 | Output Files ===
- gwas_results_ADD.txt: 主要关联分析结果 | Main association analysis results"""
        
        if self.config.correction_method in ["bonferroni", "all"]:
            report += """
- significant_bonferroni.txt: Bonferroni校正显著位点 | Bonferroni significant loci
- bonferroni_info.txt: Bonferroni校正信息 | Bonferroni correction information"""
        
        if self.config.correction_method in ["suggestive", "all"]:
            report += """
- significant_suggestive.txt: 提示性关联位点 | Suggestive association loci
- suggestive_info.txt: 提示性关联信息 | Suggestive association information"""
        
        if self.config.correction_method in ["fdr", "all"]:
            report += """
- significant_fdr.txt: FDR校正显著位点 | FDR significant loci
- fdr_info.txt: FDR校正信息 | FDR correction information"""
        
        report += """
- manhattan_plot.png: Manhattan图 | Manhattan plot
- qq_plot.png: QQ图 | QQ plot
"""
        
        # 如果有显著位点，显示每种方法的top hits | Show top hits for each method if significant loci exist
        if correction_results:
            report += "\n=== 各校正方法的前5个最显著位点 | Top 5 Most Significant Loci by Each Method ===\n"
            
            for method, count in correction_results.items():
                if count > 0:
                    method_file = f"significant_{method}.txt"
                    if Path(method_file).exists():
                        method_df = pd.read_csv(method_file, sep='\t')
                        top_hits = method_df.head(5)
                        
                        report += f"\n{method.upper()}方法 | {method.upper()} Method:\n"
                        if self.config.trait_type == "qualitative":
                            report += "CHR\tSNP\tBP\tA1\tOR\tP\n"
                            for _, row in top_hits.iterrows():
                                report += f"{row['CHR']}\t{row['SNP']}\t{row['BP']}\t{row['A1']}\t{row['OR']:.3f}\t{row['P']:.2e}\n"
                        else:
                            report += "CHR\tSNP\tBP\tA1\tBETA\tP\n"
                            for _, row in top_hits.iterrows():
                                report += f"{row['CHR']}\t{row['SNP']}\t{row['BP']}\t{row['A1']}\t{row['BETA']:.3f}\t{row['P']:.2e}\n"
        
        # 保存报告 | Save report
        with open("analysis_report.txt", "w", encoding='utf-8') as f:
            f.write(report)
        
        self.logger.info("分析报告已生成 | Analysis report generated: analysis_report.txt")

class VisualizationGenerator:
    """可视化生成器 | Visualization Generator"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def generate_plots(self):
        """生成可视化图形 | Generate visualization plots"""
        self.logger.info("生成可视化图形 | Generating visualization plots...")
        
        try:
            # 创建R脚本 | Create R script
            r_script = f'''
# 安装和载入必要的包 | Install and load required packages
if (!require("qqman", quietly=TRUE)) {{
    install.packages("qqman", repos="https://cran.rstudio.com/")
    library(qqman)
}}

# 检查结果文件是否存在 | Check if results file exists
if (!file.exists("gwas_results_ADD.txt")) {{
    stop("错误: 找不到 gwas_results_ADD.txt 文件 | Error: Cannot find gwas_results_ADD.txt file")
}}

# 读取GWAS结果 | Read GWAS results
cat("读取GWAS结果... | Reading GWAS results...\\n")
gwas <- read.table("gwas_results_ADD.txt", header=TRUE, stringsAsFactors=FALSE)

if (nrow(gwas) == 0) {{
    stop("错误: GWAS结果文件为空 | Error: GWAS results file is empty")
}}

cat("数据行数 | Data rows:", nrow(gwas), "\\n")
cat("表型类型 | Trait type: {self.config.trait_type}\\n")
cat("校正方法 | Correction method: {self.config.correction_method}\\n")

# 处理非标准染色体编号 | Handle non-standard chromosome naming
gwas$CHR_original <- gwas$CHR
gwas$CHR_numeric <- as.numeric(gsub("OV0?", "", gwas$CHR))

# 如果转换失败，使用序号 | If conversion fails, use index
if (any(is.na(gwas$CHR_numeric))) {{
    unique_chrs <- unique(gwas$CHR)
    chr_mapping <- setNames(1:length(unique_chrs), unique_chrs)
    gwas$CHR_numeric <- chr_mapping[gwas$CHR]
}}

# 准备绘图数据 | Prepare plotting data
gwas_plot <- data.frame(
    SNP = gwas$SNP,
    CHR = gwas$CHR_numeric,
    BP = gwas$BP,
    P = gwas$P
)

# 移除无效数据 | Remove invalid data
gwas_plot <- gwas_plot[!is.na(gwas_plot$P) & gwas_plot$P > 0 & gwas_plot$P <= 1, ]

if (nrow(gwas_plot) == 0) {{
    stop("错误: 没有有效的P值数据 | Error: No valid P-value data")
}}

cat("有效数据点数 | Valid data points:", nrow(gwas_plot), "\\n")

# 设置显著性阈值 | Set significance thresholds
total_snps <- nrow(gwas_plot)
bonferroni_threshold <- {self.config.bonferroni_alpha} / total_snps
suggestive_threshold <- {self.config.suggestive_threshold}

# 生成Manhattan图 | Generate Manhattan plot
cat("生成Manhattan图... | Generating Manhattan plot...\\n")
trait_label <- if ("{self.config.trait_type}" == "qualitative") "Disease Resistance (Qualitative)" else "Quantitative Trait"
correction_label <- "{self.config.correction_method}"

png("manhattan_plot.png", width=1400, height=700, res=100)
tryCatch({{
    manhattan(gwas_plot, 
              main=paste("PLINK GWAS Manhattan Plot -", trait_label, "- Correction:", correction_label), 
              genomewideline = -log10(bonferroni_threshold),
              suggestiveline = -log10(suggestive_threshold),
              col = c("blue4", "orange3"))
    
    # 添加阈值标签 | Add threshold labels
    abline(h = -log10(bonferroni_threshold), col = "red", lty = 2, lwd = 1)
    abline(h = -log10(suggestive_threshold), col = "blue", lty = 2, lwd = 1)
    
    # 添加图例 | Add legend
    legend("topright", 
           legend = c(paste("Bonferroni:", format(bonferroni_threshold, scientific=TRUE, digits=2)),
                     paste("Suggestive:", format(suggestive_threshold, scientific=TRUE, digits=2))),
           col = c("red", "blue"), lty = 2, cex = 0.8)
}}, error = function(e) {{
    plot(1, type="n", main="Manhattan Plot Error", xlab="Position", ylab="-log10(P)")
    text(1, 1, paste("Error:", e$message), cex=0.8)
}})
dev.off()

# 生成QQ图 | Generate QQ plot
cat("生成QQ图... | Generating QQ plot...\\n")
png("qq_plot.png", width=600, height=600, res=100)
tryCatch({{
    qq(gwas_plot$P, main=paste("QQ Plot -", trait_label, "- Correction:", correction_label))
}}, error = function(e) {{
    plot(1, type="n", main="QQ Plot Error", xlab="Expected", ylab="Observed")
    text(1, 1, paste("Error:", e$message), cex=0.8)
}})
dev.off()

# 保存染色体映射 | Save chromosome mapping
write.table(data.frame(
    Original_CHR = unique(gwas$CHR_original),
    Numeric_CHR = unique(gwas$CHR_numeric)
), "chromosome_mapping.txt", row.names=FALSE, quote=FALSE, sep="\\t")

# 保存阈值信息 | Save threshold information
cat("\\n=== 显著性阈值信息 | Significance Threshold Information ===\\n")
cat("Bonferroni阈值 | Bonferroni threshold:", format(bonferroni_threshold, scientific=TRUE, digits=3), "\\n")
cat("提示性阈值 | Suggestive threshold:", format(suggestive_threshold, scientific=TRUE, digits=3), "\\n")

# 统计各阈值下的显著SNP数 | Count significant SNPs under each threshold
bonferroni_count <- sum(gwas_plot$P < bonferroni_threshold)
suggestive_count <- sum(gwas_plot$P < suggestive_threshold)

cat("Bonferroni显著SNP数 | Bonferroni significant SNPs:", bonferroni_count, "\\n")
cat("提示性显著SNP数 | Suggestive significant SNPs:", suggestive_count, "\\n")

cat("\\n=== 可视化完成 | Visualization completed ===\\n")
cat("生成的文件 | Generated files:\\n")
cat("- manhattan_plot.png\\n") 
cat("- qq_plot.png\\n")
cat("- chromosome_mapping.txt\\n")

cat("\\n=== 数据统计 | Data statistics ===\\n")
cat("总SNP数 | Total SNPs:", nrow(gwas_plot), "\\n")
cat("染色体数 | Number of chromosomes:", length(unique(gwas_plot$CHR)), "\\n")
cat("最小P值 | Minimum P-value:", min(gwas_plot$P), "\\n")
'''
            
            # 保存R脚本 | Save R script
            with open("plot_results.R", "w") as f:
                f.write(r_script)
            
            # 运行R脚本 | Run R script
            result = self.cmd_runner.run(["Rscript", "plot_results.R"], 
                                       "生成可视化图形 | Generating visualization plots", 
                                       check=False)
            
            if result.returncode == 0:
                self.logger.info("可视化图形生成成功 | Visualization plots generated successfully")
            else:
                self.logger.warning("可视化图形生成失败，请手动运行 | Visualization generation failed, please run manually: Rscript plot_results.R")
        
        except Exception as e:
            self.logger.error(f"生成可视化脚本失败 | Failed to generate visualization script: {e}")
