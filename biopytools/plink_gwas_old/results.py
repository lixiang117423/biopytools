"""
PLINK GWAS结果处理模块 | PLINK GWAS Results Processing Module
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_results(self, main_result: str) -> pd.DataFrame:
        """处理分析结果 | Process analysis results"""
        self.logger.info("处理分析结果 | Processing analysis results...")
        
        result_file = f"{main_result}.assoc.logistic"
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
        add_results.to_csv("gwas_results_ADD.txt", sep='\t', index=False)
        
        # 筛选显著位点 | Filter significant loci
        sig_threshold = float(self.config.significance_threshold)
        sug_threshold = float(self.config.suggestive_threshold)
        
        significant = add_results[add_results['P'] < sig_threshold]
        suggestive = add_results[add_results['P'] < sug_threshold]
        
        significant.to_csv("significant_hits.txt", sep='\t', index=False)
        suggestive.to_csv("suggestive_hits.txt", sep='\t', index=False)
        
        self.logger.info(f"总SNP数 | Total SNPs: {len(add_results)}")
        self.logger.info(f"全基因组显著位点 | Genome-wide significant loci (P<{sig_threshold}): {len(significant)}")
        self.logger.info(f"提示性关联位点 | Suggestive association loci (P<{sug_threshold}): {len(suggestive)}")
        
        return add_results

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
            'total_snps': 0,
            'total_samples': 0,
            'cases': 0,
            'controls': 0,
            'chromosomes': 0,
            'significant_hits': 0,
            'suggestive_hits': 0
        }
        
        # 读取文件统计 | Read file statistics
        if Path("data_qc1.bim").exists():
            bim_df = pd.read_csv("data_qc1.bim", sep="\t", header=None)
            stats['total_snps'] = len(bim_df)
            stats['chromosomes'] = bim_df.iloc[:, 0].nunique()
        
        if Path("data_qc1.fam").exists():
            fam_df = pd.read_csv("data_qc1.fam", sep=r"\s+", header=None)
            stats['total_samples'] = len(fam_df)
            stats['cases'] = (fam_df.iloc[:, 5] == 2).sum()
            stats['controls'] = (fam_df.iloc[:, 5] == 1).sum()
        
        if results_df is not None:
            sig_threshold = float(self.config.significance_threshold)
            sug_threshold = float(self.config.suggestive_threshold)
            stats['significant_hits'] = (results_df['P'] < sig_threshold).sum()
            stats['suggestive_hits'] = (results_df['P'] < sug_threshold).sum()
        
        # 生成报告 | Generate report
        report = f"""
=== PLINK GWAS分析完成报告 | PLINK GWAS Analysis Report ===
分析时间 | Analysis time: {stats['analysis_time']}
主要分析 | Main analysis: {stats['main_analysis']}
染色体类型 | Chromosome type: 支持非标准编号（如OV开头） | Support non-standard naming (e.g., OV prefix)

=== 样本统计 | Sample Statistics ===
质控后总样本数 | Total samples after QC: {stats['total_samples']}
病例数 (感病) | Cases (susceptible): {stats['cases']}
对照数 (抗病) | Controls (resistant): {stats['controls']}

=== SNP统计 | SNP Statistics ===
质控后总SNP数 | Total SNPs after QC: {stats['total_snps']}
染色体数 | Number of chromosomes: {stats['chromosomes']}

=== 关联分析结果 | Association Analysis Results ===
全基因组显著位点数 | Genome-wide significant loci (P<{self.config.significance_threshold}): {stats['significant_hits']}
提示性关联位点数 | Suggestive association loci (P<{self.config.suggestive_threshold}): {stats['suggestive_hits']}

=== 分析参数 | Analysis Parameters ===
个体缺失率阈值 | Individual missing rate threshold: {self.config.mind}
SNP缺失率阈值 | SNP missing rate threshold: {self.config.geno}
最小等位基因频率 | Minor allele frequency: {self.config.maf}
Hardy-Weinberg平衡P值 | Hardy-Weinberg equilibrium p-value: {self.config.hwe}
LD r²阈值 | LD r² threshold: {self.config.ld_r2_threshold}
使用主成分数 | Number of PCs used: {self.config.pca_use}

=== 结果解释 | Results Interpretation ===
- OR > 1: 该等位基因增加感病风险 | This allele increases susceptibility risk
- OR < 1: 该等位基因具有保护作用（促进抗病） | This allele has protective effect (promotes resistance)
- 显著位点需要在独立群体中验证 | Significant loci need validation in independent populations

=== 输出文件 | Output Files ===
- gwas_results_ADD.txt: 主要关联分析结果 | Main association analysis results
- significant_hits.txt: 显著关联位点 | Significant association loci
- suggestive_hits.txt: 提示性关联位点 | Suggestive association loci
- manhattan_plot.png: Manhattan图（需运行可视化） | Manhattan plot (requires visualization)
- qq_plot.png: QQ图（需运行可视化） | QQ plot (requires visualization)
"""
        
        # 如果有显著位点，添加top hits | Add top hits if significant loci exist
        if results_df is not None and len(results_df) > 0:
            report += "\n前10个最显著的关联位点 | Top 10 most significant association loci:\n"
            report += "CHR\tSNP\tBP\tA1\tOR\tP\n"
            top_hits = results_df.nsmallest(10, 'P')[['CHR', 'SNP', 'BP', 'A1', 'OR', 'P']]
            for _, row in top_hits.iterrows():
                report += f"{row['CHR']}\t{row['SNP']}\t{row['BP']}\t{row['A1']}\t{row['OR']:.3f}\t{row['P']:.2e}\n"
        
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
            r_script = '''
# 安装和载入必要的包 | Install and load required packages
if (!require("qqman", quietly=TRUE)) {
    install.packages("qqman", repos="https://cran.rstudio.com/")
    library(qqman)
}

# 检查结果文件是否存在 | Check if results file exists
if (!file.exists("gwas_results_ADD.txt")) {
    stop("错误: 找不到 gwas_results_ADD.txt 文件 | Error: Cannot find gwas_results_ADD.txt file")
}

# 读取GWAS结果 | Read GWAS results
cat("读取GWAS结果... | Reading GWAS results...\\n")
gwas <- read.table("gwas_results_ADD.txt", header=TRUE, stringsAsFactors=FALSE)

if (nrow(gwas) == 0) {
    stop("错误: GWAS结果文件为空 | Error: GWAS results file is empty")
}

cat("数据行数 | Data rows:", nrow(gwas), "\\n")

# 处理非标准染色体编号 | Handle non-standard chromosome naming
gwas$CHR_original <- gwas$CHR
gwas$CHR_numeric <- as.numeric(gsub("OV0?", "", gwas$CHR))

# 如果转换失败，使用序号 | If conversion fails, use index
if (any(is.na(gwas$CHR_numeric))) {
    unique_chrs <- unique(gwas$CHR)
    chr_mapping <- setNames(1:length(unique_chrs), unique_chrs)
    gwas$CHR_numeric <- chr_mapping[gwas$CHR]
}

# 准备绘图数据 | Prepare plotting data
gwas_plot <- data.frame(
    SNP = gwas$SNP,
    CHR = gwas$CHR_numeric,
    BP = gwas$BP,
    P = gwas$P
)

# 移除无效数据 | Remove invalid data
gwas_plot <- gwas_plot[!is.na(gwas_plot$P) & gwas_plot$P > 0 & gwas_plot$P <= 1, ]

if (nrow(gwas_plot) == 0) {
    stop("错误: 没有有效的P值数据 | Error: No valid P-value data")
}

cat("有效数据点数 | Valid data points:", nrow(gwas_plot), "\\n")

# 生成Manhattan图 | Generate Manhattan plot
cat("生成Manhattan图... | Generating Manhattan plot...\\n")
png("manhattan_plot.png", width=1200, height=600, res=100)
tryCatch({
    manhattan(gwas_plot, 
              main="PLINK GWAS Manhattan Plot (Disease Resistance)", 
              genomewideline = -log10(5e-8),
              suggestiveline = -log10(1e-5),
              col = c("blue4", "orange3"))
}, error = function(e) {
    plot(1, type="n", main="Manhattan Plot Error", xlab="Position", ylab="-log10(P)")
    text(1, 1, paste("Error:", e$message), cex=0.8)
})
dev.off()

# 生成QQ图 | Generate QQ plot
cat("生成QQ图... | Generating QQ plot...\\n")
png("qq_plot.png", width=600, height=600, res=100)
tryCatch({
    qq(gwas_plot$P, main="QQ Plot")
}, error = function(e) {
    plot(1, type="n", main="QQ Plot Error", xlab="Expected", ylab="Observed")
    text(1, 1, paste("Error:", e$message), cex=0.8)
})
dev.off()

# 保存染色体映射 | Save chromosome mapping
write.table(data.frame(
    Original_CHR = unique(gwas$CHR_original),
    Numeric_CHR = unique(gwas$CHR_numeric)
), "chromosome_mapping.txt", row.names=FALSE, quote=FALSE, sep="\\t")

cat("\\n=== 可视化完成 | Visualization completed ===\\n")
cat("生成的文件 | Generated files:\\n")
cat("- manhattan_plot.png\\n") 
cat("- qq_plot.png\\n")
cat("- chromosome_mapping.txt\\n")

cat("\\n=== 数据统计 | Data statistics ===\\n")
cat("总SNP数 | Total SNPs:", nrow(gwas_plot), "\\n")
cat("染色体数 | Number of chromosomes:", length(unique(gwas_plot$CHR)), "\\n")
cat("最小P值 | Minimum P-value:", min(gwas_plot$P), "\\n")
if (min(gwas_plot$P) < 5e-8) {
    cat("显著位点数 | Significant loci (P<5e-8):", sum(gwas_plot$P < 5e-8), "\\n")
}
cat("提示性位点数 | Suggestive loci (P<1e-5):", sum(gwas_plot$P < 1e-5), "\\n")
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
