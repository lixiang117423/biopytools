"""
 ADMIXTURE结果处理和可视化模块|ADMIXTURE Results Processing and Visualization Module

协变量生成(GWAS covariates)、R绘图、总结报告。所有产物写04_results/;
读取Q/P从03_admixture/、fam从02_plink/。R脚本不自动装包(超算无网)。
"""

import glob
import os
import shutil
import subprocess
from shutil import which

import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage


class CovariateGenerator:
    """ 协变量生成器|Covariate Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_gwas_covariates(self, best_k: int):
        """ 生成GWAS协变量文件|Generate GWAS covariate file (写04_results/)"""
        self.logger.info(f"生成GWAS协变量文件|Generating GWAS covariate file (K={best_k})")

        # 定位Q文件(03_admixture/,兼容ADAMIXTURE双格式)|Locate Q file
        base = self.config.base_name
        q_f1 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.Q")
        q_f2 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.{best_k}.Q")
        q_file = q_f1 if os.path.exists(q_f1) else q_f2
        if not os.path.exists(q_file):
            raise FileNotFoundError(
                f"Q文件不存在|Q file not found: tried {q_f1} and {q_f2}")

        q_data = pd.read_csv(q_file, sep=r'\s+', header=None)

        # 读取个体信息(从02_plink/的fam)|Read individual info from fam
        fam_file = os.path.join(self.config.plink_dir, f"{base}.fam")
        if not os.path.exists(fam_file):
            raise FileNotFoundError(f"FAM文件不存在|FAM file not found: {fam_file}")
        fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None)

        # 用前K-1个祖先成分作为协变量(避免共线性)|First K-1 components as covariates
        covariates = pd.DataFrame()
        covariates['FID'] = fam_data.iloc[:, 0].values
        covariates['IID'] = fam_data.iloc[:, 1].values
        for i in range(best_k - 1):
            covariates[f'PC{i+1}'] = q_data.iloc[:, i].values

        covar_file = os.path.join(self.config.results_dir, "gwas_covariates.txt")
        covariates.to_csv(covar_file, sep='\t', header=True, index=False)

        self.logger.info(f"GWAS协变量文件已保存|GWAS covariate file saved: {covar_file}")
        self.logger.info("在PLINK中使用|Use in PLINK: --covar gwas_covariates.txt --covar-name PC1,PC2,...")

        return covar_file


class PlotFilesCollector:
    """ 绘图文件收集器|Plot files collector (写05_plot_files/,供下游R绘图包)

    收拢下游R绘图包的两个输入:
    - 02_plink/{base}.fam -> 05_plot_files/admixture_ready.fam (复制改名,R包固定文件名)
    - 03_admixture/*.Q -> 05_plot_files/ (按原名复制,覆盖所有K)
    复制而非移动: 02/03 目录文件是断点续传依据; 文件极小,每次运行刷新,
    保证重新跑更宽K范围后新Q也会进目录。
    """

    # 下游R包要求的fam文件名|fam filename expected by downstream R package
    FAM_DST_NAME = "admixture_ready.fam"

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def collect(self) -> dict:
        """ 复制fam与Q文件到05_plot_files/|Copy fam and Q files into 05_plot_files/

        Returns:
            dict: {'fam': 目标fam路径|dest fam path, 'q_files': 目标Q路径列表|dest Q paths}
        """
        # dry_run: 不创建目录不复制|dry-run: create nothing, copy nothing
        if self.config.dry_run:
            self.logger.info(
                "模拟运行: 跳过绘图文件收集|Dry run: skipping plot files collection")
            return {'fam': os.path.join(self.config.plot_files_dir, self.FAM_DST_NAME),
                    'q_files': []}

        os.makedirs(self.config.plot_files_dir, exist_ok=True)

        # fam: 复制并改名|fam: copy and rename
        fam_src = os.path.join(self.config.plink_dir, f"{self.config.base_name}.fam")
        if not os.path.exists(fam_src):
            raise FileNotFoundError(f"FAM文件不存在|FAM file not found: {fam_src}")
        fam_dst = os.path.join(self.config.plot_files_dir, self.FAM_DST_NAME)
        shutil.copy2(fam_src, fam_dst)
        self.logger.info(f"FAM文件已复制|FAM file copied: {fam_src} -> {fam_dst}")

        # Q: 按原名复制全部K|Q: copy all K as-is
        q_srcs = sorted(glob.glob(os.path.join(self.config.admixture_dir, "*.Q")))
        q_dsts = []
        if not q_srcs:
            self.logger.warning(
                "未找到Q文件|No .Q files in 03_admixture/ (ADMIXTURE可能未产出任何K|analysis may have produced no K)")
        for q_src in q_srcs:
            q_dst = os.path.join(self.config.plot_files_dir, os.path.basename(q_src))
            shutil.copy2(q_src, q_dst)
            q_dsts.append(q_dst)
        if q_srcs:
            self.logger.info(
                f"Q文件已复制|Q files copied: {len(q_srcs)} 个到|files to {self.config.plot_files_dir}")

        return {'fam': fam_dst, 'q_files': q_dsts}


class ClusterAssignmentGenerator:
    """ 簇归属表生成器|Cluster Assignment Table Generator (写04_results/)

    每个K下逐样品取Q行最大比例的cluster(1-based)作为归属,附该最大比例值;
    混合样品照常归属,比例值低肉眼可见,阈值由用户自行筛选。
    注意: ADMIXTURE各K的cluster编号互相独立(K=3的Pop1与K=4的Pop1无对应关系)。
    |For each K, assign each sample to the cluster (1-based) with the max Q
    proportion, plus that max value. Admixed samples are still assigned — the
    low proportion is visible for manual filtering. Cluster numbering is
    independent across K runs (Pop1@K=3 is unrelated to Pop1@K=4).
    """

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate(self):
        """ 生成簇归属宽表|Generate cluster assignment wide table

        Returns:
            str | None: 输出CSV路径; 无任何可用K或dry_run时返回None
        """
        # dry_run: 不写文件|dry-run: write nothing
        if self.config.dry_run:
            self.logger.info(
                "模拟运行: 跳过簇归属表生成|Dry run: skipping cluster assignment generation")
            return None

        # fam提供FID/IID(行序=Q行序)|fam provides FID/IID (row order = Q row order)
        fam_file = os.path.join(self.config.plink_dir, f"{self.config.base_name}.fam")
        if not os.path.exists(fam_file):
            raise FileNotFoundError(f"FAM文件不存在|FAM file not found: {fam_file}")
        fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None)

        table = pd.DataFrame()
        table['FID'] = fam_data.iloc[:, 0].values
        table['IID'] = fam_data.iloc[:, 1].values

        found_ks = []
        for k in range(self.config.min_k, self.config.max_k + 1):
            q_file = self._locate_q_file(k)
            if q_file is None:
                self.logger.warning(
                    f"K={k} 的Q文件缺失,跳过该K|Q file missing for K={k}, skipping this K")
                continue

            q_data = pd.read_csv(q_file, sep=r'\s+', header=None)
            if len(q_data) != len(table):
                # 行数不齐时FID/IID会错位,必须跳过|misaligned FID/IID hazard, must skip
                self.logger.warning(
                    f"K={k} 的Q行数({len(q_data)})与fam行数({len(table)})不一致,跳过该K"
                    f"|Q rows ({len(q_data)}) != fam rows ({len(table)}), skipping K={k}")
                continue

            # argmax归属(1-based)+最大比例值+构成相似组|argmax assignment + max proportion + composition-similarity group
            table[f'K{k}'] = (q_data.idxmax(axis=1) + 1).values
            table[f'K{k}_max'] = q_data.max(axis=1).values
            table[f'K{k}_qgroup'] = self._qgroup_labels(q_data, k)
            found_ks.append(k)

        if not found_ks:
            self.logger.warning(
                "未找到任何可用的Q文件,跳过簇归属表|No usable Q files, skipping cluster assignment table")
            return None

        out_file = os.path.join(self.config.results_dir, "cluster_assignment.csv")
        table.to_csv(out_file, index=False)
        self.logger.info(
            f"簇归属表已保存|Cluster assignment table saved: {out_file} "
            f"(覆盖K|covers K: {','.join(str(k) for k in found_ks)})")
        return out_file

    def _qgroup_labels(self, q_data: pd.DataFrame, k: int):
        """ 祖先构成相似组|Ancestry-composition similarity groups

        对Q行(样品的K维构成向量)层次聚类(Euclidean + average linkage,
        确定性无随机种子),切成至多K组;按组大小降序重编号(1=最大组,
        并列按首成员序号破平)。注意: qgroup是"构成相似样品分组",
        与K{n}列的argmax祖先群体归属语义不同——混合个体多的面板里
        qgroup可能把杂交个体单独成组。
        |Hierarchical clustering of Q rows (Euclidean + average linkage,
        deterministic), cut into at most K groups; renumbered by size
        descending (1=largest, ties by first member index). NOTE: qgroup
        groups samples by composition similarity — distinct from the argmax
        ancestral-population assignment in the K{n} column.
        """
        n = len(q_data)
        if n < 2:
            # 单样品无法聚类,全体归组1|no clustering possible, all group 1
            return [1] * n

        z = linkage(q_data.values, method='average', metric='euclidean')
        raw_labels = fcluster(z, t=k, criterion='maxclust')

        # 重编号: 组大小降序,并列按首成员序号|renumber: size desc, ties by first member
        sizes = {}
        first_idx = {}
        for i, lab in enumerate(raw_labels):
            sizes[lab] = sizes.get(lab, 0) + 1
            first_idx.setdefault(lab, i)
        order = sorted(sizes, key=lambda lab: (-sizes[lab], first_idx[lab]))
        remap = {lab: rank + 1 for rank, lab in enumerate(order)}
        return [remap[lab] for lab in raw_labels]

    def _locate_q_file(self, k: int):
        """ 定位Q文件(03_admixture/,兼容双格式)|Locate Q file (dual-format compatible)

        ADMIXTURE产出{base}.{k}.Q; ADAMIXTURE产出{base}.{k}.{k}.Q。
        """
        base = self.config.base_name
        for name in (f"{base}.{k}.Q", f"{base}.{k}.{k}.Q"):
            path = os.path.join(self.config.admixture_dir, name)
            if os.path.exists(path):
                return path
        return None


class PlotGenerator:
    """ 绘图生成器|Plot Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_plots(self, q_data: pd.DataFrame, best_k: int):
        """ 生成可视化图表(写04_results/)|Generate visualization plots (to 04_results/)"""
        self.logger.info("正在生成可视化图表|Generating visualization plots...")

        if which("Rscript") is None:
            self.logger.warning("Rscript未安装或不在PATH中，跳过图表生成|Rscript not in PATH, skipping plot generation")
            return None

        # 生成并保存R脚本到04_results/|Generate R script into 04_results/
        r_script = self._create_r_script(best_k)
        r_file = os.path.join(self.config.results_dir, "generate_plots.R")
        with open(r_file, 'w', encoding='utf-8') as f:
            f.write(r_script)

        # 记录命令(§2.2.1)|Log command
        self.logger.info(f"命令|Command: Rscript {r_file}")

        try:
            result = subprocess.run(
                ["Rscript", r_file],
                cwd=self.config.results_dir,
                capture_output=True,
                text=True,
                check=True,
            )
            self.logger.info("R图表生成成功|R plots generated successfully")
            if result.stdout:
                self.logger.info(f"R脚本输出|R script output:\n{result.stdout}")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"R图表生成失败|R plot generation failed.")
            self.logger.warning(f"R脚本错误信息|R script error message:\n{e.stderr}")
        except FileNotFoundError:
            self.logger.warning("Rscript未安装，跳过图表生成|Rscript not installed, skipping plot generation")

        return r_file

    def _create_r_script(self, best_k: int):
        """ 创建R脚本|Create R script (绝对路径指向by-step目录,pdf写cwd=04_results/)"""
        base = self.config.base_name
        q_path = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.Q")
        fam_path = os.path.join(self.config.plink_dir, f"{base}.fam")
        cv_path = os.path.join(self.config.admixture_dir, "cv_results.csv")

        return f'''#!/usr/bin/env Rscript
#  ADMIXTURE结果可视化脚本|ADMIXTURE Results Visualization Script

# 检查必要R包(缺失则跳过,不在超算自动安装)|Check packages; skip if missing (no auto-install on HPC)
needed <- c("ggplot2", "reshape2", "dplyr")
missing_pkgs <- character(0)
for (pkg in needed) {{
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {{
        missing_pkgs <- c(missing_pkgs, pkg)
    }}
}}
if (length(missing_pkgs) > 0) {{
    cat("缺少R包,跳过绘图|Missing R packages, skipping plots:", paste(missing_pkgs, collapse = ", "), "\\n")
    quit(save = "no", status = 0)
}}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

# 读取数据(绝对路径)|Read data (absolute paths)
q_file_path <- "{q_path}"
fam_file_path <- "{fam_path}"
cv_file_path <- "{cv_path}"

if (!file.exists(q_file_path) || !file.exists(fam_file_path)) {{
    stop(" Q或FAM文件不存在|Q or FAM file not found.")
}}

q_data <- read.table(q_file_path)
colnames(q_data) <- paste0("Pop", 1:{best_k})

#  添加个体信息|Add individual info
fam_data <- read.table(fam_file_path)
q_data$Individual_ID <- 1:nrow(q_data)
q_data$FID <- fam_data$V1
q_data$IID <- fam_data$V2

# 1. 个体祖先成分条形图|Individual ancestry bar plot
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
hist(admixture_level, main = "Admixture Level Distribution",
     xlab = "Admixture Level (1 - Max Ancestry)", ylab = "Number of Individuals",
     col = "lightblue", breaks = 20)
hist(max_ancestry, main = "Maximum Ancestry Distribution",
     xlab = "Maximum Ancestry Proportion", ylab = "Number of Individuals",
     col = "lightgreen", breaks = 20)
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
        """ 生成分析总结(写04_results/)|Generate analysis summary (to 04_results/)"""
        summary_file = os.path.join(self.config.results_dir, "analysis_summary.txt")

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(" ADMIXTURE分析总结报告|ADMIXTURE Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")

            f.write(" 输入文件|Input Files:\n")
            f.write(f"  -  VCF文件|VCF file: {self.config.vcf_file}\n\n")

            f.write(" 分析参数|Analysis Parameters:\n")
            f.write(f"  -  方法|Method: {self.config.method}\n")
            f.write(f"  -  K值范围|K range: {self.config.min_k} - {self.config.max_k}\n")
            f.write(f"  -  交叉验证折数|CV folds: {self.config.cv_folds}\n")
            f.write(f"  -  线程数|Threads: {self.config.threads}\n")
            f.write(f"  -  MAF阈值|MAF threshold: {self.config.maf}\n")
            f.write(f"  -  缺失率阈值|Missing rate threshold: {self.config.missing_rate}\n")
            f.write(f"  -  HWE p值阈值|HWE p-value threshold: {self.config.hwe_pvalue}\n")
            f.write(f"  -  LD剪枝|LD pruning: {self.config.ld_prune} (window={self.config.ld_window}, step={self.config.ld_step}, r2={self.config.ld_r2})\n\n")

            f.write(" 分析结果|Analysis Results:\n")
            f.write(f"  -  最优K值|Best K value (lowest CV error): {best_k}\n")
            f.write(f"  -  总个体数|Total individuals: {stats['total_individuals']}\n")
            f.write(f"  -  高度混合个体数 (max ancestry < 0.7)|Highly admixed individuals: {stats['highly_admixed']}\n")
            f.write(f"  -  纯合个体数 (max ancestry > 0.9)|Pure individuals: {stats['pure_individuals']}\n")
            f.write(f"  -  平均混合程度|Mean admixture level: {stats['mean_admixture_level']:.4f}\n")
            f.write(f"  -  平均最大祖先成分|Mean max ancestry: {stats['mean_max_ancestry']:.4f}\n\n")

            f.write("输出文件|Output Files:\n")
            f.write(f"  -  个体祖先成分|Individual ancestry proportions: 04_results/admixture_proportions.csv\n")
            f.write(f"  -  簇归属表(每K每样品)|Cluster assignment (per K per sample): 04_results/cluster_assignment.csv\n")
            f.write(f"  -  GWAS协变量|GWAS covariates: 04_results/gwas_covariates.txt\n")
            f.write(f"  -  交叉验证结果|Cross-validation results: 03_admixture/cv_results.csv\n")
            f.write(f"  -  可视化图表|Visualization plots: 04_results/*.pdf\n")
            f.write(f"  -  分析总结|Analysis summary: 04_results/analysis_summary.txt\n")
            f.write(f"  -  下游R绘图输入|Downstream R plotting inputs: 05_plot_files/\n")
            f.write(f"  -  软件版本|Software versions: 00_pipeline_info/software_versions.yml\n")
            f.write(f"输出目录|Output directory: {self.config.output_dir}\n")

        self.logger.info(f"分析总结已保存|Analysis summary saved: {summary_file}")
        return summary_file
