"""
R脚本生成器|R Script Generator

生成rMVP分析的R脚本|Generate R scripts for rMVP analysis
"""

from pathlib import Path
from typing import List


class RMVPRScriptGenerator:
    """rMVP R脚本生成器|rMVP R Script Generator"""

    def __init__(self, config):
        """
        初始化|Initialize

        Args:
            config: RMVPConfig对象|RMVPConfig object
        """
        self.config = config
        self.output_dir = Path(config.output_dir)

    def generate_data_conversion_script(self) -> str:
        """
        生成数据转换脚本|Generate data conversion script

        将VCF转换为rMVP格式|Convert VCF to rMVP format
        - LD去连锁开启时：完整VCF仅生成基因型（全部SNP），去连锁VCF计算K/PCA
        - LD pruning enabled: full VCF -> genotype only (all SNPs), pruned VCF -> K/PCA

        Returns:
            R脚本内容|R script content
        """
        # 使用绝对路径|Use absolute paths
        vcf_abs = str(Path(self.config.vcf_file).resolve())
        phe_abs = str(Path(self.config.pheno_file).resolve())
        out_abs = str(self.output_dir.resolve())

        if self.config.ld_pruning:
            # 去连锁VCF路径（由PLINK run_ld_pruning产出）|Pruned VCF path (produced by PLINK run_ld_pruning)
            pruned_vcf_abs = str((self.output_dir / f"{self.config.output_prefix}_pruned.vcf").resolve())
            mvp_data_block = f"""# 完整VCF → 基因型（全部SNP，不计算K/PC）|Full VCF -> genotype (all SNPs, no K/PC)
cat("转换完整VCF（全部SNP）|Converting full VCF (all SNPs)...\\n")
mvp_data <- MVP.Data(
    fileVCF = "{vcf_abs}",
    filePhe = "{phe_abs}",
    fileKin = FALSE,
    filePC = FALSE,
    out = file.path(out_path, output_prefix),
    ncpus = {self.config.ncpus}
)

# 去连锁VCF → K + PCA（kinship/PCA在去连锁SNP上计算）|Pruned VCF -> K + PCA (on LD-pruned SNPs)
cat("转换去连锁VCF（计算K/PCA）|Converting pruned VCF (computing K/PCA)...\\n")
mvp_data_pruned <- MVP.Data(
    fileVCF = "{pruned_vcf_abs}",
    filePhe = "{phe_abs}",
    fileKin = TRUE,
    filePC = TRUE,
    out = file.path(out_path, paste0(output_prefix, "_pruned")),
    ncpus = {self.config.ncpus}
)"""
            output_files_cat = """cat("  -", file.path(out_path, paste0(output_prefix, ".geno.desc")), " (全部SNP|all SNPs)\\n")
cat("  -", file.path(out_path, paste0(output_prefix, ".phe")), "\\n")
cat("  -", file.path(out_path, paste0(output_prefix, "_pruned.kin.desc")), " (去连锁|pruned)\\n")
cat("  -", file.path(out_path, paste0(output_prefix, "_pruned.pc.desc")), " (去连锁|pruned)\\n")"""
        else:
            mvp_data_block = f"""cat("正在转换VCF文件|Converting VCF file...\\n")
mvp_data <- MVP.Data(
    fileVCF = "{vcf_abs}",
    filePhe = "{phe_abs}",
    fileKin = TRUE,
    filePC = TRUE,
    out = file.path(out_path, output_prefix),
    ncpus = {self.config.ncpus}
)"""
            output_files_cat = """cat("  -", file.path(out_path, paste0(output_prefix, ".geno.desc")), "\\n")
cat("  -", file.path(out_path, paste0(output_prefix, ".phe")), "\\n")
cat("  -", file.path(out_path, paste0(output_prefix, ".kin.desc")), "\\n")
cat("  -", file.path(out_path, paste0(output_prefix, ".pc.desc")), "\\n")"""

        script = f"""# rMVP数据转换脚本|rMVP Data Conversion Script
# 自动生成于|Auto-generated at: {self._get_timestamp()}

# 加载rMVP包|Load rMVP package
suppressMessages({{
    library(rMVP)
}})

# 设置参数|Set parameters
file_vcf <- "{vcf_abs}"
file_phe <- "{phe_abs}"
output_prefix <- "{self.config.output_prefix}"
out_path <- "{out_abs}"

cat("\\n=== rMVP数据转换|rMVP Data Conversion ===\\n")
cat("VCF文件|VCF file:", file_vcf, "\\n")
cat("表型文件|Phenotype file:", file_phe, "\\n")
cat("输出前缀|Output prefix:", output_prefix, "\\n")
cat("输出目录|Output directory:", out_path, "\\n")
cat("LD去连锁|LD pruning:", {"TRUE" if self.config.ld_pruning else "FALSE"}, "\\n\\n")

# 转换VCF为rMVP格式|Convert VCF to rMVP format
{mvp_data_block}

cat("\\n数据转换完成|Data conversion completed\\n")
cat("输出文件|Output files:\\n")
{output_files_cat}
"""
        return script

    def _build_kinship_pca_blocks(self):
        """
        构建K/PCA读取与协方差参数的R代码块|Build K/PCA reading and covariate R blocks

        LD去连锁开启时：K与PCA从去连锁SNP读取，PCA作为CV.*协变量传入（避免内部双重PCA）
        LD pruning on: K and PCA read from pruned SNPs, PCA passed as CV.* covariates (no internal double PCA)
        关闭时：保持原逻辑（K来自完整VCF，nPC.*内部PCA）

        Returns:
            (kin_pc_block, covariate_args)|R代码块元组
        """
        if self.config.ld_pruning:
            kin_pc_block = f'''# 读取去连锁K和PCA（kinship/PCA来自去连锁SNP）|Read pruned K and PCA (from LD-pruned SNPs)
cat("正在读取去连锁Kinship矩阵|Reading pruned Kinship matrix...\\n")
K <- attach.big.matrix(file.path(out_path, paste0(output_prefix, "_pruned.kin.desc")))

cat("正在读取去连锁PCA|Reading pruned PCA...\\n")
pc <- bigmemory::as.matrix(attach.big.matrix(file.path(out_path, paste0(output_prefix, "_pruned.pc.desc"))))

# 取前n个PC作为协变量|Take first n PCs as covariates
pc_glm <- pc[, seq_len(min({self.config.n_pc_glm}, ncol(pc))), drop = FALSE]
pc_mlm <- pc[, seq_len(min({self.config.n_pc_mlm}, ncol(pc))), drop = FALSE]
pc_farmcpu <- pc[, seq_len(min({self.config.n_pc_farmcpu}, ncol(pc))), drop = FALSE]'''
            covariate_args = '''    CV.GLM = pc_glm,
    CV.MLM = pc_mlm,
    CV.FarmCPU = pc_farmcpu,'''
        else:
            kin_pc_block = '''# 读取Kinship和PCA|Read Kinship and PCA
if (file.exists(file.path(out_path, paste0(output_prefix, ".kin.desc")))) {
    cat("正在读取Kinship矩阵|Reading Kinship matrix...\\n")
    K <- attach.big.matrix(file.path(out_path, paste0(output_prefix, ".kin.desc")))
} else {
    K <- NULL
    cat("未使用Kinship矩阵|Kinship matrix not used\\n")
}

if (file.exists(file.path(out_path, paste0(output_prefix, ".pc.desc")))) {
    cat("正在读取PCA|Reading PCA...\\n")
    pc <- bigmemory::as.matrix(attach.big.matrix(file.path(out_path, paste0(output_prefix, ".pc.desc"))))
} else {
    pc <- NULL
    cat("未使用PCA|PCA not used\\n")
}'''
            covariate_args = f'''    nPC.GLM = {self.config.n_pc_glm},
    nPC.MLM = {self.config.n_pc_mlm},
    nPC.FarmCPU = {self.config.n_pc_farmcpu},'''
        return kin_pc_block, covariate_args

    def generate_gwas_script(self, trait_name: str, trait_col: int) -> str:
        """
        生成单个表型的GWAS分析脚本|Generate GWAS script for single trait

        Args:
            trait_name: 表型名称|Trait name
            trait_col: 表型列号（从2开始）|Trait column number (starts from 2)

        Returns:
            R脚本内容|R script content
        """
        # 生成输出文件名|Generate output file names
        trait_prefix = f"{self.config.output_prefix}_{trait_name}"

        # 生成模型列表|Generate model list
        models_str = ', '.join([f'"{m}"' for m in self.config.models])

        # 使用绝对路径|Use absolute paths
        out_abs = str(self.output_dir.resolve())

        # K/PCA来源与协方差参数（根据LD去连锁分支）|K/PCA source and covariate args (branched on LD pruning)
        kin_pc_block, covariate_args = self._build_kinship_pca_blocks()

        script = f"""# rMVP GWAS分析脚本|rMVP GWAS Analysis Script
# 表型|Trait: {trait_name}
# 自动生成于|Auto-generated at: {self._get_timestamp()}

# 加载rMVP包|Load rMVP package
suppressMessages({{
    library(rMVP)
    library(bigmemory)
}})

# 设置参数|Set parameters
output_prefix <- "{self.config.output_prefix}"
trait_name <- "{trait_name}"
trait_col <- {trait_col}
trait_prefix <- "{trait_prefix}"
out_path <- "{out_abs}"
ncpus <- {self.config.ncpus}
maxLine <- {self.config.maxLine}

cat("\\n=== rMVP GWAS分析|rMVP GWAS Analysis ===\\n")
cat("表型名称|Trait name:", trait_name, "\\n")
cat("表型列号|Trait column:", trait_col, "\\n")
cat("分析模型|Models:", {models_str}, "\\n")
cat("CPU核心数|CPU cores:", ncpus, "\\n")
cat("输出目录|Output directory:", out_path, "\\n\\n")

# 读取数据|Read data
cat("正在读取数据|Reading data...\\n")

genotype <- attach.big.matrix(file.path(out_path, paste0(output_prefix, ".geno.desc")))
phenotype_full <- read.table(file.path(out_path, paste0(output_prefix, ".phe")), header = TRUE)
map <- read.table(file.path(out_path, paste0(output_prefix, ".geno.map")), header = TRUE)

# 提取单个表型|Extract single trait
phenotype <- phenotype_full[, c(1, trait_col)]
colnames(phenotype)[2] <- trait_name

cat("样本数|Sample size:", nrow(phenotype), "\\n")
cat("SNP数|SNP count:", ncol(genotype), "\\n")

{kin_pc_block}

# 运行GWAS|Run GWAS
cat("\\n开始GWAS分析|Starting GWAS analysis...\\n")

start_time <- Sys.time()

imvp <- MVP(
    phe = phenotype,
    geno = genotype,
    map = map,
    K = K,
{covariate_args}
    maxLine = maxLine,
    ncpus = ncpus,
    vc.method = "{self.config.vc_method}",
    method.bin = "{self.config.method_bin}",
    maxLoop = {self.config.max_loop},
    threshold = {self.config.threshold},
    method = c({models_str}),
    file.output = c("pmap", "pmap.signal", "plot", "log"),
    file.type = "{self.config.file_type}",
    dpi = {self.config.dpi},
    memo = trait_name,
    outpath = out_path,
    verbose = TRUE
)

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "secs")

cat("\\n分析完成|Analysis completed\\n")
cat("运行时间|Runtime:", round(as.numeric(duration), 2), "秒|seconds\\n")

# 保存结果到RData文件|Save results to RData file
result_file <- file.path(out_path, paste0(trait_prefix, ".RData"))
save(imvp, file = result_file)
cat("结果已保存|Results saved to:", result_file, "\\n")
"""
        return script

    def generate_batch_script(self, trait_names: List[str]) -> str:
        """
        生成批量分析脚本（一次运行所有表型）|Generate batch analysis script

        Args:
            trait_names: 表型名称列表|List of trait names

        Returns:
            R脚本内容|R script content
        """
        models_str = ', '.join([f'"{m}"' for m in self.config.models])

        # 使用绝对路径|Use absolute paths
        out_abs = str(self.output_dir.resolve())

        # K/PCA来源与协方差参数（根据LD去连锁分支）|K/PCA source and covariate args (branched on LD pruning)
        kin_pc_block, covariate_args = self._build_kinship_pca_blocks()

        script = f"""# rMVP批量GWAS分析脚本|rMVP Batch GWAS Analysis Script
# 分析{len(trait_names)}个表型|Analyzing {len(trait_names)} traits
# 自动生成于|Auto-generated at: {self._get_timestamp()}

# 加载rMVP包|Load rMVP package
suppressMessages({{
    library(rMVP)
    library(bigmemory)
}})

# 设置参数|Set parameters
output_prefix <- "{self.config.output_prefix}"
out_path <- "{out_abs}"
ncpus <- {self.config.ncpus}
maxLine <- {self.config.maxLine}

cat("\\n=== rMVP批量GWAS分析|rMVP Batch GWAS Analysis ===\\n")
cat("表型数量|Number of traits:", {len(trait_names)}, "\\n")
cat("分析模型|Models:", {models_str}, "\\n")
cat("CPU核心数|CPU cores:", ncpus, "\\n")
cat("输出目录|Output directory:", out_path, "\\n\\n")

# 读取数据|Read data
cat("正在读取数据|Reading data...\\n")

genotype <- attach.big.matrix(file.path(out_path, paste0(output_prefix, ".geno.desc")))
phenotype_full <- read.table(file.path(out_path, paste0(output_prefix, ".phe")), header = TRUE)
map <- read.table(file.path(out_path, paste0(output_prefix, ".geno.map")), header = TRUE)

{kin_pc_block}

# 批量分析所有表型|Batch analyze all traits
trait_names <- c({', '.join([f'"{t}"' for t in trait_names])})
n_traits <- length(trait_names)

all_results <- list()

for (i in 1:n_traits) {{
    trait_name <- trait_names[i]
    trait_col <- i + 1  # 表型从第2列开始|Traits start from column 2

    cat("\\n")
    cat("=" , rep("=", 60), "\\n", sep = "")
    cat("正在分析表型", i, "/", n_traits, ":", trait_name, "\\n")
    cat("Analyzing trait", i, "/", n_traits, ":", trait_name, "\\n")
    cat("=" , rep("=", 60), "\\n", sep = "")

    # 提取单个表型|Extract single trait
    phenotype <- phenotype_full[, c(1, trait_col)]
    colnames(phenotype)[2] <- trait_name

    # 运行GWAS|Run GWAS
    start_time <- Sys.time()

    imvp <- tryCatch({{
        MVP(
            phe = phenotype,
            geno = genotype,
            map = map,
            K = K,
{covariate_args}
            maxLine = maxLine,
            ncpus = ncpus,
            vc.method = "{self.config.vc_method}",
            method.bin = "{self.config.method_bin}",
            maxLoop = {self.config.max_loop},
            threshold = {self.config.threshold},
            method = c({models_str}),
            file.output = c("pmap", "pmap.signal", "plot", "log"),
            file.type = "{self.config.file_type}",
            dpi = {self.config.dpi},
            memo = trait_name,
            outpath = out_path,
            verbose = TRUE
        )
    }}, error = function(e) {{
        cat("错误|Error:", e$message, "\\n")
        NULL
    }})

    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "secs")

    if (!is.null(imvp)) {{
        cat("\\n表型", trait_name, "分析完成|Trait", trait_name, "analysis completed\\n")
        cat("运行时间|Runtime:", round(as.numeric(duration), 2), "秒|seconds\\n")

        # 保存结果|Save results
        trait_prefix <- paste0(output_prefix, "_", trait_name)
        result_file <- file.path(out_path, paste0(trait_prefix, ".RData"))
        save(imvp, file = result_file)

        all_results[[trait_name]] <- imvp
    }} else {{
        cat("\\n表型", trait_name, "分析失败|Trait", trait_name, "analysis failed\\n")
    }}

    # 清理内存|Clean memory
    gc()
}}

cat("\\n")
cat("=" , rep("=", 60), "\\n", sep = "")
cat("所有表型分析完成|All traits analysis completed\\n")
cat("=" , rep("=", 60), "\\n", sep = "")

# 保存所有结果|Save all results
all_results_file <- file.path(out_path, paste0(output_prefix, "_all_results.RData"))
save(all_results, file = all_results_file)
cat("\\n所有结果已保存|All results saved to:", all_results_file, "\\n")
"""
        return script

    def _get_timestamp(self) -> str:
        """获取当前时间戳|Get current timestamp"""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
