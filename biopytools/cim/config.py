"""
CIM分析配置管理模块|CIM Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class CIMConfig:
    """CIM分析配置类|CIM Analysis Configuration Class"""

    # 必需参数|Required parameters
    input_file: str  # VCF文件路径|VCF file path
    pheno_file: str  # 表型文件路径|Phenotype file path
    output_dir: str  # 输出目录|Output directory

    # 群体类型参数|Cross type parameters
    cross_type: str = "f2"  # 群体类型|Cross type (f2/bc)

    # 遗传图谱参数|Genetic map parameters
    map_mode: str = "mstmap"  # cM来源模式|cM source mode (physical/estimate/mstmap)

    # MSTmap参数|MSTmap parameters
    mstmap_pvalue: float = 1e-6  # 聚类p值起始值(自动调优)|Clustering p-value start value (auto-tuned)
    mstmap_distfun: str = "kosambi"  # 距离函数|Distance function (kosambi/haldane)
    mstmap_path: str = "~/miniforge3/envs/Rqtl/bin/mstmap"  # MSTmap二进制路径|MSTmap binary path

    # 标记过滤参数|Marker filtering parameters
    maf: float = 0.05  # 最小等位基因频率|Minor allele frequency threshold
    max_missing: float = 0.1  # 最大缺失率|Maximum missing rate

    # CIM参数|CIM parameters
    n_marcovar: int = 10  # 协因子数量|Number of marker covariates
    window: float = 10.0  # 窗口大小(cM)|Window size (cM)
    method: str = "hk"  # 扫描方法|Scanning method (hk/em/imp)
    step: float = 1.0  # 伪标记步长(cM)|Pseudomarker step (cM)

    # 置换检验参数|Permutation test parameters
    n_perm: int = 1000  # 置换次数|Number of permutations

    # RF质控参数|RF QC parameters
    max_het_rate: float = 0.6  # H基因型最大比例|Max heterozygous genotype rate
    # K近邻局部平均RF最大值(规则2, 由全染色体平均重构: 全局均值在正确RF定义下
    # 无区分力, 远距离对≈0.5主导使正常标记也≈0.45-0.5)|Max local mean RF to the
    # k//2 nearest markers (rule 2, rebuilt from chromosome-wide mean: the global
    # mean cannot discriminate - distant pairs ~0.5 push good markers to ~0.45-0.5)
    max_mean_rf: float = 0.35
    rf_knn: int = 10  # 局部RF窗口邻居数(前后各rf_knn//2)|Local RF window neighbors

    # 局部重组热点过滤(双亲群体大片交换,短距离高频重组位点删除)|Local recombination
    # hotspot filter (large-block exchange: short-distance high-RF pairs are errors)
    enable_local_hotspot: bool = True   # 总开关|master switch
    local_hotspot_dist: int = 1000      # 物理距离阈值bp: 距离更近的相邻对才检查|only pairs closer than this
    local_hotspot_rf: float = 0.20      # 相邻RF软阈值: 达到即给两侧标记各记1分|soft RF (scores both neighbors)
    local_hotspot_hard_rf: float = 0.30  # 硬阈值: 单对RF达到即两侧标记都删|hard RF (delete both)
    local_hotspot_score_cut: int = 1    # 标记级删除线: 软评分>=此值删除|delete markers with score >= this
    local_hotspot_relative: float = 0.0  # 相对判据系数: >0时阈值=max(绝对值,系数*近距对RF中位数)|relative floor
    local_hotspot_min_markers: int = 100  # 只处理标记数>=此值的染色体(小scaffold跳过)|skip tiny scaffolds

    # 物理距离降采样|Physical distance thinning
    min_phys_gap: int = 0               # 最小相邻物理距离bp; 0=关闭|Min physical gap bp; 0=off

    # LD pruning参数|LD pruning parameters
    ld_window: int = 50  # LD计算窗口(SNP数)|LD window size (SNP count)
    ld_step: int = 5  # LD计算步长(SNP数)|LD step size (SNP count)
    ld_r2: float = 0.1  # LD r2阈值|LD r2 threshold
    skip_ld: bool = True  # 跳过LD降维（默认跳过，输入数据通常已预降维）|Skip LD pruning (default: input usually pre-pruned)

    # 环境参数|Environment parameters
    r_env: str = "~/miniforge3/envs/Rqtl"  # R conda环境路径或名称|R conda env path or name
    threads: int = 1  # 并行线程数|Number of parallel threads

    # 内部路径(自动生成)|Internal paths (auto-generated)
    qc_dir: str = field(init=False, default="")
    pre_rf_dir: str = field(init=False, default="")
    post_rf_dir: str = field(init=False, default="")
    pre_rf_tidy_dir: str = field(init=False, default="")
    post_rf_tidy_dir: str = field(init=False, default="")
    pre_rf_physical_dir: str = field(init=False, default="")
    pre_rf_mstmap_dir: str = field(init=False, default="")
    pre_rf_plots_dir: str = field(init=False, default="")
    post_rf_physical_dir: str = field(init=False, default="")
    post_rf_mstmap_dir: str = field(init=False, default="")
    post_rf_plots_dir: str = field(init=False, default="")
    r_script_dir: str = field(init=False, default="")
    log_dir: str = field(init=False, default="")
    info_dir: str = field(init=False, default="")
    log_file: str = field(init=False, default="")

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开MSTmap路径中的~|Expand ~ in mstmap_path
        self.mstmap_path = expand_path(self.mstmap_path)
        # 展开R环境路径|Expand R env path
        self.r_env = expand_path(self.r_env)

        # 规范化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.pheno_file = os.path.normpath(os.path.abspath(self.pheno_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 创建输出子目录|Create output subdirectories
        self.info_dir = os.path.join(self.output_dir, "00_pipeline_info")
        self.qc_dir = os.path.join(self.output_dir, "01_qc")
        self.pre_rf_dir = os.path.join(self.output_dir, "02_cim", "pre_rf")
        self.post_rf_dir = os.path.join(self.output_dir, "02_cim", "post_rf")
        self.pre_rf_tidy_dir = os.path.join(self.pre_rf_dir, "tidy_files")
        self.post_rf_tidy_dir = os.path.join(self.post_rf_dir, "tidy_files")
        self.pre_rf_physical_dir = os.path.join(self.pre_rf_dir, "physical")
        self.pre_rf_mstmap_dir = os.path.join(self.pre_rf_dir, "mstmap")
        self.pre_rf_plots_dir = os.path.join(self.pre_rf_dir, "plots")
        self.post_rf_physical_dir = os.path.join(self.post_rf_dir, "physical")
        self.post_rf_mstmap_dir = os.path.join(self.post_rf_dir, "mstmap")
        self.post_rf_plots_dir = os.path.join(self.post_rf_dir, "plots")
        self.r_script_dir = self.post_rf_dir  # 默认指向post-rf（向后兼容）
        self.log_dir = os.path.join(self.output_dir, "99_logs")

        for d in [self.info_dir, self.qc_dir,
                  self.pre_rf_dir, self.pre_rf_tidy_dir,
                  self.pre_rf_physical_dir, self.pre_rf_mstmap_dir, self.pre_rf_plots_dir,
                  self.post_rf_dir, self.post_rf_tidy_dir,
                  self.post_rf_physical_dir, self.post_rf_mstmap_dir, self.post_rf_plots_dir,
                  self.log_dir]:
            Path(d).mkdir(parents=True, exist_ok=True)

        self.log_file = os.path.join(self.log_dir, "cim.log")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        if not os.path.exists(self.pheno_file):
            errors.append(f"表型文件不存在|Phenotype file does not exist: {self.pheno_file}")

        # 检查群体类型|Check cross type
        valid_cross_types = ["f2", "bc"]
        if self.cross_type not in valid_cross_types:
            errors.append(f"群体类型无效|Invalid cross type: {self.cross_type}. "
                          f"支持|Supported: {', '.join(valid_cross_types)}")

        # 检查cM模式|Check map mode
        valid_map_modes = ["physical", "estimate", "mstmap"]
        if self.map_mode not in valid_map_modes:
            errors.append(f"cM模式无效|Invalid map mode: {self.map_mode}. "
                          f"支持|Supported: {', '.join(valid_map_modes)}")

        # 检查CIM方法|Check CIM method
        valid_methods = ["hk", "em", "imp"]
        if self.method not in valid_methods:
            errors.append(f"扫描方法无效|Invalid method: {self.method}. "
                          f"支持|Supported: {', '.join(valid_methods)}")

        # 检查参数范围|Check parameter ranges
        if not (0 < self.maf <= 0.5):
            errors.append(f"MAF阈值应在0到0.5之间|MAF threshold should be between 0 and 0.5: {self.maf}")

        if not (0 < self.max_missing <= 1.0):
            errors.append(f"最大缺失率应在0到1之间|Max missing rate should be between 0 and 1: {self.max_missing}")

        if not (0 < self.max_het_rate <= 1.0):
            errors.append(f"H基因型比例阈值应在0到1之间|Max het rate should be between 0 and 1: {self.max_het_rate}")

        if not (0 < self.max_mean_rf <= 1.0):
            errors.append(f"平均RF阈值应在0到1之间|Max mean RF should be between 0 and 1: {self.max_mean_rf}")

        if self.rf_knn < 2:
            errors.append(f"局部RF窗口邻居数应>=2|RF KNN window must be >= 2: {self.rf_knn}")

        # 检查局部热点过滤参数|Check local hotspot filter params
        if self.local_hotspot_dist <= 0:
            errors.append(f"热点距离阈值应>0|Hotspot distance threshold must be > 0: {self.local_hotspot_dist}")
        if not (0 < self.local_hotspot_rf <= 0.5):
            errors.append(f"热点软RF阈值应在(0,0.5]|Hotspot soft RF must be in (0,0.5]: {self.local_hotspot_rf}")
        if not (0 < self.local_hotspot_hard_rf <= 0.5):
            errors.append(f"热点硬RF阈值应在(0,0.5]|Hotspot hard RF must be in (0,0.5]: {self.local_hotspot_hard_rf}")
        if self.local_hotspot_hard_rf < self.local_hotspot_rf:
            errors.append(f"热点硬RF阈值应>=软RF阈值|Hotspot hard RF must be >= soft RF: "
                          f"{self.local_hotspot_hard_rf} < {self.local_hotspot_rf}")
        if self.local_hotspot_score_cut < 1:
            errors.append(f"热点评分删除线应>=1|Hotspot score cutoff must be >= 1: {self.local_hotspot_score_cut}")
        if self.local_hotspot_relative < 0:
            errors.append(f"热点相对判据系数应>=0|Hotspot relative factor must be >= 0: {self.local_hotspot_relative}")
        if self.local_hotspot_min_markers < 2:
            errors.append(f"热点过滤最小标记数应>=2|Hotspot min markers must be >= 2: {self.local_hotspot_min_markers}")
        if self.min_phys_gap < 0:
            errors.append(f"最小物理距离应>=0|Min physical gap must be >= 0: {self.min_phys_gap}")

        if self.n_marcovar < 1:
            errors.append(f"协因子数量应>=1|Number of covariates should be >= 1: {self.n_marcovar}")

        if self.window <= 0:
            errors.append(f"窗口大小应>0|Window size should be > 0: {self.window}")

        if self.step <= 0:
            errors.append(f"步长应>0|Step should be > 0: {self.step}")

        if self.n_perm < 0:
            errors.append(f"置换次数应>=0|Number of permutations should be >= 0: {self.n_perm}")

        if self.mstmap_distfun not in ["kosambi", "haldane"]:
            errors.append(f"MSTmap距离函数无效|Invalid MSTmap distance function: {self.mstmap_distfun}. "
                          f"支持|Supported: kosambi, haldane")

        if errors:
            raise ValueError("\n".join(errors))

        return True
