"""
Pixy配置管理模块|Pixy Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class PixyConfig:
    """Pixy配置类|Pixy Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str = ""  # 输入VCF文件（需用bgzip压缩并建立tabix索引）|Input VCF file (must be bgzip-compressed and tabix-indexed)
    pop_file: str = ""  # 群体文件（两列：样本ID 群体名）|Population file (two columns: sample_id population_name)
    output_dir: str = ""  # 输出目录|Output directory

    # 统计量选择|Statistics selection
    calc_pi: bool = True  # 计算pi（核苷酸多样性）|Calculate pi (nucleotide diversity)
    calc_dxy: bool = True  # 计算dxy（群体间核苷酸差异）|Calculate dxy (nucleotide divergence between populations)
    calc_fst: bool = True  # 计算fst（遗传分化系数）|Calculate fst (genetic differentiation coefficient)
    calc_watterson_theta: bool = True  # 计算Watterson's theta|Calculate Watterson's theta
    calc_tajima_d: bool = True  # 计算Tajima's D|Calculate Tajima's D

    # 窗口参数|Window parameters
    window_size: Optional[int] = None  # 窗口大小（bp），pixy要求必须指定窗口大小或BED文件|Window size in bp (pixy requires window_size or bed_file)
    bed_file: Optional[str] = None  # BED文件定义窗口（自定义大小窗口）|BED file defining windows (custom-sized windows)
    sites_file: Optional[str] = None  # 位点文件（只计算特定位点）|Sites file (calculate only specific sites)

    # 质控参数|Quality control parameters
    min_samples: int = 0  # 每个群体最小样本数（0表示不限制）|Minimum samples per population (0=no limit)
    max_missing: float = 1.0  # 最大缺失率（默认1.0=不限制）|Maximum missing rate (default 1.0=no limit)
    min_maf: float = 0.0  # 最小等位基因频率（默认0.0=不限制）|Minor allele frequency (default 0.0=no limit)
    zscore_window: Optional[int] = None  # Z-score过滤窗口大小（不设置则不过滤）|Z-score filtering window size (null=no filter)

    # 染色体参数|Chromosome parameters
    chromosomes: Optional[List[str]] = None  # 指定染色体列表（不设置则全部）|List of chromosomes (null for all)

    # 工具路径|Tool paths
    pixy_path: str = "pixy"  # pixy可执行文件路径|pixy executable path
    conda_env: str = "~/miniforge3/envs/pixy_v.2.0.0"  # conda环境路径|conda environment path

    # 其他参数|Other parameters
    threads: int = 12  # 线程数|Number of threads
    keep_intermediate: bool = False  # 保留中间文件|Keep intermediate files
    verbose: bool = False  # 详细输出|Verbose output
    bypass_invariant_check: bool = False  # 绕过不变位点检查（默认自动检测并自动启用）|Bypass invariant sites check (default: auto-detect and auto-enable)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 路径展开|Path expansion
        self.conda_env = os.path.expanduser(self.conda_env)
        if self.vcf_file:
            self.vcf_path = Path(self.vcf_file)
        if self.pop_file:
            self.pop_path = Path(self.pop_file)
        if self.output_dir:
            self.output_path = Path(self.output_dir).resolve()
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = str(self.output_path)
        if self.bed_file:
            self.bed_path = Path(self.bed_file)
        else:
            self.bed_path = None
        if self.sites_file:
            self.sites_path = Path(self.sites_file)
        else:
            self.sites_path = None

        # 染色体列表|Chromosome list
        if self.chromosomes:
            self.chromosome_list = self.chromosomes
        else:
            self.chromosome_list = None

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查VCF文件|Check VCF file
        if not self.vcf_file:
            errors.append("VCF文件不能为空|VCF file cannot be empty")
        elif not self.vcf_path.exists():
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        # 检查VCF索引文件|Check VCF index file
        if self.vcf_file and self.vcf_path.exists():
            # 索引文件应该与VCF文件同名，加上.csi或.tbi后缀
            # Index file should have the same name as VCF file, with .csi or .tbi suffix
            vcf_index = Path(str(self.vcf_path) + '.csi')
            if not vcf_index.exists():
                vcf_index_tbi = Path(str(self.vcf_path) + '.tbi')
                if not vcf_index_tbi.exists():
                    errors.append(f"VCF索引文件不存在（需要用tabix建立索引）|VCF index file does not exist (needs to be indexed with tabix): {self.vcf_path}.csi or {self.vcf_path}.tbi")

        # 检查群体文件|Check population file
        if not self.pop_file:
            errors.append("群体文件不能为空|Population file cannot be empty")
        elif not self.pop_path.exists():
            errors.append(f"群体文件不存在|Population file does not exist: {self.pop_file}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        # 检查BED文件|Check BED file
        if self.bed_file and not self.bed_path.exists():
            errors.append(f"BED文件不存在|BED file does not exist: {self.bed_file}")

        # 检查sites文件|Check sites file
        if self.sites_file and not self.sites_path.exists():
            errors.append(f"位点文件不存在|Sites file does not exist: {self.sites_file}")

        # 检查窗口参数互斥|Check mutual exclusion of window parameters
        window_params = [self.window_size is not None, self.bed_file is not None, self.sites_file is not None]
        if sum(window_params) > 1:
            errors.append("窗口大小、BED文件和位点文件只能指定其中一个|Can only specify one of window_size, bed_file, or sites_file")

        # pixy要求必须指定window_size或bed_file（sites_file除外）
        # pixy requires window_size or bed_file (except for sites_file)
        if not (self.window_size is not None or self.bed_file is not None or self.sites_file is not None):
            errors.append("pixy要求必须指定窗口大小(-w/--window-size)、BED文件(-b/--bed-file)或位点文件(-s/--sites-file)|pixy requires window_size (-w/--window-size), bed_file (-b/--bed-file), or sites_file (-s/--sites-file)")

        # 检查至少选择一个统计量|Check at least one statistic is selected
        if not any([self.calc_pi, self.calc_dxy, self.calc_fst, self.calc_watterson_theta, self.calc_tajima_d]):
            errors.append("至少需要选择一个统计量（pi、dxy、fst、watterson_theta、tajima_d）|At least one statistic must be selected (pi, dxy, fst, watterson_theta, tajima_d)")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查质控参数范围|Check QC parameter ranges
        if self.max_missing < 0 or self.max_missing > 1:
            errors.append(f"最大缺失率必须在0-1之间|max_missing must be between 0-1: {self.max_missing}")

        if self.min_maf < 0 or self.min_maf > 1:
            errors.append(f"最小等位基因频率必须在0-1之间|min_maf must be between 0-1: {self.min_maf}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
