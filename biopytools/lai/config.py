"""
LAI模块配置类|LAI Module Configuration Classes
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class LAIConfig:
    """LAI计算基础配置类|LAI Calculation Base Configuration Class"""

    # 必需参数|Required parameters
    genome: str = ""  # 基因组FASTA文件|Genome FASTA file
    output_dir: str = ""  # 输出目录|Output directory

    # 运行参数|Runtime parameters
    threads: int = 12  # 线程数|Number of threads
    mode: str = "full"  # 运行模式: full(完整流程), harvest(仅候选识别), retrieve(仅筛选), calculate(仅LAI计算)|Run mode
    skip_completed: bool = True  # 跳过已完成的步骤|Skip completed steps

    # 工具路径|Tool paths (conda环境路径)|Tool paths (conda environment paths)
    conda_env_ltr_harvest: str = "~/miniforge3/envs/ltr_harvest_parallel_v.1.2"
    conda_env_ltr_finder: str = "~/miniforge3/envs/ltr_finder_parallel_v.1.3"
    conda_env_ltr_retriever: str = "~/miniforge3/envs/ltr_retriever_v.3.0.1"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if self.genome:
            self.genome_path = Path(self.genome).resolve()  # 转换为绝对路径|Convert to absolute path
        if self.output_dir:
            self.output_path = Path(self.output_dir).resolve()  # 转换为绝对路径|Convert to absolute path
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            # 创建临时目录|Create temp directory
            self.tmp_dir = self.output_path / "tmp"
            self.tmp_dir.mkdir(parents=True, exist_ok=True)

        # 展开conda环境路径|Expand conda environment paths
        self.conda_env_ltr_harvest = expand_path(self.conda_env_ltr_harvest)
        self.conda_env_ltr_finder = expand_path(self.conda_env_ltr_finder)
        self.conda_env_ltr_retriever = expand_path(self.conda_env_ltr_retriever)

        # 设置工具路径|Set tool paths from conda environments
        harvest_conda = Path(self.conda_env_ltr_harvest)
        finder_conda = Path(self.conda_env_ltr_finder)
        retriever_conda = Path(self.conda_env_ltr_retriever)

        self.gt_path = str(harvest_conda / "bin" / "gt")
        self.ltr_finder_path = str(finder_conda / "bin" / "LTR_FINDER_parallel")
        self.ltr_retriever_path = str(retriever_conda / "bin" / "LTR_retriever")
        self.lai_path = str(retriever_conda / "bin" / "LAI")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not self.genome:
            errors.append("基因组文件不能为空|Genome file cannot be empty")
        elif not self.genome_path.exists():
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        # 检查运行模式|Check run mode
        valid_modes = ['full', 'harvest', 'retrieve', 'calculate']
        if self.mode not in valid_modes:
            errors.append(f"运行模式必须是以下之一|Run mode must be one of: {valid_modes}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查conda环境|Check conda environments
        if not Path(self.conda_env_ltr_harvest).exists():
            errors.append(f"LTR_harvest conda环境不存在|LTR_harvest conda environment does not exist: {self.conda_env_ltr_harvest}")
        if not Path(self.conda_env_ltr_finder).exists():
            errors.append(f"LTR_finder conda环境不存在|LTR_finder conda environment does not exist: {self.conda_env_ltr_finder}")
        if not Path(self.conda_env_ltr_retriever).exists():
            errors.append(f"LTR_retriever conda环境不存在|LTR_retriever conda environment does not exist: {self.conda_env_ltr_retriever}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class LTRHarvestConfig:
    """LTRharvest配置类|LTRharvest Configuration Class"""

    # LTRharvest参数|LTRharvest parameters
    minlenltr: int = 100  # 最小LTR长度|Minimum LTR length
    maxlenltr: int = 7000  # 最大LTR长度|Maximum LTR length
    mintsd: int = 4  # 最小目标位点 duplication|Minimum target site duplication
    maxtsd: int = 6  # 最大目标位点 duplication|Maximum target site duplication
    motif: str = "TGCA"  # LTR基序|LTR motif
    motifmis: int = 1  # 基序错配数|Motif mismatches
    similar: int = 85  # 相似度阈值|Similarity threshold
    vic: int = 10  # 相邻间隔距离|Vicinity distance
    seed: int = 20  # 种子长度|Seed length

    # 索引参数|Index parameters
    index_memory: int = 3500  # 索引内存限制(MB)|Index memory limit (MB)


@dataclass
class LTRFinderConfig:
    """LTR_FINDER配置类|LTR_FINDER Configuration Class"""

    # LTR_FINDER参数|LTR_FINDER parameters
    size: int = 1000000  # 分块大小|Chunk size
    time: int = 300  # 时间限制(秒)|Time limit (seconds)

    # 其他参数|Other parameters
    harvest_out: bool = True  # 使用harvest输出格式|Use harvest output format


@dataclass
class LTRRetrieverConfig:
    """LTR_retriever配置类|LTR_retriever Configuration Class"""

    # LTR_retriever参数|LTR_retriever parameters
    minlenltr: int = 100  # 最小LTR长度|Minimum LTR length
    maxlenltr: int = 7000  # 最大LTR长度|Maximum LTR length
    mintsd: int = 4  # 最小目标位点 duplication|Minimum target site duplication
    maxtsd: int = 6  # 最大目标位点 duplication|Maximum target site duplication

    # 筛选参数|Filter parameters
    seqid: bool = True  # 是否包含序列ID|Whether to include sequence IDs


@dataclass
class LAICalculateConfig:
    """LAI计算配置类|LAI Calculation Configuration Class"""

    # LAI计算参数|LAI calculation parameters
    window: int = 3000000  # 窗口大小(bp)|Window size (bp)
    step: int = 300000  # 步长(bp)|Step size (bp)

    # 估算模式|Estimation mode
    quick: bool = False  # 快速估算|Quick estimation
    qquick: bool = False  # 超快速(不估算一致性)|Super quick (no identity estimation)

    # 多倍体参数|Polyploid parameters
    mono: Optional[str] = None  # 单倍体序列列表文件|Monoploid sequence list file
    iden: Optional[float] = None  # 手动指定LTR一致性|Manually specify LTR identity
    totltr: Optional[float] = None  # 手动指定总LTR含量|Manually specify total LTR content
    genome_size: Optional[int] = None  # 手动指定基因组大小|Manually specify genome size

    # 其他参数|Other parameters
    threads: int = 4  # BLAST线程数|BLAST threads
