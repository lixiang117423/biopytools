"""
PopLDdecay配置管理模块|PopLDdecay Configuration Management Module
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path


@dataclass
class PopLDdecayConfig:
    """PopLDdecay配置类|PopLDdecay Configuration Class"""

    # 必需参数|Required parameters
    input_file: str = ""  # 输入文件（VCF或Genotype格式）|Input file (VCF or Genotype format)
    output_prefix: str = ""  # 输出文件前缀|Output file prefix

    # 输入类型|Input type
    input_type: str = "vcf"  # 输入文件类型|Input file type (vcf/genotype)

    # 过滤参数|Filter parameters
    max_dist: int = 300  # 最大距离（kb）|Max distance in kb
    min_maf: float = 0.005  # 最小等位基因频率|Min minor allele frequency
    max_het: float = 0.88  # 最大杂合率|Max ratio of heterozygous alleles
    max_miss: float = 0.25  # 最大缺失率|Max ratio of missing alleles

    # 子群体分析|Subgroup analysis
    subpop_file: str = ""  # 子群体样本列表文件|Subgroup sample list file

    # 输出类型|Output type
    out_type: int = 1  # 输出类型|Output type (1: R^2, 2: R^2 & D', 3: Pairwise LD)

    # 其他参数|Other parameters
    ehh_site: str = "NA"  # EHH区域起始位点|EHH region start site
    output_filtered_snp: bool = False  # 输出过滤后的SNP|Output filtered SNPs

    # 绘图参数|Plotting parameters
    plot: bool = True  # 是否绘图|Whether to plot
    bin1: int = 10  # 短距离bin大小|Bin size for short distance
    bin2: int = 100  # 长距离bin大小|Bin size for long distance
    break_point: int = 100  # 短/长距离分界点|Break point between short/long distance
    max_x: Optional[int] = None  # 最大X坐标|Max X coordinate
    measure: str = "r2"  # LD度量方法|LD measure method (r2/D/both)
    method: str = "MeanBin"  # 绘图方法|Plotting method (MeanBin/HW/MedianBin/PercentileBin)
    percentile: float = 0.5  # 百分位数（用于PercentileBin）|Percentile for PercentileBin

    # 阈值推荐参数|Threshold recommendation parameters
    recommend_threshold: bool = True  # 是否推荐LD阈值|Whether to recommend LD threshold (default: True)

    # 软件路径|Software paths
    poplddecay_path: str = field(
        default_factory=lambda: get_tool_path(
            'poplddecay',
            '~/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay',
            'POPLDDECAY_PATH'
        )
    )
    plot_one_pop_path: str = field(
        default_factory=lambda: expand_path(
            '~/software/PopLDdecay/PopLDdecay-3.43/bin/Plot_OnePop.pl'
        )
    )
    plot_multi_pop_path: str = field(
        default_factory=lambda: expand_path(
            '~/software/PopLDdecay/PopLDdecay-3.43/bin/Plot_MultiPop.pl'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.poplddecay_path = expand_path(self.poplddecay_path)
        self.plot_one_pop_path = expand_path(self.plot_one_pop_path)
        self.plot_multi_pop_path = expand_path(self.plot_multi_pop_path)
        if self.input_file:
            self.input_path = Path(self.input_file)
        if self.output_prefix:
            self.output_path = Path(self.output_prefix)
            # 创建输出目录|Create output directory
            self.output_path.parent.mkdir(parents=True, exist_ok=True)
        if self.subpop_file:
            self.subpop_path = Path(self.subpop_file)

    def validate(self):
        """验证配置参数|Validate configuration parameters

        Returns:
            bool: 验证通过返回True|Return True if validation passes

        Raises:
            ValueError: 参数验证失败时|When parameter validation fails
        """
        errors = []

        # 检查输入文件|Check input file
        if not self.input_file:
            errors.append("输入文件不能为空|Input file cannot be empty")
        elif not self.input_path.exists():
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        # 检查输出文件|Check output file
        if not self.output_prefix:
            errors.append("输出文件前缀不能为空|Output file prefix cannot be empty")

        # 检查输入类型|Check input type
        if self.input_type not in ["vcf", "genotype"]:
            errors.append(f"输入类型必须是vcf或genotype|Input type must be vcf or genotype: {self.input_type}")

        # 检查max_dist|Check max_dist
        if self.max_dist <= 0:
            errors.append(f"最大距离必须为正数|Max distance must be positive: {self.max_dist}")

        # 检查MAF范围|Check MAF range
        if not 0.0 <= self.min_maf <= 1.0:
            errors.append(f"最小等位基因频率必须在0-1之间|Min MAF must be between 0-1: {self.min_maf}")

        # 检查杂合率范围|Check heterozygous rate range
        if not 0.0 <= self.max_het <= 1.0:
            errors.append(f"最大杂合率必须在0-1之间|Max het rate must be between 0-1: {self.max_het}")

        # 检查缺失率范围|Check missing rate range
        if not 0.0 <= self.max_miss <= 1.0:
            errors.append(f"最大缺失率必须在0-1之间|Max missing rate must be between 0-1: {self.max_miss}")

        # 检查子群体文件|Check subpopulation file
        if self.subpop_file and not self.subpop_path.exists():
            errors.append(f"子群体文件不存在|Subpopulation file does not exist: {self.subpop_file}")

        # 检查输出类型|Check output type
        if not 1 <= self.out_type <= 8:
            errors.append(f"输出类型必须在1-8之间|Output type must be between 1-8: {self.out_type}")

        # 检查绘图参数|Check plotting parameters
        if self.bin1 <= 0:
            errors.append(f"短距离bin大小必须为正数|bin1 must be positive: {self.bin1}")
        if self.bin2 <= 0:
            errors.append(f"长距离bin大小必须为正数|bin2 must be positive: {self.bin2}")
        if self.break_point <= 0:
            errors.append(f"分界点必须为正数|break point must be positive: {self.break_point}")
        if self.measure not in ["r2", "D", "both"]:
            errors.append(f"度量方法必须是r2、D或both|Measure must be r2, D, or both: {self.measure}")
        if self.method not in ["MeanBin", "HW", "MedianBin", "PercentileBin"]:
            errors.append(f"绘图方法无效|Invalid plotting method: {self.method}")
        if not 0.0 <= self.percentile <= 1.0:
            errors.append(f"百分位数必须在0-1之间|Percentile must be between 0-1: {self.percentile}")

        # 检查软件路径|Check software paths
        plot_one_pop = Path(self.plot_one_pop_path)
        if not plot_one_pop.exists():
            errors.append(f"Plot_OnePop.pl脚本不存在|Plot_OnePop.pl script does not exist: {self.plot_one_pop_path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_output_files(self) -> dict:
        """获取所有输出文件路径|Get all output file paths

        Returns:
            dict: 输出文件路径字典|Dictionary of output file paths
        """
        output_prefix = str(self.output_path)
        return {
            'stat': f"{output_prefix}.stat.gz",
            'figure': f"{output_prefix}.png",
            'log': f"{output_prefix}.log",
            'threshold': f"{output_prefix}_threshold_recommendations.tsv"
        }
