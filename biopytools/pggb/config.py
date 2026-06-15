"""
PGGB配置类|PGGB Configuration Class
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import expand_path


@dataclass
class PGGBConfig:
    """PGGB泛基因组图构建配置类|PGGB Pangenome Graph Builder Configuration Class"""

    # 必需参数|Required parameters
    input_fa: str = ""
    output_dir: str = ""

    # conda环境|conda environment
    conda_env: str = "pggb_v.0.7.4"

    # wfmash参数|wfmash parameters
    segment_length: int = 5000       # 比对分段长度|Segment length for mapping
    block_length: int = 0            # 最小block长度过滤(0=自动=5*segment_length)|Min block length (0=auto=5*segment_length)
    map_pct_id: int = 90             # 比对一致度|Percent identity for mapping
    n_mappings: int = 1              # 每个segment的mapping数量|Number of mappings per segment
    mash_kmer: int = 19              # mash kmer大小|Mash kmer size
    mash_kmer_thres: float = 0.001   # 忽略最频繁kmer的比例|Ignore top frequent kmers fraction
    no_splits: bool = False           # 禁用序列拆分|Disable sequence splitting during mapping
    sparse_map: str = ""             # 稀疏映射比例(''=不稀疏, 'auto'=自动)|Sparse mapping fraction
    exclude_delim: str = ""          # 跳过相同前缀的映射分隔符|Delimiter for excluding self-mappings

    # seqwish参数|seqwish parameters
    min_match_length: int = 23       # 最小匹配长度|Minimum match length
    sparse_factor: float = 0.0       # 稀疏因子|Sparse factor (0=no sparsification)
    transclose_batch: str = "10M"    # 传递闭包批处理大小|Transitive closure batch size

    # smoothxg参数|smoothxg parameters
    skip_normalization: bool = False # 跳过图归一化|Skip graph normalization
    n_haplotypes: int = 0            # 单倍型数(0=自动检测)|Number of haplotypes (0=auto-detect)
    max_path_jump: int = 0           # 最大路径跳跃|Maximum path jump
    max_edge_jump: int = 0           # 最大边跳跃|Maximum edge jump
    target_poa_length: str = "700,1100"  # POA目标长度|POA target sequence length
    poa_params: str = ""             # POA参数(''=默认asm20)|POA parameters (''=default asm20)
    poa_padding: float = 0.001       # POA填充|POA padding
    pad_max_depth: int = 100         # 最大填充深度|Maximum padding depth

    # odgi/vg参数|odgi/vg parameters
    skip_viz: bool = True            # 跳过可视化(默认跳过，节省时间)|Skip visualizations
    do_stats: bool = False           # 生成统计信息|Generate statistics
    vcf_spec: str = ""               # VCF输出参考规范|VCF output reference spec

    # 通用参数|General parameters
    threads: int = 24                # 线程数|Number of compute threads
    poa_threads: int = 0             # POA线程数(0=同threads)|POA threads (0=same as threads)
    resume: bool = False             # 断点续传|Resume from existing outputs
    compress: bool = False           # 压缩输出|Compress output files
    keep_temp_files: bool = False    # 保留临时文件|Keep intermediate files
    input_paf: str = ""              # 外部PAF文件(跳过wfmash)|External PAF file (skip wfmash)
    temp_dir: str = ""               # 临时文件目录(默认=output_dir)|Temporary directory

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        if self.input_fa:
            self.input_fa = expand_path(self.input_fa)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
        if self.input_paf:
            self.input_paf = expand_path(self.input_paf)
        if self.temp_dir:
            self.temp_dir = expand_path(self.temp_dir)

        # 确保输出目录存在|Ensure output directory exists
        if self.output_dir:
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not self.input_fa:
            errors.append("输入FASTA文件不能为空|Input FASTA file cannot be empty")
        elif not Path(self.input_fa).exists():
            errors.append(f"输入FASTA文件不存在|Input FASTA file not found: {self.input_fa}")

        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        if self.threads <= 0:
            errors.append("线程数必须为正整数|Thread count must be positive")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"PGGBConfig(\n"
            f"  input_fa={self.input_fa!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  conda_env={self.conda_env!r},\n"
            f"  segment_length={self.segment_length!r},\n"
            f"  map_pct_id={self.map_pct_id!r},\n"
            f"  n_haplotypes={self.n_haplotypes!r},\n"
            f"  threads={self.threads!r},\n"
            f"  resume={self.resume!r}\n"
            f")"
        )
