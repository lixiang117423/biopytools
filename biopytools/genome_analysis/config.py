"""
Genome Analysis Configuration
基因组分析配置类
"""

import os
from dataclasses import dataclass, field
from typing import Optional

from ..common.paths import get_domain_tool_path, expand_path


@dataclass
class GenomeAnalysisConfig:
    """基因组分析配置类|Genome Analysis Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str
    output_dir: str

    # 可选参数|Optional parameters
    read_length: int = 150
    kmer_size: int = 21
    threads: int = 64
    hash_size: str = "10G"
    max_kmer_cov: int = 1000
    read1_suffix: str = "*_1.clean.fq.gz"
    skip_smudgeplot: bool = False
    ploidy: int = 2

    # 工具路径|Tool paths
    # 域环境映射见 biopytools/common/env_map.py (jellyfish/genomescope2 归 asm 域;
    # FastK/smudgeplot 不在映射表, 保持裸命令名作为 legacy 回退)
    # |Domain mapping in common/env_map.py (jellyfish/genomescope2 -> asm; FastK/
    # smudgeplot not mapped, keep bare command name as legacy fallback)
    jellyfish_path: str = field(
        default_factory=lambda: get_domain_tool_path('jellyfish', 'jellyfish', 'JELLYFISH_PATH')
    )
    genomescope2_path: str = field(
        default_factory=lambda: get_domain_tool_path('genomescope2', 'genomescope2', 'GENOMESCOPE2_PATH')
    )
    fastk_path: str = field(
        default_factory=lambda: get_domain_tool_path('FastK', 'FastK', 'FASTK_PATH')
    )
    smudgeplot_path: str = field(
        default_factory=lambda: get_domain_tool_path('smudgeplot', 'smudgeplot', 'SMUDGEPLOT_PATH')
    )

    # 运行时变量|Runtime variables
    kcov: Optional[float] = field(default=None, init=False)
    inferred_ploidy: Optional[int] = field(default=None, init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有工具路径(含~的域环境路径展开为绝对路径; 裸命令名保持原样)
        # |Expand tool paths (~ -> absolute; bare command names stay as-is)
        self.jellyfish_path = expand_path(self.jellyfish_path)
        self.genomescope2_path = expand_path(self.genomescope2_path)
        self.fastk_path = expand_path(self.fastk_path)
        self.smudgeplot_path = expand_path(self.smudgeplot_path)

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时|When configuration parameters are invalid
        """
        # 检查必需目录|Check required directories
        if not self.input_dir or not os.path.exists(self.input_dir):
            raise ValueError(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 检查参数范围|Check parameter ranges
        if self.read_length <= 0:
            raise ValueError(f"读长必须大于0|Read length must be > 0: {self.read_length}")

        if self.kmer_size <= 0:
            raise ValueError(f"K-mer大小必须大于0|K-mer size must be > 0: {self.kmer_size}")

        if self.threads <= 0:
            raise ValueError(f"线程数必须大于0|Threads must be > 0: {self.threads}")

        if self.max_kmer_cov <= 0:
            raise ValueError(f"最大覆盖度必须大于0|Max coverage must be > 0: {self.max_kmer_cov}")

        if self.ploidy is not None and (self.ploidy < 1 or self.ploidy > 6):
            raise ValueError(f"倍性必须在1-6之间|Ploidy must be between 1-6: {self.ploidy}")

        return True

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"GenomeAnalysisConfig(\n"
            f"  input_dir={self.input_dir!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  read_length={self.read_length!r},\n"
            f"  kmer_size={self.kmer_size!r},\n"
            f"  threads={self.threads!r},\n"
            f"  hash_size={self.hash_size!r},\n"
            f"  max_kmer_cov={self.max_kmer_cov!r},\n"
            f"  read1_suffix={self.read1_suffix!r},\n"
            f"  skip_smudgeplot={self.skip_smudgeplot!r},\n"
            f"  ploidy={self.ploidy!r}\n"
            f")"
        )
