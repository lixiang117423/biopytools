"""基因密度计算配置类|Gene density configuration"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path


@dataclass
class GeneDensityConfig:
    """基因密度计算配置类|Gene density calculation configuration

    按固定大小窗口统计每条染色体各区间的基因数量与基因密度(基因/Mb)
    |Count genes and gene density (genes/Mb) per fixed-size window along each chromosome
    """

    # 必需参数|Required
    gff_file: str

    # 可选参数|Optional
    output_dir: str = "./gene_density_output"
    window_size: int = 100000
    feature_type: str = "gene"
    genome_file: Optional[str] = None
    prefix: Optional[str] = None
    generate_plot: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing

        - 输入路径展开为绝对路径(§11.3.1)|Expand input paths to absolute
        - 创建输出目录|Create output directory
        - 设置by-step输出路径(§12)|Set by-step output paths
        """
        # 输入文件展开~/$VAR并解析为绝对路径(下游计算依赖绝对路径, §11.3.1)
        # |Expand ~/$VAR and resolve input files to absolute (§11.3.1)
        self.gff_file = str(Path(expand_path(self.gff_file)).resolve())
        if self.genome_file:
            self.genome_file = str(Path(expand_path(self.genome_file)).resolve())

        # 输出目录|Output directory
        self.output_dir = expand_path(self.output_dir)
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 前缀默认值(可追溯性§12.3)|Prefix default (traceability)
        if not self.prefix:
            self.prefix = Path(self.gff_file).stem

        # by-step输出路径|By-step output paths (§12)
        self.density_tsv = str(
            Path(self.output_dir) / "01_density" / f"{self.prefix}.gene_density.tsv"
        )
        self.plot_png = str(
            Path(self.output_dir) / "02_plot" / f"{self.prefix}.gene_density.png"
        )
        self.versions_yml = str(
            Path(self.output_dir) / "00_pipeline_info" / "software_versions.yml"
        )

    def validate(self):
        """验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当任何校验失败|When any validation fails
        """
        errors = []

        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file not found: {self.gff_file}")

        if self.window_size <= 0:
            errors.append(f"窗口大小必须为正数|Window size must be positive: {self.window_size}")

        if self.genome_file and not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
