"""
WGDI配置管理模块|WGDI Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import get_tool_path, expand_path


@dataclass
class WGDIConfig:
    """WGDI通用配置类|WGDI Common Configuration Class"""

    # WGDI软件路径|WGDI software path
    wgdi_path: str = field(
        default_factory=lambda: get_tool_path(
            'wgdi',
            '~/miniforge3/envs/wgdi_v.0.75/bin/wgdi',
            'WGDI_PATH'
        )
    )

    # 通用参数|Common parameters
    threads: int = 8
    working_dir: str = "."

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.wgdi_path = expand_path(self.wgdi_path)
        self.working_dir = expand_path(self.working_dir)

        # 创建工作目录|Create working directory
        Path(self.working_dir).mkdir(parents=True, exist_ok=True)


@dataclass
class DotPlotConfig(WGDIConfig):
    """DotPlot专用配置|DotPlot-specific configuration"""

    # 必需文件|Required files
    blast_file: str = ""
    gff1_file: str = ""
    gff2_file: str = ""
    lens1_file: str = ""
    lens2_file: str = ""

    # 可选参数|Optional parameters
    genome1_name: str = "Genome1"
    genome2_name: str = "Genome2"
    blast_reverse: bool = False
    multiple: int = 1
    score: int = 100
    evalue: float = 1e-5
    repeat_number: int = 10
    position: str = "order"  # order | start | end

    # 可选文件|Optional files
    ancestor_left: Optional[str] = None
    ancestor_top: Optional[str] = None

    # 输出参数|Output parameters
    markersize: float = 0.5
    figsize: str = "10,10"
    savefig: str = "dotplot.png"

    def __post_init__(self):
        super().__post_init__()
        # 展开所有文件路径|Expand all file paths
        self.blast_file = expand_path(self.blast_file)
        self.gff1_file = expand_path(self.gff1_file)
        self.gff2_file = expand_path(self.gff2_file)
        self.lens1_file = expand_path(self.lens1_file)
        self.lens2_file = expand_path(self.lens2_file)

        if self.ancestor_left:
            self.ancestor_left = expand_path(self.ancestor_left)
        if self.ancestor_top:
            self.ancestor_top = expand_path(self.ancestor_top)

    def validate(self):
        """验证配置参数|Validate configuration"""
        errors = []

        required_files = [
            ("BLAST文件|BLAST file", self.blast_file),
            ("GFF1文件|GFF1 file", self.gff1_file),
            ("GFF2文件|GFF2 file", self.gff2_file),
            ("LENS1文件|LENS1 file", self.lens1_file),
            ("LENS2文件|LENS2 file", self.lens2_file),
        ]

        for desc, path in required_files:
            if not Path(path).exists():
                errors.append(f"{desc}不存在|does not exist: {path}")

        if self.position not in ["order", "start", "end"]:
            errors.append(f"position参数错误|position must be order/start/end: {self.position}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class CollinearityConfig(WGDIConfig):
    """Collinearity专用配置|Collinearity-specific configuration"""

    # 必需文件|Required files
    blast_file: str = ""
    gff1_file: str = ""
    gff2_file: str = ""
    lens1_file: str = ""
    lens2_file: str = ""

    # 处理参数|Processing parameters
    blast_reverse: bool = False
    comparison: str = "genomes"  # genomes | chromosomes
    multiple: int = 1
    process: int = 8
    evalue: float = 1e-5
    score: int = 100

    # 评分参数|Scoring parameters
    grading: str = "50,40,25"  # red,blue,gray scores
    mg: str = "40,40"  # max_gap values

    # 其他参数|Other parameters
    pvalue: float = 1.0
    repeat_number: int = 20
    position: str = "order"

    # 输出参数|Output parameters
    savefile: str = "collinearity.txt"

    def __post_init__(self):
        super().__post_init__()
        # 展开路径|Expand paths
        self.blast_file = expand_path(self.blast_file)
        self.gff1_file = expand_path(self.gff1_file)
        self.gff2_file = expand_path(self.gff2_file)
        self.lens1_file = expand_path(self.lens1_file)
        self.lens2_file = expand_path(self.lens2_file)

    def validate(self):
        """验证配置参数|Validate configuration"""
        errors = []

        required_files = [
            ("BLAST文件|BLAST file", self.blast_file),
            ("GFF1文件|GFF1 file", self.gff1_file),
            ("GFF2文件|GFF2 file", self.gff2_file),
            ("LENS1文件|LENS1 file", self.lens1_file),
            ("LENS2文件|LENS2 file", self.lens2_file),
        ]

        for desc, path in required_files:
            if not Path(path).exists():
                errors.append(f"{desc}不存在|does not exist: {path}")

        if self.comparison not in ["genomes", "chromosomes"]:
            errors.append(f"comparison参数错误|comparison must be genomes/chromosomes")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class CalKsConfig(WGDIConfig):
    """CalKs专用配置|CalKs-specific configuration"""

    # 必需文件|Required files
    collinearity_file: str = ""
    fasta1_file: str = ""
    fasta2_file: str = ""

    # 输出参数|Output parameters
    savefile: str = "ks_results.txt"

    def __post_init__(self):
        super().__post_init__()
        # 展开路径|Expand paths
        self.collinearity_file = expand_path(self.collinearity_file)
        self.fasta1_file = expand_path(self.fasta1_file)
        self.fasta2_file = expand_path(self.fasta2_file)

    def validate(self):
        """验证配置参数|Validate configuration"""
        errors = []

        required_files = [
            ("共线性文件|Collinearity file", self.collinearity_file),
            ("FASTA1文件|FASTA1 file", self.fasta1_file),
            ("FASTA2文件|FASTA2 file", self.fasta2_file),
        ]

        for desc, path in required_files:
            if not Path(path).exists():
                errors.append(f"{desc}不存在|does not exist: {path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
