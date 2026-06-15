"""
蛋白质到基因组比对配置管理模块|Protein to Genome Alignment Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def expand_path(file_path: str) -> str:
    """展开路径中的~符号|Expand tilde in path

    Args:
        file_path: 文件路径|File path

    Returns:
        str: 展开后的绝对路径|Expanded absolute path
    """
    return os.path.normpath(os.path.abspath(os.path.expanduser(file_path)))


@dataclass
class Pep2GenomeConfig:
    """蛋白质到基因组比对配置类|Protein to Genome Alignment Configuration Class"""

    # 输入文件|Input files
    genome_fa: str  # 基因组FASTA文件|Genome FASTA file
    protein_fa: str  # 蛋白质FASTA文件|Protein FASTA file

    # 输出目录|Output directory
    output_dir: str  # 输出目录|Output directory

    # 处理参数|Processing parameters
    threads: int = 12  # 线程数|Number of threads

    # Miniprot路径配置|Miniprot path configuration
    miniprot_path: str = "miniprot"  # 默认使用工具名，让系统从PATH查找|Default uses tool name, let system find from PATH

    # 输出选项|Output options
    export_gff3: bool = True  # 是否导出GFF3格式|Whether to export GFF3 format
    export_bed: bool = True  # 是否导出BED格式|Whether to export BED format
    export_statistics: bool = True  # 是否生成统计报告|Whether to generate statistics report
    extract_sequences: bool = True  # 是否提取基因组序列|Whether to extract genome sequences

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.genome_fa = expand_path(self.genome_fa)
        self.protein_fa = expand_path(self.protein_fa)
        self.output_dir = expand_path(self.output_dir)

        # 只有在miniprot_path不是默认值时才展开|Only expand miniprot_path if not default
        if self.miniprot_path and self.miniprot_path != "miniprot":
            self.miniprot_path = expand_path(self.miniprot_path)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome_fa):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_fa}")

        # 检查蛋白质文件|Check protein file
        if not os.path.exists(self.protein_fa):
            errors.append(f"蛋白质文件不存在|Protein file does not exist: {self.protein_fa}")

        # 检查miniprot工具（仅当用户提供自定义路径时）|Check miniprot tool (only when user provides custom path)
        if self.miniprot_path != "miniprot" and not os.path.exists(self.miniprot_path):
            errors.append(f"miniprot工具不存在|miniprot tool does not exist: {self.miniprot_path}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
