"""
序列到基因组比对配置管理模块|Sequence to Genome Alignment Configuration Management Module
支持DNA和蛋白质序列|Support DNA and protein sequences
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
class Seq2GenomeConfig:
    """序列到基因组比对配置类|Sequence to Genome Alignment Configuration Class"""

    # 必需参数（无默认值）|Required parameters (no default values)
    genome_fa: str  # 基因组FASTA文件|Genome FASTA file
    query_fa: str  # 查询序列FASTA文件（DNA或蛋白质）|Query sequence FASTA file (DNA or protein)
    output_dir: str  # 输出目录|Output directory

    # 可选参数（有默认值）|Optional parameters (with default values)
    # 向后兼容：蛋白质FASTA文件|Backward compatibility: Protein FASTA file
    protein_fa: Optional[str] = None
    # 序列类型：'dna', 'protein', 或None（自动检测）|Sequence type: 'dna', 'protein', or None (auto-detect)
    query_type: Optional[str] = None

    # 处理参数|Processing parameters
    threads: int = 12  # 线程数|Number of threads

    # 工具路径配置|Tool path configuration
    miniprot_path: str = "miniprot"  # Miniprot工具路径|Miniprot tool path
    minimap2_path: str = "minimap2"  # Minimap2工具路径|Minimap2 tool path

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
        self.query_fa = expand_path(self.query_fa)
        self.output_dir = expand_path(self.output_dir)

        # 向后兼容：如果使用了旧的protein_fa参数名，自动转换为query_fa
        # Backward compatibility: if old protein_fa parameter is used, automatically convert to query_fa
        if hasattr(self, 'protein_fa'):
            # 如果直接传入了protein_fa参数|If protein_fa parameter is passed directly
            if self.protein_fa and self.protein_fa != self.query_fa:
                self.query_fa = expand_path(self.protein_fa)

        # 只有在工具路径不是默认值时才展开|Only expand tool paths if not default
        if self.miniprot_path and self.miniprot_path != "miniprot":
            self.miniprot_path = expand_path(self.miniprot_path)
        if self.minimap2_path and self.minimap2_path != "minimap2":
            self.minimap2_path = expand_path(self.minimap2_path)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome_fa):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_fa}")

        # 检查查询序列文件|Check query sequence file
        if not os.path.exists(self.query_fa):
            errors.append(f"查询序列文件不存在|Query sequence file does not exist: {self.query_fa}")

        # 检查miniprot工具（仅当用户提供自定义路径时）|Check miniprot tool (only when user provides custom path)
        if self.miniprot_path != "miniprot" and not os.path.exists(self.miniprot_path):
            errors.append(f"miniprot工具不存在|miniprot tool does not exist: {self.miniprot_path}")

        # 检查minimap2工具（仅当用户提供自定义路径时）|Check minimap2 tool (only when user provides custom path)
        if self.minimap2_path != "minimap2" and not os.path.exists(self.minimap2_path):
            errors.append(f"minimap2工具不存在|minimap2 tool does not exist: {self.minimap2_path}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        # 验证序列类型参数（如果指定了）|Validate sequence type parameter (if specified)
        if self.query_type is not None:
            if self.query_type not in ['dna', 'protein', 'auto']:
                errors.append(f"序列类型必须是'dna', 'protein'或'auto'|Sequence type must be 'dna', 'protein' or 'auto': {self.query_type}")
            if self.query_type == 'auto':
                self.query_type = None  # None表示自动检测|None means auto-detect

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_aligner_type(self) -> str:
        """
        获取比对器类型|Get aligner type

        Returns:
            str: 'minimap2' 或 'miniprot'
        """
        if self.query_type:
            return 'minimap2' if self.query_type == 'dna' else 'miniprot'
        else:
            # 如果未指定类型，将在运行时自动检测|If type not specified, will auto-detect at runtime
            return 'auto'


# 保持向后兼容的别名|Keep backward compatible alias
Pep2GenomeConfig = Seq2GenomeConfig
