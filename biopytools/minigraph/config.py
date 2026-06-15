"""
Minigraph配置管理模块|Minigraph Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class MinigraphBuildConfig:
    """Minigraph构建配置类|Minigraph Build Configuration Class"""

    # 必需参数|Required parameters
    ref_fasta: str  # 参考基因组FASTA
    sample_fastas: List[str]  # 样本基因组FASTA列表

    # 输出配置|Output configuration
    output_gfa: str = './pangenome.gfa'

    # 构建参数|Build parameters
    preset: str = 'ggs'  # 图构建预设 (g/gs/ggs)
    min_identity: float = 0.9  # 最小序列相似度
    min_aln_len: int = 100000  # 最小比对长度
    max_gap: int = 1000000  # 最大gap大小

    # 性能参数|Performance parameters
    threads: int = 16
    batch_size: Optional[int] = None  # K参数 (MB)

    # 外部工具路径|External tool paths
    minigraph_path: str = 'minigraph'
    gfatools_path: str = 'gfatools'

    # 处理选项|Processing options
    keep_intermediate: bool = False  # 保留中间文件
    append_mode: bool = False  # 追加模式

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.ref_fasta = expand_path(self.ref_fasta)
        self.sample_fastas = [expand_path(f) for f in self.sample_fastas]
        self.output_gfa = expand_path(self.output_gfa)

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_gfa)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 验证preset|Validate preset
        valid_presets = ['g', 'gs', 'ggs']
        if self.preset not in valid_presets:
            raise ValueError(f"无效的preset|Invalid preset: {self.preset}. "
                           f"必须为|Must be one of {valid_presets}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        required_files = [
            ('参考基因组FASTA文件|Reference FASTA file', self.ref_fasta),
        ]

        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查样本文件|Check sample files
        for i, fasta in enumerate(self.sample_fastas, 1):
            if not os.path.exists(fasta):
                errors.append(f"样本{i}的FASTA文件不存在|Sample {i} FASTA does not exist: {fasta}")

        # 验证参数范围|Validate parameter ranges
        if not 0 < self.min_identity <= 1:
            errors.append(f"min_identity必须在(0,1]范围内|min_identity must be in (0,1]")

        if self.min_aln_len < 0:
            errors.append(f"min_aln_len必须>=0|min_aln_len must be >= 0")

        if self.threads < 1:
            errors.append(f"线程数必须>=1|Thread count must be >= 1")

        if self.batch_size is not None and self.batch_size <= 0:
            errors.append(f"batch_size必须>0|batch_size must be > 0")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class MinigraphCallConfig:
    """Minigraph SV调用配置类|Minigraph SV Call Configuration Class"""

    # 必需参数|Required parameters
    graph_gfa: str  # 泛基因组图GFA文件
    sample_fastas: List[str]  # 样本基因组FASTA列表

    # 输出配置|Output configuration
    output_dir: str = './minigraph_call'

    # 调用参数|Call parameters
    preset: str = 'asm'  # SV调用预设 (asm)
    call_mode: bool = True  # 启用--call模式

    # 性能参数|Performance parameters
    threads: int = 16

    # 外部工具路径|External tool paths
    minigraph_path: str = 'minigraph'

    # 处理选项|Processing options
    output_bed: bool = True  # 输出BED格式
    output_gaf: bool = False  # 输出GAF格式

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.graph_gfa = expand_path(self.graph_gfa)
        self.sample_fastas = [expand_path(f) for f in self.sample_fastas]
        self.output_dir = expand_path(self.output_dir)

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 验证preset|Validate preset
        valid_presets = ['asm']
        if self.preset not in valid_presets:
            raise ValueError(f"无效的preset|Invalid preset: {self.preset}. "
                           f"必须为|Must be one of {valid_presets}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.graph_gfa):
            errors.append(f"GFA文件不存在|GFA file does not exist: {self.graph_gfa}")

        # 检查样本文件|Check sample files
        for i, fasta in enumerate(self.sample_fastas, 1):
            if not os.path.exists(fasta):
                errors.append(f"样本{i}的FASTA文件不存在|Sample {i} FASTA does not exist: {fasta}")

        # 验证参数|Validate parameters
        if self.threads < 1:
            errors.append(f"线程数必须>=1|Thread count must be >= 1")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class MinigraphBubbleConfig:
    """Minigraph bubble提取配置类|Minigraph Bubble Extraction Configuration Class"""

    # 必需参数|Required parameters
    graph_gfa: str  # 泛基因组图GFA文件

    # 输出配置|Output configuration
    output_bed: str = './sv_bubbles.bed'

    # 外部工具路径|External tool paths
    gfatools_path: str = 'gfatools'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.graph_gfa = expand_path(self.graph_gfa)
        self.output_bed = expand_path(self.output_bed)

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_bed)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.graph_gfa):
            errors.append(f"GFA文件不存在|GFA file does not exist: {self.graph_gfa}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class MinigraphMapConfig:
    """Minigraph序列映射配置类|Minigraph Sequence Mapping Configuration Class"""

    # 必需参数|Required parameters
    graph_gfa: str  # 泛基因组图GFA文件
    query_fastas: List[str]  # 查询序列FASTA列表

    # 输出配置|Output configuration
    output_gaf: str = './mapping.gaf'

    # 映射参数|Mapping parameters
    preset: str = 'lr'  # 映射预设 (sr/lr/map-pb/map-ont/asm)
    max_intron_len: Optional[int] = None  # 最大内含子长度（用于asm）

    # 性能参数|Performance parameters
    threads: int = 16
    batch_size: Optional[int] = None  # K参数 (MB)

    # 外部工具路径|External tool paths
    minigraph_path: str = 'minigraph'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.graph_gfa = expand_path(self.graph_gfa)
        self.query_fastas = [expand_path(f) for f in self.query_fastas]
        self.output_gaf = expand_path(self.output_gaf)

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_gaf)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 验证preset|Validate preset
        valid_presets = ['sr', 'lr', 'map-pb', 'map-ont', 'asm']
        if self.preset not in valid_presets:
            raise ValueError(f"无效的preset|Invalid preset: {self.preset}. "
                           f"必须为|Must be one of {valid_presets}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.graph_gfa):
            errors.append(f"GFA文件不存在|GFA file does not exist: {self.graph_gfa}")

        # 检查查询文件|Check query files
        for i, fasta in enumerate(self.query_fastas, 1):
            if not os.path.exists(fasta):
                errors.append(f"查询{i}的FASTA文件不存在|Query {i} FASTA does not exist: {fasta}")

        # 验证参数|Validate parameters
        if self.threads < 1:
            errors.append(f"线程数必须>=1|Thread count must be >= 1")

        if self.batch_size is not None and self.batch_size <= 0:
            errors.append(f"batch_size必须>0|batch_size must be > 0")

        if errors:
            raise ValueError("\n".join(errors))

        return True
