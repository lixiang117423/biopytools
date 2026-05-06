"""
GFF重命名配置管理模块|GFF Renamer Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class GFFRenamerConfig:
    """GFF重命名配置类|GFF Renamer Configuration Class"""

    # 输入输出参数|Input and output parameters
    input_file: str
    output_file: str

    # 重命名参数|Renaming parameters
    prefix: str  # ID前缀|ID prefix (e.g., CDRT)
    species: str  # 物种缩写|Species abbreviation (e.g., Ov)

    # AGAT清洗参数|AGAT cleaning parameters
    agat_path: str = field(
        default_factory=lambda: get_tool_path('agat', '~/miniforge3/envs/agat_v.1.7.0/bin/agat_convert_sp_gxf2gxf.pl', 'AGAT_PATH')
    )
    skip_gff_clean: bool = False

    # 并行处理参数|Parallel processing parameters
    threads: int = 12  # 并行线程数|Number of parallel threads

    # mRNA映射文件参数|mRNA mapping file parameters
    output_mrna_mapping: bool = False  # 是否输出mRNA映射文件|Whether to output mRNA mapping file
    mrna_mapping_file: Optional[str] = None  # mRNA映射文件路径|mRNA mapping file path

    # 新增参数|New parameters
    chr_mapping_file: Optional[str] = None  # 染色体映射文件路径|Chromosome mapping file path
    naming_format: str = "standard"  # 命名格式|Naming format (standard/simple/compact)
    include_utr: bool = False  # 是否包含UTR特征|Whether to include UTR features

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 规范化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(expand_path(self.input_file)))
        self.output_file = os.path.normpath(os.path.abspath(expand_path(self.output_file)))
        self.agat_path = os.path.normpath(os.path.abspath(expand_path(self.agat_path)))

        # 创建输出目录|Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)

        # 处理mRNA映射文件路径|Process mRNA mapping file path
        if self.output_mrna_mapping and self.mrna_mapping_file:
            self.mrna_mapping_file = os.path.normpath(os.path.abspath(self.mrna_mapping_file))
            # 创建映射文件目录|Create mapping file directory
            mapping_dir = Path(self.mrna_mapping_file).parent
            mapping_dir.mkdir(parents=True, exist_ok=True)
        elif self.output_mrna_mapping and not self.mrna_mapping_file:
            # 如果未指定路径，自动生成|Auto-generate if not specified
            base_name = Path(self.output_file).stem
            self.mrna_mapping_file = str(Path(self.output_file).parent / f"{base_name}_mrna_mapping.tsv")

        # 规范化染色体映射文件路径|Normalize chromosome mapping file path
        if self.chr_mapping_file:
            self.chr_mapping_file = os.path.normpath(os.path.abspath(self.chr_mapping_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        if not self.input_file.endswith(('.gff', '.gff3')):
            errors.append(
                f"输入文件格式不正确|Input file format incorrect: {self.input_file}. "
                f"支持|Support: .gff, .gff3"
            )

        # 检查参数|Check parameters
        if not self.prefix:
            errors.append(f"前缀不能为空|Prefix cannot be empty")

        if not self.species:
            errors.append(f"物种缩写不能为空|Species abbreviation cannot be empty")

        # 检查染色体映射文件|Check chromosome mapping file
        if self.chr_mapping_file and not os.path.exists(self.chr_mapping_file):
            errors.append(f"染色体映射文件不存在|Chromosome mapping file does not exist: {self.chr_mapping_file}")

        # 检查命名格式|Check naming format
        valid_formats = ['standard', 'simple', 'compact']
        if self.naming_format not in valid_formats:
            errors.append(f"命名格式无效|Invalid naming format: {self.naming_format}. 支持格式|Supported formats: {', '.join(valid_formats)}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
