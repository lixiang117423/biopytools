"""
NCBI分类学注释配置管理模块|NCBI Taxonomy Annotation Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class NCBITaxoConfig:
    """NCBI分类学注释配置类|NCBI Taxonomy Annotation Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_prefix: str

    # 输入类型配置|Input type configuration
    input_type: str = "auto"          # auto/blast/accession
    blast_column: int = 2              # BLAST结果中accession所在的列（从1开始）
    min_alignment_length: int = 1000   # 最小比对长度过滤（bp）
    fetch_titles: bool = False         # 是否获取accession的序列描述

    # 数据库配置|Database configuration
    taxid_db: str = "~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz"

    # 分类学配置|Taxonomy configuration
    lineage_format: str = "{k};{p};{c};{o};{f};{g};{s}"
    keep_full_lineage: bool = True     # 是否保留完整lineage

    # 统计配置|Statistics configuration
    stats_by: list = None              # 统计层级，默认['genus', 'species']
    stats_target: str = "both"         # blast_hits/unique_accessions/both
    stats_output: str = "txt"          # txt/csv

    # 性能配置|Performance configuration
    threads: int = 4

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认统计层级|Set default statistics levels
        if self.stats_by is None:
            self.stats_by = ['genus', 'species']

        # 展开数据库路径|Expand database path
        self.taxid_db = os.path.expanduser(self.taxid_db)
        if not os.path.isabs(self.taxid_db):
            self.taxid_db = os.path.abspath(self.taxid_db)

        # 标准化输入输出路径|Normalize input and output paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(os.path.abspath(self.output_prefix))
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        # 检查数据库文件|Check database file
        if not os.path.exists(self.taxid_db):
            errors.append(f"TaxID数据库不存在|TaxID database not found: {self.taxid_db}")

        # 验证输入类型|Validate input type
        if self.input_type not in ['auto', 'blast', 'accession']:
            errors.append(f"无效的输入类型|Invalid input type: {self.input_type}")

        # 验证统计目标|Validate statistics target
        if self.stats_target not in ['blast_hits', 'unique_accessions', 'both']:
            errors.append(f"无效的统计目标|Invalid statistics target: {self.stats_target}")

        # 验证统计输出格式|Validate statistics output format
        if self.stats_output not in ['txt', 'csv']:
            errors.append(f"无效的统计输出格式|Invalid statistics output format: {self.stats_output}")

        # 验证统计层级|Validate statistics levels
        valid_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        for level in self.stats_by:
            if level not in valid_levels:
                errors.append(f"无效的统计层级|Invalid statistics level: {level}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
