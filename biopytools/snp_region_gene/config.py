"""SNP区域基因提取配置类|SNP Region Gene Extractor Configuration"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class SnpRegionConfig:
    """SNP区域基因提取配置类|SNP Region Gene Extractor Configuration"""

    # 必需参数|Required parameters
    snp_file: str
    gff_file: str
    genome_file: str

    # 可选参数|Optional parameters
    left: int = 0
    right: int = 0
    promoter: int = 2000
    output_prefix: str = "./snp_region_output"
    gffread_path: str = "gffread"
    seqkit_path: str = "seqkit"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        output_dir = str(Path(self.output_prefix).parent)
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 设置输出文件路径|Set output file paths
        base_path = self.output_prefix
        self.cds_output = f"{base_path}_cds.fasta"
        self.protein_output = f"{base_path}_protein.fasta"
        self.gene_list_output = f"{base_path}_gene_list.txt"

        # 临时文件路径|Temporary file paths
        self.temp_cds = f"{base_path}_temp_all_cds.fasta"
        self.temp_protein = f"{base_path}_temp_all_protein.fasta"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查文件存在性|Check file existence
        if not os.path.exists(self.snp_file):
            errors.append(f"SNP文件不存在|SNP file not found: {self.snp_file}")

        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file not found: {self.gff_file}")

        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_file}")

        # 检查参数有效性|Check parameter validity
        if self.left < 0:
            errors.append(f"上游距离不能为负数|Left distance cannot be negative: {self.left}")

        if self.right < 0:
            errors.append(f"下游距离不能为负数|Right distance cannot be negative: {self.right}")

        if self.promoter < 0:
            errors.append(f"启动子距离不能为负数|Promoter distance cannot be negative: {self.promoter}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
