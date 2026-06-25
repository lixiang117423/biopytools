"""
连锁不平衡热图配置管理模块|LD Heatmap Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import get_tool_path, expand_path


@dataclass
class LDBlockShowConfig:
    """连锁不平衡热图配置类|LD Heatmap Configuration Class"""

    # 必需参数|Required parameters
    output_prefix: str
    region: str  # 格式|format: chr:start-end

    # 输入文件（三选一，互斥）|Input files (one of three, mutually exclusive)
    vcf_file: Optional[str] = None  # VCF格式输入|VCF format input
    in_genotype: Optional[str] = None  # Genotype格式输入|Genotype format input
    in_plink: Optional[str] = None  # Plink格式输入(bed+bim+fam或ped+map前缀)|Plink format input prefix

    # 软件路径|Software path
    ldblockshow_path: str = field(
        default_factory=lambda: get_tool_path(
            'ldblockshow', '~/software/LDBlockShow/bin/LDBlockShow', 'LDBLOCKSHOW_PATH'
        )
    )

    # LD度量选择|LD statistic selection
    sele_var: int = 1  # 1: D', 2: R², 3/4: Both

    # 过滤参数|Filter parameters
    maf: float = 0.05  # 最小次要等位基因频率|Minimum minor allele frequency
    miss: float = 0.25  # 最大缺失率|Maximum missing ratio
    hwe: float = 0.0  # Hardy-Weinberg平衡P值阈值|Hardy-Weinberg equilibrium P-value threshold
    het: float = 1.0  # 最大杂合率|Maximum heterozygosity ratio
    enable_oth_var: bool = False  # 允许indel/SV/CNV变异|Allow bi-indel bi-sv bi-cnv variants

    # Block检测参数|Block detection parameters
    block_type: int = 1  # 1: Gabriel, 2: SolidSpine, 3: BlockCut, 4: FixBlock, 5: NoBlock
    block_cut: str = "0.85:0.90"  # BlockType3的cutoff|Cutoff for BlockType3
    fix_block: Optional[str] = None  # 固定block文件|Fixed block file

    # 可视化参数|Visualization parameters
    in_gwas: Optional[str] = None  # GWAS P值文件|GWAS P-value file
    in_gff: Optional[str] = None  # GFF3注释文件|GFF3 annotation file
    mer_min_snp_num: int = 50  # 合并网格的最小SNP数|Minimum SNP number to merge grids

    # 输出格式|Output format
    out_png: bool = True  # 输出PNG格式|Output PNG format
    out_pdf: bool = False  # 输出PDF格式|Output PDF format

    # 其他参数|Other parameters
    sub_pop: Optional[str] = None  # 亚群样本文件|Subgroup sample file
    tag_snp_cut: float = 0.80  # TagSNP的LD cutoff|LD cutoff for TagSNP

    # ShowLDSVG绘图参数|ShowLDSVG drawing parameters
    cutline: float = 5.0  # GWAS P值显著性阈值(-log10)|GWAS P-value significance cutoff (-log10)
    point_size: Optional[int] = None  # GWAS散点大小|GWAS point size
    top_site: Optional[str] = None  # 指定GWAS峰值位点(chr:pos)|Specify GWAS peak site (chr:pos)
    no_log_p: bool = False  # 不对P值取-log10|Do not -log10 transform P-value
    no_gene_name: bool = False  # 不显示基因名|Do not show gene names
    show_num: bool = False  # 在热图中显示R²/D'值|Show R²/D' values in heatmap
    spe_snp_name: Optional[str] = None  # 特殊SNP名称文件|Special SNP name file
    show_gwas_spe_snp: bool = False  # 在GWAS图中显示特殊SNP名称|Show special SNP names in GWAS plot
    resize_h: Optional[int] = None  # 图像高度，宽度按比例自动调整|Image height, width auto-adjusted
    no_show_ldist: Optional[int] = None  # 超过此距离的SNP对不显示LD|NoShow pairwise LD over this distance

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        #  关键：展开所有包含~的路径（规范11.3.1要求）|CRITICAL: Expand all paths with ~ (spec 11.3.1)
        self.ldblockshow_path = expand_path(self.ldblockshow_path)

        # 从output_prefix提取输出目录|Extract output directory from output_prefix
        output_dir = os.path.dirname(os.path.abspath(self.output_prefix))
        self.output_path = Path(output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        if self.vcf_file:
            self.vcf_file = os.path.normpath(os.path.abspath(expand_path(self.vcf_file)))
        if self.in_genotype:
            self.in_genotype = os.path.normpath(os.path.abspath(expand_path(self.in_genotype)))
        if self.in_plink:
            self.in_plink = os.path.normpath(os.path.abspath(expand_path(self.in_plink)))
        self.ldblockshow_path = os.path.normpath(self.ldblockshow_path)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件（三选一）|Check input files (one required)
        input_count = sum(bool(x) for x in [self.vcf_file, self.in_genotype, self.in_plink])
        if input_count == 0:
            errors.append("必须指定输入文件|Must specify one input: --vcf-file, --in-genotype, or --in-plink")
        elif input_count > 1:
            errors.append("输入文件只能指定一个|Only one input allowed: --vcf-file, --in-genotype, or --in-plink")
        else:
            if self.vcf_file and not os.path.exists(self.vcf_file):
                errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")
            if self.in_genotype and not os.path.exists(self.in_genotype):
                errors.append(f"Genotype文件不存在|Genotype file does not exist: {self.in_genotype}")
            if self.in_plink and not os.path.exists(self.in_plink):
                errors.append(f"Plink文件不存在|Plink file does not exist: {self.in_plink}")

        if not os.path.exists(self.ldblockshow_path):
            errors.append(f"LDBlockShow可执行文件不存在|LDBlockShow executable does not exist: {self.ldblockshow_path}")

        # 检查可选文件|Check optional files
        if self.in_gwas and not os.path.exists(self.in_gwas):
            errors.append(f"GWAS文件不存在|GWAS file does not exist: {self.in_gwas}")

        if self.in_gff and not os.path.exists(self.in_gff):
            errors.append(f"GFF文件不存在|GFF file does not exist: {self.in_gff}")

        if self.fix_block and not os.path.exists(self.fix_block):
            errors.append(f"FixBlock文件不存在|FixBlock file does not exist: {self.fix_block}")

        if self.sub_pop and not os.path.exists(self.sub_pop):
            errors.append(f"SubPop文件不存在|SubPop file does not exist: {self.sub_pop}")

        # 检查参数范围|Check parameter ranges
        if self.sele_var not in [1, 2, 3, 4]:
            errors.append(f"SeleVar必须是1/2/3/4|SeleVar must be 1/2/3/4: {self.sele_var}")

        if self.block_type not in [1, 2, 3, 4, 5]:
            errors.append(f"BlockType必须是1/2/3/4/5|BlockType must be 1/2/3/4/5: {self.block_type}")

        if self.maf < 0 or self.maf > 1:
            errors.append(f"MAF必须在0-1之间|MAF must be between 0-1: {self.maf}")

        if self.miss < 0 or self.miss > 1:
            errors.append(f"Miss必须在0-1之间|Miss must be between 0-1: {self.miss}")

        if self.hwe < 0 or self.hwe > 1:
            errors.append(f"HWE必须在0-1之间|HWE must be between 0-1: {self.hwe}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
