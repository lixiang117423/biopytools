"""
VCF系统发育分析配置管理模块|VCF Phylogenetic Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path

@dataclass
class PhyloConfig:
    """系统发育分析配置类|Phylogenetic Analysis Configuration Class"""

    # 输入文件|Input files
    vcf_file: Optional[str] = None
    distance_matrix: Optional[str] = None

    # 输出文件|Output files
    output_dir: str = "."  # 输出目录|Output directory
    output_prefix: str = "phylo_analysis"  # 输出文件前缀|Output file prefix
    tree_output: Optional[str] = None

    # 工具路径|Tool paths
    vcf2dis_path: str = "VCF2Dis"
    nw_reroot_path: str = "~/miniforge3/envs/newick_utils_v.1.6/bin/nw_reroot"

    # 重根化参数|Rerooting parameters
    outgroup: Optional[str] = None  # 外群样本标签（逗号分隔）|Outgroup sample labels (comma-separated)

    # 行为控制|Behavior control
    skip_vcf2dis: bool = False

    # 工作目录|Working directory
    working_dir: str = "."

    # 内部属性|Internal attributes
    base_name: str = "vcf_phylo"
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.nw_reroot_path = expand_path(self.nw_reroot_path)

        self.working_path = Path(self.working_dir).resolve()

        # 设置路径|Setup paths
        self._setup_paths()

        # 标准化路径|Normalize paths
        if self.vcf_file:
            self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))

        if self.distance_matrix:
            self.distance_matrix = os.path.normpath(os.path.abspath(self.distance_matrix))

        # 处理外群参数|Process outgroup parameter
        self.outgroup_list = []
        if self.outgroup:
            # 支持逗号分隔的多个外群|Support comma-separated multiple outgroups
            self.outgroup_list = [og.strip() for og in self.outgroup.split(',') if og.strip()]

        # 设置重根化输出文件路径|Setup rerooted output file path
        if self.outgroup_list:
            base_name = Path(self.tree_output).stem
            output_dir = Path(self.tree_output).parent
            outgroup_suffix = "_".join(self.outgroup_list)
            self.rerooted_tree_output = str(output_dir / f"{base_name}_{outgroup_suffix}_root.nwk")
        else:
            self.rerooted_tree_output = None
        
    def _setup_paths(self):
        """设置路径|Setup paths"""
        # 创建输出目录|Create output directory
        self.output_dir_path = Path(self.output_dir).resolve()
        self.output_dir_path.mkdir(parents=True, exist_ok=True)

        # 如果没有指定距离矩阵文件，使用输出目录+前缀|If distance matrix not specified, use output dir + prefix
        if not self.distance_matrix:
            self.distance_matrix = str(self.output_dir_path / f"{self.output_prefix}.dis")

        # 如果没有指定树输出文件，使用输出目录+前缀|If tree output not specified, use output dir + prefix
        if not self.tree_output:
            self.tree_output = str(self.output_dir_path / f"{self.output_prefix}.nwk")
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not self.skip_vcf2dis:
            if not self.vcf_file:
                errors.append("VCF文件路径必须指定|VCF file path must be specified")
            elif not os.path.exists(self.vcf_file):
                errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")
        else:
            if not self.distance_matrix:
                errors.append("跳过VCF2Dis时必须指定距离矩阵文件|Distance matrix file must be specified when skipping VCF2Dis")
            elif not os.path.exists(self.distance_matrix):
                errors.append(f"距离矩阵文件不存在|Distance matrix file does not exist: {self.distance_matrix}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
