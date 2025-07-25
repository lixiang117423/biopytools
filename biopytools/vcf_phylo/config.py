"""
VCF系统发育分析配置管理模块 | VCF Phylogenetic Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PhyloConfig:
    """系统发育分析配置类 | Phylogenetic Analysis Configuration Class"""
    
    # 输入文件 | Input files
    vcf_file: Optional[str] = None
    distance_matrix: Optional[str] = None
    
    # 输出文件 | Output files
    output_prefix: str = "phylo_analysis"
    tree_output: Optional[str] = None
    
    # 工具路径 | Tool paths
    vcf2dis_path: str = "VCF2Dis"
    
    # 行为控制 | Behavior control
    skip_vcf2dis: bool = False
    
    # 工作目录 | Working directory
    working_dir: str = "."
    
    # 内部属性 | Internal attributes
    base_name: str = "vcf_phylo"
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.working_path = Path(self.working_dir).resolve()
        
        # 设置路径 | Setup paths
        self._setup_paths()
        
        # 标准化路径 | Normalize paths
        if self.vcf_file:
            self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        
        if self.distance_matrix:
            self.distance_matrix = os.path.normpath(os.path.abspath(self.distance_matrix))
        
        self.tree_output = os.path.normpath(os.path.abspath(self.tree_output))
        
    def _setup_paths(self):
        """设置路径 | Setup paths"""
        # 如果没有指定距离矩阵文件，使用输出前缀 | If distance matrix not specified, use output prefix
        if not self.distance_matrix:
            self.distance_matrix = self.output_prefix
            
        # 如果没有指定树输出文件，使用输出前缀 | If tree output not specified, use output prefix
        if not self.tree_output:
            self.tree_output = f"{self.output_prefix}.nwk"
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input files
        if not self.skip_vcf2dis:
            if not self.vcf_file:
                errors.append("VCF文件路径必须指定 | VCF file path must be specified")
            elif not os.path.exists(self.vcf_file):
                errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        else:
            if not self.distance_matrix:
                errors.append("跳过VCF2Dis时必须指定距离矩阵文件 | Distance matrix file must be specified when skipping VCF2Dis")
            elif not os.path.exists(self.distance_matrix):
                errors.append(f"距离矩阵文件不存在 | Distance matrix file does not exist: {self.distance_matrix}")
        
        # 检查输出目录 | Check output directories
        output_dir = Path(self.distance_matrix).parent
        tree_output_dir = Path(self.tree_output).parent
        
        for dir_path, desc in [(output_dir, "距离矩阵输出目录"), (tree_output_dir, "系统发育树输出目录")]:
            if not dir_path.exists():
                try:
                    dir_path.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    errors.append(f"{desc}创建失败 | Failed to create {desc}: {dir_path}, 错误: {e}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
