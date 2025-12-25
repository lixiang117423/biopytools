"""
配置管理模块 | Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

@dataclass
class FilterConfig:
    """变异提取配置类 | Variant Extraction Configuration Class"""
    
    # 输入文件 | Input files
    gff_file: str
    exonic_file: str
    all_variant_file: str
    gene_id: Optional[str] = None
    gene_list_file: Optional[str] = None
    
    # 输出参数 | Output parameters
    output_dir: str = './filter_output'
    output_format: str = 'excel'  # 'excel' or 'txt'
    
    # 分析参数 | Analysis parameters
    extend_bp: int = 5000  # 上下游扩展范围 | Upstream/downstream extension
    threads: int = 88  # 线程数 | Number of threads
    
    # 内部属性 | Internal attributes
    gene_ids: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 创建输出目录 | Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))
        self.exonic_file = os.path.normpath(os.path.abspath(self.exonic_file))
        self.all_variant_file = os.path.normpath(os.path.abspath(self.all_variant_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.gene_list_file:
            self.gene_list_file = os.path.normpath(os.path.abspath(self.gene_list_file))
        
        # 解析基因ID列表 | Parse gene ID list
        self._parse_gene_ids()
    
    def _parse_gene_ids(self):
        """解析基因ID | Parse gene IDs"""
        if self.gene_id:
            self.gene_ids = [self.gene_id]
        elif self.gene_list_file:
            with open(self.gene_list_file, 'r') as f:
                self.gene_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        else:
            raise ValueError("必须指定基因ID或基因列表文件 | Must specify gene ID or gene list file")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查GFF文件 | Check GFF file
        if not os.path.exists(self.gff_file):
            errors.append(f"❌ GFF文件不存在 | GFF file does not exist: {self.gff_file}")
        
        # 检查外显子变异文件 | Check exonic variant file
        if not os.path.exists(self.exonic_file):
            errors.append(f"❌ 外显子变异文件不存在 | Exonic variant file does not exist: {self.exonic_file}")
        
        # 检查所有变异文件 | Check all variant file
        if not os.path.exists(self.all_variant_file):
            errors.append(f"❌ 所有变异文件不存在 | All variant file does not exist: {self.all_variant_file}")
        
        # 检查基因列表文件 | Check gene list file
        if self.gene_list_file and not os.path.exists(self.gene_list_file):
            errors.append(f"❌ 基因列表文件不存在 | Gene list file does not exist: {self.gene_list_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.extend_bp < 0:
            errors.append(f"❌ 扩展范围必须为非负整数 | Extension range must be non-negative: {self.extend_bp}")
        
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.output_format not in ['excel', 'txt']:
            errors.append(f"❌ 输出格式必须是'excel'或'txt' | Output format must be 'excel' or 'txt': {self.output_format}")
        
        if not self.gene_ids:
            errors.append("❌ 未找到任何基因ID | No gene IDs found")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
