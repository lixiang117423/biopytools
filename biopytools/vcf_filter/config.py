"""
VCF文件筛选配置管理模块 | VCF File Filtering Configuration Management Module
"""

import os
from pathlib import Path
from typing import Optional, List, Union

class VCFFilterConfig:
    """VCF筛选配置类 | VCF Filtering Configuration Class"""
    
    def __init__(self, vcf_file: str, output_file: Optional[str] = None,
                 chr_name: Optional[Union[str, List[str]]] = None,
                 start: Optional[int] = None, end: Optional[int] = None,
                 convert_format: bool = False, plink_path: str = "plink",
                 allow_extra_chr: bool = True, min_maf: Optional[float] = None,
                 max_missing: Optional[float] = None, 
                 quality_threshold: Optional[float] = None,
                 min_depth: Optional[int] = None, max_depth: Optional[int] = None,
                 keep_samples: Optional[List[str]] = None,
                 remove_samples: Optional[List[str]] = None,
                 keep_ids: Optional[List[str]] = None,
                 remove_ids: Optional[List[str]] = None,
                 biallelic_only: bool = False, remove_indels: bool = False,
                 skip_validation: bool = True, verbose: bool = False):
        
        # 输入输出文件 | Input/Output files
        self.vcf_file = os.path.normpath(os.path.abspath(vcf_file))
        
        if output_file is None:
            base_name = Path(vcf_file).stem
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            output_file = f"{base_name}_filtered.vcf"
        
        self.output_file = os.path.normpath(os.path.abspath(output_file))
        self.output_dir = Path(self.output_file).parent
        self.base_name = Path(self.output_file).stem
        self.output_path = Path(self.output_file)
        
        # 位置筛选参数 | Position filtering parameters
        self.chr_name = chr_name
        self.start = start
        self.end = end
        
        # 格式转换参数 | Format conversion parameters
        self.convert_format = convert_format
        self.plink_path = plink_path
        self.allow_extra_chr = allow_extra_chr
        
        # 质量控制参数 | Quality control parameters
        self.min_maf = min_maf
        self.max_missing = max_missing
        self.quality_threshold = quality_threshold
        self.min_depth = min_depth
        self.max_depth = max_depth
        
        # 样本筛选参数 | Sample filtering parameters
        self.keep_samples = keep_samples
        self.remove_samples = remove_samples
        
        # 变异位点筛选参数 | Variant filtering parameters
        self.keep_ids = keep_ids
        self.remove_ids = remove_ids
        self.biallelic_only = biallelic_only
        self.remove_indels = remove_indels
        
        # 性能优化参数 | Performance optimization parameters
        self.skip_validation = skip_validation
        self.verbose = verbose
        
        # 创建输出目录 | Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 验证配置 | Validate configuration
        if not self.skip_validation:
            self.validate()
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查位置参数 | Check position parameters
        if self.start is not None and self.start < 0:
            errors.append(f"起始位置必须为非负整数 | Start position must be non-negative: {self.start}")
        
        if self.end is not None and self.end < 0:
            errors.append(f"结束位置必须为非负整数 | End position must be non-negative: {self.end}")
        
        if (self.start is not None and self.end is not None and 
            self.start > self.end):
            errors.append(f"起始位置不能大于结束位置 | Start position cannot be greater than end position: {self.start} > {self.end}")
        
        # 检查质控参数 | Check QC parameters
        if self.min_maf is not None and not 0 <= self.min_maf <= 0.5:
            errors.append(f"最小MAF必须在0-0.5之间 | Minimum MAF must be between 0-0.5: {self.min_maf}")
        
        if self.max_missing is not None and not 0 <= self.max_missing <= 1:
            errors.append(f"最大缺失率必须在0-1之间 | Maximum missing rate must be between 0-1: {self.max_missing}")
        
        if self.quality_threshold is not None and self.quality_threshold < 0:
            errors.append(f"质量阈值必须为非负数 | Quality threshold must be non-negative: {self.quality_threshold}")
        
        if errors:
            raise ValueError("\\n".join(errors))
        
        return True
