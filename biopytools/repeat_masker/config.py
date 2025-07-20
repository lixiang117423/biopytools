"""
重复序列分析配置管理模块 | Repeat Sequence Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class RepeatConfig:
    """重复序列分析配置类 | Repeat Sequence Analysis Configuration Class"""
    
    # 必需文件 | Required files
    genome_file: str
    output_dir: str = './repeat_output'
    
    # 物种和数据库配置 | Species and database configuration
    species: str = "human"
    custom_lib: Optional[str] = None  # 自定义重复序列库 | Custom repeat library
    
    # RepeatMasker参数 | RepeatMasker parameters
    rm_threads: int = 8  # 线程数 | Number of threads
    rm_soft_mask: bool = True  # 软屏蔽（小写） | Soft masking (lowercase)
    rm_generate_gff: bool = True  # 生成GFF文件 | Generate GFF file
    rm_exclude_low_complexity: bool = True  # 排除低复杂度序列 | Exclude low complexity sequences
    
    # RepeatModeler参数 | RepeatModeler parameters
    rm_model_threads: int = 8  # RepeatModeler线程数 | RepeatModeler threads
    rm_use_ltr_struct: bool = True  # 使用LTRStruct | Use LTRStruct
    rm_run_modeler: bool = True  # 是否运行RepeatModeler | Whether to run RepeatModeler
    
    # EDTA参数 | EDTA parameters
    edta_threads: int = 8  # EDTA线程数 | EDTA thread count
    edta_species: str = "others"  # EDTA物种类型 | EDTA species type (rice, maize, others)
    edta_step: str = "all"  # EDTA执行步骤 | EDTA execution step (all, filter, final, anno)
    edta_sensitive: int = 0  # 敏感模式 | Sensitive mode (0=default, 1=sensitive)
    edta_anno: int = 1  # 注释模式 | Annotation mode (0=no, 1=yes)
    edta_evaluate: int = 1  # 评估模式 | Evaluation mode (0=no, 1=yes)
    edta_overwrite: int = 0  # 覆盖现有结果 | Overwrite existing results (0=no, 1=yes)
    edta_run_analysis: bool = True  # 是否运行EDTA分析 | Whether to run EDTA analysis
    
    # TRF参数 | TRF parameters
    trf_match_weight: int = 2  # 匹配权重 | Match weight
    trf_mismatch_penalty: int = 7  # 错配惩罚 | Mismatch penalty
    trf_indel_penalty: int = 7  # 插入缺失惩罚 | Indel penalty
    trf_min_score: int = 80  # 最小比对得分 | Minimum alignment score
    trf_max_period: int = 10  # 最大周期长度 | Maximum period size
    trf_min_copies: int = 50  # 最小重复次数 | Minimum copy number
    trf_max_length: int = 500  # 最大期望长度 | Maximum expected length
    trf_run_analysis: bool = True  # 是否运行TRF分析 | Whether to run TRF analysis
    
    # 分析选项 | Analysis options
    run_de_novo: bool = True  # 运行从头预测 | Run de novo prediction
    run_database_search: bool = True  # 运行数据库搜索 | Run database search
    run_tandem_repeats: bool = True  # 运行串联重复分析 | Run tandem repeat analysis
    run_edta: bool = True  # 运行EDTA转座元件分析 | Run EDTA transposable element analysis
    
    # 过滤参数 | Filtering parameters
    min_repeat_length: int = 50  # 最小重复序列长度 | Minimum repeat length
    max_divergence: float = 0.25  # 最大分化度 | Maximum divergence
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.custom_lib:
            self.custom_lib = os.path.normpath(os.path.abspath(self.custom_lib))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查自定义库文件 | Check custom library file
        if self.custom_lib and not os.path.exists(self.custom_lib):
            errors.append(f"自定义重复序列库不存在 | Custom repeat library does not exist: {self.custom_lib}")
        
        # 检查参数范围 | Check parameter ranges
        if self.rm_threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.rm_threads}")
        
        if self.edta_threads <= 0:
            errors.append(f"EDTA线程数必须为正整数 | EDTA thread count must be positive integer: {self.edta_threads}")
        
        if self.min_repeat_length <= 0:
            errors.append(f"最小重复序列长度必须为正整数 | Minimum repeat length must be positive: {self.min_repeat_length}")
        
        if not 0 <= self.max_divergence <= 1:
            errors.append(f"最大分化度必须在0-1之间 | Maximum divergence must be between 0-1: {self.max_divergence}")
        
        # 检查EDTA参数 | Check EDTA parameters
        valid_edta_species = ["rice", "maize", "others"]
        if self.edta_species not in valid_edta_species:
            errors.append(f"无效的EDTA物种类型 | Invalid EDTA species type: {self.edta_species}. Valid options: {valid_edta_species}")
        
        valid_edta_steps = ["all", "filter", "final", "anno"]
        if self.edta_step not in valid_edta_steps:
            errors.append(f"无效的EDTA步骤 | Invalid EDTA step: {self.edta_step}. Valid options: {valid_edta_steps}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
