"""
GATK Joint Genotyping 配置管理模块|GATK Joint Genotyping Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import get_domain_tool_path, expand_path

@dataclass
class JointConfig:
    """Joint Genotyping 配置类|Joint Genotyping Configuration Class"""
    
    # 输入输出|Input/Output
    input_dir: str
    reference: str
    output_dir: str = './joint_genotyping_output'
    
    # 计算资源|Computing resources
    threads: int = 12
    memory: str = '100G'  # Java内存设置|Java memory setting
    
    # 区间设置|Interval settings
    intervals: Optional[str] = None  # 染色体或区间文件|Chromosome or interval file
    
    # 过滤参数 - SNP|Filtering parameters - SNP
    snp_qd: float = 2.0          # Quality by Depth
    snp_fs: float = 60.0         # Fisher Strand
    snp_mq: float = 40.0         # Mapping Quality
    snp_mqrs: float = -12.5      # MappingQualityRankSum
    snp_rprs: float = -8.0       # ReadPosRankSum
    snp_sor: float = 3.0         # StrandOddsRatio
    
    # 过滤参数 - INDEL|Filtering parameters - INDEL
    indel_qd: float = 2.0
    indel_fs: float = 200.0
    indel_rprs: float = -20.0
    indel_sor: float = 10.0
    
    # 工具路径(功能域环境自动解析, 回退裸命令名靠PATH)
    # |Tool paths (auto domain env, fallback to bare name via PATH)
    gatk_path: str = field(default_factory=lambda: get_domain_tool_path(
        'gatk', 'gatk', 'GATK_PATH'))
    bcftools_path: str = field(default_factory=lambda: get_domain_tool_path(
        'bcftools', 'bcftools', 'BCFTOOLS_PATH'))
    samtools_path: str = field(default_factory=lambda: get_domain_tool_path(
        'samtools', 'samtools', 'SAMTOOLS_PATH'))
    
    # 内部属性|Internal attributes
    base_name: str = 'joint_genotyping'
    file_type: str = 'unknown'  # 'gvcf' or 'vcf'
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 展开工具路径|Expand tool paths
        self.gatk_path = expand_path(self.gatk_path)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 展开工具路径|Expand tool paths
        self.gatk_path = expand_path(self.gatk_path)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 展开工具路径|Expand tool paths
        self.gatk_path = expand_path(self.gatk_path)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 展开工具路径|Expand tool paths
        self.gatk_path = expand_path(self.gatk_path)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.samtools_path = expand_path(self.samtools_path)
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查输入目录|Check input directory
        if not os.path.exists(self.input_dir):
            errors.append(f" 输入目录不存在|Input directory does not exist: {self.input_dir}")
        
        # 检查参考基因组|Check reference genome
        if not os.path.exists(self.reference):
            errors.append(f" 参考基因组不存在|Reference genome does not exist: {self.reference}")
        
        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f" 线程数必须为正整数|Thread count must be positive: {self.threads}")
        
        # 检查区间文件|Check interval file
        if self.intervals and not os.path.exists(self.intervals):
            errors.append(f" 区间文件不存在|Interval file does not exist: {self.intervals}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
