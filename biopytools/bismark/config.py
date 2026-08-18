"""
Bismark流程配置管理模块|Bismark Pipeline Configuration Management Module
"""
# (此文件无改动|No changes in this file)
import os
from dataclasses import dataclass, field
from pathlib import Path

from ..common.paths import get_domain_tool_path, expand_path

@dataclass
class BismarkConfig:
    """Bismark流程配置类|Bismark Pipeline Configuration Class"""
    
    raw_dir: str
    genome_fa: str
    output_dir: str
    
    threads: int = 88
    sort_buffer: str = '400G'
    no_overlap: bool = True
    pattern: str = '_1_clean.fq.gz'
    
    # 工具路径|Tool paths (功能域环境自动解析, 回退裸命令名靠PATH)
    # |Tool paths (auto domain env, fallback to bare name via PATH)
    bismark_path: str = field(default_factory=lambda: get_domain_tool_path(
        'bismark', 'bismark', 'BISMARK_PATH'))
    bismark_genome_preparation_path: str = field(default_factory=lambda: get_domain_tool_path(
        'bismark_genome_preparation', 'bismark_genome_preparation', 'BISMARK_GENOME_PREPARATION_PATH'))
    bowtie2_path: str = field(default_factory=lambda: get_domain_tool_path(
        'bowtie2', 'bowtie2', 'BOWTIE2_PATH'))
    bismark_methylation_extractor_path: str = field(default_factory=lambda: get_domain_tool_path(
        'bismark_methylation_extractor', 'bismark_methylation_extractor', 'BISMARK_METHYLATION_EXTRACTOR_PATH'))
    
    genome_dir: str = field(init=False) 
    mapping_dir: str = field(init=False)
    result_dir: str = field(init=False)
    tmp_dir: str = field(init=False)
    
    def __post_init__(self):
        self.raw_dir = os.path.normpath(os.path.abspath(self.raw_dir))
        self.genome_fa = os.path.normpath(os.path.abspath(self.genome_fa))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        self.genome_dir = os.path.dirname(self.genome_fa)
        self.mapping_dir = os.path.join(self.output_dir, "mapping")
        self.result_dir = os.path.join(self.output_dir, "result")
        self.tmp_dir = os.path.join(self.output_dir, "tmp")

        # 展开工具路径|Expand tool paths
        self.bismark_path = expand_path(self.bismark_path)
        self.bismark_genome_preparation_path = expand_path(self.bismark_genome_preparation_path)
        self.bowtie2_path = expand_path(self.bowtie2_path)
        self.bismark_methylation_extractor_path = expand_path(self.bismark_methylation_extractor_path)
    
    def validate(self):
        errors = []
        if not os.path.exists(self.raw_dir):
            errors.append(f"原始数据目录不存在|Raw data directory not found: {self.raw_dir}")
        if not os.path.exists(self.genome_fa):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_fa}")
        if self.threads <= 0:
            errors.append("线程数必须大于0|Number of threads must be greater than 0")
        if not self.pattern.endswith('.gz'):
            errors.append("模式(--pattern)必须以.gz结尾|Pattern (--pattern) must end with .gz")
        
        for dir_path in [self.output_dir, self.genome_dir]:
            try:
                Path(dir_path).mkdir(parents=True, exist_ok=True)
                test_file = Path(dir_path) / ".writable_test"
                test_file.touch()
                test_file.unlink()
            except Exception as e:
                errors.append(f"目录不可写或无法创建|Directory is not writable or cannot be created: {dir_path}, 错误|Error: {e}")
        
        if errors:
            raise ValueError("\n".join(errors))
        return True
