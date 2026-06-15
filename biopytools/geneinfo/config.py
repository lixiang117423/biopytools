"""
GFF3工具配置管理模块|GFF3 Utilities Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Set, List


GFF_EXTENSIONS = ('.gff3.gz', '.gff.gz', '.gff3', '.gff')


def _find_gff_files(directory: str) -> List[str]:
    """扫描目录下的GFF3/GFF文件|Scan directory for GFF3/GFF files"""
    gff_files = []
    for entry in sorted(os.listdir(directory)):
        full_path = os.path.join(directory, entry)
        if not os.path.isfile(full_path):
            continue
        lower = entry.lower()
        if any(lower.endswith(ext) for ext in GFF_EXTENSIONS):
            gff_files.append(os.path.normpath(os.path.abspath(full_path)))
    return gff_files


def get_sample_name(file_path: str) -> str:
    """从GFF文件路径提取样品名|Extract sample name from GFF file path"""
    basename = os.path.basename(file_path)
    lower = basename.lower()
    for ext in GFF_EXTENSIONS:
        if lower.endswith(ext):
            return basename[:-len(ext)]
    return basename


@dataclass
class GFFConfig:
    """GFF3处理配置类|GFF3 Processing Configuration Class"""

    # 必需参数|Required parameters
    gff3_file: str
    output_file: str

    # 可选参数|Optional parameters
    transcript_types: Set[str] = None
    gene_type: str = 'gene'

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    # 高级选项|Advanced options
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    threads: int = 1  # 线程数|Number of threads

    # 运行时属性|Runtime attributes
    gff3_files: List[str] = field(default_factory=list, init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if self.transcript_types is None:
            self.transcript_types = {'mRNA', 'transcript'}

        # 标准化路径|Normalize paths
        self.gff3_file = os.path.normpath(os.path.abspath(self.gff3_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # 收集GFF3文件列表|Collect GFF3 files list
        if os.path.isdir(self.gff3_file):
            self.gff3_files = _find_gff_files(self.gff3_file)
        else:
            self.gff3_files = [self.gff3_file]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.gff3_file):
            errors.append(f"GFF3文件或目录不存在|GFF3 file or directory does not exist: {self.gff3_file}")

        if not self.gff3_files:
            errors.append(f"目录中未找到GFF3/GFF文件|No GFF3/GFF files found in directory: {self.gff3_file}")

        if not self.transcript_types:
            errors.append("转录本类型不能为空|Transcript types cannot be empty")

        if not self.gene_type:
            errors.append("基因类型不能为空|Gene type cannot be empty")

        if errors:
            raise ValueError("\n".join(errors))

        return True
