"""
GFFcompare两两比较配置管理模块|GFFcompare Pairwise Comparison Configuration Module
"""

import os
import glob
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path


@dataclass
class GffCompareConfig:
    """GFFcompare两两比较配置类|GFFcompare Pairwise Comparison Configuration Class"""

    # 输入文件|Input files
    input_files: List[str] = field(default_factory=list)
    output_dir: str = './gffcompare_output'

    # GFFcompare参数|GFFcompare parameters
    exon_range: Optional[int] = None       # -e
    tss_distance: Optional[int] = None      # -d
    discard_single_exon_query: bool = False # -M
    discard_single_exon_ref: bool = False   # -N
    ref_overlap_only: bool = False          # -R
    query_overlap_only: bool = False        # -Q
    no_tmap_refmap: bool = False            # -T
    strict_match: bool = False             # --strict-match
    cds_match: bool = False                # --cds-match
    genome_seq: Optional[str] = None       # -s
    cprefix: Optional[str] = None          # -p
    verbose_mode: bool = False             # -V

    # 运行参数|Runtime parameters
    force: bool = False

    # 工具路径|Tool path
    gffcompare_path: str = field(
        default_factory=lambda: get_tool_path(
            'gffcompare',
            '~/miniforge3/envs/gffcompare_v.0.12.10/bin/gffcompare'
        )
    )

    # 内部属性|Internal attributes
    output_path: Path = field(default=None, init=False)
    resolved_files: List[str] = field(default_factory=list, init=False)
    sample_names: List[str] = field(default_factory=list, init=False)

    VALID_EXTENSIONS = ('.gff', '.gtf', '.gff3', '.gff.gz', '.gtf.gz', '.gff3.gz')
    MAX_FILES = 5
    MIN_FILES = 2

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.output_dir = expand_path(self.output_dir)
        self.gffcompare_path = expand_path(self.gffcompare_path)
        if self.genome_seq:
            self.genome_seq = expand_path(self.genome_seq)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 解析输入文件|Resolve input files
        self.resolved_files = self._resolve_input_files(self.input_files)
        self.sample_names = [self._extract_sample_name(f) for f in self.resolved_files]

        # 创建输出子目录|Create output subdirectories
        for subdir in ('00_pipeline_info', '01_gffcompare', '02_summary', '99_logs'):
            (self.output_path / subdir).mkdir(parents=True, exist_ok=True)

    def _resolve_input_files(self, input_files: List[str]) -> List[str]:
        """解析输入文件列表，支持文件和目录|Resolve input files, support files and directories"""
        resolved = []
        for item in input_files:
            item = os.path.normpath(os.path.abspath(expand_path(item)))

            if os.path.isfile(item):
                self._check_extension(item)
                resolved.append(item)
            elif os.path.isdir(item):
                found = []
                for ext in self.VALID_EXTENSIONS:
                    found.extend(glob.glob(os.path.join(item, f'*{ext}')))
                found = sorted(set(found))
                if not found:
                    raise ValueError(
                        f"目录中未找到GFF/GTF文件|No GFF/GTF files found in directory: {item}"
                    )
                for f in found:
                    self._check_extension(f)
                print(f"在目录中找到{len(found)}个GFF/GTF文件|Found {len(found)} GFF/GTF file(s) in directory: {item}")
                resolved.extend(found)
            else:
                raise ValueError(f"文件或目录不存在|File or directory does not exist: {item}")

        if not resolved:
            raise ValueError("未找到有效的GFF/GTF输入文件|No valid GFF/GTF input files found")

        return sorted(set(resolved))

    def _check_extension(self, file_path: str):
        """检查文件扩展名|Check file extension"""
        lower = file_path.lower()
        if not any(lower.endswith(ext) for ext in self.VALID_EXTENSIONS):
            raise ValueError(
                f"不支持的文件格式，仅支持GFF/GTF文件|Unsupported file format, "
                f"only GFF/GTF files supported: {file_path}"
            )

    @staticmethod
    def _extract_sample_name(file_path: str) -> str:
        """从文件路径提取样本名称|Extract sample name from file path"""
        basename = os.path.basename(file_path)
        for ext in ('.gff3.gz', '.gtf.gz', '.gff.gz', '.gff3', '.gtf', '.gff'):
            if basename.lower().endswith(ext):
                return basename[:-len(ext)]
        return Path(basename).stem

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if len(self.resolved_files) < self.MIN_FILES:
            errors.append(
                f"至少需要{self.MIN_FILES}个输入文件进行比较|"
                f"At least {self.MIN_FILES} input files required for comparison, "
                f"got {len(self.resolved_files)}"
            )

        if len(self.resolved_files) > self.MAX_FILES:
            errors.append(
                f"输入文件不能超过{self.MAX_FILES}个|"
                f"Input files cannot exceed {self.MAX_FILES}, "
                f"got {len(self.resolved_files)}"
            )

        # 检查文件存在|Check files exist
        for f in self.resolved_files:
            if not os.path.exists(f):
                errors.append(f"输入文件不存在|Input file does not exist: {f}")

        # 检查样本名唯一性|Check sample name uniqueness
        seen = set()
        for name in self.sample_names:
            if name in seen:
                errors.append(f"重复的样本名称|Duplicate sample name: {name}")
            seen.add(name)

        # 检查可选参数|Check optional parameters
        if self.exon_range is not None and self.exon_range <= 0:
            errors.append(f"exon_range必须为正整数|exon_range must be positive: {self.exon_range}")

        if self.tss_distance is not None and self.tss_distance <= 0:
            errors.append(f"tss_distance必须为正整数|tss_distance must be positive: {self.tss_distance}")

        if self.genome_seq and not os.path.exists(self.genome_seq):
            errors.append(f"基因组序列文件不存在|Genome sequence file not found: {self.genome_seq}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
