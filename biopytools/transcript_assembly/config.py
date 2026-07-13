"""
转录本组装配置管理模块|Transcript Assembly Configuration Module
"""

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path, get_samtools_path


@dataclass
class TranscriptAssemblyConfig:
    """转录本组装配置类|Transcript Assembly Configuration Class"""

    # 输入|Input (input_dir 与 bam_files 互斥|mutually exclusive)
    output_dir: str = "./transcript_output"
    genome_file: Optional[str] = None
    input_dir: Optional[str] = None
    bam_files: Optional[List[str]] = None

    # 处理参数|Processing parameters
    threads: int = 12
    fastq_pattern: str = '*_1.clean.fq.gz'
    sample_timeout: int = 43200  # 单样本超时秒|sample timeout seconds
    read_type: str = "auto"  # auto|short|long
    guide_gff: Optional[str] = None  # -G 参考注释|reference annotation
    output_transcripts: bool = False  # 额外输出 transcripts.fa|also output cDNA

    # 步骤控制|Step control
    step: Optional[int] = None  # 1-6, None=全部|None=all

    # 日志|Logging
    log_file: Optional[str] = None
    log_level: str = "INFO"

    # 高级|Advanced
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    force: bool = False

    # 内部|Internal
    samples: List[dict] = None

    # 工具路径|Tool paths (§13.6 传完整路径|pass full path)
    stringtie_bin: str = field(default_factory=lambda: get_tool_path(
        'stringtie', '~/miniforge3/envs/RNA_Seq/bin/stringtie', 'STRINGTIE_PATH'))
    gffread_bin: str = field(default_factory=lambda: get_tool_path(
        'gffread', '~/miniforge3/envs/RNA_Seq/bin/gffread', 'GFFREAD_PATH'))
    hisat2_bin: str = field(default_factory=lambda: get_tool_path(
        'hisat2', '~/miniforge3/envs/RNA_Seq/bin/hisat2', 'HISAT2_PATH'))
    hisat2_build_bin: str = field(default_factory=lambda: get_tool_path(
        'hisat2-build', '~/miniforge3/envs/RNA_Seq/bin/hisat2-build', 'HISAT2_BUILD_PATH'))
    samtools_bin: str = field(default_factory=get_samtools_path)

    def __post_init__(self):
        """初始化后处理|Post-init: expand paths, create dirs"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 展开所有工具路径 ~ (§11.3.1)|expand ~ in tool paths
        for attr in ('stringtie_bin', 'gffread_bin', 'hisat2_bin',
                     'hisat2_build_bin', 'samtools_bin'):
            setattr(self, attr, expand_path(getattr(self, attr)))

        # 展开用户路径|expand user paths
        if self.genome_file:
            self.genome_file = os.path.normpath(expand_path(self.genome_file))
        if self.input_dir:
            self.input_dir = os.path.normpath(expand_path(self.input_dir))
        if self.bam_files:
            self.bam_files = [os.path.normpath(expand_path(b)) for b in self.bam_files]
        if self.guide_gff:
            self.guide_gff = os.path.normpath(expand_path(self.guide_gff))

        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 流程元数据 + 日志目录|pipeline_info + logs dirs (§12)
        Path(self.output_dir, '00_pipeline_info').mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置|Validate configuration"""
        errors = []

        # 互斥:input_dir 与 bam_files 恰一|exactly one of input_dir/bam_files
        has_fastq = bool(self.input_dir)
        has_bam = bool(self.bam_files)
        if has_fastq and has_bam:
            errors.append("input_dir 与 bam_files 互斥|input_dir and bam_files are mutually exclusive")
        if not has_fastq and not has_bam:
            errors.append("必须提供 input_dir 或 bam_files|must provide input_dir or bam_files")

        # genome 条件必需|genome conditionally required
        need_genome = has_fastq or self.output_transcripts
        if need_genome and (not self.genome_file or not os.path.exists(self.genome_file)):
            errors.append(f"该模式需要基因组文件|genome file required for this mode: {self.genome_file}")
        elif self.genome_file and not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_file}")

        # 输入存在性|input existence
        if has_fastq and not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")
        if has_bam:
            for b in self.bam_files:
                if not os.path.exists(b):
                    errors.append(f"BAM 文件不存在|BAM file not found: {b}")
        if self.guide_gff and not os.path.exists(self.guide_gff):
            errors.append(f"参考注释文件不存在|Guide GFF not found: {self.guide_gff}")

        # 线程|threads
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        # read_type
        if self.read_type not in ('auto', 'short', 'long'):
            errors.append(f"无效 read_type|Invalid read_type: {self.read_type} (auto|short|long)")

        # step 范围|step range
        if self.step is not None and self.step not in [1, 2, 3, 4, 5, 6]:
            errors.append(f"无效步骤|Invalid step: {self.step} (1-6)")
        # BAM 模式禁止 step 1-3|BAM mode forbids steps 1-3
        if has_bam and self.step is not None and self.step in [1, 2, 3]:
            errors.append(f"BAM 模式不支持步骤|BAM mode does not support step: {self.step} (仅 4-6|only 4-6)")

        # 工具存在性|tool existence (§13.6)
        for name, path in [('stringtie', self.stringtie_bin), ('gffread', self.gffread_bin),
                           ('samtools', self.samtools_bin)]:
            if not os.path.exists(path):
                errors.append(f"工具不存在|Tool not found [{name}]: {path}")
        if has_fastq:  # FASTQ 模式还需 hisat2
            for name, path in [('hisat2', self.hisat2_bin), ('hisat2-build', self.hisat2_build_bin)]:
                if not os.path.exists(path):
                    errors.append(f"工具不存在|Tool not found [{name}]: {path}")

        if errors:
            raise ValueError("\n".join(errors))
        return True
