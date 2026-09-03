"""全长转录本分析配置模块|Full-length RNA analysis configuration

输入嗅探规则|Input sniffing rules:
- *.subreads.bam            -> pacbio_subreads (ccs -> refine -> 引擎)
- *.ccs.bam/*.hifi_reads.bam -> pacbio_ccs     (refine -> 引擎)
- *.fa/fasta/fq/fastq(.gz)  -> reads           (直接引擎,需显式 --data-type pacbio|ont)
"""

import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path, get_tool_path

# Clontech SMARTer cDNA 文库标准引物(Iso-Seq 通用默认)
# |Clontech SMARTer cDNA library standard primers (universal Iso-Seq default)
CLONTECH_PRIMERS_FA = """>primer_5p
AAGCAGTGGTATCAACGCAGAGTACATGGG
>primer_3p
GTACTCTGCGTTGATACCACTGCTT
"""

_SUBREADS_RE = re.compile(r"\.subreads\.bam(\.pbi)?$", re.IGNORECASE)
_CCS_RE = re.compile(r"\.(ccs|hifi_reads|consensusreads)\.bam(\.pbi)?$", re.IGNORECASE)
_READS_RE = re.compile(r"\.(fa|fasta|fq|fastq)(\.gz)?$", re.IGNORECASE)


def detect_input_type(reads: List[str]) -> str:
    """嗅探输入形态,混杂或未知报错|Sniff input kind, error on mixed/unknown

    Returns:
        "pacbio_subreads" | "pacbio_ccs" | "reads"
    """
    kinds = set()
    for r in reads:
        name = os.path.basename(r)
        if _SUBREADS_RE.search(name):
            kinds.add("pacbio_subreads")
        elif _CCS_RE.search(name):
            kinds.add("pacbio_ccs")
        elif _READS_RE.search(name):
            kinds.add("reads")
        else:
            raise ValueError(f"无法识别输入文件类型|Unknown input type: {name}")
    if len(kinds) != 1:
        raise ValueError(f"输入文件类型混杂,不支持混跑|Mixed input types not supported: {sorted(kinds)}")
    return kinds.pop()


@dataclass
class RnaIsoConfig:
    """全长转录本分析配置|Full-length RNA analysis configuration"""

    reads: List[str]                    # 输入文件(可多个,同一样品多个run)
    data_type: Optional[str] = None     # pacbio|ont,仅reads文件时必填
    reference: Optional[str] = None     # 参考基因组FASTA(isoquant引擎必填)
    genedb: Optional[str] = None        # 参考注释GTF/GFF(可选)
    engine: str = "isoquant"            # isoquant|isoseq3|both
    primers: Optional[str] = None       # 引物fasta(None=内置Clontech)
    min_passes: int = 1                 # ccs最小pass数(Iso-Seq官方推荐1)
    threads: int = 12                   # 线程数
    prefix: str = "rna_sample"          # 样本前缀
    output_dir: str = "./rna_iso_output"

    # 工具路径(§11.2)|Tool paths
    ccs_path: str = field(default_factory=lambda: get_tool_path(
        'ccs', '~/miniforge3/envs/isoseq_v.4.0.0/bin/ccs', 'CCS_PATH'))
    isoseq3_path: str = field(default_factory=lambda: get_tool_path(
        'isoseq3', '~/miniforge3/envs/isoseq_v.4.0.0/bin/isoseq3', 'ISOSEQ3_PATH'))
    isoquant_path: str = field(default_factory=lambda: get_tool_path(
        'isoquant', '~/miniforge3/envs/isoseq_v.4.0.0/bin/isoquant', 'ISOQUANT_PATH'))
    samtools_path: str = field(default_factory=lambda: get_tool_path(
        'samtools', '~/miniforge3/envs/isoseq_v.4.0.0/bin/samtools', 'SAMTOOLS_PATH'))

    # 派生属性(__post_init__填充)|Derived attributes (filled in __post_init__)
    input_kind: str = None              # pacbio_subreads|pacbio_ccs|reads
    isoquant_data_type: str = None      # pacbio_ccs|nanopore
    isoseq3_supported: str = None       # full|none
    needs_ccs: bool = False             # 是否跑ccs步骤
    needs_refine: bool = False          # 是否跑refine步骤
    ccs_dir: str = None
    refine_dir: str = None
    isoquant_dir: str = None
    isoseq3_dir: str = None
    stat_dir: str = None
    log_dir: str = None
    tmp_dir: str = None

    def __post_init__(self):
        # 展开~路径(§11.B)|Expand ~ paths
        for p in ("reference", "genedb", "primers", "output_dir"):
            if getattr(self, p):
                setattr(self, p, expand_path(getattr(self, p)))
        self.reads = [expand_path(r) for r in self.reads]
        for p in ("ccs_path", "isoseq3_path", "isoquant_path", "samtools_path"):
            setattr(self, p, expand_path(getattr(self, p)))
        if self.data_type:
            self.data_type = self.data_type.lower()
        self.engine = self.engine.lower()

        # 嗅探与派生|Sniff and derive
        self.input_kind = detect_input_type(self.reads)
        self.needs_ccs = self.input_kind == "pacbio_subreads"
        self.needs_refine = self.input_kind in ("pacbio_subreads", "pacbio_ccs")
        if self.input_kind == "reads":
            self.isoquant_data_type = "nanopore" if self.data_type == "ont" else "pacbio_ccs"
        else:
            self.isoquant_data_type = "pacbio_ccs"
        # isoseq3 26.2矩阵: 只要有BAM(subreads或ccs)即可 refine->cluster2 全链
        # |26.2 matrix: any PacBio BAM supports full refine->cluster2 chain
        self.isoseq3_supported = "full" if self.needs_refine else "none"

        # 输出目录(§12.2 by-step)|Output directories
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.ccs_dir = os.path.join(self.output_dir, "01_ccs")
        self.refine_dir = os.path.join(self.output_dir, "02_refine")
        self.isoquant_dir = os.path.join(self.output_dir, "03_isoquant")
        self.isoseq3_dir = os.path.join(self.output_dir, "04_isoseq3")
        self.stat_dir = os.path.join(self.output_dir, "00_pipeline_info")
        self.log_dir = os.path.join(self.output_dir, "99_logs")
        self.tmp_dir = os.path.join(self.output_dir, "tmp")

    @property
    def run_isoquant(self) -> bool:
        return self.engine in ("isoquant", "both")

    @property
    def run_isoseq3(self) -> bool:
        return self.engine in ("isoseq3", "both")

    def validate(self):
        """校验配置,一次性收集全部错误|Validate, collect all errors at once"""
        errors = []

        for r in self.reads:
            if not os.path.exists(r):
                errors.append(f"输入reads文件不存在|Reads file not found: {r}")

        if self.data_type not in (None, "pacbio", "ont"):
            errors.append(f"data-type须为pacbio或ont|data-type must be pacbio/ont: {self.data_type}")
        if self.input_kind == "reads" and self.data_type is None:
            errors.append("reads文件(fasta/fastq)必须显式指定--data-type pacbio或ont|"
                          "reads files (fasta/fastq) require explicit --data-type pacbio or ont")

        # isoseq3可用性|isoseq3 availability
        if self.run_isoseq3 and self.isoseq3_supported == "none":
            errors.append("isoseq3引擎仅支持PacBio BAM输入(subreads.bam/ccs.bam),"
                          "不支持fasta/fastq或ONT数据|isoseq3 engine requires PacBio BAM input, "
                          "not fasta/fastq or ONT data")

        if self.run_isoquant:
            if not self.reference:
                errors.append("IsoQuant引擎需要参考基因组(--reference)|IsoQuant engine requires --reference")
            elif not os.path.exists(self.reference):
                errors.append(f"参考基因组不存在|Reference not found: {self.reference}")

        if self.genedb and not os.path.exists(self.genedb):
            errors.append(f"注释文件不存在|Genedb not found: {self.genedb}")

        if self.primers:
            if not os.path.exists(self.primers):
                errors.append(f"引物文件不存在|Primers not found: {self.primers}")
            else:
                n_records = sum(1 for line in open(self.primers) if line.startswith(">"))
                if n_records != 2:
                    errors.append(f"primers fasta须恰好2条(5p/3p)|primers fasta needs exactly "
                                  f"2 records (5p/3p), got {n_records}")

        if self.engine not in ("isoquant", "isoseq3", "both"):
            errors.append(f"engine须为isoquant/isoseq3/both|engine must be isoquant/isoseq3/both: {self.engine}")
        if self.threads <= 0:
            errors.append(f"线程数须为正整数|Threads must be positive: {self.threads}")
        if self.min_passes < 1:
            errors.append(f"min-passes须>=1|min-passes must be >=1: {self.min_passes}")

        if errors:
            raise ValueError("\n".join(errors))
        return True

    def write_primers_fasta(self) -> str:
        """写出引物fasta到tmp(内置Clontech或复制用户文件)|Write primers fasta to tmp"""
        Path(self.tmp_dir).mkdir(parents=True, exist_ok=True)
        target = os.path.join(self.tmp_dir, "primers.fasta")
        if self.primers:
            shutil.copyfile(self.primers, target)
        else:
            with open(target, "w") as f:
                f.write(CLONTECH_PRIMERS_FA)
        return target
