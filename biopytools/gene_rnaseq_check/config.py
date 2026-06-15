"""候选基因RNA-seq转录验证配置模块|Candidate Gene RNA-seq Validation Configuration Module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path


@dataclass
class GeneRnaseqCheckConfig:
    """候选基因RNA-seq转录验证配置类|Candidate Gene RNA-seq Validation Config"""

    # 必需参数|Required parameters
    genome_fa: str = ''
    annotation_gff: str = ''
    gene_list: str = ''
    reads_dir: str = ''
    output_dir: str = './gene_rnaseq_check_output'

    # 可选参数|Optional parameters
    reads_pattern: str = ''
    threads: int = 12
    sample_timeout: int = 21600
    steps: str = 'all'
    force: bool = False
    verbose: bool = False
    quiet: bool = False

    # 分类阈值|Classification thresholds
    expressed_complete_cov_threshold: float = 90.0
    expressed_complete_exon_threshold: float = 80.0
    expressed_partial_low: float = 30.0
    expressed_partial_high: float = 90.0
    expressed_partial_exon_threshold: float = 50.0
    boundary_issue_min_cov: float = 60.0
    boundary_issue_flank_ratio: float = 0.3
    boundary_issue_flank_min_depth: float = 5.0
    not_expressed_max_cov: float = 10.0
    not_expressed_flank_max: float = 5.0
    strandness_confidence: float = 70.0
    junction_tolerance: int = 5
    flanking_window: int = 500

    # 工具路径|Tool paths
    hisat2_path: str = field(
        default_factory=lambda: get_tool_path(
            'hisat2', '~/miniforge3/envs/RNA_Seq/bin/hisat2', 'HISAT2_PATH'
        )
    )
    hisat2_build_path: str = field(
        default_factory=lambda: get_tool_path(
            'hisat2-build', '~/miniforge3/envs/RNA_Seq/bin/hisat2-build', 'HISAT2_BUILD_PATH'
        )
    )
    extract_splice_sites_path: str = field(
        default_factory=lambda: get_tool_path(
            'extract_splice_sites.py',
            '~/miniforge3/envs/RNA_Seq/bin/extract_splice_sites.py',
            'EXTRACT_SPLICE_SITES_PATH',
        )
    )
    extract_exons_path: str = field(
        default_factory=lambda: get_tool_path(
            'extract_exons.py',
            '~/miniforge3/envs/RNA_Seq/bin/extract_exons.py',
            'EXTRACT_EXONS_PATH',
        )
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'samtools', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools', 'SAMTOOLS_PATH'
        )
    )
    bedtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'bedtools', '~/miniforge3/envs/Population_genetics/bin/bedtools', 'BEDTOOLS_PATH'
        )
    )
    stringtie_path: str = field(
        default_factory=lambda: get_tool_path(
            'stringtie', '~/miniforge3/envs/RNA_Seq/bin/stringtie', 'STRINGTIE_PATH'
        )
    )
    gffcompare_path: str = field(
        default_factory=lambda: get_tool_path(
            'gffcompare', '~/miniforge3/envs/gffcompare_v.0.12.10/bin/gffcompare', 'GFFCOMPARE_PATH'
        )
    )
    infer_experiment_path: str = field(
        default_factory=lambda: get_tool_path(
            'infer_experiment.py',
            '~/miniforge3/envs/RSeQC_v.5.0.4/bin/infer_experiment.py',
            'INFER_EXPERIMENT_PATH',
        )
    )

    # 内部属性|Internal attributes
    _samples: list = field(default_factory=list, init=False)
    _gene_ids: set = field(default_factory=set, init=False)
    output_path: Path = field(default=None, init=False)
    log_dir: Path = field(default=None, init=False)
    strandness: str = field(default='', init=False)  # 'FR', 'RF', or 'unstranded'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.log_dir = self.output_path / '99_logs'
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 展开所有路径|Expand all paths
        path_fields = [
            'genome_fa', 'annotation_gff', 'gene_list', 'reads_dir',
            'hisat2_path', 'hisat2_build_path',
            'extract_splice_sites_path', 'extract_exons_path',
            'samtools_path', 'bedtools_path', 'stringtie_path',
            'gffcompare_path', 'infer_experiment_path',
        ]
        for f_name in path_fields:
            val = getattr(self, f_name, '')
            if val:
                setattr(self, f_name, expand_path(val))

        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.genome_fa = os.path.normpath(os.path.abspath(self.genome_fa))
        self.annotation_gff = os.path.normpath(os.path.abspath(self.annotation_gff))

        # 解析基因ID列表|Parse gene ID list
        self._parse_gene_list()

        # 检测样本|Detect samples
        self._detect_samples()

    def _parse_gene_list(self):
        """解析基因ID列表文件|Parse gene ID list file"""
        if not self.gene_list:
            return
        p = Path(expand_path(self.gene_list))
        if p.is_file():
            with open(p) as f:
                self._gene_ids = {
                    line.strip()
                    for line in f
                    if line.strip() and not line.startswith('#')
                }
            self.gene_list = os.path.normpath(os.path.abspath(str(p)))

    def _detect_samples(self):
        """检测RNA-seq样本|Detect RNA-seq samples"""
        import glob as glob_mod

        r_dir = self.reads_dir
        if not r_dir or not os.path.isdir(r_dir):
            return

        # 自定义模式|Custom pattern
        if self.reads_pattern and '*' in self.reads_pattern:
            samples = self._parse_with_pattern(r_dir, self.reads_pattern)
            if samples:
                self._samples = samples
                return

        # 默认模式|Default patterns
        default_patterns = [
            ("*_1.fq.gz", "*_2.fq.gz"),
            ("*_R1.fq.gz", "*_R2.fq.gz"),
            ("*.R1.fastq.gz", "*.R2.fastq.gz"),
            ("*_1.fastq.gz", "*_2.fastq.gz"),
            ("*_f1.fq.gz", "*_r2.fq.gz"),
        ]
        read_indicators = [
            ("R1", "R2"), ("_1", "_2"), (".1", ".2"),
            ("_f1", "_r2"), ("_F1", "_R2"),
        ]

        for p1, p2 in default_patterns:
            if '*' not in p1:
                continue
            prefix, suffix = p1.split("*", 1)
            r1_ind = r2_ind = None
            for r1, r2 in read_indicators:
                if r1 in suffix:
                    r1_ind, r2_ind = r1, r2
                    break
            if not r1_ind:
                continue

            fq1_files = sorted(glob_mod.glob(os.path.join(r_dir, p1)))
            samples = []
            for fq1 in fq1_files:
                basename = os.path.basename(fq1)
                name = basename
                if prefix:
                    name = name.replace(prefix, "", 1)
                if suffix:
                    name = name.replace(suffix, "", 1)
                r2_suffix = suffix.replace(r1_ind, r2_ind)
                fq2 = os.path.join(r_dir, prefix + name + r2_suffix)
                if os.path.exists(fq2):
                    samples.append({"name": name, "fastq1": fq1, "fastq2": fq2})

            if samples:
                self._samples = samples
                return

    def _parse_with_pattern(self, r_dir: str, pattern: str) -> list:
        """用自定义模式解析样本|Parse samples with custom pattern"""
        import glob as glob_mod

        parts = pattern.split("*")
        if len(parts) != 2:
            return []
        prefix, suffix = parts
        read_indicators = [
            ("R1", "R2"), ("_1", "_2"), (".1", ".2"),
            ("_f1", "_r2"), ("_F1", "_R2"),
        ]
        r1_ind = r2_ind = None
        for r1, r2 in read_indicators:
            if r1 in suffix:
                r1_ind, r2_ind = r1, r2
                break
        if not r1_ind:
            return []

        fq1_files = sorted(glob_mod.glob(os.path.join(r_dir, pattern)))
        samples = []
        for fq1 in fq1_files:
            basename = os.path.basename(fq1)
            name = basename
            if prefix:
                name = name.replace(prefix, "", 1)
            if suffix:
                name = name.replace(suffix, "", 1)
            r2_suffix = suffix.replace(r1_ind, r2_ind)
            fq2 = os.path.join(r_dir, prefix + name + r2_suffix)
            if os.path.exists(fq2):
                samples.append({"name": name, "fastq1": fq1, "fastq2": fq2})
        return samples

    def get_steps_list(self) -> List[str]:
        """解析步骤列表|Parse steps list"""
        valid_steps = [
            'align', 'parse_gff', 'coverage',
            'junction', 'stringtie', 'classify', 'boundary',
        ]
        if self.steps.lower() == 'all':
            return valid_steps[:]
        return [s.strip() for s in self.steps.split(',') if s.strip() in valid_steps]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not self.genome_fa or not os.path.exists(self.genome_fa):
            errors.append(f"基因组文件不存在|Genome FASTA not found: {self.genome_fa}")
        if not self.annotation_gff or not os.path.exists(self.annotation_gff):
            errors.append(f"注释文件不存在|Annotation GFF3 not found: {self.annotation_gff}")
        if not self.gene_list or not os.path.exists(self.gene_list):
            errors.append(f"基因ID列表文件不存在|Gene list file not found: {self.gene_list}")
        if not self._gene_ids:
            errors.append("基因ID列表为空|Gene ID list is empty")
        if not self.reads_dir or not os.path.isdir(self.reads_dir):
            errors.append(f"Reads目录不存在|Reads directory not found: {self.reads_dir}")
        if not self._samples:
            errors.append("未检测到RNA-seq样本|No RNA-seq samples detected")
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查工具路径|Check tool paths
        tool_checks = [
            ('HISAT2', self.hisat2_path),
            ('samtools', self.samtools_path),
            ('bedtools', self.bedtools_path),
            ('StringTie', self.stringtie_path),
            ('gffcompare', self.gffcompare_path),
            ('infer_experiment.py', self.infer_experiment_path),
        ]
        for name, path in tool_checks:
            if not os.path.exists(path):
                errors.append(f"{name} 不存在|{name} not found: {path}")

        if errors:
            raise ValueError("\n".join(errors))
        return True
