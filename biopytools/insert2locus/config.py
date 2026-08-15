"""insert2locus配置与样本发现|insert2locus config and sample discovery"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from ..common.paths import expand_path, get_tool_path

# 配对后缀候选(检测顺序即优先级)|Paired suffix candidates (order = priority)
SUFFIX_CANDIDATES = [
    ("_1.clean.fq.gz", "_2.clean.fq.gz"),
    ("_1.fq.gz", "_2.fq.gz"),
]


def detect_suffixes(directory: Path) -> Optional[Tuple[str, str]]:
    """自动检测目录内配对后缀|Auto-detect paired suffixes in directory"""
    for r1, r2 in SUFFIX_CANDIDATES:
        if list(directory.glob(f"*{r1}")):
            return r1, r2
    return None


def extract_sample_name(filename: str, read1_suffix: str) -> str:
    """文件名去后缀得样本名|Strip suffix from filename to get sample name"""
    if filename.endswith(read1_suffix):
        return filename[: -len(read1_suffix)]
    return Path(filename).stem


@dataclass
class Insert2locusConfig:
    """insert2locus配置类|insert2locus configuration class"""

    # 必需|Required
    input_path: str            # fastq目录或R1文件|fastq dir or R1 file
    insert_fasta: str          # 完整构建(载体+插入片段)|Full construct (vector+insert)
    output_dir: str
    # 单独插入序列(可选;不给则-f整体当insert)|Standalone insert seq (optional)
    tdna_fasta: Optional[str] = None

    # 样本与资源|Sample & resources
    threads: int = 12
    sort_mem: str = "2G"       # samtools sort单线程内存|samtools sort per-thread memory
    read1_suffix: Optional[str] = None

    # 步移参数(与run.sh默认一致)|Walking params (same defaults as run.sh)
    max_rounds: int = 30
    min_softclip: int = 25
    min_unmapped: int = 400
    min_growth: int = 50
    mapq_min: int = 1
    repeat_cap: int = 10000
    junction_flank: int = 50
    # LB/RB目标侧翼长度;None=尽可能走远,靠自然收敛刹住|
    # Target LB/RB flank length; None = walk as far as naturally converges
    target_flank: Optional[int] = None

    # 工具路径(优先级:env var > config.yml > 默认~展开)|
    # Tool paths (priority: env var > config.yml > ~-expanded default)
    bwa_path: str = field(default_factory=lambda: get_tool_path(
        'bwa', '~/miniforge3/envs/Population_genetics/bin/bwa', 'BWA_PATH'))
    samtools_path: str = field(default_factory=lambda: get_tool_path(
        'samtools', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools', 'SAMTOOLS_PATH'))
    seqkit_path: str = field(default_factory=lambda: get_tool_path(
        'seqkit', '~/miniforge3/envs/BioinfTools/bin/seqkit', 'SEQKIT_PATH'))
    spades_path: str = field(default_factory=lambda: get_tool_path(
        'spades', '~/miniforge3/envs/spades_v.4.3.0/bin/spades.py', 'SPADES_PATH'))

    # 执行控制|Execution control
    force: bool = False
    log_level: str = "INFO"

    def __post_init__(self):
        # ⚠️ 所有可能含~的路径必须展开|All ~ paths must be expanded
        self.input_path = expand_path(self.input_path)
        self.insert_fasta = expand_path(self.insert_fasta)
        self.output_dir = expand_path(self.output_dir)
        if self.tdna_fasta:
            self.tdna_fasta = expand_path(self.tdna_fasta)
        for tool in ("bwa_path", "samtools_path", "seqkit_path", "spades_path"):
            setattr(self, tool, expand_path(getattr(self, tool)))
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self._input = Path(self.input_path)
        self.read2_suffix: Optional[str] = None
        if self.read1_suffix is None and self._input.is_dir():
            detected = detect_suffixes(self._input)
            if detected:
                self.read1_suffix, self.read2_suffix = detected
            else:
                self.read1_suffix, self.read2_suffix = SUFFIX_CANDIDATES[1]
        elif self.read1_suffix is None:
            # 单文件模式:从文件名推断;文件不存在时回落默认后缀,由validate统一报错|
            # Single-file mode: infer from filename; fall back to default when the
            # file does not exist so validate() reports it in one place
            for r1, r2 in SUFFIX_CANDIDATES:
                if self._input.name.endswith(r1):
                    self.read1_suffix, self.read2_suffix = r1, r2
                    break
            else:
                if self._input.exists():
                    raise ValueError(
                        f"无法识别输入文件配对后缀|Cannot detect pair suffix: "
                        f"{self._input.name}")
                self.read1_suffix, self.read2_suffix = SUFFIX_CANDIDATES[1]
        else:
            # 用户显式指定read1后缀,read2按命名约定推断|
            # User-specified read1 suffix, infer read2 by naming convention
            for r1, r2 in SUFFIX_CANDIDATES:
                if r1 == self.read1_suffix:
                    self.read2_suffix = r2
                    break
            if self.read2_suffix is None:
                self.read2_suffix = self.read1_suffix.replace("_1", "_2", 1)

    def discover_samples(self) -> List[Tuple[str, Path, Path]]:
        """发现样本(R1,R2)三元组|Discover (name, R1, R2) triples"""
        errors = []
        if self._input.is_file():
            sample = extract_sample_name(self._input.name, self.read1_suffix)
            r2 = self._input.parent / f"{sample}{self.read2_suffix}"
            if not r2.exists():
                raise ValueError(f"缺少配对文件|Missing mate: {r2}")
            return [(sample, self._input, r2)]
        pairs = []
        for r1_file in sorted(self._input.glob(f"*{self.read1_suffix}")):
            sample = extract_sample_name(r1_file.name, self.read1_suffix)
            r2_file = self._input / f"{sample}{self.read2_suffix}"
            if not r2_file.exists():
                errors.append(
                    f"样本{sample}缺少配对文件|Missing mate for {sample}: {r2_file}")
                continue
            pairs.append((sample, r1_file, r2_file))
        if errors:
            raise ValueError("\n".join(errors))
        if not pairs:
            raise ValueError(
                f"目录中未发现样本|No samples found with suffix "
                f"{self.read1_suffix}: {self._input}")
        return pairs

    def validate(self) -> None:
        """收集全部错误一次性抛出|Collect all errors and raise at once"""
        errors = []
        if not Path(self.input_path).exists():
            errors.append(f"输入不存在|Input not found: {self.input_path}")
        if not Path(self.insert_fasta).exists():
            errors.append(f"插入序列不存在|Insert fasta not found: {self.insert_fasta}")
        if self.tdna_fasta and not Path(self.tdna_fasta).exists():
            errors.append(
                f"tdna序列不存在|tdna fasta not found: {self.tdna_fasta}")
        if self.threads <= 0:
            errors.append("线程数必须为正|Threads must be positive")
        if self.max_rounds < 1:
            errors.append("最大轮数必须>=1|Max rounds must be >= 1")
        if errors:
            raise ValueError("\n".join(errors))
