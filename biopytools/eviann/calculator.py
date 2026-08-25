"""
EviAnn输入数据处理与流程编排|EviAnn Input Processing & Pipeline Orchestration

核心职责|Core responsibilities:
- 转录组文件自动识别(二代/三代、配对/聚类、mix匹配)
- EviAnn -r 描述文件与样本清单生成
- 多run文件合并、输出目录编排、命令构建与执行
"""

import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from biopytools.common.conda_runner import CommandRunner, build_conda_command
from biopytools.common.paths import expand_path

# 三代长读文件名关键词|Long-read filename keywords
LONG_READ_KEYWORDS = (
    "isoseq", "iso_seq", "iso-seq", "iso",
    "ccs", "hifi", "pacbio", "nanopore", "ont", "dorado",
    "simplex", "duplex",
)

# 过滤标记(剥离样本名时与关键词同样处理)|Filter markers (stripped like keywords)
_FILTER_MARKERS = ("clean",)

# 双端配对正则|Paired-end pattern: <prefix>[._-](r?1|r?2|read1|read2)([._-]|$)
_PAIR_RE = re.compile(
    r"^(.*?)[._-](r[12]|read[12]|[12])(?:[._-]|$)", re.IGNORECASE)

# 转录组数据扩展名(允许 .gz)|Transcriptome data extensions (.gz allowed)
_READ_EXTS = (".fq", ".fastq", ".fa", ".fasta", ".fna", ".fas", ".bam")

# eviann 最终结果后缀|EviAnn final result suffixes
RESULT_SUFFIXES = (".pseudo_label.gff", ".proteins.fasta", ".transcripts.fasta")


def _strip_extensions(name: str) -> str:
    """去掉测序扩展名(.fq.gz/.fasta/.bam等)|Strip sequencing extensions"""
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
        lower = lower[:-3]
    for ext in (".fastq", ".fasta", ".fq", ".fa", ".fna", ".fas", ".bam"):
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def _keyword_re(keyword: str) -> re.Pattern:
    """词边界关键词正则(允许尾随数字编号)|Word-boundary keyword pattern
    (trailing digits allowed, e.g. S1_ccs1/S1_ccs2)"""
    return re.compile(
        r"(?:^|[._-])" + re.escape(keyword) + r"(?:\d+)?(?:[._-]|$)",
        re.IGNORECASE)


def _contains_long_keyword(base: str) -> bool:
    """文件名含三代关键词(词边界)|Contains long-read keyword (word boundary)"""
    return any(_keyword_re(k).search(base) for k in LONG_READ_KEYWORDS)


def _contains_marker(base: str, marker: str) -> bool:
    """文件名含标记(词边界)|Contains marker (word boundary)"""
    return bool(_keyword_re(marker).search(base))


def classify_file(path: str) -> str:
    """识别转录组文件类型|Classify transcriptome file type

    规则|Rules:
    - fasta 一律视为三代|fasta is always long-read
    - bam 按关键词区分 bam/bam_isoseq|bam split by keyword
    - fastq 含三代关键词 → 三代|fastq with keyword → long
    - fastq 无配对痕迹且含 .clean → 三代(用户默认三代后缀 *.clean.fq.gz)
      |fastq without pair trace but with .clean → long (default long suffix)
    - 其余 fastq → 二代|otherwise → short

    Args:
        path: 文件路径|File path

    Returns:
        'short_fastq'|'long_fastq'|'long_fasta'|'short_bam'|'long_bam'|'unknown'
    """
    lower = os.path.basename(path).lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    if lower.endswith(".bam"):
        base = lower[:-4]
        return "long_bam" if _contains_long_keyword(base) else "short_bam"
    for ext in (".fasta", ".fastq", ".fna", ".fas", ".fa", ".fq"):
        if lower.endswith(ext):
            base = lower[:-len(ext)]
            if ext in (".fa", ".fasta", ".fna", ".fas"):
                return "long_fasta"
            if _contains_long_keyword(base):
                return "long_fastq"
            if (extract_pair_prefix(base) is None
                    and _contains_marker(base, "clean")):
                return "long_fastq"
            return "short_fastq"
    return "unknown"


def extract_pair_prefix(name: str) -> Optional[Tuple[str, str]]:
    """提取双端配对前缀和方向|Extract paired-end sample prefix and mate

    支持 _R1/_R2、_1/_2、.R1.、read1/read2、_R1_001、_1.clean 等形态
    |Supports _R1/_R2, _1/_2, .R1., read1/read2, _R1_001, _1.clean etc.

    Args:
        name: 纯文件名|File name only

    Returns:
        (样本前缀, '1'|'2') 或 None|(sample prefix, mate) or None
    """
    base = _strip_extensions(name)
    m = _PAIR_RE.search(base)
    if not m or not m.group(1):
        return None
    direction = "1" if m.group(2)[-1] == "1" else "2"
    return (m.group(1), direction)


def strip_long_keywords(name: str) -> str:
    """剥离三代关键词与过滤标记得到样本名|Strip keywords/markers to get sample name

    Args:
        name: 纯文件名|File name only

    Returns:
        样本名|Sample name
    """
    base = _strip_extensions(name)
    original = base
    for kw in sorted(LONG_READ_KEYWORDS + _FILTER_MARKERS, key=len,
                     reverse=True):
        base = _keyword_re(kw).sub("_", base)
    base = re.sub(r"[._-]{2,}", "_", base).strip("._-")
    return base or original


def pair_short_reads(
        paths: List[str]) -> Tuple[List[Tuple[str, List[str], List[str]]], List[str]]:
    """双端fastq分组|Group paired-end fastq by sample

    同一样本的多个run归入同组(合并前形态)|Multiple runs of one sample grouped

    Args:
        paths: 二代fastq路径列表|Short-read fastq paths

    Returns:
        (配对组, 单端文件)|(paired groups, single-end files)
        配对组: (样本名, r1列表, r2列表)|paired: (sample, r1 list, r2 list)
    """
    groups: Dict[str, Dict[str, List[str]]] = {}
    singles: List[str] = []
    for p in paths:
        parsed = extract_pair_prefix(os.path.basename(p))
        if parsed is None:
            singles.append(p)
            continue
        sample, direction = parsed
        groups.setdefault(sample, {"1": [], "2": []})[direction].append(p)
    paired = [(s, d["1"], d["2"]) for s, d in groups.items()
              if d["1"] and d["2"]]
    for s, d in groups.items():
        if not (d["1"] and d["2"]):
            # 缺一侧配对 → 归单端|one-sided pair → single
            singles.extend(d["1"] + d["2"])
    return paired, singles


def cluster_long_reads(paths: List[str]) -> List[Tuple[str, List[str]]]:
    """长读文件按样本聚类|Cluster long-read files by sample

    Args:
        paths: 长读路径列表|Long-read paths

    Returns:
        [(样本名, 文件列表)] 保持首次出现顺序|[(sample, files)] in first-seen order
    """
    groups: Dict[str, List[str]] = {}
    for p in paths:
        sample = strip_long_keywords(os.path.basename(p))
        groups.setdefault(sample, []).append(p)
    return list(groups.items())


def match_mix(short_names: List[str],
              long_names: List[str]) -> Dict[str, Optional[str]]:
    """为短读样本匹配长读样本|Match long-read samples to short-read samples

    规则|Rules:
    - 样本名完全相等或互为前缀 → 候选|exact or prefix relation → candidate
    - 冲突时短名更长者优先,取最具体的候选|longest short name wins, most specific candidate

    Args:
        short_names: 二代样本名列表|Short-read sample names
        long_names: 三代样本名列表|Long-read sample names

    Returns:
        {短样本名: 匹配的长样本名或None}|{short name: matched long name or None}
    """
    result: Dict[str, Optional[str]] = {}
    available = set(long_names)
    for s in sorted(set(short_names), key=len, reverse=True):
        candidates = [l for l in available
                      if s == l or s.startswith(l) or l.startswith(s)]
        if candidates:
            best = max(candidates, key=len)
            result[s] = best
            available.remove(best)
        else:
            result[s] = None
    return result


@dataclass
class Experiment:
    """一个EviAnn实验(描述文件一行)|One EviAnn experiment (one -r line)"""

    sample: str
    tag: str = "fastq"  # fastq|isoseq|mix|bam|bam_isoseq|bam_mix
    # 合并后的单文件路径|Merged single-file paths
    r1: Optional[str] = None
    r2: Optional[str] = None
    long_reads: Optional[str] = None
    short_bam: Optional[str] = None
    long_bam: Optional[str] = None
    # 合并前的多文件列表|Pre-merge file lists
    r1_files: List[str] = field(default_factory=list)
    r2_files: List[str] = field(default_factory=list)
    long_files: List[str] = field(default_factory=list)
    short_bam_files: List[str] = field(default_factory=list)
    long_bam_files: List[str] = field(default_factory=list)


def format_rnaseq_line(exp: Experiment) -> str:
    """格式化为EviAnn -r 描述文件行|Format one EviAnn -r line"""
    if exp.tag == "fastq":
        return f"{exp.r1} {exp.r2} fastq" if exp.r2 else f"{exp.r1} fastq"
    if exp.tag == "isoseq":
        return f"{exp.long_reads} isoseq"
    if exp.tag == "mix":
        return f"{exp.r1} {exp.r2} {exp.long_reads} mix"
    if exp.tag == "bam":
        return f"{exp.short_bam} bam"
    if exp.tag == "bam_isoseq":
        return f"{exp.long_bam} bam_isoseq"
    if exp.tag == "bam_mix":
        return f"{exp.short_bam} {exp.long_bam} bam_mix"
    raise ValueError(f"未知标签|Unknown tag: {exp.tag}")


def write_rnaseq_list(experiments: List[Experiment], out_path: str) -> None:
    """写出EviAnn -r 描述文件|Write EviAnn -r description file

    注意:不写注释行(EviAnn按行解析)|Note: no comment lines (EviAnn parses lines)
    """
    if not experiments:
        raise ValueError("无实验可写|No experiments to write")
    lines = [format_rnaseq_line(e) for e in experiments]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _split_files(cell: str) -> List[str]:
    """逗号分隔文件列表|Comma-separated file list"""
    return [p.strip() for p in cell.split(",") if p.strip()]


def _build_experiment(sample: str, r1s: List[str], r2s: List[str],
                      longs: List[str], tag_col: str) -> Experiment:
    """按列内容推断实验|Infer experiment from columns"""
    def is_bam(p):
        return classify_file(p) in ("short_bam", "long_bam")

    if tag_col:
        exp = Experiment(sample=sample, tag=tag_col)
        if tag_col == "fastq":
            exp.r1_files, exp.r2_files = r1s, r2s
        elif tag_col == "isoseq":
            exp.long_files = longs
        elif tag_col == "mix":
            exp.r1_files, exp.r2_files, exp.long_files = r1s, r2s, longs
        elif tag_col == "bam":
            exp.short_bam_files = r1s
        elif tag_col == "bam_isoseq":
            exp.long_bam_files = longs
        elif tag_col == "bam_mix":
            exp.short_bam_files, exp.long_bam_files = r1s, longs
        else:
            raise ValueError(f"未知标签|Unknown tag: {tag_col}")
        return exp

    r1 = r1s[0] if r1s else None
    r2 = r2s[0] if r2s else None
    long = longs[0] if longs else None
    if r1 and is_bam(r1):
        if long and is_bam(long):
            return Experiment(sample, "bam_mix",
                              short_bam_files=r1s, long_bam_files=longs)
        return Experiment(sample, "bam", short_bam_files=r1s)
    if long and is_bam(long):
        return Experiment(sample, "bam_isoseq", long_bam_files=longs)
    if r1 and long:
        return Experiment(sample, "mix",
                          r1_files=r1s, r2_files=r2s, long_files=longs)
    if r1 and r2:
        return Experiment(sample, "fastq", r1_files=r1s, r2_files=r2s)
    if r1:
        return Experiment(sample, "fastq", r1_files=r1s)
    if long:
        return Experiment(sample, "isoseq", long_files=longs)
    raise ValueError(f"样本行无有效文件|No valid files for sample: {sample}")


def parse_sample_sheet(path: str, check_files: bool = True) -> List[Experiment]:
    """解析样本清单TSV|Parse sample sheet TSV

    格式|Format: sample_id\\tr1\\tr2\\tlong_reads[\\ttag]
    - 列可为空;每格逗号分隔多个文件|columns optional; comma-separated files per cell
    - # 开头为注释|# starts a comment

    Args:
        path: 清单路径|Sheet path
        check_files: 校验文件存在|Verify file existence

    Returns:
        实验列表|Experiment list
    """
    experiments: List[Experiment] = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                raise ValueError(
                    f"清单行格式错误|Bad sheet line {lineno}: {line}")
            sample = cols[0].strip()
            if not sample:
                raise ValueError(f"清单样本名为空|Empty sample name at line {lineno}")
            r1s = _split_files(cols[1]) if len(cols) > 1 else []
            r2s = _split_files(cols[2]) if len(cols) > 2 else []
            longs = _split_files(cols[3]) if len(cols) > 3 else []
            tag_col = cols[4].strip() if len(cols) > 4 else ""
            if check_files:
                for fpath in r1s + r2s + longs:
                    if not os.path.exists(expand_path(fpath)):
                        raise ValueError(
                            f"清单文件不存在|Sheet file not found: {fpath} "
                            f"(line {lineno})")
            experiments.append(
                _build_experiment(sample, r1s, r2s, longs, tag_col))
    if not experiments:
        raise ValueError(f"清单无有效行|No valid rows in sheet: {path}")
    return experiments


def write_sample_sheet_tsv(experiments: List[Experiment],
                           out_path: str) -> None:
    """写出样本清单(供人工查看/修改)|Write sample sheet for manual review"""
    def cell(files: List[str], single: Optional[str]) -> str:
        return ",".join(files) if files else (single or "")

    header = "# sample_id\tr1\tr2\tlong_reads\ttag"
    lines = [header]
    for e in experiments:
        r1 = cell(e.r1_files, e.r1) or cell(e.short_bam_files, e.short_bam)
        r2 = cell(e.r2_files, e.r2)
        long = cell(e.long_files, e.long_reads) or cell(
            e.long_bam_files, e.long_bam)
        lines.append(f"{e.sample}\t{r1}\t{r2}\t{long}\t{e.tag}")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def merge_files(paths: List[str], dest_dir: str, sample: str,
                kind: str) -> Optional[str]:
    """多run文件合并(二进制拼接)|Merge multi-run files (binary concat)

    单文件透传;目标已存在则跳过(断点续传友好)
    |Single file passed through; existing target skipped (resume friendly)

    Args:
        paths: 文件路径列表|File paths
        dest_dir: 合并输出目录|Merge destination dir
        sample: 样本名|Sample name
        kind: 类别(r1/r2/long/shortbam/longbam)|Kind

    Returns:
        最终单文件路径|Final single-file path
    """
    if not paths:
        return None
    os.makedirs(dest_dir, exist_ok=True)
    ext = "".join(Path(paths[0]).suffixes)
    merged = os.path.join(dest_dir, f"{sample}.{kind}.merged{ext}")
    if os.path.exists(merged):
        return merged
    if len(paths) == 1:
        return paths[0]
    with open(merged, "wb") as out:
        for p in paths:
            with open(p, "rb") as fh:
                shutil.copyfileobj(fh, out)
    return merged


def _merge_experiment_files(exp: Experiment, merged_dir: str) -> Experiment:
    """填充实验的单文件字段|Populate experiment single-file fields"""
    if exp.r1_files:
        exp.r1 = merge_files(exp.r1_files, merged_dir, exp.sample, "r1")
    if exp.r2_files:
        exp.r2 = merge_files(exp.r2_files, merged_dir, exp.sample, "r2")
    if exp.long_files:
        exp.long_reads = merge_files(exp.long_files, merged_dir, exp.sample,
                                     "long")
    if exp.short_bam_files:
        exp.short_bam = merge_files(exp.short_bam_files, merged_dir,
                                    exp.sample, "shortbam")
    if exp.long_bam_files:
        exp.long_bam = merge_files(exp.long_bam_files, merged_dir, exp.sample,
                                   "longbam")
    return exp


class EviAnnRunner:
    """EviAnn注释流程编排器|EviAnn Annotation Pipeline Orchestrator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        out = config.output_dir
        self.work_dir = os.path.join(out, "work")
        self.inputs_dir = os.path.join(out, "01_inputs")
        self.results_dir = os.path.join(out, "results")
        self.info_dir = os.path.join(out, "00_pipeline_info")
        self.logs_dir = os.path.join(out, "99_logs")
        self.merged_dir = os.path.join(self.work_dir, "merged")

    # ---------- 工作区|Workspace ----------

    def prepare_workspace(self) -> None:
        """创建目录结构并软链基因组|Create layout, symlink genome"""
        for d in (self.info_dir, self.inputs_dir, self.work_dir,
                  self.results_dir, self.logs_dir):
            os.makedirs(d, exist_ok=True)
        # 基因组软链(保持原文件名,EviAnn以文件名做输出前缀)
        # |genome symlink (keep name; EviAnn uses file name as output prefix)
        genome_name = os.path.basename(self.config.genome)
        link = os.path.join(self.work_dir, genome_name)
        if not os.path.exists(link):
            os.symlink(self.config.genome, link)
            self.logger.info(f"基因组软链|Genome symlinked: {link}")
        elif not os.path.islink(link) or os.path.realpath(
                link) != os.path.realpath(self.config.genome):
            self.logger.warning(
                f"work目录已存在同名文件且指向不同路径|Existing file in work "
                f"points elsewhere, keeping it for resume: {link}")

    # ---------- 输入处理|Input processing ----------

    @staticmethod
    def _is_read_file(name: str) -> bool:
        """按原始文件名判断是否测序数据|Check read extensions on raw name"""
        lower = name.lower()
        if lower.endswith(".gz"):
            lower = lower[:-3]
        return any(lower.endswith(e) for e in _READ_EXTS)

    def _collect_input_files(self) -> List[str]:
        """展开 --rnaseq-data 条目为文件列表|Expand --rnaseq-data entries"""
        collected: List[str] = []
        for entry in self.config.rnaseq_data.split(","):
            entry = entry.strip()
            if not entry:
                continue
            if os.path.isdir(entry):
                for name in sorted(os.listdir(entry)):
                    p = os.path.join(entry, name)
                    if os.path.isfile(p) and self._is_read_file(name):
                        collected.append(p)
            elif os.path.isfile(entry):
                collected.append(entry)
            else:
                raise ValueError(f"输入不存在|Input not found: {entry}")
        seen: set = set()
        uniq: List[str] = []
        for p in collected:
            if p not in seen:
                seen.add(p)
                uniq.append(p)
        if not uniq:
            raise ValueError("未找到转录组数据文件|No transcriptome files found")
        return uniq

    def _build_experiments_from_files(self,
                                      files: List[str]) -> List[Experiment]:
        """自动识别文件并构建实验|Auto-classify files into experiments"""
        from collections import defaultdict
        by_type = defaultdict(list)
        for f in files:
            by_type[classify_file(f)].append(f)
        for ftype, fs in by_type.items():
            self.logger.info(
                f"识别|Classified {ftype}: {len(fs)} 个文件|files")

        experiments: List[Experiment] = []

        # 二代配对 + 三代聚类 + mix 匹配|pair shorts, cluster longs, match mix
        paired, singles = pair_short_reads(by_type.get("short_fastq", []))
        long_groups = cluster_long_reads(
            by_type.get("long_fastq", []) + by_type.get("long_fasta", []))
        mix_map = match_mix([s for s, _, _ in paired],
                            [s for s, _ in long_groups])
        long_by_name = dict(long_groups)
        used_longs: set = set()

        for sample, r1s, r2s in paired:
            matched = mix_map.get(sample)
            if matched:
                experiments.append(
                    Experiment(sample=sample, tag="mix",
                                r1_files=r1s, r2_files=r2s,
                                long_files=long_by_name[matched]))
                used_longs.add(matched)
                self.logger.info(
                    f"混合配对|Mix matched: {sample} + {matched}")
            else:
                experiments.append(
                    Experiment(sample=sample, tag="fastq",
                                r1_files=r1s, r2_files=r2s))

        for f in singles:
            # 缺配对一侧 → 单端二代|missing mate → single-end short
            sample = strip_long_keywords(os.path.basename(f))
            experiments.append(
                Experiment(sample=sample, tag="fastq", r1_files=[f]))

        for name, fs in long_groups:
            if name not in used_longs:
                experiments.append(
                    Experiment(sample=name, tag="isoseq", long_files=fs))

        # bam 类|bam handling
        short_bams = by_type.get("short_bam", [])
        long_bams = by_type.get("long_bam", [])
        long_bam_by_name: Dict[str, str] = {}
        for b in long_bams:
            long_bam_by_name.setdefault(
                strip_long_keywords(os.path.basename(b)), b)
        short_bam_names = [
            _strip_extensions(os.path.basename(b)) for b in short_bams]
        bam_map = match_mix(short_bam_names, list(long_bam_by_name))
        used_long_bams: set = set()
        for bam_path, name in zip(short_bams, short_bam_names):
            matched = bam_map.get(name)
            if matched:
                experiments.append(
                    Experiment(sample=name, tag="bam_mix",
                                short_bam_files=[bam_path],
                                long_bam_files=[long_bam_by_name[matched]]))
                used_long_bams.add(matched)
            else:
                experiments.append(
                    Experiment(sample=name, tag="bam",
                                short_bam_files=[bam_path]))
        for name, bam_path in long_bam_by_name.items():
            if name not in used_long_bams:
                experiments.append(
                    Experiment(sample=name, tag="bam_isoseq",
                                long_bam_files=[bam_path]))

        if not experiments:
            raise ValueError("未识别出任何转录组实验|No experiments identified")
        return experiments

    def prepare_experiments(self) -> List[Experiment]:
        """按输入模式准备实验列表|Prepare experiments per input mode"""
        if self.config.sample_sheet:
            self.logger.info(
                f"使用样本清单|Using sample sheet: {self.config.sample_sheet}")
            return parse_sample_sheet(self.config.sample_sheet)
        # 自动模式|auto mode
        self.logger.info(
            f"自动识别转录组数据|Auto-classifying RNA-seq data: "
            f"{self.config.rnaseq_data}")
        experiments = self._build_experiments_from_files(
            self._collect_input_files())
        # 写出清单副本供人工查看修改|write sheet copy for manual review
        sheet_path = os.path.join(self.inputs_dir, "sample_sheet.tsv")
        write_sample_sheet_tsv(experiments, sheet_path)
        self.logger.info(
            f"样本清单已生成(可修改后以--sample-sheet重跑)|Sample sheet "
            f"written (edit & rerun with --sample-sheet): {sheet_path}")
        return experiments

    def prepare_rnaseq_list(self) -> str:
        """准备EviAnn -r 描述文件|Prepare EviAnn -r list

        Returns:
            描述文件路径|Description file path
        """
        if self.config.rnaseq:
            self.logger.info(
                f"透传原生描述文件|Passthrough -r file: {self.config.rnaseq}")
            return self.config.rnaseq
        experiments = self.prepare_experiments()
        for exp in experiments:
            _merge_experiment_files(exp, self.merged_dir)
        list_path = os.path.join(self.inputs_dir, "rnaseq_list.txt")
        write_rnaseq_list(experiments, list_path)
        for line in open(list_path):
            self.logger.info(f"  实验行|Experiment: {line.rstrip()}")
        return list_path

    # ---------- 命令构建|Command building ----------

    def build_command(self, rnaseq_list: str) -> List[str]:
        """构建EviAnn命令|Build EviAnn command

        Args:
            rnaseq_list: -r 描述文件路径|Description file path

        Returns:
            完整命令列表(conda run 包装)|Full command list (conda run wrapped)
        """
        cfg = self.config
        eviann_sh = os.path.join(cfg.eviann_path, "bin", "eviann.sh")
        # 传纯文件名,cwd=work(EviAnn以文件名做输出前缀,中间文件落work)
        # |bare genome name, cwd=work (prefix from file name, intermediates in work)
        genome_name = os.path.basename(cfg.genome)
        args: List[str] = []
        if cfg.threads != 12:
            args += ["-t", str(cfg.threads)]
        args += ["-g", genome_name]
        args += ["-r", rnaseq_list]
        if cfg.transcripts:
            args += ["-e", cfg.transcripts]
        if cfg.proteins:
            args += ["-p", cfg.proteins]
        if cfg.uniprot:
            args += ["-s", cfg.uniprot]
        if cfg.max_intron:
            args += ["-m", str(cfg.max_intron)]
        if cfg.partial:
            args.append("--partial")
        if cfg.ploidy != 2:
            args += ["-d", str(cfg.ploidy)]
        if cfg.cds_gff:
            args += ["-c", cfg.cds_gff]
        if cfg.lncrna_tpm != 1.0:
            args += ["--lncrnamintpm", str(cfg.lncrna_tpm)]
        if cfg.min_prot:
            args += ["--min_prot", str(cfg.min_prot)]
        if cfg.functional:
            args.append("--functional")
        if cfg.mito_contigs:
            args += ["--mito_contigs", cfg.mito_contigs]
        if cfg.extra_gff:
            args += ["--extra", cfg.extra_gff]
        if cfg.debug:
            args.append("--debug")
        if cfg.verbose:
            args.append("--verbose")
        return build_conda_command(eviann_sh, args)

    # ---------- 结果与版本|Results & versions ----------

    def results_complete(self) -> bool:
        """最终结果齐全|All final results present"""
        genome_name = os.path.basename(self.config.genome)
        return all(
            os.path.exists(
                os.path.join(self.results_dir, genome_name + suffix))
            for suffix in RESULT_SUFFIXES)

    def collect_results(self) -> None:
        """复制最终结果到results/(work保留续传状态)|Copy results (keep work for resume)"""
        genome_name = os.path.basename(self.config.genome)
        for suffix in RESULT_SUFFIXES:
            src = os.path.join(self.work_dir, genome_name + suffix)
            dst = os.path.join(self.results_dir, genome_name + suffix)
            if os.path.exists(src):
                shutil.copy2(src, dst)
                self.logger.info(f"结果已收集|Collected: {genome_name}{suffix}")
            else:
                self.logger.warning(
                    f"结果缺失|Missing result: {genome_name}{suffix}")

    def write_versions(self) -> None:
        """写software_versions.yml|Write software_versions.yml"""
        eviann_sh = os.path.join(self.config.eviann_path, "bin", "eviann.sh")
        version = "unknown"
        try:
            runner = CommandRunner(self.logger)
            ok, stdout, _ = runner.run(
                build_conda_command(eviann_sh, ["--version"]),
                "探测EviAnn版本|Detecting EviAnn version")
            if ok:
                m = re.search(r"version\s+([\d.]+)", stdout)
                if m:
                    version = m.group(1)
        except Exception as e:
            self.logger.warning(f"版本探测失败|Version detection failed: {e}")
        from biopytools.eviann import __version__
        content = (
            f"module: biopytools.eviann {__version__}\n"
            f"eviann: {version}\n"
        )
        out = os.path.join(self.info_dir, "software_versions.yml")
        with open(out, "w") as f:
            f.write(content)
        self.logger.info(f"版本信息|Versions written: {out}")

    # ---------- 主流程|Main flow ----------

    def run(self) -> bool:
        """运行完整流程|Run full pipeline

        Returns:
            成功|Success
        """
        try:
            self.prepare_workspace()
            rnaseq_list = self.prepare_rnaseq_list()
            if self.results_complete():
                self.logger.info(
                    "最终结果已存在,跳过EviAnn执行|Results complete, "
                    "skipping EviAnn run")
                self.write_versions()
                return True
            cmd = self.build_command(rnaseq_list)
            runner = CommandRunner(self.logger, working_dir=self.work_dir)
            ok, _, stderr = runner.run(
                cmd, "EviAnn基因组注释|EviAnn genome annotation",
                capture_output=False)
            if not ok:
                self.logger.error(
                    f"EviAnn执行失败|EviAnn failed: "
                    f"{stderr[-500:] if stderr else ''}")
                return False
            self.collect_results()
            self.write_versions()
            self.logger.info(
                f"注释完成|Annotation complete. 结果目录|Results: "
                f"{self.results_dir}")
            return True
        except Exception as e:
            self.logger.error(f"流程出错|Pipeline error: {e}")
            return False
