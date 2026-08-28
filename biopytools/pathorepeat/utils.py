"""pathorepeat工具函数模块|pathorepeat Utilities Module

日志管理器(三文件分离) + 输入收集(单文件/目录) + 样品命名
|Logger (3-file separation) + input collection (single/dir) + sample naming
"""

import logging
import os
import re
import sys
from typing import Dict, List, NamedTuple, Optional

from .config import FASTA_GZ_SUFFIXES, FASTA_SUFFIXES


class PathorepeatLogger:
    """pathorepeat日志管理器|pathorepeat Logger Manager

    stdout(<=INFO) + stderr(>=WARNING) + 三文件:
    pathorepeat.log(全量) pathorepeat.out.log(<=INFO) pathorepeat.err.log(>=WARNING)
    """

    def __init__(self, logs_dir: str, log_level: str = 'INFO'):
        self.log_file = os.path.join(logs_dir, 'pathorepeat.log')
        self.out_log_file = os.path.join(logs_dir, 'pathorepeat.out.log')
        self.err_log_file = os.path.join(logs_dir, 'pathorepeat.err.log')
        os.makedirs(logs_dir, exist_ok=True)
        self.logger = self._setup_logging(log_level)

    def _setup_logging(self, log_level: str) -> logging.Logger:
        """设置日志(named logger,不污染root)|Setup logging (named, no root pollution)"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger('biopytools.pathorepeat')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        specs = [(self.log_file, None),
                 (self.out_log_file, lambda r: r.levelno <= logging.INFO),
                 (self.err_log_file, lambda r: r.levelno >= logging.WARNING)]
        for path, level_filter in specs:
            handler = logging.FileHandler(path, encoding='utf-8')
            handler.setLevel(logging.DEBUG)
            if level_filter:
                handler.addFilter(level_filter)
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


def genome_name(path: str) -> str:
    """样品名=basename剥全部后缀层|Sample name = basename, all suffix layers stripped"""
    name = os.path.basename(path)
    while True:
        root, ext = os.path.splitext(name)
        if ext.lower() in ('.gz',) + FASTA_SUFFIXES:
            name = root
        else:
            break
    return name


def format_number(num: int) -> str:
    """大于1百万用M单位2位小数(§5.3)|Numbers >1M use M unit, 2 decimals"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


def collect_genomes(path: str) -> List[str]:
    """收集基因组(单文件/目录两种形态)|Collect genomes (single file / directory)

    目录:glob *.fa/*.fna/*.fasta(大小写不敏感,不递归),重名 stem 报错
    |Directory: glob fa/fna/fasta (case-insensitive, no recursion), dup stems error
    """
    errors: List[str] = []
    genomes: List[str] = []

    if os.path.isdir(path):
        for name in sorted(os.listdir(path)):
            full = os.path.join(path, name)
            if not os.path.isfile(full):
                continue
            low = name.lower()
            if low.endswith(FASTA_GZ_SUFFIXES):
                errors.append(f"检测到压缩FASTA,请先解压|Compressed FASTA found, "
                              f"decompress first: {name}")
            elif low.endswith(FASTA_SUFFIXES):
                if os.path.getsize(full) == 0:
                    errors.append(f"基因组文件为空|Empty genome file: {full}")
                else:
                    genomes.append(os.path.abspath(full))
        stems = {}
        for g in genomes:
            stem = genome_name(g)
            if stem in stems:
                errors.append(f"样品名重复(不同后缀)|Duplicate sample name: "
                              f"{g} vs {stems[stem]}")
            else:
                stems[stem] = g
        if not genomes:
            errors.append(f"目录中未找到FASTA({', '.join(FASTA_SUFFIXES)})"
                          f"|No FASTA files in directory: {path}")
    elif os.path.isfile(path):
        low = path.lower()
        if low.endswith(FASTA_GZ_SUFFIXES):
            errors.append(f"不支持压缩FASTA,请先解压|Compressed FASTA not supported, "
                          f"decompress first: {path}")
        elif not low.endswith(FASTA_SUFFIXES):
            errors.append(f"输入应为 {', '.join(FASTA_SUFFIXES)} FASTA"
                          f"|Input must be a FASTA: {path}")
        elif os.path.getsize(path) == 0:
            errors.append(f"基因组文件为空|Empty genome file: {path}")
        else:
            genomes.append(os.path.abspath(path))
    else:
        errors.append(f"输入路径不存在|Input path not found: {path}")

    if errors:
        raise ValueError('\n'.join(errors))
    return genomes


class RepeatHit(NamedTuple):
    """RepeatMasker 逐条命中(1-based inclusive)|One RM hit (1-based inclusive)"""
    seqid: str
    start: int
    end: int
    family: str      # 库内家族名(.out 第9列)|library family name (.out col 9)
    family_class: str  # 类/超家族(.out 第10列)|class/superfamily (.out col 10)


def parse_repeatmasker_out(path: str) -> List[RepeatHit]:
    """解析 RepeatMasker .out(首3行头,数据行首列为整数)|Parse RM .out

    "left" 括号列含空格(如 "( 998577)")会被拆成多个 token,
    解析时跳过 end 之后所有 "(" 开头的 token 再取家族名/类
    |The "(left)" column splits into extra tokens; skip '('-prefixed
    tokens after "end" before reading family/class
    """
    hits: List[RepeatHit] = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            tokens = line.split()
            # 数据行:>=10字段且首列为整数;表头/空行跳过
            # |Data rows: >=10 fields with integer col 1; skip headers/blanks
            if len(tokens) < 10 or not tokens[0].lstrip('-').isdigit():
                continue
            try:
                seqid, start, end = tokens[4], int(tokens[5]), int(tokens[6])
            except (ValueError, IndexError):
                continue
            rest = tokens[7:]
            # "left" 列真实输出为 "( 998577)":括号与数字间有空格,拆成
            # "(" 与 "998577)" 两个 token,两个都得跳过才能取到家族名
            # |Real "left" col is "( 998577)": the space splits it into
            # "(" and "998577)", both must be skipped to reach the family
            while rest and (rest[0].startswith('(') or rest[0].endswith(')')):
                rest.pop(0)
            if len(rest) < 2:
                continue
            hits.append(RepeatHit(seqid=seqid, start=start, end=end,
                                  family=rest[0], family_class=rest[1]))
    return hits


def parse_repeatmasker_tbl(path: str) -> Dict[str, Optional[float]]:
    """解析 .tbl 总体统计|Parse .tbl overall stats"""
    text = open(path, encoding='utf-8').read()
    patterns = {
        'sequences': re.compile(r'sequences:\s+(\d+)'),
        'total_length': re.compile(r'total length:\s+(\d+) bp'),
        'masked_bp': re.compile(r'letters masked:\s+(\d+) bp'),
        'masked_pct': re.compile(r'letters masked:\s+\d+ bp \(\s*([\d.]+) %\)'),
    }
    stats: Dict[str, Optional[float]] = {}
    for key, pat in patterns.items():
        m = pat.search(text)
        stats[key] = (int(m.group(1)) if key != 'masked_pct'
                      else float(m.group(1))) if m else None
    return stats


def parse_tesorter_cls_tsv(path: str) -> Dict[str, Dict[str, str]]:
    """解析 TEsorter cls.tsv(#TE Order Superfamily Clade ...)|Parse TEsorter cls.tsv"""
    result: Dict[str, Dict[str, str]] = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#TE'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            result[fields[0]] = {'order': fields[1], 'superfamily': fields[2],
                                 'clade': fields[3]}
    return result


def fasta_ids(path: str) -> List[str]:
    """FASTA 头行 id(首个空白前)|FASTA header ids (before first whitespace)"""
    ids = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def _gff_attr(attr_field: str, key: str) -> Optional[str]:
    """取 GFF 属性值|Get GFF attribute value"""
    for item in attr_field.split(';'):
        k, _, v = item.partition('=')
        if k.strip() == key:
            return v.strip()
    return None


def load_effector_regions(path: str) -> List[Dict[str, object]]:
    """读 effector 候选区(BED/GFF3,统一1-based闭区间)|Load effector regions

    BED 为 0-based 半开,此处转 1-based inclusive 与 RepeatMasker 坐标对齐
    |BED is 0-based half-open; converted to 1-based inclusive to match RM coords
    """
    is_gff = path.lower().endswith(('.gff', '.gff3'))
    regions: List[Dict[str, object]] = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            if is_gff:
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                attrs = fields[8] if len(fields) > 8 else ''
                name = _gff_attr(attrs, 'ID') or _gff_attr(attrs, 'Name') \
                    or f"{fields[0]}:{fields[3]}-{fields[4]}"
                regions.append({'id': name, 'seqid': fields[0],
                                'start': int(fields[3]), 'end': int(fields[4])})
            else:
                fields = line.split()
                if len(fields) < 3:
                    raise ValueError(f"BED至少3列|BED needs 3+ columns: {line}")
                name = fields[3] if len(fields) > 3 and fields[3] != '.' \
                    else f"{fields[0]}:{fields[1]}-{fields[2]}"
                regions.append({'id': name, 'seqid': fields[0],
                                'start': int(fields[1]) + 1, 'end': int(fields[2])})
    return regions


def find_overlaps(region: Dict[str, object], hits: List[RepeatHit]) -> List[RepeatHit]:
    """同染色体且区间相交的命中|Hits on same seqid with intersecting intervals"""
    return [h for h in hits
            if h.seqid == region['seqid']
            and h.start <= region['end'] and region['start'] <= h.end]


def overlap_bp(region: Dict[str, object], hit: RepeatHit) -> int:
    """重叠碱基数(相邻=0)|Overlap bp (adjacent=0)"""
    lo = max(int(region['start']), hit.start)
    hi = min(int(region['end']), hit.end)
    return max(0, hi - lo + 1)


def nearest_repeat_distance(region: Dict[str, object],
                            hits: List[RepeatHit]) -> Optional[int]:
    """最近 repeat 距离:重叠=0,无同染色体=None|Nearest repeat distance"""
    if find_overlaps(region, hits):
        return 0
    distances = []
    for h in hits:
        if h.seqid != region['seqid']:
            continue
        if h.start > int(region['end']):
            # 距离=两区间之间夹的碱基数(gap-1),与测试口径一致
            # |Distance = bp strictly between the two intervals (gap-1)
            distances.append(h.start - int(region['end']) - 1)
        elif h.end < int(region['start']):
            distances.append(int(region['start']) - h.end - 1)
    return min(distances) if distances else None
