"""junction reads钓取与诱饵抽取|Junction read fishing and bait extraction"""

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def read_fasta_dict(fasta) -> Dict[str, str]:
    """读fasta成{name:seq}|Read fasta into {name:seq}"""
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    """解析CIGAR为(op,len)列表|Parse CIGAR into (op,len) list"""
    if cigar == "*":
        return []
    return [(op, int(num)) for num, op in _CIGAR_RE.findall(cigar)]


def cigar_softclips(cigar: str) -> Tuple[int, int]:
    """首尾soft-clip长度|(lead,trail) soft-clip lengths"""
    ops = parse_cigar(cigar)
    if not ops:
        return (0, 0)
    lead = ops[0][1] if ops[0][0] == "S" else 0
    trail = ops[-1][1] if ops[-1][0] == "S" else 0
    return (lead, trail)


def cigar_matched(cigar: str) -> int:
    """M/=/X匹配总长|Total matched length (M/=/X)"""
    return sum(n for op, n in parse_cigar(cigar) if op in ("M", "=", "X"))


def extract_flank_bait(sam_lines: Iterable[str], min_softclip: int,
                       min_unmapped: int) -> List[Tuple[str, str]]:
    """抽纯侧翼诱饵:未比对整条+首尾soft-clip段(弃insert锚段)|
    Pure-flank bait: whole unmapped contigs + soft-clip ends (insert anchors dropped)"""
    baits: List[Tuple[str, str]] = []
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        qname, flag, cigar, seq = f[0], int(f[1]), f[5], f[9]
        if flag & 4:
            if len(seq) >= min_unmapped:
                baits.append((qname, seq))
            continue
        lead, trail = cigar_softclips(cigar)
        if lead >= min_softclip:
            baits.append((f"{qname}_L", seq[:lead]))
        if trail >= min_softclip:
            baits.append((f"{qname}_R", seq[len(seq) - trail:]))
    return baits


def is_junction_alignment(cigar: str, flag: int, min_flank: int) -> Optional[str]:
    """边界类型L/R/LR|Border type L/R/LR"""
    if flag & 4:
        return None
    lead, trail = cigar_softclips(cigar)
    l_ok, r_ok = lead >= min_flank, trail >= min_flank
    if l_ok and r_ok:
        return "LR"
    if l_ok:
        return "L"
    if r_ok:
        return "R"
    return None


def mate_suffix_names(read_names: Iterable[str], mate: int) -> List[str]:
    """补回/1或/2后缀(bwa去掉了原始fastq的后缀)|
    Append /1 or /2 suffix (bwa strips it from fastq headers)"""
    return [f"{name}/{mate}" for name in read_names]
