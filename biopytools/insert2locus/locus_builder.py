"""完整locus判定与junction标注|Complete locus determination and junction annotation"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .data_processing import cigar_matched, cigar_softclips, is_junction_alignment, parse_cigar


@dataclass
class JunctionInfo:
    """junction contig信息|Junction contig info"""
    contig: str
    border: str          # L/R/LR/unmapped
    pos_on_insert: int
    anchor_bp: int
    flank_bp: int
    lead: int
    trail: int
    anchor_class: str = "insert"   # insert/backbone 锚区归属|anchor region class


def _ref_span(cigar: str) -> int:
    """CIGAR参考跨度(M/D/N/=/X耗参考)|Reference span of CIGAR"""
    return sum(n for op, n in parse_cigar(cigar) if op in ("M", "D", "N", "=", "X"))


def classify_anchor(pos: int, span: int,
                    insert_regions: Optional[List[Tuple[int, int]]]) -> str:
    """锚段归属:≥50%落insert区判insert,否则backbone|
    Anchor class: >=50% overlapping insert regions -> insert, else backbone"""
    if not insert_regions:
        return "insert"
    start, end = pos, pos + span - 1
    overlap = 0
    for r_start, r_end in insert_regions:
        overlap += max(0, min(end, r_end) - max(start, r_start) + 1)
    return "insert" if overlap * 2 >= span else "backbone"


def insert_regions_from_alignment(sam_lines) -> List[Tuple[int, int]]:
    """tdna比对构建的SAM→insert区段(合并重叠,弃secondary)|
    tdna-vs-construct SAM -> merged insert regions (secondary dropped)"""
    spans: List[Tuple[int, int]] = []
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        start = int(f[3])
        spans.append((start, start + _ref_span(f[5]) - 1))
    if not spans:
        return []
    spans.sort()
    merged = [list(spans[0])]
    for start, end in spans[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def max_flank_extent(sam_lines, contig_seqs: Dict[str, str]) -> int:
    """当前侧翼最远延伸(softclip段与未比对contig取最大)|
    Furthest flank extension: max of softclip ends and unmapped contigs"""
    best = 0
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4:
            seq = contig_seqs.get(f[0], f[9])
            best = max(best, len(seq))
            continue
        if flag & 256 or flag & 2048:
            continue
        lead, trail = cigar_softclips(f[5])
        best = max(best, lead, trail)
    return best


def max_matched_by_query(sam_lines) -> Dict[str, int]:
    """每个query的最大匹配长度(侧翼植物性校验用)|
    Max matched length per query (for flank plant check)"""
    out: Dict[str, int] = {}
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if int(f[1]) & 4:
            out.setdefault(f[0], 0)
            continue
        out[f[0]] = max(out.get(f[0], 0), cigar_matched(f[5]))
    return out


@dataclass
class CompleteLocus:
    """完整插入位点|Complete insertion locus"""
    contig_name: str
    seq: str
    lead: int           # LB长度|LB length
    trail: int          # RB长度|RB length
    insert_record: str
    insert_len: int
    insert_start: int   # locus内insert起点(1-based)|insert start in locus (1-based)
    insert_end: int     # locus内insert终点|insert end in locus


def _iter_alignments(sam_lines):
    """逐条产出比对列|Yield alignment columns line by line"""
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        yield {"qname": f[0], "flag": int(f[1]), "rname": f[2],
               "pos": int(f[3]), "cigar": f[5], "seq": f[9]}


def find_complete_locus(sam_lines, contig_seqs: Dict[str, str],
                        insert_record_lens: Dict[str, int],
                        min_anchor_frac: float = 0.95) -> Optional[CompleteLocus]:
    """找锚段覆盖主record>=95%且双端有翼的最长contig|
    Find longest contig anchored >=95% of its record with both-side flanks"""
    best: Optional[CompleteLocus] = None
    for aln in _iter_alignments(sam_lines):
        flag = aln["flag"]
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        rec_len = insert_record_lens.get(aln["rname"])
        if not rec_len:
            continue
        anchor = cigar_matched(aln["cigar"])
        if anchor < rec_len * min_anchor_frac:
            continue
        lead, trail = cigar_softclips(aln["cigar"])
        if lead == 0 or trail == 0:
            continue
        seq = contig_seqs.get(aln["qname"], aln["seq"])
        if best is None or len(seq) > len(best.seq):
            best = CompleteLocus(
                contig_name=aln["qname"], seq=seq, lead=lead, trail=trail,
                insert_record=aln["rname"], insert_len=rec_len,
                insert_start=lead + 1, insert_end=lead + anchor)
    return best


def annotate_junctions(sam_lines, contig_seqs: Dict[str, str],
                       junction_flank: int, min_unmapped: int,
                       insert_regions: Optional[List[Tuple[int, int]]] = None
                       ) -> List[JunctionInfo]:
    """标注全部junction contig(按flank_bp倒序)|Annotate all junction contigs (desc by flank_bp)"""
    out: List[JunctionInfo] = []
    for aln in _iter_alignments(sam_lines):
        flag = aln["flag"]
        seq = contig_seqs.get(aln["qname"], aln["seq"])
        if flag & 4:
            if len(seq) >= min_unmapped:
                out.append(JunctionInfo(
                    contig=aln["qname"], border="unmapped", pos_on_insert=0,
                    anchor_bp=0, flank_bp=len(seq), lead=0, trail=0))
            continue
        if flag & 256 or flag & 2048:
            continue
        lead, trail = cigar_softclips(aln["cigar"])
        border = is_junction_alignment(aln["cigar"], flag, junction_flank)
        if border:
            out.append(JunctionInfo(
                contig=aln["qname"], border=border, pos_on_insert=aln["pos"],
                anchor_bp=cigar_matched(aln["cigar"]),
                flank_bp=lead + trail, lead=lead, trail=trail,
                anchor_class=classify_anchor(
                    aln["pos"], _ref_span(aln["cigar"]), insert_regions)))
    out.sort(key=lambda j: j.flank_bp, reverse=True)
    return out


def write_junction_outputs(junctions: List[JunctionInfo],
                           contig_seqs: Dict[str, str],
                           locus: Optional[CompleteLocus], out_dir: Path,
                           sample: str) -> Dict[str, Path]:
    """写junction/locus产物|Write junction & locus outputs"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    junc_fa = out_dir / f"{sample}_junction_contigs.fasta"
    junc_tsv = out_dir / f"{sample}_junction_report.tsv"
    locus_fa = out_dir / f"{sample}_complete_locus.fasta"
    with open(junc_fa, "w") as fh:
        for j in junctions:
            fh.write(f">{j.contig} border={j.border} pos={j.pos_on_insert} "
                     f"anchor={j.anchor_bp}bp flank={j.flank_bp}bp\n")
            fh.write(contig_seqs[j.contig] + "\n")
    with open(junc_tsv, "w") as fh:
        fh.write("contig\tborder\tpos_on_insert\tanchor_bp\tflank_bp\t"
                 "anchor_class\n")
        for j in junctions:
            fh.write(f"{j.contig}\t{j.border}\t{j.pos_on_insert}\t"
                     f"{j.anchor_bp}\t{j.flank_bp}\t{j.anchor_class}\n")
    if locus:
        with open(locus_fa, "w") as fh:
            fh.write(f">{sample}_complete_locus LB={locus.lead}bp "
                     f"insert={locus.insert_end - locus.insert_start + 1}bp "
                     f"RB={locus.trail}bp\n{locus.seq}\n")
    return {"junction_fasta": junc_fa, "junction_tsv": junc_tsv,
            "locus_fasta": locus_fa}
