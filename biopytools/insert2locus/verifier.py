"""覆盖验证与分级|Coverage verification and grading"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .data_processing import parse_cigar


@dataclass
class SegmentStat:
    """区段覆盖统计|Segment coverage statistics"""
    name: str
    start: int
    end: int
    length: int
    mean_depth: float
    min_depth: int
    zero_windows: int


def _zero_windows(depths: List[int]) -> int:
    """零覆盖窗口数(连续0算1个)|Count zero-coverage windows (consecutive 0 = 1)"""
    n = 0
    prev_zero = False
    for d in depths:
        if d == 0 and not prev_zero:
            n += 1
        prev_zero = d == 0
    return n


def segment_coverage(depth_rows, locus) -> List[SegmentStat]:
    """LB/insert/RB分段统计|Per-segment (LB/insert/RB) statistics"""
    if locus is None:
        return []
    by_pos = dict(depth_rows)
    bounds = [
        ("LB", 1, locus.lead),
        ("insert", locus.lead + 1, locus.lead + locus.insert_len),
        ("RB", locus.lead + locus.insert_len + 1,
         locus.lead + locus.insert_len + locus.trail),
    ]
    segs = []
    for name, start, end in bounds:
        depths = [by_pos.get(p, 0) for p in range(start, end + 1)]
        segs.append(SegmentStat(
            name=name, start=start, end=end, length=max(len(depths), 0),
            mean_depth=round(sum(depths) / len(depths), 1) if depths else 0.0,
            min_depth=min(depths) if depths else 0,
            zero_windows=_zero_windows(depths)))
    return segs


def classify(segments: List[SegmentStat], lb_junction_reads: int,
             rb_junction_reads: int, junction_flank: int, locus,
             flanks_plant: Optional[Tuple[bool, bool]] = None) -> str:
    """PASS/WARN/FAIL分级|Grade PASS/WARN/FAIL"""
    if locus is None or not segments:
        return "FAIL"
    if any(s.zero_windows > 0 for s in segments):
        return "WARN"
    if locus.lead < junction_flank or locus.trail < junction_flank:
        return "WARN"
    if lb_junction_reads == 0 or rb_junction_reads == 0:
        return "WARN"
    if flanks_plant is not None and not all(flanks_plant):
        return "WARN"   # 侧翼匹配构建=载体骨架来源|Flank matches construct
    return "PASS"


def region_depth(depth_rows,
                 regions: Dict[str, List[Tuple[int, int]]]) -> Dict[str, dict]:
    """分区覆盖(insert区/骨架区,backbone信号)|Region coverage (insert vs backbone)"""
    by_pos: Dict[Tuple[str, int], int] = {}
    for record, pos, depth in depth_rows:
        by_pos[(record, pos)] = depth
    out = {}
    for label, spans in regions.items():
        positions = []
        for start, end in spans:
            positions.extend(range(start, end + 1))
        if not positions:
            out[label] = {"mean": 0.0, "covered_frac": 0.0, "length": 0}
            continue
        # depth_rows可能多record;跨record按位置累加(单record构建场景为主)|
        # Sum depth across records per position (single-record construct typical)
        depth_by_pos: Dict[int, int] = {}
        for (_rec, pos), d in by_pos.items():
            depth_by_pos[pos] = depth_by_pos.get(pos, 0) + d
        depths = [depth_by_pos.get(p, 0) for p in positions]
        covered = sum(1 for d in depths if d > 0)
        out[label] = {
            "mean": round(sum(depths) / len(depths), 1),
            "covered_frac": round(covered / len(depths), 4),
            "length": len(depths),
        }
    return out


def per_record_depth(depth_rows, total_lens: Dict[str, int]) -> Dict[str, dict]:
    """逐record覆盖(backbone线索)|Per-record coverage (backbone signal)"""
    agg: Dict[str, Dict[str, int]] = {}
    for record, _pos, depth in depth_rows:
        a = agg.setdefault(record, {"sum": 0, "n": 0, "cov": 0})
        a["sum"] += depth
        a["n"] += 1
        if depth > 0:
            a["cov"] += 1
    out = {}
    for record, a in agg.items():
        total = total_lens.get(record, a["n"])
        out[record] = {
            "mean": round(a["sum"] / a["n"], 1) if a["n"] else 0.0,
            "covered_frac": round(a["cov"] / total, 4) if total else 0.0,
            "length": total,
        }
    return out


def count_junction_reads(sam_lines, boundary: int) -> int:
    """跨边界reads数:比对区间覆盖boundary且CIGAR含S|
    Reads spanning the boundary: alignment covers boundary and CIGAR has S"""
    n = 0
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4 or flag & 256:
            continue
        cigar = f[5]
        if "S" not in cigar:
            continue
        ops = parse_cigar(cigar)
        ref_start = int(f[3])
        ref_span = sum(length for op, length in ops if op in ("M", "D", "=", "X"))
        ref_end = ref_start + ref_span - 1
        if ref_start <= boundary <= ref_end:
            n += 1
    return n
