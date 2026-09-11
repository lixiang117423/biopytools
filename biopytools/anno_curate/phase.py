"""CDS相位算法|CDS phase algorithm

phase[i] = (cumLen % 3 == 0) ? 0 : 3 - cumLen % 3,cumLen=当前 CDS 之前
所有 CDS 长度之和;只信 CDS 长度不信原 phase(spec §6)
|phase from cumulative CDS length only, original phase never trusted (spec §6)
"""

from dataclasses import dataclass, field
from typing import List

from .gff_model import GffRecord


def order_cds_by_strand(cds_records: List[GffRecord]) -> List[GffRecord]:
    """链方向定 CDS 顺序|Order CDS by strand"""
    strand = cds_records[0].strand if cds_records else '+'
    return sorted(cds_records, key=lambda r: r.start,
                  reverse=(strand == '-'))


def compute_expected_phases(cds_records: List[GffRecord]) -> List[str]:
    """期望相位公式|Expected-phase formula"""
    phases = []
    cum_len = 0
    for i, rec in enumerate(cds_records):
        if i == 0:
            phases.append('0')
        else:
            m = cum_len % 3
            phases.append('0' if m == 0 else str(3 - m))
        cum_len += rec.end - rec.start + 1
    return phases


@dataclass
class PhaseResult:
    """转录本相位判定结果|Per-transcript phase verdict"""

    status: str                       # VALID | CORRECTED | INVALID_LENGTH
    expected: List[str] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """VALID/CORRECTED 均算有效|Both VALID and CORRECTED are valid"""
        return self.status in ('VALID', 'CORRECTED')

    @property
    def is_problematic(self) -> bool:
        """仅 INVALID_LENGTH 需人工修|Only INVALID_LENGTH needs manual fix"""
        return self.status == 'INVALID_LENGTH'


def validate_transcript_phases(cds_records: List[GffRecord]) -> PhaseResult:
    """三态判定(不改输入)|Three-state verdict (input untouched)"""
    ordered = order_cds_by_strand(cds_records)
    total = sum(r.end - r.start + 1 for r in ordered)
    if total % 3 != 0:
        return PhaseResult(
            'INVALID_LENGTH',
            issues=[f'CDS 总长非 3 倍数|Total CDS length not multiple of 3: '
                    f'{total}'])
    expected = compute_expected_phases(ordered)
    issues = []
    mismatch = False
    for rec, exp in zip(ordered, expected):
        if rec.phase == '.':
            mismatch = True
            issues.append(f'NO PHASE(第{rec.line_no}行 phase 缺失)'
                          f'|missing phase at line {rec.line_no}')
        elif rec.phase != exp:
            mismatch = True
            issues.append(f"INVALID PHASE 第{rec.line_no}行 phase={rec.phase} "
                          f"期望|expected {exp}")
    if mismatch:
        return PhaseResult('CORRECTED', expected=expected, issues=issues)
    return PhaseResult('VALID', expected=expected)


def apply_phases(cds_records: List[GffRecord], phases: List[str]):
    """就地覆盖 phase|Overwrite phases in place"""
    for rec, ph in zip(cds_records, phases):
        rec.phase = ph
