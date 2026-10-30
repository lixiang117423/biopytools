"""phase 单测|phase tests"""

from biopytools.anno_curate.gff_model import GffRecord
from biopytools.anno_curate.phase import (PhaseResult, apply_phases,
                                          compute_expected_phases,
                                          order_cds_by_strand,
                                          validate_transcript_phases)


def cds(start, end, strand='+', phase='.', seqid='chr1', line_no=0):
    return GffRecord(seqid, '.', 'CDS', start, end, '.', strand, phase,
                     {'Parent': 't1'}, line_no)


def test_expected_phases_basic():
    # 长度 3,4,5,6:cumLen 依次 0,3,7,12 → phase = 0,0,(3-1)=2,0
    recs = [cds(1, 3), cds(10, 13), cds(20, 24), cds(30, 35)]
    assert compute_expected_phases(recs) == ['0', '0', '2', '0']


def test_expected_phases_all_multiple_of_3():
    recs = [cds(1, 3), cds(10, 12), cds(20, 25)]
    assert compute_expected_phases(recs) == ['0', '0', '0']


def test_order_plus_strand_ascending():
    recs = [cds(50, 60), cds(1, 10)]
    assert [r.start for r in order_cds_by_strand(recs)] == [1, 50]


def test_order_minus_strand_descending():
    recs = [cds(1, 10, '-'), cds(50, 60, '-')]
    assert [r.start for r in order_cds_by_strand(recs)] == [50, 1]


def test_invalid_length():
    recs = [cds(1, 3), cds(10, 13)]          # 总长 3+4=7 非 3 倍数
    res = validate_transcript_phases(recs)
    assert res.status == 'INVALID_LENGTH'
    assert res.is_problematic and not res.is_valid


def test_valid():
    # 长度 4+5=9(3 倍数),期望 phase 0,2
    recs = [cds(1, 4, phase='0'), cds(10, 14, phase='2')]
    res = validate_transcript_phases(recs)
    assert res.status == 'VALID'
    assert res.is_valid and not res.is_problematic


def test_corrected_on_missing_or_wrong():
    for p in ('.', '1'):                      # 缺失与错误同等对待
        # 长度 4+5=9(3 倍数),期望 phase 0,2
        recs = [cds(1, 4, phase='0'), cds(10, 14, phase=p)]
        res = validate_transcript_phases(recs)
        assert res.status == 'CORRECTED'
        assert res.is_valid
        assert res.expected == ['0', '2']


def test_missing_phase_flagged_in_issues():
    # 长度 4+5=9(3 倍数),进入 CORRECTED 分支才检 phase 问题
    recs = [cds(1, 4, phase='.'), cds(10, 14, phase='2')]
    res = validate_transcript_phases(recs)
    assert any('NO PHASE' in i for i in res.issues)


def test_apply_phases_overwrites():
    recs = [cds(1, 3, phase='.'), cds(10, 13, phase='9')]
    apply_phases(recs, ['0', '2'])
    assert [r.phase for r in recs] == ['0', '2']
