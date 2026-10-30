"""genome2sv SV 序列输出 + PAV 矩阵单测|Tests for SV sequence & PAV matrix output

纯函数层:FaidxReader seek、GT→PAV、五类型序列提取、符号化 fallback、PAV 表格式化。
|Pure-function layer: FaidxReader seek, GT→PAV, per-type extraction, symbolic
fallback, PAV table formatting.
"""
from pathlib import Path

import pytest

from biopytools.genome2sv.utils import (
    FaidxReader, extract_sv_sequence, format_pav_binary, format_pav_matrix,
    gt_present, parse_info_str, revcomp, stable_sv_id)

# 测试用参考:chrA 100bp、行宽 10,序列可预测|Test ref: chrA 100bp, width 10
_CHRA = "".join("ACGT"[i % 4] for i in range(100))


def _write_fasta_with_fai(tmp_path, seqs, line_width=10):
    """写 FASTA 并手算 .fai(免外部工具)|Write FASTA and hand-build .fai."""
    fa = tmp_path / "ref.fasta"
    fai_lines = []
    offset = 0
    with open(fa, "w") as fh:
        for name, seq in seqs.items():
            header = f">{name}\n"
            fh.write(header)
            offset += len(header.encode())
            fai_lines.append(f"{name}\t{len(seq)}\t{offset}\t"
                             f"{min(line_width, len(seq))}\t{line_width + 1}")
            for i in range(0, len(seq), line_width):
                chunk = seq[i:i + line_width]
                fh.write(chunk + "\n")
                offset += len(chunk.encode()) + 1
    Path(str(fa) + ".fai").write_text("\n".join(fai_lines) + "\n")
    return str(fa)


# ---------- FaidxReader ----------

def test_faidx_reader_fetch_within_line(tmp_path):
    """单行内取序列|Fetch within one line."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    reader = FaidxReader(fa)
    assert reader.fetch("chrA", 1, 4) == "ACGT"


def test_faidx_reader_fetch_spanning_lines(tmp_path):
    """跨行取序列(换行剥离)|Fetch spanning lines (newlines stripped)."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    reader = FaidxReader(fa)
    assert reader.fetch("chrA", 8, 15) == _CHRA[7:15]


def test_faidx_reader_fetch_full_contig(tmp_path):
    """整条染色体|Full contig fetch equals source."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    reader = FaidxReader(fa)
    assert reader.fetch("chrA", 1, 100) == _CHRA


def test_faidx_reader_clamps_end_to_length(tmp_path):
    """end 超染色体长度时截到末端|Clamp end to contig length."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    reader = FaidxReader(fa)
    assert reader.fetch("chrA", 98, 500) == _CHRA[97:]


def test_faidx_reader_missing_fai_raises(tmp_path):
    """.fai 缺失时抛 FileNotFoundError|Missing .fai raises."""
    fa = tmp_path / "ref.fasta"
    fa.write_text(">chrA\nACGT\n")
    with pytest.raises(FileNotFoundError):
        FaidxReader(str(fa))


# ---------- GT → PAV ----------

@pytest.mark.parametrize("gt,expected", [
    ("1/1", 1), ("0/1", 1), ("1", 1), ("1|0", 1),
    ("./.", 0), (".", 0), ("0/0", 0),
])
def test_gt_present_variants(gt, expected):
    """GT 含 1 → 1,缺失/纯参考 → 0|GT with allele 1 → 1 else 0."""
    assert gt_present(gt) == expected


def test_gt_present_ignores_trailing_fields():
    """只看冒号前 GT 子字段|Only the colon-free GT sub-field counts."""
    assert gt_present("1/1:NA:500:0,0") == 1
    assert gt_present("./.:NaN:0") == 0


# ---------- revcomp / stable id ----------

def test_revcomp():
    assert revcomp("AACCGGTT") == "AACCGGTT"
    assert revcomp("ACGTN") == "NACGT"


def test_stable_sv_id_padding_and_dup_normalization():
    """DUP 子类型归一为 DUP,编号 5 位补零|DUP subtypes normalize; 5-digit pad."""
    assert stable_sv_id("INS", 1) == "pan_sv.INS.00001"
    assert stable_sv_id("DEL", 123) == "pan_sv.DEL.00123"
    assert stable_sv_id("DUP:TANDEM", 7) == "pan_sv.DUP.00007"
    assert stable_sv_id("DUP:INT", 8) == "pan_sv.DUP.00008"


# ---------- extract_sv_sequence ----------

def test_extract_ins_from_alt(tmp_path):
    """INS:ALT 剥 anchor 前缀|INS: strip anchor prefix from ALT."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    seq, source = extract_sv_sequence(
        "INS", "chrA", 10, "A", "AGGGTTTAAA", {"END": 10}, FaidxReader(fa))
    assert seq == "GGGTTTAAA"
    assert source == "alt"


def test_extract_ins_symbolic_returns_none(tmp_path):
    """符号化 INS 无参考来源,跳过|Symbolic INS has no source; skip."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    assert extract_sv_sequence(
        "INS", "chrA", 10, "N", "<INS>", {"END": 10, "SVLEN": 300},
        FaidxReader(fa)) is None


def test_extract_del_prefers_ref_field(tmp_path):
    """DEL:REF 字段即被删序列|DEL: REF field carries deleted sequence."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    ref_allele = "A" + _CHRA[10:15]  # anchor + 被删 5bp|anchor + 5bp deleted
    seq, source = extract_sv_sequence(
        "DEL", "chrA", 11, ref_allele, "A", {"END": 15}, FaidxReader(fa))
    assert seq == ref_allele
    assert source == "ref"


def test_extract_del_symbolic_falls_back_to_region(tmp_path):
    """符号化 DEL:按 [POS,END] 从参考提取|Symbolic DEL: fetch [POS,END]."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    seq, source = extract_sv_sequence(
        "DEL", "chrA", 11, "N", "<DEL>", {"END": 15}, FaidxReader(fa))
    assert seq == _CHRA[10:15]
    assert source == "region"


def test_extract_inv_from_alt(tmp_path):
    """INV:ALT 已是 revcomp,直接用|INV: ALT is already revcomp."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    region = _CHRA[10:15]
    seq, source = extract_sv_sequence(
        "INV", "chrA", 11, region, revcomp(region), {"END": 15}, FaidxReader(fa))
    assert seq == revcomp(region)
    assert source == "alt"


def test_extract_inv_symbolic_falls_back_to_region_revcomp(tmp_path):
    """符号化 INV:区间提取 + revcomp|Symbolic INV: region + revcomp."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    seq, source = extract_sv_sequence(
        "INV", "chrA", 11, "N", "<INV>", {"END": 15}, FaidxReader(fa))
    assert seq == revcomp(_CHRA[10:15])
    assert source == "region_revcomp"


@pytest.mark.parametrize("svtype", ["DUP", "DUP:TANDEM", "DUP:INT"])
def test_extract_dup_always_from_region(tmp_path, svtype):
    """DUP 各型:统一按 [POS,END] 提取重复单元|All DUP flavors: unit region."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    seq, source = extract_sv_sequence(
        svtype, "chrA", 11, "N", "<" + svtype + ">", {"END": 15}, FaidxReader(fa))
    assert seq == _CHRA[10:15]
    assert source == "region"


def test_extract_bnd_returns_none(tmp_path):
    """BND 无区间,跳过|BND has no interval; skip."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    assert extract_sv_sequence(
        "BND", "chrA", 11, "N", "N[chrA:99[", {}, FaidxReader(fa)) is None


def test_extract_unknown_type_returns_none(tmp_path):
    """未知类型跳过|Unknown type skipped."""
    fa = _write_fasta_with_fai(tmp_path, {"chrA": _CHRA})
    assert extract_sv_sequence(
        "TRA", "chrA", 11, "N", "<TRA>", {"END": 15}, FaidxReader(fa)) is None


def test_extract_del_region_reader_failure_returns_none(tmp_path):
    """reader 初始化失败时优雅跳过|Graceful skip when reader is None."""
    assert extract_sv_sequence(
        "DEL", "chrA", 11, "N", "<DEL>", {"END": 15}, None) is None


# ---------- parse_info_str ----------

def test_parse_info_str():
    info = parse_info_str("SVTYPE=DEL;END=1499;SVLEN=-500")
    assert info == {"SVTYPE": "DEL", "END": "1499", "SVLEN": "-500"}


# ---------- PAV 表格式化 ----------

_ROWS = [
    ("pan_sv.INS.00001", "chr1", 100, 100, "INS", 102, [1, 0]),
    ("pan_sv.DEL.00002", "chr1", 999, 1499, "DEL", -500, [1, 1]),
]


def test_format_pav_matrix_tsv():
    out = format_pav_matrix(_ROWS, ["hap1", "hap2"])
    lines = out.strip().split("\n")
    assert lines[0] == "sv_id\tchrom\tpos\tend\tsvtype\tsvlen\thap1\thap2"
    assert lines[1] == ("pan_sv.INS.00001\tchr1\t100\t100\tINS\t102\t1\t0")
    assert lines[2] == ("pan_sv.DEL.00002\tchr1\t999\t1499\tDEL\t-500\t1\t1")


def test_format_pav_binary_tsv():
    out = format_pav_binary(_ROWS, ["hap1", "hap2"])
    lines = out.strip().split("\n")
    assert lines[0] == "sv_id\thap1\thap2"
    assert lines[1] == "pan_sv.INS.00001\t1\t0"
    assert lines[2] == "pan_sv.DEL.00002\t1\t1"
