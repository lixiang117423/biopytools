"""gff_model 单测|gff_model tests"""

import pytest

from biopytools.anno_curate.gff_model import (
    GffRecord, GxfDocument, feature_rank, parse_gxf, sort_records,
    write_gff3)

GFF3_MINIMAL = """##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\texon\t100\t200\t.\t+\t.\tParent=t1
chr1\t.\tCDS\t100\t195\t.\t+\t0\tParent=t1
"""


@pytest.fixture
def minimal_doc(tmp_path):
    p = tmp_path / 'in.gff3'
    p.write_text(GFF3_MINIMAL)
    return parse_gxf(str(p))


def test_parse_minimal_gff3(minimal_doc):
    assert len(minimal_doc.records) == 4
    r = minimal_doc.records[3]
    assert (r.seqid, r.type, r.start, r.end, r.phase) == \
        ('chr1', 'CDS', 100, 195, '0')
    assert r.attributes == {'Parent': 't1'}


def test_parse_swaps_reversed_coords(tmp_path):
    p = tmp_path / 'rev.gff3'
    p.write_text("chr1\t.\tCDS\t200\t100\t.\t+\t0\tParent=t1\n")
    doc = parse_gxf(str(p))
    assert (doc.records[0].start, doc.records[0].end) == (100, 200)


def test_parse_attribute_order_preserved(tmp_path):
    p = tmp_path / 'attr.gff3'
    p.write_text("chr1\t.\texon\t1\t9\t.\t+\t.\tID=e1;Parent=t1;Note=x\n")
    doc = parse_gxf(str(p))
    assert list(doc.records[0].attributes) == ['ID', 'Parent', 'Note']


def test_parse_bad_columns_collected(tmp_path):
    p = tmp_path / 'bad.gff3'
    p.write_text("chr1\t.\tCDS\t1\t9\t.\t+\t0\n"          # 8 列
                 "chr1\t.\tCDS\tone\t9\t.\t+\t0\tParent=t\n")  # 坐标非数字
    with pytest.raises(ValueError) as ei:
        parse_gxf(str(p))
    msg = str(ei.value)
    assert '1' in msg and '2' in msg   # 两个行号都在


def test_case_variant_detected(tmp_path):
    p = tmp_path / 'lower.gff3'
    p.write_text("chr1\t.\tcds\t1\t9\t.\t+\t0\tParent=t1\n")
    doc = parse_gxf(str(p))
    assert len(doc.case_variants) == 1
    assert 'cds' in doc.case_variants[0]


def test_parse_gtf_converted(tmp_path):
    p = tmp_path / 'in.gtf'
    p.write_text("chr1\tsrc\texon\t100\t200\t.\t+\t.\t"
                 "gene_id \"g1\"; transcript_id \"t1\";\n")
    doc = parse_gxf(str(p))
    r = doc.records[0]
    assert r.type == 'exon'
    assert r.attributes['Parent'] == 't1'
    assert r.attributes['gene_id'] == 'g1'
    assert 'transcript_id' not in r.attributes


def test_parse_gtf_transcript_row_becomes_id_and_parent(tmp_path):
    # I-1 回归:GTF transcript 行自身是转录本级,transcript_id→ID、
    # gene_id→Parent;不能再按子特征规则把 transcript_id 塞进 Parent
    # (否则该行成为无 ID 孤儿,⑦ 拆假基因、all_tids 翻倍)
    # |I-1 regression: a GTF transcript row IS transcript-level;
    # transcript_id becomes ID, gene_id becomes Parent
    p = tmp_path / 'in.gtf'
    p.write_text(
        "chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\t"
        "gene_id \"g1\"; transcript_id \"t1\";\n"
        "chr1\tsrc\texon\t100\t200\t.\t+\t.\t"
        "gene_id \"g1\"; transcript_id \"t1\";\n")
    doc = parse_gxf(str(p))
    assert doc.records[0].attributes == {'ID': 't1', 'Parent': 'g1'}
    exon = doc.records[1]                     # 子特征行规则不变
    assert exon.attributes['Parent'] == 't1'
    assert exon.attributes['gene_id'] == 'g1'


def test_parse_gtf_transcript_row_without_gene_id(tmp_path):
    # 无 gene_id 的 transcript 行:只获得 ID,不造 Parent
    # |Transcript row without gene_id: ID only, no Parent invented
    p = tmp_path / 'in.gtf'
    p.write_text("chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\t"
                 "transcript_id \"t1\";\n")
    doc = parse_gxf(str(p))
    assert doc.records[0].attributes == {'ID': 't1'}


def test_parse_gtf_attrs_without_trailing_semicolon(tmp_path):
    # M-5 回归:末个属性无分号不再丢(仍要求 key "value" 引号形态)
    # |M-5 regression: a final attribute without a semicolon is kept
    p = tmp_path / 'in.gtf'
    p.write_text("chr1\tsrc\texon\t100\t200\t.\t+\t.\t"
                 "gene_id \"g1\"; transcript_id \"t1\"\n")
    doc = parse_gxf(str(p))
    r = doc.records[0]
    assert r.attributes['Parent'] == 't1'
    assert r.attributes['gene_id'] == 'g1'


def test_fasta_block_passthrough(tmp_path):
    p = tmp_path / 'fa.gff3'
    p.write_text(GFF3_MINIMAL + "##FASTA\n>chr1\nACGT\n")
    doc = parse_gxf(str(p))
    assert doc.fasta_block == '##FASTA\n>chr1\nACGT\n'
    assert len(doc.records) == 4


def test_indexes(minimal_doc):
    minimal_doc.rebuild_indexes()
    assert 't1' in minimal_doc.id_index
    assert [r.type for r in minimal_doc.children_by_parent['t1']] == \
        ['exon', 'CDS']


def test_feature_rank_transcript_level():
    assert feature_rank('gene') == 1
    assert feature_rank('mRNA') == 2
    assert feature_rank('transcript') == 3
    assert feature_rank('tRNA') == 3      # 非编码转录本级不落 rank 8
    assert feature_rank('lnc_RNA') == 3
    assert feature_rank('exon') == 4
    assert feature_rank('five_prime_UTR') == 5
    assert feature_rank('CDS') == 6
    assert feature_rank('three_prime_UTR') == 7
    assert feature_rank('bioregion') == 8


def test_sort_chr_first_seen_order():
    chr_order = ['chr10', 'chr2']
    recs = [
        GffRecord('chr2', '.', 'gene', 50, 60, '.', '+', '.', {}, 1),
        GffRecord('chr10', '.', 'gene', 500, 600, '.', '+', '.', {}, 2),
    ]
    out = sort_records(recs, chr_order)
    assert [r.seqid for r in out] == ['chr10', 'chr2']   # 首现序非字典序


def test_write_gff3_roundtrip(minimal_doc, tmp_path):
    out = tmp_path / 'out.gff3'
    write_gff3(minimal_doc, minimal_doc.records, str(out))
    text = out.read_text()
    assert text.startswith('##gff-version 3\n')
    assert '\tCDS\t100\t195\t.\t+\t0\tParent=t1\n' in text
