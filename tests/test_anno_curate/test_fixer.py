"""fixer 单测|fixer tests"""

import logging

import pytest

from biopytools.anno_curate.fixer import Gxffixer
from biopytools.anno_curate.gff_model import parse_gxf


def make_fixer(tmp_path, text, name='in.gff3'):
    p = tmp_path / name
    p.write_text(text)
    doc = parse_gxf(str(p))
    return Gxffixer(doc, logging.getLogger('biopytools.anno_curate.test'))


def test_circular_parent_rewritten(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=t1\n"
        "chr1\t.\tCDS\t100\t195\t.\t+\t0\tParent=t1\n")
    fixer.check_id_parent_conflicts()
    mrna = fixer.doc.records[0]
    assert mrna.attributes['Parent'] == 't1.gene'   # 改 mRNA 自身 Parent


def test_duplicate_transcript_id_hard_fail(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tmRNA\t300\t400\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tmRNA\t500\t600\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\t.\tmRNA\t700\t800\t.\t+\t.\tID=t2;Parent=g1\n")
    with pytest.raises(ValueError) as ei:
        fixer.check_id_parent_conflicts()
    msg = str(ei.value)
    assert 't1' in msg and 't2' in msg   # 全部重复一次列完


def test_gene_rows_removed_when_no_mrna(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1\n"
        "chr1\t.\texon\t100\t200\t.\t+\t.\tID=e1\n")
    fixer.check_usage_of_gene_feature()
    assert [r.type for r in fixer.doc.records] == ['exon']


def test_gene_rows_kept_when_mrna_exists(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1\n")
    fixer.check_usage_of_gene_feature()
    assert len(fixer.doc.records) == 2


def test_dangling_mrna_split(tmp_path):
    # 计划笔误修正|Plan typo fix:悬空的是无孩子的 t1(spec §5⑨),连同独占父
    # gene g1 拆出;带 CDS 孩子的 t2 留主文件。原断言把两对 ID 写反,且 kept
    # 列含无 ID 的 CDS 会 KeyError
    # |The childless t1 (plus its exclusive gene g1) is the dangling one; t2
    # with a CDS child stays. Original assertions swapped the two pairs, and
    # the kept list contains the ID-less CDS (KeyError with ['ID'])
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tgene\t700\t900\t.\t+\t.\tID=g2\n"
        "chr1\t.\tmRNA\t700\t800\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\t.\tCDS\t700\t795\t.\t+\t0\tParent=t2\n")
    dangling = fixer.split_dangling_mrna()
    assert sorted(r.attributes['ID'] for r in dangling) == ['g1', 't1']
    kept_ids = [r.attributes.get('ID') for r in fixer.doc.records]
    assert kept_ids == ['g2', 't2', None]   # 末位 None=t2 的 CDS(无ID)|None=CDS w/o ID


def test_dangling_gene_kept_if_has_other_normal_mrna(tmp_path):
    # g1 同时有一个悬空 t1 与一个正常 t2 → g1 留主文件
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tmRNA\t700\t800\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\t.\tCDS\t700\t795\t.\t+\t0\tParent=t2\n")
    dangling = fixer.split_dangling_mrna()
    assert [r.attributes['ID'] for r in dangling] == ['t1']
    kept_ids = [r.attributes.get('ID') for r in fixer.doc.records]
    assert 'g1' in kept_ids and 't2' in kept_ids


def test_dangling_parent_first_value_of_multi_parent(tmp_path):
    # C-1 预建索引口径:Parent 逗号多值只认首值(与 rebuild_indexes
    # 一致)——t1 Parent=g1,gx 记在 g1 名下,g1 随悬空 t1 带走
    # |C-1 prebuilt-index semantics: first value of a comma-joined
    # Parent wins (same as rebuild_indexes); g1 leaves with t1
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g1,gx\n"
        "chr1\t.\tmRNA\t700\t800\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\t.\tCDS\t700\t795\t.\t+\t0\tParent=t2\n")
    dangling = fixer.split_dangling_mrna()
    ids = sorted(r.attributes['ID'] for r in dangling
                 if r.attributes.get('ID'))
    assert ids == ['g1', 't1']


def test_recall_mrna_for_orphan_children(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tCDS\t100\t195\t.\t+\t0\tParent=t1;gene_id=g1\n"
        "chr1\t.\tCDS\t300\t395\t.\t+\t0\tParent=t1;gene_id=g1\n")
    fixer.recall_mrna_features()
    fixer.doc.rebuild_indexes()
    mrna = fixer.doc.id_index.get('t1')
    assert mrna is not None
    assert mrna.type == 'mRNA'
    assert mrna.source == 'Recall'
    assert (mrna.start, mrna.end) == (100, 395)      # 子特征 min/max
    assert mrna.attributes['Parent'] == 'g1'         # gene_id 回退
    assert 'Parent' not in mrna.attributes or \
        mrna.attributes['Parent'] == 'g1'


def test_recall_mrna_parent_fallback_id_dash_mrna(tmp_path):
    # 子特征无 Parent → 用 ID+"-mRNA" 当转录本号
    fixer = make_fixer(tmp_path,
        "chr1\t.\tCDS\t100\t195\t.\t+\t0\tID=c1;gene_id=g1\n")
    fixer.recall_mrna_features()
    fixer.doc.rebuild_indexes()
    mrna = fixer.doc.id_index.get('c1-mRNA')
    assert mrna is not None and mrna.source == 'Recall'


def test_existing_transcript_untouched(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t50\t400\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t195\t.\t+\t0\tParent=t1\n")
    fixer.recall_mrna_features()
    fixer.doc.rebuild_indexes()
    mrna = fixer.doc.id_index['t1']
    assert mrna.source != 'Recall'                   # 已存在的 mRNA 不动


def test_fix_gene_builds_missing_gene(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tmRNA\t700\t900\t.\t+\t.\tID=t2;geneID=g9\n")
    fixer.fix_gene_features()
    fixer.doc.rebuild_indexes()
    g1 = fixer.doc.id_index.get('g1')
    g9 = fixer.doc.id_index.get('g9')
    assert g1 is not None and (g1.start, g1.end) == (100, 500)
    assert g9 is not None and (g9.start, g9.end) == (700, 900)
    assert g1.source == 'anno_curate'
    # t2 的 Parent 被接到 g9
    t2 = fixer.doc.id_index['t2']
    assert t2.attributes['Parent'] == 'g9'


def test_fix_gene_new_id_dash_gene(tmp_path):
    # mRNA 无 Parent 无 geneID → 新建 <mRNA_ID>.gene 并接上
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1\n")
    fixer.fix_gene_features()
    fixer.doc.rebuild_indexes()
    assert fixer.doc.id_index.get('t1.gene') is not None
    assert fixer.doc.id_index['t1'].attributes['Parent'] == 't1.gene'


def test_recall_mrna_no_fake_parent_on_gene(tmp_path):
    # gene 行(有 ID 无 Parent)不得被写上伪 Parent(GFF3 根特征带
    # Parent 非法);同文件无 Parent 的 CDS 正常回退
    # |Gene rows must not get a fake Parent; the CDS still falls back
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t50\t600\t.\t+\t.\tID=g1\n"
        "chr1\t.\tCDS\t100\t195\t.\t+\t0\tID=c1\n")
    fixer.recall_mrna_features()
    gene = fixer.doc.records[0]
    cds = fixer.doc.records[1]
    assert 'Parent' not in gene.attributes
    assert cds.attributes['Parent'] == 'c1-mRNA'


def test_utr_from_mrna_minus_cds(tmp_path):
    # 无 exon 无 UTR → mRNA−CDS 差集:+链左5'右3'
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t200\t399\t.\t+\t0\tParent=t1\n")
    fixer.fix_utr_features()
    utrs = sorted((r for r in fixer.doc.records
                   if 'UTR' in r.type), key=lambda r: r.start)
    assert [(t.type, t.start, t.end) for t in utrs] == \
        [('five_prime_UTR', 100, 199), ('three_prime_UTR', 400, 500)]


def test_utr_minus_strand_swapped(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t-\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t200\t399\t.\t-\t0\tParent=t1\n")
    fixer.fix_utr_features()
    utrs = {r.type: r for r in fixer.doc.records if 'UTR' in r.type}
    assert utrs['five_prime_UTR'].start == 400      # -链右端=5'
    assert utrs['three_prime_UTR'].end == 199


def test_utr_from_exon_minus_cds_with_adjacent_shrink(tmp_path):
    # 有 exon → 逐 exon 拆;exonEnd==cdsStart 相切缩 1bp
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t1\t199\t.\t+\t.\tParent=t1\n"
        "chr1\t.\texon\t200\t399\t.\t+\t.\tParent=t1\n"
        "chr1\t.\texon\t400\t500\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t200\t399\t.\t+\t0\tParent=t1\n")
    fixer.fix_utr_features()
    utrs = sorted((r for r in fixer.doc.records
                   if 'UTR' in r.type), key=lambda r: r.start)
    assert [(t.type, t.start, t.end) for t in utrs] == \
        [('five_prime_UTR', 1, 199), ('three_prime_UTR', 400, 500)]


def test_utr_exon_inside_cds_skipped(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t210\t300\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t200\t399\t.\t+\t0\tParent=t1\n")
    fixer.fix_utr_features()
    assert not [r for r in fixer.doc.records if 'UTR' in r.type]


def test_exon_merged_from_cds_utr(tmp_path):
    # CDS 与 UTR gap=1bp 也合并
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t200\t399\t.\t+\t0\tParent=t1\n"
        "chr1\t.\tfive_prime_UTR\t100\t199\t.\t+\t.\tParent=t1\n")
    fixer.fix_exon_features()
    exons = [r for r in fixer.doc.records if r.type == 'exon']
    assert len(exons) == 1
    assert (exons[0].start, exons[0].end) == (100, 399)
    assert exons[0].source == 'anno_curate'


def test_exon_not_rebuilt_when_present(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t1\t250\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t200\t399\t.\t+\t0\tParent=t1\n")
    before = [r for r in fixer.doc.records if r.type == 'exon']
    fixer.fix_exon_features()
    after = [r for r in fixer.doc.records if r.type == 'exon']
    assert before == after                     # 已有 exon 不重建


def test_duplicated_child_ids_suffix(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\texon\t100\t200\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\t.\texon\t300\t400\t.\t+\t.\tID=e1;Parent=t2\n"
        "chr1\t.\texon\t500\t600\t.\t+\t.\tID=e2;Parent=t3\n")
    fixer.fix_duplicated_child_ids()
    ids = sorted(r.attributes['ID'] for r in fixer.doc.records)
    assert ids == ['e1', 'e1.u1', 'e2']         # 第一个保留,后续 .u<N>


def test_recall_mrna_no_fake_parent_on_bare_mrna(tmp_path):
    # 无 Parent 的 mRNA 行不回退,留给 ⑦ 按回退链建 t1.gene
    # |Bare mRNA rows don't fall back; ⑦ builds t1.gene instead
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1\n")
    fixer.recall_mrna_features()
    t1 = fixer.doc.records[0]
    assert 'Parent' not in t1.attributes
    fixer.fix_gene_features()
    assert t1.attributes['Parent'] == 't1.gene'


def test_run_full_pipeline(tmp_path):
    # brief 笔误修正|Brief typo fix:stats['dangling_mrna'] 仅在实际发现悬空
    # 时才写入(无悬空走提前 return,键不存在),断言改 .get 兜底 0
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t398\t.\t+\t.\tParent=t1\n")   # 相位缺失,总长299非3倍
    result = fixer.run()
    assert result.stats.get('dangling_mrna', 0) == 0
    assert 't1' in result.problematic_transcripts       # 299%3!=0


def test_phase_corrected_in_run(tmp_path):
    # brief 笔误修正|Brief typo fix:原第二段 CDS 500-599(100bp)+100-399
    # (300bp)=400 非 3 倍,会进 INVALID_LENGTH 而非 corrected;改 500-598
    # (99bp) 总长 399 为 3 倍数
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t399\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t500\t598\t.\t+\t.\tParent=t1\n")  # 总长399为3倍但phase缺失
    result = fixer.run()
    assert result.corrected_transcripts == ['t1']
    cds = [r for r in result.doc.records if r.type == 'CDS']
    assert [r.phase for r in sorted(cds, key=lambda r: r.start)] == ['0', '0']


def test_write_outputs_extra_files_conditional(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t399\t.\t+\t0\tParent=t1\n")  # 一切正常
    result = fixer.run()
    fixed = tmp_path / 'a.fixed.gff3'
    written = fixer.write_outputs(
        str(fixed),
        str(tmp_path / 'a.dangling.gff3'),
        str(tmp_path / 'a_phase_corrected.gff3'),
        str(tmp_path / 'a_phase_problematic.gff3'))
    assert fixed.exists()
    assert written == {'fixed': True, 'dangling': False,
                       'corrected': False, 'problematic': False}
    assert not (tmp_path / 'a.dangling.gff3').exists()


def test_write_outputs_corrected_contains_parent_gene(tmp_path):
    # brief 笔误修正|Brief typo fix:同 test_phase_corrected_in_run,第二段
    # CDS 599→598 使总长 399 为 3 倍数(否则进 problematic,corrected 文件不落盘)
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t600\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t600\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t399\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t500\t598\t.\t+\t.\tParent=t1\n")
    result = fixer.run()
    fixer.write_outputs(
        str(tmp_path / 'a.fixed.gff3'),
        str(tmp_path / 'a.dangling.gff3'),
        str(tmp_path / 'a_phase_corrected.gff3'),
        str(tmp_path / 'a_phase_problematic.gff3'))
    text = (tmp_path / 'a_phase_corrected.gff3').read_text()
    assert 'ID=g1' in text          # 含父 gene
    assert '\tCDS\t100\t399\t.\t+\t0\t' in text   # phase 已被修正写入


def test_dangling_skips_idless_transcript_row(tmp_path):
    fixer = make_fixer(tmp_path,
        "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tParent=g1\n"
        "chr1\t.\tgene\t700\t900\t.\t+\t.\tID=g2\n"
        "chr1\t.\tmRNA\t700\t800\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\t.\tCDS\t700\t795\t.\t+\t0\tParent=t2\n")
    dangling = fixer.split_dangling_mrna()
    assert dangling == []   # 无 ID mRNA 行跳过悬空判定,不虚报


def test_run_phase_minus_strand_applies_in_chain_order(tmp_path):
    # -链:CDS 文件序升序,链序降序;缺失 phase 须按链序修正
    # dispatch 笔误修正|Dispatch typo fix:原 100-399(300bp)+500-599(100bp)
    # =400 非 3 倍会进 INVALID_LENGTH;注释算术(cumLen=100→'2',总长399)要求
    # 第2段 CDS 长 299bp → 改 100-398
    fixer = make_fixer(tmp_path,
        "chr1\t.\tmRNA\t100\t600\t.\t-\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t100\t398\t.\t-\t.\tParent=t1\n"    # 链序第2条
        "chr1\t.\tCDS\t500\t599\t.\t-\t.\tParent=t1\n")   # 链序第1条
    fixer.run_phase_correction()
    cds_by_start = {r.start: r for r in fixer.doc.records
                    if r.type == 'CDS'}
    # 链序第1条(500-599,len100) phase='0';第2条 cumLen=100→100%3=1→'2'
    assert cds_by_start[500].phase == '0'
    assert cds_by_start[100].phase == '2'
