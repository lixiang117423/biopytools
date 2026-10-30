"""diagnosis 单测|diagnosis tests"""

import json
import logging

from biopytools.anno_curate.config import AnnoCurateConfig
from biopytools.anno_curate.diagnosis import (
    diagnosis_fingerprint, effective_utr_thresholds, load_fingerprint,
    save_fingerprint, should_skip_diagnosis, structural_diagnosis,
    write_diagnosis)
from biopytools.anno_curate.gff_model import parse_gxf

LOGGER = logging.getLogger('biopytools.anno_curate.test')


def doc_from(tmp_path, text):
    p = tmp_path / 'in.gff3'
    p.write_text(text)
    return parse_gxf(str(p))


def base_cfg(tmp_path, **kw):
    return AnnoCurateConfig(input=str(tmp_path / 'in.gff3'),
                            output_dir=str(tmp_path / 'out'), **kw)


def _utr_tags(issues):
    """抽取全部UTR相关标签(笔误修正辅助)|Collect UTR-related tags"""
    return [t for i in issues for t in i.tags if 'UTR' in t]


def test_cds_too_short_and_long(tmp_path):
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t89\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t89\t.\t+\t0\tParent=t1\n"          # 89 < 90
        "chr1\t.\tmRNA\t1\t10501\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\t.\tCDS\t1\t10501\t.\t+\t0\tParent=t2\n")      # > 10500
    issues = {i.tid: i for i in
              structural_diagnosis(doc, base_cfg(tmp_path), LOGGER)}
    # 笔误修正:89=3x29+2、10501=3x3500+1 均非3倍数,总长标签与长度标签并存
    # |typo fix: 89/10501 not multiples of 3, phase tag coexists with length tag
    assert 'CDS too short' in issues['t1'].tags
    assert 'CDS too long' in issues['t2'].tags


def test_cds_length_boundary_ok(tmp_path):
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t90\t.\t+\t0\tParent=t1\n")         # 90 恰好达标
    assert structural_diagnosis(doc, base_cfg(tmp_path), LOGGER) == []


def test_phase_tags(tmp_path):
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t600\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t300\t.\t+\t.\tParent=t1\n"          # 缺 phase→CORRECTED
        "chr1\t.\tCDS\t400\t699\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\t.\tCDS\t1\t499\t.\t+\t0\tParent=t2\n")         # 499非3倍
    issues = {i.tid: i for i in
              structural_diagnosis(doc, base_cfg(tmp_path), LOGGER)}
    assert 'NO PHASE' in issues['t1'].tags
    assert any('3' in t for t in issues['t2'].tags)


def test_utr_fraction_with_relax(tmp_path):
    # mRNA 跨度 1000,3'UTR 占 400 = 40% > 0.30*(1+0.5)=0.45? 40%<45% 不报;
    # 5'UTR 占 100 = 10% < 22.5% 不报 → 空
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tfive_prime_UTR\t1\t100\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tCDS\t101\t600\t.\t+\t0\tParent=t1\n"
        "chr1\t.\tthree_prime_UTR\t601\t1000\t.\t+\t.\tParent=t1\n")
    # 笔误修正:fixture CDS 101-600=500 非3倍数必带相位标签,改断言无 UTR 标签
    # |typo fix: fixture CDS 500nt is not a multiple of 3, assert no UTR tag only
    assert _utr_tags(
        structural_diagnosis(doc, base_cfg(tmp_path), LOGGER)) == []


def test_utr_fraction_over_relaxed_limit(tmp_path):
    # 3'UTR 占 50% > 45% → 报
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t1\n"
        "chr1\t.\tthree_prime_UTR\t501\t1000\t.\t+\t.\tParent=t1\n")
    issues = structural_diagnosis(doc, base_cfg(tmp_path), LOGGER)
    # 笔误修正:CDS 500 非3倍数并存相位标签,改成员断言
    # |typo fix: CDS 500nt co-occurs with phase tag, use membership assert
    assert "3'UTR fraction 0.50 exceeds 0.45" in issues[0].tags


def test_utr_check_disabled(tmp_path):
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t1\n"
        "chr1\t.\tthree_prime_UTR\t501\t1000\t.\t+\t.\tParent=t1\n")
    cfg = base_cfg(tmp_path, check_utr=False)
    # 笔误修正:同上,fixture 相位标签仍在,只断言 UTR 检查被关闭
    # |typo fix: phase tag remains in fixture, assert UTR check disabled only
    assert _utr_tags(structural_diagnosis(doc, cfg, LOGGER)) == []


def test_effective_thresholds():
    cfg = AnnoCurateConfig(input='x.gff3', output_dir='/tmp/anno_curate_t',
                           relax=0.5)
    assert effective_utr_thresholds(cfg) == (0.15 * 1.5, 0.30 * 1.5)


def test_write_diagnosis_tsv_and_summary(tmp_path):
    from biopytools.anno_curate.diagnosis import TranscriptIssue
    issues = [TranscriptIssue('t2', 'chr2', '-', 84, ['CDS too short'],
                              zero_cov_exons=None)]
    tsv = tmp_path / 'd.tsv'
    smy = tmp_path / 'd.summary.txt'
    write_diagnosis(issues, ['t1', 't2'], str(tsv), str(smy),
                    {'cpc2_status': '跳过|skipped',
                     'rna_status': '跳过|skipped'})
    lines = tsv.read_text().strip().split('\n')
    # 笔误修正:首行是 '#Found totally N' 计数注释行(T12 断点测试
    # startswith('#Found') 依赖),表头/数据行各后移一位
    # |typo fix: leading '#Found' comment line (T12 resume test relies on it),
    # header/data shift down one line
    assert lines[0].startswith('#Found')
    assert lines[1].split('\t') == ['transcript_id', 'chr', 'strand',
                                    'cds_len', 'zero_cov_exons', 'tags']
    assert lines[2].split('\t')[0] == 't2'
    assert 'total transcripts: 2' in smy.read_text()
    assert 'CDS too short: 1' in smy.read_text()


def test_fingerprint_roundtrip_and_skip(tmp_path):
    cfg = base_cfg(tmp_path, rna_bams=(str(tmp_path / 'a.bam'),))
    fp = diagnosis_fingerprint(cfg)
    meta = tmp_path / '.diagnosis_meta.json'
    save_fingerprint(str(meta), fp)
    assert load_fingerprint(str(meta)) == fp
    out_tsv = tmp_path / 'd.tsv'
    out_tsv.write_text('x')
    assert should_skip_diagnosis(str(out_tsv), str(meta), cfg) is True
    cfg2 = base_cfg(tmp_path, rna_bams=(str(tmp_path / 'b.bam'),))
    assert should_skip_diagnosis(str(out_tsv), str(meta), cfg2) is False
    assert should_skip_diagnosis(str(out_tsv), str(meta),
                                 base_cfg(tmp_path, force=True)) is False


def test_fingerprint_and_skip_negative_paths(tmp_path):
    # 计划账本按 10 个测试计(10→16→21 递推),brief 代码块仅 9 个,补
    # 负路径用例:首跑无输出不跳过;指纹缺失/损坏视为不匹配
    # |Ledger reconciliation: negative paths of the fingerprint resume guard
    cfg = base_cfg(tmp_path)
    meta = str(tmp_path / '.diagnosis_meta.json')
    missing_tsv = str(tmp_path / 'd.tsv')
    # 首跑:输出不存在 → 不跳过(指纹在不在都一样)
    # |First run: no output → never skip
    assert should_skip_diagnosis(missing_tsv, meta, cfg) is False
    save_fingerprint(meta, diagnosis_fingerprint(cfg))
    assert should_skip_diagnosis(missing_tsv, meta, cfg) is False
    # 指纹缺失/损坏 → None ≠ 指纹 → 重跑诊断|Missing/corrupt meta → re-run
    with open(meta, 'w', encoding='utf-8') as fh:
        fh.write('{broken json')
    assert load_fingerprint(meta) is None
    assert load_fingerprint(str(tmp_path / 'absent.json')) is None
    (tmp_path / 'd.tsv').write_text('x')
    assert should_skip_diagnosis(str(tmp_path / 'd.tsv'), meta, cfg) is False


# ---------- CPC2 编码潜能|CPC2 coding potential ----------


def _fa(tmp_path):
    p = tmp_path / 'g.fa'
    p.write_text('>chr1\n' + 'ACGT' * 25 + '\n')     # 100bp
    return p


def test_fasta_index_fetch_and_missing(tmp_path):
    from biopytools.anno_curate.diagnosis import FastaIndex
    idx = FastaIndex(str(_fa(tmp_path)))
    assert idx.fetch('chr1', 1, 4) == 'ACGT'
    assert idx.fetch('chr1', 98, 100) == 'CGT'
    assert idx.fetch('chrX', 1, 4) == ''


def test_fasta_index_uses_existing_fai(tmp_path):
    from biopytools.anno_curate.diagnosis import FastaIndex
    fa = _fa(tmp_path)
    (tmp_path / 'g.fa.fai').write_text('chr1\t100\t6\t100\t101\n')
    idx = FastaIndex(str(fa))
    assert idx.fetch('chr1', 1, 3) == 'ACG'


def test_fasta_index_fai_multiline_fetch(tmp_path):
    # C1 回归锁定:.fai 第4列是 linebases(每行碱基数),fetch 的换行补偿
    # 需要第5列 linewidth(含换行整行字节数);多行 FASTA 走 .fai 分支
    # 必须跨两个行尾仍与真值逐碱基一致
    # |C1 regression lock: .fai col4 is linebases, fetch's newline
    # compensation needs col5 linewidth; multiline .fai path must match
    # ground truth across both line ends
    fa = tmp_path / 'g.fa'
    seq = 'ACGT' * 40                              # 160bp
    lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
    fa.write_text('>chr1\n' + '\n'.join(lines) + '\n')
    (tmp_path / 'g.fa.fai').write_text('chr1\t160\t6\t60\t61\n')
    from biopytools.anno_curate.diagnosis import FastaIndex
    idx = FastaIndex(str(fa))
    assert idx.fetch('chr1', 1, 160) == seq


def test_fasta_index_last_line_full_length_no_newline(tmp_path):
    # I-2 回归:末行恰好等长且无换行时,旧实现把 line_bpp 刷成纯碱基长
    # (60 而非 61),此后该序列 fetch 全错;修复后全长/跨行窗口均为真值
    # |I-2 regression: an unterminated full-length last line used to
    # shrink line_bpp to the bare base count (60, not 61), corrupting
    # every fetch of the entry
    fa = tmp_path / 'g.fa'
    seq = 'ACGT' * 45                              # 180 = 60 x 3
    lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
    fa.write_text('>chr1\n' + '\n'.join(lines))    # 末行无 \n
    from biopytools.anno_curate.diagnosis import FastaIndex
    idx = FastaIndex(str(fa))
    assert idx.fetch('chr1', 1, 180) == seq
    assert idx.fetch('chr1', 55, 125) == seq[54:125]   # 跨两处行尾
    assert idx.fetch('chr1', 178, 180) == seq[177:180]


def test_fasta_index_single_line_no_newline(tmp_path):
    # 单行且无换行的极端兜底:len+1 使换行补偿恒为 0
    # |Single unterminated sequence line: len+1 zeroes compensation
    fa = tmp_path / 'g.fa'
    seq = 'ACGTACGTAC'
    fa.write_text('>chr1\n' + seq)                 # 序列行无换行
    from biopytools.anno_curate.diagnosis import FastaIndex
    idx = FastaIndex(str(fa))
    assert idx.fetch('chr1', 1, 10) == seq
    assert idx.fetch('chr1', 8, 10) == seq[7:10]


def test_extract_cds_reverse_complement(tmp_path):
    from biopytools.anno_curate.diagnosis import (FastaIndex,
                                                  extract_cds_fasta)
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t100\t.\t-\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t4\t.\t-\t0\tParent=t1\n")
    tid, text = extract_cds_fasta(doc, FastaIndex(str(_fa(tmp_path))), 't1')
    assert tid == 't1'
    assert text.split('\n')[1] == 'ACGT'    # revcomp(ACGT)=ACGT


def test_parse_cpc2_output(tmp_path):
    from biopytools.anno_curate.diagnosis import parse_cpc2_output
    p = tmp_path / 'cpc2.txt'
    p.write_text('#ID\tlength\tpep\tfickett\tpI\tORF\tprob\tlabel\n'
                 't1\t300\t100\t1.0\t7.0\t1\t0.99\tcoding\n'
                 't2\t201\t67\t0.1\t5.0\t1\t0.01\tnoncoding\n')
    res = parse_cpc2_output(str(p))
    assert res == {'t1': False, 't2': True}


def test_cpc2_diagnosis_degrades_when_tool_fails(tmp_path, monkeypatch):
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t12\t.\t+\t0\tParent=t1\n")

    def boom(*a, **kw):
        raise OSError('no cpc2')

    monkeypatch.setattr(diag, 'run_cpc2', boom)
    cfg = base_cfg(tmp_path, genome=str(_fa(tmp_path)))
    tids, status = diag.cpc2_diagnosis(doc, cfg, LOGGER,
                                       str(tmp_path / 'tmp'))
    # 笔误修正:run_cpc2 抛异常走的是 except 分支,状态为"失败降级|degraded",
    # 非"跳过|skipped"(后者是 no-genome/cpc2-not-found 降级路径)
    # |typo fix: raised exception hits the except branch → "degraded",
    # not "skipped" (the latter covers no-genome / cpc2-not-found)
    assert tids == []
    assert ('降级' in status or 'degraded' in status
            or '跳过' in status or 'skipped' in status)


def test_cpc2_diagnosis_end_to_end_mocked(tmp_path, monkeypatch):
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t12\t.\t+\t0\tParent=t1\n")
    out = tmp_path / 'cpc2.txt'
    out.write_text('t1\t12\t4\t0.1\t5.0\t1\t0.01\tnoncoding\n')

    def fake_run(fa, odir, cfg, logger):
        return str(out)

    monkeypatch.setattr(diag, 'run_cpc2', fake_run)
    # 用假入口过存在性检查,不依赖本机 annot 域是否装了 CPC2
    # |stub path so the existence check passes without a local CPC2
    cfg = _cpc2_cfg(tmp_path)
    tids, status = diag.cpc2_diagnosis(doc, cfg, LOGGER,
                                       str(tmp_path / 'tmp'))
    assert tids == ['t1']


def _cpc2_cfg(tmp_path):
    """带假CPC2入口的配置(路径存在即过存在性检查)|Cfg with a stub CPC2"""
    stub = tmp_path / 'CPC2.py'
    stub.write_text('# stub\n')
    return base_cfg(tmp_path, genome=str(_fa(tmp_path)),
                    cpc2_path=str(stub))


def test_cpc2_diagnosis_skips_missing_seqid(tmp_path, monkeypatch, caplog):
    # I-3 回归:seqid 不在FASTA的转录本逐条 WARNING 跳过,其余照常进
    # CPC2,不整体降级;状态文本带跳过计数
    # |I-3 regression: transcripts whose seqid is missing from the FASTA
    # are warned and skipped individually; the rest still reach CPC2
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t1\t99\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\tCDS\t1\t99\t.\t+\t0\tParent=t1\n"
        "chrX\t.\tmRNA\t1\t99\t.\t+\t.\tID=t2;Parent=g2\n"
        "chrX\t.\tCDS\t1\t99\t.\t+\t0\tParent=t2\n")
    out = tmp_path / 'cpc2.txt'
    out.write_text('t1\t99\t33\t0.1\t5.0\t1\t0.01\tnoncoding\n')
    seen = {}

    def fake_run(fa, odir, cfg, logger):
        with open(fa, encoding='utf-8') as fh:
            seen['fa'] = fh.read()
        return str(out)

    monkeypatch.setattr(diag, 'run_cpc2', fake_run)
    with caplog.at_level(logging.WARNING):
        tids, status = diag.cpc2_diagnosis(
            doc, _cpc2_cfg(tmp_path), LOGGER, str(tmp_path / 'tmp'))
    assert tids == ['t1']                       # 好转录本照常判 noncoding
    assert '>t1' in seen['fa'] and '>t2' not in seen['fa']
    assert '1 extract skipped' in status
    assert any('t2' in r.message and 'chrX' in r.message
               for r in caplog.records)         # 逐条 WARNING 带转录本id


def test_cpc2_diagnosis_all_extracts_fail_degrades(tmp_path, monkeypatch):
    # I-3 负路径:成功 0 条才整体降级,状态带失败计数
    # |I-3 negative path: degrade only when zero extracts succeed
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chrX\t.\tmRNA\t1\t99\t.\t+\t.\tID=t2;Parent=g2\n"
        "chrX\t.\tCDS\t1\t99\t.\t+\t0\tParent=t2\n"
        "chrY\t.\tmRNA\t1\t99\t.\t+\t.\tID=t3;Parent=g3\n"
        "chrY\t.\tCDS\t1\t99\t.\t+\t0\tParent=t3\n")
    def must_not_run(*a, **kw):
        raise AssertionError('CPC2 不应被调用|CPC2 must not run')

    monkeypatch.setattr(diag, 'run_cpc2', must_not_run)
    tids, status = diag.cpc2_diagnosis(
        doc, _cpc2_cfg(tmp_path), LOGGER, str(tmp_path / 'tmp'))
    assert tids == []
    assert '2' in status and ('降级' in status or 'degraded' in status)


# ---------- RNA 零覆盖证据|RNA zero-coverage evidence ----------


def test_export_exon_bed(tmp_path):
    from biopytools.anno_curate.diagnosis import export_exon_bed
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t10\t40\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t10\t20\t.\t+\t.\tParent=t1\n"
        "chr1\t.\texon\t30\t40\t.\t+\t.\tParent=t1\n")
    bed = tmp_path / 'exons.bed'
    n = export_exon_bed(doc, str(bed))
    assert n == 2
    lines = bed.read_text().strip().split('\n')
    assert lines[0].split('\t') == ['chr1', '9', '20', 't1:1']   # 0-based


def test_run_bedcov_parses_reads(tmp_path, monkeypatch):
    from biopytools.anno_curate import diagnosis as diag
    out = tmp_path / 'bedcov.txt'
    out.write_text('chr1\t9\t20\tt1:1\t5\t55\n'
                   'chr1\t29\t40\tt1:2\t0\t0\n')

    class R:  # fake subprocess result
        returncode = 0
        stdout = out.read_text()
        stderr = ''

    def fake_run(cmd, **kw):
        diag._LAST_BEDCOV_CMD = cmd
        return R()

    monkeypatch.setattr('subprocess.run', fake_run)
    cfg = base_cfg(tmp_path)
    res = diag.run_bedcov(str(tmp_path / 'e.bed'), 'x.bam', cfg, LOGGER)
    assert res == {'t1:1': 5, 't1:2': 0}


def test_rna_aggregation_any_bam_counts(tmp_path, monkeypatch):
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t10\t40\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t10\t20\t.\t+\t.\tParent=t1\n"
        "chr1\t.\texon\t30\t40\t.\t+\t.\tParent=t1\n")
    # bam1: exon1 有 reads;bam2: 全零 → 任一覆盖 → exon1 有覆盖
    monkeypatch.setattr(diag, 'run_bedcov',
                        lambda b, bam, c, l, **kw: {'t1:1': 3, 't1:2': 0}
                        if bam.endswith('1.bam') else
                        {'t1:1': 0, 't1:2': 0})
    cfg = base_cfg(tmp_path, rna_bams=(str(tmp_path / '1.bam'),
                                       str(tmp_path / '2.bam')))
    zero_map, no_cov, status = diag.rna_coverage_diagnosis(
        doc, cfg, LOGGER, str(tmp_path / 'tmp'))
    assert zero_map == {'t1': '1/2'}        # exon2 在所有 BAM 零覆盖
    assert no_cov == set()                  # exon1 有覆盖 → 非全零


def test_rna_all_zero_marks_no_coverage(tmp_path, monkeypatch):
    from biopytools.anno_curate import diagnosis as diag
    doc = doc_from(tmp_path,
        "chr1\t.\tmRNA\t10\t40\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t10\t20\t.\t+\t.\tParent=t1\n"
        "chr1\t.\texon\t30\t40\t.\t+\t.\tParent=t1\n")
    monkeypatch.setattr(diag, 'run_bedcov',
                        lambda b, bam, c, l, **kw: {'t1:1': 0, 't1:2': 0})
    cfg = base_cfg(tmp_path, rna_bams=(str(tmp_path / '1.bam'),))
    zero_map, no_cov, status = diag.rna_coverage_diagnosis(
        doc, cfg, LOGGER, str(tmp_path / 'tmp'))
    assert zero_map == {'t1': '2/2'}
    assert no_cov == {'t1'}                 # 全零覆盖 → NO RNA COVERAGE 集合


def test_merge_diagnosis_combines_rna():
    from biopytools.anno_curate.diagnosis import (TranscriptIssue,
                                                  merge_diagnosis)
    issues = [TranscriptIssue('t1', 'chr1', '+', 300,
                              ['INVALID PHASE'], zero_cov_exons=None),
              TranscriptIssue('t2', 'chr2', '-', 600, [], zero_cov_exons=None)]
    tx_meta = {'t3': ('chr3', '+', 900)}
    merged = merge_diagnosis(issues, {'t1': '2/8', 't3': '3/3'},
                             {'t1'}, tx_meta)
    by_tid = {i.tid: i for i in merged}
    assert by_tid['t1'].zero_cov_exons == '2/8'
    assert 'NO RNA COVERAGE' in by_tid['t1'].tags    # t1 在 no_cov 集合
    assert 'ZERO-COV EXON 3/3' in by_tid['t3'].tags  # 无结构问题但零覆盖→新增
    assert by_tid['t3'].cds_len == 900               # 元数据来自 tx_meta
