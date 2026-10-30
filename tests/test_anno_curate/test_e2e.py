"""anno_curate e2e(计算节点)|e2e (compute node)

真实 samtools(build/index/bedcov)+ 真实 CPC2(-g 触发,annot 域),
完整 CLI 子进程跑通并断言输出结构与内容。

跑法|Run: bsub -q c02 -K -o <log> sh -c \
    'cd <repo_root> && python -m pytest \
     tests/test_anno_curate/test_e2e.py -m e2e -v'
"""

import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.e2e

REPO_ROOT = Path(__file__).resolve().parents[2]

# 笔误核查记录(先例14处,本文件逐条手工推演过):
# - t1: CDS 200-600(401nt)+401-500(100nt)=501=3x167 → 3倍数,相位可修;
#   期望相位 0/1(cum=401,401%3=2,phase=3-2=1)
# - t2: CDS 100-399=300nt → 3倍数+phase0 合法 → 干净转录本不入诊断表
# - t3: 无CDS,exon e9 重复x2 → 子ID去重为 e9/e9.u1,转录本级无诊断
# - t4: CDS 1-80=80nt → 非3倍数(80%3=2)且<90 → 双tag 有意构造
BAD_GFF = """##gff-version 3
chr1\t.\tgene\t600\t100\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t600\t100\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t600\t200\t.\t+\t.\tParent=t1
chr1\t.\tCDS\t500\t401\t.\t+\t.\tParent=t1
chr2\t.\tmRNA\t50\t400\t.\t-\t.\tID=t2
chr2\t.\texon\t50\t400\t.\t-\t.\tParent=t2
chr2\t.\tCDS\t100\t399\t.\t-\t0\tParent=t2
chr3\t.\tgene\t10\t99\t.\t+\t.\tID=g3
chr3\t.\tmRNA\t10\t99\t.\t+\t.\tID=t3;Parent=g3
chr3\t.\texon\t10\t99\t.\t+\t.\tID=e9;Parent=t3
chr3\t.\texon\t10\t99\t.\t+\t.\tID=e9;Parent=t3
chr4\t.\tmRNA\t1\t80\t.\t+\t.\tID=t4;Parent=g4
chr4\t.\tCDS\t1\t80\t.\t+\t0\tParent=t4
"""

# CPC2 场景:2 转录本均有 CDS(600nt=3倍数,相位合法,无结构问题),
# tB 负链双外显子覆盖 revcomp+多段拼接路径
CPC2_GFF = """##gff-version 3
chr1\t.\tgene\t101\t700\t.\t+\t.\tID=gA
chr1\t.\tmRNA\t101\t700\t.\t+\t.\tID=tA;Parent=gA
chr1\t.\texon\t101\t700\t.\t+\t.\tParent=tA
chr1\t.\tCDS\t101\t700\t.\t+\t0\tParent=tA
chr2\t.\tgene\t201\t1100\t.\t-\t.\tID=gB
chr2\t.\tmRNA\t201\t1100\t.\t-\t.\tID=tB;Parent=gB
chr2\t.\texon\t201\t500\t.\t-\t.\tParent=tB
chr2\t.\texon\t801\t1100\t.\t-\t.\tParent=tB
chr2\t.\tCDS\t801\t1100\t.\t-\t0\tParent=tB
chr2\t.\tCDS\t201\t500\t.\t-\t0\tParent=tB
"""

# GTF 场景:同基因 g1 两个可变剪接体(transcript 行自身带 gene_id+
# transcript_id);修复前 transcript 行被转成"Parent=t1 的无 ID 孤儿",
# ⑦ 会拆出 t1/t2 各自的假 gene 且转录本计数翻倍
# |GTF scenario: two isoforms of one gene; before the transcript-row fix
# each row became an ID-less orphan and ⑦ split them into separate genes
GTF_ISOFORMS = (
    'chr1\tsrc\ttranscript\t100\t399\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\texon\t100\t399\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\tCDS\t100\t399\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t1";\n'
    'chr1\tsrc\ttranscript\t500\t799\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t2";\n'
    'chr1\tsrc\texon\t500\t799\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t2";\n'
    'chr1\tsrc\tCDS\t500\t799\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t2";\n')


def _samtools_base():
    """定位真实samtools|Locate real samtools(PATH→align域→conda run)"""
    found = shutil.which('samtools')
    if found:
        return [found]
    align = Path('~/miniforge3/envs/align/bin/samtools').expanduser()
    if align.exists():
        return [str(align)]
    return ['conda', 'run', '-n', 'align', '--no-capture-output', 'samtools']


def _run_cli(*args, cwd):
    """跑完整CLI子进程|Run full CLI as subprocess

    cwd 指向临时目录防 CLI 向仓库写任何中间产物;biopytools 经
    PYTHONPATH 导入(包无 __main__.py,python -m biopytools 不可用,
    单模块直跑 main.py 末尾的 __main__ 块)。
    |cwd points to a temp dir so the CLI never writes into the repo;
    biopytools imports via PYTHONPATH (the package has no __main__.py,
    so we drive the __main__ block of anno_curate/main.py directly).
    """
    env = dict(os.environ)
    env['PYTHONPATH'] = (str(REPO_ROOT) + os.pathsep +
                         env.get('PYTHONPATH', ''))
    return subprocess.run(
        [sys.executable, '-m', 'biopytools.anno_curate.main', *args],
        shell=False, capture_output=True, text=True, cwd=str(cwd), env=env)


def _read_diagnosis_rows(tsv_text):
    """诊断TSV→{tid: cols}|Diagnosis TSV → {tid: columns}"""
    rows = {}
    for line in tsv_text.splitlines():
        if not line or line.startswith('#') or \
                line.startswith('transcript_id'):
            continue
        cols = line.split('\t')
        rows[cols[0]] = cols
    return rows


@pytest.fixture
def e2e_env(tmp_path):
    gff = tmp_path / 'e2e.gff3'
    gff.write_text(BAD_GFF)
    # 构造小RNA BAM: t2 外显子区间给几条read。
    # SAM 头必须声明全部 4 条染色体(brief 笔误修正:原 fixture 只有 chr2,
    # samtools bedcov 对头外染色体的 BED 行报错 rc=2 致 RNA 分支整体降级;
    # 真实同基因组 BAM 头必含全部 seq)|Header must declare all 4 chroms.
    sam = tmp_path / 'r.sam'
    sam.write_text(
        '@HD\tVN:1.6\tSO:coordinate\n'
        '@SQ\tSN:chr1\tLN:1000\n'
        '@SQ\tSN:chr2\tLN:1000\n'
        '@SQ\tSN:chr3\tLN:100\n'
        '@SQ\tSN:chr4\tLN:100\n'
        'r1\t0\tchr2\t60\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n'
        'r2\t0\tchr2\t300\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\t*\n')
    bam = tmp_path / 'rna.bam'
    subprocess.run(_samtools_base() + ['view', '-bS', '-o', str(bam),
                                       str(sam)],
                   check=True, shell=False)
    return gff, bam, tmp_path / 'out'


def test_e2e_full(e2e_env):
    gff, bam, out = e2e_env
    r = _run_cli('-i', str(gff), '-o', str(out), '--rna-bam', str(bam),
                 cwd=out.parent)
    assert r.returncode == 0, r.stderr
    fixed = (out / '01_fixed' / 'e2e.fixed.gff3').read_text()
    # 坐标已交换+相位已修+层级完整
    assert '\tCDS\t200\t600' in fixed or '\tCDS\t401\t600' in fixed
    assert '\tmRNA\t' in fixed and '\tgene\t' in fixed
    # 相位真修:t1 两段 CDS 期望相位 0/1(401nt 后 cum%3=2→3-2=1)
    assert '\tCDS\t200\t600\t.\t+\t0\t' in fixed
    assert '\tCDS\t401\t500\t.\t+\t1\t' in fixed
    # 层级补建:t2 无 Parent→t2.gene;t4 Parent 指向缺失 g4→补建 g4
    assert 'ID=t2.gene' in fixed and 'ID=g4' in fixed
    # 相位修复/问题转录本附加文件(t1 修相位/t4 非3倍数)
    assert (out / '01_fixed' / 'e2e_phase_corrected.gff3').exists()
    assert (out / '01_fixed' / 'e2e_phase_problematic.gff3').exists()
    tsv = (out / '02_diagnosis' / 'e2e.diagnosis.tsv').read_text()
    assert 't4' in tsv                 # 80nt CDS too short
    assert 'e9' not in tsv             # ID去重不出现在转录本级诊断
    assert '#Found totally 3 transcripts' in tsv
    assert (out / '02_diagnosis' /
            '.diagnosis_meta.json').exists()

    # ---------- bedcov 列序验证(T11 遗留)|bedcov cols[4]=reads ----------
    # 手工预期:BAM 仅 chr2 有 reads(r1:60-79,r2:300-319);
    # t2:1(chr2 49-400) reads>0 → t2 不入表;其余转录本外显子全零:
    # t1:1(chr1)/t3:1+t3:2(chr3)/t4:1(chr4) → 1/1、2/2、1/1。
    # 若 cols[4] 解析错列(如取到 BED 名列/其他列),t2 会误入零覆盖表。
    rows = _read_diagnosis_rows(tsv)
    assert set(rows) == {'t1', 't3', 't4'}, rows
    assert 't2' not in rows            # 有reads外显子未误判零覆盖
    assert rows['t1'][4] == '1/1', rows['t1']
    assert rows['t3'][4] == '2/2', rows['t3']
    assert rows['t4'][4] == '1/1', rows['t4']
    assert 'NO RNA COVERAGE' in rows['t4'][5]
    # t4 双结构 tag:CDS 总长非3倍数 + too short
    assert 'CDS too short' in rows['t4'][5]
    assert 'multiple of 3' in rows['t4'][5]
    smy = (out / '02_diagnosis' / 'e2e.diagnosis.summary.txt').read_text()
    assert 'NO RNA COVERAGE: 3' in smy
    assert 'RNA coverage: 完成|done (1 BAMs, 3 transcripts with ' \
           'zero-cov exons)' in smy
    assert 'CPC2: 跳过|skipped (no --genome)' in smy


def test_e2e_duplicate_id_hard_fail(tmp_path):
    gff = tmp_path / 'dup.gff3'
    gff.write_text("chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\t.\tmRNA\t20\t29\t.\t+\t.\tID=t1;Parent=g1\n")
    r = _run_cli('-i', str(gff), '-o', str(tmp_path / 'out'),
                 cwd=tmp_path)
    assert r.returncode == 1
    assert 't1' in (r.stdout + r.stderr)


def test_e2e_gtf_transcript_rows(tmp_path):
    """GTF transcript 行转换(I-1)|GTF transcript-row conversion (I-1)"""
    gtf = tmp_path / 'iso.gtf'
    gtf.write_text(GTF_ISOFORMS)
    out = tmp_path / 'out'
    r = _run_cli('-i', str(gtf), '-o', str(out), cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    fixed = (out / '01_fixed' / 'iso.fixed.gff3').read_text()
    # 恰好一条 gene,两个转录本级同挂 g1(可变剪接不拆家)
    # |Exactly one gene; both transcripts parented to g1
    rows = [ln.split('\t') for ln in fixed.splitlines()
            if '\t' in ln and not ln.startswith('#')]
    gene_lines = [ln for ln in rows if ln[2] == 'gene']
    assert len(gene_lines) == 1 and 'ID=g1' in gene_lines[0]
    tx_lines = [ln for ln in rows if ln[2] == 'transcript']
    assert len(tx_lines) == 2
    assert all('ID=t' in ln[8] and 'Parent=g1' in ln[8]
               for ln in tx_lines)
    assert {ln[8].split(';')[0] for ln in tx_lines} == {'ID=t1', 'ID=t2'}
    # 子特征挂各自转录本;无 t1/t2 命名的假 gene 残留
    # |Children parent to their own transcript; no fake tN genes
    assert 'Parent=t1' in fixed and 'Parent=t2' in fixed
    assert 'ID=t1.gene' not in fixed and 'ID=t2.gene' not in fixed
    # all_tids 不翻倍:恰好 2 个转录本级记录且全部无结构问题
    # |Roster not doubled: exactly 2 transcripts, all clean
    smy = (out / '02_diagnosis' / 'iso.diagnosis.summary.txt').read_text()
    assert 'total transcripts: 2' in smy
    assert 'transcripts with issues: 0' in smy
    assert (out / '01_fixed' / 'iso_phase_corrected.gff3').exists()


def _write_genome(path):
    """合成小基因组(60列折行,锻炼多行FASTA偏移索引)|
    Synthesize a small genome (60-col wrapped lines)"""
    rng = random.Random(42)
    with open(path, 'w', encoding='utf-8') as fh:
        for name, n in (('chr1', 2000), ('chr2', 1500)):
            fh.write(f'>{name}\n')
            s = ''.join(rng.choices('ACGT', k=n))
            for i in range(0, n, 60):
                fh.write(s[i:i + 60] + '\n')


def test_e2e_cpc2_genome(tmp_path):
    """CPC2 真跑(annot 域)|Real CPC2 run (annot domain env)"""
    gff = tmp_path / 'cpc2.gff3'
    gff.write_text(CPC2_GFF)
    genome = tmp_path / 'genome.fa'
    _write_genome(genome)
    out = tmp_path / 'out'
    r = _run_cli('-i', str(gff), '-o', str(out), '-g', str(genome),
                 cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    diag = out / '02_diagnosis'
    assert (diag / 'cpc2.diagnosis.tsv').exists()
    assert (diag / '.diagnosis_meta.json').exists()
    smy = (diag / 'cpc2.diagnosis.summary.txt').read_text()
    # CPC2 状态必须真跑完成;降级时把模块日志尾巴带进失败信息,不许静默
    # |CPC2 must actually run; on degradation surface the module log tail
    m = re.search(r'^CPC2: (.+)$', smy, re.MULTILINE)
    assert m, smy
    log_tail = ''
    log_path = out / '99_logs' / 'anno_curate.log'
    if log_path.exists():
        log_tail = log_path.read_text()[-2000:]
    assert m.group(1).startswith('完成|done'), \
        f"CPC2 未真跑完成|CPC2 not done: {m.group(1)}\n{log_tail}"
    # done 行带 noncoding N/2 计数(2 条 CDS 转录本进了 CPC2)
    assert re.search(r'CPC2: 完成\|done \(noncoding \d/2\)', smy), smy
    # 无 --rna-bam 时 RNA 分支明确跳过
    assert 'RNA coverage: 跳过|skipped (no --rna-bam)' in smy
