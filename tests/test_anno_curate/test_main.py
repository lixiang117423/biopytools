"""main 编排单测|main orchestration tests"""

import logging

import pytest

from biopytools.anno_curate.main import run_pipeline
from biopytools.anno_curate.config import AnnoCurateConfig


# 笔误修正(brief 第12处):第二条 CDS 500-599(100bp)使总长 400 非整除3,
# 相位走 INVALID_LENGTH 不写 0,test_full_run_layout 的 phase 断言必失败;
# 改 500-598(99bp)总长 399=3x133,CORRECTED 路径生效
# |Typo fix (12th in brief): second CDS 500-599 (100bp) gives total 400,
# not a multiple of 3, so phase stays '.' via INVALID_LENGTH and the phase
# assert below fails; 500-598 (99bp) totals 399=3x133, enabling CORRECTED
GOOD_GFF = """##gff-version 3
chr1\t.\tgene\t100\t600\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t600\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t100\t399\t.\t+\t.\tParent=t1
chr1\t.\tCDS\t500\t598\t.\t+\t.\tParent=t1
"""


@pytest.fixture
def setup(tmp_path):
    gff = tmp_path / 'in.gff3'
    gff.write_text(GOOD_GFF)
    out = tmp_path / 'out'
    cfg = AnnoCurateConfig(input=str(gff), output_dir=str(out))
    return cfg, out


def test_full_run_layout(setup):
    import logging
    cfg, out = setup
    assert run_pipeline(cfg, logging.getLogger('t')) == 0
    assert (out / '01_fixed' / 'in.fixed.gff3').exists()
    assert (out / '02_diagnosis' / 'in.diagnosis.tsv').exists()
    assert (out / '02_diagnosis' / 'in.diagnosis.summary.txt').exists()
    assert (out / '00_pipeline_info' / 'software_versions.yml').exists()
    assert (out / '99_logs').exists()
    fixed = (out / '01_fixed' / 'in.fixed.gff3').read_text()
    assert '\tCDS\t100\t399\t.\t+\t0\t' in fixed   # phase 已修


def test_two_level_resume(setup):
    import logging
    cfg, out = setup
    run_pipeline(cfg, logging.getLogger('t'))
    tsv = out / '02_diagnosis' / 'in.diagnosis.tsv'
    tsv.write_text('#edited\n')
    # 再跑:参数一致 → 诊断跳过,编辑保留
    run_pipeline(cfg, logging.getLogger('t'))
    assert tsv.read_text() == '#edited\n'
    # 换参数 → 只重跑诊断
    cfg2 = AnnoCurateConfig(input=cfg.input, output_dir=str(out),
                            relax=0.2)
    run_pipeline(cfg2, logging.getLogger('t'))
    assert tsv.read_text().startswith('#Found')


def test_skip_diagnosis(setup):
    import logging
    cfg, out = setup
    cfg.skip_diagnosis = True
    run_pipeline(cfg, logging.getLogger('t'))
    assert (out / '01_fixed' / 'in.fixed.gff3').exists()
    assert not (out / '02_diagnosis').exists()


# 相位齐全版(GOOD_GFF 的 CDS phase '.'→'0',相位 VALID,corrected 不落盘)
# |All-phases-present variant (phase '.'→'0' on CDS, VALID, no corrected)
CLEAN_GFF = GOOD_GFF.replace('\t.\tParent=t1', '\t0\tParent=t1')


def test_force_purges_stale_extras(tmp_path):
    # M-1 回归:第一轮相位缺失→corrected 落盘;人工修好输入后 --force
    # 重跑,本轮没有 corrected,旧文件必须清掉不残留
    # |M-1 regression: round 1 writes _phase_corrected; after fixing the
    # input, a --force rerun must not leave the stale copy behind
    import logging
    gff = tmp_path / 'in.gff3'
    gff.write_text(GOOD_GFF)
    out = tmp_path / 'out'
    assert run_pipeline(AnnoCurateConfig(input=str(gff),
                        output_dir=str(out)), logging.getLogger('t')) == 0
    corrected = out / '01_fixed' / 'in_phase_corrected.gff3'
    assert corrected.exists()
    gff.write_text(CLEAN_GFF)
    assert run_pipeline(AnnoCurateConfig(input=str(gff), force=True,
                        output_dir=str(out)), logging.getLogger('t')) == 0
    assert not corrected.exists()
    assert not (out / '01_fixed' / 'in_phase_problematic.gff3').exists()
    assert not (out / '01_fixed' / 'in.fixed.dangling.gff3').exists()
    assert (out / '01_fixed' / 'in.fixed.gff3').exists()
    # --force 连旧诊断目录一起清后重建
    # |--force purges and rebuilds the diagnosis dir too
    assert (out / '02_diagnosis' / 'in.diagnosis.tsv').exists()


def test_refix_without_force_clears_stale_extras(tmp_path):
    # M-1 附属:主输出缺失(上一轮中途崩溃)重跑修复时,条件附加文件同样
    # 以本轮判定为准,旧副本不冒充结果
    # |M-1 companion: re-fixing after a crashed round (main output
    # missing) also clears the conditional extras
    import logging
    gff = tmp_path / 'in.gff3'
    gff.write_text(GOOD_GFF)
    out = tmp_path / 'out'
    assert run_pipeline(AnnoCurateConfig(input=str(gff),
                        output_dir=str(out)), logging.getLogger('t')) == 0
    corrected = out / '01_fixed' / 'in_phase_corrected.gff3'
    assert corrected.exists()
    (out / '01_fixed' / 'in.fixed.gff3').unlink()   # 模拟中途崩溃
    gff.write_text(CLEAN_GFF)
    assert run_pipeline(AnnoCurateConfig(input=str(gff),
                        output_dir=str(out)), logging.getLogger('t')) == 0
    assert not corrected.exists()
    assert (out / '01_fixed' / 'in.fixed.gff3').exists()


def test_tmp_dir_cleaned_when_diagnosis_raises(tmp_path, monkeypatch):
    # M-2 回归:诊断分支异常也清 tmp(finally 兜底)
    # |M-2 regression: tmp cleaned even when diagnosis raises
    import logging
    from biopytools.anno_curate import main as m
    gff = tmp_path / 'in.gff3'
    gff.write_text(GOOD_GFF)
    out = tmp_path / 'out'
    cfg = AnnoCurateConfig(input=str(gff), output_dir=str(out))

    def boom(*a, **kw):
        raise RuntimeError('diag exploded')

    monkeypatch.setattr(m, 'structural_diagnosis', boom)
    with pytest.raises(RuntimeError):
        run_pipeline(cfg, logging.getLogger('t'))
    assert (out / '01_fixed' / 'in.fixed.gff3').exists()   # 修复已完成
    assert not (out / 'tmp').exists()                      # tmp 已清
