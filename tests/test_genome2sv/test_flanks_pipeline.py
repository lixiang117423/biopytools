"""genome2sv 步骤6 flank 侧翼序列测试|genome2sv step-6 flank pipeline tests"""

import logging

from biopytools.genome2sv.config import Genome2SVConfig
from biopytools.genome2sv.pipeline import Genome2SVPipeline

MERGED_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
c1\t200\t.\tA\tATTTTTTTTTT\t.\tPASS\tSVTYPE=INS;SVLEN=10;END=200\tGT:PSV:LN:QV\t1/1:0:10:60\t0/0:0:.:.
c1\t350\t.\tN\tN[ c1:351[\t.\tPASS\tSVTYPE=TRA;SVLEN=0;END=350;CHR2=c1\tGT:PSV:LN:QV\t0/1:0:0:60\t0/1:0:0:60
c1\t300\t.\tTTTTTTTTTT\tT\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=309\tGT:PSV:LN:QV\t./.:.:.:.\t1/1:0:10:50
"""


def make_pipeline(tmp_path, flank=50):
    """参考(fai 手写)+merged VCF 的最小流程|minimal pipeline with ref+fai+merged VCF"""
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">c1\n" + "ACGT" * 100 + "\n")
    q_fa = tmp_path / "q.fa"
    q_fa.write_text(">c1\n" + "ACGT" * 100 + "\n")
    fof = tmp_path / "samples.fof"
    fof.write_text(f"REF\t{ref_fa}\nQ1\t{q_fa}\n")
    cfg = Genome2SVConfig(input_fof=str(fof), ref_sample="REF",
                          output_dir=str(tmp_path / "out"), flank=flank)
    # pipeline 读 reference/ 下的参考副本:补 fasta + fai(单行 fasta:
    # offset=4(>c1 共 4 字节),linebases=400,linewidth=401)
    ref_copy = cfg.reference_dir / "ref.fa"
    ref_copy.write_text(">c1\n" + "ACGT" * 100 + "\n")
    ref_copy.with_suffix(".fa.fai").write_text("c1\t400\t4\t400\t401\n")
    vcf = cfg.merged_dir / "pan_sv_survivor.vcf"
    vcf.write_text(MERGED_VCF)
    pipeline = Genome2SVPipeline(cfg, logging.getLogger("genome2sv"))
    return cfg, pipeline, str(vcf)


class TestGenerateSvFlanks:
    def test_ins_end_artifact_forced_to_pos(self, tmp_path):
        """SURVIVOR 合并 VCF 的 INS END 可能带伪偏移(实测 POS+65536),
        插入位点以 POS 为准|bogus INS END (POS+65536 artifact) is forced to POS"""
        cfg, pipeline, _ = make_pipeline(tmp_path)
        vcf = cfg.merged_dir / "pan_sv_survivor.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "c1\t200\t.\tA\tATTTTTTTTT\t.\tPASS\t"
            "SVTYPE=INS;SVLEN=10;END=100000\tGT:LN:QV\t1/1:10:60\n")
        assert pipeline.generate_sv_flanks(str(vcf))
        fa = (cfg.flank_dir / "sv_flank50bp.fa").read_text()
        header = fa.splitlines()[0]
        assert "flank=150-250" in header, header
        seq = fa.split("\n", 1)[1].replace("\n", "")
        assert seq == ("ACGT" * 100)[149:250]

    def test_fasta_records_skip_tra(self, tmp_path):
        """TRA 跳过,INS/DEL 出序列且 header 完整|TRA skipped; INS/DEL emitted"""
        cfg, pipeline, vcf = make_pipeline(tmp_path)
        assert pipeline.generate_sv_flanks(vcf)
        fa = (cfg.flank_dir / "sv_flank50bp.fa").read_text()
        records = {}
        for block in fa.split(">")[1:]:
            header, _, seq = block.partition("\n")
            # FASTA 60 列换行,比对前去掉换行|FASTA wrapped at 60 cols
            records[header.split()[0]] = (header, seq.replace("\n", "").strip())
        assert set(records) == {"pan_sv.INS.00001", "pan_sv.DEL.00001"}
        ref = "ACGT" * 100
        ins_header, ins_seq = records["pan_sv.INS.00001"]
        assert ins_header == ("pan_sv.INS.00001 type=INS chrom=c1 pos=200 "
                              "end=200 svlen=10 flank=150-250 len=101 "
                              "samples=S1")
        assert ins_seq == ref[149:250]
        del_header, del_seq = records["pan_sv.DEL.00001"]
        assert del_header == ("pan_sv.DEL.00001 type=DEL chrom=c1 pos=300 "
                              "end=309 svlen=-10 flank=250-359 len=110 "
                              "samples=S2")
        assert del_seq == ref[249:359]

    def test_tsv_genotypes(self, tmp_path):
        """TSV 按样本出 GT/LN/QV|TSV carries per-sample GT/LN/QV"""
        cfg, pipeline, vcf = make_pipeline(tmp_path)
        pipeline.generate_sv_flanks(vcf)
        lines = (cfg.flank_dir / "sv_flank50bp.tsv").read_text().splitlines()
        assert lines[0].split("\t") == [
            "sv_id", "chrom", "pos", "end", "svtype", "svlen",
            "flank_lo", "flank_hi",
            "S1_GT", "S1_LN", "S1_QV", "S2_GT", "S2_LN", "S2_QV", "sequence"]
        del_cells = lines[2].split("\t")
        assert del_cells[:3] == ["pan_sv.DEL.00001", "c1", "300"]
        assert del_cells[8:14] == ["./.", ".", ".", "1/1", "10", "50"]

    def test_xlsx_written(self, tmp_path):
        cfg, pipeline, vcf = make_pipeline(tmp_path)
        pipeline.generate_sv_flanks(vcf)
        assert (cfg.flank_dir / "sv_flank50bp.xlsx").exists()

    def test_checkpoint_skips(self, tmp_path):
        """三输出齐全时跳过,不覆盖既有文件|all outputs present -> skip untouched"""
        cfg, pipeline, vcf = make_pipeline(tmp_path)
        sentinel = {}
        for name in ("sv_flank50bp.fa", "sv_flank50bp.tsv", "sv_flank50bp.xlsx"):
            p = cfg.flank_dir / name
            p.write_text("sentinel")
            sentinel[name] = p
        assert pipeline.generate_sv_flanks(vcf)
        assert all(p.read_text() == "sentinel" for p in sentinel.values())

    def test_missing_vcf_graceful(self, tmp_path):
        """merged VCF 缺失优雅跳过|missing merged VCF degrades gracefully"""
        cfg, pipeline, _ = make_pipeline(tmp_path)
        assert pipeline.generate_sv_flanks(
            str(cfg.merged_dir / "no_such.vcf"))
        assert not any(cfg.flank_dir.iterdir())
