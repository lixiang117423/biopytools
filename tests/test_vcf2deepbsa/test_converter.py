"""vcf2deepbsa converter 测试|vcf2deepbsa converter tests

输出格式必须与原 deepbsa vcf2csv 逐字节一致:
无 header,逗号分隔,每行 CHROM,POS,REF,ALT,AD_ref_s1,AD_alt_s1,...
|Output must be byte-identical to the original deepbsa vcf2csv:
no header, comma-separated, each row CHROM,POS,REF,ALT,AD_ref_s1,AD_alt_s1,...
"""

from pathlib import Path

from biopytools.vcf2deepbsa.config import Vcf2DeepBsaConfig
from biopytools.vcf2deepbsa.converter import Vcf2DeepBsaConverter


def run_convert(tmp_path, vcf, logger, out_name="out"):
    """跑一次转换,返回 (csv路径, converter)|Run once, return (csv path, converter)"""
    cfg = Vcf2DeepBsaConfig(input_vcf=vcf, output_dir=str(tmp_path / out_name))
    conv = Vcf2DeepBsaConverter(cfg, logger)
    csv_path = Path(conv.run())
    return csv_path, conv


class TestBasicConversion:
    """基本转换|Basic conversion"""

    def test_csv_content_exact(self, tmp_path, make_vcf, logger):
        csv_path, _ = run_convert(tmp_path, make_vcf(), logger)
        assert csv_path.name == "x_deepbsa.csv"
        text = open(csv_path).read()
        assert text == (
            "chr1,100,A,T,10,5,20,8\n"
            "chr1,200,G,C,0,30,15,3\n"
            "chr2,300,AT,A,7,9,2,22\n"
        )

    def test_no_header_line(self, tmp_path, make_vcf, logger):
        csv_path, _ = run_convert(tmp_path, make_vcf(), logger)
        first_line = open(csv_path).readline()
        assert not first_line.startswith("CHROM")
        assert not first_line.startswith(",")

    def test_stats_reported(self, tmp_path, make_vcf, logger):
        _, conv = run_convert(tmp_path, make_vcf(), logger)
        assert conv.stats["total"] == 3
        assert conv.stats["written"] == 3
        assert conv.stats["sample_count"] == 2

    def test_gz_input_same_output(self, tmp_path, make_vcf, logger):
        csv_plain, _ = run_convert(tmp_path, make_vcf(name="a.vcf"), logger,
                                   out_name="o1")
        csv_gz, _ = run_convert(tmp_path, make_vcf(name="a.vcf.gz"), logger,
                                out_name="o2")
        assert open(csv_plain).read() == open(csv_gz).read()
        assert csv_gz.name == "a_deepbsa.csv"

    def test_multiallelic_takes_first_two_ad(self, tmp_path, make_vcf,
                                             logger):
        vcf = make_vcf(records=[
            ("chr1", 50, "A", "T,G", ["0/1:8,3,4", "0/0:9,1,0"]),
        ])
        csv_path, _ = run_convert(tmp_path, vcf, logger)
        # ALT 含逗号 → 最小引号规则加引号(与原 pandas to_csv 行为一致)
        # |ALT contains comma → quoted per minimal rule (matches pandas to_csv)
        assert open(csv_path).read() == 'chr1,50,A,"T,G",8,3,9,1\n'

    def test_output_name_variants(self, tmp_path, make_vcf, logger):
        for name, expect in [("sample1.vcf", "sample1_deepbsa.csv"),
                             ("sample1.vcf.gz", "sample1_deepbsa.csv"),
                             ("weird.txt", "weird_deepbsa.csv")]:
            csv_path, _ = run_convert(
                tmp_path, make_vcf(name=name), logger,
                out_name="o_" + expect)
            assert csv_path.name == expect


class TestSkipReasons:
    """跳过原因分类|Skip reason categories"""

    def test_format_without_ad_skipped(self, tmp_path, make_vcf, logger):
        vcf = make_vcf(records=[
            ("chr1", 10, "A", "T", ["0/0:11", "0/1:22"]),
        ], fmt="GT:DP")
        csv_path, conv = run_convert(tmp_path, vcf, logger)
        assert open(csv_path).read() == ""
        assert conv.stats["skipped_no_ad"] == 1

    def test_dot_ad_skipped(self, tmp_path, make_vcf, logger):
        vcf = make_vcf(records=[
            ("chr1", 10, "A", "T", ["0/0:.,.", "0/1:20,8"]),
        ])
        csv_path, conv = run_convert(tmp_path, vcf, logger)
        assert open(csv_path).read() == ""
        assert conv.stats["skipped_bad_value"] == 1

    def test_single_value_ad_skipped(self, tmp_path, make_vcf, logger):
        vcf = make_vcf(records=[
            ("chr1", 10, "A", "T", ["0/0:5", "0/1:20,8"]),
        ])
        csv_path, conv = run_convert(tmp_path, vcf, logger)
        assert open(csv_path).read() == ""
        assert conv.stats["skipped_incomplete_ad"] == 1

    def test_one_sample_bad_whole_record_skipped(self, tmp_path, make_vcf,
                                                 logger):
        vcf = make_vcf(records=[
            ("chr1", 10, "A", "T", ["0/0:5,1", "0/1:20,."]),
            ("chr1", 20, "A", "T", ["0/0:6,2", "0/1:7,3"]),
        ])
        csv_path, conv = run_convert(tmp_path, vcf, logger)
        assert open(csv_path).read() == "chr1,20,A,T,6,2,7,3\n"
        assert conv.stats["written"] == 1
        assert conv.stats["skipped_bad_value"] == 1

    def test_column_mismatch_skipped(self, tmp_path, make_vcf, logger):
        # 第2条记录少1个样本列,整行跳过|2nd record lacks one sample column
        lines = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR\tS",
            "chr1\t10\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:5,1\t0/1:6,2",
            "chr1\t20\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:7,3",
        ]
        p = tmp_path / "mis.vcf"
        p.write_text("\n".join(lines) + "\n")
        csv_path, conv = run_convert(tmp_path, str(p), logger)
        assert open(csv_path).read() == "chr1,10,A,T,5,1,6,2\n"
        assert conv.stats["skipped_column_mismatch"] == 1

    def test_header_only_vcf(self, tmp_path, make_vcf, logger):
        csv_path, conv = run_convert(tmp_path, make_vcf(records=[]), logger)
        assert open(csv_path).read() == ""
        assert conv.stats["total"] == 0
        assert conv.stats["sample_count"] == 0

    def test_all_records_skipped_still_exit_ok(self, tmp_path, make_vcf,
                                               logger):
        vcf = make_vcf(records=[
            ("chr1", 10, "A", "T", ["0/0:.", "0/1:."]),
            ("chr1", 20, "A", "T", ["0/0:1", "0/1:2"]),
        ])
        csv_path, conv = run_convert(tmp_path, vcf, logger)
        assert open(csv_path).read() == ""
        assert conv.stats["total"] == 2
        assert conv.stats["written"] == 0
