"""vcf2deepbsa utils 测试|vcf2deepbsa utils tests"""

from biopytools.vcf2deepbsa.utils import ModuleLogger, VcfRecord, VcfReader


class TestVcfRecord:
    """单行解析|Single line parsing"""

    def test_parse_basic_fields(self):
        line = ("chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:AD\t"
                "0/0:10,5\t0/1:20,8")
        r = VcfRecord(line)
        assert r.CHROM == "chr1"
        assert r.POS == "100"
        assert r.REF == "A"
        assert r.ALT == "T"
        assert r.FORMAT == ["GT", "AD"]
        assert len(r.GT) == 2
        assert r.GT[0]["AD"] == "10,5"
        assert r.GT[1]["GT"] == "0/1"

    def test_sample_count(self):
        line = ("chr1\t1\t.\tA\tT\t.\t.\t.\tGT:AD:DP\t"
                "0/0:5,1:6\t0/1:4,2:6\t1/1:0,9:9")
        r = VcfRecord(line)
        assert len(r.GT) == 3
        assert r.GT[2]["DP"] == "9"


class TestVcfReader:
    """VCF 迭代器|VCF iterator"""

    def test_iterates_records(self, make_vcf):
        with VcfReader(make_vcf()) as reader:
            records = list(reader)
        assert len(records) == 3
        assert records[0].CHROM == "chr1"
        assert records[-1].CHROM == "chr2"

    def test_samples_from_header(self, make_vcf):
        with VcfReader(make_vcf()) as reader:
            assert reader.samples == ["pool_R", "pool_S"]
            list(reader)

    def test_gz_input(self, make_vcf):
        with VcfReader(make_vcf(name="x.vcf.gz")) as reader:
            assert len(list(reader)) == 3

    def test_header_only_yields_nothing(self, make_vcf):
        with VcfReader(make_vcf(records=[])) as reader:
            assert list(reader) == []


class TestModuleLogger:
    """日志管理器|Logger manager"""

    def test_three_handlers_and_format(self, tmp_path):
        mgr = ModuleLogger(log_file=str(tmp_path / "test.log"))
        log = mgr.get_logger()
        kinds = {type(h).__name__ for h in log.handlers}
        assert kinds == {"StreamHandler", "FileHandler"}
        # stdout(INFO+)+stderr(WARNING+)+file = 3 handler
        assert len(log.handlers) == 3
        for h in log.handlers:
            fmt = h.formatter._fmt
            assert "%(msecs)03d" in fmt
            assert fmt.count("[") == 0

    def test_no_log_file_still_works(self):
        log = ModuleLogger().get_logger()
        log.info("无文件日志|No file logging")
