"""vcf2deepbsa main/CLI 测试|vcf2deepbsa main/CLI tests"""

from biopytools.vcf2deepbsa.main import main, parse_arguments


class TestParseArguments:
    """参数解析|Argument parsing"""

    def test_required_and_defaults(self):
        args = parse_arguments(["-i", "x.vcf"])
        assert args.input_vcf == "x.vcf"
        assert args.output_dir == "./vcf2deepbsa_output"
        assert args.log_level == "INFO"

    def test_full_form(self):
        args = parse_arguments(["-i", "x.vcf", "-o", "/tmp/o",
                                "--log-level", "DEBUG"])
        assert args.output_dir == "/tmp/o"
        assert args.log_level == "DEBUG"


class TestMain:
    """端到端入口|End-to-end entry"""

    def test_main_converts_and_exits_zero(self, tmp_path, make_vcf):
        vcf = make_vcf()
        out = str(tmp_path / "out")
        rc = main(["-i", vcf, "-o", out])
        assert rc == 0
        assert (tmp_path / "out" / "x_deepbsa.csv").is_file()
        assert (tmp_path / "out" / "vcf2deepbsa.log").is_file()

    def test_main_missing_input_returns_one(self, tmp_path):
        rc = main(["-i", str(tmp_path / "none.vcf"),
                   "-o", str(tmp_path / "out")])
        assert rc == 1
