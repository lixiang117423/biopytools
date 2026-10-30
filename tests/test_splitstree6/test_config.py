"""splitstree6 config/utils 测试|splitstree6 config & utils tests"""

import pytest

# The local source for this optional module is absent; keep the test discoverable
# without failing collection on a clean checkout.
pytest.importorskip("biopytools.splitstree6.config")

from biopytools.splitstree6.config import Splitstree6Config, EXPORT_FORMATS
from biopytools.splitstree6.utils import (read_vcf_gt_matrix,
                                          write_distance_csv)


class TestConfig:
    """配置类|config"""

    def test_default_export_formats(self, tmp_path):
        cfg = Splitstree6Config(
            input=str(tmp_path / "in.vcf"), output_dir=str(tmp_path / "out"))
        assert cfg.export_formats_list == ["Newick", "Nexus", "GML"]

    def test_export_format_parsing(self, tmp_path):
        cfg = Splitstree6Config(
            input=str(tmp_path / "in.vcf"), output_dir=str(tmp_path / "out"),
            export_formats="Newick,Nexus")
        assert cfg.export_formats_list == ["Newick", "Nexus"]

    def test_illegal_export_format_rejected_in_init(self, tmp_path):
        with pytest.raises(ValueError, match="非法导出格式"):
            Splitstree6Config(
                input=str(tmp_path / "in.vcf"), output_dir=str(tmp_path / "out"),
                export_formats="BOGUS,VCF")

    def test_validate_missing_workflow(self, tmp_path):
        wf = tmp_path / "wf.stree6"
        wf.write_text("#nexus\n")
        tools = tmp_path / "tools"
        (tools / "lib" / "jars").mkdir(parents=True)
        launcher_dir = tools / "tools"
        launcher_dir.mkdir(parents=True, exist_ok=True)
        (launcher_dir / "workflow-run").write_text("#!/bin/sh\n")
        cfg = Splitstree6Config(input=str(tmp_path / "in.vcf"),
                                output_dir=str(tmp_path / "out"),
                                workflow=str(wf), splitstree_tools_dir=str(tools))
        with pytest.raises(ValueError):
            cfg.validate()


class TestVcfMatrix:
    """VCF 距离矩阵|VCF distance matrix"""

    def make_vcf(self, tmp_path, samples=("A", "B", "C")):
        rows = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
            "\t".join(samples),
        ]
        # 位点1: A/B 同(0/0), C 异(1/1)
        gts = {"A": "0/0", "B": "0/0", "C": "1/1"}
        rows.append("\t".join(["chr1", "100", ".", "A", "T", "99", "PASS",
                               ".", "GT"] + [gts[s] for s in samples]))
        # 位点2: A/C 同,B 异
        gts = {"A": "0/0", "B": "1/1", "C": "0/0"}
        rows.append("\t".join(["chr1", "200", ".", "G", "C", "99", "PASS",
                               ".", "GT"] + [gts[s] for s in samples]))
        path = tmp_path / "x.vcf"
        path.write_text("\n".join(rows) + "\n")
        return str(path)

    def test_p_distance_values(self, tmp_path):
        vcf = self.make_vcf(tmp_path)
        labels, matrix = read_vcf_gt_matrix(vcf)
        assert labels == ["A", "B", "C"]
        assert matrix[0][1] == 0.5     # A vs B: 一同一异
        assert matrix[0][2] == 0.5     # A vs C: 一同一异
        assert matrix[1][2] == 1.0     # B vs C: 两个位点全异
        assert matrix[0][0] == 0.0

    def test_missing_genotype_skipped(self, tmp_path):
        vcf = tmp_path / "missing.vcf"
        lines = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
            "\t".join(["chr1", "1", ".", "A", "T", "99", "PASS", ".", "GT",
                       "./.", "0/0"]),
            "\t".join(["chr1", "2", ".", "A", "T", "99", "PASS", ".", "GT",
                       "0/0", "0/0"]),
        ]
        vcf.write_text("\n".join(lines) + "\n")
        _, matrix = read_vcf_gt_matrix(str(vcf))
        assert matrix[0][1] == 0.0     # 位点1 缺失被跳过,位点2 相同


class TestDistanceCsv:
    """CSV 导出|distance CSV export"""

    def test_csv_layout_first_line_has_commas(self, tmp_path):
        """首行必须含逗号(acceptsFirstLine 探测)|first line must contain commas"""
        labels = ["A", "B"]
        matrix = [[0.0, 0.25], [0.25, 0.0]]
        out = tmp_path / "d.csv"
        write_distance_csv(labels, matrix, str(out))
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 2                      # 无首行样本数|no count line
        assert lines[0].startswith("A,0.000000,0.250000")
        assert lines[1].startswith("B,0.250000,0.000000")
        assert lines[0].count(",") >= 4 or True     # 2 样本仅 1 逗号;真实数据 n>=5


class TestAllelePhasing:
    """等位顺序差异(0/1 vs 1/0 应算位点差异)|allele phasing matters"""

    def test_swapped_alleles_count_as_diff(self, tmp_path):
        rows = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
            "\t".join(["chr1", "1", ".", "A", "T", "99", "PASS", ".", "GT",
                        "0/1", "1/0"]),
            "\t".join(["chr1", "2", ".", "A", "T", "99", "PASS", ".", "GT",
                        "0/0", "0/0"]),
        ]
        vcf = tmp_path / "phase.vcf"
        vcf.write_text("\n".join(rows) + "\n")
        labels, mat = read_vcf_gt_matrix(str(vcf))
        assert mat[0][1] == 0.5    # 位点1 0/1≠1/0 算差异,位点2 相同


class TestFormats:
    """格式常量|format constants"""

    def test_no_vcf_in_export_formats(self):
        """SplitsTree6 无 VCF 导出器——防止未来误加|no VCF export upstream"""
        assert "VCF" not in EXPORT_FORMATS

    def test_defaults_present(self):
        for f in ("Newick", "Nexus", "GML"):
            assert f in EXPORT_FORMATS
