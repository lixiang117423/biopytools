"""vcf2splitstree config/utils 测试|vcf2splitstree config & utils tests"""

import pytest

from biopytools.vcf2splitstree.config import Vcf2SplitstreeConfig
from biopytools.vcf2splitstree.utils import (read_vcf_gt_matrix,
                                             write_distance_csv)


def make_vcf(tmp_path, samples=("A", "B", "C"), name="x.vcf"):
    """构造小 VCF:3 样本 × 5 位点(含缺失/相位符/杂合)|Build a small test VCF"""
    rows = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
        "\t".join(samples),
    ]
    data = [
        (1, ["0/0", "0/0", "1/1"]),      # A=B, C 异
        (2, ["0/0", "1/1", "0/0"]),      # A=C 同
        (3, ["./.", "0/0", "0/0"]),      # A 缺失
        (4, ["1|1", "0/0", "0|1"]),      # 相位符
        (5, ["0/1", "1/1", "0/0"]),      # 杂合
    ]
    for pos, gts in data:
        rows.append("\t".join(["chr1", str(pos), ".", "A", "T", "99",
                               "PASS", ".", "GT"] + gts))
    path = tmp_path / name
    path.write_text("\n".join(rows) + "\n")
    return str(path)


class TestConfig:
    """配置类|config"""

    def test_output_csv_naming(self, tmp_path):
        cfg = Vcf2SplitstreeConfig(input=make_vcf(tmp_path),
                                   output_dir=str(tmp_path / "out"))
        assert cfg.output_csv.name == "x_distances.csv"

    def test_gz_naming_strips_gz(self, tmp_path):
        cfg = Vcf2SplitstreeConfig(input=make_vcf(tmp_path, name="y.vcf.gz"),
                                   output_dir=str(tmp_path / "out"))
        assert cfg.output_csv.name == "y_distances.csv"

    def test_non_vcf_rejected(self, tmp_path):
        p = tmp_path / "x.fasta"
        p.write_text(">a\nACGT\n")
        cfg = Vcf2SplitstreeConfig(input=str(p),
                                   output_dir=str(tmp_path / "out"))
        with pytest.raises(ValueError, match="VCF"):
            cfg.validate()

    def test_missing_input_rejected(self, tmp_path):
        cfg = Vcf2SplitstreeConfig(input=str(tmp_path / "nope.vcf"),
                                   output_dir=str(tmp_path / "out"))
        with pytest.raises(ValueError, match="输入不存在"):
            cfg.validate()


class TestVcfMatrix:
    """VCF 距离矩阵|VCF distance matrix"""

    def test_p_distance_values(self, tmp_path):
        labels, matrix = read_vcf_gt_matrix(make_vcf(tmp_path))
        assert labels == ["A", "B", "C"]
        # 位点1: A=B同 C异;位点2: A=C同 B异;位点3: A缺失跳过
        # 位点4: A(1|1) B(0/0) C(0|1);位点5: A(0/1) B(1/1) C(0/0)
        # A-B: 1同 2异 3跳 4异 5异 → 3/4=0.75
        # A-C: 1异 2同 3跳 4异 5异 → 3/4=0.75
        # B-C: 1异 2异 3同 4异 5异 → 4/5=0.80
        assert abs(matrix[0][1] - 0.75) < 1e-9
        assert abs(matrix[0][2] - 0.75) < 1e-9
        assert abs(matrix[1][2] - 0.80) < 1e-9
        assert matrix[0][0] == 0.0

    def test_swapped_alleles_count_as_diff(self, tmp_path):
        """0/1 vs 1/0 算位点差异|allele phasing matters"""
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
        assert mat[0][1] == 0.5

    def test_missing_genotype_skipped(self, tmp_path):
        rows = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
            "\t".join(["chr1", "1", ".", "A", "T", "99", "PASS", ".", "GT",
                       "./.", "0/0"]),
            "\t".join(["chr1", "2", ".", "A", "T", "99", "PASS", ".", "GT",
                       "0/0", "0/0"]),
        ]
        vcf = tmp_path / "miss.vcf"
        vcf.write_text("\n".join(rows) + "\n")
        _, matrix = read_vcf_gt_matrix(str(vcf))
        assert matrix[0][1] == 0.0


class TestDistanceCsv:
    """CSV 写出|distance CSV output"""

    def test_csv_layout_first_line_has_commas(self, tmp_path):
        """首行必须含逗号(SplitsTree6 自动探测)|first line must contain commas"""
        labels = ["A", "B", "C"]
        matrix = [[0.0, 0.25, 0.5], [0.25, 0.0, 0.75], [0.5, 0.75, 0.0]]
        out = tmp_path / "d.csv"
        write_distance_csv(labels, matrix, str(out))
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 3                    # 无首行样本数|no count line
        assert lines[0].startswith("A,0.000000,0.250000,0.500000")
        assert lines[1].startswith("B,0.250000,0.000000,0.750000")
        assert lines[0].count(",") == 3   # label后1 + 3值间2;真实 n≥5 时 ≥4
