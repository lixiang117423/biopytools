"""输入文件多分隔符解析测试|Multi-delimiter input file parsing tests

ID文件分隔符不可控:可能是Tab/分号/逗号/空格甚至混合(Excel导出常见)。
解析器应逐行自动嗅探,而不是只认Tab两列。
|ID files come with arbitrary delimiters: tab/semicolon/comma/space or even
mixed (common with Excel exports). The parser should sniff per line, not
require tab-only two columns.
"""

import logging

from biopytools.get_link_from_CNCB.utils import InputFileParser


def _write_ids(tmp_path, content: str) -> str:
    """写入临时ID文件并返回路径|Write an ID file to tmp and return its path"""
    id_file = tmp_path / "ids.txt"
    id_file.write_text(content, encoding="utf-8")
    return str(id_file)


class TestReadAndGroupByProjectDelimiters:
    """read_and_group_by_project 的多分隔符嗅探|Delimiter sniffing in read_and_group_by_project"""

    def test_tab_separated_regression(self, tmp_path):
        """Tab两列是既有行为,必须保持|Tab-separated two columns is existing behavior, must keep"""
        f = _write_ids(tmp_path, "PRJNA123\tSRR111\nPRJNA123\tSRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111", "SRR222"]}

    def test_semicolon_separated(self, tmp_path):
        f = _write_ids(tmp_path, "PRJNA123;SRR111\nPRJNA456;SRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"], "PRJNA456": ["SRR222"]}

    def test_comma_separated(self, tmp_path):
        f = _write_ids(tmp_path, "PRJNA123,SRR111\nPRJNA456,SRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"], "PRJNA456": ["SRR222"]}

    def test_space_separated(self, tmp_path):
        f = _write_ids(tmp_path, "PRJNA123 SRR111\nPRJNA456 SRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"], "PRJNA456": ["SRR222"]}

    def test_comma_space_mix(self, tmp_path):
        """Excel风格"逗号+空格"混排|Excel-style comma+space mix"""
        f = _write_ids(tmp_path, "PRJNA123, SRR111\nPRJNA456 , SRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"], "PRJNA456": ["SRR222"]}

    def test_mixed_delimiters_across_lines(self, tmp_path):
        """同一文件内各行分隔符不同|Different delimiters on different lines in one file"""
        f = _write_ids(tmp_path, "PRJNA123\tSRR111\nPRJNA456,SRR222\nPRJNA789;SRR333\nPRJNA023 SRR444\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {
            "PRJNA123": ["SRR111"],
            "PRJNA456": ["SRR222"],
            "PRJNA789": ["SRR333"],
            "PRJNA023": ["SRR444"],
        }

    def test_header_row_skipped(self, tmp_path):
        """纯字母两列视为表头跳过(真实accession必含数字)|All-alphabetic columns mean header, skipped"""
        f = _write_ids(tmp_path, "ProjectID,RunID\nPRJNA123,SRR111\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"]}

    def test_extra_columns_take_first_two(self, tmp_path):
        """超过两列取前两列,其余忽略|More than two columns: keep the first two, ignore the rest"""
        f = _write_ids(tmp_path, "PRJNA123,SRR111,sample_note\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"]}

    def test_single_field_skipped(self, tmp_path):
        """只有一个字段的行跳过,不进结果|Single-field lines are skipped"""
        f = _write_ids(tmp_path, "SRR111\nPRJNA123,SRR222\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR222"]}

    def test_blank_and_comment_lines_ignored(self, tmp_path):
        f = _write_ids(tmp_path, "\n# comment line\nPRJNA123\tSRR111\n   \n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111"]}

    def test_dedup_and_sort_kept(self, tmp_path):
        """同项目Run ID去重后排序的既有行为保持|Dedup-then-sort per project is kept"""
        f = _write_ids(tmp_path, "PRJNA123\tSRR222\nPRJNA123\tSRR111\nPRJNA123\tSRR111\n")
        result = InputFileParser.read_and_group_by_project(f)
        assert result == {"PRJNA123": ["SRR111", "SRR222"]}

    def test_missing_file_returns_none(self, tmp_path):
        result = InputFileParser.read_and_group_by_project(str(tmp_path / "no_such.txt"))
        assert result is None


class TestValidateInputFileDelimiters:
    """validate_input_file 与解析口径同步|validate_input_file stays consistent with parsing"""

    def test_validate_tab_passes(self, tmp_path):
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "PRJNA123\tSRR111\n"))
        assert is_valid

    def test_validate_semicolon_passes(self, tmp_path):
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "PRJNA123;SRR111\n"))
        assert is_valid

    def test_validate_comma_space_passes(self, tmp_path):
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "PRJNA123, SRR111\n"))
        assert is_valid

    def test_validate_extra_columns_passes(self, tmp_path):
        """解析时>2列取前两列,验证同样放行|Parsing keeps the first two of >2 columns; validation also passes"""
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "PRJNA123,SRR111,note\n"))
        assert is_valid

    def test_validate_header_row_not_counted(self, tmp_path):
        """表头行不算有效数据行|Header rows do not count as valid data lines"""
        is_valid, message = InputFileParser.validate_input_file(_write_ids(tmp_path, "ProjectID,RunID\nPRJNA123,SRR111\n"))
        assert is_valid
        assert "1" in message  # 只有1个有效数据行|Only one valid data line

    def test_validate_header_only_fails(self, tmp_path):
        """仅表头无数据行应失败|Header only, no data lines: should fail"""
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "ProjectID,RunID\n"))
        assert not is_valid

    def test_validate_single_field_fails(self, tmp_path):
        is_valid, _ = InputFileParser.validate_input_file(_write_ids(tmp_path, "SRR111\n"))
        assert not is_valid

    def test_validate_missing_file_fails(self, tmp_path):
        is_valid, _ = InputFileParser.validate_input_file(str(tmp_path / "no_such.txt"))
        assert not is_valid
