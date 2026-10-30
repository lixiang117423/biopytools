"""fastani utils 测试|fastani utils tests"""

import logging
import os

import pytest

from biopytools.fastani.utils import (
    FastaniLogger, collect_genome_files, genome_name, is_fasta_file)


class TestGenomeName:
    def test_plain_fasta(self):
        assert genome_name('/data/R0590-6.fna') == 'R0590-6'

    def test_gz_stripped_layered(self):
        assert genome_name('/data/R0590-6.fa.gz') == 'R0590-6'

    def test_fasta_gz(self):
        assert genome_name('/data/sample.fasta.gz') == 'sample'

    def test_upper_suffix(self):
        assert genome_name('/data/S1.FNA') == 'S1'

    def test_unknown_suffix_kept(self):
        assert genome_name('/data/readme.txt') == 'readme.txt'


class TestIsFastaFile:
    def test_supported(self):
        assert is_fasta_file('/a/x.fa') and is_fasta_file('/a/x.fna.gz')
        assert is_fasta_file('/a/x.FASTA')

    def test_unsupported(self):
        assert not is_fasta_file('/a/x.txt')
        assert not is_fasta_file('/a/ql.txt')


class TestCollectGenomeFiles:
    def test_directory_glob_sorted(self, tmp_path):
        for name in ('b.fa', 'a.fna', 'c.fasta', 'd.fa.gz', 'note.txt', 'log.out'):
            (tmp_path / name).write_text('>c\nACGT\n')
        result = collect_genome_files(str(tmp_path))
        assert [os.path.basename(p) for p in result] == ['a.fna', 'b.fa', 'c.fasta', 'd.fa.gz']
        assert all(os.path.isabs(p) for p in result)

    def test_single_fasta_file(self, tmp_path):
        g = tmp_path / 'x.fa'; g.write_text('>c\nACGT\n')
        assert len(collect_genome_files(str(g))) == 1

    def test_list_file_with_comments_and_blanks(self, tmp_path):
        g1 = tmp_path / 'a.fa'; g1.write_text('>c\nACGT\n')
        g2 = tmp_path / 'b.fa'; g2.write_text('>c\nACGT\n')
        lst = tmp_path / 'list.txt'
        lst.write_text(f'# 注释|comment\n\n{g1}\n  {g2}  \n')
        result = collect_genome_files(str(lst))
        assert len(result) == 2

    def test_list_line_missing_file_raises(self, tmp_path):
        lst = tmp_path / 'list.txt'
        lst.write_text(str(tmp_path / 'ghost.fa') + '\n')
        with pytest.raises(ValueError, match='ghost'):
            collect_genome_files(str(lst))

    def test_list_line_not_fasta_raises(self, tmp_path):
        other = tmp_path / 'x.txt'; other.write_text('x')
        lst = tmp_path / 'list.txt'
        lst.write_text(str(other) + '\n')
        with pytest.raises(ValueError, match='FASTA'):
            collect_genome_files(str(lst))

    def test_empty_directory_raises(self, tmp_path):
        empty = tmp_path / 'empty'; empty.mkdir()
        with pytest.raises(ValueError, match='未找到'):
            collect_genome_files(str(empty))

    def test_empty_fasta_in_directory_raises(self, tmp_path):
        (tmp_path / 'a.fna').write_text('>c\nACGT\n')
        (tmp_path / 'empty.fna').write_text('')
        with pytest.raises(ValueError, match='基因组文件为空'):
            collect_genome_files(str(tmp_path))

    def test_empty_file_in_list_raises(self, tmp_path):
        g = tmp_path / 'a.fa'; g.write_text('>c\nACGT\n')
        e = tmp_path / 'e.fa'; e.write_text('')
        lst = tmp_path / 'list.txt'
        lst.write_text(f'{g}\n{e}\n')
        with pytest.raises(ValueError, match='Empty genome file'):
            collect_genome_files(str(lst))

    def test_empty_single_fasta_raises(self, tmp_path):
        g = tmp_path / 'x.fa'; g.write_text('')
        with pytest.raises(ValueError, match='基因组文件为空'):
            collect_genome_files(str(g))

    def test_missing_path_raises(self, tmp_path):
        with pytest.raises(ValueError, match='不存在'):
            collect_genome_files(str(tmp_path / 'nope'))


class TestFastaniLogger:
    def test_three_log_files_and_streams(self, tmp_path, capsys):
        logs_dir = tmp_path / '99_logs'
        manager = FastaniLogger(str(logs_dir), log_level='INFO')
        logger = manager.get_logger()
        logger.info("开始|Starting")
        logger.warning("警告|Warning")
        for fh in logging.getLogger('biopytools.fastani').handlers:
            fh.flush()
        assert (logs_dir / 'fastani.log').exists()
        assert (logs_dir / 'fastani_out.log').exists()
        assert (logs_dir / 'fastani_err.log').exists()
        all_log = (logs_dir / 'fastani.log').read_text()
        out_log = (logs_dir / 'fastani_out.log').read_text()
        err_log = (logs_dir / 'fastani_err.log').read_text()
        assert '开始|Starting' in all_log and '警告|Warning' in all_log
        assert '开始|Starting' in out_log and '警告|Warning' not in out_log
        assert '警告|Warning' in err_log and '开始|Starting' not in err_log
        captured = capsys.readouterr()
        assert '开始|Starting' in captured.out
        assert '警告|Warning' in captured.err
