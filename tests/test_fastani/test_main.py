"""fastANI main入口测试|fastANI main entry tests"""

from unittest.mock import MagicMock, patch

import pytest

from biopytools.fastani.main import parse_arguments


class TestParseArguments:
    def test_all_vs_all_args(self):
        args = parse_arguments(['-i', 'genomes/', '-o', 'outdir/'])
        assert args.input == 'genomes/' and args.output_dir == 'outdir/'
        assert args.threads == 12 and args.kmer == 16
        assert args.frag_len == 3000 and args.min_fraction == 0.2

    def test_directional_args(self):
        args = parse_arguments(['-q', 'q.fa', '-r', 'refs/'])
        assert args.query == 'q.fa' and args.reference == 'refs/'
        assert args.input is None

    def test_custom_params(self):
        args = parse_arguments(['-i', 'g/', '-k', '15', '--frag-len', '2000',
                                '--min-fraction', '0.1', '-t', '8'])
        assert (args.kmer, args.frag_len, args.min_fraction, args.threads) == \
            (15, 2000, 0.1, 8)


class TestMain:
    def test_main_runs_calculator_and_exits_zero(self, tmp_path):
        gdir = tmp_path / 'genomes'; gdir.mkdir()
        (gdir / 'a.fna').write_text('>c\nACGT\n')
        (gdir / 'b.fna').write_text('>c\nACGT\n')
        out = tmp_path / 'out'
        fake_bin = tmp_path / 'fastANI'
        fake_bin.write_text('')
        with patch('biopytools.fastani.main.FastaniCalculator') as mock_calc, \
             patch('biopytools.fastani.main.FastaniLogger') as mock_log, \
             patch('biopytools.fastani.config.get_domain_tool_path', return_value=str(fake_bin)):
            mock_calc.return_value.run.return_value = True
            mock_log.return_value.get_logger.return_value = MagicMock()
            with patch('sys.argv', ['fastani', '-i', str(gdir), '-o', str(out)]):
                from biopytools.fastani.main import main
                with pytest.raises(SystemExit) as ei:
                    main()
        assert ei.value.code == 0
        mock_calc.assert_called_once()

    def test_main_validation_error_exits_one(self, tmp_path, capsys):
        with patch('sys.argv', ['fastani', '-o', str(tmp_path / 'out')]):
            from biopytools.fastani.main import main
            with pytest.raises(SystemExit) as ei:
                main()
        assert ei.value.code == 1
        assert '错误|Error' in capsys.readouterr().err
