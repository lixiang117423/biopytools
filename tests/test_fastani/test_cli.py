"""fastANI Click CLI测试|fastANI Click CLI tests"""

from unittest.mock import patch

from click.testing import CliRunner


def load_cli():
    from biopytools.cli.commands.fastani import fastani
    return fastani


class TestCliHelp:
    def test_help_lists_options(self):
        runner = CliRunner()
        result = runner.invoke(load_cli(), ['-h'])
        assert result.exit_code == 0
        for opt in ('--input', '--query', '--reference', '--output-dir',
                    '--threads', '--kmer', '--frag-len', '--min-fraction'):
            assert opt in result.output

    def test_registered_in_top_level(self):
        from biopytools.cli.main import COMMAND_REGISTRY
        entry = [e for e in COMMAND_REGISTRY if e[1] == 'fastani']
        assert entry and entry[0][0] == 'fastani'


class TestCliPassthrough:
    def test_args_forwarded_to_main(self, tmp_path):
        gdir = tmp_path / 'genomes'; gdir.mkdir()
        (gdir / 'a.fna').write_text('>c\nACGT\n')
        (gdir / 'b.fna').write_text('>c\nACGT\n')
        runner = CliRunner()
        captured = {}

        def fake_module_main():
            import sys
            captured['argv'] = list(sys.argv)

        with patch('biopytools.fastani.main.main', side_effect=fake_module_main):
            result = runner.invoke(
                load_cli(),
                ['-i', str(gdir), '-o', str(tmp_path / 'out'),
                 '-t', '6', '-k', '15', '--frag-len', '2500',
                 '--min-fraction', '0.15'])
        assert result.exit_code == 0
        argv = captured['argv']
        assert argv[argv.index('-i') + 1] == str(gdir)
        assert argv[argv.index('-t') + 1] == '6'
        assert argv[argv.index('-k') + 1] == '15'
        assert argv[argv.index('--frag-len') + 1] == '2500'
        assert argv[argv.index('--min-fraction') + 1] == '0.15'
        assert argv[argv.index('-o') + 1] == str(tmp_path / 'out')

    def test_defaults_forwarded(self, tmp_path):
        gdir = tmp_path / 'genomes'; gdir.mkdir()
        (gdir / 'a.fna').write_text('>c\nACGT\n')
        runner = CliRunner()
        captured = {}

        def fake_module_main():
            import sys
            captured['argv'] = list(sys.argv)

        with patch('biopytools.fastani.main.main', side_effect=fake_module_main):
            result = runner.invoke(load_cli(), ['-i', str(gdir)])
        assert result.exit_code == 0
        argv = captured['argv']
        assert argv[argv.index('-t') + 1] == '12'
        assert argv[argv.index('-o') + 1] == './fastani_output'
        assert argv[argv.index('-k') + 1] == '16'
        assert argv[argv.index('--frag-len') + 1] == '3000'
        assert argv[argv.index('--min-fraction') + 1] == '0.2'
        assert argv[argv.index('--log-level') + 1] == 'INFO'

    def test_missing_input_path_rejected(self, tmp_path):
        runner = CliRunner()
        result = runner.invoke(load_cli(), ['-i', str(tmp_path / 'nope')])
        assert result.exit_code != 0
        assert '不存在' in result.output or 'does not exist' in result.output
