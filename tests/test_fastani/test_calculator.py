"""fastANI calculator 测试|fastANI calculator tests"""

import os
from unittest.mock import MagicMock, patch

import pytest

from biopytools.fastani.calculator import FastaniCalculator
from biopytools.fastani.config import FastaniConfig


def make_setup(tmp_path, n_genomes=3, mode='all'):
    """构造calculator+假基因组|Build calculator with fake genomes"""
    genomes_dir = tmp_path / 'genomes'
    genomes_dir.mkdir()
    for i in range(n_genomes):
        (genomes_dir / f'g{i}.fna').write_text(f'>c{i}\nACGTACGTAC\n')
    kwargs = {'output_dir': str(tmp_path / 'out')}
    if mode == 'all':
        kwargs['input'] = str(genomes_dir)
    else:
        kwargs['query'] = str(genomes_dir / 'g0.fna')
        kwargs['reference'] = str(genomes_dir)
    cfg = FastaniConfig(**kwargs)
    cfg.fastani_path = str(
        tmp_path / 'miniforge3' / 'envs' / 'pop' / 'bin' / 'fastANI')
    os.makedirs(os.path.dirname(cfg.fastani_path), exist_ok=True)
    with open(cfg.fastani_path, 'w') as fh:
        fh.write('#!/bin/sh\n')
    conda_exe = tmp_path / 'miniforge3' / 'bin' / 'conda'
    conda_exe.parent.mkdir(parents=True, exist_ok=True)
    conda_exe.write_text('#!/bin/sh\nexit 0\n')
    conda_exe.chmod(0o755)
    logger = MagicMock()
    return FastaniCalculator(cfg, logger), cfg, genomes_dir


class TestCollectInputs:
    def test_all_vs_all_same_list(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path, n_genomes=3)
        q, r = calc.collect_inputs()
        assert q == r and len(q) == 3

    def test_directional_mode(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path, mode='dir')
        q, r = calc.collect_inputs()
        assert len(q) == 1 and len(r) == 3

    def test_all_vs_all_needs_two(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path, n_genomes=1)
        with pytest.raises(ValueError, match='2'):
            calc.collect_inputs()

    def test_duplicate_names_rejected(self, tmp_path):
        # 同名不同后缀→同名基因组|same stem different suffix→duplicate name
        gdir = tmp_path / 'genomes'; gdir.mkdir()
        (gdir / 'x.fa').write_text('>c\nACGT\n')
        (gdir / 'x.fna').write_text('>c\nACGT\n')
        cfg = FastaniConfig(input=str(gdir), output_dir=str(tmp_path / 'out'))
        cfg.fastani_path = str(tmp_path / 'bin' / 'fastANI')
        os.makedirs(tmp_path / 'bin', exist_ok=True)
        (tmp_path / 'bin' / 'fastANI').write_text('')
        calc = FastaniCalculator(cfg, MagicMock())
        with pytest.raises(ValueError, match='重复'):
            calc.collect_inputs()


class TestWriteListFiles:
    def test_lists_written_absolute(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        ql, rl = calc.write_list_files(q, r)
        assert os.path.basename(ql) == 'genome_list_ql.txt'
        lines = open(ql).read().strip().split('\n')
        assert len(lines) == 3 and all(os.path.isabs(p) for p in lines)


class TestRunFastani:
    def test_command_construction_and_call(self, tmp_path, monkeypatch):
        calc, cfg, genomes_dir = make_setup(tmp_path)
        monkeypatch.setenv(
            "CONDA_EXE", str(tmp_path / "miniforge3" / "bin" / "conda"))
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)

        def mock_run_with_file(cmd, **kwargs):
            # 写伪输出文件模拟fastANI成功|Write fake output to simulate success
            os.makedirs(calc.fastani_dir, exist_ok=True)
            with open(calc.out_file, 'w') as fh:
                fh.write('dummy\toutput\t99.0\t10\t10\n')
            return MagicMock(returncode=0, stdout='', stderr='')

        with patch('biopytools.fastani.calculator.subprocess.run', side_effect=mock_run_with_file) as mock_run:
            assert calc.run_fastani() is True
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        # conda包装契约(2026-09修复): 绝对conda + 'run -p <环境绝对前缀>',
        # 不再裸调'conda -n'(作业PATH上的系统conda解析不到环境)|
        # Wrapping contract: absolute conda + 'run -p <env prefix>'; bare
        # 'conda -n' is forbidden (a job's PATH conda may miss the env)
        assert os.path.basename(cmd[0]) == 'conda', "必须用conda可执行文件绝对路径"
        assert cmd[1:3] == ['run', '-p']
        assert cmd[3].rstrip('/').endswith('/envs/pop')
        assert '--no-capture-output' in cmd[0:6]
        assert cmd[cmd.index('--ql') + 1].endswith('genome_list_ql.txt')
        assert cmd[cmd.index('--rl') + 1].endswith('genome_list_rl.txt')
        assert '--matrix' in cmd
        assert cmd[cmd.index('-k') + 1] == '16'
        assert cmd[cmd.index('--fragLen') + 1] == '3000'
        assert cmd[cmd.index('--minFraction') + 1] == '0.2'
        assert mock_run.call_args[1]['shell'] is False

    def test_nonzero_exit_returns_false(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)
        with patch('biopytools.fastani.calculator.subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=255, stdout='', stderr='boom')
            assert calc.run_fastani() is False
        assert calc.logger.error.called

    def test_checkpoint_skips_existing_out(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)  # 写真实.checkpoint|writes real .checkpoint
        os.makedirs(calc.fastani_dir, exist_ok=True)
        with open(os.path.join(calc.fastani_dir, 'fastani.out'), 'w') as fh:
            fh.write('a.fna\tb.fna\t99.0\t10\t10\n')
        with patch('biopytools.fastani.calculator.subprocess.run') as mock_run:
            assert calc.run_fastani() is True
        mock_run.assert_not_called()
        assert calc.logger.info.called

    def test_checkpoint_invalidated_on_input_change(self, tmp_path):
        # 换输入集合后旧.out必须重跑|Stale .out must re-run after input change
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)
        with open(calc.out_file, 'w') as fh:
            fh.write('stale\toutput\t99.0\t10\t10\n')

        (genomes_dir / 'g9.fna').write_text('>c9\nACGTACGTAC\n')
        calc.write_list_files(q + [str(genomes_dir / 'g9.fna')], r)

        def fake_run(cmd, **kwargs):
            # 运行前旧产物应已被删除|Stale output must be removed before run
            assert not os.path.exists(calc.out_file)
            with open(calc.out_file, 'w') as fh:
                fh.write('fresh\toutput\t99.0\t10\t10\n')
            return MagicMock(returncode=0, stdout='', stderr='')

        with patch('biopytools.fastani.calculator.subprocess.run',
                   side_effect=fake_run) as mock_run:
            assert calc.run_fastani() is True
        mock_run.assert_called_once()
        assert calc.logger.warning.called
        assert open(calc.out_file).read().startswith('fresh')

    def test_empty_out_with_zero_exit_degrades(self, tmp_path):
        # 全部配对<80%时.out为空但rc=0→WARNING继续|empty .out with rc 0→WARNING, continue
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)
        os.makedirs(calc.fastani_dir, exist_ok=True)
        open(calc.out_file, 'w').close()  # 预置空文件|pre-create empty file
        with patch('biopytools.fastani.calculator.subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout='', stderr='')
            assert calc.run_fastani() is True
        assert calc.logger.warning.called

    def test_missing_out_with_zero_exit_errors(self, tmp_path):
        # rc=0但输出文件缺失→ERROR|exit 0 but output file missing→ERROR
        calc, cfg, genomes_dir = make_setup(tmp_path)
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)
        with patch('biopytools.fastani.calculator.subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout='', stderr='')
            assert calc.run_fastani() is False
        assert calc.logger.error.called


from biopytools.fastani.calculator import (
    build_ani_matrix, build_nearest_table, parse_fastani_out)


def write_out(tmp_path, lines):
    out = tmp_path / 'fastani.out'
    out.write_text('\n'.join(lines) + '\n')
    return str(out)


class TestParseFastaniOut:
    def test_tab_delimited(self, tmp_path):
        f = write_out(tmp_path, [
            '/g/A.fna\t/g/B.fna\t97.7507\t1303\t1608',
            '/g/B.fna\t/g/A.fna\t97.7601\t1310\t1600'])
        records, malformed = parse_fastani_out(f)
        assert malformed == 0
        assert records[0] == {'query': 'A', 'reference': 'B', 'ani': 97.7507,
                              'matched_frags': 1303, 'total_frags': 1608}

    def test_space_fallback(self, tmp_path):
        f = write_out(tmp_path, ['/g/A.fna /g/B.fna 97.75 100 110'])
        records, malformed = parse_fastani_out(f)
        assert malformed == 0 and records[0]['ani'] == 97.75

    def test_malformed_counted_skipped(self, tmp_path):
        f = write_out(tmp_path, ['/g/A.fna\t/g/B.fna\t97.7\t10\t11', 'garbage line'])
        records, malformed = parse_fastani_out(f)
        assert malformed == 1 and len(records) == 1

    def test_empty_file(self, tmp_path):
        f = write_out(tmp_path, [])
        records, malformed = parse_fastani_out(f)
        assert records == [] and malformed == 0


class TestBuildAniMatrix:
    def test_bidirectional_average_and_diagonal(self, tmp_path):
        records = [
            {'query': 'A', 'reference': 'B', 'ani': 97.0, 'matched_frags': 10, 'total_frags': 10},
            {'query': 'B', 'reference': 'A', 'ani': 99.0, 'matched_frags': 10, 'total_frags': 10},
        ]
        out = str(tmp_path / 'ani_matrix.tsv')
        build_ani_matrix(records, ['A', 'B'], ['A', 'B'], True, out)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[0] == ['genome', 'A', 'B']
        assert rows[1] == ['A', '100.0000', '98.0000']
        assert rows[2] == ['B', '98.0000', '100.0000']

    def test_one_direction_and_na(self, tmp_path):
        records = [{'query': 'A', 'reference': 'C', 'ani': 90.5,
                    'matched_frags': 5, 'total_frags': 10}]
        out = str(tmp_path / 'ani_matrix.tsv')
        build_ani_matrix(records, ['A', 'B'], ['A', 'C'], True, out)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[1] == ['A', '100.0000', '90.5000']
        assert rows[2] == ['B', 'NA', 'NA']

    def test_rectangular_no_special_diagonal(self, tmp_path):
        records = []
        out = str(tmp_path / 'ani_matrix.tsv')
        build_ani_matrix(records, ['A'], ['A', 'B'], False, out)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[1] == ['A', 'NA', 'NA']


class TestBuildNearestTable:
    def test_best_hit_with_af_sorted_desc(self, tmp_path):
        records = [
            {'query': 'A', 'reference': 'B', 'ani': 96.0, 'matched_frags': 8, 'total_frags': 10},
            {'query': 'A', 'reference': 'C', 'ani': 98.0, 'matched_frags': 3, 'total_frags': 10},
            {'query': 'B', 'reference': 'A', 'ani': 97.0, 'matched_frags': 10, 'total_frags': 10},
        ]
        out = str(tmp_path / 'nearest_genome.tsv')
        build_nearest_table(records, ['A', 'B'], ['A', 'B', 'C'], True, out)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[0] == ['genome', 'nearest_genome', 'ani', 'aligned_fraction',
                           'matched_frags', 'total_frags']
        # ANI降序:A(98)在前B(97)在后|ANI descending: A (98) before B (97)
        assert rows[1] == ['A', 'C', '98.0000', '0.3000', '3', '10']
        assert rows[2] == ['B', 'A', '97.0000', '1.0000', '10', '10']

    def test_self_pair_never_nearest_in_all_vs_all(self, tmp_path):
        # all-vs-all自身行不得成为最近邻|self rows must not become nearest
        records = [
            {'query': 'A', 'reference': 'A', 'ani': 100.0,
             'matched_frags': 10, 'total_frags': 10},
            {'query': 'A', 'reference': 'B', 'ani': 97.0,
             'matched_frags': 9, 'total_frags': 10},
        ]
        out = str(tmp_path / 'nearest_genome.tsv')
        build_nearest_table(records, ['A'], ['A', 'B'], True, out)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[1][1] == 'B' and rows[1][2] == '97.0000'

    def test_all_na_row_last_with_warning(self, tmp_path):
        out = str(tmp_path / 'nearest_genome.tsv')
        logger = MagicMock()
        build_nearest_table([], ['A'], ['B'], True, out, logger=logger)
        rows = [line.split('\t') for line in open(out).read().strip().split('\n')]
        assert rows[1] == ['A', 'NA', 'NA', 'NA', 'NA', 'NA']
        assert logger.warning.called


class TestIteratedMode:
    """遍历 1-vs-all 模式|iterated 1-vs-all mode"""

    def make_calc(self, tmp_path, n_genomes=3, threshold=2):
        """构造超过阈值的 all-vs-all 计算器|calculator above threshold"""
        from biopytools.fastani.config import FastaniConfig
        genomes_dir = tmp_path / 'genomes'
        genomes_dir.mkdir()
        for i in range(n_genomes):
            (genomes_dir / f'g{i}.fna').write_text(f'>c{i}\nACGTACGTAC\n')
        cfg = FastaniConfig(input=str(genomes_dir),
                            output_dir=str(tmp_path / 'out'),
                            iterated_threshold=threshold)
        cfg.fastani_path = str(tmp_path / 'fastANI')
        calc = FastaniCalculator(cfg, MagicMock())
        calc._prepare_dirs()
        return calc, genomes_dir

    def test_under_threshold_uses_single_run(self, tmp_path):
        """≤阈值不切遍历|below threshold keeps single run"""
        calc, _ = self.make_calc(tmp_path, n_genomes=2, threshold=5)
        def fake_single(ql, rl, out):
            open(out, 'w').write('g0\tg1\t99.0\t9\t10\n')
            return True
        with patch.object(calc, '_run_single', side_effect=fake_single) as single, \
             patch.object(calc, '_run_iterated', return_value=True) as it:
            calc.collect_inputs()
            q, r = calc.collect_inputs()
            calc.write_list_files(q, r)
            calc._n_query = 2
            assert calc.run_fastani(q) is True
            single.assert_called_once()
            it.assert_not_called()

    def test_over_threshold_uses_iterated(self, tmp_path):
        """超过阈值切遍历|over threshold switches to iterated"""
        calc, _ = self.make_calc(tmp_path, n_genomes=3, threshold=2)
        calc.collect_inputs()
        q, r = calc.collect_inputs()
        calc.write_list_files(q, r)
        calc._n_query = 3
        with patch.object(calc, '_run_single', return_value=True) as single, \
             patch.object(calc, '_run_iterated', return_value=True) as it:
            assert calc.run_fastani(q) is True
            it.assert_called_once()
            single.assert_not_called()

    def test_iterated_batches_and_aggregate(self, tmp_path):
        """遍历:每批 1 个 query,汇总到 fastani.out|one query per batch, aggregated"""
        calc, genomes_dir = self.make_calc(tmp_path, n_genomes=3, threshold=2)
        q, r = calc.collect_inputs()
        calc._n_query = 3

        def fake_single(ql, rl, out_batch):
            qname = open(ql).read().strip().split('/')[-1]
            refs = [x.split('/')[-1] for x in open(rl).read().strip().split('\n')]
            with open(out_batch, 'w') as fh:
                for ref in refs:
                    fh.write(f'{qname}\t{ref}\t99.0\t9\t10\n')
            return True

        with patch.object(calc, '_run_single', side_effect=fake_single):
            assert calc._run_iterated(q) is True
        # 汇总文件存在且含全部配对|aggregate has all pairs
        lines = open(calc.out_file).read().strip().split('\n')
        assert len(lines) == 9   # 3 query × 3 ref(含自比对)
        # 批次目录有 3 个 ql + 3 个 rl|3 ql + 3 rl files
        import glob
        assert len(glob.glob(os.path.join(calc.fastani_dir, 'batches', 'ql_*.txt'))) == 3
        assert len(glob.glob(os.path.join(calc.fastani_dir, 'batches', 'rl_*.txt'))) == 3

    def test_iterated_resume_skips_done_batches(self, tmp_path):
        """断点续传:已完成批次跳过|resume skips completed batches"""
        calc, genomes_dir = self.make_calc(tmp_path, n_genomes=3, threshold=2)
        q, r = calc.collect_inputs()
        calc._n_query = 3
        # 预写 batch_0000_0000.out 模拟已完成|pre-write first batch
        batch_dir = os.path.join(calc.fastani_dir, 'batches')
        os.makedirs(batch_dir, exist_ok=True)
        with open(os.path.join(batch_dir, 'batch_0000_0000.out'), 'w') as fh:
            fh.write('g0\tg0\t100.0\t10\t10\n')

        calls = []

        def fake_single(ql, rl, out_batch):
            calls.append(os.path.basename(out_batch))
            qname = open(ql).read().strip().split('/')[-1]
            with open(out_batch, 'w') as fh:
                fh.write(f'{qname}\t{qname}\t100.0\t10\t10\n')
            return True

        with patch.object(calc, '_run_single', side_effect=fake_single):
            assert calc._run_iterated(q) is True
        assert calls == ['batch_0001_0000.out', 'batch_0002_0000.out']  # 0000 跳过


class TestRefBatching:
    """reference 分批|reference batching"""

    def test_small_ref_batch_splits(self, tmp_path):
        """ref_batch_size=2 时 3 个 ref 拆成 2+1 两批|2+1 batches for 3 refs"""
        from biopytools.fastani.config import FastaniConfig
        genomes_dir = tmp_path / 'genomes'
        genomes_dir.mkdir()
        for i in range(3):
            (genomes_dir / f'g{i}.fna').write_text(f'>c{i}\nACGT\n')
        cfg = FastaniConfig(input=str(genomes_dir),
                            output_dir=str(tmp_path / 'out'),
                            iterated_threshold=2, ref_batch_size=2)
        cfg.fastani_path = str(tmp_path / 'fastANI')
        calc = FastaniCalculator(cfg, MagicMock())
        calc._prepare_dirs()
        q, _ = calc.collect_inputs()
        calc._n_query = 3
        batch_dir = os.path.join(calc.fastani_dir, 'batches')
        os.makedirs(batch_dir, exist_ok=True)

        def fake_single(ql, rl, out_batch):
            qname = open(ql).read().strip().split('/')[-1]
            refs = [x.split('/')[-1] for x in open(rl).read().strip().split('\n')]
            with open(out_batch, 'w') as fh:
                for ref in refs:
                    fh.write(f'{qname}\t{ref}\t99.0\t9\t10\n')
            return True

        with patch.object(calc, '_run_single', side_effect=fake_single):
            assert calc._run_iterated(q) is True
        # 3 query × 2 批 = 6 个批次输出|3 queries × 2 batches = 6 outputs
        import glob
        outs = sorted(glob.glob(os.path.join(batch_dir, 'batch_*.out')))
        assert len(outs) == 6
        # 汇总 3×3=9 行|aggregate 9 rows
        assert len(open(calc.out_file).read().strip().split('\n')) == 9


class TestRunOrchestration:
    def test_run_end_to_end_mocked(self, tmp_path):
        calc, cfg, genomes_dir = make_setup(tmp_path, n_genomes=2)
        calc.logger = MagicMock()

        def fake_run(cmd, **kwargs):
            # 写伪fastANI输出|Write fake fastANI output
            os.makedirs(calc.fastani_dir, exist_ok=True)
            g0 = os.path.join(str(genomes_dir), 'g0.fna')
            g1 = os.path.join(str(genomes_dir), 'g1.fna')
            with open(calc.out_file, 'w') as fh:
                fh.write(f'{g0}\t{g1}\t99.0000\t9\t10\n{g1}\t{g0}\t99.3580\t9\t10\n')
            with open(calc.out_file + '.matrix', 'w') as fh:
                fh.write('2\n')
            return MagicMock(returncode=0, stdout='', stderr='')

        with patch('biopytools.fastani.calculator.subprocess.run', side_effect=fake_run):
            assert calc.run() is True

        matrix = open(os.path.join(calc.results_dir, 'ani_matrix.tsv')).read()
        assert 'g0' in matrix and '99.1790' in matrix  # (99.0000+99.3580)/2=99.1790
        nearest = open(os.path.join(calc.results_dir, 'nearest_genome.tsv')).read()
        assert 'nearest_genome' in nearest
        versions = os.path.join(calc.pipeline_info_dir, 'software_versions.yml')
        assert os.path.exists(versions)
