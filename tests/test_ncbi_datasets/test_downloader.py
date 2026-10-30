"""downloader 核心模块测试|downloader core tests"""

import json
import subprocess
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from biopytools.ncbi_datasets.config import NCBIDatasetsConfig
from biopytools.ncbi_datasets.downloader import NCBIDatasetsDownloader
from biopytools.ncbi_datasets.utils import ModuleLogger


def make_downloader(tmp_path, **kwargs):
    """构造下载器(带临时目录与日志)|Build a downloader with tmp dirs and logger"""
    kwargs.setdefault('output_dir', str(tmp_path / 'out'))
    config = NCBIDatasetsConfig(taxon=67593, **kwargs)
    config.validate()
    logger = ModuleLogger(str(config.log_file)).get_logger()
    return NCBIDatasetsDownloader(config, logger)


def fake_run(returncode=0, stdout='', stderr=''):
    """构造假 CompletedProcess|Build a fake CompletedProcess"""
    return subprocess.CompletedProcess([], returncode, stdout=stdout, stderr=stderr)


SAMPLE_JSON_LINES = '\n'.join([
    json.dumps({'accession': 'GCF_000149865.2',
                'assembly_info': {'assembly_level': 'Complete Genome',
                                  'assembly_status': 'current'},
                'organism': {'organism_name': 'Phytophthora sojae'}}),
    json.dumps({'accession': 'GCA_000149865.3',
                'assembly_info': {'assembly_level': 'Scaffold',
                                  'assembly_status': 'current'},
                'organism': {'organism_name': 'Phytophthora sojae'}}),
    json.dumps({'total_count': 2}),  # 非 assembly 记录|non-assembly record
    'not-json-line',
]) + '\n'


class TestArgs:
    """参数构建|argument building"""

    def test_filter_args_empty(self, tmp_path):
        dl = make_downloader(tmp_path)
        assert dl.build_filter_args() == []

    def test_filter_args_full(self, tmp_path):
        dl = make_downloader(tmp_path, assembly_source='refseq',
                             assembly_level='complete,chromosome',
                             reference=True, annotated=True)
        assert dl.build_filter_args() == [
            '--assembly-source', 'RefSeq',
            '--assembly-level', 'complete,chromosome',
            '--reference', '--annotated',
        ]

    def test_include_default_genome_only(self, tmp_path):
        dl = make_downloader(tmp_path)
        assert dl.build_include_list() == 'genome'

    def test_include_all_flags(self, tmp_path):
        dl = make_downloader(tmp_path, include_gff3=True, include_protein=True,
                             include_cds=True, include_seq_report=True)
        assert dl.build_include_list() == 'genome,gff3,protein,cds,seq-report'


class TestSummary:
    """清单查询|summary query"""

    def test_extract_row_nested_fields(self, tmp_path):
        record = json.loads(SAMPLE_JSON_LINES.splitlines()[0])
        row = NCBIDatasetsDownloader._extract_row(record)
        assert row == {'accession': 'GCF_000149865.2', 'organism': 'Phytophthora sojae',
                       'assembly_level': 'Complete Genome', 'assembly_status': 'current'}

    def test_parse_skips_non_assembly_records(self, tmp_path):
        dl = make_downloader(tmp_path)
        rows = dl._parse_summary_output(SAMPLE_JSON_LINES)
        assert len(rows) == 2

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_run_summary_writes_manifest(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        mock_run.return_value = fake_run(stdout=SAMPLE_JSON_LINES)
        assert dl.run_summary() == 2
        lines = dl.config.manifest_file.read_text().strip().splitlines()
        assert lines[0] == 'accession\torganism\tassembly_level\tassembly_status'
        assert len(lines) == 3
        assert 'GCF_000149865.2' in lines[1]

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_run_summary_query_failure(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        mock_run.return_value = fake_run(returncode=1, stderr='api error')
        assert dl.run_summary() == -1

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_run_summary_no_assemblies(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        mock_run.return_value = fake_run(stdout='')
        assert dl.run_summary() == 0


class TestDownload:
    """下载|download"""

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_skips_when_zip_exists(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        with zipfile.ZipFile(dl.config.zip_file, 'w') as zf:
            zf.writestr('ncbi_dataset/data/x/f.fna', 'ACGT')
        assert dl.download() is True
        mock_run.assert_not_called()

    def test_is_valid_zip(self, tmp_path):
        dl = make_downloader(tmp_path)
        assert not dl._is_valid_zip(dl.config.zip_file)  # 不存在|missing
        dl.config.zip_file.write_bytes(b'')  # 空文件|empty
        assert not dl._is_valid_zip(dl.config.zip_file)
        with zipfile.ZipFile(dl.config.zip_file, 'w') as zf:
            zf.writestr('a.txt', 'x')
        assert dl._is_valid_zip(dl.config.zip_file)
        # 截断:去掉结尾(丢掉 EOCD)|truncated: drop the tail (EOCD gone)
        data = dl.config.zip_file.read_bytes()
        dl.config.zip_file.write_bytes(data[:-22])
        assert not dl._is_valid_zip(dl.config.zip_file)

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_corrupt_zip_redownloaded(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        dl.config.zip_file.write_bytes(b'PK\x03\x04corrupted-truncated')

        def fake_download(cmd, logger, capture=False):
            # 模拟下载重写一个完整 zip|simulate download writing a valid zip
            with zipfile.ZipFile(dl.config.zip_file, 'w') as zf:
                zf.writestr('ncbi_dataset/data/x/f.fna', 'ACGT')
            return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

        mock_run.side_effect = fake_download
        assert dl.download() is True
        mock_run.assert_called_once()  # 损坏 zip 触发重下|corrupt zip triggers re-download
        assert dl._is_valid_zip(dl.config.zip_file)

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_download_checks_resulting_zip(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        mock_run.return_value = fake_run()  # 成功但没写文件|success but no file written
        assert dl.download() is False

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_download_success(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)

        def fake_download(cmd, logger, capture=False):
            with zipfile.ZipFile(dl.config.zip_file, 'w') as zf:
                zf.writestr('ncbi_dataset/data/x/f.fna', 'ACGT')
            return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

        mock_run.side_effect = fake_download
        assert dl.download() is True
        cmd = mock_run.call_args[0][0]
        assert 'download' in cmd and '--include' in cmd and 'genome' in cmd

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_download_failure(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        mock_run.return_value = fake_run(returncode=1, stderr='boom')
        assert dl.download() is False


class TestExtract:
    """解压|extract"""

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_skips_when_dataset_dir_exists(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        dl.config.dataset_dir.mkdir(parents=True)
        assert dl.extract() is True
        mock_run.assert_not_called()

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_missing_zip_rejected(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        assert dl.extract() is False

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_extract_renames_dataset_dir(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        with zipfile.ZipFile(dl.config.zip_file, 'w') as zf:
            zf.writestr('ncbi_dataset/data/GCF_1/f.fna', 'ACGT')
        (dl.config.download_dir / 'ncbi_dataset' / 'data').mkdir(parents=True)
        mock_run.return_value = fake_run()
        assert dl.extract() is True
        assert dl.config.dataset_dir.is_dir()
        assert not (dl.config.download_dir / 'ncbi_dataset').exists()

    @patch('biopytools.ncbi_datasets.downloader.run_command')
    def test_extract_rejects_corrupt_zip(self, mock_run, tmp_path):
        dl = make_downloader(tmp_path)
        dl.config.zip_file.write_bytes(b'PK\x03\x04truncated')
        assert dl.extract() is False
        mock_run.assert_not_called()


class TestVersions:
    """版本信息|software versions"""

    def test_writes_versions_yml(self, tmp_path):
        dl = make_downloader(tmp_path)
        dl.datasets_version = 'v16.21.0'
        dl.write_software_versions()
        content = dl.config.versions_file.read_text()
        assert 'v16.21.0' in content and '1.0.0' in content


def make_fake_dataset(dl, accessions=('GCF_000864765.1',), with_gz=False):
    """构造假 NCBI 解压目录|Build a fake extracted NCBI dataset dir"""
    data_dir = dl.config.dataset_dir / 'data'
    for acc in accessions:
        d = data_dir / acc
        d.mkdir(parents=True)
        (d / f'{acc}_genomic.fna').write_text('>seq\nACGT\n')
        gff_name = 'genomic.gff.gz' if with_gz else 'genomic.gff'
        (d / gff_name).write_text('gff')
        (d / 'protein.faa').write_text('>p\nMKT\n')
        (d / 'cds_from_genomic.fna').write_text('>c\nACGT\n')


class TestOrganize:
    """整理步骤|organize step"""

    def test_unified_symlink_naming(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl)
        assert dl.organize() is True
        genomes = dl.config.organized_dir / 'genomes'
        gff = dl.config.organized_dir / 'gff'
        protein = dl.config.organized_dir / 'protein'
        cds = dl.config.organized_dir / 'cds'
        # 统一后缀: .fa 而非 .fna|unified suffix: .fa not .fna
        assert (genomes / 'GCF_000864765.1.fa').is_symlink()
        assert (gff / 'GCF_000864765.1.gff').is_symlink()
        assert (protein / 'GCF_000864765.1.faa').is_symlink()
        assert (cds / 'GCF_000864765.1.cds.fa').is_symlink()
        assert not (genomes / 'GCF_000864765.1.fna').exists()
        # 软链指向原始文件|symlinks point to original files
        assert (genomes / 'GCF_000864765.1.fa').resolve().name.endswith('_genomic.fna')

    def test_index_file_content(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl)
        dl.organize()
        lines = dl.config.index_file.read_text().strip().splitlines()
        assert lines[0] == 'accession\ttype\tpath'
        assert len(lines) == 5  # genome+gff+protein+cds
        assert any(l.startswith('GCF_000864765.1\tgenome\t') for l in lines)
        assert Path(lines[1].split('\t')[2]).exists()

    def test_gz_gff_keeps_suffix(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl, with_gz=True)
        dl.organize()
        assert (dl.config.organized_dir / 'gff' / 'GCF_000864765.1.gff.gz').is_symlink()

    def test_missing_type_skipped(self, tmp_path):
        dl = make_downloader(tmp_path, include_protein=True)
        make_fake_dataset(dl)
        (dl.config.dataset_dir / 'data' / 'GCF_000864765.1' / 'protein.faa').unlink()
        assert dl.organize() is True
        assert not (dl.config.organized_dir / 'protein' / 'GCF_000864765.1.faa').exists()
        # 索引只含实际整理出的|index only has organized entries
        lines = dl.config.index_file.read_text().strip().splitlines()
        assert len(lines) == 4

    def test_organize_disabled(self, tmp_path):
        dl = make_downloader(tmp_path, organize=False)
        make_fake_dataset(dl)
        assert dl.organize() is True
        assert not (dl.config.organized_dir / 'genomes').exists()

    def test_checkpoint_skips_second_run(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl)
        assert dl.organize() is True
        genomes = dl.config.organized_dir / 'genomes'
        n_before = len(list(genomes.iterdir()))
        assert dl.organize() is True  # 断点续传:直接跳过|checkpoint: skip
        assert len(list(genomes.iterdir())) == n_before

    def test_missing_target_recreated(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl)
        dl.organize()
        target = dl.config.organized_dir / 'genomes' / 'GCF_000864765.1.fa'
        target.unlink()
        assert dl.organize() is True
        assert target.is_symlink()

    def test_existing_file_not_overwritten(self, tmp_path):
        dl = make_downloader(tmp_path)
        make_fake_dataset(dl)
        genomes = dl.config.organized_dir / 'genomes'
        genomes.mkdir(parents=True)
        clash = genomes / 'GCF_000864765.1.fa'
        clash.write_text('user data')
        assert dl.organize() is True
        assert clash.read_text() == 'user data'  # 不覆盖用户文件|keep user file


class TestRun:
    """主流程编排|pipeline orchestration"""

    @patch('biopytools.ncbi_datasets.downloader.ensure_datasets_tool')
    def test_tool_unavailable_stops(self, mock_tool, tmp_path):
        dl = make_downloader(tmp_path)
        mock_tool.return_value = None
        assert dl.run() is False

    @patch.object(NCBIDatasetsDownloader, 'extract')
    @patch.object(NCBIDatasetsDownloader, 'download')
    @patch.object(NCBIDatasetsDownloader, 'run_summary')
    @patch('biopytools.ncbi_datasets.downloader.ensure_datasets_tool')
    def test_dry_run_stops_before_download(self, mock_tool, mock_summary,
                                           mock_download, mock_extract, tmp_path):
        dl = make_downloader(tmp_path, dry_run=True)
        mock_tool.return_value = '/fake/datasets'
        mock_summary.return_value = 2
        assert dl.run() is True
        mock_download.assert_not_called()
        mock_extract.assert_not_called()

    @patch.object(NCBIDatasetsDownloader, 'extract')
    @patch.object(NCBIDatasetsDownloader, 'download')
    @patch.object(NCBIDatasetsDownloader, 'run_summary')
    @patch('biopytools.ncbi_datasets.downloader.ensure_datasets_tool')
    def test_no_assemblies_finishes_without_download(self, mock_tool, mock_summary,
                                                     mock_download, mock_extract, tmp_path):
        dl = make_downloader(tmp_path)
        mock_tool.return_value = '/fake/datasets'
        mock_summary.return_value = 0
        assert dl.run() is True
        mock_download.assert_not_called()
        mock_extract.assert_not_called()

    @patch.object(NCBIDatasetsDownloader, 'write_software_versions')
    @patch.object(NCBIDatasetsDownloader, 'organize')
    @patch.object(NCBIDatasetsDownloader, 'extract')
    @patch.object(NCBIDatasetsDownloader, 'download')
    @patch.object(NCBIDatasetsDownloader, 'run_summary')
    @patch('biopytools.ncbi_datasets.downloader.ensure_datasets_tool')
    def test_full_pipeline(self, mock_tool, mock_summary, mock_download,
                           mock_extract, mock_organize, mock_versions, tmp_path):
        dl = make_downloader(tmp_path)
        mock_tool.return_value = '/fake/datasets'
        mock_summary.return_value = 2
        mock_download.return_value = True
        mock_extract.return_value = True
        mock_organize.return_value = True
        assert dl.run() is True
        mock_download.assert_called_once()
        mock_extract.assert_called_once()
        mock_organize.assert_called_once()
        mock_versions.assert_called()
