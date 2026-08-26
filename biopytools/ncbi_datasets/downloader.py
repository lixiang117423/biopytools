"""NCBI datasets 下载核心模块|NCBI datasets download core module

流程: summary 查清单 → 下载 zip → 解压并按 taxon 命名 → 版本信息
|Pipeline: summary manifest → download zip → extract & rename → software versions
"""

import json
import os
from pathlib import Path

from .config import ASSEMBLY_SOURCE_MAP, INCLUDE_OPTIONS, SUMMARY_FIELDS
from .utils import ensure_datasets_tool, format_number, run_command

DEFAULT_INCLUDE = 'genome'

# 整理规则: (类型, 子目录, ((源文件glob, 目标后缀), ...))
# |Organize rules: (type, subdir, ((source glob, target suffix), ...))
# 目标统一 {accession}{suffix}: 基因组 .fa / gff .gff / 蛋白 .faa / cds .cds.fa
# |Targets use {accession}{suffix}: genome .fa, gff .gff, protein .faa, cds .cds.fa
ORGANIZE_SPECS = [
    ('genome', 'genomes', (('*_genomic.fna', '.fa'), ('*_genomic.fa', '.fa'))),
    ('gff', 'gff', (('genomic.gff', '.gff'), ('genomic.gff.gz', '.gff.gz'))),
    ('protein', 'protein', (('protein.faa', '.faa'),)),
    ('cds', 'cds', (('cds_from_genomic.fna', '.cds.fa'),)),
]
# 与 include 参数的对应(用于缺失告警)|Mapping to include flags (for missing warnings)
INCLUDE_TYPE_MAP = {'include_gff3': 'gff', 'include_protein': 'protein', 'include_cds': 'cds'}


class NCBIDatasetsDownloader:
    """NCBI taxon 基因组批量下载器|NCBI taxon genome batch downloader"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.datasets_version = 'unknown'

    # ---- 参数构建|Argument building ----

    def build_filter_args(self) -> list:
        """构建筛选参数列表|Build filter argument list"""
        args = []
        if self.config.assembly_source:
            args += ['--assembly-source', ASSEMBLY_SOURCE_MAP[self.config.assembly_source]]
        if self.config.assembly_level:
            args += ['--assembly-level', self.config.assembly_level]
        if self.config.reference:
            args += ['--reference']
        if self.config.annotated:
            args += ['--annotated']
        return args

    def build_include_list(self) -> str:
        """构建下载内容列表(默认只下 genome)|Build include list (genome only by default)"""
        includes = [DEFAULT_INCLUDE]
        for attr, name in INCLUDE_OPTIONS.items():
            if getattr(self.config, attr):
                includes.append(name)
        return ','.join(includes)

    # ---- 工具版本|Tool version ----

    def get_version(self) -> str:
        """获取 datasets 版本|Get datasets CLI version"""
        result = run_command(
            [self.config.resolved_datasets_path, '--version'], self.logger, capture=True)
        if result and result.returncode == 0 and result.stdout.strip():
            # 输出形如|Output looks like: "datasets version: 18.23.0"
            version_line = result.stdout.strip().splitlines()[0]
            return version_line.replace('datasets version:', '').strip()
        return 'unknown'

    # ---- 步骤 1: 查清单|Step 1: summary manifest ----

    def run_summary(self) -> int:
        """查询 taxon 下所有 assembly 并写入清单 TSV|Query all assemblies and write manifest TSV

        Returns:
            int: assembly 数量|assembly count, -1 表示查询失败|query failed
        """
        cmd = [self.config.resolved_datasets_path, 'summary', 'genome', 'taxon',
               str(self.config.taxon)]
        cmd += self.build_filter_args()
        cmd += ['--as-json-lines']
        self.logger.info("执行|Executing: NCBI assembly 清单查询|assembly summary query")
        result = run_command(cmd, self.logger, capture=True)
        if result is None or result.returncode != 0:
            stderr = (result.stderr or '').strip() if result is not None else ''
            self.logger.error(f"summary 查询失败|summary query failed: {stderr or '未知错误|unknown error'}")
            return -1

        rows = self._parse_summary_output(result.stdout)
        if not rows:
            self.logger.warning(
                f"未找到任何 assembly|No assemblies found for taxon {self.config.taxon}"
                f" (可能筛选过严|filters may be too strict)")
            self._write_manifest(rows)
            return 0

        self._write_manifest(rows)
        self.logger.info(f"共找到|Total assemblies: {format_number(len(rows))}")
        self.logger.info(f"清单已写入|Manifest written: {self.config.manifest_file}")
        self._log_manifest_preview(rows)
        return len(rows)

    def _parse_summary_output(self, stdout: str) -> list:
        """解析 --as-json-lines 输出|Parse --as-json-lines output"""
        rows = []
        for line in stdout.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if 'accession' not in record:
                continue  # 跳过非 assembly 记录(如 total_count)|skip non-assembly records
            rows.append(self._extract_row(record))
        return rows

    @staticmethod
    def _extract_row(record: dict) -> dict:
        """从 JSON 记录提取清单字段|Extract manifest fields from a JSON record"""
        assembly_info = record.get('assembly_info', {}) or {}
        organism = record.get('organism', {}) or {}
        return {
            'accession': record.get('accession', ''),
            'organism': organism.get('organism_name', ''),
            'assembly_level': assembly_info.get('assembly_level', ''),
            'assembly_status': assembly_info.get('assembly_status', ''),
        }

    def _write_manifest(self, rows: list) -> None:
        """写入清单 TSV|Write manifest TSV"""
        with open(self.config.manifest_file, 'w', encoding='utf-8') as fh:
            fh.write('\t'.join(SUMMARY_FIELDS) + '\n')
            for row in rows:
                fh.write('\t'.join(str(row[f]) for f in SUMMARY_FIELDS) + '\n')

    def _log_manifest_preview(self, rows: list) -> None:
        """日志输出清单预览(前10行)|Log manifest preview (first 10 rows)"""
        preview = '\t'.join(SUMMARY_FIELDS) + '\n'
        for row in rows[:10]:
            preview += '\t'.join(str(row[f]) for f in SUMMARY_FIELDS) + '\n'
        self.logger.info(f"清单预览|Manifest preview:\n{preview}")

    # ---- 步骤 2: 下载|Step 2: download ----

    @staticmethod
    def _is_valid_zip(path: Path) -> bool:
        """快速校验 zip 完整性(EOCD 签名,读文件尾)|Quick zip integrity check (EOCD signature)

        截断的下载包文件头仍是 PK 魔数,必须查结尾的 End-of-central-directory
        签名才能发现|Truncated zips keep the PK header; only the trailing EOCD
        signature reveals truncation
        """
        try:
            size = path.stat().st_size
            if size < 22:  # EOCD 最小 22 字节|EOCD minimum size
                return False
            with open(path, 'rb') as fh:
                fh.seek(max(0, size - 65536))
                tail = fh.read()
            # End-of-central-directory 签名|EOCD signature
            return b'PK\x05\x06' in tail
        except OSError:
            return False

    def download(self) -> bool:
        """下载基因组 zip(存在且完整则跳过)|Download genome zip (skip if present and valid)

        已存在的 zip 会先校验完整性:损坏(下载被中断)则删除重下
        |Existing zips are integrity-checked: corrupt ones are removed and re-downloaded
        """
        zip_file = self.config.zip_file
        if zip_file.exists():
            if self._is_valid_zip(zip_file):
                self.logger.info(
                    f"跳过已完成步骤|Skipping completed step: 下载已存在|zip already "
                    f"exists: {zip_file}")
                return True
            self.logger.warning(f"检测到损坏的 zip,删除并重新下载|Corrupt zip found, "
                                f"removing and re-downloading: {zip_file}")
            try:
                zip_file.unlink()
            except OSError as e:
                self.logger.error(f"删除损坏 zip 失败|Failed to remove corrupt zip: {e}")
                return False

        cmd = [self.config.resolved_datasets_path, 'download', 'genome', 'taxon',
               str(self.config.taxon)]
        cmd += self.build_filter_args()
        cmd += ['--include', self.build_include_list(), '--filename', str(zip_file)]
        self.logger.info("执行|Executing: NCBI 基因组批量下载|genome batch download")
        result = run_command(cmd, self.logger, capture=False)
        if result is None or result.returncode != 0:
            stderr = (result.stderr or '').strip() if result is not None else ''
            self.logger.error(f"下载失败|download failed: {stderr or '未知错误|unknown error'}")
            return False
        if not self._is_valid_zip(zip_file):
            self.logger.error(f"下载产物 zip 损坏|Downloaded zip is corrupt: {zip_file}")
            return False
        self.logger.info(f"下载完成|Download completed: {zip_file}")
        return True

    # ---- 步骤 3: 解压|Step 3: extract ----

    def extract(self) -> bool:
        """解压 zip 并按 taxon 命名目录(已存在则跳过)|Extract zip and rename (skip if exists)"""
        if self.config.dataset_dir.exists():
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: 解压目录已存在|dataset dir already "
                f"exists: {self.config.dataset_dir}")
            return True
        if not self.config.zip_file.exists():
            self.logger.error(f"zip 不存在,无法解压|zip not found: {self.config.zip_file}")
            return False
        if not self._is_valid_zip(self.config.zip_file):
            self.logger.error(f"zip 文件损坏|Corrupt zip file: {self.config.zip_file}")
            return False

        cmd = ['unzip', '-o', '-q', str(self.config.zip_file), '-d', str(self.config.download_dir)]
        self.logger.info("执行|Executing: 解压基因组 zip|unzip genome archive")
        result = run_command(cmd, self.logger, capture=False)
        if result is None or result.returncode != 0:
            stderr = (result.stderr or '').strip() if result is not None else ''
            self.logger.error(f"解压失败|extraction failed: {stderr or '未知错误|unknown error'}")
            return False

        # zip 内顶层目录 ncbi_dataset → {taxon}.ncbi_dataset(§12 命名规范)
        # |rename zip top-level dir ncbi_dataset -> {taxon}.ncbi_dataset
        extracted = self.config.download_dir / 'ncbi_dataset'
        if extracted.exists() and not self.config.dataset_dir.exists():
            os.rename(str(extracted), str(self.config.dataset_dir))
        if not self.config.dataset_dir.exists():
            self.logger.error(f"解压后未找到数据目录|data dir not found after extraction: "
                              f"{self.config.dataset_dir}")
            return False
        self.logger.info(f"解压完成|Extraction completed: {self.config.dataset_dir}")
        return True

    # ---- 步骤 4: 整理|Step 4: organize ----

    def organize(self) -> bool:
        """整理下载产物到 02_organized(软链 + 统一命名)|Organize downloads into 02_organized

        每个 accession 的序列/注释统一为 {accession}.fa/.gff/.faa/.cds.fa 软链,
        按类型分目录,生成 files.tsv 索引(绝对路径,可直接喂下游 fof)
        |Per-accession symlinks with unified names, grouped by type, plus files.tsv index
        """
        if not self.config.organize:
            self.logger.info("跳过整理|Skipping organize (--no-organize)")
            return True
        if not self.config.dataset_dir.is_dir():
            self.logger.error(f"数据目录不存在,无法整理|dataset dir not found: "
                              f"{self.config.dataset_dir}")
            return False
        # 断点续传:索引完整则跳过|Checkpoint: skip when index is complete
        if self.config.index_file.exists() and self._index_targets_complete():
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 整理已存在|"
                             f"organized: {self.config.organized_dir}")
            return True

        data_dir = self.config.dataset_dir / 'data'
        if not data_dir.is_dir():
            self.logger.error(f"数据目录不存在|data dir not found: {data_dir}")
            return False

        rows = []
        type_counts = {ftype: 0 for ftype, _, _ in ORGANIZE_SPECS}
        for accession_dir in sorted(p for p in data_dir.iterdir() if p.is_dir()):
            for accession, ftype, path in self._organize_one(accession_dir, accession_dir.name):
                rows.append((accession, ftype, path))
                type_counts[ftype] += 1

        self._write_index(rows)
        self._log_organize_summary(type_counts)
        return True

    def _organize_one(self, accession_dir, accession):
        """整理单个 accession 目录,返回 (row, type) 列表|Organize one accession dir"""
        results = []
        for ftype, subdir, patterns in ORGANIZE_SPECS:
            source, suffix = None, None
            for glob_pat, target_suffix in patterns:
                matches = sorted(accession_dir.glob(glob_pat))
                if matches:
                    source, suffix = matches[0], target_suffix
                    break
            if source is None:
                continue  # 该类型未下载或缺失|type absent (not included or missing)
            target_dir = self.config.organized_dir / subdir
            target_dir.mkdir(parents=True, exist_ok=True)
            link = target_dir / f'{accession}{suffix}'
            if self._make_link(source, link):
                results.append((accession, ftype, str(link.resolve())))
        return results

    @staticmethod
    def _make_link(source, link) -> bool:
        """建立软链(已存在且指向正确则跳过,不覆盖用户文件)|Create symlink idempotently"""
        if link.is_symlink() and link.resolve() == source.resolve():
            return True
        if link.exists() or link.is_symlink():
            return False  # 已存在但指向别处,不覆盖|exists but points elsewhere, keep
        try:
            os.symlink(source, link)
            return True
        except OSError:
            return False

    def _index_targets_complete(self) -> bool:
        """检查索引中所有目标文件均存在|Check all indexed targets exist"""
        if not self.config.index_file.exists():
            return False
        try:
            with open(self.config.index_file, encoding='utf-8') as fh:
                lines = [l.strip() for l in fh if l.strip()]
        except OSError:
            return False
        if len(lines) <= 1:  # 只有表头|header only
            return False
        for line in lines[1:]:
            parts = line.split('\t')
            # 软链本身必须存在(不跟随目标)|the symlink itself must exist
            if len(parts) < 3 or not os.path.islink(parts[2]) or not Path(parts[2]).exists():
                return False
        return True

    def _write_index(self, rows: list) -> None:
        """写入整理索引 files.tsv|Write organize index files.tsv"""
        with open(self.config.index_file, 'w', encoding='utf-8') as fh:
            fh.write('accession\ttype\tpath\n')
            for accession, ftype, path in sorted(rows):
                fh.write(f'{accession}\t{ftype}\t{path}\n')
        self.logger.info(f"整理索引已写入|Index written: {self.config.index_file}")

    def _log_organize_summary(self, type_counts: dict) -> None:
        """输出整理摘要与缺失告警|Log organize summary and missing warnings"""
        parts = [f'{name}={format_number(n)}' for name, n in type_counts.items() if n]
        self.logger.info(f"整理完成|Organized ({', '.join(parts)}): {self.config.organized_dir}")
        # 用户显式要求的内容类型缺失时告警|Warn when a requested include type is missing
        for attr, ftype in INCLUDE_TYPE_MAP.items():
            if getattr(self.config, attr) and type_counts[ftype] == 0:
                self.logger.warning(f"未找到 {ftype} 文件|No {ftype} files found "
                                    f"(对应 --include-{attr.replace('_', '-')})")

    # ---- 步骤 5: 版本信息|Step 5: software versions ----

    def write_software_versions(self) -> None:
        """写入软件版本信息(§12.5)|Write software versions YAML"""
        import yaml

        from . import __version__
        versions = {
            'ncbi_datasets': {'version': __version__},
            'datasets': {'version': self.datasets_version},
        }
        try:
            with open(self.config.versions_file, 'w', encoding='utf-8') as fh:
                yaml.safe_dump(versions, fh, allow_unicode=True, sort_keys=False)
            self.logger.info(f"版本信息已写入|Versions written: {self.config.versions_file}")
        except OSError as e:
            self.logger.error(f"版本信息写入失败|Failed to write versions file: {e}")

    # ---- 主流程|Main pipeline ----

    def run(self) -> bool:
        """运行完整下载流程|Run the full download pipeline"""
        self.logger.info("开始 NCBI 基因组批量下载|Starting NCBI genome batch download")
        self.logger.info(f"taxon: {self.config.taxon} | 输出目录|output: "
                         f"{self.config.output_path}")

        # 0. 确保工具可用|Ensure tool is ready
        tool = ensure_datasets_tool(self.config.resolved_datasets_path, self.logger)
        if not tool:
            self.logger.error("datasets 工具不可用,流程终止|datasets tool unavailable, "
                              "pipeline terminated")
            return False
        self.config.resolved_datasets_path = tool
        self.datasets_version = self.get_version()
        self.logger.info(f"datasets 版本|datasets version: {self.datasets_version}")

        # 1. 查清单|Summary
        count = self.run_summary()
        if count < 0:
            return False
        if count == 0:
            self.logger.warning("无 assembly 可下载,流程结束|No assemblies to download, "
                                "pipeline finished")
            self.write_software_versions()
            return True

        # 2. 下载|Download
        if self.config.dry_run:
            self.logger.info("dry-run 模式,跳过下载|dry-run mode, skipping download")
            self.write_software_versions()
            return True
        if not self.download():
            return False

        # 3. 解压|Extract
        if not self.extract():
            return False

        # 4. 整理|Organize
        if not self.organize():
            return False

        # 5. 版本信息|Software versions
        self.write_software_versions()
        self._log_result_summary()
        self.logger.info("NCBI 基因组批量下载完成|NCBI genome batch download completed")
        return True

    def _log_result_summary(self) -> None:
        """输出结果摘要|Log result summary"""
        data_dir = self.config.dataset_dir / 'data'
        accession_count = 0
        if data_dir.is_dir():
            accession_count = sum(1 for p in data_dir.iterdir() if p.is_dir())
        self.logger.info(f"结果目录|Result dir: {self.config.dataset_dir}")
        self.logger.info(f"整理目录|Organized dir: {self.config.organized_dir}")
        self.logger.info(f"下载 assembly 数|Assemblies downloaded: {format_number(accession_count)}")
