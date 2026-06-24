"""
QIIME2分类器训练与缓存模块|QIIME2 Classifier Training and Caching Module

将原始参考数据库(SILVA/UNITE)导入并训练为sklearn分类器(.qza),训练结果缓存复用
|Import raw reference DBs (SILVA/UNITE) and train sklearn classifier (.qza), with caching
"""

import gzip
import hashlib
import json
import os
import tarfile
from typing import Optional, Tuple

from .config import Qiime2Config
from .utils import CommandRunner, build_conda_command
from .config import (
    SILVA_NR99_FNAME, UNITE_TGZ_FNAME, UNITE_FASTA_FNAME, UNITE_TAXONOMY_FNAME
)


class ClassifierTrainer:
    """分类器训练器|Classifier Trainer"""

    def __init__(self, config: Qiime2Config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd = cmd_runner

    # ------------------------------------------------------------------
    # 原始库解析|Raw database parsing
    # ------------------------------------------------------------------

    def parse_silva(self, gz_path: str, work_dir: str) -> Tuple[str, str]:
        """
        解析SILVA gz,拆分为序列文件和分类文件|Parse SILVA gz into seqs + taxonomy

        SILVA标题行格式: >ID 分号分隔8级分类(最后一级物种名可含空格)
        |SILVA header: >ID semicolon-separated 8-level taxonomy (species may contain spaces)

        Args:
            gz_path: SILVA gz文件路径|SILVA gz file path
            work_dir: 工作目录|working directory

        Returns:
            (seqs.fasta, taxonomy.tsv)路径|(seqs.fasta, taxonomy.tsv) paths
        """
        seqs_path = os.path.join(work_dir, 'silva_seqs.fasta')
        tax_path = os.path.join(work_dir, 'silva_taxonomy.tsv')

        self.logger.info(f"解析SILVA参考库|Parsing SILVA reference: {gz_path}")

        seq_count = 0
        with gzip.open(gz_path, 'rt') as fin, \
             open(seqs_path, 'w') as f_seq, \
             open(tax_path, 'w') as f_tax:
            for line in fin:
                if line.startswith('>'):
                    # 标题行: >ID taxonomy|header line: >ID taxonomy
                    header = line[1:].rstrip('\n')
                    parts = header.split(None, 1)
                    seq_id = parts[0]
                    taxonomy = parts[1] if len(parts) > 1 else 'Unknown'
                    f_seq.write(f">{seq_id}\n")
                    f_tax.write(f"{seq_id}\t{taxonomy}\n")
                    seq_count += 1
                else:
                    # 序列行,原样写入|sequence line, write as-is
                    f_seq.write(line)

        self.logger.info(f"SILVA解析完成|SILVA parsed: {seq_count} 条序列|sequences")
        return seqs_path, tax_path

    def extract_unite(self, tgz_path: str, work_dir: str) -> Tuple[str, str]:
        """
        解压UNITE tarball,提取99%阈值的fasta和taxonomy|Extract UNITE 99% fasta and taxonomy

        Args:
            tgz_path: UNITE tgz文件路径|UNITE tgz file path
            work_dir: 工作目录|working directory

        Returns:
            (fasta, taxonomy.txt)路径
        """
        self.logger.info(f"解压UNITE参考库|Extracting UNITE reference: {tgz_path}")

        fasta_out = os.path.join(work_dir, 'unite_seqs.fasta')
        tax_out = os.path.join(work_dir, 'unite_taxonomy.txt')

        with tarfile.open(tgz_path, 'r:gz') as tar:
            # 找到目标成员|Find target members
            fasta_member = None
            tax_member = None
            for member in tar.getnames():
                basename = os.path.basename(member)
                if basename == UNITE_FASTA_FNAME:
                    fasta_member = member
                elif basename == UNITE_TAXONOMY_FNAME:
                    tax_member = member

            if not fasta_member or not tax_member:
                raise RuntimeError(
                    f"UNITE tarball缺少必需文件|UNITE tarball missing required files: "
                    f"{UNITE_FASTA_FNAME} / {UNITE_TAXONOMY_FNAME}"
                )

            # 提取到工作目录(去掉子目录)|Extract to work dir (strip subdirectory)
            for member, out_path in [(fasta_member, fasta_out), (tax_member, tax_out)]:
                fobj = tar.extractfile(member)
                if fobj is None:
                    raise RuntimeError(f"无法提取UNITE文件|Cannot extract UNITE file: {member}")
                with open(out_path, 'wb') as f:
                    f.write(fobj.read())

        self.logger.info("UNITE提取完成|UNITE extracted")
        return fasta_out, tax_out

    # ------------------------------------------------------------------
    # 导入与训练命令|Import and training commands
    # ------------------------------------------------------------------

    def _import_sequences(self, fasta_path: str, output_qza: str) -> bool:
        """导入参考序列|Import reference sequences"""
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'import',
            '--type', 'FeatureData[Sequence]',
            '--input-path', fasta_path,
            '--output-path', output_qza,
            '--input-format', 'DNAFASTAFormat',
            '--validate-level', self.config.validate_level,
        ])
        return self.cmd.run(cmd, f"导入参考序列|Import reference sequences: {os.path.basename(fasta_path)}")[0]

    def _import_taxonomy(self, tax_path: str, output_qza: str,
                         has_header: bool) -> bool:
        """导入分类信息|Import taxonomy"""
        input_format = 'TSVTaxonomyFormat' if has_header else 'HeaderlessTSVTaxonomyFormat'
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'import',
            '--type', 'FeatureData[Taxonomy]',
            '--input-path', tax_path,
            '--output-path', output_qza,
            '--input-format', input_format,
        ])
        return self.cmd.run(cmd, f"导入分类信息|Import taxonomy ({input_format})")[0]

    def _extract_reads_region(self, ref_seqs_qza: str, output_qza: str) -> bool:
        """
        从参考序列中提取目标区域(in silico PCR)|Extract target region from reference (in silico PCR)

        仅16S使用|16S only
        """
        args = [
            'feature-classifier', 'extract-reads',
            '--i-sequences', ref_seqs_qza,
            '--p-f-primer', self.config.fwd_primer,
            '--p-r-primer', self.config.rev_primer,
            '--p-min-length', str(self.config.min_length),
            '--p-n-jobs', str(self.config.threads),
            '--o-reads', output_qza,
        ]
        if self.config.max_length > 0:
            args.extend(['--p-max-length', str(self.config.max_length)])
        cmd = build_conda_command(self.config.qiime_bin, args)
        return self.cmd.run(cmd, "提取目标区域|Extract target region (in silico PCR)")[0]

    def _fit_classifier(self, ref_reads_qza: str, ref_tax_qza: str,
                        output_qza: str) -> bool:
        """训练sklearn朴素贝叶斯分类器|Train sklearn naive Bayes classifier"""
        cmd = build_conda_command(self.config.qiime_bin, [
            'feature-classifier', 'fit-classifier-naive-bayes',
            '--i-reference-reads', ref_reads_qza,
            '--i-reference-taxonomy', ref_tax_qza,
            '--o-classifier', output_qza,
        ])
        return self.cmd.run(cmd, "训练分类器|Train classifier (fit-classifier-naive-bayes)")[0]

    # ------------------------------------------------------------------
    # 缓存与对外入口|Cache and public entry
    # ------------------------------------------------------------------

    def _db_tag(self) -> str:
        """数据库标签|Database tag"""
        return 'silva138.2' if self.config.is_16s else 'unite10_99'

    def _cache_path(self) -> str:
        """分类器缓存路径(按amplicon+引物+库区分)|Classifier cache path"""
        primer_hash = hashlib.md5(
            f"{self.config.fwd_primer}|{self.config.rev_primer}".encode()
        ).hexdigest()[:8]
        db_tag = self._db_tag()
        # 16S区分区域(引物相关),ITS为全长|16S keyed by region (primers), ITS full-length
        region = f"{self.config.fwd_primer[:6]}F-{self.config.rev_primer[:6]}R" if self.config.is_16s else 'fulllen'
        fname = f"classifier__{db_tag}__{self.config.amplicon}-{region}__{primer_hash}__sklearn.qza"
        return os.path.join(self.config.classifier_cache_dir, fname)

    def _write_cache_meta(self, classifier_qza: str) -> None:
        """写缓存元数据JSON(审计用)|Write cache metadata JSON (for audit)"""
        meta = {
            'amplicon': self.config.amplicon,
            'database_tag': self._db_tag(),
            'database_path': self.config.ref_database_path,
            'fwd_primer': self.config.fwd_primer,
            'rev_primer': self.config.rev_primer,
            'region_extracted': self.config.is_16s,
            'min_length': self.config.min_length,
            'max_length': self.config.max_length,
            'threads': self.config.threads,
        }
        meta_path = classifier_qza.replace('.qza', '.json')
        try:
            with open(meta_path, 'w') as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
        except Exception:
            pass

    def train(self) -> Tuple[str, str]:
        """
        训练分类器并写入缓存|Train classifier and write to cache

        Returns:
            (classifier_qza路径, 来源描述)|(classifier_qza path, source description)
        """
        cache_path = self._cache_path()
        os.makedirs(self.config.classifier_cache_dir, exist_ok=True)

        # 断点续传: 缓存已存在则跳过|Checkpoint: skip if cache exists
        if os.path.exists(cache_path):
            self.logger.info(f"分类器缓存命中,跳过训练|Classifier cache hit, skip training: {cache_path}")
            return cache_path, f"缓存|cached: {cache_path}"

        self.logger.info("=" * 60)
        self.logger.info(f"开始训练分类器|Start training classifier ({self.config.amplicon})")
        self.logger.info(f"数据库|Database: {self.config.ref_database_path}")
        self.logger.info(f"缓存输出|Cache output: {cache_path}")
        self.logger.warning(
            "分类器训练耗时较长(16S约30-60分钟),仅首次运行|"
            "Classifier training is slow (~30-60min for 16S), first run only"
        )

        # 训练工作目录(置于输出work下,便于排查)|Training work dir under output/work
        train_work = os.path.join(self.config.work_dir, 'classifier_train')
        os.makedirs(train_work, exist_ok=True)

        ref_seqs_qza = os.path.join(train_work, 'ref_seqs.qza')
        ref_tax_qza = os.path.join(train_work, 'ref_taxonomy.qza')

        if self.config.is_16s:
            # 16S: 解析SILVA → 导入 → 提取区域 → 训练|16S: parse SILVA → import → extract region → train
            seqs_fasta, tax_tsv = self.parse_silva(self.config.ref_database_path, train_work)

            if not self._import_sequences(seqs_fasta, ref_seqs_qza):
                raise RuntimeError("参考序列导入失败|Reference sequence import failed")
            if not self._import_taxonomy(tax_tsv, ref_tax_qza, has_header=False):
                raise RuntimeError("分类信息导入失败|Taxonomy import failed")

            ref_reads_qza = os.path.join(train_work, 'ref_reads_trimmed.qza')
            if not self._extract_reads_region(ref_seqs_qza, ref_reads_qza):
                raise RuntimeError("目标区域提取失败|Region extraction failed")
        else:
            # ITS: 解压UNITE → 导入 → 直接训练(不提取区域)|
            # ITS: extract UNITE → import → train directly (no region extraction)
            seqs_fasta, tax_txt = self.extract_unite(self.config.ref_database_path, train_work)

            if not self._import_sequences(seqs_fasta, ref_seqs_qza):
                raise RuntimeError("参考序列导入失败|Reference sequence import failed")
            if not self._import_taxonomy(tax_txt, ref_tax_qza, has_header=True):
                raise RuntimeError("分类信息导入失败|Taxonomy import failed")

            ref_reads_qza = ref_seqs_qza  # ITS直接用全长|ITS uses full-length directly

        if not self._fit_classifier(ref_reads_qza, ref_tax_qza, cache_path):
            raise RuntimeError("分类器训练失败|Classifier training failed")

        self._write_cache_meta(cache_path)
        self.logger.info(f"分类器训练完成并已缓存|Classifier trained and cached: {cache_path}")
        return cache_path, f"新训练|newly trained: {cache_path}"

    def get_or_train(self) -> Tuple[str, str]:
        """
        获取分类器: 优先用户指定→缓存→训练|Get classifier: user-specified → cache → train

        Returns:
            (classifier_qza路径, 来源描述)|(classifier_qza path, source description)
        """
        # 1. 用户显式指定|User-specified
        if self.config.classifier:
            self.logger.info(f"使用用户指定分类器|Using user-specified classifier: {self.config.classifier}")
            return self.config.classifier, f"用户指定|user-specified: {self.config.classifier}"

        # 2. 缓存或训练|Cache or train
        return self.train()
