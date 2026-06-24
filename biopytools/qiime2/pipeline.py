"""
QIIME2流程编排模块|QIIME2 Pipeline Orchestration Module

双端扩增子全流程: 导入→切引物→去噪(ASV)→OTU聚类→分类注释→系统发育→多样性+抽平→导出
|Paired-end amplicon pipeline: import→trim→denoise(ASV)→OTU→taxonomy→phylogeny→diversity+rarefaction→export
"""

import os
import shutil
from datetime import datetime
from typing import Callable, Dict, List, Optional, Tuple

from .config import Qiime2Config
from .utils import (
    CommandRunner, PairFinder, build_conda_command, generate_manifest,
    reverse_complement, generate_software_versions_yml, format_number
)


class Qiime2Pipeline:
    """QIIME2流程编排器|QIIME2 Pipeline Orchestrator"""

    def __init__(self, config: Qiime2Config, logger,
                 classifier_qza: str, classifier_source: str):
        self.config = config
        self.logger = logger
        self.classifier_qza = classifier_qza
        self.classifier_source = classifier_source

        # 命令执行器工作目录设为输出目录|Command runner working dir = output dir
        self.cmd = CommandRunner(logger, config.output_dir)

        # 流程中传递的中间产物路径|Intermediate artifact paths threaded through pipeline
        self.demux_qza: Optional[str] = None          # 导入后的双端序列|imported paired sequences
        self.denoise_input: Optional[str] = None      # 去噪输入(demux或trimmed)|denoise input
        self.feature_table: Optional[str] = None      # 当前活动特征表(ASV或OTU)|active feature table
        self.rep_seqs: Optional[str] = None           # 当前活动代表序列|active representative seqs
        self.rooted_tree: Optional[str] = None        # 有根树|rooted tree
        self.sampling_depth_used: int = 0             # 实际使用的抽样深度|sampling depth used

    # ------------------------------------------------------------------
    # 断点续传辅助|Checkpoint helpers
    # ------------------------------------------------------------------

    def _is_completed(self, output_file: str) -> bool:
        """检查步骤是否已完成(输出文件存在性)|Check if step is done (output existence)"""
        return os.path.exists(output_file)

    def _run_step(self, step_name: str, output_file: str,
                  run_func: Callable, *args, **kwargs) -> bool:
        """带断点续传的步骤执行器|Step runner with checkpoint resume (§10.2)"""
        if self._is_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: {step_name}")
            return True

        self.logger.info(f"开始步骤|Starting step: {step_name}")
        return run_func(*args, **kwargs)

    def _ensure_metadata(self, pairs: Dict[str, Tuple[str, str]]) -> str:
        """确保有元数据文件(用户提供或自动生成)|Ensure metadata file (user or auto-generated)"""
        if self.config.metadata_file and os.path.exists(self.config.metadata_file):
            return self.config.metadata_file

        # 生成最小元数据(满足core-metrics对--m-metadata-file的要求)|
        # Generate minimal metadata (satisfies core-metrics --m-metadata-file requirement)
        meta_path = os.path.join(self.config.info_dir, 'sample_metadata.tsv')
        with open(meta_path, 'w') as f:
            f.write('sample-id\tgroup\n')
            f.write('#q2:types\tcategorical\n')
            for sample_name in sorted(pairs.keys()):
                f.write(f'{sample_name}\tall\n')
        self.logger.info(f"自动生成样品元数据|Auto-generated sample metadata: {meta_path}")
        return meta_path

    # ------------------------------------------------------------------
    # 步骤1: 导入|Step 1: Import
    # ------------------------------------------------------------------

    def step_import(self, pairs: Dict[str, Tuple[str, str]]) -> bool:
        """导入双端FASTQ|Import paired-end FASTQ"""
        manifest_path = os.path.join(self.config.info_dir, 'sample_manifest.tsv')
        generate_manifest(pairs, manifest_path)
        self.logger.info(f"样品清单|Sample manifest: {len(pairs)} 个样品|samples")
        for sample_name in sorted(pairs.keys()):
            self.logger.info(f"  {sample_name}")

        self.demux_qza = os.path.join(self.config.import_dir, 'demux.qza')

        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'import',
            '--type', 'SampleData[PairedEndSequencesWithQuality]',
            '--input-path', manifest_path,
            '--output-path', self.demux_qza,
            '--input-format', 'PairedEndFastqManifestPhred33V2',
        ])
        if not self.cmd.run(cmd, "导入双端FASTQ|Import paired-end FASTQ"):
            return False

        # 质量汇总(可视化,失败不致命)|Quality summary (visualization, non-fatal)
        demux_qzv = os.path.join(self.config.import_dir, 'demux.qzv')
        cmd = build_conda_command(self.config.qiime_bin, [
            'demux', 'summarize',
            '--i-data', self.demux_qza,
            '--o-visualization', demux_qzv,
        ])
        self.cmd.run(cmd, "序列质量汇总|Demux quality summary")

        self.denoise_input = self.demux_qza
        return True

    # ------------------------------------------------------------------
    # 步骤2: 切引物|Step 2: Trim primers
    # ------------------------------------------------------------------

    def step_trim(self) -> bool:
        """切除引物(cutadapt)|Trim primers (cutadapt)"""
        trimmed_qza = os.path.join(self.config.trim_dir, 'demux_trimmed.qza')

        # --p-front-r需传反向引物的反向互补序列|--p-front-r needs reverse complement of reverse primer
        rev_rc = reverse_complement(self.config.rev_primer)
        self.logger.info(f"反向引物反向互补|Reverse primer RC: {self.config.rev_primer} -> {rev_rc}")

        cmd = build_conda_command(self.config.qiime_bin, [
            'cutadapt', 'trim-paired',
            '--i-demultiplexed-sequences', self.denoise_input,
            '--p-front-f', self.config.fwd_primer,
            '--p-front-r', rev_rc,
            '--p-discard-untrimmed',
            '--p-error-rate', '0.1',
            '--p-cores', str(self.config.threads),
            '--o-trimmed-sequences', trimmed_qza,
        ])
        if not self.cmd.run(cmd, "切除引物|Trim primers (cutadapt)"):
            return False

        self.denoise_input = trimmed_qza
        return True

    # ------------------------------------------------------------------
    # 步骤3: 去噪(ASV)|Step 3: Denoise (ASV)
    # ------------------------------------------------------------------

    def step_denoise(self) -> bool:
        """DADA2双端去噪,产出ASV表/代表序列/统计|DADA2 paired denoise → ASV table/seqs/stats"""
        table_qza = os.path.join(self.config.denoise_dir, 'table.qza')
        rep_seqs_qza = os.path.join(self.config.denoise_dir, 'rep_seqs.qza')
        stats_qza = os.path.join(self.config.denoise_dir, 'denoising_stats.qza')

        cmd = build_conda_command(self.config.qiime_bin, [
            'dada2', 'denoise-paired',
            '--i-demultiplexed-seqs', self.denoise_input,
            '--p-trunc-len-f', str(self.config.trunc_len_f),
            '--p-trunc-len-r', str(self.config.trunc_len_r),
            '--p-trim-left-f', str(self.config.trim_left_f),
            '--p-trim-left-r', str(self.config.trim_left_r),
            '--p-n-threads', str(self.config.threads),
            '--p-chimera-method', 'consensus',
            '--o-table', table_qza,
            '--o-representative-sequences', rep_seqs_qza,
            '--o-denoising-stats', stats_qza,
        ])
        if not self.cmd.run(cmd, "DADA2双端去噪(ASV)|DADA2 paired denoise (ASV)"):
            return False

        self.feature_table = table_qza
        self.rep_seqs = rep_seqs_qza

        # 特征表汇总(失败不致命)|Feature table summary (non-fatal)
        table_qzv = os.path.join(self.config.denoise_dir, 'table.qzv')
        cmd = build_conda_command(self.config.qiime_bin, [
            'feature-table', 'summarize',
            '--i-table', table_qza,
            '--o-visualization', table_qzv,
        ])
        self.cmd.run(cmd, "特征表汇总|Feature table summary")

        # 代表序列制表|Tabulate representative sequences
        rep_seqs_qzv = os.path.join(self.config.denoise_dir, 'rep_seqs.qzv')
        cmd = build_conda_command(self.config.qiime_bin, [
            'feature-table', 'tabulate-seqs',
            '--i-data', rep_seqs_qza,
            '--o-visualization', rep_seqs_qzv,
        ])
        self.cmd.run(cmd, "代表序列制表|Tabulate representative sequences")

        return True

    # ------------------------------------------------------------------
    # 步骤4: OTU聚类(可选)|Step 4: OTU clustering (optional)
    # ------------------------------------------------------------------

    def step_otu(self) -> bool:
        """vsearch从头聚类OTU,替换下游使用的特征表/序列|vsearch de novo OTU clustering"""
        otu_table_qza = os.path.join(self.config.otu_dir, 'otu_table.qza')
        otu_seqs_qza = os.path.join(self.config.otu_dir, 'otu_seqs.qza')

        cmd = build_conda_command(self.config.qiime_bin, [
            'vsearch', 'cluster-features-de-novo',
            '--i-sequences', self.rep_seqs,
            '--i-table', self.feature_table,
            '--p-perc-identity', str(self.config.perc_identity),
            '--p-strand', 'plus',
            '--p-threads', str(self.config.threads),
            '--o-clustered-table', otu_table_qza,
            '--o-clustered-sequences', otu_seqs_qza,
        ])
        if not self.cmd.run(cmd, f"vsearch OTU聚类({self.config.perc_identity})|vsearch OTU clustering"):
            return False

        # 下游改用OTU表/序列|Downstream uses OTU table/seqs
        self.feature_table = otu_table_qza
        self.rep_seqs = otu_seqs_qza
        return True

    # ------------------------------------------------------------------
    # 步骤5: 分类注释|Step 5: Taxonomy classification
    # ------------------------------------------------------------------

    def step_taxonomy(self) -> bool:
        """sklearn分类注释|sklearn taxonomy classification"""
        taxonomy_qza = os.path.join(self.config.taxonomy_dir, 'taxonomy.qza')

        cmd = build_conda_command(self.config.qiime_bin, [
            'feature-classifier', 'classify-sklearn',
            '--i-classifier', self.classifier_qza,
            '--i-reads', self.rep_seqs,
            '--p-confidence', str(self.config.confidence),
            '--p-n-jobs', str(self.config.threads),
            '--o-classification', taxonomy_qza,
        ])
        if not self.cmd.run(cmd, "sklearn分类注释|sklearn classification"):
            return False

        # 分类制表|Tabulate taxonomy
        taxonomy_qzv = os.path.join(self.config.taxonomy_dir, 'taxonomy.qzv')
        cmd = build_conda_command(self.config.qiime_bin, [
            'metadata', 'tabulate',
            '--m-input-file', taxonomy_qza,
            '--o-visualization', taxonomy_qzv,
        ])
        self.cmd.run(cmd, "分类注释制表|Tabulate taxonomy")

        return True

    # ------------------------------------------------------------------
    # 步骤6: 系统发育建树|Step 6: Phylogeny
    # ------------------------------------------------------------------

    def step_phylogeny(self) -> bool:
        """mafft比对+fasttree建树|mafft alignment + fasttree tree"""
        aligned_qza = os.path.join(self.config.phylogeny_dir, 'aligned_rep_seqs.qza')
        masked_qza = os.path.join(self.config.phylogeny_dir, 'masked_aligned_rep_seqs.qza')
        unrooted_qza = os.path.join(self.config.phylogeny_dir, 'unrooted_tree.qza')
        rooted_qza = os.path.join(self.config.phylogeny_dir, 'rooted_tree.qza')

        cmd = build_conda_command(self.config.qiime_bin, [
            'phylogeny', 'align-to-tree-mafft-fasttree',
            '--i-sequences', self.rep_seqs,
            '--p-n-threads', str(self.config.threads),
            '--o-alignment', aligned_qza,
            '--o-masked-alignment', masked_qza,
            '--o-tree', unrooted_qza,
            '--o-rooted-tree', rooted_qza,
        ])
        if not self.cmd.run(cmd, "mafft比对+fasttree建树|mafft + fasttree phylogeny"):
            return False

        self.rooted_tree = rooted_qza
        return True

    # ------------------------------------------------------------------
    # 抽样深度自动计算|Auto sampling depth
    # ------------------------------------------------------------------

    @staticmethod
    def _percentile(values: List[float], p: float) -> float:
        """计算百分位(不依赖numpy)|Compute percentile (no numpy dependency)"""
        if not values:
            return 0.0
        s = sorted(values)
        if len(s) == 1:
            return s[0]
        rank = (p / 100.0) * (len(s) - 1)
        lo = int(rank)
        hi = min(lo + 1, len(s) - 1)
        frac = rank - lo
        return s[lo] * (1 - frac) + s[hi] * frac

    def _auto_sampling_depth(self) -> Tuple[int, Dict[str, int]]:
        """导出特征表→biom→tsv,统计每样本reads总数,取第10百分位|
        Export table→biom→tsv, compute per-sample totals, take 10th percentile"""
        biom_path = os.path.join(self.config.work_dir, 'feature_table.biom')
        tsv_path = os.path.join(self.config.work_dir, 'feature_table.tsv')

        # 导出为biom|Export to biom
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'export',
            '--input-path', self.feature_table,
            '--output-path', biom_path,
        ])
        if not self.cmd.run(cmd, "导出特征表(biom)|Export feature table (biom)"):
            return 1, {}

        # biom转tsv|biom to tsv
        biom_bin = os.path.join(os.path.dirname(self.config.qiime_bin), 'biom')
        cmd = build_conda_command(biom_bin, [
            'convert', '-i', biom_path, '-o', tsv_path, '--to-tsv',
        ])
        if not self.cmd.run(cmd, "biom转tsv|biom to tsv"):
            return 1, {}

        # 解析tsv,统计每样本总数|Parse tsv, compute per-sample totals
        per_sample: Dict[str, int] = {}
        sample_names: List[str] = []
        with open(tsv_path) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('#OTU ID'):
                    sample_names = line.split('\t')[1:]
                    continue
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 2:
                    continue
                counts = fields[1:]
                for i, name in enumerate(sample_names):
                    if i >= len(counts):
                        break
                    try:
                        per_sample[name] = per_sample.get(name, 0) + int(float(counts[i]))
                    except ValueError:
                        pass

        if not per_sample:
            return 1, {}

        totals = list(per_sample.values())
        depth = max(1, int(self._percentile(totals, 10)))

        self.logger.info("每样本reads分布|Per-sample reads distribution:")
        for name in sorted(per_sample.keys()):
            self.logger.info(f"  {name}: {format_number(per_sample[name])}")
        self.logger.info(
            f"自动抽样深度(第10百分位)|Auto sampling depth (10th percentile): {depth}"
        )
        return depth, per_sample

    # ------------------------------------------------------------------
    # 步骤7: 多样性+抽平|Step 7: Diversity + rarefaction
    # ------------------------------------------------------------------

    def step_diversity(self, metadata_path: str) -> bool:
        """核心多样性指标(alpha/beta/抽平)|Core diversity metrics (alpha/beta/rarefaction)"""
        # 确定抽样深度|Determine sampling depth
        if self.config.sampling_depth > 0:
            self.sampling_depth_used = self.config.sampling_depth
            self.logger.info(f"使用用户指定抽样深度|Using user sampling depth: {self.sampling_depth_used}")
        else:
            depth, _ = self._auto_sampling_depth()
            self.sampling_depth_used = depth

        if self.sampling_depth_used < 1:
            self.logger.error("抽样深度无效|Invalid sampling depth")
            return False

        if self.config.skip_phylogeny or self.rooted_tree is None:
            # ITS或跳过系统发育: 非系统发育核心指标|ITS or skip phylogeny: non-phylogenetic core metrics
            cm_dir = os.path.join(self.config.diversity_dir, 'core_metrics')
            cmd = build_conda_command(self.config.qiime_bin, [
                'diversity', 'core-metrics',
                '--i-table', self.feature_table,
                '--p-sampling-depth', str(self.sampling_depth_used),
                '--m-metadata-file', metadata_path,
                '--output-dir', cm_dir,
                '--p-n-jobs', str(self.config.threads),
            ])
        else:
            # 16S: 系统发育核心指标(含UniFrac/Faith)|16S: phylogenetic core metrics (UniFrac/Faith)
            cm_dir = os.path.join(self.config.diversity_dir, 'core_metrics_phylogenetic')
            cmd = build_conda_command(self.config.qiime_bin, [
                'diversity', 'core-metrics-phylogenetic',
                '--i-table', self.feature_table,
                '--i-phylogeny', self.rooted_tree,
                '--p-sampling-depth', str(self.sampling_depth_used),
                '--m-metadata-file', metadata_path,
                '--output-dir', cm_dir,
                '--p-n-jobs-or-threads', str(self.config.threads),
            ])
        if not self.cmd.run(cmd, "核心多样性指标|Core diversity metrics"):
            return False

        # alpha稀释曲线|Alpha rarefaction curve
        alpha_qzv = os.path.join(self.config.diversity_dir, 'alpha_rarefaction.qzv')
        max_depth = max(self.sampling_depth_used * 2, self.sampling_depth_used + 1)
        args = [
            'diversity', 'alpha-rarefaction',
            '--i-table', self.feature_table,
            '--p-min-depth', '1',
            '--p-max-depth', str(max_depth),
            '--p-steps', '10',
            '--m-metadata-file', metadata_path,
            '--o-visualization', alpha_qzv,
        ]
        if self.rooted_tree is not None:
            args.extend(['--i-phylogeny', self.rooted_tree])
        cmd = build_conda_command(self.config.qiime_bin, args)
        self.cmd.run(cmd, "alpha稀释曲线|Alpha rarefaction")

        return True

    # ------------------------------------------------------------------
    # 步骤8: 导出|Step 8: Export
    # ------------------------------------------------------------------

    def step_export(self) -> bool:
        """导出用户最终产物(序列/丰度表/分类)|Export user-facing products"""
        prefix = 'asv' if self.config.method == 'asv' else 'otu'

        # 代表序列→FASTA|Representative sequences → FASTA
        seqs_fasta = os.path.join(self.config.export_dir, f'{prefix}_sequences.fasta')
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'export',
            '--input-path', self.rep_seqs,
            '--output-path', seqs_fasta,
            '--output-format', 'DNAFASTAFormat',
        ])
        if not self.cmd.run(cmd, f"导出{prefix.upper()}序列|Export {prefix.upper()} sequences"):
            return False

        # 特征表→BIOM|Feature table → BIOM
        biom_out = os.path.join(self.config.export_dir, 'feature_table.biom')
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'export',
            '--input-path', self.feature_table,
            '--output-path', biom_out,
        ])
        if not self.cmd.run(cmd, "导出丰度表(biom)|Export feature table (biom)"):
            return False

        # BIOM→TSV|BIOM → TSV
        tsv_out = os.path.join(self.config.export_dir, 'feature_table.tsv')
        biom_bin = os.path.join(os.path.dirname(self.config.qiime_bin), 'biom')
        cmd = build_conda_command(biom_bin, [
            'convert', '-i', biom_out, '-o', tsv_out, '--to-tsv',
        ])
        if not self.cmd.run(cmd, "biom转tsv|biom to tsv"):
            return False

        # 分类→TSV|Taxonomy → TSV
        taxonomy_qza = os.path.join(self.config.taxonomy_dir, 'taxonomy.qza')
        taxonomy_out = os.path.join(self.config.export_dir, 'taxonomy.tsv')
        cmd = build_conda_command(self.config.qiime_bin, [
            'tools', 'export',
            '--input-path', taxonomy_qza,
            '--output-path', taxonomy_out,
            '--output-format', 'TSVTaxonomyFormat',
        ])
        self.cmd.run(cmd, "导出分类表|Export taxonomy (tsv)")

        # 合并分类进丰度表|Merge taxonomy into feature table
        merged_out = os.path.join(self.config.export_dir, 'feature_table_w_taxonomy.tsv')
        self._merge_table_taxonomy(tsv_out, taxonomy_out, merged_out)

        return True

    def _merge_table_taxonomy(self, table_tsv: str, taxonomy_tsv: str,
                              out_path: str) -> None:
        """合并特征表与分类信息|Merge feature table with taxonomy"""
        # 读取分类: feature_id → taxon|Read taxonomy: feature_id → taxon
        taxonomy: Dict[str, str] = {}
        if os.path.exists(taxonomy_tsv):
            with open(taxonomy_tsv) as f:
                header_seen = False
                for line in f:
                    line = line.rstrip('\n')
                    if not line or line.startswith('#q2:types'):
                        continue
                    if not header_seen:
                        header_seen = True  # 跳过表头|skip header
                        continue
                    fields = line.split('\t')
                    if len(fields) >= 2:
                        taxonomy[fields[0]] = fields[1]

        # 读取特征表并合并|Read feature table and merge
        try:
            with open(table_tsv) as f, open(out_path, 'w') as out:
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('#OTU ID'):
                        out.write(f'{line}\ttaxonomy\n')
                        continue
                    if line.startswith('#') or not line:
                        out.write(line + '\n')
                        continue
                    fields = line.split('\t')
                    fid = fields[0]
                    taxon = taxonomy.get(fid, 'Unassigned')
                    out.write(f'{line}\t{taxon}\n')
            self.logger.info(f"合并表已生成|Merged table generated: {out_path}")
        except Exception as e:
            self.logger.warning(f"合并分类失败|Failed to merge taxonomy: {e}")

    # ------------------------------------------------------------------
    # 主流程|Main pipeline
    # ------------------------------------------------------------------

    def _reset_step_dirs_if_forced(self) -> None:
        """force模式下清空步骤目录(保留日志与信息)|Clear step dirs if force (keep logs/info)"""
        if not self.config.force:
            return
        self.logger.warning("force模式: 清空已有步骤输出|force: clearing existing step outputs")
        step_dirs = [self.config.import_dir, self.config.trim_dir,
                     self.config.denoise_dir, self.config.otu_dir,
                     self.config.taxonomy_dir, self.config.phylogeny_dir,
                     self.config.diversity_dir, self.config.export_dir,
                     self.config.work_dir]
        for d in step_dirs:
            if os.path.exists(d):
                shutil.rmtree(d)
                os.makedirs(d, exist_ok=True)

    def run_pipeline(self) -> bool:
        """运行完整流程|Run complete pipeline"""
        start_time = datetime.now()

        self.logger.info("=" * 60)
        self.logger.info("QIIME2流程开始|QIIME2 pipeline started")
        self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"扩增子|Amplicon: {self.config.amplicon} | 方法|Method: {self.config.method}")
        self.logger.info(f"分类器|Classifier: {self.classifier_source}")

        # 检测双端样品|Detect paired-end samples
        finder = PairFinder(self.config.input_dir, self.config.r1_suffix, self.config.r2_suffix)
        pairs = finder.find_pairs()
        if not pairs:
            self.logger.error(
                f"未检测到任何双端样品|No paired-end samples detected "
                f"(r1_suffix={self.config.r1_suffix}, r2_suffix={self.config.r2_suffix})"
            )
            return False
        self.logger.info(f"检测到样品|Detected samples: {len(pairs)}")

        self._reset_step_dirs_if_forced()
        metadata_path = self._ensure_metadata(pairs)

        # 步骤1: 导入(标记=导入的demux.qza)|Step 1: Import (marker=imported demux.qza)
        if not self._run_step("01_import|Import",
                              os.path.join(self.config.import_dir, 'demux.qza'),
                              self.step_import, pairs):
            return False

        # 步骤2: 切引物(可选,标记=trimmed.qza)|Step 2: Trim (optional, marker=trimmed.qza)
        if not self.config.skip_cutadapt:
            if not self._run_step("02_trim|Trim",
                                  os.path.join(self.config.trim_dir, 'demux_trimmed.qza'),
                                  self.step_trim):
                return False
        else:
            self.logger.info("跳过引物切除(--skip-cutadapt)|Skipping primer trimming")

        # 步骤3: 去噪(ASV)|Step 3: Denoise
        if not self._run_step("03_denoise|Denoise",
                              os.path.join(self.config.denoise_dir, 'table.qza'),
                              self.step_denoise):
            return False

        # 步骤4: OTU聚类(可选)|Step 4: OTU (optional)
        if self.config.method == 'otu':
            if not self._run_step("04_otu|OTU",
                                  os.path.join(self.config.otu_dir, 'otu_table.qza'),
                                  self.step_otu):
                return False

        # 步骤5: 分类注释|Step 5: Taxonomy
        if not self._run_step("05_taxonomy|Taxonomy",
                              os.path.join(self.config.taxonomy_dir, 'taxonomy.qza'),
                              self.step_taxonomy):
            return False

        # 步骤6: 系统发育(可选)|Step 6: Phylogeny (optional)
        if not self.config.skip_phylogeny:
            if not self._run_step("06_phylogeny|Phylogeny",
                                  os.path.join(self.config.phylogeny_dir, 'rooted_tree.qza'),
                                  self.step_phylogeny):
                return False
        else:
            self.logger.info("跳过系统发育建树(ITS或--skip-phylogeny)|Skipping phylogeny")

        # 步骤7: 多样性+抽平(标记=core_metrics输出目录)|Step 7: Diversity (marker=core_metrics dir)
        cm_subdir = 'core_metrics' if (self.config.skip_phylogeny or self.rooted_tree is None) \
                    else 'core_metrics_phylogenetic'
        if not self._run_step("07_diversity|Diversity",
                              os.path.join(self.config.diversity_dir, cm_subdir),
                              self.step_diversity, metadata_path):
            return False

        # 步骤8: 导出(标记=序列fasta)|Step 8: Export (marker=sequences fasta)
        export_prefix = 'asv' if self.config.method == 'asv' else 'otu'
        if not self._run_step("08_export|Export",
                              os.path.join(self.config.export_dir, f'{export_prefix}_sequences.fasta'),
                              self.step_export):
            return False

        # 生成版本信息|Generate version info
        generate_software_versions_yml(
            self.config.output_dir, self.config,
            self.classifier_source, self.sampling_depth_used, start_time
        )

        elapsed = (datetime.now() - start_time).total_seconds()
        self.logger.info("=" * 60)
        self.logger.info("QIIME2流程完成|QIIME2 pipeline completed")
        self.logger.info(f"总耗时|Total time: {elapsed:.1f}s")
        self.logger.info(f"抽样深度|Sampling depth: {self.sampling_depth_used}")
        self.logger.info(f"最终产物|Final products: {self.config.export_dir}/")
        self.logger.info(f"日志|Log: {self.config.log_dir}/")

        return True
