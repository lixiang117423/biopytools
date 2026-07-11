"""疫霉菌基因组注释流程|Oomycete Genome Annotation Pipeline

复刻 T2T 证据驱动 Augustus 手工流程(Phase1 + Phase2 全量):
01 重复屏蔽 -> 02 RNA-seq 比对 -> 03 三代转录本比对(可选) ->
04 蛋白比对(可选) -> 05 hints -> 06 训练集 -> 07 Augustus 预测 -> 08 LTR(可选,正交)
|Replicates the T2T evidence-driven Augustus pipeline (Phase 1 + 2):
01 repeat mask -> 02 RNA-seq align -> 03 long-read align (opt) ->
04 protein align (opt) -> 05 hints -> 06 training set -> 07 Augustus -> 08 LTR (opt)
"""

import glob
import os
import shutil
from pathlib import Path

from .config import OomyceteAnnoConfig
from .utils import (
    CommandRunner,
    OomyceteAnnoLogger,
    build_conda_command,
    build_genemark_command,
    find_rnaseq_pairs,
    miniprot_gff_to_protein_hints,
)


class OomyceteAnnoRunner:
    """疫霉菌注释流程运行器|Oomycete annotation pipeline runner."""

    def __init__(self, config: OomyceteAnnoConfig, logger=None):
        if not isinstance(config, OomyceteAnnoConfig):
            raise TypeError("config 必须是 OomyceteAnnoConfig|config must be OomyceteAnnoConfig")
        self.config = config

        if not logger:
            logger = OomyceteAnnoLogger(config.log_dir).get_logger()
        self.logger = logger
        self.cmd = CommandRunner(logger, config.working_dir)

    def _ascii_workdir(self, tag: str) -> Path:
        """~/tmp 下 ASCII 临时工作目录(避开中文路径), 每基因组独立且持久(供续传)
        |ASCII temp workdir under ~/tmp; per-genome, persistent for resume.

        GeneMark probuild / Augustus 等 C 二进制及 gmes_petap.pl 的 abs_path 会把中文
        路径解成乱码。把输入真实拷贝进此 ASCII 目录, cwd 设为此目录运行, 完成后拷回产物。
        |C binaries and gmes_petap.pl's abs_path mangle CJK paths. Copy inputs into this
        ASCII dir, run with cwd=here, then copy outputs back to the real output dir.
        """
        import hashlib
        d = Path.home() / "tmp" / f"oomycete_{tag}_{hashlib.md5(self.config.genome.encode()).hexdigest()[:8]}"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def _copy_in(self, src: str, dst_path: Path) -> Path:
        """把 src 真实拷贝到 dst_path(未拷过或大小不符才拷); 返回 dst_path
        |Real-copy src to dst_path if absent/size-mismatch; return dst_path."""
        src = os.path.abspath(src)
        if not dst_path.exists() or os.path.getsize(dst_path) != os.path.getsize(src):
            shutil.copy2(src, str(dst_path))
        return dst_path

    # ----------------------------------------------------------
    # 主流程|Main pipeline
    # ----------------------------------------------------------
    def run(self) -> bool:
        """运行完整流程(graceful degradation)|Run full pipeline."""
        c = self.config
        self.logger.info("=" * 70)
        self.logger.info("开始疫霉菌基因组注释|Starting oomycete genome annotation")
        self.logger.info("=" * 70)
        self.logger.info(f"基因组|Genome: {c.genome}")
        self.logger.info(f"物种|Species: {c.species}")
        self.logger.info(
            f"证据|Evidence: RNA-seq={bool(c.rnaseq_dirs)} 三代|iso-seq={bool(c.isoseq)} "
            f"蛋白|protein={bool(c.prot_seq)}"
        )

        try:
            # 01 重复屏蔽|repeat masking
            if c.skip_repeat:
                self.logger.info("跳过重复屏蔽|Skipping repeat masking")
                masked_genome = c.genome
            else:
                masked_genome = self._step_repeat_masking()
                if not masked_genome:
                    self.logger.warning("重复屏蔽失败, 使用原基因组|Masking failed, using raw genome")
                    masked_genome = c.genome

            # 02 RNA-seq 比对 -> BAM|RNA-seq align
            rna_bam = None
            if c.rnaseq_dirs and not c.skip_rna:
                rna_bam = self._step_rna_align(masked_genome)

            # 03 三代转录本比对 -> GFF3(供 TransDecoder 训练集)|long-read align
            iso_gff = None
            if c.isoseq and not c.skip_iso:
                iso_gff = self._step_iso_align(masked_genome)

            # 04 蛋白比对 -> protein hints|protein align
            protein_hints = None
            if c.prot_seq and not c.skip_protein:
                protein_hints = self._step_protein_hints(masked_genome)

            # 05 hints|hints
            intron_gff = None
            if rna_bam:
                intron_gff = self._step_intron_hints(rna_bam)
            hintsfile = self._step_merge_hints(intron_gff, protein_hints)

            # 06 训练集(三代->TransDecoder, 否则 GeneMark ES/ET/EP+)|training set
            training_gtf = self._step_training_set(masked_genome, intron_gff, iso_gff)
            if not training_gtf:
                self.logger.error("训练集生成失败|Training set generation failed")
                return False

            # 07 Augustus 训练 + 预测|Augustus train + predict
            gff = self._step_augustus(masked_genome, training_gtf, hintsfile)
            if not gff:
                self.logger.error("Augustus 预测失败|Augustus prediction failed")
                return False

            self.logger.info("=" * 70)
            self.logger.info(f"注释完成|Annotation completed: {gff}")
            self.logger.info("=" * 70)
            self._write_versions(gff)

            # 08 LTR(正交 TE, 失败不阻断)|LTR (orthogonal, non-fatal)
            if not c.skip_ltr:
                try:
                    self._step_ltr(masked_genome)
                except Exception as e:
                    self.logger.warning(f"LTR 注释失败(不影响基因注释结果)|LTR failed (non-fatal): {e}")
            return True

        except Exception as e:
            self.logger.error(f"流程失败|Pipeline failed: {e}")
            raise

    # ----------------------------------------------------------
    # 步骤1: 重复屏蔽|Repeat masking
    # ----------------------------------------------------------
    def _step_repeat_masking(self) -> str:
        """RepeatModeler + RepeatMasker 软屏蔽|Soft-mask repeats."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤1: 重复序列屏蔽|Step 1: Repeat masking")
        self.logger.info("-" * 70)

        masked_genome = self._find_masked_genome()
        if masked_genome:
            self.logger.info(f"跳过(已完成)|Skipping (done): {masked_genome}")
            return masked_genome

        genome_abs = os.path.abspath(c.genome)
        repeat_dir_abs = os.path.abspath(c.repeat_dir)
        db_name = f"{c.species}_db"

        cmd = build_conda_command(c.build_database_bin, ["-name", db_name, genome_abs])
        if not self.cmd.run(cmd, "构建 RepeatModeler 数据库|Build RepeatModeler DB", cwd=c.repeat_dir):
            raise RuntimeError("BuildDatabase 失败|BuildDatabase failed")

        cmd = build_conda_command(
            c.repeatmodeler_bin, ["-database", db_name, "-threads", str(c.threads)]
        )
        if not self.cmd.run(cmd, "运行 RepeatModeler|Run RepeatModeler", cwd=c.repeat_dir):
            raise RuntimeError("RepeatModeler 失败|RepeatModeler failed")

        consensi_list = glob.glob(os.path.join(repeat_dir_abs, "RM_*", "consensi.fa.classified"))
        if not consensi_list:
            consensi_list = glob.glob(os.path.join(repeat_dir_abs, "**", "consensi.fa.classified"), recursive=True)
        if not consensi_list:
            raise RuntimeError("找不到 RepeatModeler 的 consensi.fa.classified|consensi not found")
        consensi_abs = os.path.abspath(consensi_list[0])
        self.logger.info(f"重复序列库|Repeat library: {consensi_abs}")

        mask_opt = ["-xsmall"] if c.soft_masking else []
        cmd = build_conda_command(
            c.repeatmasker_bin,
            ["-lib", consensi_abs, "-dir", repeat_dir_abs, "-pa", str(c.threads)]
            + mask_opt + [genome_abs],
        )
        if not self.cmd.run(cmd, "运行 RepeatMasker|Run RepeatMasker"):
            raise RuntimeError("RepeatMasker 失败|RepeatMasker failed")

        masked_genome = self._find_masked_genome()
        if not masked_genome:
            raise RuntimeError(f"找不到屏蔽后的基因组|Masked genome not found in {repeat_dir_abs}")
        self.logger.info(f"屏蔽完成|Masking completed: {masked_genome}")
        return masked_genome

    def _find_masked_genome(self) -> str:
        """定位屏蔽后的基因组|Locate masked genome."""
        c = self.config
        genome_name = os.path.basename(c.genome)
        for cand in (
            os.path.join(c.repeat_dir, f"{genome_name}.masked"),
            os.path.join(c.repeat_dir, f"{os.path.splitext(genome_name)[0]}.masked.fa"),
        ):
            if os.path.exists(cand):
                return cand
        masked_files = glob.glob(os.path.join(os.path.abspath(c.repeat_dir), "*.masked"))
        return masked_files[0] if masked_files else ""

    # ----------------------------------------------------------
    # 步骤2: RNA-seq 比对|RNA-seq alignment (HISAT2)
    # ----------------------------------------------------------
    def _step_rna_align(self, genome: str) -> str:
        """HISAT2 比对 -> 排序+索引 BAM|HISAT2 align to sorted+indexed BAM."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤2: RNA-seq 比对|Step 2: RNA-seq alignment (HISAT2)")
        self.logger.info("-" * 70)

        sorted_bam = os.path.join(c.rna_align_dir, "rnaseq.sorted.bam")
        if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
            self.logger.info(f"跳过(已完成)|Skipping (done): {sorted_bam}")
            return sorted_bam

        genome_abs = os.path.abspath(genome)
        index_prefix = os.path.join(c.rna_align_dir, c.species)

        cmd = build_conda_command(c.hisat2_build_bin, ["-p", str(c.threads), genome_abs, index_prefix])
        if not self.cmd.run(cmd, "构建 HISAT2 索引|Build HISAT2 index"):
            raise RuntimeError("HISAT2 索引构建失败|HISAT2 index build failed")

        pairs = find_rnaseq_pairs(c.rnaseq_dirs, c.read1_pattern, c.read2_pattern, self.logger)
        if not pairs:
            raise RuntimeError("找不到 RNA-seq 配对文件|No RNA-seq pairs found")
        self.logger.info(f"找到 {len(pairs)} 对 RNA-seq|Found {len(pairs)} RNA-seq pairs")

        if len(pairs) > 1:
            r1_in = os.path.join(c.rna_align_dir, "merged_R1.fq.gz")
            r2_in = os.path.join(c.rna_align_dir, "merged_R2.fq.gz")
            if not (os.path.exists(r1_in) and os.path.exists(r2_in)):
                self._cat_files([p[0] for p in pairs], r1_in, "合并 R1|Merge R1")
                self._cat_files([p[1] for p in pairs], r2_in, "合并 R2|Merge R2")
        else:
            r1_in, r2_in = pairs[0]

        sam_file = os.path.join(c.rna_align_dir, "rnaseq.sam")
        hisat_args = ["-x", index_prefix, "-1", r1_in, "-2", r2_in,
                      "-p", str(c.threads), "-S", sam_file]
        if c.rna_strandness:
            hisat_args.extend(["--rna-strandness", c.rna_strandness])
        cmd = build_conda_command(c.hisat2_bin, hisat_args)
        if not self.cmd.run(cmd, "HISAT2 比对|HISAT2 alignment"):
            raise RuntimeError("HISAT2 比对失败|HISAT2 alignment failed")

        cmd = build_conda_command(c.samtools_bin, ["sort", "-@", str(c.threads), "-o", sorted_bam, sam_file])
        if not self.cmd.run(cmd, "排序 BAM|Sort BAM"):
            raise RuntimeError("BAM 排序失败|BAM sort failed")
        cmd = build_conda_command(c.samtools_bin, ["index", sorted_bam])
        if not self.cmd.run(cmd, "索引 BAM|Index BAM"):
            raise RuntimeError("BAM 索引失败|BAM index failed")

        if os.path.exists(sam_file):
            os.remove(sam_file)
        self.logger.info(f"RNA-seq 比对完成|RNA-seq alignment done: {sorted_bam}")
        return sorted_bam

    def _cat_files(self, srcs, dst, desc):
        """合并文件(直接 cat, 不走 conda)|concatenate files (plain cat)."""
        self.logger.info(f"命令|Command: cat {' '.join(srcs)} > {dst}")
        with open(dst, "wb") as out:
            for s in srcs:
                with open(s, "rb") as f:
                    shutil.copyfileobj(f, out)

    # ----------------------------------------------------------
    # 步骤3: 三代转录本比对(Phase2)|Long-read align (GMAP -> genome GFF3)
    # ----------------------------------------------------------
    def _step_iso_align(self, genome: str) -> str:
        """GMAP 比对三代转录本 -> 基因组坐标 GFF3(供 TransDecoder 映射)
        |GMAP align long-read transcripts to genome-coordinate GFF3 (for TransDecoder)."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤3: 三代转录本比对|Step 3: Long-read alignment (GMAP)")
        self.logger.info("-" * 70)

        iso_gff = os.path.join(c.iso_align_dir, "iso.gff3")
        if os.path.exists(iso_gff):
            self.logger.info(f"跳过(已完成)|Skipping (done): {iso_gff}")
            return iso_gff

        genome_abs = os.path.abspath(genome)
        isoseq_abs = os.path.abspath(c.isoseq)

        # 3.1 GMAP 建库|build index
        cmd = build_conda_command(
            c.gmap_build_bin, ["-d", c.species, "-D", os.path.abspath(c.iso_align_dir), genome_abs]
        )
        if not self.cmd.run(cmd, "GMAP 建库|GMAP build index"):
            raise RuntimeError("GMAP 建库失败|GMAP build failed")

        # 3.2 GMAP 比对 -> gff3_gene(基因组坐标 mRNA/exon)|align to gff3_gene
        cmd = build_conda_command(
            c.gmap_bin,
            ["-d", c.species, "-D", os.path.abspath(c.iso_align_dir),
             "-t", str(c.threads), "-f", "gff3_gene", isoseq_abs],
        )
        # gmap 默认输出到 stdout, 重定向到文件|gmap prints to stdout, redirect
        with open(iso_gff, "w") as out_fh:
            ok = self.cmd.run(cmd, "GMAP 比对三代转录本|GMAP align iso-seq", stdout=out_fh)
        if not ok or not os.path.exists(iso_gff):
            raise RuntimeError("GMAP 比对失败|GMAP align failed")
        self.logger.info(f"三代比对完成|Long-read alignment done: {iso_gff}")
        return iso_gff

    # ----------------------------------------------------------
    # 步骤4: 蛋白 hints(Phase2)|Protein hints (miniprot)
    # ----------------------------------------------------------
    def _step_protein_hints(self, genome: str) -> str:
        """miniprot 比对蛋白 -> Augustus CDSpart hints|miniprot align to CDSpart hints."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤4: 蛋白证据 hints|Step 4: Protein hints (miniprot)")
        self.logger.info("-" * 70)

        protein_hints = os.path.join(c.hints_dir, "protein.gff")
        miniprot_gff = os.path.join(c.protein_align_dir, "miniprot.gff3")

        # 转换层: protein.gff 非空才跳过(0 字节 = 上次转换失败需重转)
        # |conversion: skip only if protein.gff non-empty (0-byte = prior failure)
        if os.path.exists(protein_hints) and os.path.getsize(protein_hints) > 0:
            self.logger.info(f"跳过(已完成)|Skipping (done): {protein_hints}")
            return protein_hints

        # miniprot 比对层: miniprot.gff3 存在则复用, 否则重跑|reuse gff3 if exists, else run
        if not os.path.exists(miniprot_gff):
            genome_abs = os.path.abspath(genome)
            prot_abs = os.path.abspath(c.prot_seq)
            cmd = build_conda_command(
                c.miniprot_bin, ["--gff", "-t", str(c.threads), genome_abs, prot_abs]
            )
            with open(miniprot_gff, "w") as out_fh:
                ok = self.cmd.run(cmd, "miniprot 蛋白比对|miniprot align proteins", stdout=out_fh)
            if not ok or not os.path.exists(miniprot_gff):
                raise RuntimeError("miniprot 比对失败|miniprot align failed")
        else:
            self.logger.info(f"复用已有 miniprot GFF|Reusing existing miniprot GFF: {miniprot_gff}")

        n = miniprot_gff_to_protein_hints(miniprot_gff, protein_hints)
        self.logger.info(f"蛋白 hints 生成|Protein hints: {protein_hints} ({n} CDSpart)")
        if n == 0:
            self.logger.warning("蛋白 hints 为空(miniprot 无 CDS 比对?)|Protein hints empty (no CDS hits?)")
        return protein_hints

    # ----------------------------------------------------------
    # 步骤5: intron hints + 合并 hintsfile|Intron hints + merged hintsfile
    # ----------------------------------------------------------
    def _step_intron_hints(self, bam: str) -> str:
        """bam2hints 抽 intron hints|Extract intron hints via bam2hints."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤5a: intron hints|Step 5a: Intron hints (bam2hints)")
        self.logger.info("-" * 70)

        intron_gff = os.path.join(c.hints_dir, "intron.gff")
        # 非空才视为完成(0 字节 = 上次失败需重跑)|non-empty = done (0-byte = prior failure)
        if os.path.exists(intron_gff) and os.path.getsize(intron_gff) > 0:
            self.logger.info(f"跳过(已完成)|Skipping (done): {intron_gff}")
            return intron_gff

        cmd = build_conda_command(c.bam2hints_bin, ["--in=" + bam, "--out=" + intron_gff])
        # bam2hints 链接 libbamtools.so.2.5.2, Augustus env 缺该库 -> 注入 LD_LIBRARY_PATH 兜底
        # |bam2hints links libbamtools.so.2.5.2 absent from env; inject LD_LIBRARY_PATH
        env = os.environ.copy()
        ld = self._libbamtools_ld_path()
        if ld:
            env["LD_LIBRARY_PATH"] = ld + ":" + env.get("LD_LIBRARY_PATH", "")
        if not self.cmd.run(cmd, "bam2hints 抽剪接位点|bam2hints introns", env=env):
            raise RuntimeError("bam2hints 失败|bam2hints failed")
        self.logger.info(f"intron hints 完成|Intron hints done: {intron_gff}")
        return intron_gff

    def _libbamtools_ld_path(self) -> str:
        """为 bam2hints 准备 libbamtools.so.2.5.2(Augustus env 缺库, 用 2.5.x 软链兜底)
        |Provide libbamtools.so.2.5.2 (env lacks it; symlink from a 2.5.x).

        bam2hints 二进制链接 BamTools soname 2.5.2, 但 Augustus_v.3.5.0 env 没装该库。
        BamTools 2.5.x 为 patch 版 ABI 兼容, 用现成的 2.5.3(braker env) 软链兜底。
        |bam2hints links BamTools soname 2.5.2 but the env lacks it; 2.5.x patch
        versions are ABI-compatible, so symlink from an available 2.5.x (braker env).
        """
        cache = Path(os.path.expanduser("~/.cache/biopytools/lib"))
        cache.mkdir(parents=True, exist_ok=True)
        link = cache / "libbamtools.so.2.5.2"
        if link.is_symlink() or link.exists():
            return str(cache)
        # 找现成 libbamtools.so.2.*(排除 2.5.2 自身)|find existing 2.5.x (exclude self)
        cands = [
            x for x in glob.glob(os.path.expanduser("~/miniforge3/envs/*/lib/libbamtools.so.2.*"))
            if not x.endswith("libbamtools.so.2.5.2")
        ]
        if not cands:
            self.logger.warning(
                "找不到 libbamtools.so.2.*, bam2hints 可能失败|No libbamtools.so.2.* found"
            )
            return ""
        src = cands[0]
        try:
            os.symlink(src, str(link))
            self.logger.info(f"创建 libbamtools 软链兜底|Created libbamtools symlink: {link} -> {src}")
        except OSError as e:
            self.logger.warning(f"创建软链失败|Symlink failed: {e}")
            return ""
        return str(cache)

    def _step_merge_hints(self, intron_gff: str, protein_hints: str) -> str:
        """合并 intron + protein hints -> Augustus hintsfile|Merge into hintsfile."""
        c = self.config
        parts = [p for p in (intron_gff, protein_hints) if p and os.path.exists(p)]
        if not parts:
            self.logger.info("无 hints 证据, Augustus 将 ab initio 预测|No hints: ab initio Augustus")
            return ""
        hintsfile = os.path.join(c.hints_dir, "hintsfile.gff")
        self.logger.info(f"合并 hints|Merging hints: {[os.path.basename(p) for p in parts]} -> {hintsfile}")
        with open(hintsfile, "w") as out:
            for p in parts:
                with open(p) as f:
                    shutil.copyfileobj(f, out)
        return hintsfile

    # ----------------------------------------------------------
    # 步骤6: 训练集(三代->TransDecoder / 否则 GeneMark)|Training set
    # ----------------------------------------------------------
    def _step_training_set(self, genome: str, intron_gff: str, iso_gff: str) -> str:
        """训练集来源分流:三代->TransDecoder; 否则 GeneMark(ET/EP+/ES)
        |Branch: long-read -> TransDecoder; else GeneMark (ET/EP+/ES)."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤6: 训练集生成|Step 6: Training set")
        self.logger.info("-" * 70)

        # 有三代 -> T2T 法: GMAP+TransDecoder 基因组坐标 gene models|long-read -> T2T method
        if iso_gff:
            try:
                gtf = self._step_training_from_isoseq(genome, iso_gff)
                if gtf:
                    return gtf
            except Exception as e:
                self.logger.warning(
                    f"三代训练集失败, 回退 GeneMark|Iso-seq training failed, fallback GeneMark: {e}"
                )
        return self._step_genemark(genome, intron_gff)

    def _step_training_from_isoseq(self, genome: str, iso_gff: str) -> str:
        """三代 GMAP GFF3 + TransDecoder -> 基因组坐标 gene models
        |Long-read GMAP GFF3 + TransDecoder -> genome-coordinate gene models."""
        c = self.config
        genome_gtf = os.path.join(c.training_dir, "iso_training.gff3")
        if os.path.exists(genome_gtf):
            self.logger.info(f"跳过(已完成)|Skipping (done): {genome_gtf}")
            return genome_gtf

        isoseq_abs = os.path.abspath(c.isoseq)
        td_gff = os.path.join(c.training_dir, os.path.basename(isoseq_abs) + ".transdecoder.gff3")

        # 6a TransDecoder.LongOrfs + Predict(在转录本序列上)|on transcript seqs, cwd=training_dir
        cmd = build_conda_command(c.transdecoder_longorfs_bin, ["-t", isoseq_abs])
        if not self.cmd.run(cmd, "TransDecoder.LongOrfs", cwd=c.training_dir):
            raise RuntimeError("TransDecoder.LongOrfs 失败|TransDecoder.LongOrfs failed")
        cmd = build_conda_command(c.transdecoder_predict_bin, ["-t", isoseq_abs])
        if not self.cmd.run(cmd, "TransDecoder.Predict", cwd=c.training_dir):
            raise RuntimeError("TransDecoder.Predict 失败|TransDecoder.Predict failed")

        if not os.path.exists(td_gff):
            raise RuntimeError(f"TransDecoder 未产出 {td_gff}|TransDecoder output missing")

        # 6b cdna_alignment_orf_to_genome_orf.pl: ORF(转录本坐标)->基因组坐标
        # |map ORF (transcript coords) to genome via GMAP alignment
        # 参数顺序: <orf.gff3> <alignment.gff3> <transcript.fa>|arg order per TransDecoder util
        cmd = build_conda_command(
            c.cdna_orf_to_genome_bin, [td_gff, os.path.abspath(iso_gff), isoseq_abs]
        )
        with open(genome_gtf, "w") as out_fh:
            ok = self.cmd.run(cmd, "ORF 映射到基因组|Map ORFs to genome", cwd=c.training_dir, stdout=out_fh)
        if not ok or not os.path.exists(genome_gtf) or os.path.getsize(genome_gtf) == 0:
            raise RuntimeError("ORF 基因组映射失败或空|ORF-to-genome mapping failed or empty")
        self.logger.info(f"三代训练集完成|Iso-seq training set done: {genome_gtf}")
        return genome_gtf

    def _step_genemark(self, genome: str, intron_gff: str) -> str:
        """GeneMark 自训练(ET 有 intron / EP+ 有蛋白 / ES 纯基因组)|GeneMark auto-mode."""
        c = self.config
        genemark_gtf = os.path.join(c.training_dir, "genemark.gtf")
        if os.path.exists(genemark_gtf) and os.path.getsize(genemark_gtf) > 0:
            self.logger.info(f"跳过(已完成)|Skipping (done): {genemark_gtf}")
            return genemark_gtf

        # GeneMark probuild(C) + gmes_petap.pl 的 abs_path 不支持中文路径
        # -> ~/tmp 下 ASCII 临时目录: 真实拷贝输入进去, cwd 运行, 完成后拷回 genemark.gtf
        # |probuild + abs_path mangle CJK -> ASCII temp dir under ~/tmp: copy inputs in,
        # run with cwd=here, copy genemark.gtf back
        work = self._ascii_workdir("gm")
        seq_in = self._copy_in(os.path.abspath(genome), work / "genome.fa")

        args = ["--sequence", str(seq_in), "--cores", str(c.threads)]
        # 始终用 GeneMark-ES(基因组自训练): 最稳健。
        # GeneMark-ET 需带 read-support 分数 + 链信息的 intron; 非链特异性 RNA-seq 的
        # bam2hints 输出(score=0/strand=.)会致 parse_ET.pl 除零崩溃(已实测)。
        # intron + 蛋白证据全部改走 Augustus hints(见 _step_merge_hints/_step_augustus),
        # 不依赖 GeneMark-ET/EP+, 证据不浪费。
        # |Always GeneMark-ES (genome self-training): most robust. GeneMark-ET needs introns
        # with score+strand; unstranded bam2hints output crashes parse_ET.pl (verified).
        # Intron + protein evidence flow to Augustus hints instead, so nothing is wasted.
        self.logger.info("使用 GeneMark-ES(基因组自训练)|Using GeneMark-ES (genome self-training)")
        args.append("--ES")
        # intron_gff 不再喂给 GeneMark(改由 Augustus hints 使用)|introns go to Augustus, not GeneMark
        _ = intron_gff

        cmd = build_genemark_command(c.gmes_petap_path, c.genemark_perl_env, args)
        if not self.cmd.run(cmd, "运行 GeneMark|Run GeneMark", cwd=str(work)):
            raise RuntimeError("GeneMark 失败|GeneMark failed")

        gtf_src = work / "genemark.gtf"
        if not gtf_src.exists():
            raise RuntimeError(f"GeneMark 未产出 genemark.gtf|genemark.gtf not found: {gtf_src}")
        shutil.copy2(gtf_src, genemark_gtf)
        self.logger.info(f"GeneMark 训练集完成(临时目录|tmp dir: {work})|GeneMark training set done: {genemark_gtf}")
        return genemark_gtf

    # ----------------------------------------------------------
    # 步骤7: Augustus 训练 + 预测|Augustus train + predict
    # ----------------------------------------------------------
    def _step_augustus(self, genome: str, training_gtf: str, hintsfile: str) -> str:
        """etraining + augustus 预测|Augustus etraining + prediction."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤7: Augustus 训练 + 预测|Step 7: Augustus train + predict")
        self.logger.info("-" * 70)

        # Augustus(C)+etraining: ~/tmp ASCII 临时目录 + 直接调二进制。
        # conda run 的 activation 会把 AUGUSTUS_CONFIG_PATH 覆盖成 env config, 故必须直接调
        # (不走 build_conda_command)。转换脚本(Perl)可直接调、能读中文路径。
        # |run via ~/tmp ASCII dir + direct binary call (conda run activation overrides
        # AUGUSTUS_CONFIG_PATH to env config, so direct call is required).
        work = self._ascii_workdir("aug")
        augustus_config = self._prepare_augustus_config()
        aug_env = os.environ.copy()
        aug_env["AUGUSTUS_CONFIG_PATH"] = augustus_config

        species_dir = os.path.join(augustus_config, "species", c.species)
        params_exist = (
            os.path.isdir(species_dir)
            and any(f.endswith("_exon_probs.pbl") for f in os.listdir(species_dir))
        ) if os.path.isdir(species_dir) else False

        if params_exist:
            self.logger.info(f"跳过 etraining(species 参数已存在)|Skip etraining (params exist)")
        else:
            # etraining 要 GenBank 位置参数(不是 --trainfile)。GTF->GFF->GenBank 两步转换:
            # gtf2gff.pl 读 GTF(Perl 处理中文路径, via stdin) -> gff2gbSmallDNA.pl(+基因组)
            # ->GenBank。输出到 ASCII work 目录供 C 二进制读取。
            # |etraining takes GenBank positional (not --trainfile). GTF->GFF->GenBank;
            # Perl converters handle CJK input paths; outputs land in ASCII work dir.
            training_gb = work / "training.gb"
            if not training_gb.exists():
                gff_tmp = work / "genemark.gff"
                with open(training_gtf, "rb") as gtf_fh:
                    if not self.cmd.run([c.gtf2gff_bin, "--out=" + str(gff_tmp)],
                                        "gtf2gff GTF->GFF", stdin=gtf_fh):
                        raise RuntimeError("gtf2gff 失败|gtf2gff failed")
                if not self.cmd.run(
                    [c.gff2gb_bin, str(gff_tmp), os.path.abspath(genome), "1000", str(training_gb)],
                    "gff2gbSmallDNA -> GenBank",
                ):
                    raise RuntimeError("gff2gbSmallDNA 失败|gff2gbSmallDNA failed")
            # etraining: 直接调, GenBank 作位置参数|direct call, GenBank positional
            if not self.cmd.run([c.etraining_bin, "--species=" + c.species, str(training_gb)],
                                "Augustus etraining", env=aug_env):
                raise RuntimeError("etraining 失败|etraining failed")
            self.logger.info(f"etraining 完成|etraining done: {species_dir}")

        out_gff = os.path.join(c.augustus_dir, f"{c.species}.augustus.gff")
        if os.path.exists(out_gff) and os.path.getsize(out_gff) > 0:
            self.logger.info(f"跳过预测(已完成)|Skip predict (done): {out_gff}")
            return out_gff

        # predict: 直接调, --gff3=on(非 --gff3), 基因组+hintsfile 用 ASCII 拷贝
        # |predict: direct call, --gff3=on (not --gff3), ASCII copies for genome+hints
        pred_args = ["--species=" + c.species, "--gff3=on"]
        if hintsfile and os.path.exists(hintsfile):
            hints_in = self._copy_in(os.path.abspath(hintsfile), work / "hintsfile.gff")
            pred_args.extend(["--hintsfile=" + str(hints_in),
                              "--extrinsicCfgFile=" + os.path.abspath(c.extrinsic_cfg)])
        genome_in = self._copy_in(os.path.abspath(genome), work / "genome.fa")
        pred_args.append(str(genome_in))

        with open(out_gff, "w") as out_fh:
            ok = self.cmd.run([c.augustus_bin] + pred_args, "Augustus 预测|Augustus prediction",
                              env=aug_env, stdout=out_fh)
        if not ok or not os.path.exists(out_gff) or os.path.getsize(out_gff) == 0:
            raise RuntimeError("Augustus 预测失败或空输出|Augustus prediction failed or empty")
        # 拷回 species 参数(可见性 + 自包含)|copy species params back to output dir
        self._copyback_species_params(augustus_config, c.augustus_dir, c.species)
        self.logger.info(f"预测完成(临时目录|tmp dir: {work})|Prediction done: {out_gff}")
        return out_gff

    def _copyback_species_params(self, aug_config: str, target_dir, species: str):
        """把 ASCII 临时 config 里的 species 参数拷回输出目录(失败不影响结果)
        |Copy species params from ASCII temp config back to output (non-fatal)."""
        src_species = os.path.join(aug_config, "species", species)
        if not os.path.isdir(src_species):
            return
        dst_species = os.path.join(str(target_dir), "augustus_config", "species", species)
        try:
            if os.path.exists(dst_species):
                shutil.rmtree(dst_species)
            os.makedirs(os.path.dirname(dst_species), exist_ok=True)
            shutil.copytree(src_species, dst_species)
            self.logger.info(f"species 参数已拷回|Species params copied back: {dst_species}")
        except Exception as e:
            self.logger.warning(f"species 参数拷回失败(不影响结果)|Copyback failed (non-fatal): {e}")

    def _prepare_augustus_config(self) -> str:
        """拷贝 Augustus config 到 ASCII 临时目录 + 从 generic seed species
        |Copy config to ASCII temp dir + seed species from generic.

        etraining 不自动创建新 species, 需 seed: 把 generic 参数文件按原名拷进
        species/<species>/(内部引用 generic_weightmatrix.txt 等需原名), 并建
        <species>_parameters.cfg(从 generic_parameters.cfg, etraining 按此名查找)。
        |etraining does not auto-create a new species; seed by copying generic param
        files (keeping names so internal refs resolve) + <species>_parameters.cfg.
        """
        c = self.config
        dst = self._ascii_workdir("aug") / "augustus_config"
        species_dst = dst / "species" / c.species
        # 已拷贝 + seed(species 目录存在)则跳过|skip if copied + seeded (species dir exists)
        if species_dst.is_dir():
            return str(dst)
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(c.augustus_config_src, str(dst))
        generic_dir = dst / "species" / "generic"
        species_dst.mkdir(parents=True, exist_ok=True)
        for f in generic_dir.glob("generic_*"):
            shutil.copy2(str(f), str(species_dst / f.name))
        shutil.copy2(
            str(generic_dir / "generic_parameters.cfg"),
            str(species_dst / f"{c.species}_parameters.cfg"),
        )
        self.logger.info(f"已 seed species {c.species}(from generic)|Seeded species {c.species}")
        return str(dst)

    # ----------------------------------------------------------
    # 步骤8: LTR 注释(Phase2, 正交 TE)|LTR annotation (orthogonal)
    # ----------------------------------------------------------
    def _step_ltr(self, genome: str):
        """gt ltrharvest + LTR_retriever(失败由上层捕获)|gt ltrharvest + LTR_retriever."""
        c = self.config
        self.logger.info("-" * 70)
        self.logger.info("步骤8: LTR 注释|Step 8: LTR annotation")
        self.logger.info("-" * 70)

        genome_base = os.path.basename(c.genome)
        ltr_marker = os.path.join(c.ltr_dir, genome_base + ".pass.list")
        if os.path.exists(ltr_marker) and os.path.getsize(ltr_marker) > 0:
            self.logger.info(f"跳过(已完成)|Skipping (done): {ltr_marker}")
            return

        # gt(C)不支持中文路径 -> ~/tmp 下 ASCII 临时目录运行, 完成后拷回关键产物
        # |gt (C) can't do CJK -> run in ASCII temp dir, copy key outputs back
        work = self._ascii_workdir("ltr")
        genome_in = self._copy_in(os.path.abspath(genome), work / "genome.fa")
        genome_name = "genome.fa"
        index_name = str(work / (genome_name + ".suf"))
        harvest_out = str(work / (genome_name + ".harvest.scn"))

        # 8a suffixerator 建 suffix array 索引(-lcp 生成 lcp 表, ltrharvest 需要; 旧版用 -lis)
        # |build suffix array index (-lcp makes lcp table required by ltrharvest; old flag was -lis)
        cmd = build_conda_command(
            c.gt_bin,
            ["suffixerator", "-db", str(genome_in), "-indexname", index_name,
             "-dna", "-suf", "-lcp", "-des", "-ssp", "-md5"],
        )
        if not self.cmd.run(cmd, "gt suffixerator 建索引|gt suffixerator", cwd=str(work)):
            raise RuntimeError("gt suffixerator 失败|gt suffixerator failed")

        # 8b ltrharvest 候选 LTR|harvest candidate LTRs
        cmd = build_conda_command(
            c.gt_bin,
            ["ltrharvest", "-index", index_name, "-seqids", "yes",
             "-minlenltr", "100", "-maxlenltr", "7000", "-mintsd", "4", "-maxtsd", "6",
             "-similar", "85", "-vic", "10", "-seed", "20"],
        )
        with open(harvest_out, "w") as out_fh:
            ok = self.cmd.run(cmd, "gt ltrharvest", cwd=str(work), stdout=out_fh)
        if not ok:
            raise RuntimeError("gt ltrharvest 失败|gt ltrharvest failed")

        # 8c LTR_retriever 过滤+注释|filter + annotate
        cmd = build_conda_command(
            c.ltr_retriever_bin,
            ["-genome", str(genome_in), "-inharvest", harvest_out, "-threads", str(c.threads)],
        )
        if not self.cmd.run(cmd, "LTR_retriever|LTR_retriever", cwd=str(work)):
            raise RuntimeError("LTR_retriever 失败|LTR_retriever failed")

        # 拷回 LTR_retriever 关键产物|copy key outputs back
        n = 0
        for f in work.glob(genome_name + "*"):
            try:
                shutil.copy2(str(f), os.path.join(c.ltr_dir, os.path.basename(f)))
                n += 1
            except Exception as e:
                self.logger.debug(f"拷回失败|copyback failed: {f}: {e}")
        self.logger.info(f"LTR 注释完成(拷回 {n} 个文件)|LTR done ({n} files back): {c.ltr_dir}")

    # ----------------------------------------------------------
    # 版本记录|Version record
    # ----------------------------------------------------------
    def _write_versions(self, final_gff: str):
        """写 software_versions.yml(失败不阻断)|Write versions (non-blocking)."""
        c = self.config
        out = c.pipeline_info_dir / "software_versions.yml"
        tools = {
            "RepeatModeler": c.repeatmodeler_bin, "RepeatMasker": c.repeatmasker_bin,
            "HISAT2": c.hisat2_bin, "samtools": c.samtools_bin,
            "GMAP": c.gmap_bin, "miniprot": c.miniprot_bin,
            "bam2hints": c.bam2hints_bin, "augustus": c.augustus_bin,
            "etraining": c.etraining_bin, "GeneMark-ES": c.gmes_petap_path,
        }
        if c.gt_bin and os.path.exists("/" + c.gt_bin.lstrip("/")):
            tools["gt"] = c.gt_bin
        if c.ltr_retriever_bin:
            tools["LTR_retriever"] = c.ltr_retriever_bin
        lines = ["pipeline:", "  name: biopytools oomycete_anno", "  version: '1.0.0'", "tools:"]
        for name, path in tools.items():
            ok, ver, _ = self.cmd.run_capture(build_conda_command(path, ["--version"]), timeout=30)
            ver = (ver or "").strip().splitlines()[0] if ok and (ver or "").strip() else "unknown"
            lines.append(f"  {name}:\n    version: '{ver}'\n    path: '{path}'")
        lines.append("parameters:")
        lines.append(f"  species: {c.species}")
        lines.append(f"  threads: {c.threads}")
        lines.append(f"  rna_strandness: '{c.rna_strandness}'")
        lines.append(f"  evidence: RNA-seq={bool(c.rnaseq_dirs)} iso-seq={bool(c.isoseq)} protein={bool(c.prot_seq)}")
        lines.append(f"  final_gff: {final_gff}")
        try:
            out.write_text("\n".join(lines) + "\n", encoding="utf-8")
            self.logger.info(f"版本信息已保存|Versions saved: {out}")
        except Exception as e:
            self.logger.warning(f"版本信息写入失败|Failed to write versions: {e}")
