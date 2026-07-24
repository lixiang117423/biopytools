"""rnaseq2vcf 数据处理步骤|pipeline steps (index, qc, align, call, filter)"""

import os
import shlex
from typing import List, Tuple

from .utils import build_conda_command, get_conda_env, CommandRunner


class GenomeIndexer:
    """构建共享基因组索引(faidx/dict/hisat2 索引/剪接位点)|Build shared genome index"""

    def __init__(self, config, logger, runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.runner = runner

    def _hisat2_sibling(self, name: str) -> str:
        """从 hisat2_path 推导同 env 下工具路径|Derive sibling tool path in same env"""
        return self.config.hisat2_path.replace('/hisat2', f'/{name}')

    def _index_commands(self) -> List[List[str]]:
        """faidx + dict + hisat2-build(始终)|always"""
        cfg = self.config
        ht2_prefix = os.path.join(cfg.genome_index_dir, cfg.genome_name)
        return [
            build_conda_command(cfg.samtools_path, ['faidx', cfg.ref_genome_fa]),
            build_conda_command(cfg.gatk_path, ['CreateSequenceDictionary', '-R', cfg.ref_genome_fa]),
            build_conda_command(self._hisat2_sibling('hisat2-build'), [cfg.ref_genome_fa, ht2_prefix]),
        ]

    def _splice_sites_command(self) -> List[str]:
        """hisat2_extract_splice_sites(仅当提供 gff3)|only when gff3 provided"""
        cfg = self.config
        return build_conda_command(self._hisat2_sibling('hisat2_extract_splice_sites.py'), [cfg.gff3_file])

    def run(self) -> bool:
        cfg = self.config
        idx_dir = cfg.genome_index_dir
        ht2_prefix = os.path.join(idx_dir, cfg.genome_name)
        ss_file = os.path.join(idx_dir, f"{cfg.genome_name}.ss")
        has_ss = bool(cfg.gff3_file)
        done = (os.path.exists(cfg.ref_genome_fa + '.fai')
                and os.path.exists(ht2_prefix + '.1.ht2')
                and (os.path.exists(ss_file) if has_ss else True))
        if done and not cfg.force:
            self.logger.info("跳过已完成步骤|Skipping completed step: genome_index")
            return True
        self.logger.info("开始构建基因组索引|Building genome index")
        cmds = self._index_commands()
        descs = ["samtools faidx", "gatk CreateSequenceDictionary", "hisat2-build"]
        for i in range(3):
            if not self.runner.run(cmds[i], f"索引|index: {descs[i]}"):
                return False
        if has_ss:
            ss_cmd = ' '.join(self._splice_sites_command()) + f" > '{ss_file}'"
            if not self.runner.run(ss_cmd, "提取剪接位点|extract splice sites"):
                return False
        else:
            self.logger.info("未提供 GFF3,跳过剪接位点(HISAT2 将 de novo 发现 junction)|"
                             "No GFF3, skipping splice-site extraction (HISAT2 de novo)")
        return True


class QualityController:
    """fastp 质控(自递归 biopytools fastp)|QC via biopytools fastp"""

    def __init__(self, config, logger, runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.runner = runner

    def run_sample(self, sample: str, r1: str, r2: str) -> Tuple[str, str]:
        cfg = self.config
        qc_dir = os.path.join(cfg.output_dir, sample, "01_qc")
        os.makedirs(qc_dir, exist_ok=True)
        r1_clean = os.path.join(qc_dir, f"{sample}_1.clean.fq.gz")
        r2_clean = os.path.join(qc_dir, f"{sample}_2.clean.fq.gz")
        if os.path.exists(r1_clean) and os.path.exists(r2_clean) and not cfg.force:
            self.logger.info(f"跳过已完成 QC|Skipping QC: {sample}")
            return r1_clean, r2_clean
        cmd = ['biopytools', 'fastp', '-i', r1, '-I', r2,
               '-o', r1_clean, '-O', r2_clean, '-w', str(cfg.threads)]
        if not self.runner.run(cmd, f"fastp 质控|fastp QC: {sample}"):
            raise RuntimeError(f"QC 失败|QC failed: {sample}")
        return r1_clean, r2_clean


class Aligner:
    """HISAT2 剪接感知比对 → 排序 BAM|HISAT2 splice-aware align → sorted BAM"""

    def __init__(self, config, logger, runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.runner = runner

    def _align_command_str(self, sample: str, r1: str, r2: str, out_bam: str, log_file: str) -> str:
        """构建比对命令(单 conda run bash -c 包裹管道)|
        Build align command (single conda run bash -c wrapping the pipe).

        hisat2 与 samtools 同在 RNA_Seq env,用单个 conda run 包裹管道(避免 conda run|conda run,§13.2.1)。
        hisat2 stderr 重定向到日志文件;samtools 写到 -o 文件(stdout 为空,无二进制污染)。
        """
        cfg = self.config
        ht2_prefix = os.path.join(cfg.genome_index_dir, cfg.genome_name)
        t = str(cfg.threads)
        env = get_conda_env(cfg.hisat2_path)  # RNA_Seq
        ss_opt = ""
        if cfg.gff3_file:
            ss_file = os.path.join(cfg.genome_index_dir, f"{cfg.genome_name}.ss")
            ss_opt = f"--known-splicesite-infile {ss_file} "
        # env 内 hisat2/samtools 在 PATH,用裸名调用|bare names resolve inside conda env
        pipeline = (f"hisat2 -x {ht2_prefix} {ss_opt}-p {t} --dta -1 {r1} -2 {r2} 2> {log_file} | "
                    f"samtools sort -@ {t} -o {out_bam} -")
        return f"conda run -n {env} --no-capture-output bash -c '{pipeline}'"

    def run_sample(self, sample: str, r1: str, r2: str) -> str:
        cfg = self.config
        align_dir = os.path.join(cfg.output_dir, sample, "02_align")
        os.makedirs(align_dir, exist_ok=True)
        out_bam = os.path.join(align_dir, f"{sample}.sorted.bam")
        log_file = os.path.join(align_dir, f"{sample}.hisat2.log")
        if os.path.exists(out_bam) and not cfg.force:
            self.logger.info(f"跳过已完成比对|Skipping alignment: {sample}")
            return out_bam
        cmd_str = self._align_command_str(sample, r1, r2, out_bam, log_file)
        desc = f"HISAT2 比对|HISAT2 align: {sample} (日志|log: {log_file})"
        if not self.runner.run_with_progress(cmd_str, desc):
            raise RuntimeError(f"比对失败|Alignment failed: {sample}")
        return out_bam


class Caller:
    """GATK RNA-seq calling: RG→MarkDup→SplitNCigar→HaplotypeCaller|GATK RNA-seq calling"""

    def __init__(self, config, logger, runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.runner = runner

    def _gatk(self, args: List[str]) -> List[str]:
        return build_conda_command(self.config.gatk_path, args)

    def _add_rg_command(self, sample, in_bam, out_bam):
        return self._gatk(['AddOrReplaceReadGroups', '-I', in_bam, '-O', out_bam,
                           '--RGID', sample, '--RGLB', 'lib1', '--RGPL', 'illumina',
                           '--RGPU', 'unit1', '--RGSM', sample])

    def _markdup_command(self, in_bam, out_bam, metrics):
        return self._gatk(['MarkDuplicates', '-I', in_bam, '-O', out_bam,
                           '-M', metrics, '--CREATE_INDEX'])

    def _splitncigar_command(self, sample, in_bam, out_bam):
        cfg = self.config
        return self._gatk(['SplitNCigarReads', '-R', cfg.ref_genome_fa, '-I', in_bam, '-O', out_bam])

    def _haplotype_command(self, sample, in_bam, out_vcf):
        cfg = self.config
        return self._gatk(['HaplotypeCaller', '-R', cfg.ref_genome_fa, '-I', in_bam,
                           '--dont-use-soft-clipped-bases',
                           '--standard-min-confidence-threshold-for-calling', str(cfg.min_conf),
                           '-O', out_vcf])

    def run_sample(self, sample: str, in_bam: str) -> str:
        cfg = self.config
        call_dir = os.path.join(cfg.output_dir, sample, "03_calling")
        os.makedirs(call_dir, exist_ok=True)
        rg = os.path.join(call_dir, f"{sample}.rg.bam")
        dedup = os.path.join(call_dir, f"{sample}.dedup.bam")
        split = os.path.join(call_dir, f"{sample}.split.bam")
        raw_vcf = os.path.join(call_dir, f"{sample}.raw.vcf.gz")
        if os.path.exists(raw_vcf) and not cfg.force:
            self.logger.info(f"跳过已完成 calling|Skipping calling: {sample}")
            return raw_vcf
        steps = [
            (self._add_rg_command(sample, in_bam, rg), f"加读组|AddOrReplaceReadGroups: {sample}"),
            (self._markdup_command(rg, dedup, os.path.join(call_dir, f"{sample}.metrics")),
             f"去重|MarkDuplicates: {sample}"),
            (self._splitncigar_command(sample, dedup, split), f"拆分 N reads|SplitNCigarReads: {sample}"),
            (self._haplotype_command(sample, split, raw_vcf), f"变异检测|HaplotypeCaller: {sample}"),
        ]
        for cmd, desc in steps:
            # shlex.join 正确引用含 shell 特殊字符(如 FS>30.0 的 >)的参数|properly quote shell-special args
            if not self.runner.run_with_progress(shlex.join(cmd), desc):
                raise RuntimeError(f"calling 失败|calling failed: {sample} ({desc})")
        return raw_vcf


class VariantFilter:
    """GATK VariantFiltration(FS/QD/cluster)+ bcftools PASS|GATK filter + PASS"""

    def __init__(self, config, logger, runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.runner = runner

    def _gatk_wrap(self, args):
        return build_conda_command(self.config.gatk_path, args)

    def _bcftools_wrap(self, args):
        return build_conda_command(self.config.bcftools_path, args)

    def _variant_filtration_command(self, in_vcf, out_vcf):
        cfg = self.config
        return self._gatk_wrap(['VariantFiltration', '-R', cfg.ref_genome_fa, '-V', in_vcf,
                                '--window', str(cfg.cluster_window), '--cluster', str(cfg.cluster_size),
                                '--filter-name', 'FS', '--filter-expression', f"FS > {cfg.fs_threshold}",
                                '--filter-name', 'QD', '--filter-expression', f"QD < {cfg.qd_threshold}",
                                '-O', out_vcf])

    def run_sample(self, sample: str, raw_vcf: str) -> str:
        cfg = self.config
        f_dir = os.path.join(cfg.output_dir, sample, "04_filter")
        os.makedirs(f_dir, exist_ok=True)
        filtered = os.path.join(f_dir, f"{sample}.filtered.vcf.gz")
        pass_vcf = os.path.join(f_dir, f"{sample}.pass.vcf.gz")
        if os.path.exists(pass_vcf) and not cfg.force:
            self.logger.info(f"跳过已完成过滤|Skipping filtering: {sample}")
            return pass_vcf
        if not self.runner.run_with_progress(shlex.join(self._variant_filtration_command(raw_vcf, filtered)),
                                             f"GATK 过滤|VariantFiltration: {sample}"):
            raise RuntimeError(f"过滤失败|filter failed: {sample}")
        # bcftools 从 VariantFiltration 产物(filtered)提取 PASS|extract PASS from filtered
        view_cmd = self._bcftools_wrap(['view', '-f', 'PASS', '-Oz', '-o', pass_vcf, filtered])
        if not self.runner.run(shlex.join(view_cmd), f"提取 PASS|bcftools PASS: {sample}"):
            raise RuntimeError(f"PASS 提取失败|PASS extract failed: {sample}")
        if not self.runner.run(shlex.join(self._bcftools_wrap(['index', pass_vcf])), f"索引|index: {pass_vcf}"):
            self.logger.warning(f"VCF 索引失败(可继续)|VCF index failed (continuing): {pass_vcf}")
        return pass_vcf
