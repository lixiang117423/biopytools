"""mixrace GTX 后端编排|mixrace GTX backend orchestration.

调 biopytools fastq2vcf-gtx(clean fastq → 比对+gVCF+联合calling 一步)产出
joint VCF 与 per-sample BAM;从 BAM 提取 mapped reads 供 k-mer 分析。
|Call fastq2vcf-gtx (clean fastq -> mapping + joint calling in one module);
extract mapped reads from BAMs for the k-mer stage.
"""
import os
import shlex
from pathlib import Path
from typing import Optional, Tuple

from .pipeline import _done
from .utils import get_conda_env

_EXCL = "0x904"          # unmapped+secondary+supplementary
_EXCL_INCL_UNMAPPED = "0x900"


def run_gtx(config, runner, ckpt, fastq_dir: str) -> Tuple[Optional[str], Optional[str]]:
    """step2: fastq2vcf-gtx 一步比对+联合calling|one-step mapping + joint calling.

    Returns:
        (bam_dir, joint_vcf);失败 (None, None) 不建断点。|(None, None) on failure.
    """
    out = Path(config.output_dir) / "03_gtx"
    joint_vcf = out / "04_joint_calling" / "gtx_joint_raw.vcf.gz"
    bam_dir = out / "03_mapping" / "bam"
    if config.enable_checkpoint and _done(ckpt, "gtx", joint_vcf):
        runner.logger.info("跳过已完成步骤|Skipping completed step: gtx")
        return str(bam_dir), str(joint_vcf)
    runner.logger.info("开始步骤|Starting step: GTX mapping+joint calling")
    ok, _, _ = runner.run(
        f"biopytools fastq2vcf-gtx --clean-fastq-dir {fastq_dir} -g {config.genome} "
        f"-o {out} -t {config.threads}",
        "GTX 比对+联合calling|GTX mapping + joint calling")
    if not ok:
        runner.logger.error("GTX 失败,未建断点|GTX failed, no checkpoint")
        return None, None
    if not runner.dry_run and not joint_vcf.exists():
        runner.logger.error(f"GTX 未产出 joint VCF: {joint_vcf}|GTX produced no joint VCF")
        return None, None
    if config.enable_checkpoint:
        ckpt.create("gtx")
    return str(bam_dir), str(joint_vcf)


def extract_mapped_fastq(config, runner, ckpt, sample: str,
                        bam: str) -> Optional[Tuple[str, str]]:
    """step4 前置:从 BAM 提取达标 mapped reads 成对 fastq|mapped-read pair extraction.

    samtools view -F 0x904 [-q Q] 过滤后 samtools fastq;输出
    05_kmer/mapped_fastq/{sample}_1/2.mapped.fq.gz。
    """
    out_dir = Path(config.output_dir) / "05_kmer" / "mapped_fastq"
    out_dir.mkdir(parents=True, exist_ok=True)
    r1 = out_dir / f"{sample}_1.mapped.fq.gz"
    r2 = out_dir / f"{sample}_2.mapped.fq.gz"
    if config.enable_checkpoint and _done(ckpt, f"mapped_{sample}", r1):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: mapped_{sample}")
        return str(r1), str(r2)
    q = config.min_mapq
    qpart = f" -q {q}" if q > 0 else ""
    env = get_conda_env(config.samtools_path)
    st_path = shlex.quote(config.samtools_path)
    runner.logger.info(f"开始步骤|Starting step: extract mapped reads {sample}")
    ok, _, _ = runner.run(
        f"conda run -n {env} --no-capture-output bash -c "
        f"'{st_path} view -b -F {_EXCL}{qpart} {bam} - | "
        f"{st_path} fastq -1 {r1} -2 {r2} -'",
        f"提取mapped reads|extract mapped reads {sample} (MAPQ>={q if q > 0 else 0})")
    if not ok:
        runner.logger.error(f"mapped reads 提取失败: {sample}|extraction failed")
        return None
    if config.enable_checkpoint and (runner.dry_run or r1.exists()):
        ckpt.create(f"mapped_{sample}")
    return str(r1), str(r2)


def count_mapped(runner, config, bam: str) -> Tuple[Optional[int], Optional[int]]:
    """(总primary含unmapped, 达标mapped) 计数|(total primary, qualified mapped)."""
    st = config.samtools_path
    ok_t, out_t, _ = runner.run_conda(st, ["view", "-c", "-F", _EXCL_INCL_UNMAPPED, bam],
                                      f"计数total|count total {os.path.basename(bam)}")
    args_m = ["view", "-c", "-F", _EXCL]
    if config.min_mapq > 0:
        args_m += ["-q", str(config.min_mapq)]
    args_m.append(bam)
    ok_m, out_m, _ = runner.run_conda(st, args_m,
                                      f"计数mapped|count mapped {os.path.basename(bam)}")

    def _parse(ok, out):
        return int(out.strip()) if ok and out.strip().isdigit() else None

    return _parse(ok_t, out_t), _parse(ok_m, out_m)
