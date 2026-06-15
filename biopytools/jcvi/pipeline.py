"""
JCVI共享分析管道|JCVI Shared Analysis Pipeline

mcscan和allelic共享的步骤1-6:
1. 发现样本
2. 提取蛋白质/CDS序列 (gffread)
3. GFF→BED + uniq
4. LAST比对
5. Blastfilter过滤
6. Synteny scan → .anchors + .lifted.anchors
"""

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Tuple
from itertools import combinations

from .config import JcviBaseConfig
from .utils import (
    JcviLogger,
    build_jcvi_command,
    discover_samples,
    get_sample_name,
    get_jcvi_stem,
)


@dataclass
class PipelineResult:
    """管道运行结果|Pipeline run result"""
    samples: List[str] = field(default_factory=list)
    pair_list: List[Tuple[str, str]] = field(default_factory=list)
    pair_dirs: dict = field(default_factory=dict)  # {(name_a, name_b): Path}
    pep_dir: Optional[Path] = None
    bed_dir: Optional[Path] = None


class JcviPipeline:
    """JCVI共享分析管道|JCVI Shared Analysis Pipeline"""

    def __init__(self, config: JcviBaseConfig, logger: Optional[JcviLogger] = None):
        self.config = config
        self.logger_obj = logger
        self.logger = logger.get_logger() if logger else None

    def run(self) -> PipelineResult:
        """运行共享管道步骤1-6|Run shared pipeline steps 1-6"""
        # 步骤1: 发现样本
        self._log_section("步骤1/6: 发现样本|Step 1/6: Discover samples")
        samples = discover_samples(
            self.config.input_dir,
            self.config.gff_ext,
            self.config.fa_ext,
            self.logger,
        )

        if len(samples) < 2:
            self.logger.error(
                f"至少需要2个有效样本|At least 2 valid samples required, found {len(samples)}"
            )
            return PipelineResult()

        self.logger.info(f"共发现 {len(samples)} 个样本|Found {len(samples)} samples total")

        if self.config.pairs:
            pair_list = self._resolve_pairs(samples)
        else:
            pair_list = list(combinations(samples, 2))

        self.logger.info(f"共 {len(pair_list)} 组两两比较|{len(pair_list)} pairwise comparisons")

        # 步骤2: 提取序列
        self._log_section("步骤2/6: 提取序列|Step 2/6: Extract sequences")
        pep_dir = Path(self.config.output_dir) / "01_pep"
        pep_dir.mkdir(parents=True, exist_ok=True)
        for prefix in samples:
            self._step_extract_pep(prefix, pep_dir)

        # 步骤3: GFF→BED
        self._log_section("步骤3/6: GFF转BED|Step 3/6: Convert GFF to BED")
        bed_dir = Path(self.config.output_dir) / "02_bed"
        bed_dir.mkdir(parents=True, exist_ok=True)
        for prefix in samples:
            self._step_gff2bed(prefix, bed_dir)
            self._step_bed_uniq(prefix, bed_dir)

        # 步骤4-6: 两两比较
        self._log_section("步骤4/6: 两两共线性分析|Step 4/6: Pairwise collinearity analysis")
        pair_dirs = {}
        for idx, (prefix_a, prefix_b) in enumerate(pair_list, 1):
            name_a = get_sample_name(prefix_a)
            name_b = get_sample_name(prefix_b)
            self.logger.info(f"\n[{idx}/{len(pair_list)}] {name_a} vs {name_b}")

            pair_dir = Path(self.config.output_dir) / "03_pairwise" / f"{name_a}_vs_{name_b}"
            pair_dir.mkdir(parents=True, exist_ok=True)

            success = self._run_pairwise(prefix_a, prefix_b, bed_dir, pep_dir, pair_dir)
            if success:
                pair_dirs[(name_a, name_b)] = pair_dir

        return PipelineResult(
            samples=samples,
            pair_list=pair_list,
            pair_dirs=pair_dirs,
            pep_dir=pep_dir,
            bed_dir=bed_dir,
        )

    # ---- 断点续传 + 命令执行 ----

    def _is_step_completed(self, output_file: str,
                           cwd: Optional[str] = None) -> bool:
        if cwd:
            output_file = os.path.join(cwd, output_file)
        return Path(output_file).exists() and Path(output_file).stat().st_size > 0

    def _run_step(self, step_name: str, output_file: str, cmd: List[str],
                  cwd: Optional[str] = None) -> bool:
        if self._is_step_completed(output_file, cwd):
            self.logger.info(f"  跳过已完成|Skipping completed: {step_name}")
            return True

        self.logger.info(f"  执行|Executing: {step_name}")
        self.logger.info(f"  命令|Command: {' '.join(cmd)}")

        try:
            run_env = None
            if cmd and '/envs/' in cmd[0] and cmd[0].endswith('python'):
                env_bin = os.path.dirname(cmd[0])
                run_env = os.environ.copy()
                run_env['PATH'] = env_bin + os.pathsep + run_env.get('PATH', '')

            result = subprocess.run(cmd, cwd=cwd, env=run_env, timeout=None)
            if result.returncode != 0:
                self.logger.error(
                    f"  {step_name} 失败 (exit code: {result.returncode})|{step_name} failed")
                return False

            if not self._is_step_completed(output_file, cwd):
                self.logger.error(
                    f"  {step_name} 完成但输出文件不存在|{step_name} completed but output missing: {output_file}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"  {step_name} 出错|{step_name} error: {e}")
            return False

    # ---- 子步骤 ----

    def _step_extract_pep(self, prefix: str, pep_dir: Path) -> bool:
        sample_name = get_sample_name(prefix)
        gff_file = prefix + self.config.gff_ext
        fa_file = prefix + self.config.fa_ext

        if self.config.dbtype == "nucl":
            out_file = str(pep_dir / f"{sample_name}.cds")
            cmd = ['gffread', gff_file, '-g', fa_file, '-x', out_file]
            step_desc = f"提取CDS序列 ({sample_name})"
        else:
            out_file = str(pep_dir / f"{sample_name}.pep")
            cmd = ['gffread', gff_file, '-g', fa_file, '-y', out_file]
            step_desc = f"提取蛋白质序列 ({sample_name})"

        return self._run_step(step_desc, out_file, cmd)

    def _step_gff2bed(self, prefix: str, bed_dir: Path) -> bool:
        sample_name = get_sample_name(prefix)
        gff_file = prefix + self.config.gff_ext
        bed_file = str(bed_dir / f"{sample_name}.bed")

        cmd = build_jcvi_command('jcvi.formats.gff', [
            'bed',
            '--type=' + self.config.gff_type,
            '--key=' + self.config.gff_key,
            gff_file,
            '-o', bed_file,
        ], self.config.conda_env)

        return self._run_step(f"GFF转BED ({sample_name})", bed_file, cmd)

    def _step_bed_uniq(self, prefix: str, bed_dir: Path) -> bool:
        sample_name = get_sample_name(prefix)
        bed_file = f"{sample_name}.bed"
        uniq_file = f"{get_jcvi_stem(sample_name)}.uniq.bed"

        cmd = build_jcvi_command('jcvi.formats.bed', ['uniq', bed_file],
                                 self.config.conda_env)

        return self._run_step(f"BED去重 ({sample_name})", uniq_file, cmd,
                              cwd=str(bed_dir))

    def _run_pairwise(self, prefix_a: str, prefix_b: str,
                      bed_dir: Path, pep_dir: Path, pair_dir: Path) -> bool:
        """
        运行两两共线性分析(步骤4-6)|Run pairwise collinearity analysis (steps 4-6)

        流程:
            1. 序列比对 (LAST/diamond)
            2. Blastfilter过滤
            3. Synteny scan → .anchors + .lifted.anchors
        """
        name_a = get_sample_name(prefix_a)
        name_b = get_sample_name(prefix_b)
        stem_a = get_jcvi_stem(name_a)
        stem_b = get_jcvi_stem(name_b)
        suffix = ".cds" if self.config.dbtype == "nucl" else ".pep"
        pprefix = f"{stem_a}.{stem_b}"

        # 创建符号链接
        for name in [name_a, name_b]:
            stem = get_jcvi_stem(name)
            bed_src = str((bed_dir / f"{get_jcvi_stem(name)}.uniq.bed").resolve())
            pep_src = str((pep_dir / f"{name}{suffix}").resolve())
            bed_link = str(pair_dir / f"{stem}.bed")
            pep_link = str(pair_dir / f"{name}{suffix}")

            for src, link in [(bed_src, bed_link), (pep_src, pep_link)]:
                if Path(link).exists() or Path(link).is_symlink():
                    Path(link).unlink()
                os.symlink(src, link)

        # 子步骤1: 序列比对
        last_file = f"{pprefix}.last"
        if self.config.align_soft == "last":
            cmd = build_jcvi_command('jcvi.apps.align', [
                'last',
                f"{name_b}{suffix}", f"{name_a}{suffix}",
                f"--dbtype={self.config.dbtype}",
                f"--cpus={self.config.threads}",
            ], self.config.conda_env)
            if not self._run_step(f"LAST比对|LAST alignment ({name_a} vs {name_b})",
                                  last_file, cmd, cwd=str(pair_dir)):
                return False
        elif self.config.align_soft == "diamond_blastp":
            cmd = build_jcvi_command('jcvi.apps.align', [
                'blast',
                f"{name_b}{suffix}", f"{name_a}{suffix}",
                f"--dbtype={self.config.dbtype}",
                f"--cpus={self.config.threads}",
            ], self.config.conda_env)
            if not self._run_step(f"Blast比对|Blast alignment ({name_a} vs {name_b})",
                                  last_file, cmd, cwd=str(pair_dir)):
                return False

        # 子步骤2: Blastfilter
        filtered_last = f"{last_file}.filtered"
        filter_args = [last_file, f"--cscore={self.config.cscore}"]
        if self.config.no_strip_names:
            filter_args.append("--no_strip_names")
        cmd = build_jcvi_command('jcvi.compara.blastfilter', filter_args,
                                 self.config.conda_env)
        if not self._run_step(f"过滤比对结果|Filter alignment ({name_a} vs {name_b})",
                              filtered_last, cmd, cwd=str(pair_dir)):
            return False

        # 子步骤3: Synteny scan
        anchors_file = f"{pprefix}.anchors"
        use_liftover = getattr(self.config, 'liftover', True)
        lifted_anchors_file = f"{pprefix}.lifted.anchors" if use_liftover else f"{pprefix}.anchors"
        scan_args = [
            'scan',
            filtered_last,
            anchors_file,
            f"--min_size={self.config.min_size}",
            "--dist=20",
        ]
        if use_liftover:
            scan_args.append(f"--liftover={last_file}")
        if self.config.no_strip_names:
            scan_args.append("--no_strip_names")
        cmd = build_jcvi_command('jcvi.compara.synteny', scan_args, self.config.conda_env)
        if not self._run_step(f"共线性扫描|Synteny scan ({name_a} vs {name_b})",
                              lifted_anchors_file, cmd, cwd=str(pair_dir)):
            return False

        return True

    # ---- 工具方法 ----

    def _resolve_pairs(self, samples: List[str]) -> List[tuple]:
        sample_map = {get_sample_name(s): s for s in samples}
        result = []
        for pair_str in self.config.pairs:
            parts = pair_str.split(',')
            if len(parts) != 2:
                self.logger.warning(f"  无效配对格式, 跳过|Invalid pair format, skip: {pair_str}")
                continue
            name_a, name_b = parts[0].strip(), parts[1].strip()
            if name_a not in sample_map or name_b not in sample_map:
                self.logger.warning(f"  样本不存在, 跳过|Sample not found, skip: {name_a},{name_b}")
                continue
            result.append((sample_map[name_a], sample_map[name_b]))
        return result

    def _log_section(self, title: str):
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(title)
        self.logger.info("=" * 60)
