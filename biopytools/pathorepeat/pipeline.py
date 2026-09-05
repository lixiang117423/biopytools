"""pathorepeat 流程编排模块|pathorepeat Pipeline Module

四步:建库建模→屏蔽→TEsorter分类→汇总;步骤级断点续传;批量样品循环
|4 steps: model→mask→classify→summarize; step-level resume; batch loop
"""

import os
import re
import shutil
import subprocess
from typing import Dict, List, Optional

from biopytools.common.conda_runner import build_conda_command

from . import __version__
from .config import MASKING_FLAGS, PathorepeatConfig
from .report import (
    write_batch_summary, write_effector_overlap, write_families_classified,
    write_repeat_summary,
)
from .utils import (
    fasta_ids, load_effector_regions, parse_repeatmasker_out,
    parse_repeatmasker_tbl, parse_tesorter_cls_tsv,
)

STEP_DIRS = {
    'modeler': '01_modeler',
    'masker': '02_masker',
    'tesorter': '03_tesorter',
    'summary': '04_summary',
    'info': '00_pipeline_info',
    'logs': '99_logs',
}

# 版本探测正则(须匹配真实二进制输出,如 "RepeatModeler version 2.0.9")
# |Version probe patterns (must match real binary outputs, e.g. "RepeatModeler version 2.0.9")
VERSION_PATTERNS = {
    'repeatmodeler': r'RepeatModeler\s+(?:v(?:ersion)?\s+)?([\d.]+)',
    'repeatmasker': r'RepeatMasker\s+version\s+([\d.]+)',
    'tesorter': r'([\d.]+)',
}


def run_command(cmd: List[str], logger, description: str,
                working_dir: Optional[str] = None,
                extra_env: Optional[Dict[str, str]] = None) -> bool:
    """执行外部命令(INFO 先记录完整命令,§2.2.1)|Run external command

    所有命令已经 build_conda_command 包装(含 --no-capture-output);
    extra_env 以继承父环境方式注入(如 FAMDB_DIR)
    |extra_env is merged over the inherited environment (e.g. FAMDB_DIR)
    """
    logger.info(f"执行|Executing: {description}")
    logger.info(f"命令|Command: {' '.join(cmd)}")
    env = {**os.environ, **extra_env} if extra_env else None
    try:
        result = subprocess.run(cmd, shell=False, cwd=working_dir, env=env,
                                capture_output=True, text=True)
    except FileNotFoundError as e:
        logger.error(f"工具不存在|Tool not found: {e}")
        return False
    if result.returncode != 0:
        tail = (result.stderr or '')[-2000:]
        logger.error(f"命令失败(rc={result.returncode})|Command failed "
                     f"(rc={result.returncode}): {tail}")
        return False
    return True


class PathorepeatPipeline:
    """pathorepeat 流程编排器|pathorepeat pipeline orchestrator"""

    def __init__(self, config: PathorepeatConfig, logger):
        self.config = config
        self.logger = logger
        self.dirs = {k: os.path.join(config.output_dir, v)
                     for k, v in STEP_DIRS.items()}
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)

    # ---------- 通用|Common ----------

    def _step_done(self, marker: str) -> bool:
        """断点续传:marker 存在且非空即跳过(§10.2)|Resume via done-marker"""
        return (self.config.skip_completed and os.path.exists(marker)
                and os.path.getsize(marker) > 0)

    def _modeler_run_dir(self, name: str) -> str:
        """RepeatModeler 工作目录(按样品隔离)|RM working dir (per-sample)"""
        return os.path.join(self.dirs['modeler'], f'{name}_rm_run')

    def _lib_path(self, name: str) -> str:
        """de novo 库路径(done-marker)|De novo library path (done-marker)"""
        return os.path.join(self._modeler_run_dir(name),
                            f'{name}_db-families.fa')

    def probe_versions(self) -> Dict[str, str]:
        """探测三工具版本(登录节点安全)|Probe tool versions"""
        probes = {
            'repeatmodeler': (self.config.repeatmodeler_path, ['-version']),
            'repeatmasker': (self.config.repeatmasker_path, ['-v']),
            'tesorter': (self.config.tesorter_path, ['-v']),
        }
        versions = {}
        for key, (path, args) in probes.items():
            cmd = build_conda_command(path, args)
            self.logger.info(f"执行|Executing: 探测版本|version probe ({key})")
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            try:
                result = subprocess.run(
                    cmd, shell=False,
                    capture_output=True, text=True, timeout=120)
                text = (result.stdout or '') + (result.stderr or '')
                m = re.search(VERSION_PATTERNS[key], text)
                versions[key] = m.group(1) if m else 'unknown'
            except (subprocess.SubprocessError, OSError):
                versions[key] = 'unknown'
        return versions

    def write_software_versions(self, versions: Dict[str, str]) -> None:
        """写 00_pipeline_info/software_versions.yml(§12.5)|Write versions yml"""
        path = os.path.join(self.dirs['info'], 'software_versions.yml')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(f"pathorepeat_module: {__version__}\n")
            fh.write(f"repeatmodeler: {versions.get('repeatmodeler', 'unknown')}\n")
            fh.write(f"repeatmasker: {versions.get('repeatmasker', 'unknown')}\n")
            fh.write(f"tesorter: {versions.get('tesorter', 'unknown')}\n")

    # ---------- Step 1|建库+建模 ----------

    def _modeler_fallback(self, name: str, run_dir: str) -> Optional[str]:
        """分类已尝试失败时降级采用未分类 consensi.fa|Fall back to unclassified consensi.fa

        判据:consensi.fa 非空 且 rmod.log 中出现过 RepeatClassifier
        (RM2 主流程已跑完、仅分类失败——典型原因 FAMDB_DIR/Dfam 数据未配置)。
        --no-skip-completed(强制重建)时不降级,保持重建语义。
        |Criteria: non-empty consensi.fa + RepeatClassifier seen in rmod.log
        (RM2 finished modeling, only classification failed). Skipped on force
        rerun (--no-skip-completed).
        """
        if not self.config.skip_completed:
            return None
        cons = os.path.join(run_dir, 'consensi.fa')
        rmod_log = os.path.join(run_dir, f'{name}_db-rmod.log')
        if not (os.path.exists(cons) and os.path.getsize(cons) > 0):
            return None
        classifier_attempted = False
        if os.path.exists(rmod_log):
            with open(rmod_log, encoding='utf-8',
                      errors='ignore') as fh:
                classifier_attempted = 'RepeatClassifier' in fh.read()
        if not classifier_attempted:
            return None
        marker = self._lib_path(name)
        shutil.copyfile(cons, marker)
        self.logger.warning(
            f"RepeatModeler 分类失败(常见原因:FAMDB_DIR/Dfam 数据未配置),"
            f"降级采用未分类库 consensi.fa 继续流程,家族分类由 TEsorter 承担;"
            f"配置 Dfam 后可删除 {marker} 重跑以获得 RM2 自带分类"
            f"|RepeatModeler classification failed (usually FAMDB_DIR/Dfam data "
            f"missing); falling back to unclassified consensi.fa, classification "
            f"delegated to TEsorter; configure Dfam and remove {marker} to rerun "
            f"with RM2 built-in classification")
        return marker

    def _run_modeler(self, genome: str, name: str) -> Optional[str]:
        """BuildDatabase + RepeatModeler;返回库路径或 None|Returns lib path or None"""
        marker = self._lib_path(name)
        if self._step_done(marker):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: "
                             f"modeler ({name})")
            return marker
        run_dir = self._modeler_run_dir(name)
        os.makedirs(run_dir, exist_ok=True)
        fallback = self._modeler_fallback(name, run_dir)
        if fallback:
            return fallback
        db_name = f'{name}_db'

        cmd = build_conda_command(
            self.config.build_database_path, ['-name', db_name, genome])
        if not run_command(cmd, self.logger, f"BuildDatabase 建库|Building "
                                             f"database ({name})",
                           working_dir=run_dir):
            return None

        args = ['-database', db_name, '-threads', str(self.config.threads),
                '-dir', run_dir]
        if self.config.ltr_struct:
            args.append('-LTRStruct')
        cmd = build_conda_command(self.config.repeatmodeler_path, args)
        desc = f"RepeatModeler 建模|RepeatModeler de novo ({name})"
        rm_env = ({'FAMDB_DIR': self.config.famdb_dir}
                  if self.config.famdb_dir else None)
        if not run_command(cmd, self.logger, desc, working_dir=run_dir,
                           extra_env=rm_env):
            return None
        if not os.path.exists(marker):
            self.logger.error(f"RepeatModeler 未产出库文件|No library produced: "
                              f"{marker}")
            return None
        return marker

    # ---------- Step 2|屏蔽 ----------

    def _rename_masker_outputs(self, genome: str, name: str) -> None:
        """RepeatMasker 原生输出名改为 {name}.* 前缀|Rename RM outputs"""
        base = os.path.basename(genome)
        mapping = {
            f'{base}.masked': f'{name}_masked.fa',
            f'{base}.out': f'{name}.out',
            f'{base}.gff': f'{name}.gff',
            f'{base}.tbl': f'{name}.tbl',
        }
        for src_name, dst_name in mapping.items():
            src = os.path.join(self.dirs['masker'], src_name)
            dst = os.path.join(self.dirs['masker'], dst_name)
            if os.path.exists(src):
                os.replace(src, dst)

    def _run_masker(self, genome: str, name: str, lib: str) -> Optional[str]:
        """RepeatMasker 软屏蔽;返回 masked 路径或 None|Returns masked path"""
        masker_dir = self.dirs['masker']
        marker = os.path.join(masker_dir, f'{name}_masked.fa')
        if self._step_done(marker):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: "
                             f"masker ({name})")
            return marker
        args = ['-pa', str(self.config.threads), '-lib', lib]
        args += MASKING_FLAGS[self.config.masking_mode]
        args += ['-gff', '-dir', masker_dir, genome]
        cmd = build_conda_command(self.config.repeatmasker_path, args)
        if not run_command(cmd, self.logger,
                           f"RepeatMasker 屏蔽|RepeatMasker masking ({name})"):
            return None
        self._rename_masker_outputs(genome, name)
        if not os.path.exists(marker):
            self.logger.error(f"RepeatMasker 未产出 masked 文件|No masked output: "
                              f"{marker}")
            return None
        return marker

    # ---------- Step 3|TEsorter 分类(失败降级) ----------

    def _run_tesorter(self, lib: str, name: str) -> Optional[str]:
        """TEsorter 分类;失败返回 None(降级继续)|Returns cls.tsv path or None"""
        db = self.config.tesorter_db
        prefix = os.path.join(self.dirs['tesorter'],
                              f'{name}_db-families_{db}')
        marker = f'{prefix}.cls.tsv'
        if self._step_done(marker):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: "
                             f"tesorter ({name})")
            return marker
        args = [lib, '-db', db, '-p', str(self.config.threads),
                '-pre', prefix, '-tmp',
                os.path.join(self.config.output_dir, 'tmp')]
        if self.config.db_hmm:
            args += ['--db-hmm', self.config.db_hmm]
        cmd = build_conda_command(self.config.tesorter_path, args)
        if not run_command(cmd, self.logger,
                           f"TEsorter 分类|TEsorter classification ({name})"):
            self.logger.warning(
                f"TEsorter 失败,降级为全 unknown 继续汇总|TEsorter failed; "
                f"degrading to all-unknown and continuing ({name})")
            return None
        if not os.path.exists(marker):
            self.logger.warning(
                f"TEsorter 未产出分类表,降级继续|No cls.tsv produced; degrading "
                f"({name})")
            return None
        return marker

    # ---------- Step 4|汇总 ----------

    def _run_summary(self, genome: str, name: str, cls_path: Optional[str],
                     lib: str) -> Dict:
        """汇总 + 可选 effector 交叉检查;返回指标|Summary and metrics"""
        summary_dir = self.dirs['summary']
        tbl = parse_repeatmasker_tbl(os.path.join(self.dirs['masker'],
                                                  f'{name}.tbl'))
        hits = parse_repeatmasker_out(os.path.join(self.dirs['masker'],
                                                   f'{name}.out'))
        cls_map = (parse_tesorter_cls_tsv(cls_path) if cls_path else {})
        write_repeat_summary(
            os.path.join(summary_dir, f'{name}_repeat_summary.tsv'),
            tbl, hits, cls_map, lib, genome_fasta=genome)
        write_families_classified(
            os.path.join(summary_dir, f'{name}_families_classified.tsv'),
            lib, cls_map)
        if self.config.effector_bed:
            regions = load_effector_regions(self.config.effector_bed)
            write_effector_overlap(
                os.path.join(summary_dir, f'{name}_effector_overlap.tsv'),
                regions, hits)
        elif self.config.effector_gff:
            regions = load_effector_regions(self.config.effector_gff)
            write_effector_overlap(
                os.path.join(summary_dir, f'{name}_effector_overlap.tsv'),
                regions, hits)

        families = fasta_ids(lib)
        classified = [f for f in families if f in cls_map]
        masked_pct = tbl.get('masked_pct')
        return {
            'sample': name, 'status': 'ok' if cls_path else 'degraded',
            'masked_pct': masked_pct,
            'n_families': len(families),
            'classified_pct': (len(classified) / len(families) * 100)
                              if families else 0.0,
            'reason': '' if cls_path else
                      'TEsorter分类失败,unknown处理|TEsorter degraded to unknown',
        }

    def _cleanup_tmp(self) -> None:
        """清理 TEsorter 临时目录(§12.4,docs 承诺)|Remove tmp dir (§12.4)

        ignore_errors 保证目录不存在时幂等|ignore_errors makes it idempotent
        """
        self.logger.info("清理临时目录|Cleaning tmp directory")
        shutil.rmtree(os.path.join(self.config.output_dir, 'tmp'),
                      ignore_errors=True)

    # ---------- 入口|Entry ----------

    def run_for_genome(self, genome: str, name: str) -> Dict:
        """单样品完整四步;失败不抛错返回 failed|One sample, 4 steps, no raise"""
        try:
            lib = self._run_modeler(genome, name)
            if lib is None:
                return {'sample': name, 'status': 'failed',
                        'masked_pct': None, 'n_families': None,
                        'classified_pct': None,
                        'reason': 'RepeatModeler失败|RepeatModeler failed'}
            masked = self._run_masker(genome, name, lib)
            if masked is None:
                return {'sample': name, 'status': 'failed',
                        'masked_pct': None, 'n_families': None,
                        'classified_pct': None,
                        'reason': 'RepeatMasker失败|RepeatMasker failed'}
            cls_path = self._run_tesorter(lib, name)
            return self._run_summary(genome, name, cls_path, lib)
        except (OSError, ValueError) as e:
            self.logger.error(f"样品处理异常|Sample error ({name}): {e}")
            return {'sample': name, 'status': 'failed', 'masked_pct': None,
                    'n_families': None, 'classified_pct': None, 'reason': str(e)}

    def run(self, genomes: List[str]) -> bool:
        """批量入口:逐样品循环+batch_summary;全部成功返回 True|Batch entry"""
        versions = self.probe_versions()
        self.logger.info(
            f"版本|Versions: RepeatModeler {versions['repeatmodeler']} / "
            f"RepeatMasker {versions['repeatmasker']} / "
            f"TEsorter {versions['tesorter']}")
        self.write_software_versions(versions)

        results: List[Dict] = []
        for i, genome in enumerate(genomes, 1):
            name = self.config.sample_name(genome)
            self.logger.info(f"开始样品|Starting sample ({i}/{len(genomes)}): "
                             f"{name}")
            results.append(self.run_for_genome(genome, name))
            self.logger.info(f"样品完成|Sample done ({i}/{len(genomes)}): "
                             f"{name} → {results[-1]['status']}")

        if len(genomes) > 1:
            write_batch_summary(
                os.path.join(self.dirs['summary'], 'batch_summary.tsv'), results)
        self._cleanup_tmp()
        ok = all(r['status'] in ('ok', 'degraded') for r in results)
        if not ok:
            failed = [r['sample'] for r in results if r['status'] == 'failed']
            self.logger.error(f"存在失败样品|Failed samples: "
                              f"{', '.join(failed)}")
        return ok
