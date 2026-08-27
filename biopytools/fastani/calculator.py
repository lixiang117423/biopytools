"""fastANI核心计算模块|fastANI Calculator Module

编排: 收集输入 → 写列表 → 运行fastANI → 解析.out → 矩阵+最近邻 → 版本信息
|Orchestrate: collect → write lists → run fastANI → parse .out → matrix + nearest → versions
"""

import os
import subprocess
from typing import Dict, List, Tuple

from ..common.conda_runner import build_conda_command

from .utils import collect_genome_files, genome_name


class FastaniCalculator:
    """fastANI计算器|fastANI Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.pipeline_info_dir = os.path.join(config.output_dir, '00_pipeline_info')
        self.fastani_dir = os.path.join(config.output_dir, '01_fastani')
        self.results_dir = os.path.join(config.output_dir, '02_results')
        self.logs_dir = os.path.join(config.output_dir, '99_logs')
        self.out_file = os.path.join(self.fastani_dir, 'fastani.out')
        self._ql_file = None
        self._rl_file = None
        self._n_query = 0
        self._n_reference = 0

    def _prepare_dirs(self):
        """创建输出子目录|Create output subdirectories"""
        for d in (self.pipeline_info_dir, self.fastani_dir,
                  self.results_dir, self.logs_dir):
            os.makedirs(d, exist_ok=True)

    # ---------- 输入编排|Input orchestration ----------

    @staticmethod
    def _check_no_duplicates(paths: List[str], side: str):
        """侧内基因组名查重|Duplicate genome-name check within one side"""
        seen = {}
        for p in paths:
            name = genome_name(p)
            if name in seen:
                raise ValueError(
                    f"基因组名重复|Duplicate genome name: {name} "
                    f"({seen[name]} 与|and {p});同侧文件主名必须唯一|"
                    f"stems must be unique within {side}")
            seen[name] = p

    def collect_inputs(self) -> Tuple[List[str], List[str]]:
        """收集双侧基因组|Collect genomes for both sides

        Returns: (query_paths, reference_paths);all-vs-all两侧相同|both sides identical
        """
        if self.config.all_vs_all:
            genomes = collect_genome_files(self.config.input)
            if len(genomes) < 2:
                raise ValueError("all-vs-all至少需要2个基因组|all-vs-all needs at least "
                                 f"2 genomes, got {len(genomes)}")
            self._check_no_duplicates(genomes, 'input')
            return genomes, genomes

        query = collect_genome_files(self.config.query)
        reference = collect_genome_files(self.config.reference)
        if not query:
            raise ValueError("query侧未收集到基因组|No genomes collected for query")
        if not reference:
            raise ValueError("reference侧未收集到基因组|No genomes collected for reference")
        self._check_no_duplicates(query, 'query')
        self._check_no_duplicates(reference, 'reference')
        # 跨侧同名仅提示(自比对真实运行约100%)|Cross-side overlap: warn only (self pairs run ~100%)
        q_names = {genome_name(p) for p in query}
        overlap = q_names & {genome_name(p) for p in reference}
        if overlap:
            self.logger.warning(
                f"query与reference存在同名基因组|Overlapping names across sides: "
                f"{sorted(overlap)};自比对将真实运行(约100%),计入矩阵|self pairs "
                f"will actually run (~100%) and are included in the matrix")
        return query, reference

    def write_list_files(self, query_paths: List[str],
                         reference_paths: List[str]) -> Tuple[str, str]:
        """写fastANI --ql/--rl列表|Write fastANI --ql/--rl list files"""
        os.makedirs(self.fastani_dir, exist_ok=True)
        ql_file = os.path.join(self.fastani_dir, 'genome_list_ql.txt')
        rl_file = os.path.join(self.fastani_dir, 'genome_list_rl.txt')
        for path_list, out_file in ((query_paths, ql_file),
                                    (reference_paths, rl_file)):
            with open(out_file, 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(path_list) + '\n')
        self._ql_file = ql_file
        self._rl_file = rl_file
        self.logger.info(
            f"基因组列表已写|Genome lists written: query={len(query_paths)}, "
            f"reference={len(reference_paths)}")
        return ql_file, rl_file

    # ---------- 运行|Run ----------

    def _input_lists_newer(self) -> bool:
        """输入列表是否比结果新(mtime)|Whether any input list is newer than result"""
        if not os.path.exists(self.out_file):
            return False
        out_mtime = os.path.getmtime(self.out_file)
        for lst in (self._ql_file, self._rl_file):
            if lst and os.path.exists(lst) and os.path.getmtime(lst) > out_mtime:
                return True
        return False

    def _fastani_cmd(self, ql_file, rl_file, out_file) -> List[str]:
        """构建 fastANI 命令|Build the fastANI command"""
        return build_conda_command(self.config.fastani_path, [
            '--ql', ql_file, '--rl', rl_file,
            '-o', out_file, '--matrix',
            '-t', str(self.config.threads),
            '-k', str(self.config.kmer),
            '--fragLen', str(self.config.frag_len),
            '--minFraction', str(self.config.min_fraction),
        ])

    def _run_single(self, ql_file, rl_file, out_file) -> bool:
        """运行单次 fastANI|Run one fastANI invocation"""
        cmd = self._fastani_cmd(ql_file, rl_file, out_file)
        self.logger.info("执行|Executing: fastANI全基因组ANI计算|fastANI whole-genome ANI")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False,
                                    capture_output=True, text=True)
        except FileNotFoundError as e:
            self.logger.error(f"fastANI未找到|fastANI not found: {e}")
            return False
        if result.returncode != 0:
            self.logger.error(
                f"fastANI运行失败|fastANI failed (exit {result.returncode}): "
                f"{(result.stderr or '').strip()[:2000]}")
            return False
        if not os.path.exists(out_file):
            self.logger.error(
                "fastANI成功退出但输出文件缺失|fastANI exited 0 but output "
                f"file missing: {out_file}")
            return False
        return True

    def _run_iterated(self, query_paths) -> bool:
        """逐轮 1-vs-all,reference 也分批:内存 = 1 query + 一批 reference 草图
        |Iterated 1-vs-all with reference batching: memory = 1 query + batch
        of reference sketches (resume per batch)

        fastANI 会把 reference 列表全部 sketch 驻留内存(线性累积),
        因此把 reference 拆成小批(默认 50),每轮只加载一批,内存可控。
        |fastANI keeps all reference sketches in RAM (linear growth), so we
        split references into small batches (default 50) per round.

        Returns:
            bool: 全部批次完成|all batches completed
        """
        n = len(query_paths)
        ref_batch = getattr(self.config, 'ref_batch_size', 50)
        self.logger.warning(
            f"基因组数 {n} > 阈值 {self.config.iterated_threshold},"
            f"使用逐轮 1-vs-all + reference 分批(每批 {ref_batch},内存友好)|"
            f"{n} genomes over threshold; iterated 1-vs-all with reference "
            f"batching (batch={ref_batch}, low memory)")

        batch_dir = os.path.join(self.fastani_dir, 'batches')
        os.makedirs(batch_dir, exist_ok=True)

        # 汇总文件(先存在则整体跳过)|aggregate output (skip if complete)
        if os.path.exists(self.out_file) and os.path.getsize(self.out_file) > 0:
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: fastANI "
                f"({self.out_file} 已存在|exists)")
            return True
        # 汇总文件存在但可能残缺(上次汇总中断)→ 重建|stale aggregate → rebuild
        if os.path.exists(self.out_file):
            os.remove(self.out_file)

        total_batches = 0
        done = 0
        for i, q in enumerate(query_paths):
            ql = os.path.join(batch_dir, f'ql_{i:04d}.txt')
            with open(ql, 'w', encoding='utf-8') as fh:
                fh.write(q + '\n')
            # reference 分批|reference batching
            for j in range(0, n, ref_batch):
                refs = query_paths[j:j + ref_batch]
                rl = os.path.join(batch_dir, f'rl_{i:04d}_{j:04d}.txt')
                with open(rl, 'w', encoding='utf-8') as fh:
                    fh.write('\n'.join(refs) + '\n')
                out_batch = os.path.join(
                    batch_dir, f'batch_{i:04d}_{j:04d}.out')
                total_batches += 1
                if os.path.exists(out_batch) and os.path.getsize(out_batch) > 0:
                    done += 1
                    continue
                if not self._run_single(ql, rl, out_batch):
                    self.logger.error(
                        f"批次失败|Batch ({i},{j}) failed: {q}; 重跑将从中断处"
                        f"继续|rerun resumes from here")
                    return False
                done += 1
                if done % max(1, total_batches // 10) == 0 or done == total_batches:
                    self.logger.info(
                        f"批次进度|Batch progress: {done}/{total_batches}")

        # 汇总所有批次 → fastani.out|aggregate batches
        with open(self.out_file, 'w', encoding='utf-8') as fh:
            for i in range(n):
                for j in range(0, n, ref_batch):
                    out_batch = os.path.join(
                        batch_dir, f'batch_{i:04d}_{j:04d}.out')
                    if not os.path.exists(out_batch):
                        self.logger.error(f"批次缺失|Missing batch: {out_batch}")
                        return False
                    with open(out_batch, encoding='utf-8') as bfh:
                        fh.write(bfh.read())
        self.logger.info(
            f"遍历完成,已汇总|Iterated runs done, aggregated "
            f"({done}/{total_batches} batches)")
        return True


    def run_fastani(self, query_paths: List[str] = None) -> bool:
        """运行fastANI(断点续传)|Run fastANI (checkpoint resume)

        须先调用write_list_files(内部依赖_ql_file/_rl_file)
        |Requires write_list_files first (uses _ql_file/_rl_file)
        """
        # 遍历模式分流|iterated mode branch
        if self.config.all_vs_all and self.config.iterated \
                and self._n_query > self.config.iterated_threshold:
            return self._run_iterated(query_paths)

        # 输入变化(任一列表比结果新)→ 作废旧结果重跑|input changed (any list newer
        # than result) → invalidate stale output and re-run
        if os.path.exists(self.out_file) and self._input_lists_newer():
            self.logger.warning(
                "输入基因组列表已变化,作废旧结果重跑|Input lists changed, "
                "invalidating stale fastANI output")
            os.remove(self.out_file)

        if os.path.exists(self.out_file) and os.path.getsize(self.out_file) > 0:
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: fastANI "
                f"({self.out_file} 已存在|exists)")
            return True

        if not self._run_single(self._ql_file, self._rl_file, self.out_file):
            return False

        if os.path.getsize(self.out_file) == 0:
            # 空输出+rc0=全部配对<80%|Empty output with rc 0 = all pairs below 80%
            self.logger.warning(
                "fastANI输出为空(所有配对ANI<80%或无可比对片段),矩阵将记NA|"
                "fastANI output empty (all pairs <80% ANI); matrix cells will be NA")
            return True

        self.logger.info("fastANI运行完成|fastANI run completed")
        return True

    # ---------- 编排|Orchestration ----------

    def _generate_software_versions(self):
        """生成software_versions.yml(§12.5)|Generate software_versions.yml"""
        import yaml

        versions = {}
        try:
            cmd = build_conda_command(self.config.fastani_path, ['-v'])
            self.logger.info("执行|Executing: fastANI版本探测|fastANI version probe")
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            ver = (result.stdout or result.stderr or '').strip().split('\n')[0]
            versions['fastANI'] = {'version': ver or 'unknown',
                                   'path': self.config.fastani_path}
        except Exception as e:
            self.logger.warning(f"获取fastANI版本失败|Failed to get fastANI version: {e}")
            versions['fastANI'] = {'version': 'unknown',
                                   'path': self.config.fastani_path}

        from . import __version__
        info = {
            'pipeline': {'name': 'biopytools fastani', 'version': __version__},
            'tools': versions,
            'parameters': {
                'mode': 'all-vs-all' if self.config.all_vs_all else 'query-vs-ref',
                'threads': self.config.threads,
                'kmer': self.config.kmer,
                'frag_len': self.config.frag_len,
                'min_fraction': self.config.min_fraction,
                'n_query': self._n_query,
                'n_reference': self._n_reference,
            },
        }
        out = os.path.join(self.pipeline_info_dir, 'software_versions.yml')
        try:
            with open(out, 'w', encoding='utf-8') as fh:
                yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
            self.logger.info(f"软件版本信息已保存|Software versions saved: {out}")
        except Exception as e:
            self.logger.warning(f"保存软件版本信息失败|Failed to save versions: {e}")

    def run(self) -> bool:
        """主流程|Main pipeline"""
        self._prepare_dirs()
        query_paths, ref_paths = self.collect_inputs()
        self._n_query, self._n_reference = len(query_paths), len(ref_paths)
        query_names = [genome_name(p) for p in query_paths]
        ref_names = [genome_name(p) for p in ref_paths]
        self.write_list_files(query_paths, ref_paths)

        if not self.run_fastani(query_paths):
            return False

        records, malformed = parse_fastani_out(self.out_file)
        if malformed:
            self.logger.warning(
                f"跳过{malformed}行格式异常记录|Skipped {malformed} malformed lines")

        build_ani_matrix(records, query_names, ref_names,
                         self.config.all_vs_all,
                         os.path.join(self.results_dir, 'ani_matrix.tsv'))
        build_nearest_table(records, query_names, ref_names,
                            self.config.all_vs_all,
                            os.path.join(self.results_dir, 'nearest_genome.tsv'),
                            logger=self.logger)
        self._generate_software_versions()

        n_pairs = len(query_names) * len(ref_names)
        if self.config.all_vs_all:
            n_pairs -= len(query_names)  # 排除自身|exclude self
        na_pairs = max(n_pairs - len(records), 0)
        self.logger.info(
            f"完成|Completed: {len(records)}对报告|pairs reported, "
            f"约{na_pairs}对<80%记NA|~{na_pairs} pairs <80% marked NA")
        return True


# ---------- 模块级纯函数|Module-level pure functions ----------


def parse_fastani_out(out_file: str) -> Tuple[List[Dict], int]:
    """解析fastANI .out(5列)|Parse fastANI .out (5 columns)

    列|Columns: query reference ANI matched_frags total_frags
    tab分隔优先,空白分隔回退;坏行计数不中断|tab first, whitespace fallback;
    malformed lines counted, not fatal
    """
    records: List[Dict] = []
    malformed = 0
    with open(out_file, encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 5:
                parts = line.split()
            if len(parts) != 5:
                malformed += 1
                continue
            try:
                records.append({
                    'query': genome_name(parts[0]),
                    'reference': genome_name(parts[1]),
                    'ani': float(parts[2]),
                    'matched_frags': int(parts[3]),
                    'total_frags': int(parts[4]),
                })
            except ValueError:
                malformed += 1
    return records, malformed


def build_ani_matrix(records: List[Dict], query_names: List[str],
                     ref_names: List[str], all_vs_all: bool, out_file: str) -> None:
    """生成ANI矩阵TSV|Write ANI matrix TSV

    双向平均(同官方--matrix);单向取单向;<80%记NA;all-vs-all对角100
    |Bidirectional average (as official --matrix); one-way kept; <80% NA;
    all-vs-all diagonal 100
    """
    pair = {(r['query'], r['reference']): r for r in records}
    lines = ['genome\t' + '\t'.join(ref_names)]
    for q in query_names:
        cells = []
        for r in ref_names:
            if all_vs_all and q == r:
                cells.append('100.0000')
                continue
            fwd = pair.get((q, r))
            rev = pair.get((r, q))
            if fwd and rev:
                cells.append(f"{(fwd['ani'] + rev['ani']) / 2:.4f}")
            elif fwd or rev:
                rec = fwd if fwd else rev
                cells.append(f"{rec['ani']:.4f}")
            else:
                cells.append('NA')
        lines.append(q + '\t' + '\t'.join(cells))
    with open(out_file, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')


def build_nearest_table(records: List[Dict], query_names: List[str],
                        ref_names: List[str], all_vs_all: bool,
                        out_file: str, logger=None) -> None:
    """生成最近邻表|Write nearest-neighbor table

    每query取ANI最高的reference(all-vs-all排除自身);全NA行告警且排最后;
    按ANI降序|Best ANI reference per query (self excluded in all-vs-all);
    all-NA rows warned and sorted last; ANI descending
    """
    by_query: Dict[str, List[Dict]] = {}
    for r in records:
        by_query.setdefault(r['query'], []).append(r)

    rows = []
    na_genomes = []
    for q in query_names:
        candidates = [r for r in by_query.get(q, [])
                      if not (all_vs_all and r['reference'] == q)]
        if not candidates:
            rows.append((q, 'NA', 'NA', 'NA', 'NA', 'NA'))
            na_genomes.append(q)
            continue
        best = max(candidates, key=lambda r: (r['ani'], r['matched_frags']))
        af = best['matched_frags'] / best['total_frags'] if best['total_frags'] else 0.0
        rows.append((q, best['reference'], f"{best['ani']:.4f}", f"{af:.4f}",
                     str(best['matched_frags']), str(best['total_frags'])))

    header = ('genome\tnearest_genome\tani\taligned_fraction\t'
              'matched_frags\ttotal_frags')
    rows.sort(key=lambda row: (row[2] == 'NA',
                               -float(row[2]) if row[2] != 'NA' else 0.0))
    with open(out_file, 'w', encoding='utf-8') as fh:
        fh.write(header + '\n')
        for row in rows:
            fh.write('\t'.join(row) + '\n')

    if na_genomes and logger:
        logger.warning(
            f"{len(na_genomes)}个基因组无任何>=80%的配对|Genomes with no pair >=80% "
            f"ANI: {na_genomes}")

