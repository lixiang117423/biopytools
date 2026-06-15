"""泛基因组Block构建 - 两两比对模块|Pan-Blocks Construction - Pairwise Alignment Module"""

import os
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Tuple
import logging

from .config import PanBlocksConfig
from .utils import build_conda_command, CommandRunner


class PairwiseAligner:
    """两两基因组比对执行器|Pairwise Genome Aligner"""

    def __init__(self, config: PanBlocksConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger)
        self.coords_dir = Path(config.coords_dir)

    def run_all_alignments(self):
        """运行所有两两基因组比对|Run all pairwise genome alignments"""
        genome_order = self.config.genome_order_list
        genome_map = {g['name']: g['path'] for g in self.config.genomes}

        pairs = self._generate_all_pairs(genome_order, genome_map)
        pending = self._filter_pending(pairs)
        total = len(pairs)
        total_all = len(pairs) + (len(genome_order) * (len(genome_order) - 1) - total)

        self.logger.info(f"总比对数|Total pairwise alignments: {total_all}")
        self.logger.info(f"已完成|Completed: {total_all - total}, 待执行|Pending: {total}")

        if total == 0:
            self.logger.info("所有比对已完成|All alignments completed")
            return True

        if not pending:
            return True

        self.logger.info(f"开始并行比对|Starting parallel alignments: {self.config.parallel_alignments} 并发|concurrent")
        success_count = 0
        fail_count = 0

        with ThreadPoolExecutor(max_workers=self.config.parallel_alignments) as executor:
            futures = {executor.submit(self._run_single_alignment, *pair): pair for pair in pending}
            for future in as_completed(futures):
                ref_name, query_name = futures[future]
                try:
                    success = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    self.logger.error(f"比对异常|Alignment exception: {ref_name} vs {query_name}: {e}")
                    fail_count += 1

        self.logger.info(f"比对完成|Alignments completed: 成功|success={success_count}, 失败|fail={fail_count}")
        return fail_count == 0

    def _generate_all_pairs(self, genome_order: List[str],
                            genome_map: dict) -> List[Tuple[str, str, str, str]]:
        """生成所有有向基因组对|Generate all directed genome pairs"""
        pairs = []
        for ref_name in genome_order:
            for query_name in genome_order:
                if ref_name == query_name:
                    continue
                pairs.append((ref_name, query_name, genome_map[ref_name], genome_map[query_name]))
        return pairs

    def _filter_pending(self, pairs: List[Tuple]) -> List[Tuple]:
        """过滤已完成的比对|Filter out completed alignments"""
        pending = []
        for pair in pairs:
            ref_name, query_name = pair[0], pair[1]
            coords_file = self.coords_dir / f"{ref_name}.vs.{query_name}.filtered.coords"
            if coords_file.exists() and coords_file.stat().st_size > 0:
                self.logger.info(f"跳过已完成比对|Skipping completed: {ref_name} vs {query_name}")
            else:
                pending.append(pair)
        return pending

    def _run_single_alignment(self, ref_name: str, query_name: str,
                              ref_path: str, query_path: str) -> bool:
        """运行单个两两比对|Run single pairwise alignment"""
        prefix = self.coords_dir / f"{ref_name}.vs.{query_name}"
        delta_file = f"{prefix}.delta"
        filtered_delta_file = f"{prefix}.filtered.delta"
        coords_file = f"{prefix}.filtered.coords"

        nucmer_threads = max(1, self.config.threads // self.config.parallel_alignments)

        # Step 1: nucmer
        nucmer_cmd = build_conda_command(
            self.config.nucmer_path,
            ['--mum', '-t', str(nucmer_threads), ref_path, query_path, '--prefix', str(prefix)]
        )
        success, _, stderr = self.runner.run(nucmer_cmd, description=f"nucmer: {ref_name} vs {query_name}")
        if not success:
            self.logger.error(f"nucmer失败|nucmer failed: {ref_name} vs {query_name}: {stderr}")
            return False

        # Step 2: delta-filter (输出到stdout，需重定向到文件)
        delta_filter_cmd = build_conda_command(
            self.config.delta_filter_path,
            ['-l', str(self.config.min_alignment_length), '-r', '-q', str(delta_file)]
        )
        if not self._redirect_stdout_to_file(delta_filter_cmd, filtered_delta_file):
            self.logger.error(f"delta-filter失败|delta-filter failed: {ref_name} vs {query_name}")
            self._cleanup_intermediates(delta_file)
            return False

        # Step 3: show-coords
        show_coords_cmd = build_conda_command(
            self.config.show_coords_path,
            ['-TrHcl', str(filtered_delta_file)]
        )
        if not self._redirect_stdout_to_file(show_coords_cmd, coords_file):
            self._cleanup_intermediates(delta_file, filtered_delta_file)
            return False

        # Step 4: 清理中间文件
        self._cleanup_intermediates(delta_file, filtered_delta_file)

        # 验证输出
        if not os.path.exists(coords_file) or os.path.getsize(coords_file) == 0:
            self.logger.warning(f"比对结果为空|Alignment result is empty: {ref_name} vs {query_name}")
            return True

        return True

    def _redirect_stdout_to_file(self, cmd: List[str], output_file: str) -> bool:
        """执行命令并将stdout重定向到文件|Execute command and redirect stdout to file"""
        import subprocess
        self.logger.info(f"命令|Command: {' '.join(cmd)} > {output_file}")
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, shell=False)
            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed, stderr: {result.stderr.decode('utf-8', errors='ignore')}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {e}")
            return False

    def _cleanup_intermediates(self, delta_file: str, filtered_delta_file: str):
        """清理中间文件|Clean up intermediate files"""
        for f in [delta_file, filtered_delta_file]:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass
