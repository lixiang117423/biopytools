"""check_reads 流程编排器|check_reads orchestrator

并行 gz 完整性校验 → 0 字节/配对完整性 → 报告与退出码
|Parallel gzip integrity check → empty/pairing checks → report & exit code
"""

import gzip
import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

from .utils import collect_fastq_files, detect_read_pair


class CheckReadsCalculator:
    """fastq 完整性检查器|FASTQ integrity checker"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    # ---------- 单文件检查|per-file checks ----------

    @staticmethod
    def check_gzip(path: str) -> bool:
        """流式校验 gz 完整性|Stream-verify gzip integrity"""
        try:
            with gzip.open(path, "rb") as fh:
                while fh.read(1 << 20):
                    pass
            return True
        except (OSError, EOFError):
            return False

    def check_file(self, path: str) -> str:
        """单个文件状态|Per-file status

        Returns:
            状态|status: OK / CORRUPT / EMPTY / PLAIN(明文,跳过压缩校验)
        """
        if os.path.getsize(path) == 0:
            return "EMPTY"
        if not path.lower().endswith(".gz"):
            return "PLAIN"
        return "OK" if self.check_gzip(path) else "CORRUPT"

    # ---------- 配对完整性|pairing check ----------

    @staticmethod
    def check_pairing(fastq_paths) -> list:
        """检测缺配对的样本|Detect samples missing a mate

        Returns:
            [(样本名|sample, 缺|missing 'R1'|'R2'|'mate')]
        """
        issues = []
        by_sample = {}
        for p in fastq_paths:
            parsed = detect_read_pair(os.path.basename(p))
            if parsed is None:
                continue  # 单端|single-end
            sample, direction = parsed
            by_sample.setdefault(sample, set()).add(direction)
        for sample in sorted(by_sample):
            dirs = by_sample[sample]
            if len(dirs) == 1:
                missing = "R2" if "1" in dirs else "R1"
                issues.append((sample, missing))
        return issues

    # ---------- 报告|report ----------

    def _write_report(self, statuses, pairing_issues, singles,
                      ignored, start_time) -> None:
        """写报告文件|Write report file"""
        lines = [
            f"=== check_reads 报告|report | {start_time:%Y-%m-%d %H:%M:%S} ===",
            f"输入|Input: {self.config.input_dir} | 线程|threads: {self.config.threads}",
            "-- 文件状态|file status --",
        ]
        counts = {}
        for path in sorted(statuses):
            st = statuses[path]
            counts[st] = counts.get(st, 0) + 1
            lines.append(f"{st:7s} {path}")
        lines.append("-- 配对检查|pairing --")
        if pairing_issues:
            for sample, missing in pairing_issues:
                lines.append(f"缺{missing} missing_{missing}: {sample}")
        else:
            lines.append("(无|none)")
        lines.append("-- 汇总|summary --")
        total = len(statuses)
        lines.append(f"文件总数|files: {total}")
        for st in ("OK", "CORRUPT", "EMPTY", "PLAIN"):
            lines.append(f"{st}: {counts.get(st, 0)}")
        lines.append(f"配对样本|paired samples: "
                     f"{len({os.path.basename(p) for p in statuses}) - len(singles)}")
        lines.append(f"单端|single-end: {len(singles)}")
        lines.append(f"忽略非fastq|ignored: {len(ignored)}")
        if ignored[:5]:
            lines.append(f"  示例|e.g.: {', '.join(ignored[:5])}")
        problem = counts.get("CORRUPT", 0) + counts.get("EMPTY", 0) + \
            len(pairing_issues)
        lines.append(f"结果|RESULT: "
                     f"{'发现问题(见上方)|issues found' if problem else '全部通过|all clean'}")
        with open(self.config.report_file, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        self.logger.info(f"报告已写入|Report written: {self.config.report_file}")

    # ---------- 主流程|main run ----------

    def run(self) -> int:
        """端到端检查|Run end-to-end. Returns exit code 0/1."""
        start = datetime.now()
        self.logger.info("开始 reads 完整性检查|Starting reads integrity check")
        self.logger.info(f"输入|Input: {self.config.input_dir} | "
                         f"线程|threads: {self.config.threads}")

        fastq_paths, ignored = collect_fastq_files(
            self.config.input_dirs, self.logger)
        if not fastq_paths:
            self.logger.error("未找到 fastq 文件|No fastq files found")
            return 1
        self.logger.info(f"共找到|Total fastq files: {len(fastq_paths)}")

        # 并行完整性校验|parallel integrity check
        statuses = {}
        with ThreadPoolExecutor(max_workers=self.config.threads) as pool:
            for path, status in zip(fastq_paths, pool.map(self.check_file, fastq_paths)):
                statuses[path] = status
        n_ok = sum(1 for s in statuses.values() if s == "OK")
        n_corrupt = sum(1 for s in statuses.values() if s == "CORRUPT")
        n_empty = sum(1 for s in statuses.values() if s == "EMPTY")
        self.logger.info(f"完整性|Integrity: OK={n_ok} CORRUPT={n_corrupt} "
                         f"EMPTY={n_empty}")

        # 配对检查|pairing
        pairing_issues = self.check_pairing(fastq_paths)
        singles = [p for p in fastq_paths
                   if detect_read_pair(os.path.basename(p)) is None]
        if pairing_issues:
            self.logger.warning(
                f"缺配对样本|Samples missing a mate: "
                f"{', '.join(f'{s}->{m}' for s, m in pairing_issues)}")
        self.logger.info(f"配对样本|Paired: "
                         f"{len(fastq_paths) - len(singles) - 2 * len(pairing_issues)} 对|pairs, "
                         f"单端|single-end: {len(singles)}")

        self._write_report(statuses, pairing_issues, singles, ignored, start)

        problems = n_corrupt + n_empty + len(pairing_issues)
        if problems:
            self.logger.error(
                f"发现问题|Issues found: CORRUPT={n_corrupt} EMPTY={n_empty} "
                f"缺配对|missing mate={len(pairing_issues)}")
            return 1
        self.logger.info("全部通过|All clean")
        return 0
