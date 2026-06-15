"""
统计汇总报告模块|Summary Report Module

GFFcompare 统计、校正决策汇总
"""

import os
import logging
from typing import Dict, List, Optional, TYPE_CHECKING

from .utils import CommandRunner, build_conda_command

if TYPE_CHECKING:
    from .config import RnaseqValConfig


class SummaryReporter:
    """统计汇总报告生成器|Summary Report Generator"""

    def __init__(self, config: 'RnaseqValConfig', logger: logging.Logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: logger 实例|Logger instance
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def generate_summary(
        self,
        gffcmp_prefix: Optional[str] = None,
        correction_dir: Optional[str] = None,
    ) -> Optional[str]:
        """生成汇总报告|Generate summary report

        Args:
            gffcmp_prefix: GFFcompare 输出前缀|GFFcompare output prefix
            correction_dir: 校正结果目录|Correction results directory

        Returns:
            Optional[str]: 输出 TSV 文件路径，失败返回 None|Output TSV path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "07_correction")
        os.makedirs(out_dir, exist_ok=True)
        summary_file = os.path.join(out_dir, "summary.tsv")

        self.logger.step("步骤: 生成统计汇总|Step: Generating summary report")

        lines = []

        # 1. GFFcompare 统计|GFFcompare statistics
        if gffcmp_prefix:
            stats_file = f"{gffcmp_prefix}.stats"
            if os.path.exists(stats_file):
                lines.append("# GFFcompare 统计|GFFcompare Statistics")
                lines.append(f"{'指标|Metric':<40}\t{'值|Value'}")
                lines.append("-" * 60)
                self._parse_gffcmp_stats(stats_file, lines)
                lines.append("")

        # 2. 校正决策统计|Correction decision statistics
        if correction_dir and os.path.isdir(correction_dir):
            lines.append("# 校正决策统计|Correction Decision Statistics")
            lines.append(f"{'类别|Category':<40}\t{'数量|Count'}\t{'占比|Percentage'}")
            lines.append("-" * 60)
            self._count_correction_categories(correction_dir, lines)
            lines.append("")

        # 3. 数据来源统计|Data source statistics
        lines.append("# 数据来源统计|Data Source Statistics")
        sr_count = len(self.config.samples_sr) if self.config.samples_sr else 0
        lr_count = len(self.config.samples_lr) if self.config.samples_lr else 0
        lines.append(f"{'二代样本数|SR sample count':<40}\t{sr_count}")
        lines.append(f"{'三代样本数|LR sample count':<40}\t{lr_count}")
        lines.append(f"{'三代平台|LR platform':<40}\t{self.config.lr_platform}")
        lines.append(f"{'链特异性|Strandness':<40}\t{self.config.strandness}")
        lines.append("")

        # 写入文件|Write to file
        with open(summary_file, "w") as f:
            f.write("\n".join(lines) + "\n")

        self.logger.info(f"汇总报告已写入|Summary report written: {summary_file}")
        return summary_file

    def _parse_gffcmp_stats(self, stats_file: str, lines: List[str]):
        """解析 GFFcompare .stats 文件|Parse GFFcompare .stats file

        Args:
            stats_file: stats 文件路径|stats file path
            lines: 输出行列表|Output lines list
        """
        # 提取关键指标：Sensitivity, Precision
        interesting_metrics = [
            "Sensitive", "Precision",
            "Novel loci", "Novel isoforms",
        ]

        with open(stats_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                for metric in interesting_metrics:
                    if metric in line:
                        parts = line.split()
                        if len(parts) >= 2:
                            lines.append(f"{parts[0]:<40}\t{parts[1]}")
                        break

    def _count_correction_categories(self, correction_dir: str, lines: List[str]):
        """统计各类别转录本数量|Count transcripts in each correction category

        Args:
            correction_dir: 校正结果目录|Correction directory
            lines: 输出行列表|Output lines list
        """
        import glob

        total = 0
        category_files = {
            "KEEP_HIGH_CONFIDENCE": "keep_high_confidence.tsv",
            "KEEP_LR_SUPPORTED": "keep_lr_supported.tsv",
            "KEEP_SR_SUPPORTED": "keep_sr_supported.tsv",
            "KEEP_WITH_FLAG": "keep_with_flag.tsv",
            "MANUAL_REVIEW": "manual_review.tsv",
        }

        counts = {}
        for cat, filename in category_files.items():
            filepath = os.path.join(correction_dir, filename)
            if os.path.exists(filepath):
                with open(filepath, "r") as f:
                    line_count = sum(1 for l in f if l.strip() and not l.startswith("#"))
                # 减去表头行|Subtract header line
                count = max(line_count - 1, 0)
            else:
                count = 0
            counts[cat] = count
            total += count

        for cat, count in counts.items():
            pct = (count / total * 100) if total > 0 else 0
            lines.append(f"{cat:<40}\t{count}\t{pct:.1f}%")

        lines.append(f"{'总计|Total':<40}\t{total}\t100.0%")

    def run_multiqc(self, qc_dirs: Optional[List[str]] = None) -> bool:
        """运行 MultiQC 汇总 QC 报告|Run MultiQC to aggregate QC reports

        Args:
            qc_dirs: QC 报告目录列表|QC report directory list

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        if not qc_dirs:
            self.logger.info("无 QC 目录，跳过 MultiQC|No QC directories, skipping MultiQC")
            return True

        out_dir = os.path.join(self.config.output_dir, "07_correction")
        os.makedirs(out_dir, exist_ok=True)

        dirs_str = " ".join(qc_dirs)
        cmd = f"multiqc {dirs_str} -o {out_dir}"
        cmd = build_conda_command(cmd, self.config.conda_env)

        success = self.cmd_runner.run(cmd, "MultiQC 汇总|MultiQC aggregation")
        return success
