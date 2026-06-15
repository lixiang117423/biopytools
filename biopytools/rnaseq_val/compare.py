"""
转录本结构比较模块|Transcript Structure Comparison Module

GFFcompare 多来源 GTF 与参考注释的比较
"""

import os
import logging
from typing import List, Optional, TYPE_CHECKING

from .utils import CommandRunner, FileValidator, build_conda_command

if TYPE_CHECKING:
    from .config import RnaseqValConfig


class GFFCompareRunner:
    """GFFcompare 运行器|GFFcompare Runner"""

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
        self.file_validator = FileValidator(logger)

    def run(
        self,
        query_gtfs: List[str],
        label: str = "merged",
    ) -> Optional[str]:
        """运行 GFFcompare|Run GFFcompare

        Args:
            query_gtfs: query GTF 文件列表（二代合并 GTF + 三代 GTF）|Query GTF file list
            label: 输出前缀标签|Output prefix label

        Returns:
            Optional[str]: 输出前缀路径（不含扩展名），失败返回 None|Output prefix path, None on failure
        """
        if not query_gtfs:
            self.logger.error("无 query GTF 文件|No query GTF files provided")
            return None

        out_dir = os.path.join(self.config.output_dir, "06_compare")
        os.makedirs(out_dir, exist_ok=True)

        prefix = os.path.join(out_dir, label)
        annotated_gtf = f"{prefix}.annotated.gtf"

        # 断点续传|Checkpoint
        if os.path.exists(annotated_gtf) and os.path.getsize(annotated_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"GFFcompare 已完成，跳过|GFFcompare done, skipping")
                return prefix

        self.logger.step("步骤: 转录本结构比较 (GFFcompare)|Step: Transcript comparison (GFFcompare)")
        self.logger.info(f"参考注释|Reference annotation: {self.config.annotation_gtf}")
        self.logger.info(f"query GTF 数量|Query GTF count: {len(query_gtfs)}")

        # 构建 GFFcompare 命令|Build GFFcompare command
        cmd_parts = [
            "gffcompare",
            "-r", self.config.annotation_gtf,
            "-o", prefix,
        ]
        for gtf in query_gtfs:
            cmd_parts.append(gtf)

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "GFFcompare")

        if not success:
            self.logger.error("GFFcompare 运行失败|GFFcompare failed")
            return None

        # 检查关键输出文件|Check key output files
        stats_file = f"{prefix}.stats"
        tracking_file = f"{prefix}.tracking"

        if os.path.exists(stats_file):
            self.logger.info(f"GFFcompare 统计文件|GFFcompare stats: {stats_file}")
        if os.path.exists(tracking_file):
            self.logger.info(f"GFFcompare tracking 文件|GFFcompare tracking: {tracking_file}")

        self.logger.info(f"GFFcompare 完成|GFFcompare done: {annotated_gtf}")
        return prefix


class TrackingParser:
    """GFFcompare .tracking 文件解析器|GFFcompare .tracking file parser"""

    def __init__(self, logger: logging.Logger):
        """初始化|Initialize

        Args:
            logger: logger 实例|Logger instance
        """
        self.logger = logger

    def parse(self, tracking_file: str) -> List[dict]:
        """解析 .tracking 文件|Parse .tracking file

        tracking 文件格式：
        列1: query 转录本 ID
        列2: 参考转录本 ID (class_code 赋值来源)
        列3+: 各输入 GTF 中对应的转录本 ID

        Args:
            tracking_file: .tracking 文件路径|.tracking file path

        Returns:
            List[dict]: 每行解析为字典，包含 transcript_id, class_code, ref_id, source_ids 等
                Each row parsed as a dict with transcript_id, class_code, ref_id, source_ids, etc.
        """
        results = []

        if not os.path.exists(tracking_file):
            self.logger.error(f"tracking 文件不存在|tracking file not found: {tracking_file}")
            return results

        # 读取注释 GTF 获取 class_code 信息
        # tracking 文件本身不包含 class_code 列，需要从 annotated GTF 获取
        # 这里先解析 tracking 结构，class_code 在 correct.py 中从 annotated.gtf 获取
        with open(tracking_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue

                record = {
                    "transcript_id": parts[0],
                    "ref_id": parts[1] if parts[1] != "-" else None,
                }
                results.append(record)

        self.logger.info(f"解析 tracking 文件|Parsed tracking file: {len(results)} 条记录|records")
        return results
