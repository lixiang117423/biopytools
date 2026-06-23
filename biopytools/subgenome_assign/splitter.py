"""
FASTA 拆分模块|FASTA Splitter

按归属表把目标基因组拆成每个亲本一个 FASTA（+ 可选未归属 FASTA）
|Split target genome by assignment into one FASTA per parent (+ optional unassigned)
"""

import os
from pathlib import Path
from typing import List


class FastaSplitter:
    """FASTA 拆分器|FASTA Splitter"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def split(
        self,
        rows: List[dict],
        target_fasta: str,
        out_dir: Path,
    ) -> dict:
        """
        按归属拆分 FASTA|Split FASTA by assignment

        Args:
            rows: 归属结果列表|list of assignment dicts
            target_fasta: 目标基因组 FASTA|target genome FASTA
            out_dir: 输出目录|output directory

        Returns:
            {parent_name: out_path} 字典|dict of parent_name to output FASTA path
        """
        if not self.config.split_fasta:
            self.logger.info("未启用 FASTA 拆分|FASTA splitting disabled")
            return {}

        # 按归属分组染色体|Group chroms by assignment
        from collections import defaultdict
        groups = defaultdict(list)
        for r in rows:
            groups[r['assigned']].append(r['chrom'])

        # 确保 samtools faidx 索引存在|Ensure .fai exists
        fai = Path(target_fasta).with_suffix(target_fasta.endswith('.gz') and '.fai' or Path(target_fasta).suffix + '.fai')
        if not fai.exists():
            self.logger.info(f"构建 .fai 索引|Building .fai: {fai.name}")
            args = ['faidx', target_fasta]
            cmd = [self.config.samtools_path] + args
            if not self.cmd_runner.run(cmd, f"samtools faidx|samtools faidx"):
                return {}

        # 每个亲本一份 FASTA（包括 UNASSIGNED 如果开启）|One FASTA per parent
        results = {}
        for parent_name, chroms in groups.items():
            if parent_name == 'UNASSIGNED' and not self.config.keep_unassigned:
                self.logger.info(f"跳过 UNASSIGNED（共 {len(chroms)} 条）|Skipping UNASSIGNED")
                continue

            out_file = out_dir / f"subgenome_{parent_name}.fa"
            self._extract(target_fasta, chroms, out_file)
            results[parent_name] = str(out_file)
            self.logger.info(
                f"  {parent_name}: {len(chroms)} 条染色体|chroms -> {out_file.name}"
            )

        return results

    def _extract(self, target_fasta: str, chroms: list, out_file: Path) -> bool:
        """samtools faidx 提取多个染色体|Extract multiple chroms via samtools faidx"""
        # 列表文件（避免命令行过长）|List file (avoid long cmdline)
        list_file = out_file.with_suffix('.list')
        with open(list_file, 'w') as fh:
            for c in chroms:
                fh.write(c + "\n")

        args = ['faidx', target_fasta]
        # samtools faidx 支持 -r <region_file> 批量提取
        # |samtools faidx supports -r <region_file> for batch extraction
        args.extend(['-r', str(list_file)])

        cmd = [self.config.samtools_path] + args

        # 输出重定向，用 run_pipeline|Redirect output, use run_pipeline
        pipeline = f"{' '.join(cmd)} > {out_file}"
        return self.cmd_runner.run_pipeline(
            pipeline,
            f"samtools faidx|samtools faidx: {out_file.name}"
        )
