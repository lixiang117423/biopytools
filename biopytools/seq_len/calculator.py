"""seq_len 核心计算逻辑|seq_len core calculator (orchestration, no CLI)"""

import os
from typing import List, Tuple

from .utils import compute_summary, format_number, iter_seq_lengths

# 单条记录:(来源文件 basename, 序列 id, 长度)|record: (source basename, seq_id, length)
Record = Tuple[str, str, int]

_SUMMARY_HEADER = ("source_file", "num_seqs", "total_length", "n50",
                    "min_length", "max_length", "mean_length")


class SeqLenCalculator:
    """收集长度 → 过滤/排序 → 写主表 + 汇总|Collect lengths, filter/sort, write tables"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run(self) -> bool:
        """执行全流程|Run the full pipeline"""
        cfg = self.config

        # 1. 流式收集记录(不载入整条序列)|stream-collect records
        records: List[Record] = []
        for path in cfg.input_files:
            source = os.path.basename(path)
            for seq_id, length in iter_seq_lengths(path):
                records.append((source, seq_id, length))

        # 2. 最小长度过滤|min-length filter
        if cfg.min_length > 0:
            records = [r for r in records if r[2] >= cfg.min_length]

        # 3. 按长度降序(可选)|optional length-descending sort
        if cfg.sort:
            records.sort(key=lambda r: r[2], reverse=True)

        # 4. 写主表|write main table
        self._write_main(records)
        self.logger.info(
            f"写入主表|Wrote main table: {cfg.tsv_path} "
            f"({format_number(len(records))} 条序列|sequences)")

        # 5. 写汇总(可选)|write summary (optional)
        if cfg.summary:
            self._write_summary(records)

        return True

    def _write_main(self, records: List[Record]) -> None:
        """写主表 TSV(文件夹模式带 source_file 列)|Write main TSV (source_file col in folder mode)"""
        cfg = self.config
        with open(cfg.tsv_path, 'w') as f:
            if cfg.is_folder:
                f.write("source_file\tsequence_id\tlength\n")
                for source, seq_id, length in records:
                    f.write(f"{source}\t{seq_id}\t{length}\n")
            else:
                f.write("sequence_id\tlength\n")
                for _source, seq_id, length in records:
                    f.write(f"{seq_id}\t{length}\n")

    def _write_summary(self, records: List[Record]) -> None:
        """写汇总 TSV:每文件一行 + 文件夹模式追加 ALL 行|Write summary TSV (per-file + ALL in folder mode)"""
        cfg = self.config
        # 按来源分组(顺序无关,统计是基于全量)|group by source
        by_source = {}
        for source, _seq_id, length in records:
            by_source.setdefault(source, []).append(length)

        with open(cfg.summary_path, 'w') as f:
            f.write('\t'.join(_SUMMARY_HEADER) + '\n')
            # 每文件一行,按输入文件顺序|one row per input file, in input order
            for path in cfg.input_files:
                source = os.path.basename(path)
                self._write_summary_row(f, source, by_source.get(source, []))
            # 文件夹模式(>1 文件)追加全局 ALL 行|folder mode (>1 file): append ALL row
            if len(cfg.input_files) > 1:
                all_lengths = [length for _s, _i, length in records]
                self._write_summary_row(f, 'ALL', all_lengths)

        self.logger.info(f"写入汇总|Wrote summary: {cfg.summary_path}")

    def _write_summary_row(self, f, source: str, lengths: List[int]) -> None:
        """写一行汇总统计|Write one summary row"""
        s = compute_summary(lengths)
        f.write(f"{source}\t{s['num_seqs']}\t{s['total_length']}\t{s['n50']}\t"
                f"{s['min_length']}\t{s['max_length']}\t{s['mean_length']}\n")
        total = format_number(s['total_length'])
        self.logger.info(
            f"{source}: {s['num_seqs']} 条|seqs, 总长|total {total}, "
            f"N50 {format_number(s['n50'])}, "
            f"最短|min {s['min_length']}, 最长|max {s['max_length']}")
