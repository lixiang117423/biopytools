"""
归属判定模块|Assignment Module

解析 PAF 文件，按染色体聚合匹配碱基数，判定每条染色体的亚基因组归属
|Parse PAF files, aggregate per-chromosome matching bases, and assign each
chromosome to its most likely subgenome
"""

import collections
from pathlib import Path
from typing import Dict, List

from .utils import format_number


class SubgenomeAssigner:
    """亚基因组归属判定器|Subgenome Assigner"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    @staticmethod
    def _parse_paf(path: str) -> Dict[str, Dict[str, int]]:
        """
        解析单个 PAF 文件，按 query 聚合匹配碱基数
        |Parse one PAF, aggregate matching bases per query

        PAF 列定义（前 12 列必须）|PAF columns (first 12 required):
            0 qname  1 qlen  2 qstart  3 qend
            4 +/-    5 tname 6 tlen    7 tstart 8 tend
            9 nmatch 10 blocklen 11 mapq

        Returns:
            {qname: {'match': int, 'block': int}}
        """
        agg = collections.defaultdict(lambda: {"match": 0, "block": 0})
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                qname = f[0]
                nmatch = int(f[9])
                blocklen = int(f[10])
                agg[qname]["match"] += nmatch
                agg[qname]["block"] += blocklen
        return agg

    def assign(
        self,
        target_fai: Path,
        paf_paths: Dict[str, str],
    ) -> List[dict]:
        """
        对目标基因组每条染色体判定亚基因组归属
        |Assign each chromosome in target to a subgenome

        Args:
            target_fai: 目标基因组 .fai 路径|Path to target genome .fai
            paf_paths: {parent_name: paf_path}|dict of parent_name to PAF path

        Returns:
            归属结果列表（每个 dict 一行）|list of result dicts
        """
        self.logger.info("解析 PAF 并判定归属|Parsing PAFs and assigning subgenomes")

        # 每个 PAF 解析一次|Parse each PAF once
        per_parent_stats: Dict[str, Dict[str, Dict[str, int]]] = {}
        for parent_name, paf in paf_paths.items():
            self.logger.info(f"解析|Parse: {paf}")
            per_parent_stats[parent_name] = self._parse_paf(paf)

        # 从 .fai 读取所有染色体（防止漏掉未比对上的）|Read all chroms from .fai
        all_chroms = []
        with open(target_fai) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts:
                    all_chroms.append(parts[0])

        # 每条染色体计算归属|Compute assignment per chromosome
        parent_names = list(paf_paths.keys())
        rows = []
        for chrom in all_chroms:
            row = {'chrom': chrom}
            scores = {}
            for parent_name in parent_names:
                stats = per_parent_stats[parent_name].get(chrom, {"match": 0, "block": 0})
                row[f"{parent_name}_match"] = stats["match"]
                row[f"{parent_name}_block"] = stats["block"]
                scores[parent_name] = stats["match"]

            # 找最高分亲本|Find best-scoring parent
            total = sum(scores.values())
            if total == 0:
                row['assigned'] = 'UNASSIGNED'
                row['confidence'] = 0.0
                row['best_score'] = 0
                row['total_score'] = 0
            else:
                best_parent = max(scores, key=scores.get)
                best_score = scores[best_parent]
                row['assigned'] = best_parent
                row['confidence'] = round(best_score / total, 4)
                row['best_score'] = best_score
                row['total_score'] = total

            # 低置信度标记|Mark low confidence
            row['status'] = (
                'LOW_CONFIDENCE' if (
                    row['assigned'] != 'UNASSIGNED'
                    and row['confidence'] < self.config.min_conf
                ) else 'OK'
            )

            rows.append(row)

        return rows

    def write_tsv(self, rows: List[dict], out_file: Path) -> bool:
        """写出归属 TSV|Write assignment TSV"""
        if not rows:
            self.logger.error("无归属数据|No assignment data")
            return False

        # 动态列：chrom + <parent>_match + <parent>_block + best/total/conf/status
        parent_cols = []
        for k in rows[0]:
            if k.endswith('_match') or k.endswith('_block'):
                parent_cols.append(k)

        header = (
            ['chrom']
            + parent_cols
            + ['best_score', 'total_score', 'assigned', 'confidence', 'status']
        )

        with open(out_file, 'w') as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                line = [str(r.get(col, '')) for col in header]
                fh.write("\t".join(line) + "\n")

        self.logger.info(f"归属表已写|Assignment TSV written: {out_file}")
        return True

    def summarize(self, rows: List[dict]) -> None:
        """打印归属汇总|Print assignment summary"""
        if not rows:
            return

        # 按亲本统计|Count per parent
        per_parent = collections.Counter(r['assigned'] for r in rows)
        low_conf = sum(1 for r in rows if r['status'] == 'LOW_CONFIDENCE')
        unassigned = per_parent.pop('UNASSIGNED', 0)

        self.logger.info("=" * 60)
        self.logger.info(f"归属汇总|Assignment summary: {len(rows)} 条染色体|chromosomes")
        self.logger.info("-" * 60)
        for parent, count in sorted(per_parent.items()):
            self.logger.info(
                f"  {parent}: {count} 条|chroms "
                f"({format_number(sum(r['best_score'] for r in rows if r['assigned']==parent))} bp match)"
            )
        if unassigned:
            self.logger.warning(f"  UNASSIGNED: {unassigned} 条|chroms (无任何亲本比对)")
        if low_conf:
            self.logger.warning(
                f"  LOW_CONFIDENCE: {low_conf} 条|chroms "
                f"(confidence < {self.config.min_conf}，需手工复核)"
            )
            # 列出低置信度染色体|List low-confidence chroms
            for r in rows:
                if r['status'] == 'LOW_CONFIDENCE':
                    self.logger.warning(
                        f"    {r['chrom']}: confidence={r['confidence']}, "
                        f"best={r['assigned']}({r['best_score']}), "
                        f"total={r['total_score']}"
                    )
        self.logger.info("=" * 60)
