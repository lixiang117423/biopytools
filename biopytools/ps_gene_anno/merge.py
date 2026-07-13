"""
ps-gene-anno 合并输出(GFF3)|Merge output (GFF3)
braker.gff3(移除错误合并基因+其mRNA+子) + gap GFF3 行 → merged.gff3
|braker.gff3 (minus merged genes) + gap GFF3 lines → merged.gff3
"""

import os
from typing import List, Set

from .gap_analysis import _parse_gff3_attr


def merge_gff3(braker_gff3: str, gap_lines: List[str],
               merged_gene_ids: Set[str], out_path: str) -> str:
    """
    合并 braker.gff3 与 gap 基因 → merged.gff3|Merge braker.gff3 + gap genes

    移除 merged_gene_ids 中的基因(gene + 其 mRNA + 其 exon/CDS 等子, 按 Parent 链)
    |Remove merged genes and their children (mRNA/exon/CDS) via Parent chain
    其余 braker 行原样保留, 追加 gap_lines。
    """
    # 第一遍: 找 merged gene 的所有 mRNA IDs(用于移除其子)|find merged mRNAs
    merged_mrna_ids: Set[str] = set()
    if os.path.exists(braker_gff3):
        with open(braker_gff3) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                feat, attr = cols[2], cols[8]
                attrs = _parse_gff3_attr(attr)
                if feat in ('mRNA', 'transcript') and attrs.get('Parent') in merged_gene_ids:
                    mid = attrs.get('ID', '')
                    if mid:
                        merged_mrna_ids.add(mid)

    # 第二遍: 移除 merged gene/mRNA/子, 保留其余|filter
    kept_lines = []
    if os.path.exists(braker_gff3):
        with open(braker_gff3) as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    kept_lines.append(line.rstrip('\n'))
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    kept_lines.append(line.rstrip('\n'))
                    continue
                feat, attr = cols[2], cols[8]
                attrs = _parse_gff3_attr(attr)
                gid = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                # 移除: merged gene 本身 / 其 mRNA / 其 mRNA 的子
                if feat == 'gene' and gid in merged_gene_ids:
                    continue
                if feat in ('mRNA', 'transcript') and parent in merged_gene_ids:
                    continue
                if parent in merged_mrna_ids:
                    # exon/CDS/start_codon/stop_codon/intron 等 mRNA 的子
                    continue
                kept_lines.append(line.rstrip('\n'))

    with open(out_path, 'w') as out:
        out.write("##gff-version 3\n")
        out.write("\n".join(kept_lines))
        if kept_lines:
            out.write("\n")
        if gap_lines:
            out.write("\n".join(gap_lines))
            out.write("\n")
    return out_path
