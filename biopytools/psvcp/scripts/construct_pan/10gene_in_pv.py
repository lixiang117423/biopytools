#!/usr/bin/env python3
"""step10: 找完全落在 PAV 插入区间内的 query feature
|Find query features fully inside PAV insertion intervals

R 原版 10gene_in_pv_from_gff_parLapply2.R 的 Python 向量化重写。
在 **query 原坐标系** 下判断: gff feature [start,end] 是否完全落在 bed3 的
query 插入区间 [bed3.col8, bed3.col8 + len(bed3.col5) - 1] 内; 命中则 cbind
(gff 行 + bed3 行) 输出。无命中不写文件(对齐 R: if nrow>0 才 write)。
|Python vectorized rewrite. Tests in **query coord** whether gff feature
fully inside bed3 query insertion interval; cbind matched rows. No file
on 0 hits (mirrors R).

bed3 列结构(8 列)|bed3 columns (8):
  [1]pan_chr [2]pan_start [3]pan_end [4]type [5]序列 [6]source [7]原chr [8]原起点

用法|Usage: python3 10gene_in_pv.py <query_gff> <bed3> <out>
"""

import sys
import numpy as np
import pandas as pd


def gene_in_pv(gff_path: str, bed3_path: str, out_path: str) -> int:
    """返回命中行数(0 则不写文件)|Return hit count; write file only if >0"""
    gff = pd.read_csv(gff_path, sep='\t', header=None, dtype=str, comment='#')
    bed3 = pd.read_csv(bed3_path, sep='\t', header=None, dtype=str, comment='#')

    g_chr = gff[0].values
    g_start = gff[3].values.astype(np.int64)
    g_end = gff[4].values.astype(np.int64)
    # R 语义: chr 循环用 bed3[,1](pan 染色体名集合); 选 bed3 行用 bed3[,7](query_chr)==chr;
    # gff 用 [,1]==chr → 实际匹配 gff_chr==bed3 query_chr(均在 query 坐标系)
    # |R semantics: loop chr over bed3[,1]; pick bed3 by bed3[,7](query_chr)==chr;
    # gff by [,1]==chr → effectively gff_chr==bed3 query_chr (both in query coord)
    b_pan = bed3[0].values                              # R bed3[,1] chr 循环来源
    b_qchr = bed3[6].values                             # R bed3[,7] query_chr
    b_qstart = bed3[7].values.astype(np.int64)          # R bed3[,8] query 起点
    b_len = bed3[4].map(len).values.astype(np.int64)    # R nchar(bed3[,5])
    b_qend = b_qstart + b_len - 1

    out_rows = []
    for chr_ in pd.unique(b_pan):   # 对齐 R unique(bed3[,1])
        gm = np.where(g_chr == chr_)[0]
        bm = np.where(b_qchr == chr_)[0]   # 对齐 R bed3[bed3[,7]==chr]
        if gm.size == 0 or bm.size == 0:
            continue
        # bed3 按 query 起点排序(插入位点不重叠 → qend 亦递增)
        order = np.argsort(b_qstart[bm], kind='stable')
        bs = b_qstart[bm][order]
        be = b_qend[bm][order]
        bm_ord = bm[order]
        gs = g_start[gm]
        ge = g_end[gm]
        # 每 gff 找 b_qstart <= g_start 的最近一个 bed3, 验证 qend >= g_end
        k = np.searchsorted(bs, gs, side='right') - 1
        kc = np.clip(k, 0, bs.size - 1)
        hit = (k >= 0) & (bs[kc] <= gs) & (be[kc] >= ge)
        if hit.any():
            g_hit = gm[hit]
            b_hit = bm_ord[kc[hit]]
            g_arr = gff.values
            b_arr = bed3.values
            for gi, bi in zip(g_hit, b_hit):
                out_rows.append('\t'.join(g_arr[gi].tolist() + b_arr[bi].tolist()))

    if out_rows:
        with open(out_path, 'w', encoding='utf-8') as out:
            out.write('\n'.join(out_rows) + '\n')
    return len(out_rows)


def main():
    if len(sys.argv) != 4:
        sys.exit("用法|Usage: python3 10gene_in_pv.py <query_gff> <bed3> <out>")
    n = gene_in_pv(sys.argv[1], sys.argv[2], sys.argv[3])
    print(f"10gene_in_pv: {n} 行命中|rows in PAV → {sys.argv[3]}", file=sys.stderr)


if __name__ == '__main__':
    main()
