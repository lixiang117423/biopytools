#!/usr/bin/env python3
"""step8: 用 bed2 平移 gff 坐标(插入位点累积偏移)
|Shift gff coords by accumulated bed2 insertion offsets

R 原版 8gff_update_by_bed2info_parLapply.R 的 Python 向量化重写。
对每条 gff 的 start/end, 累加所有 bed2.start < coord 的 (bed2.end-bed2.start),
即插入位点造成的坐标偏移。bed2 须按 (chr, start) 升序(R break 依赖, 已验证)。
|Python vectorized rewrite. For each gff start/end, accumulate
(bed2.end-bed2.start) over all bed2.start < coord. bed2 must be sorted
by (chr, start) — R's break relies on it (verified).

用法|Usage: python3 8gff_update.py <gff> <bed2> <out>
"""

import sys
import numpy as np
import pandas as pd


def gff_update_by_bed2(gff_path: str, bed2_path: str, out_path: str):
    gff = pd.read_csv(gff_path, sep='\t', header=None, dtype=str, comment='#')
    bed2 = pd.read_csv(bed2_path, sep='\t', header=None, dtype=str, comment='#')
    g_chr = gff[0].values
    g_start = gff[3].values.astype(np.int64)
    g_end = gff[4].values.astype(np.int64)
    b_chr = bed2[0].values
    b_start = bed2[1].values.astype(np.int64)
    b_len = bed2[2].values.astype(np.int64) - b_start   # R bed2[,3]-bed2[,2]

    new_start = g_start.copy()
    new_end = g_end.copy()
    for chr_ in pd.unique(b_chr):   # R parLapply over chr_in_bed2
        gm = np.where(g_chr == chr_)[0]
        bm = np.where(b_chr == chr_)[0]
        if gm.size == 0 or bm.size == 0:
            continue
        bs = b_start[bm]                                   # 已升序(全局排序子集)
        cum = np.concatenate([[0], np.cumsum(b_len[bm])])  # cum[i] = 前 i 个 bed2 长度和
        # R: gff[,4] > bed2[,2]  <=>  bed2.start < gff.start (严格小于)
        idx_s = np.searchsorted(bs, g_start[gm], side='left')
        new_start[gm] = g_start[gm] + cum[idx_s]
        idx_e = np.searchsorted(bs, g_end[gm], side='left')
        new_end[gm] = g_end[gm] + cum[idx_e]

    gff[3] = new_start.astype(str)
    gff[4] = new_end.astype(str)
    # R 输出序: chr_in_bed2 顺序 + chr_no_in_bed2, 每 chr 内 gff 原序
    # |R order: chr_in_bed2 then chr_no_in_bed2, gff original order within each chr
    bset = set(pd.unique(b_chr))
    chr_no_bed2 = [c for c in pd.unique(g_chr) if c not in bset]
    order = list(pd.unique(b_chr)) + chr_no_bed2
    parts = [gff[g_chr == c] for c in order if (g_chr == c).any()]
    pd.concat(parts, ignore_index=True).to_csv(
        out_path, sep='\t', header=False, index=False)


def main():
    if len(sys.argv) != 4:
        sys.exit("用法|Usage: python3 8gff_update.py <gff> <bed2> <out>")
    gff_update_by_bed2(sys.argv[1], sys.argv[2], sys.argv[3])
    print(f"8gff_update → {sys.argv[3]}", file=sys.stderr)


if __name__ == '__main__':
    main()
