#!/usr/bin/env python3
"""step9: 按 bed3 插入点切分被穿透的 feature
|Split features penetrated by bed3 insertion points

R 原版 9gff_split_by_bed3_6.R 的 Python 重写。对非 gene feature, 若 bed3 插入点
bed3[,2] 落在其 [start, end) 内 (start <= p < end), 切成两条 [s, p] 与
[bed3[,3]+1, e], 第9列 ID 加 _1/_2 后缀。gene 行不切(R 原版 gene 分支 pv_split
是死代码, 函数未 return)。多点切同一 feature 时后缀累积。bed3[,2]/bed3[,3] 与
gff(update1.gff)同在 pan 坐标系。
|Python rewrite. Non-gene feature penetrated by bed3[,2] (start<=p<end) is
split into [s,p] and [bed3[,3]+1,e] with ID suffix _1/_2. Gene branch is
dead code in R. Multi-point cuts accumulate suffixes.

用法|Usage: python3 9gff_split.py <gff> <bed3> <out>
"""

import re
import sys
import numpy as np
import pandas as pd

_ID_RE = re.compile(r'(ID=[^;]+);')


def _add_id_suffix(attr: str, suffix: str) -> str:
    """复现 R gsub('(ID=.+?);','\\1<suffix>;'): 首个 ID=...; 加后缀
    |Mirror R gsub: append suffix to first ID=...;"""
    return _ID_RE.sub(r'\1' + suffix + ';', attr)


def gff_split_by_bed3(gff_path: str, bed3_path: str, out_path: str):
    # 空 gff 守卫:注释/空文件会触发 EmptyDataError,写空文件返回(对齐 R 原版空输出)
    # |Empty-gff guard: comment-only/empty file raises EmptyDataError; write empty & return
    try:
        gff = pd.read_csv(gff_path, sep='\t', header=None, dtype=str, comment='#')
    except pd.errors.EmptyDataError:
        open(out_path, 'w').close()
        return
    if gff.empty:
        open(out_path, 'w').close()
        return
    bed3 = pd.read_csv(bed3_path, sep='\t', header=None, dtype=str, comment='#')
    g_chr = gff[0].values
    b_chr = bed3[0].values
    b_start = bed3[1].values.astype(np.int64)   # R bed3[,2] 切点
    b_end = bed3[2].values.astype(np.int64)     # R bed3[,3]

    out_dfs = []
    bset = set(pd.unique(b_chr))
    for chr_ in pd.unique(b_chr):   # R parLapply unique(bed3[,1])
        gm = np.where(g_chr == chr_)[0]
        bm = np.where(b_chr == chr_)[0]
        if gm.size == 0:            # 守卫: gff 0 行 → 空(对齐 R return 空 df)
            continue
        if bm.size == 0:            # 守卫: bed3 0 行 → gff 原样
            out_dfs.append(gff.iloc[gm])
            continue
        rows = gff.iloc[gm].values.tolist()
        starts = np.array([int(r[3]) for r in rows])
        ends = np.array([int(r[4]) for r in rows])
        types = np.array([r[2] for r in rows])
        for j in bm:
            p = int(b_start[j])
            pe = int(b_end[j])
            mask = (starts <= p) & (p < ends) & (types != 'gene')
            hit = np.where(mask)[0]
            if hit.size == 0:
                continue
            # 从后往前替换: 前面索引不变; a/b 坐标已不含 p 在内部, 同点不再切
            for i in reversed(hit):
                orig = rows[i]
                a = orig.copy()
                a[4] = str(p)
                a[8] = _add_id_suffix(orig[8], '_1')
                b = orig.copy()
                b[3] = str(pe + 1)
                b[8] = _add_id_suffix(orig[8], '_2')
                rows[i:i + 1] = [a, b]
            starts = np.array([int(r[3]) for r in rows])
            ends = np.array([int(r[4]) for r in rows])
            types = np.array([r[2] for r in rows])
        out_dfs.append(pd.DataFrame(rows, columns=gff.columns))

    # chr_no_in_bed3 的 gff 行(R setdiff(chr_in_gff, chr_in_bed3), 原样保留)
    for c in [c for c in pd.unique(g_chr) if c not in bset]:
        out_dfs.append(gff.iloc[np.where(g_chr == c)[0]])

    # 空 gff 守卫:无任何输出行时写空文件,避免 pd.concat([]) 报错(对齐 R 原版空输出)
    # |Empty-gff guard: write empty file when nothing to emit, avoid pd.concat([]) error
    if not out_dfs:
        open(out_path, 'w').close()
        return

    pd.concat(out_dfs, ignore_index=True).to_csv(
        out_path, sep='\t', header=False, index=False)


def main():
    if len(sys.argv) != 4:
        sys.exit("用法|Usage: python3 9gff_split.py <gff> <bed3> <out>")
    gff_split_by_bed3(sys.argv[1], sys.argv[2], sys.argv[3])
    print(f"9gff_split → {sys.argv[3]}", file=sys.stderr)


if __name__ == '__main__':
    main()
