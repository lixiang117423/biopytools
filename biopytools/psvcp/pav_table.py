"""PSVCP PAV 信息表与矩阵生成|PAV info table & matrix generation

从 pan.pav.gff(Insertion PAV 清单)生成两张表:
- 信息表: 每个 PAV 的 pan 区间 + 来源样本 + 原基因组位置 + 长度
- 矩阵:   样本 × PAV 的 0/1 矩阵(来源样本语义)
|Generate two tables from pan.pav.gff (Insertion PAV list):
- info: per-PAV pan interval + source sample + orig-genome position + length
- matrix: sample × PAV 0/1 matrix (source-sample semantics)

第9列 ID 格式 ID=样本_染色体_起点(1-based); 已用真实数据验证 pan 段与
原基因组段逐碱基一致、起点偏差 0(见 pav_probe 验证)。
|Column-9 ID format ID=sample_chr_start(1-based); verified base-for-base
identical between pan segment and orig-genome segment, start offset 0.

矩阵语义: 1 = 该样本贡献此插入(source), 0 = 未标注(含 ref 及未参与比对的
query)。注意这是「来源标注」矩阵,非样本间真实共享性 —— PSVCP 原流程只做
query vs ref 比对,未做样本两两比较。
|Matrix semantics: 1 = source sample of this insertion, 0 = otherwise
(incl. ref and queries not compared). This is a source-label matrix,
NOT true inter-sample presence/absence.
"""

import re
from typing import Tuple


def _split_chr_start(id_val: str, source: str) -> Tuple[str, int]:
    """从 'S42_Chr06_54016560' + source='S42' 提取 (染色体, 起点)
    用 rpartition 分最后一个下划线,染色体名含下划线也正确
    |Extract (chr, start) via last-underscore split; works even if chr name has underscores"""
    prefix = source + '_'
    rest = id_val[len(prefix):] if id_val.startswith(prefix) else id_val
    if '_' not in rest:
        return rest, 0
    chr_part, _, start_part = rest.rpartition('_')
    try:
        return chr_part, int(start_part)
    except ValueError:
        return chr_part, 0   # 起点非数字:染色体仍取 rpartition 左部|non-numeric start: keep chr_part


def _parse_pav_line(fields):
    """解析一行 PAV gff(9 列)→ dict;字段不足/坐标非法返回 None
    |Parse one PAV gff line (9 cols) → dict; None on malformed line"""
    if len(fields) < 9:
        return None
    try:
        pan_start = int(fields[3])
        pan_end = int(fields[4])
    except ValueError:
        return None
    pan_chr, source = fields[0], fields[1]
    length = pan_end - pan_start + 1
    m = re.search(r'ID=([^;]+)', fields[8])
    id_val = m.group(1) if m else ''
    orig_chr, orig_start = _split_chr_start(id_val, source)
    return {
        'pav_id': f"{pan_chr}:{pan_start}-{pan_end}",
        'pan_chr': pan_chr, 'pan_start': pan_start, 'pan_end': pan_end,
        'length_bp': length, 'source': source,
        'orig_chr': orig_chr, 'orig_start': orig_start,
        'orig_end': orig_start + length - 1 if orig_start else 0,
    }


def iter_pav_records(pav_gff: str):
    """逐行读取 PAV gff(跳过注释/空行),yield 解析后的 dict
    |Yield parsed dict per PAV line (skip comments/blanks)"""
    with open(pav_gff, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            rec = _parse_pav_line(line.split('\t'))
            if rec:
                yield rec


_INFO_COLS = ['pav_id', 'pan_chr', 'pan_start', 'pan_end', 'length_bp',
              'source', 'orig_chr', 'orig_start', 'orig_end']


def generate_pav_info(pav_gff: str, out_path: str) -> int:
    """生成 PAV 信息表 tsv,返回写出行数|Write PAV info table, return row count"""
    n = 0
    with open(out_path, 'w', encoding='utf-8') as out:
        out.write('#' + '\t'.join(_INFO_COLS) + '\n')
        for rec in iter_pav_records(pav_gff):
            out.write('\t'.join(str(rec[c]) for c in _INFO_COLS) + '\n')
            n += 1
    return n


def _read_samples(genome_list_path: str):
    """读 genome_list,返回去 .fa 后缀的样本名列表(保持顺序)|Read sample names (strip .fa)"""
    samples = []
    with open(genome_list_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                samples.append(line[:-3] if line.endswith('.fa') else line)
    return samples


def generate_pav_matrix(pav_gff: str, genome_list_path: str, out_path: str) -> int:
    """生成样本×PAV 0/1 矩阵(来源样本=1,其余=0)|Write sample×PAV 0/1 matrix"""
    samples = _read_samples(genome_list_path)
    n = 0
    with open(out_path, 'w', encoding='utf-8') as out:
        out.write('# PAV 矩阵(来源样本语义): 1=该样本贡献此插入, 0=其余(含 ref 及未比对 query)\n')
        out.write('# 注:仅来源样本标1,非样本间真实共享性(原流程只做 query vs ref)\n')
        out.write('\t'.join(['pav_id'] + samples) + '\n')
        for rec in iter_pav_records(pav_gff):
            row = [rec['pav_id']] + ['1' if s == rec['source'] else '0' for s in samples]
            out.write('\t'.join(row) + '\n')
            n += 1
    return n


def main():
    """CLI: python pav_table.py <pav_gff> -l <genome_list> --info <out> --matrix <out>"""
    import argparse
    p = argparse.ArgumentParser(description='PSVCP PAV 信息表/矩阵生成|PAV info/matrix generator')
    p.add_argument('pav_gff', help='pan.pav.gff 或 pan.pav.sorted.gff')
    p.add_argument('-l', '--genome-list', help='genome_list 文本(矩阵需要)')
    p.add_argument('--info', help='输出信息表 tsv 路径')
    p.add_argument('--matrix', help='输出矩阵 tsv 路径')
    args = p.parse_args()
    if args.info:
        n = generate_pav_info(args.pav_gff, args.info)
        print(f"信息表|info: {n} 行|rows → {args.info}")
    if args.matrix:
        if not args.genome_list:
            p.error('--matrix 需要 -l/--genome-list |requires genome-list')
        n = generate_pav_matrix(args.pav_gff, args.genome_list, args.matrix)
        print(f"矩阵|matrix: {n} 行|rows → {args.matrix}")


if __name__ == '__main__':
    main()
