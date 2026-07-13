"""
ps-gene-anno 证据扫描|miniprot 蛋白→基因组(不做去冗余)+ GFF3 解析
|Evidence scan: miniprot protein-to-genome (no dedup) + GFF3 parsing
"""

import os
import re
from dataclasses import dataclass, field
from typing import List, Tuple, Optional


@dataclass
class MiniprotHit:
    """miniprot 蛋白命中|miniprot protein hit"""
    query_id: str
    chrom: str
    start: int                              # 1-based inclusive
    end: int                                # 1-based inclusive
    strand: str                             # '+' or '-'
    cds_exons: List[Tuple[int, int, int]] = field(default_factory=list)  # [(start,end,phase)] 1-based
    identity: float = 0.0                   # %
    coverage: float = 0.0                   # %


def _parse_attr(attr_str: str) -> dict:
    """解析 GFF3 第9列 attributes(ID=...;Parent=...)|Parse GFF3 column 9"""
    attrs = {}
    for kv in attr_str.split(';'):
        kv = kv.strip()
        if '=' in kv:
            k, v = kv.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _get_float_attr(attrs: dict, names: List[str]) -> float:
    """多字段名容错取 float(miniprot 版本差异)|Tolerant float lookup"""
    for n in names:
        if n in attrs:
            try:
                return float(attrs[n])
            except ValueError:
                pass
    return 0.0


def parse_miniprot_gff3(gff3_path: str,
                        min_identity: float,
                        min_coverage: float) -> List[MiniprotHit]:
    """
    解析 miniprot --gff 输出 → List[MiniprotHit]|Parse miniprot GFF3 to hits

    真实 miniprot GFF3 格式要点(已用实测数据校准):
    |Real miniprot GFF3 notes (calibrated with real data):
    - Identity=0~1 小数(需 ×100 转百分比)|Identity is 0-1 fraction, ×100 to percent
    - 无 Coverage 字段; 从 ##PAF 行 query 全长 + Target 范围计算
    |No Coverage field; compute from ##PAF query length + Target span
    - feature: mRNA + CDS(+ 可选 stop_codon)

    过滤 identity < min_identity 或 coverage < min_coverage 的命中。
    |Filter hits below identity/coverage thresholds.
    """
    if not os.path.exists(gff3_path):
        return []

    query_len: dict = {}            # query_id → 蛋白全长(from ##PAF)|prot full len
    mrna_meta: dict = {}            # tid → meta
    cds_by_parent: dict = {}        # parent_tid → [(start,end,phase)]

    with open(gff3_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('##PAF'):
                # ##PAF query qlen qstart qend strand target ...
                cols = line.split()
                if len(cols) >= 3:
                    try:
                        query_len[cols[1]] = int(cols[2])
                    except ValueError:
                        pass
                continue
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom, src, feat, start, end, score, strand, phase, attr = cols[:9]
            attrs = _parse_attr(attr)
            if feat == 'mRNA':
                tid = attrs.get('ID', '')
                query_id = ''
                tstart = tend = 0
                if 'Target' in attrs:
                    tparts = attrs['Target'].split()
                    query_id = tparts[0] if tparts else ''
                    if len(tparts) >= 3:
                        try:
                            tstart, tend = int(tparts[1]), int(tparts[2])
                        except ValueError:
                            pass
                # Identity 0-1 小数 → ×100 百分比|Identity 0-1 → ×100
                identity = _get_float_attr(attrs, ['Identity', 'identity']) * 100.0
                # coverage = Target 覆盖残基 / query 全长|coverage from Target span
                qlen = query_len.get(query_id, 0)
                if qlen > 0 and tend >= tstart:
                    coverage = (tend - tstart + 1) / qlen * 100.0
                else:
                    coverage = 0.0
                mrna_meta[tid] = {
                    'chrom': chrom, 'start': int(start), 'end': int(end),
                    'strand': strand, 'identity': identity,
                    'coverage': coverage, 'query_id': query_id,
                }
            elif feat == 'CDS':
                parent = attrs.get('Parent', '')
                try:
                    ph = int(phase) if phase != '.' else 0
                except ValueError:
                    ph = 0
                cds_by_parent.setdefault(parent, []).append(
                    (int(start), int(end), ph))

    # 组装 MiniprotHit + 过滤|Assemble + filter
    hits = []
    for tid, meta in mrna_meta.items():
        cds = sorted(cds_by_parent.get(tid, []), key=lambda x: x[0])
        hit = MiniprotHit(
            query_id=meta['query_id'], chrom=meta['chrom'],
            start=meta['start'], end=meta['end'], strand=meta['strand'],
            cds_exons=cds, identity=meta['identity'], coverage=meta['coverage'])
        if hit.identity >= min_identity and hit.coverage >= min_coverage:
            hits.append(hit)
    return hits


def run_miniprot(genome: str, prot_seq: str, out_gff3: str,
                 config, cmd_runner, logger) -> bool:
    """
    运行 miniprot 蛋白→基因组比对(--gff, 不做去冗余)
    |Run miniprot protein-to-genome (--gff, no dedup)

    Args:
        genome: 未 mask 原始基因组|Unmasked raw genome
        prot_seq: 蛋白序列|Protein sequences
        out_gff3: 输出 GFF3 路径|Output GFF3 path
        config/miniprot_bin/threads 提供|provides miniprot_bin/threads
        cmd_runner: 命令执行器(自动 conda 包装)|command runner (auto conda wrap)
        logger: 日志器|logger

    Returns:
        bool: 是否成功|success
    """
    logger.info(f"miniprot 蛋白扫描(不做去冗余)|miniprot scan (no dedup): "
                f"{os.path.basename(prot_seq)} -> {os.path.basename(genome)}")
    cmd = (f"{config.miniprot_bin} --gff -t {config.threads} "
           f"{genome} {prot_seq} > {out_gff3}")
    if not cmd_runner.run_command(cmd, "miniprot 蛋白→基因组比对|miniprot protein alignment"):
        logger.error("miniprot 失败|miniprot failed")
        return False
    logger.info(f"miniprot 输出|miniprot output: {out_gff3}")
    return True
