"""群体基因型判定核心（纯函数）|Population genotype classification (pure functions)

本模块仅包含纯函数/纯数据类，无文件 IO、无子进程调用。
被 vcf_extractor 构建 IndelRecord，被 main 调用 classify_indel。
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class IndelRecord:
    """INDEL记录|INDEL record

    genotypes: sample_id -> 原始GT字符串（如 '1/1', '0|0', './.', '1/1:DP:10'）
    统计时由 parse_gt 解析为分类
    """
    chrom: str
    pos: int
    end: int
    ref: str
    alt: str
    indel_type: str        # insertion / deletion
    indel_size: int
    genotypes: Dict[str, str] = field(default_factory=dict)


@dataclass
class Candidate:
    """候选分子标记|Candidate marker"""
    indel: IndelRecord
    direction: str          # resistant_specific / susceptible_specific
    present_group: str      # resistant / susceptible
    r_stats: Dict[str, float] = field(default_factory=dict)
    s_stats: Dict[str, float] = field(default_factory=dict)


def parse_gt(gt_str: str) -> str:
    """
    解析GT字段为分类|Parse GT field into category

    支持形如 '0/0', '0|0', '1/1', '2/2', '0/1', '1|0', '1/2', './.' 的 GT，
    以及带 DP/AD 等后缀的字符串（仅取冒号前的 GT 部分）。

    Returns: 'hom_ref' | 'hom_alt' | 'het' | 'missing'
    """
    # 只取GT（去掉DP/AD等）|take GT part only
    gt = gt_str.split(':')[0]
    alleles = re.split(r'[|/]', gt)
    if '.' in alleles or '' in alleles:
        return 'missing'
    try:
        ints = [int(a) for a in alleles]
    except ValueError:
        return 'missing'
    if all(a == 0 for a in ints):
        return 'hom_ref'
    if len(set(ints)) == 1 and ints[0] > 0:
        return 'hom_alt'
    return 'het'


def compute_group_stats(genotypes: Dict[str, str], group_samples: List[str]) -> Dict[str, float]:
    """
    计算组内基因型统计|Compute within-group genotype stats

    分母为组内全部样品数（缺失样品也计入）。

    Returns: dict 含 hom_alt_rate / hom_ref_rate / het_rate / missing_rate / n
    """
    n = len(group_samples)
    if n == 0:
        return {'hom_alt_rate': 0.0, 'hom_ref_rate': 0.0, 'het_rate': 0.0,
                'missing_rate': 1.0, 'n': 0}
    counts = {'hom_alt': 0, 'hom_ref': 0, 'het': 0, 'missing': 0}
    for s in group_samples:
        raw = genotypes.get(s)
        if raw is None:
            counts['missing'] += 1
        else:
            counts[parse_gt(raw)] += 1
    return {
        'hom_alt_rate': counts['hom_alt'] / n,
        'hom_ref_rate': counts['hom_ref'] / n,
        'het_rate': counts['het'] / n,
        'missing_rate': counts['missing'] / n,
        'n': n,
    }


def classify_indel(indel: IndelRecord, r_samples: List[str], s_samples: List[str],
                   config) -> Optional[Candidate]:
    """
    共显性准宽松判定|Codominant loose classification

    判定规则：
      - "present" 方向：该组 hom_alt_rate >= min_group_consistency
        且 missing_rate <= max_missing_rate；
      - "absent" 方向：该组 hom_ref_rate >= min_group_consistency
        且 missing_rate <= max_missing_rate。
      - 仅当恰好一组 present + 另一组 absent 时返回 Candidate；否则返回 None。
        （两组同时满足 == 数据矛盾；两组都不满足 == 不分离；均跳过。）
    """
    r_stats = compute_group_stats(indel.genotypes, r_samples)
    s_stats = compute_group_stats(indel.genotypes, s_samples)
    mc = config.min_group_consistency
    mm = config.max_missing_rate

    def present_ok(stats):
        return stats['hom_alt_rate'] >= mc and stats['missing_rate'] <= mm

    def absent_ok(stats):
        return stats['hom_ref_rate'] >= mc and stats['missing_rate'] <= mm

    r_present = present_ok(r_stats) and absent_ok(s_stats)
    s_present = present_ok(s_stats) and absent_ok(r_stats)

    # 两组都"present"或都不"present"时跳过（前者为矛盾，后者为不分离）
    if r_present == s_present:
        return None
    if r_present:
        return Candidate(indel=indel, direction='resistant_specific',
                         present_group='resistant', r_stats=r_stats, s_stats=s_stats)
    return Candidate(indel=indel, direction='susceptible_specific',
                     present_group='susceptible', r_stats=r_stats, s_stats=s_stats)
