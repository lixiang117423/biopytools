"""
CIM分析数据转换模块|CIM Analysis Data Conversion Module

负责VCF解析、标记过滤、LD降维、tidy文件生成和R脚本生成
Handles VCF parsing, marker filtering, LD pruning, tidy file generation and R script generation
"""

import gzip
import logging
import os
import shlex
import shutil
import sys
import tempfile
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import CIMConfig
from .utils import CommandRunner, build_conda_command


def parse_vcf_genotypes(vcf_path: str, pheno_file: str,
                        logger: logging.Logger) -> Tuple[np.ndarray, pd.DataFrame, List[str], np.ndarray]:
    """
    解析VCF文件提取基因型|Parse VCF file to extract genotypes

    Args:
        vcf_path: VCF文件路径|VCF file path
        pheno_file: 表型文件路径|Phenotype file path
        logger: 日志器|Logger

    Returns:
        tuple: (genotype_matrix, marker_info, samples, phenotype_array)
            genotype_matrix: numpy array (n_markers x n_samples), 1=AA, 2=AB, 3=BB, NA=missing
            marker_info: DataFrame with chr, pos, id columns
            samples: list of sample names
            phenotype_array: numpy array of phenotype values
    """
    logger.info("步骤1: 解析VCF和表型文件|Step 1: Parsing VCF and phenotype files")

    # 读取表型文件|Read phenotype file
    pheno_df = pd.read_csv(pheno_file, sep='\t', usecols=range(2), header=0, on_bad_lines='skip')
    if 'sample' not in pheno_df.columns or 'value' not in pheno_df.columns:
        logger.error("表型文件格式错误，需要sample和value列|Phenotype file format error, need 'sample' and 'value' columns")
        raise ValueError("表型文件必须包含sample和value列|Phenotype file must have 'sample' and 'value' columns")

    pheno_samples = pheno_df['sample'].astype(str).tolist()
    pheno_values = pheno_df['value'].values
    pheno_index = {s: i for i, s in enumerate(pheno_samples)}
    n_samples = len(pheno_samples)

    logger.info(f"表型样本数|Phenotype samples: {n_samples}")

    # 解析VCF|Parse VCF
    opener = gzip.open if vcf_path.endswith('.gz') else open
    markers_chr = []
    markers_pos = []
    markers_id = []
    genotype_rows = []
    vcf_samples = None
    sample_indices = None
    total_lines = 0

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                vcf_samples = fields[9:]
                # 建立VCF样本到表型样本的索引映射
                sample_indices = []
                missing_samples = []
                for s in vcf_samples:
                    if s in pheno_index:
                        sample_indices.append(pheno_index[s])
                    else:
                        sample_indices.append(None)
                        missing_samples.append(s)
                if missing_samples:
                    logger.warning(f"VCF中{len(missing_samples)}个样本不在表型文件中，已跳过|"
                                   f"{len(missing_samples)} VCF samples not in phenotype file, skipped: "
                                   f"{', '.join(missing_samples[:10])}"
                                   f"{'...' if len(missing_samples) > 10 else ''}")
                continue

            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue

            total_lines += 1
            chrom = fields[0]
            pos = int(fields[1])
            snp_id = fields[2] if fields[2] != '.' else f"{chrom}_{pos}"

            # 提取GT字段|Extract GT field
            genotypes = np.full(n_samples, np.nan, dtype=np.float64)
            for j, gt_field in enumerate(fields[9:]):
                if j >= len(sample_indices):
                    break
                idx = sample_indices[j]
                if idx is None:
                    continue
                gt = gt_field.split(':')[0]
                if gt == '0/0' or gt == '0|0':
                    genotypes[idx] = 1  # AA
                elif gt == '0/1' or gt == '0|1' or gt == '1/0' or gt == '1|0':
                    genotypes[idx] = 2  # AB
                elif gt == '1/1' or gt == '1|1':
                    genotypes[idx] = 3  # BB
                # else: ./. or other → NaN (missing)

            markers_chr.append(chrom)
            markers_pos.append(pos)
            markers_id.append(snp_id)
            genotype_rows.append(genotypes)

    genotype_matrix = np.vstack(genotype_rows)  # (n_markers x n_samples)
    marker_info = pd.DataFrame({
        'chr': markers_chr,
        'pos': markers_pos,
        'id': markers_id
    })

    logger.info(f"VCF解析完成|VCF parsed: {total_lines} markers, {n_samples} samples")
    logger.info(f"基因型缺失率|Genotype missing rate: {np.isnan(genotype_matrix).mean():.4f}")
    logger.info(f"染色体列表|Chromosomes: {sorted(marker_info['chr'].unique())}")

    return genotype_matrix, marker_info, pheno_samples, pheno_values


def filter_markers(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                   maf: float, max_missing: float,
                   logger: logging.Logger) -> Tuple[np.ndarray, pd.DataFrame, Dict]:
    """
    基于MAF和缺失率过滤标记|Filter markers by MAF and missing rate

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        maf: 最小等位基因频率|Minimum allele frequency
        max_missing: 最大缺失率|Maximum missing rate
        logger: 日志器|Logger

    Returns:
        tuple: (filtered_genotype_matrix, filtered_marker_info, stats_dict)
    """
    logger.info(f"步骤2: 标记过滤 (MAF>={maf}, missing<={max_missing})|Step 2: Marker filtering")

    n_total = genotype_matrix.shape[0]
    stats = {'total_markers': n_total}

    # 缺失率过滤|Missing rate filter
    missing_rate = np.isnan(genotype_matrix).mean(axis=1)
    keep_missing = missing_rate <= max_missing
    n_after_missing = keep_missing.sum()
    stats['after_missing_filter'] = int(n_after_missing)
    stats['removed_by_missing'] = int(n_total - n_after_missing)
    logger.info(f"缺失率过滤: {n_total} → {n_after_missing} (移除|removed {stats['removed_by_missing']})")

    # MAF过滤|MAF filter
    # 在已过滤缺失率的矩阵上计算MAF
    # Calculate MAF on the missing-filtered matrix
    filtered_geno = genotype_matrix[keep_missing]

    n_called = (~np.isnan(filtered_geno)).sum(axis=1)
    n_called = np.maximum(n_called, 1)  # 避免除零|Avoid division by zero

    b_count = np.zeros(filtered_geno.shape[0])
    b_count += np.nansum(filtered_geno == 2, axis=1)  # AB贡献1个B
    b_count += 2 * np.nansum(filtered_geno == 3, axis=1)  # BB贡献2个B

    b_freq = b_count / (2 * n_called)
    maf_val = np.minimum(b_freq, 1 - b_freq)

    keep_maf_in_filtered = maf_val >= maf

    n_after_maf = keep_maf_in_filtered.sum()
    stats['after_maf_filter'] = int(n_after_maf)
    stats['removed_by_maf'] = int(n_after_missing - n_after_maf)
    logger.info(f"MAF过滤: {n_after_missing} → {n_after_maf} (移除|removed {stats['removed_by_maf']})")

    # 应用过滤: 先missing再MAF|Apply filters: missing then MAF
    filtered_genotype = filtered_geno[keep_maf_in_filtered]
    filtered_marker = marker_info[keep_missing].reset_index(drop=True)
    filtered_marker = filtered_marker[keep_maf_in_filtered].reset_index(drop=True)

    # 移除全NA的样本列|Remove completely NA sample columns
    sample_missing = np.isnan(filtered_genotype).all(axis=0)
    keep_sample_mask = ~sample_missing
    if sample_missing.any():
        n_rm_samples = sample_missing.sum()
        logger.warning(f"移除{int(n_rm_samples)}个全缺失样本|Removing {int(n_rm_samples)} samples with all missing genotypes")
        filtered_genotype = filtered_genotype[:, keep_sample_mask]
    stats['keep_sample_mask'] = keep_sample_mask

    stats['final_markers'] = filtered_genotype.shape[0]
    stats['final_samples'] = filtered_genotype.shape[1]
    logger.info(f"过滤完成|Filtering done: {stats['final_markers']} markers, {stats['final_samples']} samples")

    return filtered_genotype, filtered_marker, stats


def _remove_local_hotspots(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                           dist_thr: int, soft_rf: float, hard_rf: float,
                           relative_factor: float, score_cut: int, min_markers: int,
                           logger: logging.Logger,
                           qc_dir: str = "") -> Tuple[np.ndarray, Dict]:
    """短距离重组热点过滤|Remove short-distance recombination hotspots.

    双亲群体染色体交换是大片交换,物理距离很近(<dist_thr bp)的相邻标记间若出现高频
    重组(RF>=阈值),几乎必为基因型错误/旁系同源比对,删除。两轮:
    轮1 硬阈值: 单对RF>=hard_rf → 两侧标记都删(无法判断哪个错,且<thr内标记
              对QTL分辨率完全冗余); 轮2 软评分: 基于原始相邻对评分, 每个未删标记
              与原始左右邻居中RF>=soft_rf且距离<thr的对数>=score_cut → 删
              (hard标记自己不参与计分, 但仍是邻居关系的供分方)。
    |Large-block exchange in biparental crosses: high RF between physically close
    neighbors is almost surely a genotyping/paralog error. Round 1 (hard):
    RF>=hard_rf deletes BOTH markers; round 2 (soft): scores computed on original
    adjacency pairs, delete markers scoring >= score_cut (hard-flagged markers
    receive no score but still feed scores to their neighbors).

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息 (chr/pos/id)|Marker info
        dist_thr: 物理距离阈值bp, 距离更近的相邻对才检查|Only pairs closer than this
        soft_rf: 相邻RF软阈值, 达到即给两侧标记各记1分|Soft RF (scores both neighbors)
        hard_rf: 相邻RF硬阈值, 达到即两侧标记都删|Hard RF (delete both)
        relative_factor: 相对判据系数; >0时阈值=max(绝对值, 系数*近距对RF中位数)|Relative floor
        score_cut: 标记级删除线, 软评分>=此值删除|Delete markers with score >= this
        min_markers: 只处理标记数>=此值的染色体(小scaffold跳过)|Skip tiny scaffolds
        logger: 日志器|Logger
        qc_dir: QC输出目录, 非空且有删除时导出诊断TSV|QC dir for diagnostics TSV

    Returns:
        tuple: (remove掩码 True=删除, stats)|(boolean removal mask, stats dict)
    """
    n_total = genotype_matrix.shape[0]
    stats = {'markers_before_hotspot': n_total}
    hard_flags = np.zeros(n_total, dtype=bool)
    score = np.zeros(n_total, dtype=int)
    reasons = {}   # idx -> 'hard_pair' | 'score'
    max_rf_seen = {}  # idx -> 最高相邻RF(诊断用)|max adjacent RF (diagnostics)

    for chr_val in marker_info['chr'].unique():
        chr_idx = np.where(marker_info['chr'].values == chr_val)[0]
        if len(chr_idx) < max(2, min_markers):
            continue   # 小scaffold不处理|skip tiny scaffolds
        pos = marker_info['pos'].values[chr_idx].astype(np.int64)
        order = np.argsort(pos, kind='stable')   # 防御性排序|defensive sort
        chr_idx = chr_idx[order]
        pos = pos[order]
        n = len(chr_idx)

        adj_dist = pos[1:] - pos[:-1]
        near = adj_dist < dist_thr
        adj_rf = np.full(n - 1, np.nan)
        for k in range(n - 1):
            rf = _calc_single_rf(genotype_matrix[chr_idx[k]],
                                 genotype_matrix[chr_idx[k + 1]])
            adj_rf[k] = rf if rf is not None else np.nan

        # 相对判据: 阈值 = max(绝对值, 系数*近距对RF中位数)|relative: floor = factor*median
        med = float(np.nanmedian(adj_rf[near])) if near.any() else 0.0   # 空切片守卫|empty-slice guard
        soft, hard = soft_rf, hard_rf
        if relative_factor > 0:
            soft = max(soft_rf, relative_factor * med)
            hard = max(hard_rf, soft)
        logger.info(f"热点过滤 {chr_val}: {n} markers, 近距对{int(near.sum())}, "
                    f"soft={soft:.3f} hard={hard:.3f} (近距对RF中位数={med:.3f})"
                    f"|hotspot {chr_val}: {n} markers, {int(near.sum())} close pairs")

        # 轮1: 硬阈值对 → 两侧都删|round 1: hard pairs delete both
        for k in np.where(near & (adj_rf >= hard))[0]:
            hard_flags[chr_idx[k]] = True
            hard_flags[chr_idx[k + 1]] = True
            reasons.setdefault(chr_idx[k], 'hard_pair')
            reasons.setdefault(chr_idx[k + 1], 'hard_pair')
            max_rf_seen[chr_idx[k]] = max(max_rf_seen.get(chr_idx[k], 0), adj_rf[k])
            max_rf_seen[chr_idx[k + 1]] = max(max_rf_seen.get(chr_idx[k + 1], 0), adj_rf[k])

        # 轮2: 软评分(基于原始相邻对; hard标记不参与计分)|round 2: soft scoring on original adjacency
        for k in np.where(near & (adj_rf >= soft) & ~(adj_rf >= hard))[0]:
            a, b = chr_idx[k], chr_idx[k + 1]
            if not hard_flags[a]:
                score[a] += 1
            if not hard_flags[b]:
                score[b] += 1
            max_rf_seen[a] = max(max_rf_seen.get(a, 0), adj_rf[k])
            max_rf_seen[b] = max(max_rf_seen.get(b, 0), adj_rf[k])

    remove = hard_flags | (score >= score_cut)
    for i in np.where(remove)[0]:
        reasons.setdefault(int(i), 'score')
    stats['removed_by_hotspot'] = int(remove.sum())
    stats['removed_by_hard_pair'] = int(hard_flags.sum())
    stats['removed_by_soft_score'] = int(remove.sum() - hard_flags.sum())

    # 诊断TSV|diagnostics TSV
    if remove.any() and qc_dir:
        removed_idx = np.where(remove)[0]
        diag = pd.DataFrame({
            'marker_id': marker_info.loc[remove, 'id'].values,
            'chr': marker_info.loc[remove, 'chr'].values,
            'pos': marker_info.loc[remove, 'pos'].values,
            'reason': [reasons.get(int(i), 'score') for i in removed_idx],
            'score': score[remove],
            'max_adj_rf': [round(float(max_rf_seen.get(int(i), np.nan)), 4)
                            for i in removed_idx],
        })
        diag_path = os.path.join(qc_dir, 'local_hotspot_removed.tsv')
        diag.to_csv(diag_path, sep='\t', index=False)
        logger.info(f"热点删除清单已保存|Hotspot removal report saved: {diag_path}")

    stats['markers_after_hotspot'] = int((~remove).sum())
    logger.info(f"热点过滤|Hotspot filter: {n_total} → {stats['markers_after_hotspot']} "
                f"(硬对删{int(hard_flags.sum())}, 评分删{int(remove.sum() - hard_flags.sum())})")
    # remove: True=删除的布尔掩码|boolean mask, True = remove
    return remove, stats


def thin_markers_by_distance(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                             min_gap: int, logger: logging.Logger
                             ) -> Tuple[np.ndarray, pd.DataFrame, Dict]:
    """物理距离降采样: 相邻距离<min_gap的密集簇仅保留簇头。

    66万SNP相邻中位距离仅37-55bp,对99个体的QTL分辨率(~1cM)完全冗余;
    稀疏化后MSTmap建图与CIM置换检验耗时/内存降一个量级。
    |Physical thinning: keep only the first marker of any cluster denser than
    min_gap. Adjacent median is 37-55bp - fully redundant for ~1cM resolution.

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息 (chr/pos/id)|Marker info
        min_gap: 最小相邻物理距离bp|Min physical gap in bp
        logger: 日志器|Logger

    Returns:
        tuple: (降采样后矩阵, 降采样后marker_info, stats)|(thinned matrix, marker info, stats)
    """
    n_total = genotype_matrix.shape[0]
    keep = np.ones(n_total, dtype=bool)
    per_chr = {}
    for chr_val in marker_info['chr'].unique():
        chr_idx = np.where(marker_info['chr'].values == chr_val)[0]
        pos = marker_info['pos'].values[chr_idx].astype(np.int64)
        order = np.argsort(pos, kind='stable')
        chr_idx = chr_idx[order]
        pos = pos[order]
        last_kept = -10 ** 15
        for i, p in zip(chr_idx, pos):
            if p - last_kept >= min_gap:
                last_kept = int(p)   # 保留簇头|keep cluster head
            else:
                keep[i] = False
        per_chr[chr_val] = {'before': len(chr_idx),
                            'after': int(keep[chr_idx].sum())}
    filtered_genotype = genotype_matrix[keep]
    filtered_marker = marker_info[keep].reset_index(drop=True)
    stats = {'markers_before_thin': n_total,
             'markers_after_thin': filtered_genotype.shape[0],
             'min_gap': min_gap,
             'per_chromosome': per_chr}
    logger.info(f"物理降采样|Physical thinning (>= {min_gap}bp): "
                f"{n_total} → {filtered_genotype.shape[0]} markers")
    for c, d in per_chr.items():
        if d['before'] - d['after'] > 0:
            logger.info(f"  {c}: {d['before']} → {d['after']}")
    return filtered_genotype, filtered_marker, stats


def filter_rf_markers(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                      max_het_rate: float, max_mean_rf: float,
                      logger: logging.Logger,
                      qc_dir: str = "",
                      rf_knn: int = 10) -> Tuple[np.ndarray, pd.DataFrame, Dict]:
    """
    基于重组频率和杂合率过滤标记|Filter markers by recombination frequency and heterozygous rate

    三种过滤规则（按顺序执行）|Three filtering rules (executed in order):
    1. H比例过滤：杂合基因型比例过高的标记|Het rate filter: markers with too many heterozygous genotypes
    2. K近邻平均RF过滤：与物理位置前后各rf_knn//2个标记平均RF过高的标记|KNN mean RF
       filter (rf_knn//2 nearest markers on each side)
    3. 孤立重组检测：与相邻标记RF异常高（仅日志警告，不删除）|Singleton crossover detection (warning only)

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        max_het_rate: H基因型最大比例|Max heterozygous genotype rate threshold
        max_mean_rf: K近邻平均RF最大值|Max local mean RF to k nearest markers
        logger: 日志器|Logger
        qc_dir: QC输出目录|QC output directory
        rf_knn: 局部RF窗口邻居数(前后各rf_knn//2)|Neighbor count for local RF window

    Returns:
        tuple: (filtered_genotype_matrix, filtered_marker_info, stats_dict)
    """
    logger.info(f"步骤RF: 重组频率质控 (max_het_rate={max_het_rate}, "
                f"max_mean_rf={max_mean_rf}, rf_knn={rf_knn})|"
                f"Step RF: Recombination frequency QC")

    n_total = genotype_matrix.shape[0]
    stats = {'markers_before_rf': n_total}
    keep_mask = np.ones(n_total, dtype=bool)

    # --- 规则1: H比例过滤|Rule 1: Heterozygous rate filter ---
    n_valid = (~np.isnan(genotype_matrix)).sum(axis=1)
    n_valid = np.maximum(n_valid, 1)
    n_het = (genotype_matrix == 2).sum(axis=1).astype(float)
    het_rate = n_het / n_valid

    het_keep = het_rate <= max_het_rate
    n_removed_het = (~het_keep).sum()
    keep_mask &= het_keep
    stats['removed_by_het_rate'] = int(n_removed_het)
    stats['het_rate_range'] = f"{het_rate.min():.3f} - {het_rate.max():.3f}"
    stats['het_rate_median'] = round(float(np.median(het_rate)), 4)
    logger.info(f"H比例过滤|Het rate filter: {n_total} → {keep_mask.sum()} "
                f"(移除|removed {int(n_removed_het)})")
    logger.info(f"H比例分布|Het rate distribution: median={np.median(het_rate):.3f}, "
                f"range=[{het_rate.min():.3f}, {het_rate.max():.3f}]")

    # --- 规则2: 同染色体K近邻平均RF过滤|Rule 2: K-nearest-neighbor mean RF filter ---
    # 全染色体平均RF在正确RF定义下无区分力(远距离对≈0.5主导, 正常标记全局均值也
    # ≈0.45-0.5); 改为局部判据: 标记与其物理位置前后各k//2个标记的平均RF,
    # 错误标记(基因型随机)近邻RF≈0.5被删, 好标记近邻在LD内RF≈0保留
    # |Chromosome-wide mean RF has no discriminative power under the correct RF
    # definition (distant pairs ~0.5 dominate); local criterion instead: mean RF
    # to the k//2 position-nearest markers on each side. Error markers score
    # ~0.5 vs every neighbor; good markers sit in LD neighborhoods scoring ~0.
    mean_rf = _calc_local_mean_rf(genotype_matrix, marker_info, rf_knn)
    logger.info(f"局部RF窗口: 前后各{max(1, rf_knn // 2)}个标记|"
                f"Local RF window: {max(1, rf_knn // 2)} markers each side")

    rf_keep = mean_rf <= max_mean_rf
    n_removed_rf = (~rf_keep).sum()
    keep_mask &= rf_keep
    stats['removed_by_mean_rf'] = int(n_removed_rf)
    stats['mean_rf_range'] = f"{np.nanmin(mean_rf):.3f} - {np.nanmax(mean_rf):.3f}"
    stats['mean_rf_median'] = round(float(np.nanmedian(mean_rf)), 4)
    logger.info(f"平均RF过滤|Mean RF filter: {keep_mask.sum() + int(n_removed_rf) if n_removed_rf > 0 else n_total} → {keep_mask.sum()} "
                f"(移除|removed {int(n_removed_rf)})")
    logger.info(f"平均RF分布|Mean RF distribution: median={np.nanmedian(mean_rf):.3f}, "
                f"range=[{np.nanmin(mean_rf):.3f}, {np.nanmax(mean_rf):.3f}]")
    p25, p75, p95 = np.nanpercentile(mean_rf, [25, 75, 95])
    logger.info(f"平均RF分位数|Mean RF percentiles: P25={p25:.3f}, P75={p75:.3f}, P95={p95:.3f}")
    logger.info(f"过滤参考|Filter reference: "
                f"阈值0.35删{int((mean_rf > 0.35).sum())}个, "
                f"阈值0.40删{int((mean_rf > 0.40).sum())}个, "
                f"阈值0.45删{int((mean_rf > 0.45).sum())}个 "
                f"(当前阈值{max_mean_rf}删{int((mean_rf > max_mean_rf).sum())}个)")

    # --- 规则3: 孤立重组检测（仅警告，不删除）|Rule 3: Singleton detection (warning only) ---
    singleton_warnings = _detect_singletons(genotype_matrix, marker_info, keep_mask, logger)
    stats['singleton_warnings'] = len(singleton_warnings)

    # 导出孤立重组标记的H比例诊断文件|Export singleton marker H-rate diagnostics
    if singleton_warnings:
        singleton_mask = marker_info['id'].isin(singleton_warnings)
        singleton_idx = np.where(singleton_mask)[0]
        n_valid = (~np.isnan(genotype_matrix[singleton_idx])).sum(axis=1)
        n_valid = np.maximum(n_valid, 1)
        n_het = (genotype_matrix[singleton_idx] == 2).sum(axis=1).astype(float)
        het_rates = n_het / n_valid
        singleton_df = pd.DataFrame({
            'marker_id': marker_info.iloc[singleton_idx]['id'].values,
            'chr':       marker_info.iloc[singleton_idx]['chr'].values,
            'pos':       marker_info.iloc[singleton_idx]['pos'].values,
            'het_rate':  het_rates.round(4)
        }).sort_values('het_rate', ascending=False)
        if qc_dir:
            tsv_path = os.path.join(qc_dir, "singleton_het_report.tsv")
            singleton_df.to_csv(tsv_path, sep='\t', index=False)
            logger.info(f"孤立重组诊断文件已保存|Singleton report saved: {tsv_path}")
        logger.info(f"孤立重组标记H比例|Singleton H-rate: median={singleton_df['het_rate'].median():.3f}, "
                    f"高H(>0.6)共{int((singleton_df['het_rate'] > 0.6).sum())}个")

    # 应用过滤|Apply filters
    filtered_genotype = genotype_matrix[keep_mask]
    filtered_marker = marker_info[keep_mask].reset_index(drop=True)

    stats['markers_after_rf'] = filtered_genotype.shape[0]
    logger.info(f"RF质控完成|RF QC done: {n_total} → {stats['markers_after_rf']} markers")

    return filtered_genotype, filtered_marker, stats


def _calc_local_mean_rf(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                        k: int) -> np.ndarray:
    """
    每标记与物理位置最近邻居的局部平均RF|Per-marker mean RF to position-nearest neighbors

    判据(2026-08-15重构): 每标记与其**同染色体物理位置前后各 k//2 个**标记的平均RF。
    全染色体平均RF在正确RF定义下无区分力——远距离对RF≈0.5占主导, 正常标记的
    全局均值也≈0.45-0.5, 阈值无从区分好坏; 局部判据下错误标记(基因型随机)对
    所有近邻RF≈0.5被删, 好标记的近邻在LD内RF≈0保留。
    |Criterion (rebuilt 2026-08-15): mean RF to the k//2 nearest markers (by
    position) on each side, within chromosome. The chromosome-wide mean has no
    discriminative power under the correct RF definition (distant pairs ~0.5
    dominate, good markers also average ~0.45-0.5); locally, error markers
    (random genotypes) score ~0.5 vs every neighbor while good markers sit in
    LD neighborhoods scoring ~0.

    RF口径: 双纯合样本对(见_calc_single_rf), 含H/NaN样本对排除。
    |RF definition: double-homozygous pairs only (see _calc_single_rf).

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息 (chr/pos/id), 按chr分组、pos排序后取窗口|
            Marker info; window taken over pos-sorted order within each chr
        k: 邻居数, 前后各 k//2 个(染色体两端不足时用现有邻居)|Neighbor count,
            k//2 on each side (clamped at chromosome ends)

    Returns:
        ndarray: 每标记的局部平均RF(无有效邻居的标记为0)|Per-marker local mean
        RF (0 for markers without any valid neighbor)
    """
    n_total = genotype_matrix.shape[0]
    mean_rf = np.zeros(n_total)
    half = max(1, k // 2)

    for chr_val in marker_info['chr'].unique():
        chr_idx = np.where(marker_info['chr'].values == chr_val)[0]
        n = len(chr_idx)
        if n < 2:
            continue   # 单标记染色体无邻居|no neighbors on single-marker chr
        pos = marker_info['pos'].values[chr_idx].astype(np.int64)
        order = np.argsort(pos, kind='stable')   # 防御性排序|defensive sort
        chr_idx = chr_idx[order]
        g = genotype_matrix[chr_idx]

        a = (g == 1).astype(np.float64)   # AA
        b = (g == 3).astype(np.float64)   # BB
        hom = a + b                       # 纯合(NaN/H处为0)|homozygous indicator

        sum_rf = np.zeros(n)
        cnt = np.zeros(n)
        # 偏移j: 对(i, i+j)逐对累计, 覆盖前后各half个邻居|offset j covers pairs
        # (i, i+j); offsets 1..half give half neighbors on each side
        for j in range(1, min(half, n - 1) + 1):
            m = n - j
            hom_both = (hom[:m] * hom[j:]).sum(axis=1)   # 双纯合样本对数|denominator
            rec = ((a[:m] * b[j:]).sum(axis=1) +
                   (b[:m] * a[j:]).sum(axis=1))          # 重组样本数|recomb count
            with np.errstate(divide='ignore', invalid='ignore'):
                rf = rec / hom_both
            valid = hom_both > 0
            rf = np.where(valid, rf, 0.0)
            # 该对的两侧标记都累计|accumulate to both markers of each pair
            sum_rf[:m] += rf
            cnt[:m] += valid
            sum_rf[j:] += rf
            cnt[j:] += valid

        with np.errstate(divide='ignore', invalid='ignore'):
            local = np.where(cnt > 0, sum_rf / np.maximum(cnt, 1), 0.0)
        mean_rf[chr_idx] = local

    return mean_rf


def _detect_singletons(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                       keep_mask: np.ndarray, logger: logging.Logger) -> List[str]:
    """
    检测孤立重组标记（仅警告，不删除）|Detect singleton crossover markers (warning only)

    孤立重组：某标记与物理位置相邻的两个标记RF都很高（>0.5），
    说明该标记可能是基因型错误|Singleton crossover: RF > 0.5 with both neighbors

    Args:
        genotype_matrix: 基因型矩阵|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        keep_mask: 当前保留的标记掩码|Current keep mask
        logger: 日志器|Logger

    Returns:
        list: 警告的标记ID列表|List of warned marker IDs
    """
    warnings_list = []
    kept_indices = np.where(keep_mask)[0]

    for chr_val in marker_info['chr'].unique():
        chr_idx = np.where((marker_info['chr'].values == chr_val) & keep_mask)[0]
        if len(chr_idx) < 3:
            continue

        chr_geno = genotype_matrix[chr_idx]

        for k in range(1, len(chr_idx) - 1):
            g_prev = chr_geno[k - 1]
            g_curr = chr_geno[k]
            g_next = chr_geno[k + 1]

            rf_prev = _calc_single_rf(g_prev, g_curr)
            rf_next = _calc_single_rf(g_curr, g_next)

            if (rf_prev is not None and rf_prev > 0.5 and
                    rf_next is not None and rf_next > 0.5):
                marker_id = marker_info.iloc[chr_idx[k]]['id']
                warnings_list.append(marker_id)

    if warnings_list:
        logger.warning(f"孤立重组检测|Singleton detection: {len(warnings_list)} markers "
                       f"with RF>0.5 to both neighbors")
        if len(warnings_list) <= 20:
            for mid in warnings_list:
                logger.warning(f"  孤立重组标记|Singleton marker: {mid}")
        else:
            for mid in warnings_list[:10]:
                logger.warning(f"  孤立重组标记|Singleton marker: {mid}")
            logger.warning(f"  ... 省略{len(warnings_list) - 10}个|... {len(warnings_list) - 10} omitted")

    return warnings_list


def _calc_single_rf(g1: np.ndarray, g2: np.ndarray) -> Optional[float]:
    """
    计算两个标记间的重组频率|Calculate RF between two markers

    双纯合口径(2026-08-15修复): 只统计双纯合样本对; 含H的样本不提供确定的重组信息,
    原"含H记0.5"会把完全连锁标记对的RF基线抬到~0.25-0.3(H比例50%时), 阈值全部失效。
    |Double-homozygous definition (fixed 2026-08-15): count homozygous-homozygous
    sample pairs only; het samples carry no definite recombination info. The old
    "het counts 0.5" rule inflated fully-linked pairs to RF~0.25-0.3 (at 50% het),
    invalidating every threshold.

    Args:
        g1: 标记1基因型|Marker 1 genotypes
        g2: 标记2基因型|Marker 2 genotypes

    Returns:
        float or None: 重组频率，无任何双纯合样本对返回None|RF value, None if no
        double-homozygous sample pair exists
    """
    valid = ~(np.isnan(g1) | np.isnan(g2))
    both_hom = valid & (g1 != 2) & (g2 != 2)
    n = int(both_hom.sum())
    if n == 0:
        return None

    a1 = g1[both_hom] == 1
    b1 = g1[both_hom] == 3
    a2 = g2[both_hom] == 1
    b2 = g2[both_hom] == 3
    rec = ((a1 & b2) | (b1 & a2)).sum()
    return float(rec / n)


def ld_prune_markers(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                     output_dir: str, ld_window: int, ld_step: int, ld_r2: float,
                     cmd_runner: CommandRunner, logger: logging.Logger) -> Tuple[np.ndarray, pd.DataFrame, Dict]:
    """
    使用PLINK进行LD降维|LD-based marker pruning using PLINK

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        output_dir: 临时文件目录|Temporary file directory
        ld_window: LD窗口(SNP数)|LD window (SNP count)
        ld_step: LD步长(SNP数)|LD step (SNP count)
        ld_r2: r2阈值|r2 threshold
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        tuple: (pruned_genotype_matrix, pruned_marker_info, stats_dict)
    """
    logger.info(f"步骤3: LD降维 (window={ld_window}, step={ld_step}, r2={ld_r2})|Step 3: LD pruning")

    n_before = genotype_matrix.shape[0]
    stats = {'markers_before_prune': int(n_before)}

    # 创建临时工作目录(落在output_dir/tmp下，避免占满系统/tmp)|Create temp working dir under output_dir/tmp
    # mkdtemp生成唯一子目录，rmtree只删子目录，共享tmp根不会被误删|mkdtemp makes unique subdir; rmtree deletes only the subdir
    tmp_root = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_root, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix="cim_ld_prune_", dir=tmp_root)
    plink_prefix = os.path.join(tmp_dir, "plink_input")

    try:
        # 转换为PLINK格式: 缺失用NA表示 (PLINK的0在某些版本不被接受)
        # Convert to PLINK format: 1=AA, 2=AB(het), 3=BB, NA=missing
        geno_int = genotype_matrix.copy().astype(np.float64)
        geno_int[geno_int == 1] = 1  # AA → 1
        geno_int[geno_int == 2] = 2  # AB → 2
        geno_int[geno_int == 3] = 3  # BB → 3
        # NaN remains NaN → written as "0" in PED but --allow-no-sex handles this

        # 过滤掉有缺失值的标记行（PLINK不接受half-missing）
        # Filter out markers with any missing genotypes for PLINK
        valid_marker_mask = ~np.isnan(geno_int).any(axis=1)
        if valid_marker_mask.sum() < geno_int.shape[0]:
            logger.info(f"PLINK: 过滤有缺失的标记|Filtering markers with missing genotypes {geno_int.shape[0]} → {valid_marker_mask.sum()}")
            geno_int = geno_int[valid_marker_mask]
            marker_info_plink = marker_info[valid_marker_mask].reset_index(drop=True)
        else:
            marker_info_plink = marker_info

        # 写PED/MAP文件|Write PED/MAP files
        with open(f"{plink_prefix}.map", 'w') as f:
            for i in range(geno_int.shape[0]):
                f.write(f"{marker_info_plink.iloc[i]['chr']}\t{marker_info_plink.iloc[i]['id']}\t0\t{marker_info_plink.iloc[i]['pos']}\n")

        with open(f"{plink_prefix}.ped", 'w') as f:
            n_samples = geno_int.shape[1]
            for j in range(n_samples):
                sample_id = f"S{j+1}"
                geno_parts = []
                for g in geno_int[:, j]:
                    if np.isnan(g):
                        geno_parts.append("0 0")
                    elif g == 1:
                        geno_parts.append("1 1")
                    elif g == 2:
                        geno_parts.append("1 2")
                    elif g == 3:
                        geno_parts.append("2 2")
                    else:
                        geno_parts.append("0 0")
                geno_str = ' '.join(geno_parts)
                f.write(f"{sample_id}\t{sample_id}\t0\t0\t0\t-9\t{geno_str}\n")

        # 运行PLINK LD pruning|Run PLINK LD pruning
        # 经build_conda_command包装: 若plink在conda env则conda run, 否则直调(向后兼容)|
        # wrapped via build_conda_command: conda run if plink is in an env, else direct call
        plink_args = [
            '--file', plink_prefix,
            '--allow-no-sex', '--allow-extra-chr',
            '--indep-pairwise', str(ld_window), str(ld_step), str(ld_r2),
            '--out', f'{plink_prefix}_prune',
            '--silent',
        ]
        prune_cmd = ' '.join(shlex.quote(str(x)) for x in build_conda_command('plink', plink_args))
        success, stdout, stderr = cmd_runner.run(prune_cmd, "PLINK LD pruning")
        if not success:
            logger.error("PLINK LD pruning失败，回退到均匀抽样|PLINK LD pruning failed, falling back to interval sampling")
            return _fallback_interval_sampling(genotype_matrix, marker_info, 5000, stats, logger)

        # 读取pruned标记列表|Read pruned marker list
        prune_in_file = f"{plink_prefix}_prune.prune.in"
        if not os.path.exists(prune_in_file):
            logger.error("PLINK未生成prune.in文件，回退到均匀抽样|PLINK prune.in not found, falling back")
            return _fallback_interval_sampling(genotype_matrix, marker_info, 5000, stats, logger)

        prune_markers = set()
        with open(prune_in_file) as f:
            for line in f:
                prune_markers.add(line.strip())

        # 过滤到原始标记|Filter to original markers
        marker_ids = marker_info['id'].tolist()
        keep_mask = np.array([m in prune_markers for m in marker_ids])
        n_after = keep_mask.sum()
        stats['markers_after_prune'] = int(n_after)
        stats['pruning_method'] = 'LD-based (PLINK)'

        logger.info(f"LD降维完成|LD pruning done: {n_before} → {n_after} markers")

        pruned_genotype = genotype_matrix[keep_mask]
        pruned_marker = marker_info[keep_mask].reset_index(drop=True)

        return pruned_genotype, pruned_marker, stats

    finally:
        # 清理临时文件|Clean up temp files
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _fallback_interval_sampling(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                                target_count: int, stats: Dict,
                                logger: logging.Logger) -> Tuple[np.ndarray, pd.DataFrame, Dict]:
    """
    PLINK失败时的回退方案：均匀间隔抽样|Fallback: uniform interval sampling when PLINK fails

    Args:
        genotype_matrix: 基因型矩阵|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        target_count: 目标标记数|Target marker count
        stats: 统计字典|Statistics dict
        logger: 日志器|Logger

    Returns:
        tuple: (sampled_genotype, sampled_marker_info, stats)
    """
    n_total = genotype_matrix.shape[0]
    step = max(1, n_total // target_count)
    indices = list(range(0, n_total, step))

    stats['markers_after_prune'] = len(indices)
    stats['pruning_method'] = 'interval_sampling (fallback)'

    logger.info(f"均匀抽样回退|Interval sampling fallback: {n_total} → {len(indices)} markers (step={step})")

    return genotype_matrix[indices], marker_info.iloc[indices].reset_index(drop=True), stats


def _run_mstmap_single_pvalue(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                               samples: List[str], p_value: float, config: CIMConfig,
                               cmd_runner: CommandRunner,
                               logger: logging.Logger) -> List[Tuple[str, str, float]]:
    """
    以指定p_value执行MSTmap，返回所有连锁群结果|Run MSTmap with given p_value, return all LG results

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息 DataFrame|Marker info
        samples: 样本名列表|Sample names
        p_value: MSTmap cut_off_p_value|MSTmap cut_off_p_value
        config: CIM配置|CIM configuration
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        list: [(marker_id, lg_label, position), ...] 所有连锁群标记|All LG marker entries
    """
    geno_codes = {1.0: 'A', 2.0: 'X', 3.0: 'B'}
    pop_type_map = {"f2": "RIL2", "bc": "DH"}
    mstmap_pop = pop_type_map.get(config.cross_type, "RIL2")

    # 临时目录落在config.output_dir/tmp下，避免占满系统/tmp|tmp under config.output_dir/tmp
    # mkdtemp生成唯一子目录，rmtree只删子目录，共享tmp根不会被误删|mkdtemp makes unique subdir; rmtree deletes only the subdir
    tmp_root = os.path.join(config.output_dir, "tmp")
    os.makedirs(tmp_root, exist_ok=True)
    chrs = marker_info['chr'].unique()
    tmp_dir = tempfile.mkdtemp(prefix="cim_mstmap_", dir=tmp_root)
    all_map_rows = []
    lg_counter = 1

    try:
        for chr_val in chrs:
            chr_idx = np.where(marker_info['chr'].values == chr_val)[0]
            n_chr = len(chr_idx)
            if n_chr < 3:
                continue

            n_samples = genotype_matrix.shape[1]
            input_file = os.path.join(tmp_dir, f"chr_{chr_val}.txt")
            output_file = os.path.join(tmp_dir, f"chr_{chr_val}_map.txt")

            with open(input_file, 'w') as f:
                f.write(f"population_type {mstmap_pop}\n")
                f.write("population_name cross\n")
                f.write(f"distance_function {config.mstmap_distfun}\n")
                f.write(f"cut_off_p_value {p_value}\n")
                f.write("no_map_dist 20.0\n")
                f.write("no_map_size 2\n")
                f.write("missing_threshold 0.15\n")
                f.write("estimation_before_clustering yes\n")
                f.write("detect_bad_data yes\n")
                f.write("objective_function ML\n")
                f.write(f"number_of_loci {n_chr}\n")
                f.write(f"number_of_individual {n_samples}\n")
                f.write("locus_name " + " ".join(samples) + "\n")
                for idx in chr_idx:
                    marker_id = marker_info.iloc[idx]['id']
                    genos = []
                    for g in genotype_matrix[idx]:
                        if np.isnan(g):
                            genos.append('U')
                        else:
                            genos.append(geno_codes.get(g, 'U'))
                    f.write(marker_id + " " + " ".join(genos) + "\n")

            cmd = f"{config.mstmap_path} {input_file} {output_file}"
            success, stdout, stderr = cmd_runner.run(cmd, f"MSTmap {chr_val}", timeout=None)
            if not success:
                logger.warning(f"  MSTmap {chr_val} 失败，跳过|failed, skipping: {stderr[:200]}")
                continue

            chr_lgs = _parse_mstmap_output(output_file)
            n_lgs = len(chr_lgs)
            n_markers_placed = sum(len(rows) for rows in chr_lgs.values())

            for lg_name_raw, rows in chr_lgs.items():
                lg_label = f"LG{lg_counter}"
                lg_counter += 1
                for marker_id, position in rows:
                    all_map_rows.append((marker_id, lg_label, position))

            logger.info(f"    -> {n_lgs} 个连锁群|LG(s), {n_markers_placed} 个标记|markers")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return all_map_rows


def _filter_largest_lg_per_chr(all_map_rows: List[Tuple[str, str, float]],
                               marker_info: pd.DataFrame
                               ) -> Tuple[List[Tuple[str, str, float]], Set[str]]:
    """
    每染色体仅保留标记数最多的LG|Keep only the largest LG per chromosome

    MSTmap会把一条染色体拆成大LG+碎片LG(未连锁的小簇)。碎片LG标记稀疏、
    物理跨度大, R/qtl hk方法在碎片LG上会产生数值不稳定的LOD伪峰
    (实测longliuxing Chr11碎片19标记LOD=383.8, 同设置em=0.37, 且其标记
    与表型最大|corr|=0.24无真实效应)。碎片LG不进CIM, 每染色体只保留最大LG。
    MSTmap splits a chromosome into a big LG plus fragment LGs (small unlinked
    clusters). Fragments are sparse with huge physical spans; R/qtl's hk method
    produces numerically unstable bogus LOD peaks on them (measured LOD 383.8 on
    a 19-marker Chr11 fragment vs em=0.37 on the identical setup; the fragment's
    markers show max |corr|=0.24 with the phenotype). Drop fragment LGs from CIM.

    Args:
        all_map_rows: [(marker_id, lg_label, position), ...]|All LG marker entries
        marker_info: 标记信息(chr, pos, id)|Marker info DataFrame

    Returns:
        (filtered_rows, keep_lgs): 过滤后的行与保留的LG集合|Filtered rows and kept LG set
    """
    id_to_chr = marker_info.set_index('id')['chr'].to_dict()

    # 统计每个(染色体, LG)的标记数|Count markers per (chr, LG)
    chr_lg_sizes: Dict[Tuple[str, str], int] = {}
    for marker_id, lg_label, _pos in all_map_rows:
        chr_name = id_to_chr.get(marker_id)
        if chr_name is None:
            continue
        key = (chr_name, lg_label)
        chr_lg_sizes[key] = chr_lg_sizes.get(key, 0) + 1

    # 每染色体保留标记数最多的LG|Keep the biggest LG per chromosome
    biggest: Dict[str, Tuple[str, int]] = {}
    for (chr_name, lg_label), size in chr_lg_sizes.items():
        if chr_name not in biggest or size > biggest[chr_name][1]:
            biggest[chr_name] = (lg_label, size)
    keep_lgs = {lg_label for (lg_label, _size) in biggest.values()}

    # 保留LG内的标记; marker_info无染色体的标记一并剔除|Keep markers in kept LGs;
    # also drop markers whose chromosome is unknown in marker_info
    filtered_rows = [row for row in all_map_rows
                     if row[0] in id_to_chr and row[1] in keep_lgs]
    return filtered_rows, keep_lgs


def _save_map_outputs(all_map_rows: List[Tuple[str, str, float]],
                      genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                      samples: List[str], pheno_values: np.ndarray,
                      save_dir: str) -> None:
    """
    将连锁群结果保存为mstmap_map.csv, linkage_map.csv等文件|Save LG results to map files

    Args:
        all_map_rows: [(marker_id, lg_label, position), ...]|All LG marker entries
        genotype_matrix: 基因型矩阵|Genotype matrix
        marker_info: 标记信息|Marker info
        samples: 样本名列表|Sample names
        pheno_values: 表型值数组|Phenotype values
        save_dir: 保存目录|Save directory
    """
    geno_codes = {1.0: 'A', 2.0: 'H', 3.0: 'B'}

    mstmap_map_file = os.path.join(save_dir, 'mstmap_map.csv')
    linkage_map_csv = os.path.join(save_dir, 'linkage_map.csv')
    mstmap_gen_file = os.path.join(save_dir, 'mstmap_gen.csv')
    mstmap_phe_file = os.path.join(save_dir, 'mstmap_phe.csv')

    mstmap_map_df = pd.DataFrame(
        [(row[1], row[2]) for row in all_map_rows],
        columns=['chr', 'pos'],
        index=[row[0] for row in all_map_rows]
    )
    mstmap_map_df.to_csv(mstmap_map_file, index=True)

    linkage_df = mstmap_map_df.reset_index()
    linkage_df.columns = ['marker', 'chr', 'pos']
    linkage_df.to_csv(linkage_map_csv, index=False)

    # 生成mstmap_gen.csv: 行=个体, 列=retained标记 (csvs格式)
    retained_markers = mstmap_map_df.index.tolist()
    marker_id_to_row = {marker_info.iloc[i]['id']: i for i in range(len(marker_info))}
    retained_indices = [marker_id_to_row[mid] for mid in retained_markers if mid in marker_id_to_row]
    retained_geno = genotype_matrix[retained_indices]

    gen_data = []
    for i in range(len(samples)):
        row = []
        for j in range(retained_geno.shape[0]):
            g = retained_geno[j, i]
            if np.isnan(g):
                row.append('NA')
            else:
                row.append(geno_codes.get(g, 'NA'))
        gen_data.append(row)

    gen_df = pd.DataFrame(gen_data, index=samples, columns=retained_markers)
    gen_df.to_csv(mstmap_gen_file)

    phe_df = pd.DataFrame({'trait1': pheno_values[:len(samples)]}, index=samples)
    phe_df.to_csv(mstmap_phe_file)

    # 生成marker坐标对照表: 遗传坐标 ↔ 物理坐标|Generate marker coordinate index
    marker_index_file = os.path.join(save_dir, 'marker_map_index.tsv')
    marker_phys = marker_info.set_index('id')[['chr', 'pos']].copy()
    marker_phys.columns = ['chr_bp', 'pos_bp']
    index_df = mstmap_map_df.copy()
    index_df.columns = ['LG', 'pos_cM']
    index_df = index_df.join(marker_phys, how='left')
    index_df = index_df.reset_index()
    index_df.columns = ['marker', 'LG', 'pos_cM', 'chr_bp', 'pos_bp']
    index_df = index_df[['marker', 'chr_bp', 'pos_bp', 'LG', 'pos_cM']]
    index_df.to_csv(marker_index_file, sep='\t', index=False)


def build_mstmap_linkage_map(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                             samples: List[str], pheno_values: np.ndarray,
                             tidy_files: Dict[str, str],
                             output_dir: str, config: CIMConfig,
                             cmd_runner: CommandRunner,
                             logger: logging.Logger) -> Dict[str, str]:
    """
    使用C++ MSTmap构建遗传图谱，自动调优cut_off_p_value|Build genetic map with auto-tuned MSTmap

    自动迭代调优p.value，目标LG数 <= n_chr x 3。每次迭代结果保存到独立子目录。
    Auto-tune p.value with target LG count <= n_chr x 3. Each iteration saves to separate subdir.

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix (markers x samples)
        marker_info: 标记信息 DataFrame (chr, pos, id)|Marker info
        samples: 样本名列表|Sample names
        pheno_values: 表型值数组|Phenotype values
        tidy_files: tidy文件路径字典|Tidy file paths {gen_file, phe_file, map_file}
        output_dir: 输出目录|Output directory
        config: CIM配置|CIM configuration
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        dict: {mstmap_map, mstmap_gen, mstmap_phe, linkage_map} 文件路径|File paths
    """
    mstmap_map_file = os.path.join(output_dir, 'mstmap_map.csv')
    mstmap_gen_file = os.path.join(output_dir, 'mstmap_gen.csv')
    mstmap_phe_file = os.path.join(output_dir, 'mstmap_phe.csv')
    linkage_map_csv = os.path.join(output_dir, 'linkage_map.csv')

    n_chromosomes = len(marker_info['chr'].unique())
    max_lg_target = n_chromosomes * 3
    # 核心LG阈值: 标记数 > MSTmap no_map_size(_run_mstmap_single_pvalue中=2)才算真正的连锁群。
    # 孤立单标记LG是MSTmap踢出的坏标记, 不应计入"是否找回染色体"的判据。
    # Core LG threshold: only LGs with >MSTmap no_map_size(=2) markers are real linkage groups.
    # Singleton LGs are bad markers MSTmap isolated; they must NOT count against chromosome recovery.
    core_lg_min = 3

    # p.value尝试序列: 用户值, 1e-5, 1e-4, 1e-3, 1e-2
    pvalue_sequence = [config.mstmap_pvalue, 1e-5, 1e-4, 1e-3, 1e-2]
    # 去重并按从小到大排序（去掉用户值与序列中重复的）|Deduplicate and sort
    pvalue_sequence = sorted(set(pvalue_sequence))

    logger.info(f"MSTmap自动调优|Auto-tuning MSTmap: "
                f"{n_chromosomes}条染色体|chromosomes, "
                f"目标核心LG数|target core LG <= {max_lg_target}")
    logger.info(f"  p.value尝试序列|p.value sequence: "
                f"{', '.join(f'{p:.0e}' for p in pvalue_sequence)} "
                f"(核心LG阈值|core LG min markers={core_lg_min})")

    accepted_rows = None
    accepted_pvalue = None
    # 兜底最优: 选核心LG数最少(聚类最好)的迭代, 而非盲目取最后一个p.value。
    # Fallback best: pick the iteration with fewest core LGs (best clustering),
    # not blindly the last/loosest p.value.
    best_rows = None
    best_pvalue = None
    best_core_lgs = None

    for p_val in pvalue_sequence:
        lg_label = f"mstmap_pvalue_{p_val:.0e}"
        iter_dir = os.path.join(output_dir, lg_label)
        os.makedirs(iter_dir, exist_ok=True)

        logger.info(f"MSTmap自动调优|Auto-tuning: p_value={p_val:.0e} ...")

        map_rows = _run_mstmap_single_pvalue(
            genotype_matrix, marker_info, samples, p_val,
            config, cmd_runner, logger
        )

        if not map_rows:
            logger.warning(f"MSTmap自动调优|Auto-tuning: p_value={p_val:.0e} -> 0个标记|0 markers, 跳过")
            continue

        # 统计LG数: 区分核心LG(可作图)与孤立单标记LG(坏标记)|
        # Count LGs: split core LGs (mappable) from isolated singletons (bad markers)
        lg_sizes = Counter(row[1] for row in map_rows)
        n_lgs = len(lg_sizes)
        n_core_lgs = sum(1 for s in lg_sizes.values() if s >= core_lg_min)
        n_markers = len(map_rows)

        logger.info(f"MSTmap自动调优|Auto-tuning: p_value={p_val:.0e} -> "
                    f"{n_core_lgs}个核心LG|core LGs "
                    f"(共|total {n_lgs}个LG|LGs), {n_markers}个标记|markers")

        # 保存本次迭代结果|Save this iteration's results
        _save_map_outputs(map_rows, genotype_matrix, marker_info,
                          samples, pheno_values, iter_dir)

        # 更新兜底最优: 核心LG最少优先, 平手取标记多|update fallback best
        if (best_core_lgs is None or n_core_lgs < best_core_lgs
                or (n_core_lgs == best_core_lgs and n_markers > len(best_rows))):
            best_rows = map_rows
            best_pvalue = p_val
            best_core_lgs = n_core_lgs

        if n_core_lgs <= max_lg_target:
            logger.info(f"MSTmap自动调优|Auto-tuning: p_value={p_val:.0e} -> "
                        f"{n_core_lgs}个核心LG(<= {max_lg_target}), 接受|accepted")
            accepted_rows = map_rows
            accepted_pvalue = p_val
            break
        else:
            logger.info(f"MSTmap自动调优|Auto-tuning: p_value={p_val:.0e} -> "
                        f"{n_core_lgs}个核心LG(> {max_lg_target}), "
                        f"放宽至下一个p.value重试|relaxing to next p.value")

    if accepted_rows is None:
        if best_rows is None:
            # 所有p.value均聚类出0标记: MSTmap完全无法把任何标记归入连锁群。
            # 旧版在此处 f"{best_pvalue:.0e}" 对 None 格式化 → TypeError 崩溃(回归),
            # 这里改为清晰失败退出, 与 main.py 中"RF过滤后无标记"的处理保持一致。
            # |Every p.value placed 0 markers: MSTmap could not assign any marker to a
            # linkage group. The old code crashed here with TypeError on f"{None:.0e}";
            # fail cleanly instead, consistent with main.py's no-markers handling.
            logger.error(
                f"MSTmap自动调优失败: 尝试的所有p.value"
                f"({', '.join(f'{p:.0e}' for p in pvalue_sequence)}) "
                f"均未聚类出任何标记|MSTmap auto-tuning failed: every p.value placed 0 "
                f"markers. 请检查输入标记质量, 或改用 physical/estimate 模式|"
                f"Check input marker quality or switch to physical/estimate map_mode"
            )
            sys.exit(1)
        # 所有p.value的核心LG都超目标, 用聚类最好(核心LG最少)的结果, 而非最宽松p.value。
        # No p.value met the target core LG count; use the best-clustering result.
        logger.warning(f"MSTmap自动调优: 所有p.value核心LG数均超过{max_lg_target}|"
                       f"All p.values exceeded target core LG count; "
                       f"使用聚类最好结果|using best-clustering result "
                       f"(p_value={best_pvalue:.0e}, {best_core_lgs}个核心LG|core LGs)")
        accepted_rows = best_rows
        accepted_pvalue = best_pvalue

    # 每染色体仅保留最大LG进CIM|Keep only the largest LG per chromosome for CIM
    # 碎片LG标记稀疏, R/qtl hk方法在其上产生数值不稳定LOD伪峰(_filter_largest_lg_per_chr
    # 注释有实测数据), 故只把每染色体最大LG写进根目录供R脚本使用; 迭代子目录保留原始聚类结果。
    # Fragment LGs are sparse and trigger numerically unstable hk LOD artifacts (see
    # _filter_largest_lg_per_chr docstring); write only the largest LG per chromosome
    # to the root dir for the R script. Per-iteration subdirs keep the raw clustering.
    filtered_rows, keep_lgs = _filter_largest_lg_per_chr(accepted_rows, marker_info)
    all_lgs_final = set(row[1] for row in accepted_rows)
    n_dropped_lgs = len(all_lgs_final) - len(keep_lgs)
    if n_dropped_lgs > 0:
        logger.warning(f"MSTmap碎片LG过滤|Fragment LG filtering: "
                       f"每染色体仅保留最大LG进CIM|keeping only the largest LG per "
                       f"chromosome: 剔除|dropped {n_dropped_lgs}个碎片LG|fragment "
                       f"LGs ({len(accepted_rows) - len(filtered_rows)}个标记|markers), "
                       f"保留|kept {len(keep_lgs)}个LG|LGs")
    if not filtered_rows:
        logger.error(f"MSTmap碎片LG过滤后无标记|No markers left after fragment LG "
                     f"filtering. 请检查标记信息|Check marker info")
        sys.exit(1)
    _save_map_outputs(filtered_rows, genotype_matrix, marker_info,
                      samples, pheno_values, output_dir)

    lg_sizes_final = Counter(row[1] for row in accepted_rows)
    n_core_final = sum(1 for s in lg_sizes_final.values() if s >= core_lg_min)
    logger.info(f"MSTmap自动调优完成|Auto-tuning complete: "
                f"实际使用p_value={accepted_pvalue:.0e}, "
                f"{n_core_final}个核心LG|core LGs "
                f"(共|total {len(lg_sizes_final)}个LG|LGs), {len(accepted_rows)}个标记|markers")

    return {
        'mstmap_map': mstmap_map_file,
        'mstmap_gen': mstmap_gen_file,
        'mstmap_phe': mstmap_phe_file,
        'linkage_map': linkage_map_csv,
    }


def _parse_mstmap_output(output_file: str) -> Dict[str, List[Tuple[str, float]]]:
    """
    解析MSTmap输出文件|Parse MSTmap output file

    Args:
        output_file: MSTmap输出文件路径|MSTmap output file path

    Returns:
        dict: {lg_name: [(marker_id, position), ...]}
    """
    lg_data = {}
    current_lg = None
    in_group = False

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(';BEGINOFGROUP'):
                in_group = True
                continue
            if line.startswith(';ENDOFGROUP'):
                in_group = False
                current_lg = None
                continue
            if line.startswith('group '):
                current_lg = line.split(None, 1)[1]
                lg_data[current_lg] = []
                continue
            if line.startswith(';'):
                continue
            if in_group and current_lg is not None:
                parts = line.split('\t')
                if len(parts) >= 2:
                    marker_id = parts[0]
                    try:
                        position = float(parts[1])
                        lg_data[current_lg].append((marker_id, position))
                    except ValueError:
                        pass

    return lg_data


def save_filtered_vcf(original_vcf: str, marker_info: pd.DataFrame, retained_samples: List[str],
                      output_path: str, cmd_runner: CommandRunner,
                      logger: logging.Logger) -> str:
    """
    用bcftools从原始VCF中提取过滤后保留的标记和样本|Extract retained markers and samples from original VCF using bcftools

    Args:
        original_vcf: 原始VCF路径|Original VCF file path
        marker_info: 保留的标记信息|Retained marker info DataFrame (chr, pos, id)
        retained_samples: 保留的样本名列表|Retained sample name list
        output_path: 输出VCF路径|Output VCF path (.vcf.gz)
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        str: 输出文件路径，失败返回空字符串|Output file path, empty string on failure
    """
    logger.info(f"保存过滤后的VCF文件|Saving filtered VCF: {output_path}")

    n_markers = len(marker_info)
    n_samples = len(retained_samples)

    # 写入标记位置列表（bcftools -T格式：tab分隔的chr\tpos）| Write marker target list
    targets_file = output_path + ".targets.tmp"
    with open(targets_file, 'w') as f:
        for _, row in marker_info.iterrows():
            f.write(f"{row['chr']}\t{row['pos']}\n")

    # 写入样本列表（bcftools -S格式：每行一个样本名）| Write sample list
    samples_file = output_path + ".samples.tmp"
    with open(samples_file, 'w') as f:
        for s in retained_samples:
            f.write(f"{s}\n")

    try:
        # 经build_conda_command包装bcftools(向后兼容: 不在conda env则直调)|
        # bcftools wrapped via build_conda_command (backward compatible: direct call if not in env)
        view_args = ['view', '-T', targets_file, '-S', samples_file, '--force-samples',
                     '-Oz', '-o', output_path, original_vcf]
        cmd = ' '.join(shlex.quote(str(x)) for x in build_conda_command('bcftools', view_args))
        success, stdout, stderr = cmd_runner.run(cmd, "提取过滤后VCF|Extract filtered VCF")

        if not success:
            logger.warning("bcftools提取失败，跳过VCF输出|bcftools extraction failed, skipping filtered VCF output")
            return ""

        index_cmd = ' '.join(shlex.quote(str(x)) for x in build_conda_command('bcftools', ['index', output_path]))
        cmd_runner.run(index_cmd, "索引过滤后VCF|Index filtered VCF")

        logger.info(f"过滤后VCF已保存|Filtered VCF saved: {output_path} ({n_markers} markers, {n_samples} samples)")
        return output_path

    finally:
        for tmp_file in [targets_file, samples_file]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)


def build_tidy_files(genotype_matrix: np.ndarray, marker_info: pd.DataFrame,
                     samples: List[str], pheno_values: np.ndarray,
                     output_dir: str, map_mode: str,
                     logger: logging.Logger) -> Dict[str, str]:
    """
    生成R/qtl tidy格式文件|Generate R/qtl tidy format files

    Args:
        genotype_matrix: 基因型矩阵 (n_markers x n_samples)|Genotype matrix
        marker_info: 标记信息|Marker info DataFrame
        samples: 样本名列表|Sample name list
        pheno_values: 表型值数组|Phenotype values array
        output_dir: 输出目录|Output directory
        map_mode: cM模式|cM source mode (physical/estimate)
        logger: 日志器|Logger

    Returns:
        dict: 文件路径字典|File path dict {gen_file, phe_file, map_file}
    """
    logger.info("步骤4: 生成csvs格式文件|Step 4: Generating csvs format files")

    n_markers = genotype_matrix.shape[0]
    n_samples = genotype_matrix.shape[1]

    # gen.csv: 行=标记, 列=个体 (R/qtl csvs格式要求, read.cross内部处理)
    # gen.csv: rows=markers, cols=individuals (R/qtl csvs format, read.cross handles internally)
    geno_codes = {1.0: 'A', 2.0: 'H', 3.0: 'B'}
    gen_data = []
    for i in range(n_markers):
        row = [geno_codes.get(g, 'NA') if not np.isnan(g) else 'NA' for g in genotype_matrix[i, :n_samples]]
        gen_data.append(row)

    # 标记×个体 矩阵 (与旧版一致)
    gen_df = pd.DataFrame(gen_data,
                          index=marker_info['id'].tolist(),
                          columns=samples[:n_samples])
    gen_file = os.path.join(output_dir, "gen.csv")
    gen_df.to_csv(gen_file)

    # 个体×性状 (csvs格式)
    phe_df = pd.DataFrame({'trait1': pheno_values[:n_samples]},
                          index=samples[:n_samples])
    phe_file = os.path.join(output_dir, "phe.csv")
    phe_df.to_csv(phe_file)

    # map.csv: marker, chr, position(cM)
    if map_mode == "physical":
        positions = marker_info['pos'].values / 1e6  # bp → Mb (pseudo-cM)
    else:
        positions = marker_info['pos'].values / 1e6  # Mb占位; 仅estimate模式由est.map覆盖, mstmap模式的物理block沿用此Mb(不可est.map,见_build_cim_block注释)

    map_df = pd.DataFrame({
        'chr': marker_info['chr'].tolist(),
        'pos': positions
    }, index=marker_info['id'].tolist())
    map_file = os.path.join(output_dir, "map.csv")
    map_df.to_csv(map_file, index=True)

    logger.info(f"gen.csv: {n_markers} markers x {n_samples} samples")
    logger.info(f"phe.csv: {n_samples} samples x 1 trait")
    logger.info(f"map.csv: {n_markers} markers, map_mode={map_mode}")

    return {'gen_file': gen_file, 'phe_file': phe_file, 'map_file': map_file}


def _build_cim_block(config: CIMConfig, label: str, output_dir: str) -> List[str]:
    """
    生成单次CIM分析的R代码块|Generate R code block for a single CIM analysis run

    Args:
        config: CIM配置|CIM configuration
        label: 标签（如 "physical" 或 "mstmap"）|Label for output files
        output_dir: 输出根目录|Output root directory (physical/mstmap/plots 的父目录)

    Returns:
        list: R代码行列表|List of R code lines
    """
    window = config.window
    step = config.step

    plots_dir = os.path.join(output_dir, 'plots')
    label_dir = os.path.join(output_dir, label)

    results_rds = os.path.join(label_dir, 'cim_results.rds')
    perm_rds = os.path.join(label_dir, 'cim_perm_results.rds')
    lod_data_tsv = os.path.join(plots_dir, f'cim_lod_data_{label}.tsv')
    threshold_txt = os.path.join(plots_dir, f'cim_threshold_{label}.txt')
    peaks_tsv = os.path.join(plots_dir, f'cim_peaks_{label}.tsv')
    genome_pdf = os.path.join(label_dir, 'cim_genome_plot.pdf')
    chr_pdf_prefix = os.path.join(label_dir, f'cim_chr_')

    lines = []

    if label:
        lines.append(f'cat("\\n========== CIM分析 ({label}) ==========\\n")')
        lines.append('')
        if label == "physical":
            # 显式标注尺度差异:物理block图谱坐标为Mb(pseudo-cM),window/step 与
            # mstmap block(真实cM)不同尺度,两block的LOD/阈值仅在同block内可比。
            # |Mark scale mismatch explicitly: physical block map is Mb (pseudo-cM);
            # window/step differ from the mstmap block (true cM). Compare LOD/thresholds
            # only within the same block.
            lines.append(
                'cat("NOTE: 物理block图谱坐标为Mb(pseudo-cM),window/step与mstmap '
                'block(真实cM)尺度不同;两block的LOD/阈值仅在同block内可比'
                '|Physical block map is Mb (pseudo-cM); window/step differ from mstmap '
                'block (true cM); compare LOD/thresholds only within the same block\\n")')
            lines.append('')

    # 剔除孤立单标记LG(<3 markers),避免抬高置换检验阈值|
    # Remove singleton LGs (<3 markers) to avoid inflating permutation threshold
    lines.extend([
        '# 剔除孤立单标记LG(<3 markers)|Remove singleton LGs (<3 markers)',
        'cat("剔除孤立单标记LG(<3 markers)...\\n")',
        'nmar_all <- nmar(cross)',
        'singleton_lgs <- names(nmar_all)[nmar_all < 3]',
        'if (length(singleton_lgs) > 0) {',
        '    cat("  移除孤立LG|Removing singleton LGs:", length(singleton_lgs), "out of", length(nmar_all), "\\n")',
        '    cat("  移除的LG|Removed LGs:", paste(singleton_lgs, collapse=", "), "\\n")',
        '    cross <- subset(cross, chr = names(nmar_all)[nmar_all >= 3])',
        '    cat("  过滤后|After filtering:", nmar(cross), "markers in", nchr(cross), "LGs\\n")',
        '} else {',
        '    cat("  无孤立LG,无需过滤|No singleton LGs to remove\\n")',
        '}',
    ])

    # NOTE: 不对物理block调est.map()。LD降维(r2=0.1)刻意保留低LD标记,相邻标记近乎不连锁,
    # est.map会把重组率估成~0.5→Kosambi距离→∞→图谱膨胀到数百万cM→calc.genoprob生成数百万
    # pseudomarker→cim置换检验(×1000)永不完成(实测物理block 13h+未跑完,LOD畸高至305)。
    # DO NOT est.map() the physical block: LD pruning keeps markers nearly unlinked, so est.map
    # inflates the map to millions of cM -> calc.genoprob explodes -> permutation never finishes.

    lines.extend([
        '# 计算基因型概率|Calculate genotype probabilities',
        'cat("计算基因型概率...\\n")',
        'cross <- calc.genoprob(cross, step=' + str(step) + ', error.prob=0.0001)',
        'cat("基因型概率计算完成\\n")',
        '',
        '# CIM分析|CIM analysis',
        'cat("开始CIM分析...\\n")',
        'out <- cim(cross, pheno.col="trait1",',
        '           n.marcovar=' + str(config.n_marcovar) + ',',
        '           window=' + str(window) + ',',
        '           method="' + config.method + '")',
        '',
        'saveRDS(out, file="' + results_rds + '")',
        'lod_df <- as.data.frame(out)',
        'colnames(lod_df) <- c("chr", "pos", "lod")',
        'write.table(lod_df, file="' + lod_data_tsv + '",',
        '            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
        'cat("LOD扫描数据已保存|LOD scan data saved\\n")',
        'cat("CIM分析完成\\n")',
        '',
    ])

    # 置换检验|Permutation test
    threshold_lines = []
    if config.n_perm > 0:
        lines.extend([
            '# 置换检验|Permutation test',
            'cat("开始置换检验 (n=' + str(config.n_perm) + ')...\\n")',
            'set.seed(42)',
            'perm <- cim(cross, pheno.col="trait1", n.marcovar=' + str(config.n_marcovar) + ',',
            '            window=' + str(window) + ', method="' + config.method + '", n.perm=' + str(config.n_perm) + ')',
            'thr_obj <- summary(perm, alpha=0.05)',
            'threshold <- max(thr_obj)',
            'cat("显著性阈值(LOD):", threshold, "\\n")',
            'saveRDS(perm, file="' + perm_rds + '")',
            'writeLines(as.character(threshold), "' + threshold_txt + '")',
            ''
        ])
        # LOD TSV添加threshold列（方便ggplot2可视化）
        lines.extend([
            'lod_df$threshold <- threshold',
            'write.table(lod_df, file="' + lod_data_tsv + '",',
            '            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
            ''
        ])
        threshold_lines = [
            'abline(h=threshold, col="red", lty=2)',
            'mtext(paste0("Threshold: ", round(threshold, 2)), side=3, line=0, col="red", cex=0.9)',
        ]

    # 峰值表|Peak table
    if config.n_perm > 0:
        lines.extend([
            '# 提取QTL峰值|Extract QTL peaks',
            'summary_out <- summary(out)',
            'peaks <- as.data.frame(summary_out)',
            'peaks_sig <- peaks[peaks$lod >= threshold, , drop=FALSE]',
            'if (nrow(peaks_sig) > 0) {',
            '    peaks_sig <- peaks_sig[order(-peaks_sig$lod), , drop=FALSE]',
            '}',
            'write.table(peaks_sig, file="' + peaks_tsv + '",',
            '            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
            'cat("显著QTL峰值已保存|Significant QTL peaks saved\\n")',
            ''
        ])
    else:
        lines.extend([
            '# 提取QTL峰值|Extract QTL peaks',
            'summary_out <- summary(out)',
            'peaks <- as.data.frame(summary_out)',
            'peaks_sorted <- peaks[order(-peaks$lod), , drop=FALSE]',
            'write.table(peaks_sorted, file="' + peaks_tsv + '",',
            '            sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
            'cat("QTL峰值已保存(未做置换检验)|QTL peaks saved (no permutation test)\\n")',
            ''
        ])

    # 全基因组LOD图|Genome-wide LOD plot
    plot_title = f'Composite Interval Mapping (CIM) - {label}' if label else 'Composite Interval Mapping (CIM)'
    lines.extend([
        '# 全基因组LOD曲线图|Genome-wide LOD plot',
        'pdf("' + genome_pdf + '", width=12, height=6)',
        'plot(out, main="' + plot_title + '",',
        '     ylab="LOD score")',
    ])
    lines.extend(threshold_lines)
    lines.extend([
        '# 标注协因子位置|Annotate cofactor positions',
        'tryCatch({ add.cim.covar(out, pch=16, col="blue") },',
        '         error = function(e) { cat("add.cim.covar skipped:", e$message, "\\n") })',
        'dev.off()',
        'cat("全基因组LOD曲线图已保存|Genome-wide LOD plot saved\\n")',
        '',
        '# 逐染色体LOD曲线图|Per-chromosome LOD plots',
        '# 跳过孤立单标记LG(<3 markers), 不生成无意义的单标记PDF|',
        '# skip isolated singleton LGs (<3 markers); avoid junk single-marker PDFs',
        'chrs <- unique(out[, 1])',
        'nmar_vec <- nmar(cross)',
        'for (chr_i in chrs) {',
        '    chr_label <- as.character(chr_i)',
        '    if (is.na(nmar_vec[chr_label]) || nmar_vec[chr_label] < 3) next',
        '    pdf_file <- paste0("' + chr_pdf_prefix + '", chr_label, ".pdf")',
        '    pdf(pdf_file, width=8, height=5)',
        '    plot(out, chr=chr_i, main=paste("CIM -", chr_label),',
        '         ylab="LOD score")',
    ])
    lines.extend(['    ' + l for l in threshold_lines])
    lines.extend([
        '    tryCatch({ add.cim.covar(out, chr=chr_i, pch=16, col="blue") },',
        '             error = function(e) { cat("add.cim.covar skipped:", e$message, "\\n") })',
        '    dev.off()',
        '}',
        'cat("逐染色体LOD曲线图已保存|Per-chromosome LOD plots saved\\n")',
        '',
    ])

    return lines


def generate_r_cim_script(config: CIMConfig, tidy_files: Dict[str, str],
                          output_dir: str) -> str:
    """
    生成CIM分析的R脚本|Generate R script for CIM analysis

    Args:
        config: CIM配置|CIM configuration
        tidy_files: tidy文件路径字典|Tidy file path dict
        output_dir: 输出目录|Output directory

    Returns:
        str: R脚本路径|R script file path
    """
    r_script = os.path.join(output_dir, "cim_analysis.R")

    # 路径变量|Path variables
    gen_file = tidy_files['gen_file']
    phe_file = tidy_files['phe_file']
    map_file = tidy_files['map_file']

    # 公共头部：读取数据 + 去除高缺失个体|Common header: read data + remove high-missing individuals
    lines = [
        '#!/usr/bin/env Rscript',
        '# CIM Analysis Script - Generated by biopytools cim',
        '# ' + config.map_mode + ' map mode, ' + config.cross_type + ' cross',
        '',
        'library(qtl)',
    ]

    if config.map_mode == "mstmap":
        pass  # ASMap replaced by MSTmap (Python), no R library needed

    lines.extend([
        '',
        '# 构建R/qtl csvs格式genfile|Build R/qtl csvs format genfile',
        '# csvs格式要求: 行1=标记名, 行2=染色体, 行3=位置, 行4+=个体基因型',
        'cat("构建csvs格式文件...\\n")',
        '',
        'geno_raw <- read.csv("' + gen_file + '", row.names = 1, check.names = FALSE)',
        'map_df <- read.csv("' + map_file + '", header = TRUE)',
        'marker_names <- rownames(geno_raw)',
        'm_idx <- match(marker_names, as.character(map_df[,1]))',
        'chr_vec <- as.character(map_df[m_idx, 2])',
        'pos_vec <- as.character(map_df[m_idx, 3])',
        'geno_t <- t(as.matrix(geno_raw))',
        'geno_t <- cbind(rownames(geno_t), geno_t)',
        '',
        '# 拼装: header行 + chr行 + pos行 + 个体基因型行',
        'chr_row <- c("", chr_vec)',
        'pos_row <- c("", pos_vec)',
        'gen_df <- rbind(chr_row, pos_row, geno_t)',
        'colnames(gen_df) <- c("", marker_names)',
        '',
        'tmp_gen <- tempfile(fileext = ".csv")',
        'write.table(gen_df, file = tmp_gen, sep = ",", row.names = FALSE,',
        '            col.names = TRUE, quote = FALSE)',
        '',
        '# 构建cross对象|Build cross object from csvs files',
        'cross <- read.cross(format="csvs",',
        '                    genfile = tmp_gen,',
        '                    phefile = "' + phe_file + '",',
        '                    genotypes=c("A","H","B"))',
        '',
        'summary(cross)',
        '',
        '# 去除缺失基因型的个体|Remove individuals with excessive missing data',
        'n_miss <- nmissing(cross, "ind")',
        'keep <- n_miss < 0.5 * sum(nmar(cross))',
        'cross <- subset(cross, ind=keep)',
        'cat("去除高缺失个体后|After removing high-missing individuals:", nind(cross), "individuals\\n")',
        '',
    ])

    if config.map_mode == "mstmap":
        # CIM #1 with physical positions
        lines.extend(_build_cim_block(config, "physical", output_dir))

        # MSTmap: use pre-generated map files (built by Python build_mstmap_linkage_map)
        mstmap_subdir = os.path.join(output_dir, 'mstmap')
        mstmap_map_file = os.path.join(mstmap_subdir, 'mstmap_map.csv')
        mstmap_gen_file = os.path.join(mstmap_subdir, 'mstmap_gen.csv')
        mstmap_phe_file = os.path.join(mstmap_subdir, 'mstmap_phe.csv')

        lines.extend([
            'cat("\\n========== MSTmap遗传图谱 (pre-built) ==========\\n")',
            '',
            '# 用Python/MSTmap预生成的图谱文件构建csvs genfile|Build csvs genfile from MSTmap results',
            'mstmap_gen_raw <- read.csv("' + mstmap_gen_file + '", row.names = 1, check.names = FALSE)',
            'mstmap_map_df <- read.csv("' + mstmap_map_file + '", header = TRUE)',
            'mstmap_markers <- colnames(mstmap_gen_raw)',
            'am_idx <- match(mstmap_markers, as.character(mstmap_map_df[,1]))',
            'am_chr <- as.character(mstmap_map_df[am_idx, 2])',
            'am_pos <- as.character(mstmap_map_df[am_idx, 3])',
            'am_geno_t <- as.matrix(mstmap_gen_raw)',
            'am_geno_t <- cbind(rownames(am_geno_t), am_geno_t)',
            'am_chr_row <- c("", am_chr)',
            'am_pos_row <- c("", am_pos)',
            'am_df <- rbind(am_chr_row, am_pos_row, am_geno_t)',
            'colnames(am_df) <- c("", mstmap_markers)',
            'tmp_mstmap_gen <- tempfile(fileext = ".csv")',
            'write.table(am_df, file = tmp_mstmap_gen, sep = ",", row.names = FALSE,',
            '            col.names = TRUE, quote = FALSE)',
            'cross <- read.cross(format = "csvs",',
            '                    genfile = tmp_mstmap_gen,',
            '                    phefile = "' + mstmap_phe_file + '",',
            '                    genotypes = c("A", "H", "B"))',
            'n_miss <- nmissing(cross, "ind")',
            'keep <- n_miss < 0.5 * sum(nmar(cross))',
            'cross <- subset(cross, ind = keep)',
            'cat("MSTmap cross ready:", nmar(cross), "markers,", nind(cross), "individuals\\n")',
            '',
        ])

        # CIM #2 with MSTmap map
        # MSTmap块LOD系统性偏高(遗传图cM膨胀所致), 与physical块不可比量级;
        # 只用于峰物理位置互证, 最终结果以physical块为准。
        # |The mstmap block's LODs are systematically elevated (inflated cM);
        # they are not comparable in magnitude with the physical block. Use the
        # mstmap block only to cross-check peak physical locations; the physical
        # block is the final reference.
        lines.extend([
            'cat("NOTE: mstmap块LOD系统性偏高(遗传图cM膨胀), 与physical块不可直接比较量级;',
            '以physical块为最终结果, mstmap块仅用于峰物理位置互证|',
            'mstmap block LODs are systematically elevated (inflated cM), not directly',
            'comparable in magnitude with the physical block; use the physical block as',
            'the final result and the mstmap block only for cross-checking peak locations\\n")',
            '',
        ])
        lines.extend(_build_cim_block(config, "mstmap", output_dir))

        # 峰值物理坐标标注|Annotate MSTmap peaks with physical positions
        plots_dir = os.path.join(output_dir, 'plots')
        marker_index_file = os.path.join(mstmap_subdir, 'marker_map_index.tsv')
        peaks_mstmap_tsv = os.path.join(plots_dir, 'cim_peaks_mstmap.tsv')
        peaks_mstmap_phys_tsv = os.path.join(plots_dir, 'cim_peaks_mstmap_physical.tsv')
        lines.extend([
            '# MSTmap峰值物理坐标标注|Annotate MSTmap peaks with physical positions',
            'if (file.exists("' + marker_index_file + '") && file.exists("' + peaks_mstmap_tsv + '")) {',
            '    marker_idx <- read.delim("' + marker_index_file + '", stringsAsFactors = FALSE)',
            '    peaks_mst <- read.delim("' + peaks_mstmap_tsv + '", stringsAsFactors = FALSE)',
            '    if (nrow(peaks_mst) > 0) {',
            '        # 找每个峰值所在cM位置最近的marker|Find nearest marker for each peak',
            '        annotate_peak <- function(chr_lg, pos_cm, marker_idx) {',
            '            sub_idx <- marker_idx[marker_idx$LG == chr_lg, ]',
            '            if (nrow(sub_idx) == 0) return(data.frame(chr_bp=NA, pos_bp=NA))',
            '            nearest_idx <- which.min(abs(sub_idx$pos_cM - pos_cm))',
            '            data.frame(chr_bp=sub_idx$chr_bp[nearest_idx], pos_bp=sub_idx$pos_bp[nearest_idx])',
            '        }',
            '        phys <- do.call(rbind, lapply(1:nrow(peaks_mst), function(i) {',
            '            annotate_peak(as.character(peaks_mst$chr[i]), peaks_mst$pos[i], marker_idx)',
            '        }))',
            '        peaks_annotated <- cbind(peaks_mst, chr_bp=phys$chr_bp, pos_bp=phys$pos_bp)',
            '        write.table(peaks_annotated, file="' + peaks_mstmap_phys_tsv + '",',
            '                    sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
            '        cat("MSTmap峰值物理坐标已标注|Peaks annotated with physical positions\\n")',
            '    }',
            '}',
            '',
        ])

        # LOD扫描数据物理坐标标注|Annotate LOD scan data with physical positions
        # LOD行是稠密cM网格,最近标记查找会产生阶梯状坐标,故用两侧标记线性插值
        # (approx);chr_bp取最近标记的染色体(与峰值标注规则一致)。单标记LG直接取
        # 该标记坐标。写回原文件追加chr_bp/pos_bp两列,前4列保持不变。
        # LOD rows form a dense cM grid — nearest-marker lookup would produce
        # stair-step coordinates, so interpolate linearly between flanking markers
        # (approx); chr_bp follows the nearest marker (same rule as peaks). A
        # singleton LG reuses its only marker. Rewrite in place appending chr_bp/
        # pos_bp columns; the first columns stay untouched.
        lod_mstmap_tsv = os.path.join(plots_dir, 'cim_lod_data_mstmap.tsv')
        lines.extend([
            '# LOD数据物理坐标标注|Annotate LOD data with physical positions',
            'lod_mstmap_tsv <- "' + lod_mstmap_tsv + '"',
            'if (file.exists("' + marker_index_file + '") && file.exists(lod_mstmap_tsv)) {',
            '    lod_mst <- read.delim(lod_mstmap_tsv, stringsAsFactors = FALSE)',
            '    marker_idx <- read.delim("' + marker_index_file + '", stringsAsFactors = FALSE)',
            '    if (nrow(lod_mst) > 0) {',
            # 幂等守卫: 文件已含chr_bp列(如人工重跑R脚本)则跳过, 防止重复追加重名列|
            # Idempotency guard: skip if chr_bp already present (e.g. manual re-run),
            # so columns are never appended twice
            '        if (!("chr_bp" %in% colnames(lod_mst))) {',
            '            annotate_lod <- function(lg, pos_cm_vec, marker_idx) {',
            '                sub_idx <- marker_idx[marker_idx$LG == lg & !is.na(marker_idx$pos_cM) & !is.na(marker_idx$pos_bp), ]',
            '                sub_idx <- sub_idx[order(sub_idx$pos_cM), ]',
            '                sub_idx <- sub_idx[!duplicated(sub_idx$pos_cM), ]',
            '                if (nrow(sub_idx) == 0) {',
            '                    cat("  警告: LG不在marker索引中, 填NA|Warning: LG missing from marker index:", lg, "\\n")',
            '                    return(data.frame(chr_bp = rep(NA_character_, length(pos_cm_vec)),',
            '                                       pos_bp = rep(NA_real_, length(pos_cm_vec))))',
            '                }',
            '                if (nrow(sub_idx) == 1) {',
            '                    return(data.frame(chr_bp = rep(sub_idx$chr_bp[1], length(pos_cm_vec)),',
            '                                       pos_bp = rep(sub_idx$pos_bp[1], length(pos_cm_vec))))',
            '                }',
            '                pos_bp <- approx(x = sub_idx$pos_cM, y = sub_idx$pos_bp, xout = pos_cm_vec, rule = 2)$y',
            '                chr_bp <- sub_idx$chr_bp[sapply(pos_cm_vec, function(p) which.min(abs(sub_idx$pos_cM - p)))]',
            '                data.frame(chr_bp = chr_bp, pos_bp = pos_bp)',
            '            }',
            '            ann_list <- lapply(unique(lod_mst$chr), function(lg) {',
            '                rows <- which(lod_mst$chr == lg)',
            '                annotate_lod(lg, lod_mst$pos[rows], marker_idx)',
            '            })',
            '            lod_annotated <- cbind(lod_mst, do.call(rbind, ann_list))',
            '            write.table(lod_annotated, file = lod_mstmap_tsv,',
            '                        sep="\\t", row.names=FALSE, col.names=TRUE, quote=FALSE)',
            '            cat("MSTmap LOD数据物理坐标已标注|LOD data annotated with physical positions\\n")',
            '        } else {',
            '            cat("MSTmap LOD数据已含物理坐标, 跳过标注|LOD data already annotated, skipping\\n")',
            '        }',
            '    }',
            '}',
            '',
        ])

    elif config.map_mode == "estimate":
        lines.extend([
            '# 估算遗传图谱|Estimate genetic map',
            '# est.map()返回map(list),不是cross对象;必须replace.map()写回cross,',
            '# Est.map() returns a map (list), NOT a cross; replace.map() writes it back',
            '# 否则后续calc.genoprob/cim会因cross被覆盖成list而崩溃|',
            '# else calc.genoprob/cim crash because cross got overwritten by a list',
            'cat("正在估算遗传图谱...\\n")',
            'newmap <- est.map(cross, error.prob=0.0001)',
            'cross <- replace.map(cross, newmap)',
            'cat("遗传图谱估算完成\\n")',
            '',
        ])
        lines.extend(_build_cim_block(config, "", output_dir))

    else:  # physical
        lines.extend(_build_cim_block(config, "", output_dir))

    lines.append('cat("分析完成!|Analysis complete!\\n")')

    r_content = '\n'.join(lines) + '\n'

    with open(r_script, 'w') as f:
        f.write(r_content)

    return r_script


def run_r_script(script_path: str, r_env: str,
                 cmd_runner: CommandRunner, logger: logging.Logger,
                 timeout: Optional[int] = None) -> bool:
    """
    执行R脚本|Execute R script

    Args:
        script_path: R脚本路径|R script path
        r_env: conda环境名|Conda environment name
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger
        timeout: 超时时间(秒)|Timeout in seconds

    Returns:
        bool: 是否成功|Whether successful
    """
    logger.info(f"步骤5: 执行CIM R脚本|Step 5: Running CIM R script")

    success, stdout, stderr = cmd_runner.run_conda(
        r_env,
        f"Rscript {script_path}",
        "R/qtl CIM分析|R/qtl CIM analysis",
        timeout=timeout
    )

    if success:
        logger.info("CIM R脚本执行成功|CIM R script executed successfully")
        # 输出R脚本关键信息|R script key output
        for line in stdout.strip().split('\n')[-30:]:
            logger.info(f"  [R] {line}")
    else:
        logger.error("CIM R脚本执行失败|CIM R script failed")

    return success
