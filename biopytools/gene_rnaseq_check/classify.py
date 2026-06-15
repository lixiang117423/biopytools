"""候选基因分类与报告生成模块|Candidate Gene Classification and Report Generation Module"""

import os
import logging
from dataclasses import dataclass, field
from typing import Dict, List

from .config import GeneRnaseqCheckConfig
from .coverage import GeneCoverage
from .junction import GeneJunctions
from .parse_gff import EffectorRecord


# 分类标签常量|Classification tag constants
EXPRESSED_COMPLETE = 'EXPRESSED_COMPLETE'
EXPRESSED_PARTIAL = 'EXPRESSED_PARTIAL'
BOUNDARY_ISSUE = 'BOUNDARY_ISSUE'
NOT_EXPRESSED = 'NOT_EXPRESSED'
AMBIGUOUS = 'AMBIGUOUS'

# 优先级顺序（从高到低）|Priority order (high to low)
CATEGORY_PRIORITY = [
    EXPRESSED_COMPLETE, EXPRESSED_PARTIAL,
    BOUNDARY_ISSUE, NOT_EXPRESSED, AMBIGUOUS,
]


@dataclass
class ClassificationResult:
    """分类结果|Classification result"""
    gene_id: str
    chrom: str = ''
    start: int = 0
    end: int = 0
    strand: str = ''
    gene_length: int = 0
    num_exons: int = 0
    overall_coverage_pct: float = 0.0
    mean_exon_coverage_depth: float = 0.0
    min_exon_coverage_pct: float = 0.0
    upstream_500bp_mean_depth: float = 0.0
    downstream_500bp_mean_depth: float = 0.0
    junction_support: str = 'NA'
    novel_junctions: str = 'NA'
    stringtie_class_code: str = 'u'
    classification: str = ''
    needs_review: bool = False
    notes: str = ''


def classify_effector(
    gene_id: str,
    coverage: GeneCoverage,
    junctions: GeneJunctions,
    class_code: str,
    config: GeneRnaseqCheckConfig,
) -> ClassificationResult:
    """对单个基因进行分类|Classify a single gene

    Tag-based分类：所有满足条件的标签都加上，primary category按优先级取最高
    """
    tags = []

    overall_cov = coverage.overall_breadth * 100
    min_exon_cov = coverage.min_exon_breadth * 100
    gene_mean = coverage.overall_mean_depth

    upstream_ratio = (
        coverage.upstream_mean_depth / gene_mean
        if gene_mean > 0 else 0.0
    )
    downstream_ratio = (
        coverage.downstream_mean_depth / gene_mean
        if gene_mean > 0 else 0.0
    )

    novel_junc_count = (
        junctions.novel_junction_count
        if junctions.is_multi_exon else 0
    )

    # EXPRESSED_COMPLETE
    if (overall_cov >= config.expressed_complete_cov_threshold and
            min_exon_cov >= config.expressed_complete_exon_threshold and
            class_code in ('=', 'c')):
        tags.append(EXPRESSED_COMPLETE)

    # EXPRESSED_PARTIAL
    any_low_exon = any(e.breadth < config.expressed_partial_exon_threshold / 100.0
                       for e in coverage.exons)
    if (config.expressed_partial_low <= overall_cov < config.expressed_partial_high or
            any_low_exon):
        tags.append(EXPRESSED_PARTIAL)

    # BOUNDARY_ISSUE
    if overall_cov >= config.boundary_issue_min_cov:
        upstream_flag = (
            upstream_ratio > config.boundary_issue_flank_ratio and
            coverage.upstream_mean_depth > config.boundary_issue_flank_min_depth
        )
        downstream_flag = (
            downstream_ratio > config.boundary_issue_flank_ratio and
            coverage.downstream_mean_depth > config.boundary_issue_flank_min_depth
        )
        if upstream_flag or downstream_flag or novel_junc_count > 0:
            tags.append(BOUNDARY_ISSUE)

    # NOT_EXPRESSED
    if (overall_cov < config.not_expressed_max_cov and
            coverage.upstream_mean_depth <= config.not_expressed_flank_max and
            coverage.downstream_mean_depth <= config.not_expressed_flank_max):
        tags.append(NOT_EXPRESSED)

    # AMBIGUOUS: 以上都不满足
    if not tags:
        tags.append(AMBIGUOUS)

    # 确定primary category（优先级最高）|Determine primary category
    category = AMBIGUOUS
    for cat in CATEGORY_PRIORITY:
        if cat in tags:
            category = cat
            break

    # needs_review
    needs_review = category in (EXPRESSED_PARTIAL, BOUNDARY_ISSUE, AMBIGUOUS)

    # Junction字段
    if junctions.is_multi_exon:
        j_support = str(junctions.matched_junction_count)
        n_junc = str(novel_junc_count)
    else:
        j_support = 'NA'
        n_junc = 'NA'

    # Notes
    notes_parts = []
    if EXPRESSED_COMPLETE in tags:
        notes_parts.append('coverage and class_code support expression')
    if EXPRESSED_PARTIAL in tags:
        if any_low_exon:
            notes_parts.append('some exons with low coverage')
        else:
            notes_parts.append('partial coverage')
    if BOUNDARY_ISSUE in tags:
        reasons = []
        if upstream_flag:
            reasons.append('upstream coverage extends')
        if downstream_flag:
            reasons.append('downstream coverage extends')
        if novel_junc_count > 0:
            reasons.append(f'{novel_junc_count} novel junction(s)')
        notes_parts.append('boundary issue: ' + ', '.join(reasons))
    if NOT_EXPRESSED in tags:
        notes_parts.append('no expression evidence')
    if AMBIGUOUS in tags:
        notes_parts.append('insufficient evidence for clear classification')

    return ClassificationResult(
        gene_id=gene_id,
        overall_coverage_pct=round(overall_cov, 2),
        mean_exon_coverage_depth=round(
            sum(e.mean_depth for e in coverage.exons) / max(len(coverage.exons), 1), 2
        ),
        min_exon_coverage_pct=round(min_exon_cov, 2),
        upstream_500bp_mean_depth=round(coverage.upstream_mean_depth, 2),
        downstream_500bp_mean_depth=round(coverage.downstream_mean_depth, 2),
        junction_support=j_support,
        novel_junctions=n_junc,
        stringtie_class_code=class_code,
        classification=','.join(tags),
        needs_review=needs_review,
        notes='; '.join(notes_parts),
    )


def classify_all_genes(
    records: Dict[str, EffectorRecord],
    coverage_data: Dict[str, GeneCoverage],
    junction_data: Dict[str, GeneJunctions],
    class_codes: Dict[str, str],
    config: GeneRnaseqCheckConfig,
    logger: logging.Logger,
) -> List[ClassificationResult]:
    """分类所有目标基因|Classify all target genes"""
    results = []

    for gene_id, rec in records.items():
        cov = coverage_data.get(gene_id)
        junc = junction_data.get(gene_id)
        cc = class_codes.get(gene_id, 'u')

        if cov is None:
            result = ClassificationResult(
                gene_id=gene_id, classification=NOT_EXPRESSED,
                notes='no coverage data',
            )
        else:
            result = classify_effector(gene_id, cov, junc, cc, config)

        # 填充注释信息|Fill annotation info
        result.chrom = rec.chrom
        result.start = rec.gene_start
        result.end = rec.gene_end
        result.strand = rec.strand
        result.gene_length = rec.gene_length
        result.num_exons = rec.num_exons

        results.append(result)

    # 统计|Statistics
    from collections import Counter
    cat_counts = Counter(r.classification.split(',')[0] for r in results)
    review_count = sum(1 for r in results if r.needs_review)
    logger.info(f"分类统计|Classification summary:")
    for cat in CATEGORY_PRIORITY:
        logger.info(f"  {cat}: {cat_counts.get(cat, 0)}")
    logger.info(f"  需人工复核|needs_review: {review_count}")

    return results


def write_validation_report(results: List[ClassificationResult], output_path: str):
    """写入主结果表格|Write main validation report"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    headers = [
        'gene_id', 'chrom', 'start', 'end', 'strand', 'gene_length', 'num_exons',
        'overall_coverage_pct', 'mean_exon_coverage_depth', 'min_exon_coverage_pct',
        'upstream_500bp_mean_depth', 'downstream_500bp_mean_depth',
        'junction_support', 'novel_junctions', 'stringtie_class_code',
        'classification', 'needs_review', 'notes',
    ]

    with open(output_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for r in results:
            row = [
                r.gene_id, r.chrom, str(r.start), str(r.end), r.strand,
                str(r.gene_length), str(r.num_exons),
                f"{r.overall_coverage_pct:.2f}", f"{r.mean_exon_coverage_depth:.2f}",
                f"{r.min_exon_coverage_pct:.2f}",
                f"{r.upstream_500bp_mean_depth:.2f}", f"{r.downstream_500bp_mean_depth:.2f}",
                r.junction_support, r.novel_junctions, r.stringtie_class_code,
                r.classification, str(r.needs_review), r.notes,
            ]
            f.write('\t'.join(row) + '\n')


def write_flagged_bed(
    results: List[ClassificationResult],
    output_dir: str,
):
    """生成需人工复核的BED文件|Generate BED files for genes needing review

    输出按分类类别分文件，0-based half-open BED 格式
    """
    flag_dir = os.path.join(output_dir, 'flagged_for_review')
    os.makedirs(flag_dir, exist_ok=True)

    flag_categories = {
        'EXPRESSED_PARTIAL': [],
        'BOUNDARY_ISSUE': [],
        'AMBIGUOUS': [],
    }

    for r in results:
        if not r.needs_review:
            continue
        # 取primary category
        primary = r.classification.split(',')[0]
        if primary in flag_categories:
            flag_categories[primary].append(r)

    for cat_name, gene_list in flag_categories.items():
        if not gene_list:
            continue
        bed_file = os.path.join(flag_dir, f"{cat_name}.bed")
        with open(bed_file, 'w') as f:
            for r in gene_list:
                # BED: 0-based, half-open
                f.write(f"{r.chrom}\t{r.start - 1}\t{r.end}\t{r.gene_id}\t0\t{r.strand}\n")
