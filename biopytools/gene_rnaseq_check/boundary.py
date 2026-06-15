"""候选基因注释边界建议模块|Candidate Gene Annotation Boundary Suggestion Module"""

import os
import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

from .config import GeneRnaseqCheckConfig
from .coverage import GeneCoverage
from .junction import GeneJunctions, JunctionInfo
from .parse_gff import EffectorRecord
from .classify import ClassificationResult, BOUNDARY_ISSUE


@dataclass
class BoundarySuggestion:
    """边界修改建议|Boundary modification suggestion"""
    gene_id: str
    original_start: int
    original_end: int
    suggested_start: int
    suggested_end: int
    evidence: str     # 'upstream_coverage' / 'downstream_coverage' / 'novel_junction'
    confidence: str  # 'high' / 'medium' / 'low'
    notes: str


def _parse_gtf_attributes(attr_string: str) -> dict:
    """解析GTF属性列|Parse GTF attribute column"""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key.strip()] = val.strip().strip('"')
        else:
            parts = item.split(None, 1)
            if len(parts) == 2:
                attrs[parts[0].strip()] = parts[1].strip().strip('"')
    return attrs


def _parse_stringtie_gtf(gtf_file: str, chrom: str, gene_start: int,
                         gene_end: int) -> List[Dict]:
    """从StringTie GTF提取覆盖指定区域的转录本|Extract transcripts overlapping region from StringTie GTF

    Returns:
        转录本信息列表，每个包含 start, end, strand, transcript_id
    """
    transcripts = []
    if not os.path.exists(gtf_file):
        return transcripts

    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            if cols[2] != 'transcript':
                continue
            t_chrom = cols[0]
            t_start = int(cols[3])
            t_end = int(cols[4])
            t_strand = cols[6]

            # 检查是否重叠|Check overlap
            if t_chrom != chrom:
                continue
            if t_start > gene_end or t_end < gene_start:
                continue

            attrs = _parse_gtf_attributes(cols[8])
            transcripts.append({
                'start': t_start, 'end': t_end,
                'strand': t_strand,
                'transcript_id': attrs.get('transcript_id', ''),
            })

    return transcripts


class BoundaryAnalyzer:
    """边界分析器|Boundary Analyzer"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def suggest_boundaries(
        self,
        classifications: List[ClassificationResult],
        records: Dict[str, EffectorRecord],
        coverage_data: Dict[str, GeneCoverage],
        junction_data: Dict[str, GeneJunctions],
        stringtie_gtf: str,
    ) -> List[BoundarySuggestion]:
        """为BOUNDARY_ISSUE类基因生成边界建议|Generate boundary suggestions for BOUNDARY_ISSUE genes

        Args:
            classifications: 分类结果列表|Classification results
            records: 基因注释记录|Gene annotation records
            coverage_data: 覆盖度数据|Coverage data
            junction_data: Junction数据|Junction data
            stringtie_gtf: StringTie输出的GTF文件|StringTie GTF output

        Returns:
            边界建议列表|Boundary suggestion list
        """
        boundary_genes = [
            r for r in classifications
            if BOUNDARY_ISSUE in r.classification.split(',')
        ]

        if not boundary_genes:
            self.logger.info("无需边界分析|No boundary issues to analyze")
            return []

        self.logger.info(f"开始边界分析，共 {len(boundary_genes)} 个基因|"
                         f"Starting boundary analysis for {len(boundary_genes)} genes")

        suggestions = []

        for result in boundary_genes:
            gid = result.gene_id
            rec = records.get(gid)
            cov = coverage_data.get(gid)
            junc = junction_data.get(gid)

            if not rec or not cov:
                continue

            suggestion = self._analyze_gene(gid, rec, cov, junc, stringtie_gtf)
            if suggestion:
                suggestions.append(suggestion)

        return suggestions

    def _analyze_gene(
        self,
        gene_id: str,
        rec: EffectorRecord,
        cov: GeneCoverage,
        junc: Optional[GeneJunctions],
        stringtie_gtf: str,
    ) -> Optional[BoundarySuggestion]:
        """分析单个基因的边界|Analyze boundary for a single gene"""
        gene_mean = cov.overall_mean_depth

        new_start = rec.gene_start
        new_end = rec.gene_end
        evidence_parts = []
        confidence_scores = []

        # 1. 上游延伸分析|Upstream extension analysis
        if cov.upstream_mean_depth > self.config.boundary_issue_flank_min_depth:
            upstream_ratio = cov.upstream_mean_depth / gene_mean if gene_mean > 0 else 0
            if upstream_ratio > self.config.boundary_issue_flank_ratio:
                # 检查StringTie是否有更长的转录本|Check StringTie for longer transcript
                transcripts = _parse_stringtie_gtf(
                    stringtie_gtf, rec.chrom,
                    rec.gene_start - self.config.flanking_window,
                    rec.gene_end,
                )
                min_start = rec.gene_start
                for t in transcripts:
                    if t['start'] < min_start:
                        min_start = t['start']

                if min_start < rec.gene_start:
                    # StringTie转录本支持延伸|StringTie supports extension
                    new_start = min_start
                    evidence_parts.append('upstream_coverage')
                    confidence_scores.append(2)  # high
                else:
                    # 有覆盖延伸但StringTie不支持|Coverage extends but no StringTie support
                    evidence_parts.append('upstream_coverage')
                    confidence_scores.append(1)  # medium

        # 2. 下游延伸分析|Downstream extension analysis
        if cov.downstream_mean_depth > self.config.boundary_issue_flank_min_depth:
            downstream_ratio = cov.downstream_mean_depth / gene_mean if gene_mean > 0 else 0
            if downstream_ratio > self.config.boundary_issue_flank_ratio:
                transcripts = _parse_stringtie_gtf(
                    stringtie_gtf, rec.chrom,
                    rec.gene_start,
                    rec.gene_end + self.config.flanking_window,
                )
                max_end = rec.gene_end
                for t in transcripts:
                    if t['end'] > max_end:
                        max_end = t['end']

                if max_end > rec.gene_end:
                    new_end = max_end
                    evidence_parts.append('downstream_coverage')
                    confidence_scores.append(2)
                else:
                    evidence_parts.append('downstream_coverage')
                    confidence_scores.append(1)

        # 3. Novel junction分析|Novel junction analysis
        if junc and junc.is_multi_exon and junc.novel_junctions:
            # 检查novel junction是否暗示边界扩展|Check if novel junctions imply boundary extension
            extends_start = False
            extends_end = False

            for jn in junc.novel_junctions:
                # 如果novel junction的donor在当前gene start之前
                # 或acceptor在当前gene end之后
                jn_donor_1 = jn.donor + 1  # 转为1-based
                jn_acceptor_1 = jn.acceptor  # 0-based exclusive = 1-based inclusive
                if jn_donor_1 < rec.gene_start:
                    extends_start = True
                if jn_acceptor_1 > rec.gene_end:
                    extends_end = True

            if extends_start or extends_end:
                evidence_parts.append('novel_junction')
                confidence_scores.append(1)

        # 汇总|Summarize
        if not evidence_parts:
            return None

        evidence = ','.join(evidence_parts)

        # 置信度|Confidence
        if any(s >= 2 for s in confidence_scores):
            confidence = 'high'
        elif any(s >= 1 for s in confidence_scores):
            confidence = 'medium'
        else:
            confidence = 'low'

        # Notes
        notes_parts = []
        if new_start < rec.gene_start:
            notes_parts.append(f"建议上游扩展{rec.gene_start - new_start}bp|"
                               f"suggest upstream extension {rec.gene_start - new_start}bp")
        if new_end > rec.gene_end:
            notes_parts.append(f"建议下游扩展{new_end - rec.gene_end}bp|"
                               f"suggest downstream extension {new_end - rec.gene_end}bp")
        if junc and junc.novel_junctions:
            notes_parts.append(f"{junc.novel_junction_count}个novel junction|"
                               f"{junc.novel_junction_count} novel junction(s)")
        notes = '; '.join(notes_parts)

        return BoundarySuggestion(
            gene_id=gene_id,
            original_start=rec.gene_start,
            original_end=rec.gene_end,
            suggested_start=new_start,
            suggested_end=new_end,
            evidence=evidence,
            confidence=confidence,
            notes=notes,
        )


def write_boundary_suggestions(suggestions: List[BoundarySuggestion], output_path: str):
    """写入边界建议表格|Write boundary suggestions table"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    headers = [
        'gene_id', 'original_start', 'original_end',
        'suggested_start', 'suggested_end',
        'evidence', 'confidence', 'notes',
    ]

    with open(output_path, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for s in suggestions:
            row = [
                s.gene_id, str(s.original_start), str(s.original_end),
                str(s.suggested_start), str(s.suggested_end),
                s.evidence, s.confidence, s.notes,
            ]
            f.write('\t'.join(row) + '\n')
