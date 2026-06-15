"""候选基因GFF3解析模块|Candidate Gene GFF3 Parsing Module"""

import os
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional


@dataclass
class ExonInfo:
    """Exon信息|Exon information"""
    start: int   # GFF3 1-based
    end: int     # GFF3 1-based


@dataclass
class EffectorRecord:
    """候选基因注释记录|Candidate gene annotation record"""
    gene_id: str
    chrom: str
    gene_start: int        # GFF3 1-based
    gene_end: int
    strand: str
    num_exons: int
    exons: List[ExonInfo]
    introns: List[Tuple[int, int]]  # (donor_1based, acceptor_1based)
    is_multi_exon: bool
    rep_mrna_id: str
    gene_length: int = field(init=False)

    def __post_init__(self):
        self.gene_length = self.gene_end - self.gene_start + 1
        self.is_multi_exon = len(self.exons) > 1
        self.num_exons = len(self.exons)


def _parse_gff3_attributes(attr_string: str) -> Dict[str, str]:
    """解析GFF3属性列|Parse GFF3 attribute column (column 9)

    Format: key=value;key2=value2
    """
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs


def _sort_exons(exons: List[ExonInfo], strand: str) -> List[ExonInfo]:
    """按基因组坐标排序exon|Sort exons by genomic coordinates"""
    return sorted(exons, key=lambda e: e.start)


def _compute_introns(exons: List[ExonInfo]) -> List[Tuple[int, int]]:
    """从排序后的exon列表计算intron|Compute introns from sorted exon list

    Intron: 两个相邻exon之间的间隔（1-based坐标）
    intron_start = exon[i].end + 1
    intron_end = exon[i+1].start - 1
    """
    introns = []
    for i in range(len(exons) - 1):
        donor = exons[i].end + 1      # 1-based donor
        acceptor = exons[i + 1].start - 1  # 1-based acceptor
        if donor < acceptor:
            introns.append((donor, acceptor))
    return introns


def parse_gff3_effector_genes(
    gff3_file: str,
    target_gene_ids: set,
    logger: logging.Logger,
) -> Dict[str, EffectorRecord]:
    """从GFF3文件提取目标基因的注释信息|Extract target gene annotation from GFF3

    Args:
        gff3_file: GFF3文件路径|GFF3 file path
        target_gene_ids: 目标基因ID集合|Target gene ID set
        logger: 日志器|Logger

    Returns:
        基因ID到EffectorRecord的映射|Gene ID to EffectorRecord mapping
    """
    # 第一遍：收集所有feature|Pass 1: collect all features
    # gene_id -> {gene info}
    # mrna_id -> {gene_id, exon_list, cds_length}
    genes = {}
    mrnas = {}

    logger.info(f"开始解析GFF3文件|Parsing GFF3 file: {gff3_file}")

    with open(gff3_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) < 9:
                continue

            seqid = cols[0]
            feat_type = cols[2]
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            strand = cols[6]
            attrs = _parse_gff3_attributes(cols[8])

            feat_id = attrs.get('ID', '')
            parent = attrs.get('Parent', '')

            if feat_type == 'gene':
                genes[feat_id] = {
                    'chrom': seqid, 'start': start, 'end': end,
                    'strand': strand, 'mrna_ids': [],
                }
                # 兼容不同GFF3的ID格式|Handle different ID formats
                for gid in target_gene_ids:
                    if gid == feat_id:
                        break
                    # 可能 ID=gene:xxx 格式
                    if ':' in feat_id and feat_id.split(':', 1)[1] == gid:
                        genes[gid] = genes.pop(feat_id)
                        genes[gid]['_original_id'] = feat_id
                        break

            elif feat_type in ('mRNA', 'transcript'):
                mrnas[feat_id] = {
                    'gene_id': parent, 'exons': [], 'cds_length': 0,
                    'chrom': seqid, 'strand': strand,
                }
                if parent in genes:
                    genes[parent]['mrna_ids'].append(feat_id)

            elif feat_type == 'exon':
                if parent in mrnas:
                    mrnas[parent]['exons'].append(ExonInfo(start=start, end=end))

            elif feat_type == 'CDS':
                if parent in mrnas:
                    mrnas[parent]['cds_length'] += (end - start + 1)

    # 匹配目标基因ID|Match target gene IDs
    matched = {}
    missing = set()

    for gid in target_gene_ids:
        if gid not in genes:
            missing.add(gid)
            continue

        gene_info = genes[gid]
        if not gene_info['mrna_ids']:
            logger.warning(
                f"基因 {gid} 无mRNA/transcript子特征，跳过|Gene {gid} has no mRNA/transcript, skipping"
            )
            missing.add(gid)
            continue

        # 选择代表isoform：exon数最多，tie-break为CDS最长
        best_mrna = None
        best_exon_count = -1
        best_cds_length = -1

        for mrna_id in gene_info['mrna_ids']:
            if mrna_id not in mrnas:
                continue
            mrna = mrnas[mrna_id]
            exon_count = len(mrna['exons'])
            cds_length = mrna['cds_length']
            if (exon_count > best_exon_count or
                    (exon_count == best_exon_count and cds_length > best_cds_length)):
                best_exon_count = exon_count
                best_cds_length = cds_length
                best_mrna = mrna

        if best_mrna is None or not best_mrna['exons']:
            logger.warning(
                f"基因 {gid} 无有效exon，跳过|Gene {gid} has no valid exons, skipping"
            )
            missing.add(gid)
            continue

        exons = _sort_exons(best_mrna['exons'], gene_info['strand'])
        introns = _compute_introns(exons)

        record = EffectorRecord(
            gene_id=gid,
            chrom=gene_info['chrom'],
            gene_start=gene_info['start'],
            gene_end=gene_info['end'],
            strand=gene_info['strand'],
            num_exons=len(exons),
            exons=exons,
            introns=introns,
            is_multi_exon=len(exons) > 1,
            rep_mrna_id=best_mrna,
        )
        matched[gid] = record

    logger.info(
        f"匹配到 {len(matched)} 个目标基因|Matched {len(matched)} target genes, "
        f"缺失 {len(missing)} 个|missing {len(missing)}"
    )
    if missing:
        for gid in sorted(missing):
            logger.warning(f"未在GFF3中找到|Not found in GFF3: {gid}")

    return matched


def build_exon_bed(records: Dict[str, EffectorRecord], output_path: str):
    """生成exon BED文件（0-based half-open）|Generate exon BED file (0-based half-open)

    每行格式: chrom  start-1  end  gene_id  exon_index  strand
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        for gid, rec in records.items():
            for idx, exon in enumerate(rec.exons):
                f.write(f"{rec.chrom}\t{exon.start - 1}\t{exon.end}\t{gid}\t{idx}\t{rec.strand}\n")


def build_gene_bed(records: Dict[str, EffectorRecord], output_path: str):
    """生成gene区域BED文件（0-based half-open）|Generate gene region BED file"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        for gid, rec in records.items():
            f.write(f"{rec.chrom}\t{rec.gene_start - 1}\t{rec.gene_end}\t{gid}\t0\t{rec.strand}\n")


def build_flanking_bed(
    records: Dict[str, EffectorRecord],
    output_path: str,
    window: int,
    chrom_lengths: Dict[str, int],
    direction: str = 'both',
):
    """生成上下游窗口BED文件|Generate upstream/downstream window BED file

    Args:
        direction: 'upstream', 'downstream', or 'both'
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    regions = []

    for gid, rec in records.items():
        chrom_len = chrom_lengths.get(rec.chrom, float('inf'))

        if direction in ('upstream', 'both'):
            if rec.strand == '+':
                up_end = rec.gene_start - 1  # 0-based end = 1-based start - 1
                up_start = max(0, up_end - window)
            else:
                up_start = rec.gene_end  # 0-based start = 1-based end
                up_end = min(chrom_len, up_start + window)
            if up_start < up_end:
                regions.append(f"{rec.chrom}\t{up_start}\t{up_end}\t{gid}\t0\t{rec.strand}")

        if direction in ('downstream', 'both'):
            if rec.strand == '+':
                dn_start = rec.gene_end  # 0-based start = 1-based end
                dn_end = min(chrom_len, dn_start + window)
            else:
                dn_end = rec.gene_start - 1
                dn_start = max(0, dn_end - window)
            if dn_start < dn_end:
                regions.append(f"{rec.chrom}\t{dn_start}\t{dn_end}\t{gid}\t0\t{rec.strand}")

    with open(output_path, 'w') as f:
        for line in regions:
            f.write(line + '\n')


def gff3_line_to_gtf(line: str) -> Optional[str]:
    """将单行GFF3转为GTF格式|Convert a single GFF3 line to GTF format

    用于 HISAT2 的 extract_splice_sites.py / extract_exons.py
    """
    cols = line.strip().split('\t')
    if len(cols) < 9:
        return None

    attrs = _parse_gff3_attributes(cols[8])
    gene_id = attrs.get('gene_id', attrs.get('Parent', attrs.get('ID', '')))
    transcript_id = attrs.get('transcript_id', attrs.get('ID', ''))

    if cols[2] in ('gene',):
        gene_id = attrs.get('ID', '')
        transcript_id = gene_id
    elif cols[2] in ('mRNA', 'transcript'):
        gene_id = attrs.get('Parent', '')
        transcript_id = attrs.get('ID', '')
    elif cols[2] == 'exon':
        gene_id = attrs.get('gene_id', attrs.get('Parent', ''))
        transcript_id = attrs.get('transcript_id', attrs.get('Parent', ''))

    if not gene_id:
        return None

    gt_attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"'
    return '\t'.join(cols[:8] + [gt_attrs])


def convert_gff3_to_gtf(gff3_file: str, gtf_file: str, logger: logging.Logger):
    """将GFF3文件转为GTF格式|Convert GFF3 file to GTF format

    仅转换 gene/mRNA/transcript/exon/CDS 行
    """
    os.makedirs(os.path.dirname(gtf_file), exist_ok=True)
    convert_types = {'gene', 'mRNA', 'transcript', 'exon', 'CDS'}
    count = 0

    with open(gff3_file) as fin, open(gtf_file, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) >= 3 and cols[2] in convert_types:
                gtf_line = gff3_line_to_gtf(line)
                if gtf_line:
                    fout.write(gtf_line + '\n')
                    count += 1

    logger.info(f"转换了 {count} 行GFF3到GTF|Converted {count} GFF3 lines to GTF: {gtf_file}")
    return gtf_file
