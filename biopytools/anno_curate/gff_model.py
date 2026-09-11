"""GXF数据模型与解析|GXF data model and parsing

GFF3/GTF 单遍解析为内存对象;GTF 行内转 GFF3 属性
|Single-pass parse into in-memory records; GTF rows converted to GFF3 attrs
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# 转录本级类型识别(补建/悬空/排序共用,大小写不敏感)
# |Transcript-level type recognition (shared, case-insensitive)
TRANSCRIPT_RE = re.compile(
    r'mRNA|transcript\b|ncRNA|miRNA|rRNA|snoRNA|snRNA|tRNA|lnc_RNA|lncRNA',
    re.IGNORECASE)

# 显式给 rank 3 的非 mRNA 转录本级(spec §5⑪ 修正)
# |Non-mRNA transcript-level types explicitly ranked 3 (spec §5⑪)
TRANSCRIPT_LEVEL_TYPES = {
    'tRNA', 'rRNA', 'ncRNA', 'miRNA', 'snoRNA', 'snRNA',
    'lnc_RNA', 'lncRNA'}

FEATURE_RANKS = {
    'gene': 1, 'mRNA': 2, 'transcript': 3, 'exon': 4,
    'five_prime_UTR': 5, 'CDS': 6, 'three_prime_UTR': 7}

# UTR 类型识别(诊断用,小写比较)|UTR type recognition (lowercase compare)
UTR5_TYPES = {'five_prime_utr', 'five_prime_utr_region', "5'utr"}
UTR3_TYPES = {'three_prime_utr', 'three_prime_utr_region', "3'utr"}

# 大小写敏感的规范类型(变体检测)|Canonical types for case-variant detection
CANONICAL_TYPES = {'gene', 'mrna', 'transcript', 'exon', 'cds',
                   'five_prime_utr', 'three_prime_utr'}


def feature_rank(feature_type: str) -> int:
    """特征排序等级(§5⑪)|Feature sort rank (§5⑪)"""
    if feature_type in FEATURE_RANKS:
        return FEATURE_RANKS[feature_type]
    if feature_type in TRANSCRIPT_LEVEL_TYPES:
        return 3
    return 8


@dataclass
class GffRecord:
    """GFF3 单条记录|One GFF3 record"""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: Dict[str, str] = field(default_factory=dict)
    line_no: int = 0

    def to_line(self) -> str:
        """序列化为 9 列行|Serialize to a 9-column line"""
        attrs = ';'.join(f'{k}={v}' for k, v in self.attributes.items())
        return '\t'.join((self.seqid, self.source, self.type,
                          str(self.start), str(self.end), self.score,
                          self.strand, self.phase, attrs))


class GxfDocument:
    """内存 GXF 文档|In-memory GXF document"""

    def __init__(self):
        self.records: List[GffRecord] = []
        self.fasta_block: Optional[str] = None
        self.chr_order: List[str] = []
        self.case_variants: List[str] = []
        self.id_index: Dict[str, GffRecord] = {}
        self.children_by_parent: Dict[str, List[GffRecord]] = {}
        # 染色体首现去重集合(append_chr 用,须随实例初始化避免外部注入)
        # |Chromosome dedup set for append_chr (per-instance, no external
        # injection)
        self._chr_seen: set = set()

    def rebuild_indexes(self):
        """重建 ID/Parent 索引|Rebuild ID and Parent indexes"""
        self.id_index = {}
        self.children_by_parent = {}
        for r in self.records:
            rec_id = r.attributes.get('ID')
            if rec_id and rec_id not in self.id_index:
                self.id_index[rec_id] = r
            parent = r.attributes.get('Parent')
            if parent:
                # 只认第一个 Parent 值(spec §6.4)|First Parent value only
                first = parent.split(',')[0]
                self.children_by_parent.setdefault(first, []).append(r)

    def append_chr(self, seqid: str):
        """记录染色体首现序|Track chromosome first-seen order"""
        if seqid not in self._chr_seen:
            self._chr_seen.add(seqid)
            self.chr_order.append(seqid)


def _parse_attributes_gff3(attr_text: str) -> Dict[str, str]:
    """GFF3 属性列→有序dict|GFF3 attribute column → ordered dict"""
    attrs: Dict[str, str] = {}
    for part in attr_text.split(';'):
        part = part.strip()
        if not part or '=' not in part:
            continue
        key, _, val = part.partition('=')
        attrs[key.strip()] = val.strip()
    return attrs


def _parse_attributes_gtf(attr_text: str) -> Dict[str, str]:
    """GTF 属性列→GFF3 风格dict|GTF attributes → GFF3-style dict

    transcript_id→Parent(子特征);gene_id 保留;去引号;行尾无分号也收
    (仍要求 key "value" 引号形态)
    |transcript_id→Parent (child features); gene_id kept; quotes stripped;
    trailing semicolon optional (key "value" quoting still required)
    """
    attrs: Dict[str, str] = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"\s*;?', attr_text):
        attrs[m.group(1)] = m.group(2)
    transcript_id = attrs.pop('transcript_id', None)
    out: Dict[str, str] = {}
    if transcript_id:
        out['Parent'] = transcript_id
    for k, v in attrs.items():
        out[k] = v
    return out


def _looks_like_gtf(attr_text: str) -> bool:
    """GTF 识别:带引号键值或 transcript_id|GTF detection

    引号键值允许行尾无分号(与 _parse_attributes_gtf 同口径)
    |Quoted pairs allow a missing trailing semicolon (same rule as
    _parse_attributes_gtf)
    """
    return 'transcript_id' in attr_text or bool(
        re.search(r'\S+\s+"[^"]*"\s*;?', attr_text))


def parse_gxf(path: str) -> GxfDocument:
    """解析 GFF3/GTF|Parse GFF3 or GTF

    列数≠9 / 坐标非数字:收集全部问题后 ValueError(§10 硬失败)
    |Column count != 9 / non-numeric coords: collect all, raise ValueError
    """
    doc = GxfDocument()
    errors = []
    variants = {}
    with open(path, encoding='utf-8') as fh:
        in_fasta = False
        fasta_lines = []
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if in_fasta:
                fasta_lines.append(line)
                continue
            if line.startswith('##FASTA'):
                in_fasta = True
                fasta_lines.append(line)
                continue
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) != 9:
                errors.append(f"第{line_no}行列数={len(cols)}(应为9)"
                              f"|Line {line_no}: {len(cols)} columns (need 9)")
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                errors.append(f"第{line_no}行坐标非数字|Line {line_no}: "
                              f"non-numeric coords: {cols[3]},{cols[4]}")
                continue
            if start > end:
                start, end = end, start   # 坐标交换(spec §5①)
            if _looks_like_gtf(cols[8]):
                attrs = _parse_attributes_gtf(cols[8])
                # GTF 的 transcript 行自身就是转录本级:transcript_id 须转为
                # 本行 ID(有 gene_id 则作 Parent),不能按子特征规则塞进
                # Parent——否则该行变成"以自身为父的无 ID 孤儿",⑦ 会把每个
                # 可变剪接体拆成独立基因、诊断的转录本计数翻倍
                # |A GTF transcript row IS transcript-level: its
                # transcript_id becomes the row ID (gene_id its Parent);
                # the child rule would leave an ID-less orphan parenting
                # itself, so ⑦ would split isoforms into separate genes
                # and double the transcript roster
                if TRANSCRIPT_RE.search(cols[2]) and 'Parent' in attrs:
                    gene_id = attrs.get('gene_id')
                    attrs = {'ID': attrs['Parent']}
                    if gene_id:
                        attrs['Parent'] = gene_id
            else:
                attrs = _parse_attributes_gff3(cols[8])
            rec = GffRecord(cols[0], cols[1], cols[2], start, end, cols[5],
                            cols[6], cols[7], attrs, line_no)
            doc.records.append(rec)
            doc.append_chr(rec.seqid)
            low = rec.type.lower()
            if low in CANONICAL_TYPES and rec.type not in (
                    'gene', 'mRNA', 'transcript', 'exon', 'CDS',
                    'five_prime_UTR', 'three_prime_UTR'):
                variants.setdefault(low, []).append(line_no)
    if errors:
        raise ValueError('\n'.join(errors))
    for low, lines in variants.items():
        doc.case_variants.append(
            f"feature类型大小写变体|case-variant feature type: '{low}' "
            f"x{len(lines)} (行|lines: {lines[:5]}{'...' if len(lines) > 5 else ''})"
            f"——这些行不参与相位/UTR/exon逻辑"
            f"|skipped by phase/UTR/exon logic")
    if fasta_lines:
        doc.fasta_block = '\n'.join(fasta_lines) + '\n'
    doc.rebuild_indexes()
    return doc


def sort_records(records: List[GffRecord],
                 chr_order: List[str]) -> List[GffRecord]:
    """排序:染色体首现序+start+rank+line_no|Sort by chr-first-seen+start+rank"""
    rank_of = {c: i for i, c in enumerate(chr_order)}

    def key(r: GffRecord):
        return (rank_of.get(r.seqid, len(chr_order)), r.start,
                feature_rank(r.type), r.line_no)

    return sorted(records, key=key)


def write_gff3(doc: GxfDocument, records: List[GffRecord], path: str):
    """写规范 GFF3|Write canonical GFF3"""
    ordered = sort_records(records, doc.chr_order)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('##gff-version 3\n')
        for r in ordered:
            fh.write(r.to_line() + '\n')
        if doc.fasta_block:
            fh.write(doc.fasta_block)
