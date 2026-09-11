"""注释修复流水线|Annotation fix pipeline

11 步映射 GSAman GXF Fix(spec §5);步骤在内存对象上依次执行
|11 steps mirroring GSAman GXF Fix (spec §5), executed in memory
"""

import logging
import re
from dataclasses import dataclass, field
from typing import Dict, List

from .gff_model import (TRANSCRIPT_RE, GffRecord, GxfDocument, write_gff3)

# 有孩子的特征类型(悬空判定)|Child-bearing types (dangling criterion)
CHILD_TYPES = ('exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR')


def _is_transcript(rec: GffRecord) -> bool:
    """转录本级判定(大小写不敏感)|Transcript-level check (case-insensitive)"""
    return bool(TRANSCRIPT_RE.search(rec.type))


@dataclass
class GxffixResult:
    """修复结果|Fix result"""

    doc: GxfDocument
    dangling_records: List[GffRecord] = field(default_factory=list)
    corrected_transcripts: List[str] = field(default_factory=list)
    problematic_transcripts: List[str] = field(default_factory=list)
    stats: Dict[str, int] = field(default_factory=dict)


class Gxffixer:
    """修复器|Fixer"""

    def __init__(self, doc: GxfDocument, logger: logging.Logger):
        self.doc = doc
        self.logger = logger
        self.stats: Dict[str, int] = {}

    # ---------- ② ID/Parent 冲突|ID/Parent conflicts ----------

    def check_id_parent_conflicts(self):
        """环状 Parent 解开+重复转录本ID硬失败|Circular fix + dup ID fail"""
        self.doc.rebuild_indexes()
        seen: Dict[str, int] = {}
        dups: List[str] = []
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            rec_id = rec.attributes.get('ID')
            if not rec_id:
                continue
            parent = rec.attributes.get('Parent', '')
            if parent and rec_id in parent.split(','):
                # 改 mRNA 自身 Parent(spec §5②);gene 由 ⑦ 补建接住
                rec.attributes['Parent'] = f'{rec_id}.gene'
                self.stats['circular_parent_fixed'] = \
                    self.stats.get('circular_parent_fixed', 0) + 1
            seen[rec_id] = seen.get(rec_id, 0) + 1
        for rec_id, n in seen.items():
            if n > 1:
                dups.append(f"重复转录本ID|duplicate transcript ID: "
                            f"{rec_id} x{n}")
        if dups:
            raise ValueError(
                "存在重复转录本ID,须先人工解决|Duplicate transcript IDs, "
                "resolve manually first:\n" + '\n'.join(dups))

    # ---------- ③ gene 行清理|gene row cleanup ----------

    def check_usage_of_gene_feature(self):
        """无 mRNA 行时删 gene|Drop gene rows when no mRNA exists"""
        has_mrna = any(r.type == 'mRNA' for r in self.doc.records)
        if has_mrna:
            return
        before = len(self.doc.records)
        self.doc.records = [r for r in self.doc.records if r.type != 'gene']
        self.stats['gene_rows_removed'] = before - len(self.doc.records)

    # ---------- ⑨ 悬空 mRNA 拆分|dangling mRNA split ----------

    def split_dangling_mrna(self) -> List[GffRecord]:
        """隔离无子特征 mRNA(连同独占的父gene)|Isolate childless mRNAs"""
        self.doc.rebuild_indexes()
        dangling_ids = set()
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            rec_id = rec.attributes.get('ID')
            if not rec_id:
                continue    # 无 ID 转录行(GTF 畸形)不参与悬空判定
            children = self.doc.children_by_parent.get(rec_id, [])
            if not any(c.type in CHILD_TYPES or
                       c.type.lower().endswith('_utr') for c in children):
                dangling_ids.add(rec_id)
        if not dangling_ids:
            return []
        # 悬空转录本的原有子特征(非exon/CDS/UTR类型)也带走
        # |Carry the dangling transcript's odd children too
        take_ids = set(dangling_ids)
        for tid in dangling_ids:
            for child in self.doc.children_by_parent.get(tid, []):
                cid = child.attributes.get('ID')
                if cid:
                    take_ids.add(cid)
        # 预建 父ID→转录本 索引(一次 O(n)),免每个 gene 全量扫 records
        # (26k records 已 7.5s,50万级小时级);Parent 多值只认首值,
        # 与 rebuild_indexes 口径一致
        # |Pre-build parent-ID→transcripts map in one O(n) pass; the old
        # per-gene full scan was O(genes x records). Comma-joined Parent
        # values use the first entry, same rule as rebuild_indexes
        kids_by_parent: Dict[str, List[GffRecord]] = {}
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            parent = rec.attributes.get('Parent')
            if parent:
                kids_by_parent.setdefault(parent.split(',')[0], []).append(rec)
        kept = []
        dangling = []
        for rec in self.doc.records:
            rec_id = rec.attributes.get('ID', '')
            if rec_id in take_ids:
                dangling.append(rec)
                continue
            if rec.type == 'gene':
                # 独占父gene带走;仍有其他保留mRNA则留下
                # |Take the gene only if no kept transcript remains
                kids = kids_by_parent.get(rec_id, [])
                if kids and not [k for k in kids
                                 if k.attributes.get('ID') not in dangling_ids]:
                    dangling.append(rec)
                    continue
            kept.append(rec)
        self.doc.records = kept
        self.stats['dangling_mrna'] = len(dangling_ids)
        if dangling_ids:
            self.logger.warning(
                f"发现{len(dangling_ids)}个悬空mRNA,已隔离|Found "
                f"{len(dangling_ids)} dangling mRNAs, isolated")
        return dangling

    # ---------- ④ mRNA 补建|mRNA recall ----------

    def recall_mrna_features(self):
        """为有孩子没妈妈的转录本补建mRNA|Build mRNA for orphan children"""
        self.doc.rebuild_indexes()
        # Parent 缺失的子特征归到 ID+"-mRNA";gene/转录本行是根层级,
        # 写伪 Parent 在 GFF3 非法,故只认子特征类型(与 kids 过滤同谓词)
        # |Parentless child features assigned to ID+"-mRNA"; gene/transcript
        # rows are root-level (fake Parent illegal in GFF3), so only child
        # types qualify (same predicate as the kids filter)
        for rec in self.doc.records:
            if 'Parent' not in rec.attributes and rec.attributes.get('ID'):
                if rec.type in CHILD_TYPES or \
                        rec.type.lower().endswith('_utr'):
                    rec.attributes['Parent'] = \
                        f"{rec.attributes['ID']}-mRNA"
        self.doc.rebuild_indexes()
        existing = {r.attributes.get('ID') for r in self.doc.records
                    if _is_transcript(r)}
        built = []
        for parent_id, children in self.doc.children_by_parent.items():
            if parent_id in existing:
                continue
            kids = [c for c in children
                    if c.type in CHILD_TYPES or
                    c.type.lower().endswith('_utr') or c.type == 'CDS']
            if not kids:
                continue
            attrs = {}
            for c in kids:                       # 属性全量继承(去 Parent)
                for k, v in c.attributes.items():
                    if k != 'Parent':
                        attrs.setdefault(k, v)
            gene_parent = attrs.pop('gene_id', None) or \
                attrs.pop('gene_name', None)
            if gene_parent:
                attrs['Parent'] = gene_parent
            mrna = GffRecord(
                kids[0].seqid, 'Recall', 'mRNA',
                min(c.start for c in kids), max(c.end for c in kids),
                '.', kids[0].strand, '.', dict(attrs, ID=parent_id),
                line_no=kids[0].line_no)
            built.append(mrna)
        self.doc.records.extend(built)
        self.stats['mrna_recalled'] = len(built)
        if built:
            self.logger.info(
                f"补建{len(built)}个mRNA(source=Recall)|Recalled "
                f"{len(built)} mRNAs (source=Recall)")

    # ---------- ⑦ gene 补建|gene rebuild ----------

    def fix_gene_features(self):
        """为无gene的mRNA补建gene|Build gene for orphan mRNAs"""
        self.doc.rebuild_indexes()
        genes = {r.attributes.get('ID') for r in self.doc.records
                 if r.type == 'gene'}
        # 第一遍:决定每个 mRNA 的 gene_id(Parent 非空即用;否则 geneID/
        # gene_id 属性;否则 <ID>.gene)——Parent 指向不存在 gene 时不改写,
        # 交给第二遍按该 Parent 补建 gene(GSAman 回退链语义)
        # |Pass 1: resolve each mRNA's gene_id (keep non-empty Parent; else
        # geneID/gene_id attribute; else <ID>.gene). A Parent pointing to a
        # missing gene is NOT rewritten — pass 2 builds that gene instead.
        assignments = {}
        for rec in self.doc.records:
            if not _is_transcript(rec) or rec.type == 'gene':
                continue
            parent = rec.attributes.get('Parent', '')
            gene_id = parent or rec.attributes.get('geneID') or \
                rec.attributes.get('gene_id')
            if not gene_id:
                gene_id = f"{rec.attributes.get('ID', 'tx')}.gene"
            assignments[id(rec)] = gene_id
            rec.attributes['Parent'] = gene_id
        # 第二遍:新建缺失 gene,范围=该 gene 名下全部 mRNA 聚合 min/max
        # |Pass 2: build missing genes, span aggregated over member mRNAs
        by_gene = {}
        for rec in self.doc.records:
            gid = assignments.get(id(rec))
            if gid and gid not in genes:
                by_gene.setdefault(gid, []).append(rec)
        new_genes = []
        for gid, members in by_gene.items():
            new_genes.append(GffRecord(
                members[0].seqid, 'anno_curate', 'gene',
                min(m.start for m in members), max(m.end for m in members),
                '.', members[0].strand, '.', {'ID': gid},
                line_no=min(m.line_no for m in members)))
        self.doc.records.extend(new_genes)
        self.stats['genes_built'] = len(new_genes)

    # ---------- ⑤ UTR 补建|UTR rebuild ----------

    @staticmethod
    def _utr_type(is_left: bool, strand: str) -> str:
        """链方向定UTR类型|Strand decides UTR type"""
        five_is_left = strand != '-'
        if is_left:
            return 'five_prime_UTR' if five_is_left else 'three_prime_UTR'
        return 'three_prime_UTR' if five_is_left else 'five_prime_UTR'

    def fix_utr_features(self):
        """补建UTR(spec §5⑤判定顺序)|Rebuild UTRs (spec §5⑤ order)"""
        self.doc.rebuild_indexes()
        new_utrs = []
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            tid = rec.attributes.get('ID', '')
            children = self.doc.children_by_parent.get(tid, [])
            cds = [c for c in children if c.type == 'CDS']
            if not cds:
                continue
            exons = [c for c in children if c.type == 'exon']
            has_utr = any('UTR' in c.type.upper() or
                          c.type.lower().endswith('_utr') for c in children)
            cds_start = min(c.start for c in cds)
            cds_end = max(c.end for c in cds)
            strand = rec.strand
            if exons:
                # 分支一:exon−CDS 拆分|branch: exon minus CDS
                for ex in exons:
                    if ex.end < cds_start:                 # exon 全在 CDS 左
                        new_utrs.append(GffRecord(
                            ex.seqid, 'anno_curate',
                            self._utr_type(True, strand),
                            ex.start, ex.end, '.', strand, '.',
                            {'Parent': tid}, ex.line_no))
                    elif ex.start > cds_end:                # 全在右
                        new_utrs.append(GffRecord(
                            ex.seqid, 'anno_curate',
                            self._utr_type(False, strand),
                            ex.start, ex.end, '.', strand, '.',
                            {'Parent': tid}, ex.line_no))
                    else:                                   # 两端超出
                        left = ex.start < cds_start
                        right = ex.end > cds_end
                        if left:
                            new_utrs.append(GffRecord(
                                ex.seqid, 'anno_curate',
                                self._utr_type(True, strand),
                                ex.start, cds_start - 1, '.', strand, '.',
                                {'Parent': tid}, ex.line_no))
                        if right:
                            new_utrs.append(GffRecord(
                                ex.seqid, 'anno_curate',
                                self._utr_type(False, strand),
                                cds_end + 1, ex.end, '.', strand, '.',
                                {'Parent': tid}, ex.line_no))
            elif not has_utr:
                # 分支二:mRNA−CDS 差集|branch: mRNA minus CDS
                if rec.start < cds_start:
                    new_utrs.append(GffRecord(
                        rec.seqid, 'anno_curate',
                        self._utr_type(True, strand),
                        rec.start, cds_start - 1, '.', strand, '.',
                        {'Parent': tid}, rec.line_no))
                if rec.end > cds_end:
                    new_utrs.append(GffRecord(
                        rec.seqid, 'anno_curate',
                        self._utr_type(False, strand),
                        cds_end + 1, rec.end, '.', strand, '.',
                        {'Parent': tid}, rec.line_no))
        self.doc.records.extend(new_utrs)
        self.stats['utrs_built'] = len(new_utrs)

    # ---------- ⑥ exon 补建|exon rebuild ----------

    def fix_exon_features(self):
        """仅无exon行的转录本补建exon|Rebuild exons only when absent"""
        self.doc.rebuild_indexes()
        new_exons = []
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            tid = rec.attributes.get('ID', '')
            children = self.doc.children_by_parent.get(tid, [])
            if any(c.type == 'exon' for c in children):
                continue
            parts = [c for c in children
                     if c.type == 'CDS' or 'UTR' in c.type.upper() or
                     c.type.lower().endswith('_utr')]
            if not parts:
                continue
            parts.sort(key=lambda c: c.start)
            merged = [[parts[0].start, parts[0].end]]
            for c in parts[1:]:
                if merged[-1][1] + 1 >= c.start:      # gap≤1bp 合并
                    merged[-1][1] = max(merged[-1][1], c.end)
                else:
                    merged.append([c.start, c.end])
            for start, end in merged:
                new_exons.append(GffRecord(
                    rec.seqid, 'anno_curate', 'exon', start, end, '.',
                    rec.strand, '.', {'Parent': tid}, rec.line_no))
        self.doc.records.extend(new_exons)
        self.stats['exons_built'] = len(new_exons)

    # ---------- ⑧ 子特征 ID 去重|child ID dedup ----------

    def fix_duplicated_child_ids(self):
        """重复子特征ID加.u<N>|Suffix duplicated child IDs"""
        counts: Dict[str, int] = {}
        for rec in self.doc.records:
            cid = rec.attributes.get('ID')
            if cid and (rec.type in CHILD_TYPES or
                        rec.type.lower().endswith('_utr')):
                counts[cid] = counts.get(cid, 0) + 1
        uniq = 0
        seen: Dict[str, int] = {}
        for rec in self.doc.records:
            cid = rec.attributes.get('ID')
            if not cid or counts.get(cid, 0) <= 1:
                continue
            if not (rec.type in CHILD_TYPES or
                    rec.type.lower().endswith('_utr')):
                continue
            seen[cid] = seen.get(cid, 0) + 1
            if seen[cid] > 1:
                uniq += 1
                rec.attributes['ID'] = f'{cid}.u{uniq}'
        self.stats['child_ids_deduped'] = uniq

    # ---------- ⑩ 相位内嵌|inline phase correction ----------

    def run_phase_correction(self) -> tuple:
        """逐转录本相位三态处理|Per-transcript phase handling"""
        from .phase import (apply_phases, order_cds_by_strand,
                            validate_transcript_phases)
        self.doc.rebuild_indexes()
        corrected, problematic = [], []
        for rec in self.doc.records:
            if not _is_transcript(rec):
                continue
            tid = rec.attributes.get('ID', '')
            cds = [c for c in self.doc.children_by_parent.get(tid, [])
                   if c.type == 'CDS']
            if not cds:
                continue
            ordered = order_cds_by_strand(cds)
            res = validate_transcript_phases(ordered)
            if res.status == 'CORRECTED':
                apply_phases(ordered, res.expected)
                corrected.append(tid)
            elif res.status == 'INVALID_LENGTH':
                problematic.append(tid)
        self.stats['phase_corrected'] = len(corrected)
        self.stats['phase_problematic'] = len(problematic)
        return corrected, problematic

    # ---------- 总编排|orchestration ----------

    def run(self) -> GxffixResult:
        """执行全部修复步骤|Run all fix steps"""
        self.check_id_parent_conflicts()      # ②
        self.check_usage_of_gene_feature()    # ③
        self.recall_mrna_features()           # ④
        self.fix_utr_features()               # ⑤
        self.fix_exon_features()              # ⑥
        self.fix_gene_features()              # ⑦
        self.fix_duplicated_child_ids()       # ⑧
        dangling = self.split_dangling_mrna()  # ⑨
        corrected, problematic = self.run_phase_correction()  # ⑩
        # write_outputs 需要这三份产物|kept for write_outputs
        self._last_dangling = dangling
        self._last_corrected = corrected
        self._last_problematic = problematic
        return GxffixResult(
            doc=self.doc, dangling_records=dangling,
            corrected_transcripts=corrected,
            problematic_transcripts=problematic, stats=self.stats)

    def write_outputs(self, fixed_path: str, dangling_path: str,
                      corrected_path: str, problematic_path: str):
        """落盘:主文件+条件附加文件|Write main + conditional extras"""
        result_doc = self.doc
        write_gff3(result_doc, result_doc.records, fixed_path)
        written = {'fixed': True, 'dangling': False, 'corrected': False,
                   'problematic': False}
        if self._last_dangling:
            write_gff3(result_doc, self._last_dangling, dangling_path)
            written['dangling'] = True
        for ids, path, key in (
                (self._last_corrected, corrected_path, 'corrected'),
                (self._last_problematic, problematic_path, 'problematic')):
            if ids:
                recs = self._collect_with_genes(ids)
                write_gff3(result_doc, recs, path)
                written[key] = True
        return written

    def _collect_with_genes(self, transcript_ids) -> List[GffRecord]:
        """转录本全部记录+父gene|All records of transcripts + parent genes"""
        self.doc.rebuild_indexes()
        tx_ids = set(transcript_ids)
        wanted = set(tx_ids)
        for tid in transcript_ids:
            tx = self.doc.id_index.get(tid)
            if tx and tx.attributes.get('Parent'):
                for p in tx.attributes['Parent'].split(','):
                    g = self.doc.id_index.get(p)
                    if g is not None and g.type == 'gene':
                        wanted.add(p)
        collected = []
        for r in self.doc.records:
            if r.attributes.get('ID') in wanted:
                collected.append(r)
                continue
            parent = r.attributes.get('Parent', '')
            # 子行常无 ID(exon/CDS/UTR),按 Parent 归属命中转录本,否则
            # corrected/problematic 文件只剩 mRNA 行,无法人工核对相位
            # |Children are usually ID-less; claim them by Parent so the
            # extras carry the full transcript, not the mRNA row alone
            if parent and parent.split(',')[0] in tx_ids:
                collected.append(r)
        return collected
