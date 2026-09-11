"""质量诊断|Quality diagnosis

结构诊断(相位/长度/UTR)+ CPC2 + RNA 零覆盖;输出 TSV + summary
|Structural (phase/length/UTR) + CPC2 + RNA zero-coverage; TSV + summary
"""

import json
import logging
import os
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .config import AnnoCurateConfig
from .gff_model import (TRANSCRIPT_RE, UTR3_TYPES, UTR5_TYPES, GxfDocument)
from .phase import validate_transcript_phases


def _is_transcript(rec):
    """转录本级判定|Transcript-level check"""
    return bool(TRANSCRIPT_RE.search(rec.type))


@dataclass
class TranscriptIssue:
    """单转录本问题|Per-transcript issues"""

    tid: str
    seqid: str
    strand: str
    cds_len: int
    tags: List[str] = field(default_factory=list)
    zero_cov_exons: Optional[str] = None    # 'n/m',无RNA校验时None


def effective_utr_thresholds(config: AnnoCurateConfig) -> Tuple[float, float]:
    """有效UTR阈值=基准×(1+relax)|Effective thresholds (spec §7)"""
    factor = 1 + config.relax
    return config.utr5_frac_max * factor, config.utr3_frac_max * factor


def structural_diagnosis(doc: GxfDocument, config: AnnoCurateConfig,
                         logger: logging.Logger) -> List[TranscriptIssue]:
    """结构诊断:相位/长度/UTR|Structural diagnosis"""
    doc.rebuild_indexes()
    utr5_lim, utr3_lim = effective_utr_thresholds(config)
    issues: Dict[str, TranscriptIssue] = {}
    for rec in doc.records:
        if not _is_transcript(rec):
            continue
        tid = rec.attributes.get('ID', '')
        children = doc.children_by_parent.get(tid, [])
        cds = [c for c in children if c.type == 'CDS']
        cds_len = sum(c.end - c.start + 1 for c in cds)
        tags: List[str] = []
        if cds:
            phase_res = validate_transcript_phases(cds)
            # 修复后 CORRECTED 已就地覆盖;仍见 NO/INVALID 说明修复未跑或
            # 诊断对象是修复前文档|post-fix leftovers signal real issues
            if phase_res.status == 'INVALID_LENGTH':
                tags.append(f'CDS 总长非3倍数|total CDS length not '
                            f'multiple of 3: {cds_len}')
            elif phase_res.status == 'CORRECTED':
                tags.extend('NO PHASE' if 'NO PHASE' in i else 'INVALID PHASE'
                            for i in phase_res.issues)
            if cds_len < config.min_cds_nt:
                tags.append('CDS too short')
            elif cds_len > config.max_cds_nt:
                tags.append('CDS too long')
            if config.check_utr:
                span = rec.end - rec.start + 1
                utr5 = sum(c.end - c.start + 1 for c in children
                           if c.type.lower() in UTR5_TYPES)
                utr3 = sum(c.end - c.start + 1 for c in children
                           if c.type.lower() in UTR3_TYPES)
                if span > 0:
                    if utr5 / span > utr5_lim:
                        tags.append(f"5'UTR fraction {utr5/span:.2f} "
                                    f"exceeds {utr5_lim:.2f}")
                    if utr3 / span > utr3_lim:
                        tags.append(f"3'UTR fraction {utr3/span:.2f} "
                                    f"exceeds {utr3_lim:.2f}")
        if tags:
            issues[tid] = TranscriptIssue(tid, rec.seqid, rec.strand,
                                          cds_len, tags)
    return list(issues.values())


def write_diagnosis(issues: List[TranscriptIssue], all_tids: List[str],
                    out_tsv: str, out_summary: str, extras: dict):
    """写诊断表与摘要|Write diagnosis TSV and summary"""
    with open(out_tsv, 'w', encoding='utf-8') as fh:
        fh.write('#Found totally {} transcripts with potential issues.\n'
                 .format(len(issues)))
        fh.write('\t'.join(('transcript_id', 'chr', 'strand', 'cds_len',
                            'zero_cov_exons', 'tags')) + '\n')
        for i in issues:
            fh.write('\t'.join((i.tid, i.seqid, i.strand, str(i.cds_len),
                                i.zero_cov_exons or '',
                                ';'.join(i.tags))) + '\n')
    counts: Dict[str, int] = {}
    for i in issues:
        for t in i.tags:
            key = ('NO PHASE' if t == 'NO PHASE' else
                   'INVALID PHASE' if t == 'INVALID PHASE' else
                   'CDS too short' if t == 'CDS too short' else
                   'CDS too long' if t == 'CDS too long' else
                   'NO RNA COVERAGE' if t == 'NO RNA COVERAGE' else
                   t.split(' ')[0] if t.startswith('ZERO-COV') else t)
            counts[key] = counts.get(key, 0) + 1
    with open(out_summary, 'w', encoding='utf-8') as fh:
        fh.write(f'total transcripts: {len(all_tids)}\n')
        fh.write(f'transcripts with issues: {len(issues)} '
                 f'({len(issues)/len(all_tids)*100:.1f}%)\n'
                 if all_tids else 'transcripts with issues: 0\n')
        for k in sorted(counts):
            fh.write(f'{k}: {counts[k]}\n')
        fh.write('CPC2: ' + str(extras.get('cpc2_status',
                                           'skipped')) + '\n')
        fh.write('RNA coverage: ' + str(extras.get('rna_status',
                                                   'skipped')) + '\n')


def diagnosis_fingerprint(config: AnnoCurateConfig) -> dict:
    """诊断参数指纹|Diagnosis parameter fingerprint"""
    return {
        'check_utr': config.check_utr,
        'min_cds_nt': config.min_cds_nt,
        'max_cds_nt': config.max_cds_nt,
        'utr5_frac_max': config.utr5_frac_max,
        'utr3_frac_max': config.utr3_frac_max,
        'relax': config.relax,
        'genome': os.path.abspath(config.genome) if config.genome else None,
        'rna_bams': [os.path.abspath(b) for b in (config.rna_bams or ())],
    }


def save_fingerprint(path: str, fp: dict):
    """存指纹|Save fingerprint"""
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(fp, fh, indent=2, sort_keys=True)


def load_fingerprint(path: str) -> Optional[dict]:
    """读指纹|Load fingerprint"""
    try:
        with open(path, encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def should_skip_diagnosis(out_tsv: str, meta_path: str,
                          config: AnnoCurateConfig) -> bool:
    """两级断点之诊断级|Diagnosis-level resume check"""
    if config.force:
        return False
    if not os.path.exists(out_tsv):
        return False
    return load_fingerprint(meta_path) == diagnosis_fingerprint(config)


# ---------- CPC2 编码潜能|CPC2 coding potential ----------

_COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


class FastaIndex:
    """纯Python FASTA偏移索引|Pure-Python FASTA offset index"""

    def __init__(self, fasta_path: str):
        self.path = fasta_path
        self.entries = {}       # seqid -> (offset, line_bpp, seq_len)
        fai = fasta_path + '.fai'
        if os.path.exists(fai):
            self._load_fai(fai)
        else:
            self._scan()

    def _load_fai(self, fai: str):
        """解析现成.fai|Parse existing .fai

        fetch 的换行补偿需要 linewidth(col5,含换行整行字节数),
        不是 linebases(col4,每行碱基数)——用错会让多行 FASTA 行尾串位
        |fetch's newline compensation needs linewidth (col5, bytes incl.
        newline), not linebases (col4); mixing them corrupts multiline fetch
        """
        with open(fai, encoding='utf-8') as fh:
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 5:
                    self.entries[cols[0]] = (int(cols[2]), int(cols[4]),
                                             int(cols[1]))

    def _scan(self):
        """扫描建索引|Scan and index"""
        with open(self.path, 'rb') as fh:
            offset = 0
            current = None
            first_len = None
            line_bpp = None
            seq_len = 0
            for raw in fh:
                n = len(raw)
                if raw.startswith(b'>'):
                    if current:
                        self.entries[current] = (offset, line_bpp, seq_len)
                    name = raw[1:].split()[0].decode() \
                        if raw[1:].split() else ''
                    current = name
                    offset = fh.tell()
                    first_len = None
                    seq_len = 0
                elif current is not None:
                    stripped = raw.rstrip(b'\r\n')
                    seq_len += len(stripped)
                    if first_len is None:
                        first_len = len(stripped)
                    # 仅真实带换行的整行可刷新 line_bpp:末行无 \n 时
                    # n=纯碱基长(如 60 而非 61),采纳会把每行步长少算
                    # 1 字节,此后该序列 fetch 全错|Refresh line_bpp only
                    # for newline-terminated full lines; an unterminated
                    # last line would undercount the per-line stride and
                    # corrupt every fetch of this entry
                    if len(stripped) == first_len and n > len(stripped):
                        line_bpp = n
                    elif line_bpp is None:
                        # 单行且无换行的兜底:len+1 使换行补偿恒为 0
                        # |Single unterminated line: len+1 zeroes the
                        # newline compensation
                        line_bpp = len(stripped) + 1
            if current:
                self.entries[current] = (offset, line_bpp, seq_len)

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """取子序列(1-based含端)|Fetch subsequence (1-based inclusive)"""
        if seqid not in self.entries:
            return ''
        offset, line_bpp, _ = self.entries[seqid]
        span = end - start + 1
        if span <= 0:
            return ''
        # 按行边界整段读取:seek 到起点字节后按"碱基数+换行上界"一次
        # read 覆盖整个区间,剥掉换行再切片,免逐碱基 seek/read
        # (千万碱基级=千万次系统调用)|Bulk read by line boundaries: seek
        # once to the start byte, read span bytes plus a newline upper
        # bound in one go, strip newlines, slice; avoids the per-base
        # seek/read syscall storm
        stride = max(1, line_bpp - 1)
        zero = start - 1
        parts = []
        got = 0
        with open(self.path, 'rb') as fh:
            fh.seek(offset + zero + zero // stride)
            while got < span:
                chunk = fh.read(span - got + (span - got) // stride + 8)
                if not chunk:
                    break
                text = chunk.decode('ascii', 'replace')
                for nl in ('\n', '\r'):
                    text = text.replace(nl, '')
                parts.append(text)
                got += len(text)
        return ''.join(parts)[:span]


def _revcomp(seq: str) -> str:
    """反向互补|Reverse complement"""
    return seq.translate(_COMPLEMENT)[::-1]


def extract_cds_fasta(doc: GxfDocument, index: FastaIndex,
                      tid: str) -> Tuple[str, str]:
    """提取单转录本CDS序列|Extract one transcript's CDS sequence

    调用前须已 rebuild_indexes|caller must rebuild indexes first
    """
    tx = doc.id_index.get(tid)
    strand = tx.strand if tx else '+'
    cds = [c for c in doc.children_by_parent.get(tid, []) if c.type == 'CDS']
    cds = sorted(cds, key=lambda c: c.start, reverse=(strand == '-'))
    parts = []
    for c in cds:
        s = index.fetch(c.seqid, c.start, c.end).upper()
        if not s:
            raise KeyError(f"{tid}: {c.seqid} 不在FASTA|not in FASTA")
        parts.append(_revcomp(s) if strand == '-' else s)
    return tid, f'>{tid}\n' + ''.join(parts) + '\n'


def run_cpc2(cds_fasta_path: str, out_dir: str, config: AnnoCurateConfig,
             logger: logging.Logger) -> Optional[str]:
    """运行CPC2,返回输出txt路径|Run CPC2, return output txt path"""
    from biopytools.common.conda_runner import build_conda_command
    prefix = os.path.join(out_dir, 'cpc2')
    cmd = build_conda_command(config.cpc2_path, ['-i', cds_fasta_path,
                                                 '-o', prefix])
    logger.info("执行|Executing: CPC2 编码潜能|coding potential")
    logger.info(f"命令|Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    out_txt = prefix + '.txt'
    if result.returncode != 0 or not os.path.exists(out_txt):
        logger.warning(f"CPC2失败,降级跳过|CPC2 failed, degraded: "
                       f"{(result.stderr or '')[:200]}")
        return None
    return out_txt


def parse_cpc2_output(txt_path: str) -> Dict[str, bool]:
    """CPC2输出→tid是否noncoding|CPC2 output → noncoding flags"""
    out: Dict[str, bool] = {}
    with open(txt_path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 8:
                out[cols[0]] = cols[-1].strip().lower() == 'noncoding'
    return out


def cpc2_diagnosis(doc: GxfDocument, config: AnnoCurateConfig,
                   logger: logging.Logger,
                   tmp_dir: str) -> Tuple[List[str], str]:
    """CPC2 诊断总入口|CPC2 diagnosis entry"""
    if not config.genome:
        return [], '跳过|skipped (no --genome)'
    if not os.path.exists(config.cpc2_path):
        logger.warning(f"CPC2未找到,降级|CPC2 not found, degraded: "
                       f"{config.cpc2_path}")
        return [], '跳过|skipped (cpc2 not found)'
    os.makedirs(tmp_dir, exist_ok=True)
    try:
        index = FastaIndex(config.genome)
        records = []
        extract_failed = 0
        doc.rebuild_indexes()
        for rec in doc.records:
            if not _is_transcript(rec):
                continue
            tid = rec.attributes.get('ID', '')
            # 只收有CDS的转录本,免得空记录混进CPC2输入
            # |CDS-bearing transcripts only, no empty records into CPC2
            cds = [c for c in doc.children_by_parent.get(tid, [])
                   if c.type == 'CDS']
            if not cds:
                continue
            try:
                records.append(extract_cds_fasta(doc, index, tid)[1])
            except KeyError as e:
                # seqid 不在FASTA:逐条 WARNING 跳过,其余转录本照常进
                # CPC2(spec §10),不整体降级|Missing seqid: warn per
                # transcript and skip; the rest still go through CPC2
                extract_failed += 1
                logger.warning(
                    f"CPC2提取失败,跳过该转录本|CPC2 extract failed, "
                    f"transcript skipped: {e}")
        if not records:
            if extract_failed:
                return [], (f'失败降级|degraded (all {extract_failed} CDS '
                            f'extracts failed)')
            return [], '跳过|skipped (no CDS extracted)'
        fa = os.path.join(tmp_dir, 'cds.fa')
        with open(fa, 'w', encoding='utf-8') as fh:
            fh.writelines(records)
        out = run_cpc2(fa, tmp_dir, config, logger)
        if out is None:
            return [], '失败降级|failed, degraded'
        flags = parse_cpc2_output(out)
        nc = [t for t, is_nc in flags.items() if is_nc]
        status = f'完成|done (noncoding {len(nc)}/{len(flags)}'
        if extract_failed:
            status += f', {extract_failed} extract skipped'
        return nc, status + ')'
    except Exception as e:                       # noqa: BLE001 降级不崩主流程
        logger.warning(f"CPC2异常,降级|CPC2 error, degraded: {e}")
        return [], f'失败降级|degraded ({e})'


# ---------- RNA 零覆盖证据|RNA zero-coverage evidence ----------


def export_exon_bed(doc: GxfDocument, bed_path: str) -> int:
    """外显子导出BED(0-based半开)|Export exons as BED (0-based half-open)"""
    rows = 0
    with open(bed_path, 'w', encoding='utf-8') as fh:
        for rec in doc.records:
            if not _is_transcript(rec):
                continue
            tid = rec.attributes.get('ID', '')
            exons = [c for c in doc.children_by_parent.get(tid, [])
                     if c.type == 'exon']
            for k, ex in enumerate(exons, 1):
                fh.write(f'{ex.seqid}\t{ex.start - 1}\t{ex.end}\t'
                         f'{tid}:{k}\n')
                rows += 1
    return rows


def run_bedcov(bed_path: str, bam_path: str, config: AnnoCurateConfig,
               logger: logging.Logger) -> Optional[Dict[str, int]]:
    """单BAM bedcov→exon reads数|bedcov for one BAM → per-exon reads"""
    import subprocess
    from biopytools.common.conda_runner import build_conda_command
    if not (os.path.exists(bam_path + '.bai') or
            os.path.exists(bam_path + '.csi')):
        cmd_idx = build_conda_command(config.samtools_path,
                                      ['index', '-b', bam_path])
        logger.info(f"命令|Command: {' '.join(cmd_idx)}")
        r = subprocess.run(cmd_idx, shell=False, capture_output=True,
                           text=True)
        if r.returncode != 0:
            logger.warning(f"索引失败,跳过该BAM|index failed, skip: "
                           f"{bam_path}")
            return None
    cmd = build_conda_command(config.samtools_path,
                              ['bedcov', bed_path, bam_path])
    logger.info(f"命令|Command: {' '.join(cmd)}")
    r = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    if r.returncode != 0:
        logger.warning(f"bedcov失败,跳过该BAM|bedcov failed, skip: "
                       f"{bam_path}")
        return None
    out: Dict[str, int] = {}
    for line in (r.stdout or '').splitlines():
        cols = line.split('\t')
        if len(cols) >= 5:
            out[cols[3]] = int(cols[4])
    return out


def rna_coverage_diagnosis(doc: GxfDocument, config: AnnoCurateConfig,
                           logger: logging.Logger, tmp_dir: str
                           ) -> Tuple[Dict[str, str], set, str]:
    """RNA零覆盖两级判定|Two-level zero-coverage verdict

    返回 (zero_map, no_cov_tids, status):zero_map = tid -> 'n/m';
    no_cov_tids = 全部外显子零覆盖的转录本集合
    |Returns (zero_map, no_cov_tids, status)
    """
    if not config.rna_bams:
        return {}, set(), '跳过|skipped (no --rna-bam)'
    os.makedirs(tmp_dir, exist_ok=True)
    doc.rebuild_indexes()
    bed = os.path.join(tmp_dir, 'exons.bed')
    export_exon_bed(doc, bed)
    per_exon_any: Dict[str, int] = {}
    ok_bams = 0
    for bam in config.rna_bams:
        res = run_bedcov(bed, bam, config, logger)
        if res is None:
            continue
        ok_bams += 1
        for key, reads in res.items():
            # 任一 BAM 有 reads 即算覆盖|covered in any BAM
            per_exon_any[key] = max(per_exon_any.get(key, 0), reads)
    if ok_bams == 0:
        return {}, set(), '失败降级|degraded (all BAMs failed)'
    zero_map: Dict[str, str] = {}
    for rec in doc.records:
        if not _is_transcript(rec):
            continue
        tid = rec.attributes.get('ID', '')
        exons = [c for c in doc.children_by_parent.get(tid, [])
                 if c.type == 'exon']
        if not exons:
            continue
        zero = sum(1 for k in range(1, len(exons) + 1)
                   if per_exon_any.get(f'{tid}:{k}', 0) == 0)
        if zero:
            zero_map[tid] = f'{zero}/{len(exons)}'
    no_cov_tids = {tid for tid, ratio in zero_map.items()
                   if ratio.split('/')[0] == ratio.split('/')[1]}
    return zero_map, no_cov_tids, (
        f'完成|done ({ok_bams} BAMs, {len(zero_map)} transcripts with '
        f'zero-cov exons)')


def merge_diagnosis(issues: List[TranscriptIssue],
                    zero_map: Dict[str, str],
                    no_cov_tids,
                    tx_meta: Dict[str, tuple]) -> List[TranscriptIssue]:
    """RNA结果并入结构诊断|Merge RNA results into structural issues

    tx_meta: tid -> (seqid, strand, cds_len),新入表的转录本填元数据
    |tx_meta: tid -> (seqid, strand, cds_len) for newly-added rows
    """
    no_cov = set(no_cov_tids)
    info: Dict[str, TranscriptIssue] = {i.tid: i for i in issues}
    for tid, ratio in zero_map.items():
        issue = info.get(tid)
        if issue is None:
            seqid, strand, cds_len = tx_meta.get(tid, ('.', '.', 0))
            issue = TranscriptIssue(tid, seqid, strand, cds_len, [])
            info[tid] = issue
        issue.zero_cov_exons = ratio
        if tid in no_cov:
            issue.tags.append('NO RNA COVERAGE')
        else:
            issue.tags.append(f'ZERO-COV EXON {ratio}')
    return list(info.values())
