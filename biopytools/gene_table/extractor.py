"""基因信息+序列合并表核心|gene_table core: GFF parse, gene-DNA slice, gffread, merge"""

import csv
import subprocess
from typing import Dict, List

from .utils import (_open_text, build_conda_command, extract_sequence_region,
                    read_fasta_to_dict, write_fasta, format_number)

# 输出表列顺序(固定)|Fixed output column order
COLUMNS = ['Sample', 'Gene_ID', 'Transcript_ID', 'Chromosome', 'Strand',
           'Gene_Start', 'Gene_End', 'Transcript_Start', 'Transcript_End',
           'Gene_DNA', 'CDS', 'Protein']


def _parse_attrs(attr_field: str) -> Dict[str, str]:
    """解析 GFF 第9列属性|Parse GFF column-9 attributes"""
    out = {}
    for piece in attr_field.split(';'):
        if '=' in piece:
            k, v = piece.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def _load_genome(genome_file: str) -> Dict[str, str]:
    """读基因组 FASTA→{seqid: seq}(大写)|Load genome FASTA into dict (uppercase)"""
    seqs: Dict[str, str] = {}
    cur_id, chunks = None, []
    with _open_text(genome_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(chunks)
                cur_id = line[1:].split()[0] if len(line) > 1 else ''
                chunks = []
            elif line:
                chunks.append(line.strip().upper())
    if cur_id is not None:
        seqs[cur_id] = ''.join(chunks)
    return seqs


class GeneTableExtractor:
    """基因信息+序列合并表提取器|Gene info + sequence merged table extractor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    # ---- GFF 解析|GFF parse ----
    def parse_gff(self):
        """一遍扫描 GFF,收集 gene / transcript / cds_len|Single-pass GFF scan

        Returns:
            (genes, transcripts, cds_len)
            genes: {gene_id: {chr, start, end, strand}}
            transcripts: [{gene_id, transcript_id, chr, strand, t_start, t_end}]
            cds_len: {transcript_id: 总CDS长度} (用于 --longest-only)
        """
        self.logger.info(f"解析 GFF|Parsing GFF: {self.config.gff_file}")
        genes: Dict[str, Dict] = {}
        transcripts: List[Dict] = []
        cds_len: Dict[str, int] = {}
        transcript_types = {t.lower() for t in self.config.transcript_types}
        gene_type = self.config.gene_type.lower()

        with _open_text(self.config.gff_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) != 9:
                    continue
                seqid, _src, ftype, start, end, _score, strand, _phase, attr_field = fields
                ftype_l = ftype.lower()

                if ftype_l == gene_type:
                    attrs = _parse_attrs(attr_field)
                    gid = attrs.get('ID')
                    if gid:
                        genes[gid] = {'chr': seqid, 'start': int(start),
                                      'end': int(end), 'strand': strand}
                elif ftype_l in transcript_types:
                    attrs = _parse_attrs(attr_field)
                    tid = attrs.get('ID')
                    parent = attrs.get('Parent')
                    if tid and parent:
                        # Parent 可能逗号分隔,取第一个|Parent may be comma-separated; take first
                        gene_id = parent.split(',')[0]
                        transcripts.append({
                            'gene_id': gene_id, 'transcript_id': tid,
                            'chr': seqid, 'strand': strand,
                            't_start': int(start), 't_end': int(end)})
                elif ftype_l == 'cds':
                    attrs = _parse_attrs(attr_field)
                    parent = attrs.get('Parent')
                    if parent:
                        tid = parent.split(',')[0]
                        cds_len[tid] = cds_len.get(tid, 0) + (int(end) - int(start) + 1)

        self.logger.info(
            f"GFF 统计|GFF stats: 基因|genes={format_number(len(genes))} "
            f"转录本|transcripts={format_number(len(transcripts))} "
            f"有CDS|with_CDS={format_number(len(cds_len))}")
        return genes, transcripts, cds_len

    # ---- 基因 DNA 切片|gene-DNA slice ----
    def slice_gene_dna(self, genes: Dict[str, Dict]) -> Dict[str, str]:
        """加载基因组,切每个基因全长 DNA(含内含子/UTR)|Load genome, slice full-length gene DNA"""
        self.logger.info(f"加载基因组|Loading genome: {self.config.genome_file}")
        genome = _load_genome(self.config.genome_file)
        self.logger.info(f"基因组载入|Genome loaded: {format_number(len(genome))} 条序列|sequences")

        dna: Dict[str, str] = {}
        skipped = 0
        for gid, g in genes.items():
            chrom = genome.get(g['chr'])
            if chrom is None:
                # 优雅降级:染色体缺失则跳过并告警|graceful: skip + warn on missing chrom
                skipped += 1
                self.logger.warning(
                    f"基因|gene {gid} 染色体|chrom {g['chr']} 不在基因组|not in genome")
                continue
            seq = extract_sequence_region(chrom, g['start'], g['end'], g['strand'])
            if len(seq) < self.config.min_length:
                skipped += 1
                continue
            dna[gid] = seq
        if skipped:
            self.logger.warning(
                f"跳过|Skipped {skipped} 个基因(染色体缺失/长度不足)|genes (missing chrom / too short)")
        self.logger.info(f"切出基因 DNA|Sliced gene DNA: {format_number(len(dna))}")
        return dna

    # ---- 最长转录本选择|longest-transcript selection ----
    def _select_longest(self, transcripts: List[Dict], cds_len: Dict[str, int]) -> List[Dict]:
        """每基因保留最长转录本(CDS 长度优先,无 CDS 按转录本长度)|One longest transcript per gene"""
        best: Dict[str, Dict] = {}
        for t in transcripts:
            gid = t['gene_id']
            score = cds_len.get(t['transcript_id'], 0)
            # 无 CDS 时按转录本长度兜底|fall back to transcript span when no CDS
            score = score if score > 0 else (t['t_end'] - t['t_start'])
            cur = best.get(gid)
            if cur is None or score > cur['_score']:
                t2 = dict(t)
                t2['_score'] = score
                best[gid] = t2
        result = [{k: v for k, v in t.items() if k != '_score'} for t in best.values()]
        self.logger.info(f"最长转录本模式|Longest-only: {format_number(len(result))} 个基因|genes")
        return result

    # ---- gffread 出 CDS/蛋白|gffread for CDS/protein ----
    def run_gffread(self, mode: str) -> Dict[str, str]:
        """运行 gffread 出 CDS(-x) 或蛋白(-y),返回 {transcript_id: seq}|Run gffread for CDS/pep"""
        assert mode in ('cds', 'pep')
        flag = '-x' if mode == 'cds' else '-y'
        out_fa = self.config.cds_fa if mode == 'cds' else self.config.pep_fa
        label = 'CDS' if mode == 'cds' else '蛋白|Protein'
        args = [flag, out_fa, '-g', self.config.genome_file, self.config.gff_file]
        cmd = build_conda_command(self.config.gffread_path, args)
        self.logger.info(f"执行|Executing: gffread {label}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.error(
                f"gffread {label} 失败|failed (rc={result.returncode}): {result.stderr}")
            raise RuntimeError(f"gffread {mode} 失败|failed")
        seqs = read_fasta_to_dict(out_fa)
        self.logger.info(f"{label} 产出|produced: {format_number(len(seqs))} 条|records")
        return seqs

    # ---- 合并写表|merge + write ----
    def merge_and_write(self, genes, transcripts, gene_dna, cds, pep):
        """按转录本合并并写 TSV + gene.fa|Merge by transcript, write TSV + gene.fa"""
        sample = self.config.prefix
        na_cds = na_pep = na_dna = 0
        rows = []
        for t in transcripts:
            gid, tid = t['gene_id'], t['transcript_id']
            g = genes.get(gid, {})
            dna_seq = gene_dna.get(gid)
            cds_seq = cds.get(tid)
            pep_seq = pep.get(tid)
            if dna_seq is None:
                na_dna += 1
                dna_seq = 'NA'
            if cds_seq is None:
                na_cds += 1
                cds_seq = 'NA'
            if pep_seq is None:
                na_pep += 1
                pep_seq = 'NA'
            rows.append({
                'Sample': sample, 'Gene_ID': gid, 'Transcript_ID': tid,
                'Chromosome': t['chr'], 'Strand': t['strand'],
                'Gene_Start': g.get('start', 'NA'), 'Gene_End': g.get('end', 'NA'),
                'Transcript_Start': t['t_start'], 'Transcript_End': t['t_end'],
                'Gene_DNA': dna_seq, 'CDS': cds_seq, 'Protein': pep_seq})

        # lineterminator='\n' 强制 LF(csv 默认 \r\n 会使末列带 \r)|force LF; csv default \r\n leaves CR on last column
        with open(self.config.tsv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t',
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

        # gene.fa 由我们写出(基因全长DNA);cds/pep 已由 gffread 写出
        # |gene.fa written by us; cds/pep already produced by gffread
        write_fasta(sorted(gene_dna.items()), self.config.gene_fa)

        self.logger.info(
            f"写出合并表|Wrote merged table: {format_number(len(rows))} 行|rows → {self.config.tsv_path}")
        self.logger.info(
            f"缺失统计|Missing: Gene_DNA={na_dna} CDS={na_cds} Protein={na_pep}")
        return self.config.tsv_path

    # ---- 主流程|main pipeline ----
    def run(self):
        """主流程:解析→(可选最长)→切基因DNA→gffread→合并写表|Main pipeline"""
        self.logger.info("=== gene-table 开始|start ===")
        genes, transcripts, cds_len = self.parse_gff()
        if self.config.longest_only:
            transcripts = self._select_longest(transcripts, cds_len)
        gene_dna = self.slice_gene_dna(genes)
        cds = self.run_gffread('cds')
        pep = self.run_gffread('pep')
        out = self.merge_and_write(genes, transcripts, gene_dna, cds, pep)
        self.logger.info("=== gene-table 完成|done ===")
        return out
