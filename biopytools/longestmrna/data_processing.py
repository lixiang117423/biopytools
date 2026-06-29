"""
最长转录本提取数据处理模块|Longest mRNA Extraction Data Processing Module
"""

import os
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Set


class GFFGenomeAligner:
    """GFF与基因组序列名对齐器|Align GFF seqids with genome FASTA seqids

    解决场景|Solves: GFF的序列名与基因组FASTA不一致时，gffread会因为找不到序列而整体报错退出。
    本类在流程最前端检测不匹配，只保留基因组中存在的序列对应的注释行，跳过其余部分，
    从而"能提取哪些就提取哪些"，而非整体失败|
    When GFF seqids mismatch the genome FASTA, gffread aborts entirely because it cannot
    find the referenced sequences. This class detects the mismatch up front and keeps only
    the annotation lines whose seqid exists in the genome, skipping the rest, so we extract
    whatever is resolvable instead of failing the whole run.
    """

    def __init__(self, genome_file: str, gff3_file: str, logger):
        self.genome_file = genome_file
        self.gff3_file = gff3_file
        self.logger = logger

    def load_genome_seqids(self) -> Set[str]:
        """读取基因组FASTA的序列名集合|Load seqid set from genome FASTA

        取每个 > 头行第一个空白分隔的token作为序列名|
        Takes the first whitespace-delimited token of each > header as the seqid
        """
        seqids: Set[str] = set()
        with open(self.genome_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    if header:
                        seqids.add(header.split()[0])
        return seqids

    def collect_gff_seqids(self) -> Set[str]:
        """收集GFF特征行中出现的序列名|Collect seqids appearing in GFF feature lines"""
        seqids: Set[str] = set()
        with open(self.gff3_file, 'r', encoding='utf-8') as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    continue
                fields = line.split('\t')
                if fields and fields[0].strip():
                    seqids.add(fields[0].strip())
        return seqids

    def _write_filtered_gff(self, keep_seqids: Set[str], missing_seqids: Set[str],
                            output_path: str) -> Tuple[int, int]:
        """写出过滤后的GFF|Write filtered GFF

        保留：空行；普通注释行(单个#)；不引用被丢弃序列的 ## 指令(如 ##gff-version)；
        以及第一列序列名在keep_seqids中的特征行。丢弃：引用被丢弃序列的 ## 指令
        (如 ##sequence-region <missing_seqid>)，以及序列名不在keep_seqids中的特征行|
        Keeps blank lines, plain comments, ## directives that don't reference a dropped seqid,
        and feature lines whose column-1 seqid is in keep_seqids. Drops ## directives that
        reference a dropped seqid (e.g. ##sequence-region <missing_seqid>) and feature lines
        whose seqid is not in keep_seqids.
        因GFF3中同一基因的所有特征共享同一seqid，按seqid过滤可完整保留/丢弃整条基因|
        Since all features of a GFF3 gene share one seqid, seqid-based filtering keeps/drops
        whole genes intact.
        """
        kept = 0
        skipped = 0
        with open(self.gff3_file, 'r', encoding='utf-8') as fin, \
                open(output_path, 'w', encoding='utf-8') as fout:
            for line in fin:
                stripped = line.strip()
                if not stripped:
                    fout.write(line)
                    continue
                if stripped.startswith('##'):
                    # ##指令：第二token若是被丢弃的序列名则丢弃该指令，否则保留|
                    # ##directive: drop if its 2nd token is a dropped seqid, else keep
                    tokens = stripped.split()
                    if len(tokens) >= 2 and tokens[1] in missing_seqids:
                        skipped += 1
                        continue
                    fout.write(line)
                    continue
                if stripped.startswith('#'):
                    fout.write(line)  # 普通注释行原样保留|plain comment kept verbatim
                    continue
                fields = line.split('\t')
                seqid = fields[0].strip() if fields else ''
                if seqid in keep_seqids:
                    fout.write(line)
                    kept += 1
                else:
                    skipped += 1
        return kept, skipped

    def align(self, temp_manager) -> str:
        """对齐GFF与基因组序列名，返回应使用的GFF路径|Align and return the effective GFF path

        返回值|Returns:
            原 GFF 路径（完全匹配，无临时文件）或过滤后的临时 GFF 路径|
            The original GFF path (full match, no temp file) or a filtered temp GFF path

        异常|Raises:
            ValueError: GFF与基因组没有任何共有序列名，无法提取任何序列|
            ValueError: no common seqid between GFF and genome, nothing extractable
        """
        genome_seqids = self.load_genome_seqids()
        gff_seqids = self.collect_gff_seqids()

        self.logger.info(f"基因组序列数|Genome sequences: {len(genome_seqids)}")
        self.logger.info(f"GFF序列数|GFF sequences: {len(gff_seqids)}")

        missing = gff_seqids - genome_seqids
        if not missing:
            self.logger.info("GFF序列名与基因组完全匹配|All GFF seqids match the genome FASTA")
            return self.gff3_file

        matched = gff_seqids & genome_seqids
        self.logger.warning(
            f"检测到 {len(missing)} 个GFF序列在基因组FASTA中找不到，将跳过这些序列的注释|"
            f"{len(missing)} GFF seqid(s) not found in genome FASTA, their annotations will be skipped"
        )
        self.logger.warning(f"跳过的序列|Skipped seqids: {sorted(missing)}")
        self.logger.warning(f"保留的序列数|Kept seqids: {len(matched)}")

        if not matched:
            raise ValueError(
                f"GFF与基因组FASTA没有任何共有的序列名，无法提取序列|"
                f"No common seqid between GFF and genome FASTA, cannot extract sequences. "
                f"请检查 -g 与 -f 是否为配套文件|Please verify -g and -f are a matching pair"
            )

        # 写过滤后的临时GFF（由temp_manager统一清理）|Filtered temp GFF (cleaned up by temp_manager)
        temp_handle = temp_manager.create_temp_file(mode='w', suffix='.gff')
        temp_path = temp_handle.name
        temp_handle.close()
        kept, skipped = self._write_filtered_gff(matched, missing, temp_path)
        self.logger.warning(
            f"已写出过滤后的GFF：保留 {kept} 行，跳过 {skipped} 行|"
            f"Wrote filtered GFF: kept {kept} lines, skipped {skipped} lines"
        )
        return temp_path


class GFF3Parser:
    """GFF3文件解析器|GFF3 File Parser"""

    def __init__(self, gff3_path: str, logger):
        self.gff3_path = gff3_path
        self.logger = logger
        self.transcripts = defaultdict(list)

    def parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """解析GFF3属性字符串|Parse GFF3 attributes string"""
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes

    def parse(self) -> Dict[str, List[Dict[str, Any]]]:
        """解析GFF3文件|Parse GFF3 file"""
        self.logger.info(f"解析GFF3文件|Parsing GFF3 file: {self.gff3_path}")

        with open(self.gff3_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    self.logger.warning(f"行 {line_num} 格式不正确，跳过|Line {line_num} has incorrect format, skipping")
                    continue

                feature_type = fields[2]
                attributes = self.parse_attributes(fields[8])

                if feature_type == 'mRNA':
                    parent_gene = attributes.get('Parent', '').split(':')[-1]
                    transcript_info = {
                        'transcript_id': attributes.get('ID', '').split(':')[-1],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'chrom': fields[0],
                        'exons': []
                    }
                    self.transcripts[parent_gene].append(transcript_info)

                elif feature_type == 'exon':
                    parent_transcript = attributes.get('Parent', '').split(':')[-1]
                    exon_info = (int(fields[3]), int(fields[4]))

                    # 找到对应的转录本并添加外显子信息|Find corresponding transcript and add exon info
                    for gene_transcripts in self.transcripts.values():
                        for transcript in gene_transcripts:
                            if transcript['transcript_id'] == parent_transcript:
                                transcript['exons'].append(exon_info)
                                break

        self.logger.info(f"解析完成，找到 {len(self.transcripts)} 个基因|Parsing completed, found {len(self.transcripts)} genes")
        return dict(self.transcripts)

class CDSCalculator:
    """CDS长度计算器|CDS Length Calculator"""

    def __init__(self, logger):
        self.logger = logger
        self.transcript_lengths = defaultdict(int)
        self.gene_metadata = defaultdict(dict)

    def calculate_from_gff(self, gff_path: str) -> Dict[str, Dict[str, Any]]:
        """从GFF文件计算CDS长度|Calculate CDS length from GFF file"""
        self.logger.info(f"计算CDS长度|Calculating CDS lengths from: {gff_path}")

        gene_transcripts = defaultdict(list)

        with open(gff_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    continue

                chrom = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                attributes = {}
                for item in fields[8].split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        attributes[key] = value

                if feature_type == 'gene':
                    gene_id = attributes.get('ID', '').split(':')[-1]
                    self.gene_metadata[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }

                elif feature_type == 'mRNA':
                    transcript_id = attributes.get('ID', '').split(':')[-1]
                    parent_gene = attributes.get('Parent', '').split(':')[-1]
                    gene_transcripts[parent_gene].append({
                        'id': transcript_id,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'chrom': chrom
                    })

                elif feature_type == 'CDS':
                    transcript_id = attributes.get('Parent', '').split(':')[-1]
                    length = end - start + 1
                    self.transcript_lengths[transcript_id] += length

        # 找到每个基因的最长转录本|Find longest transcript for each gene
        longest_transcripts = {}
        for gene_id, transcripts in gene_transcripts.items():
            if transcripts:
                longest_transcript = max(
                    transcripts,
                    key=lambda x: self.transcript_lengths.get(x['id'], 0)
                )
                longest_transcripts[gene_id] = longest_transcript

        self.logger.info(f"找到 {len(longest_transcripts)} 个基因的最长转录本|Found longest transcripts for {len(longest_transcripts)} genes")
        return longest_transcripts

class TranscriptProcessor:
    """转录本处理器|Transcript Processor"""

    def __init__(self, logger):
        self.logger = logger

    @staticmethod
    def calculate_transcript_length(transcript: Dict[str, Any]) -> int:
        """计算转录本长度|Calculate transcript length"""
        if 'exons' in transcript and transcript['exons']:
            return sum(end - start + 1 for start, end in transcript['exons'])
        else:
            return transcript.get('end', 0) - transcript.get('start', 0) + 1

    def get_longest_transcript(self, gene_transcripts: List[Dict[str, Any]]) -> Dict[str, Any]:
        """获取最长转录本|Get longest transcript"""
        if not gene_transcripts:
            return {}

        return max(gene_transcripts, key=self.calculate_transcript_length)

    def process_all_genes(self, gene_transcripts: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Dict[str, Any]]:
        """处理所有基因，获取最长转录本|Process all genes to get longest transcripts"""
        self.logger.info("处理基因转录本，选择最长的|Processing gene transcripts, selecting longest")

        longest_transcripts = {}
        for gene_id, transcripts in gene_transcripts.items():
            if transcripts:
                longest = self.get_longest_transcript(transcripts)
                longest_transcripts[gene_id] = longest

        return longest_transcripts
