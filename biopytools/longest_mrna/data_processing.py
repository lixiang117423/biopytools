"""
最长转录本提取数据处理模块 | Longest mRNA Extraction Data Processing Module
"""

import os
from collections import defaultdict
from typing import Dict, List, Tuple, Any

class GFF3Parser:
    """GFF3文件解析器 | GFF3 File Parser"""
    
    def __init__(self, gff3_path: str, logger):
        self.gff3_path = gff3_path
        self.logger = logger
        self.transcripts = defaultdict(list)
    
    def parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """解析GFF3属性字符串 | Parse GFF3 attributes string"""
        attributes = {}
        for item in attr_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attributes[key] = value
        return attributes
    
    def parse(self) -> Dict[str, List[Dict[str, Any]]]:
        """解析GFF3文件 | Parse GFF3 file"""
        self.logger.info(f"解析GFF3文件 | Parsing GFF3 file: {self.gff3_path}")
        
        with open(self.gff3_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 9:
                    self.logger.warning(f"行 {line_num} 格式不正确，跳过 | Line {line_num} has incorrect format, skipping")
                    continue
                
                feature_type = fields[2]
                attributes = self.parse_attributes(fields[8])
                
                if feature_type == 'mRNA':
                    parent_gene = attributes.get('Parent', '').split(':')[-1]
                    transcript_info = {
                        'transcript_id': attributes.get('ID', ''),
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
                    
                    # 找到对应的转录本并添加外显子信息 | Find corresponding transcript and add exon info
                    for gene_transcripts in self.transcripts.values():
                        for transcript in gene_transcripts:
                            if transcript['transcript_id'] == parent_transcript:
                                transcript['exons'].append(exon_info)
                                break
        
        self.logger.info(f"解析完成，找到 {len(self.transcripts)} 个基因 | Parsing completed, found {len(self.transcripts)} genes")
        return dict(self.transcripts)

class CDSCalculator:
    """CDS长度计算器 | CDS Length Calculator"""
    
    def __init__(self, logger):
        self.logger = logger
        self.transcript_lengths = defaultdict(int)
        self.gene_metadata = defaultdict(dict)
    
    def calculate_from_gff(self, gff_path: str) -> Dict[str, Dict[str, Any]]:
        """从GFF文件计算CDS长度 | Calculate CDS length from GFF file"""
        self.logger.info(f"计算CDS长度 | Calculating CDS lengths from: {gff_path}")
        
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
                    gene_id = attributes.get('ID', '')
                    self.gene_metadata[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }
                
                elif feature_type == 'mRNA':
                    transcript_id = attributes.get('ID', '')
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
        
        # 找到每个基因的最长转录本 | Find longest transcript for each gene
        longest_transcripts = {}
        for gene_id, transcripts in gene_transcripts.items():
            if transcripts:
                longest_transcript = max(
                    transcripts, 
                    key=lambda x: self.transcript_lengths.get(x['id'], 0)
                )
                longest_transcripts[gene_id] = longest_transcript
        
        self.logger.info(f"找到 {len(longest_transcripts)} 个基因的最长转录本 | Found longest transcripts for {len(longest_transcripts)} genes")
        return longest_transcripts

class TranscriptProcessor:
    """转录本处理器 | Transcript Processor"""
    
    def __init__(self, logger):
        self.logger = logger
    
    @staticmethod
    def calculate_transcript_length(transcript: Dict[str, Any]) -> int:
        """计算转录本长度 | Calculate transcript length"""
        if 'exons' in transcript and transcript['exons']:
            return sum(end - start + 1 for start, end in transcript['exons'])
        else:
            return transcript.get('end', 0) - transcript.get('start', 0) + 1
    
    def get_longest_transcript(self, gene_transcripts: List[Dict[str, Any]]) -> Dict[str, Any]:
        """获取最长转录本 | Get longest transcript"""
        if not gene_transcripts:
            return {}
        
        return max(gene_transcripts, key=self.calculate_transcript_length)
    
    def process_all_genes(self, gene_transcripts: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Dict[str, Any]]:
        """处理所有基因，获取最长转录本 | Process all genes to get longest transcripts"""
        self.logger.info("处理基因转录本，选择最长的 | Processing gene transcripts, selecting longest")
        
        longest_transcripts = {}
        for gene_id, transcripts in gene_transcripts.items():
            if transcripts:
                longest = self.get_longest_transcript(transcripts)
                longest_transcripts[gene_id] = longest
        
        return longest_transcripts
