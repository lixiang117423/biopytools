"""
位置信息追踪模块 | Position Information Tracking Module
"""

import os
import gzip
from typing import Dict, List, Tuple, Optional, NamedTuple
from dataclasses import dataclass
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class KmerPosition(NamedTuple):
    """K-mer位置信息"""
    seq_name: str
    start_pos: int
    end_pos: int
    strand: str  # '+' or '-'

@dataclass
class KmerInfo:
    """完整的K-mer信息"""
    sequence: str
    positions: List[KmerPosition]
    total_count: int
    source_type: str  # 'fasta' or 'fastq'
    sample_name: str

class PositionTracker:
    """位置追踪器"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.kmer_positions: Dict[str, List[KmerPosition]] = {}
    
    def extract_fasta_kmers(self, fasta_file: str, sample_name: str) -> Dict[str, KmerInfo]:
        """从FASTA文件提取k-mer和位置信息"""
        self.logger.info(f"Extracting k-mers from FASTA: {fasta_file}")
        
        kmers = {}
        opener = gzip.open if fasta_file.endswith('.gz') else open
        
        try:
            with opener(fasta_file, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    sequence = str(record.seq).upper()
                    seq_name = record.id
                    
                    # 提取正向k-mer
                    for i in range(len(sequence) - self.config.kmer_size + 1):
                        kmer = sequence[i:i + self.config.kmer_size]
                        
                        # 跳过包含N的k-mer
                        if 'N' in kmer:
                            continue
                        
                        # 计算canonical k-mer
                        canonical_kmer = self._get_canonical_kmer(kmer)
                        strand = '+' if canonical_kmer == kmer else '-'
                        
                        position = KmerPosition(
                            seq_name=seq_name,
                            start_pos=i + 1,  # 1-based position
                            end_pos=i + self.config.kmer_size,
                            strand=strand
                        )
                        
                        if canonical_kmer not in kmers:
                            kmers[canonical_kmer] = KmerInfo(
                                sequence=canonical_kmer,
                                positions=[],
                                total_count=0,
                                source_type='fasta',
                                sample_name=sample_name
                            )
                        
                        kmers[canonical_kmer].positions.append(position)
                        kmers[canonical_kmer].total_count += 1
        
        except Exception as e:
            self.logger.error(f"Error processing FASTA file {fasta_file}: {e}")
            raise
        
        self.logger.info(f"Extracted {len(kmers)} unique k-mers from {fasta_file}")
        return kmers
    
    def extract_fastq_kmers(self, fastq_file: str, sample_name: str) -> Dict[str, KmerInfo]:
        """从FASTQ文件提取k-mer信息（不保存位置）"""
        self.logger.info(f"Extracting k-mers from FASTQ: {fastq_file}")
        
        kmers = {}
        opener = gzip.open if fastq_file.endswith('.gz') else open
        kmer_counter = 0
        
        try:
            with opener(fastq_file, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fastq'):
                    sequence = str(record.seq).upper()
                    
                    # 提取k-mer
                    for i in range(len(sequence) - self.config.kmer_size + 1):
                        kmer = sequence[i:i + self.config.kmer_size]
                        
                        # 跳过包含N的k-mer
                        if 'N' in kmer:
                            continue
                        
                        # 计算canonical k-mer
                        canonical_kmer = self._get_canonical_kmer(kmer)
                        
                        if canonical_kmer not in kmers:
                            kmer_counter += 1
                            # FASTQ k-mer不保存具体位置，用序号代替
                            position = KmerPosition(
                                seq_name=f"{sample_name}_kmer_{kmer_counter:06d}",
                                start_pos=0,
                                end_pos=0,
                                strand='+'
                            )
                            
                            kmers[canonical_kmer] = KmerInfo(
                                sequence=canonical_kmer,
                                positions=[position],
                                total_count=1,
                                source_type='fastq',
                                sample_name=sample_name
                            )
                        else:
                            kmers[canonical_kmer].total_count += 1
        
        except Exception as e:
            self.logger.error(f"Error processing FASTQ file {fastq_file}: {e}")
            raise
        
        self.logger.info(f"Extracted {len(kmers)} unique k-mers from {fastq_file}")
        return kmers
    
    def _get_canonical_kmer(self, kmer: str) -> str:
        """获取k-mer的canonical形式（字典序较小的正向或反向互补）"""
        if not self.config.canonical_form:
            return kmer
        
        reverse_complement = self._reverse_complement(kmer)
        return min(kmer, reverse_complement)
    
    def _reverse_complement(self, seq: str) -> str:
        """计算反向互补序列"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in seq[::-1])
    
    def merge_kmer_info(self, kmer_dicts: List[Dict[str, KmerInfo]]) -> Dict[str, KmerInfo]:
        """合并多个k-mer信息字典"""
        merged = {}
        
        for kmer_dict in kmer_dicts:
            for kmer, info in kmer_dict.items():
                if kmer not in merged:
                    merged[kmer] = KmerInfo(
                        sequence=kmer,
                        positions=[],
                        total_count=0,
                        source_type=info.source_type,
                        sample_name=info.sample_name
                    )
                
                merged[kmer].positions.extend(info.positions)
                merged[kmer].total_count += info.total_count
        
        return merged