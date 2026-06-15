"""
Resistify序列提取模块|Resistify Sequence Extractor Module
"""

from Bio import SeqIO
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd


class SequenceExtractor:
    """Resistify序列提取器|Resistify Sequence Extractor"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def extract_nlr_sequences(self, filtered_df: pd.DataFrame,
                              resistify_output_dir: Optional[str] = None) -> Dict[str, str]:
        """
        从nlr.fasta提取序列|Extract sequences from nlr.fasta

        Args:
            filtered_df: 过滤后的DataFrame|Filtered DataFrame
            resistify_output_dir: Resistify输出目录(默认使用config)|Resistify output directory

        Returns:
            Dict[str, str]: 序列字典{sequence_id: sequence}|Sequence dictionary
        """
        self.logger.info("提取NLR序列|Extracting NLR sequences")

        if resistify_output_dir is None:
            resistify_output_dir = self.config.resistify_output_dir

        nlr_fasta = Path(resistify_output_dir) / 'nlr.fasta'

        if not nlr_fasta.exists():
            self.logger.warning(f"NLR序列文件不存在|NLR fasta file not found: {nlr_fasta}")
            return {}

        seq_dict = SeqIO.to_dict(SeqIO.parse(str(nlr_fasta), 'fasta'))
        target_ids = set(filtered_df['Sequence'].values)
        extracted = {seq_id: str(seq_dict[seq_id].seq) for seq_id in target_ids if seq_id in seq_dict}

        self.logger.info(f"提取了{len(extracted)}条NLR序列|Extracted {len(extracted)} NLR sequences")
        return extracted

    def extract_nbarc_sequences(self, filtered_df: pd.DataFrame,
                                resistify_output_dir: Optional[str] = None) -> Dict[str, str]:
        """
        从nbarc.fasta提取序列|Extract sequences from nbarc.fasta

        Args:
            filtered_df: 过滤后的DataFrame|Filtered DataFrame
            resistify_output_dir: Resistify输出目录(默认使用config)|Resistify output directory

        Returns:
            Dict[str, str]: 序列字典{sequence_id: sequence}|Sequence dictionary
        """
        self.logger.info("提取NB-ARC序列|Extracting NB-ARC sequences")

        if resistify_output_dir is None:
            resistify_output_dir = self.config.resistify_output_dir

        nbarc_fasta = Path(resistify_output_dir) / 'nbarc.fasta'

        if not nbarc_fasta.exists():
            self.logger.warning(f"NB-ARC序列文件不存在|NB-ARC fasta file not found: {nbarc_fasta}")
            return {}

        seq_dict = SeqIO.to_dict(SeqIO.parse(str(nbarc_fasta), 'fasta'))
        target_ids = set(filtered_df['Sequence'].values)
        extracted = {seq_id: str(seq_dict[seq_id].seq) for seq_id in target_ids if seq_id in seq_dict}

        self.logger.info(f"提取了{len(extracted)}条NB-ARC序列|Extracted {len(extracted)} NB-ARC sequences")
        return extracted

    def write_fasta(self, sequences: Dict[str, str], output_file: str,
                    output_dir: Optional[str] = None):
        """
        写入FASTA文件|Write sequences to FASTA file

        Args:
            sequences: 序列字典|Sequence dictionary
            output_file: 输出文件名|Output file name
            output_dir: 输出目录(默认使用config)|Output directory
        """
        if not sequences:
            self.logger.warning("没有序列需要写入|No sequences to write")
            return

        if output_dir is None:
            output_dir = self.config.output_dir

        output_path = Path(output_dir) / output_file

        with open(output_path, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n{seq}\n")

        self.logger.info(f"写入{len(sequences)}条序列到|Wrote {len(sequences)} sequences to {output_file}")

    def extract_and_write(self, filtered_df: pd.DataFrame,
                          resistify_output_dir: Optional[str] = None,
                          output_dir: Optional[str] = None):
        """
        根据配置提取并写入序列|Extract and write sequences based on config

        Args:
            filtered_df: 过滤后的DataFrame|Filtered DataFrame
            resistify_output_dir: Resistify输出目录(默认使用config)|Resistify output directory
            output_dir: 输出目录(默认使用config)|Output directory
        """
        if resistify_output_dir is None:
            resistify_output_dir = self.config.resistify_output_dir
        if output_dir is None:
            output_dir = self.config.output_dir

        if self.config.extract_nlr_sequences:
            nlr_seqs = self.extract_nlr_sequences(filtered_df, resistify_output_dir)
            if nlr_seqs:
                output_file = f"{self.config.output_prefix}_nlr.fasta"
                self.write_fasta(nlr_seqs, output_file, output_dir)

        if self.config.extract_nbarc_sequences:
            nbarc_seqs = self.extract_nbarc_sequences(filtered_df, resistify_output_dir)
            if nbarc_seqs:
                output_file = f"{self.config.output_prefix}_nbarc.fasta"
                self.write_fasta(nbarc_seqs, output_file, output_dir)
