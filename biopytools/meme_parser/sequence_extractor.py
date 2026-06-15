"""
MEME序列提取模块|MEME Sequence Extractor Module
"""

from Bio import SeqIO
from pathlib import Path
import pandas as pd


class MemeSequenceExtractor:
    """MEME序列提取器|MEME Sequence Extractor"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def extract_motif_sequences(self, motif_df: pd.DataFrame, input_fasta: str) -> dict:
        """
        提取motif序列|Extract motif sequences

        Args:
            motif_df: motif信息DataFrame|Motif information DataFrame
            input_fasta: 输入FASTA文件路径|Input FASTA file path

        Returns:
            dict: motif序列字典{(seq_id, motif_id): sequence}|Motif sequence dictionary
        """
        self.logger.info("提取motif序列|Extracting motif sequences")

        # 读取输入FASTA|Read input FASTA
        seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))

        motif_seqs = {}

        # 遍历每个motif位点|Iterate through each motif site
        for idx, row in motif_df.iterrows():
            seq_id = row['sequence_id']
            motif_id = row['motif_id']
            start = row['start']
            end = row['end']

            # 检查序列是否存在|Check if sequence exists
            if seq_id not in seq_dict:
                self.logger.warning(f"序列未找到|Sequence not found: {seq_id}")
                continue

            # 提取motif序列|Extract motif sequence
            full_seq = seq_dict[seq_id]
            motif_seq = full_seq.seq[start-1:end]  # FASTA是1-based，Python是0-based

            # 使用(seq_id, motif_id)作为key|Use (seq_id, motif_id) as key
            key = (seq_id, motif_id)
            motif_seqs[key] = str(motif_seq)

        self.logger.info(f"提取了{len(motif_seqs)}条motif序列|Extracted {len(motif_seqs)} motif sequences")

        return motif_seqs

    def write_motif_fasta(self, motif_seqs: dict, output_file: str):
        """
        写入motif序列FASTA文件|Write motif sequences to FASTA file

        Args:
            motif_seqs: motif序列字典|Motif sequence dictionary
            output_file: 输出文件路径|Output file path
        """
        if not motif_seqs:
            self.logger.warning("没有motif序列需要写入|No motif sequences to write")
            return

        output_path = Path(self.config.output_dir) / output_file

        with open(output_path, 'w') as f:
            for key, seq in motif_seqs.items():
                # key是(seq_id, motif_id) tuple
                if isinstance(key, tuple):
                    seq_id, motif_id = key
                    header = f"{seq_id}_{motif_id}"
                else:
                    header = key
                f.write(f">{header}\n{seq}\n")

        self.logger.info(f"写入{len(motif_seqs)}条motif序列到|Wrote {len(motif_seqs)} motif sequences to {output_file}")
