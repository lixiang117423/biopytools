"""
序列提取模块|Sequence Extraction Module
"""

from Bio import SeqIO
from pathlib import Path
from typing import List, Dict
from .parser import DomainHit


class SequenceExtractor:
    """序列提取器|Sequence Extractor"""

    def __init__(self, logger, protein_fastas):
        self.logger = logger
        self.protein_fastas = protein_fastas if isinstance(protein_fastas, list) else [protein_fastas]
        self.sequence_sources = {}
        self.sequences = self._load_sequences()

    def _load_sequences(self) -> Dict[str, str]:
        """
        加载蛋白序列|Load protein sequences

        Returns:
            Dict[str, str]: 序列ID到序列的映射|Mapping from sequence ID to sequence
        """
        sequences = {}
        total = len(self.protein_fastas)

        for i, fasta_path in enumerate(self.protein_fastas, 1):
            self.logger.info(f"加载蛋白序列文件|Loading protein sequences [{i}/{total}]: {fasta_path}")
            source_name = Path(fasta_path).stem

            try:
                count_before = len(sequences)
                for record in SeqIO.parse(fasta_path, "fasta"):
                    sequences[record.id] = str(record.seq)
                    self.sequence_sources[record.id] = source_name

                loaded = len(sequences) - count_before
                self.logger.info(f"加载{loaded}条序列，累计{len(sequences)}条|Loaded {loaded} sequences, total {len(sequences)}")
            except Exception as e:
                self.logger.error(f"序列加载失败|Failed to load sequences: {fasta_path}: {e}")
                raise

        return sequences

    def extract_protein_sequences(self, hits: List[DomainHit], output_file: str):
        """
        提取匹配的蛋白序列|Extract matched protein sequences

        Args:
            hits: Domain命中列表|List of domain hits
            output_file: 输出文件路径|Output file path
        """
        self.logger.info(f"开始提取蛋白序列|Extracting protein sequences to: {output_file}")

        # 获取唯一的目标基因|Get unique target genes
        unique_targets = set(hit.target_name for hit in hits)

        extracted = 0
        with open(output_file, 'w') as f_out:
            for target_id in unique_targets:
                if target_id in self.sequences:
                    f_out.write(f">{target_id}\n")
                    # 每行60个字符|60 characters per line
                    seq = self.sequences[target_id]
                    for i in range(0, len(seq), 60):
                        f_out.write(seq[i:i+60] + "\n")
                    extracted += 1
                else:
                    self.logger.warning(f"未找到序列|Sequence not found: {target_id}")

        self.logger.info(f"蛋白序列提取完成|Protein sequences extracted: {extracted}/{len(unique_targets)}")

    def extract_domain_sequences(self, hits: List[DomainHit], output_file: str):
        """
        提取domain序列|Extract domain sequences

        Args:
            hits: Domain命中列表|List of domain hits
            output_file: 输出文件路径|Output file path
        """
        self.logger.info(f"开始提取domain序列|Extracting domain sequences to: {output_file}")

        extracted = 0
        with open(output_file, 'w') as f_out:
            for hit in hits:
                if hit.target_name in self.sequences:
                    seq = self.sequences[hit.target_name]
                    seq_len = len(seq)

                    # 检查坐标是否有效|Check if coordinates are valid
                    if hit.ali_from < 1 or hit.ali_to > seq_len or hit.ali_from > hit.ali_to:
                        self.logger.warning(
                            f"无效坐标|Invalid coordinates for {hit.target_name}: "
                            f"{hit.ali_from}-{hit.ali_to} (seq length: {seq_len})"
                        )
                        continue

                    # 提取domain序列(Python从0开始，所以需要-1)|Extract domain sequence
                    domain_seq = seq[hit.ali_from - 1:hit.ali_to]

                    # 写入序列|Write sequence
                    domain_id = f"{hit.target_name}_domain{hit.domain_number}"
                    f_out.write(f">{domain_id} {hit.query_name}|{hit.ali_from}-{hit.ali_to}|E-value:{hit.domain_evalue:.2e}\n")

                    for i in range(0, len(domain_seq), 60):
                        f_out.write(domain_seq[i:i+60] + "\n")

                    extracted += 1
                else:
                    self.logger.warning(f"未找到序列|Sequence not found: {hit.target_name}")

        self.logger.info(f"Domain序列提取完成|Domain sequences extracted: {extracted}/{len(hits)}")
