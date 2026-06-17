"""
FASTA序列ID重命名核心处理逻辑|FASTA ID Renamer Core Processing Logic
"""

from pathlib import Path
from Bio import SeqIO
from .config import FastaIDRenamerConfig
from .utils import format_number_with_padding


class FastaIDRenamer:
    """FASTA序列ID重命名器|FASTA Sequence ID Renamer"""

    def __init__(self, config, logger):
        """初始化重命名器|Initialize renamer

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # ID映射表|ID mapping table
        self.id_mapping = {}

        # 序列列表|Sequence list
        self.sequences = []

    def process(self):
        """处理FASTA文件|Process FASTA file"""
        self.logger.info("开始处理FASTA文件|Starting FASTA file processing")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")

        # 读取序列|Read sequences
        self._read_sequences()

        # 重命名序列|Rename sequences
        self._rename_sequences()

        # 写入输出文件|Write output file
        self._write_output()

        # 保存ID映射|Save ID mapping
        if self.config.save_mapping:
            self._save_id_mapping()

        # 输出统计信息|Output statistics
        self._log_statistics()

        self.logger.info("处理完成|Processing completed")

    def _read_sequences(self):
        """读取FASTA序列|Read FASTA sequences"""
        self.logger.info("读取序列|Reading sequences")

        try:
            record_count = 0
            for record in SeqIO.parse(self.config.input_file, "fasta"):
                # 解析序列头|Parse sequence header
                seq_id = record.id
                description = record.description

                # 存储序列信息|Store sequence information
                self.sequences.append({
                    'record': record,
                    'original_id': seq_id,
                    'description': description
                })

                record_count += 1

            self.logger.info(f"共读取 {record_count} 条序列|Total {record_count} sequences read")

        except Exception as e:
            self.logger.error(f"读取FASTA文件失败|Failed to read FASTA file: {e}")
            raise

    def _rename_sequences(self):
        """重命名序列|Rename sequences"""
        self.logger.info("重命名序列|Renaming sequences")

        # 按顺序重命名所有序列|Rename all sequences in order
        for idx, seq_info in enumerate(self.sequences, start=1):
            # 生成新ID|Generate new ID
            num_str = format_number_with_padding(idx, self.config.padding_width) \
                if self.config.use_zero_padding else str(idx)

            new_id = f"{self.config.prefix}{num_str}"
            seq_info['new_id'] = new_id
            self.id_mapping[seq_info['original_id']] = new_id

        self.logger.info(f"已重命名 {len(self.sequences)} 条序列|Renamed {len(self.sequences)} sequences")

    def _write_output(self):
        """写入输出文件|Write output file"""
        self.logger.info(f"写入输出文件|Writing output file: {self.config.output_file}")

        try:
            # 确定要写入的序列数量|Determine number of sequences to write
            if self.config.chr_count > 0:
                write_count = min(self.config.chr_count, len(self.sequences))
                self.logger.info(f"只输出前 {write_count} 条序列|Only outputting first {write_count} sequences")
            else:
                write_count = len(self.sequences)

            with open(self.config.output_file, 'w') as output_handle:
                for i in range(write_count):
                    seq_info = self.sequences[i]
                    # 修改序列ID|Modify sequence ID
                    record = seq_info['record']
                    record.id = seq_info['new_id']
                    record.description = seq_info['new_id']

                    # 写入序列|Write sequence
                    SeqIO.write(record, output_handle, "fasta")

            self.logger.info(f"成功写入 {write_count} 条序列|Successfully wrote {write_count} sequences")

        except Exception as e:
            self.logger.error(f"写入输出文件失败|Failed to write output file: {e}")
            raise

    def _save_id_mapping(self):
        """保存ID映射文件|Save ID mapping file"""
        mapping_file = self.config.mapping_file
        self.logger.info(f"保存ID映射文件|Saving ID mapping file: {mapping_file}")

        try:
            with open(mapping_file, 'w') as f:
                f.write("# Original_ID\tNew_ID\tDescription\n")
                for seq_info in self.sequences:
                    f.write(f"{seq_info['original_id']}\t{seq_info['new_id']}\t{seq_info['description']}\n")

            self.logger.info(f"ID映射文件已保存|ID mapping file saved: {mapping_file}")

        except Exception as e:
            self.logger.error(f"保存ID映射文件失败|Failed to save ID mapping file: {e}")
            raise

    def _log_statistics(self):
        """输出统计信息|Output statistics"""
        self.logger.info("=" * 60)
        self.logger.info("重命名统计|Renaming Statistics:")
        self.logger.info(f"  总序列数|Total sequences: {len(self.sequences)}")

        if self.config.chr_count > 0:
            chr_count = min(self.config.chr_count, len(self.sequences))
            self.logger.info(f"  输出序列数|Output sequences (first {self.config.chr_count}): {chr_count}")
            other_count = len(self.sequences) - chr_count
            self.logger.info(f"  未输出序列|Skipped sequences: {other_count}")

        self.logger.info("=" * 60)
