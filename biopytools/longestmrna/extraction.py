"""
最长转录本提取核心模块|Longest mRNA Extraction Core Module
"""

import os
import tempfile
from typing import Dict, Any
from .utils import CommandRunner, TempFileManager, build_conda_command

class SequenceExtractor:
    """序列提取器|Sequence Extractor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        # 临时文件重定向到 output/tmp,避免超算系统 /tmp 爆满|Redirect temp files to output/tmp to avoid filling system /tmp
        tmp_dir = os.path.join(os.path.dirname(config.output_file), "tmp")
        self.temp_manager = TempFileManager(logger, base_dir=tmp_dir)

    def _extract_sequences(self, longest_transcripts: Dict[str, Dict[str, Any]],
                           output_file: str, gffread_flag: str, seq_type: str) -> bool:
        """提取指定类型序列|Extract sequences of specified type

        Args:
            longest_transcripts: 最长转录本字典|Longest transcripts dict
            output_file: 输出文件路径|Output file path
            gffread_flag: gffread输出参数，'-y'蛋白/'-x'CDS|gffread output flag
            seq_type: 序列类型描述|Sequence type description
        """
        try:
            with self.temp_manager.create_temp_file(mode='w+', delete=False, encoding='utf-8') as id_file:
                for gene_id, transcript_info in longest_transcripts.items():
                    transcript_id = transcript_info.get('id', '')
                    if transcript_id:
                        id_file.write(f"{transcript_id}\n")
                temp_id_path = id_file.name

            with self.temp_manager.create_temp_file(mode='w+', delete=False, suffix='.fa') as seq_file:
                gffread_args = [
                    self.config.gff3_file,
                    '-g', self.config.genome_file,
                    gffread_flag, seq_file.name
                ]
                gffread_cmd = build_conda_command(self.config.gffread_path, gffread_args)
                self.cmd_runner.run(gffread_cmd, f"使用gffread生成{seq_type}序列|Generate {seq_type} sequences using gffread")

            seqkit_args = ['grep', '-f', temp_id_path, seq_file.name, '-o', output_file]
            seqkit_cmd = build_conda_command(self.config.seqkit_path, seqkit_args)
            result = self.cmd_runner.run(seqkit_cmd, f"使用seqkit筛选最长转录本{seq_type}序列|Filter longest transcript {seq_type} sequences using seqkit")

            if result.returncode == 0:
                self.logger.info(f"成功提取 {len(longest_transcripts)} 个最长转录本{seq_type}序列|Successfully extracted {len(longest_transcripts)} longest transcript {seq_type} sequences")
                self.logger.info(f"输出文件|Output file: {output_file}")
                return True
            else:
                self.logger.error(f"{seq_type}序列提取失败|{seq_type} sequence extraction failed")
                return False

        except Exception as e:
            self.logger.error(f"{seq_type}序列提取过程中发生错误|Error during {seq_type} sequence extraction: {e}")
            return False
        finally:
            self.temp_manager.cleanup()

    def extract_protein_sequences(self, longest_transcripts: Dict[str, Dict[str, Any]]) -> bool:
        """提取蛋白质序列|Extract protein sequences"""
        return self._extract_sequences(longest_transcripts, self.config.output_file, '-y', '蛋白')

    def extract_cds_sequences(self, longest_transcripts: Dict[str, Dict[str, Any]]) -> bool:
        """提取CDS核苷酸序列|Extract CDS nucleotide sequences"""
        return self._extract_sequences(longest_transcripts, self.config.cds_output_file, '-x', 'CDS')
