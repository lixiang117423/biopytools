"""
最长转录本提取核心模块 | Longest mRNA Extraction Core Module
"""

import os
import tempfile
from typing import Dict, Any
from .utils import CommandRunner, TempFileManager

class SequenceExtractor:
    """序列提取器 | Sequence Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.temp_manager = TempFileManager(logger)
    
    def extract_protein_sequences(self, longest_transcripts: Dict[str, Dict[str, Any]]) -> bool:
        """提取蛋白质序列 | Extract protein sequences"""
        try:
            # 创建转录本ID临时文件 | Create temporary file for transcript IDs
            with self.temp_manager.create_temp_file(mode='w+', delete=False, encoding='utf-8') as id_file:
                for gene_id, transcript_info in longest_transcripts.items():
                    transcript_id = transcript_info.get('id', '')
                    if transcript_id:
                        id_file.write(f"{transcript_id}\n")
                temp_id_path = id_file.name
            
            # 使用gffread生成蛋白质序列 | Use gffread to generate protein sequences
            with self.temp_manager.create_temp_file(mode='w+', delete=False, suffix='.fa') as protein_file:
                gffread_cmd = f'gffread "{self.config.gff3_file}" -g "{self.config.genome_file}" -y "{protein_file.name}"'
                self.cmd_runner.run(gffread_cmd, "使用gffread生成蛋白质序列 | Generate protein sequences using gffread")
            
            # 使用seqkit筛选最长转录本序列 | Use seqkit to filter longest transcript sequences
            seqkit_cmd = f'seqkit grep -f "{temp_id_path}" "{protein_file.name}" -o "{self.config.output_file}"'
            result = self.cmd_runner.run(seqkit_cmd, "使用seqkit筛选最长转录本序列 | Filter longest transcript sequences using seqkit")
            
            if result.returncode == 0:
                self.logger.info(f"✓ 成功提取 {len(longest_transcripts)} 个最长转录本序列 | Successfully extracted {len(longest_transcripts)} longest transcript sequences")
                self.logger.info(f"✓ 输出文件 | Output file: {self.config.output_file}")
                return True
            else:
                self.logger.error("✗ 序列提取失败 | Sequence extraction failed")
                return False
                
        except Exception as e:
            self.logger.error(f"✗ 序列提取过程中发生错误 | Error during sequence extraction: {e}")
            return False
        finally:
            # 清理临时文件 | Cleanup temporary files
            self.temp_manager.cleanup()
