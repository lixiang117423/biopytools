"""
FASTA ID分割处理模块 | FASTA ID Splitting Processing Module
"""

import os
from pathlib import Path
from typing import Generator, Tuple
from .utils import FastaParser

class FastaProcessor:
    """FASTA文件处理器 | FASTA File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.parser = FastaParser(logger)
        self.processed_count = 0
        self.skipped_count = 0
    
    def read_fasta_sequences(self, file_path: str) -> Generator[Tuple[str, str], None, None]:
        """读取FASTA序列 | Read FASTA sequences"""
        current_header = ""
        current_sequence = ""
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # 如果有之前的序列，先输出 | If there's a previous sequence, yield it first
                        if current_header:
                            yield (current_header, current_sequence)
                        
                        current_header = line
                        current_sequence = ""
                    else:
                        current_sequence += line
                
                # 输出最后一个序列 | Yield the last sequence
                if current_header:
                    yield (current_header, current_sequence)
                    
        except Exception as e:
            self.logger.error(f"❌ 读取FASTA文件失败 | Failed to read FASTA file: {e}")
            raise
    
    def process_header(self, header: str) -> str:
        """处理序列名称行 | Process sequence header"""
        try:
            # 检查是否跳过空行 | Check if skip empty lines
            if self.config.skip_empty and not header.strip():
                self.skipped_count += 1
                return header
            
            # 分割并提取指定位置的元素 | Split and extract element at specified position
            new_id = self.parser.split_header(header, self.config.delimiter, self.config.position)
            
            # 添加 > 符号 | Add > symbol
            processed_header = f">{new_id}"
            
            self.processed_count += 1
            
            return processed_header
            
        except Exception as e:
            self.logger.warning(f"⚠️ 处理序列名称失败，保留原格式 | Failed to process header, keeping original: {header[:50]}... Error: {e}")
            self.skipped_count += 1
            return header
    
    def create_backup(self):
        """创建原文件备份 | Create backup of original file"""
        if self.config.keep_original:
            backup_file = f"{self.config.input_file}.backup"
            try:
                import shutil
                shutil.copy2(self.config.input_file, backup_file)
                self.logger.info(f"💾 备份文件已创建 | Backup file created: {backup_file}")
            except Exception as e:
                self.logger.error(f"❌ 创建备份文件失败 | Failed to create backup: {e}")
    
    def split_fasta_ids(self):
        """执行FASTA ID分割 | Execute FASTA ID splitting"""
        self.logger.info(f"🚀 开始处理FASTA文件 | Starting to process FASTA file: {self.config.input_file}")
        
        # 获取文件统计信息 | Get file statistics
        stats = self.parser.get_fasta_stats(self.config.input_file)
        self.logger.info(f"📊 文件统计 | File statistics: {stats['total_sequences']} 序列 | sequences, {stats['file_size']} 字节 | bytes")
        
        # 自动检测分隔符 | Auto detect delimiter
        if self.config.delimiter == "auto":
            detected_delimiter = self.parser.detect_delimiter(stats['header_lines'])
            self.config.delimiter = detected_delimiter
        
        # 创建备份 | Create backup
        if self.config.keep_original:
            self.create_backup()
        
        # 处理文件 | Process file
        try:
            with open(self.config.output_file, 'w', encoding='utf-8') as out_f:
                for header, sequence in self.read_fasta_sequences(self.config.input_file):
                    # 处理序列名称 | Process header
                    new_header = self.process_header(header)
                    
                    # 写入输出文件 | Write to output file
                    out_f.write(f"{new_header}\n")
                    out_f.write(f"{sequence}\n")
                    
                    # 进度报告 | Progress report
                    if (self.processed_count + self.skipped_count) % 1000 == 0:
                        self.logger.info(f"📈 已处理 | Processed: {self.processed_count + self.skipped_count} 序列 | sequences")
            
            self.logger.info(f"✨ 处理完成 | Processing completed!")
            self.logger.info(f"✅ 成功处理 | Successfully processed: {self.processed_count} 序列 | sequences")
            self.logger.info(f"⚠️ 跳过 | Skipped: {self.skipped_count} 序列 | sequences")
            self.logger.info(f"📁 输出文件 | Output file: {self.config.output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 处理文件时发生错误 | Error occurred during file processing: {e}")
            raise
