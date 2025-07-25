"""
基因组组装数据处理模块 | Genome Assembly Data Processing Module
"""

import os
from typing import List, Dict
from .utils import CommandRunner

class ReadProcessor:
    """测序数据预处理器 | Sequencing Read Preprocessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def validate_reads(self) -> bool:
        """验证reads文件格式和质量 | Validate reads file format and quality"""
        self.logger.info("验证reads文件 | Validating reads files")
        
        reads_files = [self.config.hifi_reads]
        if self.config.ont_reads:
            reads_files.append(self.config.ont_reads)
        
        if self.config.trio_mode:
            reads_files.extend([self.config.parent1_reads, self.config.parent2_reads, self.config.child_reads])
        
        for reads_file in reads_files:
            if not self._check_file_format(reads_file):
                return False
        
        return True
    
    def _check_file_format(self, reads_file: str) -> bool:
        """检查单个reads文件格式 | Check individual reads file format"""
        self.logger.info(f"检查文件格式 | Checking file format: {reads_file}")
        
        # 检查文件扩展名 | Check file extension
        valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fasta', '.fa', '.fasta.gz', '.fa.gz']
        if not any(reads_file.endswith(ext) for ext in valid_extensions):
            self.logger.error(f"不支持的文件格式 | Unsupported file format: {reads_file}")
            return False
        
        # 检查文件大小 | Check file size
        file_size = os.path.getsize(reads_file)
        if file_size == 0:
            self.logger.error(f"文件为空 | File is empty: {reads_file}")
            return False
        
        self.logger.info(f"文件格式验证通过 | File format validation passed: {reads_file} ({file_size/1024/1024/1024:.2f} GB)")
        return True
