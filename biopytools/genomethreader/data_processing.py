"""
GenomeThreader 数据处理模块 | GenomeThreader Data Processing Module
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List
from .utils import CommandRunner, get_fasta_stats

class SequenceProcessor:
    """序列文件处理器 | Sequence File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def validate_fasta_files(self) -> bool:
        """验证FASTA文件格式 | Validate FASTA file format"""
        self.logger.info("🧬 验证FASTA文件格式 | Validating FASTA file format")
        
        files_to_check = [
            (self.config.genomic_file, "基因组文件 | Genomic file"),
        ]
        
        if self.config.cdna_file:
            files_to_check.append((self.config.cdna_file, "cDNA文件 | cDNA file"))
        if self.config.protein_file:
            files_to_check.append((self.config.protein_file, "蛋白质文件 | Protein file"))
        if self.config.est_file:
            files_to_check.append((self.config.est_file, "EST文件 | EST file"))
        
        for file_path, file_desc in files_to_check:
            if not self._validate_single_fasta(file_path, file_desc):
                return False
                
        return True
    
    def _validate_single_fasta(self, file_path: str, file_desc: str) -> bool:
        """验证单个FASTA文件 | Validate single FASTA file"""
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"❌ {file_desc}格式错误，不是有效的FASTA文件 | Invalid FASTA format: {file_path}")
                    return False
            
            # 获取文件统计信息
            stats = get_fasta_stats(file_path, self.logger)
            if stats['sequences'] == 0:
                self.logger.error(f"❌ {file_desc}为空或无有效序列 | Empty or no valid sequences: {file_path}")
                return False
            
            self.logger.info(f"✅ {file_desc}验证通过 | Validation passed: {file_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 读取{file_desc}时出错 | Error reading file: {e}")
            return False
    
    def get_sequence_statistics(self) -> Dict[str, dict]:
        """获取所有序列文件的统计信息 | Get statistics for all sequence files"""
        self.logger.info("📊 收集序列文件统计信息 | Collecting sequence file statistics")
        
        statistics = {}
        
        # 基因组文件统计
        statistics['genomic'] = get_fasta_stats(self.config.genomic_file, self.logger)
        
        # 其他序列文件统计
        if self.config.cdna_file:
            statistics['cdna'] = get_fasta_stats(self.config.cdna_file, self.logger)
        if self.config.protein_file:
            statistics['protein'] = get_fasta_stats(self.config.protein_file, self.logger)
        if self.config.est_file:
            statistics['est'] = get_fasta_stats(self.config.est_file, self.logger)
        
        return statistics

class FileManager:
    """文件管理器 | File Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def setup_output_structure(self):
        """设置输出目录结构 | Setup output directory structure"""
        self.logger.info("📁 设置输出目录结构 | Setting up output directory structure")
        
        # 创建主输出目录
        self.config.output_path.mkdir(parents=True, exist_ok=True)
        
        # 创建子目录
        subdirs = ['alignments', 'predictions', 'intermediate', 'logs']
        for subdir in subdirs:
            (self.config.output_path / subdir).mkdir(exist_ok=True)
        
        self.logger.info(f"✅ 输出目录结构创建完成 | Output directory structure created: {self.config.output_dir}")
    
    def cleanup_temp_files(self):
        """清理临时文件 | Cleanup temporary files"""
        self.logger.info("🧹 清理临时文件 | Cleaning up temporary files")
        
        temp_patterns = ['*.tmp', '*.temp', '*.intermediate']
        temp_count = 0
        
        for pattern in temp_patterns:
            for temp_file in self.config.output_path.glob(pattern):
                try:
                    temp_file.unlink()
                    temp_count += 1
                except Exception as e:
                    self.logger.warning(f"⚠️ 无法删除临时文件 | Cannot delete temp file {temp_file}: {e}")
        
        if temp_count > 0:
            self.logger.info(f"✅ 清理了 {temp_count} 个临时文件 | Cleaned up {temp_count} temporary files")
        else:
            self.logger.info("✅ 没有发现需要清理的临时文件 | No temporary files to clean up")
