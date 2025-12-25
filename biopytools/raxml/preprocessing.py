"""
🌳 RAxML数据预处理模块 | RAxML Data Preprocessing Module
"""

import os
import shutil
from pathlib import Path
from typing import Dict, Tuple
from .utils import CommandRunner

class SequenceProcessor:
    """序列文件预处理器 | Sequence File Preprocessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def validate_sequence_format(self) -> bool:
        """验证序列文件格式 | Validate sequence file format"""
        self.logger.info("🔍 验证序列文件格式 | Validating sequence file format")
        
        try:
            with open(self.config.sequence_file, 'r') as f:
                first_line = f.readline().strip()
                
                # 检查是否是PHYLIP格式 | Check if it's PHYLIP format
                if first_line and len(first_line.split()) >= 2:
                    try:
                        num_taxa = int(first_line.split()[0])
                        seq_length = int(first_line.split()[1])
                        
                        self.logger.info(f"📊 检测到PHYLIP格式 | Detected PHYLIP format: {num_taxa} taxa, {seq_length} bp")
                        return True
                    except ValueError:
                        pass
                
                # 检查是否是FASTA格式 | Check if it's FASTA format
                if first_line.startswith('>'):
                    self.logger.info("📄 检测到FASTA格式 | Detected FASTA format")
                    self.logger.warning("⚠️ RAxML需要PHYLIP格式，请转换序列格式 | RAxML requires PHYLIP format, please convert sequence format")
                    return False
                
                self.logger.error("❌ 无法识别序列文件格式 | Cannot recognize sequence file format")
                return False
                
        except Exception as e:
            self.logger.error(f"❌ 读取序列文件失败 | Failed to read sequence file: {e}")
            return False
    
    def get_sequence_statistics(self) -> Dict[str, int]:
        """获取序列统计信息 | Get sequence statistics"""
        self.logger.info("📊 获取序列统计信息 | Getting sequence statistics")
        
        stats = {'taxa': 0, 'length': 0, 'variable_sites': 0}
        
        try:
            with open(self.config.sequence_file, 'r') as f:
                first_line = f.readline().strip()
                
                if first_line and len(first_line.split()) >= 2:
                    stats['taxa'] = int(first_line.split()[0])
                    stats['length'] = int(first_line.split()[1])
                    
                    self.logger.info(f"🧬 序列统计 | Sequence statistics:")
                    self.logger.info(f"  - 序列数量 | Number of taxa: {stats['taxa']}")
                    self.logger.info(f"  - 序列长度 | Sequence length: {stats['length']} bp")
                    
                    # 估算可变位点 | Estimate variable sites (rough estimation)
                    if stats['length'] > 0:
                        # 粗略估算：假设20-30%的位点是可变的 | Rough estimate: assume 20-30% sites are variable
                        estimated_variable = int(stats['length'] * 0.25)
                        stats['variable_sites'] = estimated_variable
                        self.logger.info(f"  - 估计可变位点 | Estimated variable sites: ~{estimated_variable}")
        
        except Exception as e:
            self.logger.warning(f"⚠️ 无法获取序列统计信息 | Cannot get sequence statistics: {e}")
        
        return stats
    
    def prepare_working_directory(self):
        """准备工作目录 | Prepare working directory"""
        self.logger.info("📁 准备工作目录 | Preparing working directory")
        
        # 确保输出目录存在 | Ensure output directory exists
        self.config.output_path.mkdir(parents=True, exist_ok=True)
        
        # 如果序列文件不在输出目录中，复制一份 | Copy sequence file if not in output directory
        seq_file_path = Path(self.config.sequence_file)
        if not str(seq_file_path).startswith(str(self.config.output_path)):
            target_path = self.config.output_path / seq_file_path.name
            if not target_path.exists():
                shutil.copy2(self.config.sequence_file, target_path)
                self.logger.info(f"📋 复制序列文件到工作目录 | Copied sequence file to working directory: {target_path.name}")
        
        # 复制其他相关文件 | Copy other related files
        for file_attr in ['starting_tree', 'constraint_tree']:
            file_path = getattr(self.config, file_attr)
            if file_path and os.path.exists(file_path):
                source_path = Path(file_path)
                target_path = self.config.output_path / source_path.name
                if not target_path.exists():
                    shutil.copy2(file_path, target_path)
                    self.logger.info(f"📋 复制{file_attr}到工作目录 | Copied {file_attr} to working directory: {target_path.name}")
