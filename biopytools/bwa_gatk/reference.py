"""
📚 参考基因组索引管理模块 | Reference Genome Index Management Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class ReferenceIndexer:
    """参考基因组索引管理器 | Reference Genome Index Manager"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.reference = Path(config.reference)
    
    def check_and_build_indices(self):
        """检查并构建所有必需的索引 | Check and build all required indices"""
        self.logger.info("=" * 80)
        self.logger.info("📚 检查参考基因组索引 | Checking reference genome indices")
        self.logger.info("=" * 80)
        
        # BWA索引 | BWA indices
        self._check_bwa_index()
        
        # SAMtools索引 | SAMtools index
        self._check_samtools_index()
        
        # GATK字典 | GATK dictionary
        self._check_gatk_dict()
        
        self.logger.info("✅ 参考基因组索引检查完成 | Reference genome indices check completed")
    
    def _check_bwa_index(self):
        """检查BWA索引 | Check BWA index"""
        self.logger.info("🔍 检查BWA索引 | Checking BWA index")
        
        bwa_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        missing_indices = []
        
        for ext in bwa_extensions:
            index_file = Path(str(self.reference) + ext)
            if not index_file.exists():
                missing_indices.append(ext)
        
        if missing_indices:
            self.logger.warning(f"⚠️  缺少BWA索引文件 | Missing BWA index files: {', '.join(missing_indices)}")
            self._build_bwa_index()
        else:
            self.logger.info("  ✅ BWA索引完整 | BWA index complete")
    
    def _build_bwa_index(self):
        """构建BWA索引 | Build BWA index"""
        self.logger.info("🔨 构建BWA索引 | Building BWA index")
        cmd = f"{self.config.bwa_path} index {self.reference}"
        self.cmd_runner.run(cmd, "构建BWA索引 | Building BWA index")
    
    def _check_samtools_index(self):
        """检查SAMtools索引 | Check SAMtools index"""
        self.logger.info("🔍 检查SAMtools fai索引 | Checking SAMtools fai index")
        
        fai_file = Path(str(self.reference) + '.fai')
        
        if not fai_file.exists():
            self.logger.warning("⚠️  缺少.fai索引文件 | Missing .fai index file")
            self._build_samtools_index()
        else:
            self.logger.info("  ✅ SAMtools索引存在 | SAMtools index exists")
    
    def _build_samtools_index(self):
        """构建SAMtools索引 | Build SAMtools index"""
        self.logger.info("🔨 构建SAMtools索引 | Building SAMtools index")
        cmd = f"{self.config.samtools_path} faidx {self.reference}"
        self.cmd_runner.run(cmd, "构建SAMtools索引 | Building SAMtools index")
    
    def _check_gatk_dict(self):
        """检查GATK字典 | Check GATK dictionary"""
        self.logger.info("🔍 检查GATK字典文件 | Checking GATK dictionary")
        
        dict_file = self.reference.with_suffix('.dict')
        
        if not dict_file.exists():
            self.logger.warning("⚠️  缺少.dict字典文件 | Missing .dict dictionary file")
            self._build_gatk_dict()
        else:
            self.logger.info("  ✅ GATK字典存在 | GATK dictionary exists")
    
    def _build_gatk_dict(self):
        """构建GATK字典 | Build GATK dictionary"""
        self.logger.info("🔨 构建GATK字典 | Building GATK dictionary")
        dict_file = self.reference.with_suffix('.dict')
        cmd = f"{self.config.gatk_path} CreateSequenceDictionary -R {self.reference} -O {dict_file}"
        self.cmd_runner.run(cmd, "构建GATK字典 | Building GATK dictionary")
