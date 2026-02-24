"""
基因组索引管理模块|Genome Index Management Module
"""

import os
from pathlib import Path

class GenomeIndexer:
    """基因组索引器|Genome Indexer"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def check_and_build_index(self):
        """检查并构建基因组索引|Check and build genome index"""
        self.logger.info("检查基因组索引|Checking genome index")
        
        # BWA索引文件后缀|BWA index file extensions
        index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        genome_path = Path(self.config.genome)
        
        # 检查所有索引文件是否存在|Check if all index files exist
        missing_indices = []
        for ext in index_extensions:
            index_file = genome_path.parent / f"{genome_path.name}{ext}"
            if not index_file.exists():
                missing_indices.append(ext)
        
        if missing_indices:
            self.logger.info(f"缺少索引文件|Missing index files: {', '.join(missing_indices)}")
            self.logger.info("开始构建BWA索引|Building BWA index")
            return self.build_index()
        else:
            self.logger.info("基因组索引已存在|Genome index already exists")
            return True
    
    def build_index(self):
        """构建BWA索引|Build BWA index"""
        cmd = f"{self.config.bwa_path} index {self.config.genome}"
        
        success = self.cmd_runner.run(
            cmd,
            "构建BWA索引|Building BWA index"
        )
        
        if success:
            self.logger.info("BWA索引构建成功|BWA index built successfully")
        else:
            self.logger.error("BWA索引构建失败|BWA index building failed")
        
        return success
