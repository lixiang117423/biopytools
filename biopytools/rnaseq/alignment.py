"""
RNA-seq比对模块 | RNA-seq Alignment Module
"""

import os
from .utils import CommandRunner, FileValidator

class HISAT2Indexer:
    """HISAT2索引构建器 | HISAT2 Index Builder"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def build_hisat2_index(self) -> str:
        """构建HISAT2基因组索引 | Build HISAT2 genome index"""
        genome_path = self.config.genome_file
        threads = self.config.threads
        
        genome_dir = os.path.dirname(genome_path)
        genome_name = os.path.splitext(os.path.basename(genome_path))[0]
        index_prefix = os.path.join(genome_dir, f"{genome_name}.hisat2.index")

        # 检查索引是否已存在 | Check if index already exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2索引已存在 | HISAT2 index already exists: {index_prefix}")
            return index_prefix

        cmd = f"hisat2-build -p {threads} {genome_path} {index_prefix}"
        self.cmd_runner.run(cmd, "构建HISAT2索引 | Building HISAT2 index")
        
        self.logger.info(f"HISAT2索引构建完成 | HISAT2 index building completed: {index_prefix}")
        return index_prefix

class HISAT2Aligner:
    """HISAT2比对器 | HISAT2 Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def run_hisat2_mapping(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str) -> bool:
        """运行HISAT2比对和排序 | Run HISAT2 alignment and sorting"""
        threads = self.config.threads
        
        # 确保输出目录存在 | Ensure output directory exists
        output_dir = os.path.dirname(output_bam)
        os.makedirs(output_dir, exist_ok=True)

        # 检查BAM文件是否已存在 | Check if BAM file already exists
        if self.file_validator.check_file_exists(output_bam, "BAM文件 | BAM file"):
            return True

        cmd = (
            f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads} | "
            f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
        )

        success = self.cmd_runner.run(cmd, f"HISAT2比对和排序 | HISAT2 alignment and sorting -> {output_bam}")
        
        if success:
            self.logger.info(f"比对完成 | Alignment completed: {output_bam}")
        
        return success
