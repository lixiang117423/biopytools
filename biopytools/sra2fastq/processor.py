"""
SRA文件处理模块 |SRA File Processing Module
"""

import os
from pathlib import Path
from typing import List
from .utils import CommandRunner

class SRAProcessor:
    """SRA文件处理器|SRA File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def _build_parallel_fastq_dump_cmd(self, sra_file: str) -> str:
        """构建parallel-fastq-dump命令|Build parallel-fastq-dump command"""
        cmd_parts = [self.config.tool_path]
        
        # SRA文件 (使用-s或--sra-id)|SRA file
        cmd_parts.append(f"--sra-id {sra_file}")
        
        # 线程数|Threads
        cmd_parts.append(f"--threads {self.config.threads}")
        
        # 输出目录|Output directory
        cmd_parts.append(f"--outdir {self.config.output_dir}")
        
        # 临时目录|Temporary directory
        if self.config.tmpdir:
            cmd_parts.append(f"--tmpdir {self.config.tmpdir}")
        
        # 以下参数会传递给底层的fastq-dump|Following params pass to underlying fastq-dump
        
        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            cmd_parts.append("--split-files")
        
        # 压缩输出|Compress output
        if self.config.compress:
            cmd_parts.append("--gzip")
        
        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            cmd_parts.append("--skip-technical")
        
        # 剪切adapters|Clip adapters
        if self.config.clip:
            cmd_parts.append("--clip")
        
        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            cmd_parts.append(f"--minReadLen {self.config.min_read_len}")
        
        return " ".join(cmd_parts)
    
    def _build_fastq_dump_cmd(self, sra_file: str) -> str:
        """构建fastq-dump命令 (备选方案)|Build fastq-dump command (fallback)"""
        cmd_parts = [self.config.tool_path]
        
        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            cmd_parts.append("--split-3")
        
        # 压缩输出|Compress output
        if self.config.compress:
            cmd_parts.append("--gzip")
        
        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            cmd_parts.append("--skip-technical")
        
        # 剪切adapters|Clip adapters
        if self.config.clip:
            cmd_parts.append("--clip")
        
        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            cmd_parts.append(f"--minReadLen {self.config.min_read_len}")
        
        # 输出目录|Output directory
        cmd_parts.append(f"--outdir {self.config.output_dir}")
        
        # 输入文件|Input file
        cmd_parts.append(sra_file)
        
        return " ".join(cmd_parts)
    
    def convert_single_file(self, sra_file: str) -> bool:
        """转换单个SRA文件|Convert single SRA file"""
        base_name = Path(sra_file).stem
        self.logger.info(f"{'='*60}")
        self.logger.info(f"处理文件|Processing file: {base_name}")
        self.logger.info(f"{'='*60}")
        
        # 根据工具类型构建命令|Build command based on tool type
        if self.config.use_parallel:
            cmd = self._build_parallel_fastq_dump_cmd(sra_file)
            tool_name = "parallel-fastq-dump (多线程加速)|(multi-threaded)"
        else:
            cmd = self._build_fastq_dump_cmd(sra_file)
            tool_name = "fastq-dump (单线程)|(single-threaded)"
        
        self.logger.info(f"使用工具|Using tool: {tool_name}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        
        # 执行转换|Execute conversion
        success = self.cmd_runner.run(
            cmd, 
            f"转换SRA文件|Converting SRA file: {base_name}"
        )
        
        if success:
            self.logger.info(f"完成|Completed: {base_name}")
        else:
            self.logger.error(f"失败|Failed: {base_name}")
        
        return success
    
    def convert_all_files(self) -> dict:
        """转换所有SRA文件|Convert all SRA files"""
        total = len(self.config.input_files)
        self.logger.info(f" 共找到 {total} 个SRA文件|Found {total} SRA files")
        
        results = {
            'success': [],
            'failed': [],
            'total': total
        }
        
        for idx, sra_file in enumerate(self.config.input_files, 1):
            self.logger.info(f" 总进度|Overall Progress: [{idx}/{total}]")
            
            if self.convert_single_file(sra_file):
                results['success'].append(sra_file)
            else:
                results['failed'].append(sra_file)
        
        return results
