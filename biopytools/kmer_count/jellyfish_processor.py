"""
🐙 Jellyfish处理模块 | Jellyfish Processing Module
"""

from pathlib import Path
from typing import List

class JellyfishProcessor:
    """🐙 Jellyfish处理器 | Jellyfish Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def count_kmers(self, sample_name: str, fastq_files: List[str]) -> Path:
        """🧮 使用Jellyfish计数k-mer | Count k-mers using Jellyfish"""
        self.logger.info(f"🧮 开始k-mer计数 | Starting k-mer counting for sample: {sample_name}")
        
        # 输出文件 | Output file
        jf_file = self.config.temp_dir / f"{sample_name}.jf"
        
        # 构建命令 | Build command
        cmd_parts = [
            self.config.jellyfish_path,
            "count",
            f"-m {self.config.kmer_size}",
            f"-s {self.config.hash_size}",
            f"-t {self.config.threads}",
        ]
        
        if self.config.canonical:
            cmd_parts.append("-C")
        
        # 添加输入文件 | Add input files
        input_files = [f for f in fastq_files if f is not None]
        cmd_parts.extend(input_files)
        cmd_parts.extend(["-o", str(jf_file)])
        
        cmd = " ".join(cmd_parts)
        
        # 执行命令 | Execute command
        self.cmd_runner.run(cmd, f"🧮 K-mer计数 | K-mer counting for {sample_name}")
        
        if not jf_file.exists():
            raise FileNotFoundError(f"❌ Jellyfish输出文件未生成 | Jellyfish output file not generated: {jf_file}")
        
        self.logger.info(f"✅ K-mer计数完成 | K-mer counting completed: {jf_file}")
        return jf_file
    
    def query_kmers(self, sample_name: str, jf_file: Path, kmer_lib: Path) -> Path:
        """🔍 查询k-mer丰度 | Query k-mer abundance"""
        self.logger.info(f"🔍 查询k-mer丰度 | Querying k-mer abundance for sample: {sample_name}")
        
        # 输出文件 | Output file
        count_file = self.config.temp_dir / f"{sample_name}.kmers.count"
        
        # 构建命令 | Build command
        cmd = f"{self.config.jellyfish_path} query -s {kmer_lib} {jf_file} -o {count_file}"
        
        # 执行命令 | Execute command
        self.cmd_runner.run(cmd, f"🔍 查询k-mer丰度 | Querying k-mer abundance for {sample_name}")
        
        if not count_file.exists():
            raise FileNotFoundError(f"❌ 查询输出文件未生成 | Query output file not generated: {count_file}")
        
        self.logger.info(f"✅ K-mer丰度查询完成 | K-mer abundance query completed: {count_file}")
        return count_file
