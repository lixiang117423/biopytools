"""
🔍 K-mer查询模块 | K-mer Query Module
"""

from pathlib import Path

class KmerQueryAnalyzer:
    """K-mer查询分析器 | K-mer Query Analyzer"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def count_kmers(self) -> Path:
        """
        统计查询文件中的k-mer丰度 | Count k-mer abundance in query file
        
        Returns:
            计数文件路径 | Count file path
        """
        self.logger.info("📊 步骤3: 统计查询文件中的k-mer丰度 | Step 3: Counting k-mer abundance in query file")
        
        # 输出文件 | Output files
        kmc_db = self.config.tmp_dir / "query_kmc"
        counts_file = self.config.output_path / "query_counts.txt"
        
        # 使用绝对路径 | Use absolute paths
        kmc_db_abs = kmc_db.resolve()
        tmp_dir_abs = self.config.tmp_dir.resolve()
        counts_file_abs = counts_file.resolve()
        query_abs = Path(self.config.query).resolve()
        
        # 步骤1: 运行KMC统计 | Step 1: Run KMC counting
        cmd = (
            f"{self.config.kmc_path} "
            f"-k{self.config.kmer_length} "
            f"-ci{self.config.min_count} "
            f"-t{self.config.threads} "
            f"-m{self.config.memory_gb} "
            f"{query_abs} "
            f"{kmc_db_abs} "
            f"{tmp_dir_abs}"
        )
        
        if not self.cmd_runner.run(cmd, "使用KMC统计查询文件k-mer"):
            raise RuntimeError("❌ KMC统计失败 | KMC counting failed")
        
        # 步骤2: 导出计数结果 | Step 2: Export count results
        cmd = (
            f"{self.config.kmc_tools_path} transform "
            f"{kmc_db_abs} "
            f"dump "
            f"{counts_file_abs}"
        )
        
        if not self.cmd_runner.run(cmd, "导出k-mer计数"):
            raise RuntimeError("❌ k-mer计数导出失败 | k-mer count export failed")
        
        # 统计结果 | Count statistics
        with open(counts_file, 'r') as f:
            line_count = sum(1 for _ in f)
        
        self.logger.info(f"✅ 查询文件中发现 {line_count:,} 个不同的k-mer | Found {line_count:,} unique k-mers in query")
        self.logger.info(f"📝 计数文件已保存 | Count file saved: {counts_file}")
        
        return counts_file
