"""
K-mer分析结果处理模块 | K-mer Analysis Results Processing Module
"""

import time
from pathlib import Path

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_report(self, start_time: float, gene_count: int = 0, kmer_count: int = 0):
        """生成分析报告 | Generate analysis report"""
        self.logger.info("生成分析报告 | Generating analysis report...")
        
        # 计算运行时间
        end_time = time.time()
        runtime = int(end_time - start_time)
        runtime_str = f"{runtime // 3600}小时{(runtime % 3600) // 60}分{runtime % 60}秒"
        
        report_content = f"""
=================================================
高性能k-mer数据库分析报告 | High-Performance K-mer Database Analysis Report
=================================================

分析时间 | Analysis Time: {time.strftime('%Y-%m-%d %H:%M:%S')}
项目名称 | Project Name: {self.config.project_name}
总运行时间 | Total Runtime: {runtime_str}

输入参数 | Input Parameters:
- 基因FASTA文件 | Gene FASTA file: {self.config.gene_fasta}
- FASTQ文件目录 | FASTQ directory: {self.config.fastq_dir}
- 输出目录 | Output directory: {self.config.output_dir}
- k-mer大小 | k-mer size: {self.config.kmer_size}
- 线程数 | Thread count: {self.config.threads}
- 最小频次 | Minimum frequency: {self.config.hard_min}

分析统计 | Analysis Statistics:
- 处理基因数 | Genes processed: {gene_count}
- 提取k-mer数 | K-mers extracted: {kmer_count}

数据库文件 | Database Files:
- kmtricks索引 | kmtricks index: {self.config.kmtricks_run_dir}
- k-mer矩阵 | k-mer matrix: {self.config.kmer_matrix_file}
- RocksDB数据库 | RocksDB database: {self.config.rocksdb_dir}

查询结果 | Query Results:
- 基因k-mer列表 | Gene k-mer list: {self.config.gene_kmer_file}
- 查询结果 | Query results: {self.config.query_result_file}
- 最终k-mer矩阵 | Final k-mer matrix: {self.config.kmer_matrix_final}

性能优势 | Performance Advantages:
- 一次构建，多次查询 | Build once, query multiple times
- 查询速度: 秒级到分钟级 | Query speed: seconds to minutes
- 支持大规模数据集 (数千样本) | Support large datasets (thousands of samples)

后续使用 | Future Usage:
- 数据库已构建，可使用 --skip-build 参数快速查询新基因
- Database built, use --skip-build flag for fast queries of new genes
- 查询新基因只需1-5分钟 | Querying new genes takes only 1-5 minutes

=================================================
"""
        
        report_file = self.config.output_path / "analysis_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.logger.info(f"✓ 分析报告已生成 | Analysis report generated: {report_file}")
        return report_file

class FileManager:
    """文件管理器 | File Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def cleanup_temp_files(self):
        """清理临时文件 | Clean up temporary files"""
        self.logger.info("清理临时文件 | Cleaning up temporary files...")
        
        temp_files = [
            self.config.output_path / "import_kmer_rocksdb.py",
            self.config.output_path / "search_kmer_rocksdb.py",
            self.config.header_file
        ]
        
        for temp_file in temp_files:
            if temp_file.exists():
                try:
                    temp_file.unlink()
                    self.logger.debug(f"已删除临时文件 | Deleted temp file: {temp_file}")
                except Exception as e:
                    self.logger.warning(f"删除临时文件失败 | Failed to delete temp file {temp_file}: {e}")
    
    def check_output_completeness(self) -> bool:
        """检查输出文件完整性 | Check output file completeness"""
        self.logger.info("检查输出文件完整性 | Checking output file completeness...")
        
        required_files = [
            (self.config.fof_file, "FOF文件 | FOF file"),
            (self.config.gene_kmer_file, "基因k-mer文件 | Gene k-mer file"),
            (self.config.kmer_matrix_final, "最终k-mer矩阵 | Final k-mer matrix"),
        ]
        
        if not self.config.skip_build:
            required_files.extend([
                (self.config.kmtricks_run_dir, "kmtricks索引目录 | kmtricks index directory"),
                (self.config.kmer_matrix_file, "k-mer矩阵文件 | k-mer matrix file"),
                (self.config.rocksdb_dir, "RocksDB数据库目录 | RocksDB database directory"),
            ])
        
        missing_files = []
        for file_path, description in required_files:
            if not file_path.exists():
                missing_files.append(f"{description}: {file_path}")
        
        if missing_files:
            self.logger.error("缺少必要的输出文件 | Missing required output files:")
            for missing in missing_files:
                self.logger.error(f"  - {missing}")
            return False
        
        self.logger.info("✓ 所有必要文件检查完成 | All required files check completed")
        return True
