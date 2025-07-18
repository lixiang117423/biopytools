"""
K-mer分析主程序模块 | K-mer Analysis Main Module
"""

import argparse
import sys
import time
from .config import KmerConfig
from .utils import KmerLogger, CommandRunner, FileValidator, SequenceUtils, check_dependencies
from .data_processing import FastqFileDetector, FOFFileGenerator, GeneKmerExtractor
from .database import KmtricksRunner, RocksDBManager
from .analysis import KmerMatrixProcessor, HaplotypeAnalyzer
from .results import SummaryGenerator, FileManager

class KmerAnalyzer:
    """K-mer分析主类 | Main K-mer Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 检查依赖 | Check dependencies
        if not check_dependencies():
            sys.exit(1)
        
        # 初始化配置 | Initialize configuration
        self.config = KmerConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = KmerLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化工具类 | Initialize utility classes
        self.cmd_runner = CommandRunner(self.logger)
        self.file_validator = FileValidator(self.logger)
        self.seq_utils = SequenceUtils()
        
        # 初始化各个处理器 | Initialize processors
        self.fastq_detector = FastqFileDetector(self.logger)
        self.fof_generator = FOFFileGenerator(self.config, self.logger, self.fastq_detector)
        self.gene_extractor = GeneKmerExtractor(self.config, self.logger, self.seq_utils)
        
        self.kmtricks_runner = KmtricksRunner(self.config, self.logger, self.cmd_runner, self.file_validator)
        self.rocksdb_manager = RocksDBManager(self.config, self.logger, self.cmd_runner, self.file_validator)
        
        self.matrix_processor = KmerMatrixProcessor(self.config, self.logger)
        self.haplotype_analyzer = HaplotypeAnalyzer(self.config, self.logger)
        
        self.summary_generator = SummaryGenerator(self.config, self.logger)
        self.file_manager = FileManager(self.config, self.logger)
    
    def step1_create_fof_file(self):
        """步骤1: 创建FOF文件 | Step 1: Create FOF file"""
        return self.fof_generator.create_fof_file()
    
    def step2_build_kmtricks_database(self):
        """步骤2: 构建kmtricks数据库 | Step 2: Build kmtricks database"""
        if not self.kmtricks_runner.run_kmtricks_pipeline():
            return False
        return self.kmtricks_runner.run_kmtricks_aggregate()
    
    def step3_create_rocksdb(self):
        """步骤3: 创建RocksDB数据库 | Step 3: Create RocksDB database"""
        return self.rocksdb_manager.create_rocksdb_database()
    
    def step4_extract_gene_kmers(self):
        """步骤4: 提取基因k-mer | Step 4: Extract gene k-mers"""
        return self.gene_extractor.extract_gene_kmers()
    
    def step5_query_kmers(self):
        """步骤5: 查询k-mer | Step 5: Query k-mers"""
        return self.rocksdb_manager.query_kmers_from_database()
    
    def step6_process_results(self, gene_kmers):
        """步骤6: 处理查询结果 | Step 6: Process query results"""
        return self.matrix_processor.process_query_results(gene_kmers)
    
    def step7_haplotype_analysis(self, kmer_matrix):
        """步骤7: 单倍型分析 | Step 7: Haplotype analysis"""
        return self.haplotype_analyzer.run_haplotype_analysis(kmer_matrix)
    
    def run_single_step(self, step_num: int):
        """运行单个步骤 | Run single step"""
        step_functions = {
            1: (self.step1_create_fof_file, "创建FOF文件 | Create FOF file"),
            2: (self.step2_build_kmtricks_database, "构建kmtricks数据库 | Build kmtricks database"),
            3: (self.step3_create_rocksdb, "创建RocksDB数据库 | Create RocksDB database"),
            4: (self.step4_extract_gene_kmers, "提取基因k-mer | Extract gene k-mers"),
            5: (self.step5_query_kmers, "查询k-mer | Query k-mers"),
        }
        
        if step_num not in step_functions:
            self.logger.error(f"无效的步骤编号 | Invalid step number: {step_num}")
            return False
        
        step_func, step_name = step_functions[step_num]
        self.logger.info(f"执行步骤 {step_num} | Executing step {step_num}: {step_name}")
        
        success = step_func()
        if success:
            self.logger.info(f"步骤 {step_num} 完成 | Step {step_num} completed: {step_name}")
        else:
            self.logger.error(f"步骤 {step_num} 失败 | Step {step_num} failed: {step_name}")
        
        return success
    
    def run_analysis(self):
        """运行完整的分析流程 | Run complete analysis pipeline"""
        self.logger.info("开始高性能k-mer数据库分析流水线 | Starting high-performance k-mer database analysis pipeline...")
        
        start_time = time.time()
        gene_count = 0
        kmer_count = 0
        
        try:
            # 步骤1: 创建FOF文件 | Step 1: Create FOF file
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤1: 创建FOF文件 | Step 1: Create FOF file")
            self.logger.info("="*60)
            self.step1_create_fof_file()
            
            # 步骤2: 构建kmtricks数据库 | Step 2: Build kmtricks database
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤2: 构建kmtricks数据库 | Step 2: Build kmtricks database")
            self.logger.info("="*60)
            if not self.step2_build_kmtricks_database():
                raise Exception("kmtricks数据库构建失败 | kmtricks database build failed")
            
            # 步骤3: 创建RocksDB数据库 | Step 3: Create RocksDB database
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤3: 创建RocksDB数据库 | Step 3: Create RocksDB database")
            self.logger.info("="*60)
            if not self.step3_create_rocksdb():
                raise Exception("RocksDB数据库创建失败 | RocksDB database creation failed")
            
            # 步骤4: 提取基因k-mer | Step 4: Extract gene k-mers
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤4: 提取基因k-mer | Step 4: Extract gene k-mers")
            self.logger.info("="*60)
            gene_kmers = self.step4_extract_gene_kmers()
            gene_count = len(gene_kmers)
            kmer_count = sum(len(kmers) for kmers in gene_kmers.values())
            
            # 步骤5: 查询k-mer | Step 5: Query k-mers
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤5: 查询k-mer | Step 5: Query k-mers")
            self.logger.info("="*60)
            if not self.step5_query_kmers():
                raise Exception("k-mer查询失败 | k-mer query failed")
            
            # 步骤6: 处理查询结果 | Step 6: Process query results
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤6: 处理查询结果 | Step 6: Process query results")
            self.logger.info("="*60)
            kmer_matrix = self.step6_process_results(gene_kmers)
            
            # 步骤7: 单倍型分析 (可选) | Step 7: Haplotype analysis (optional)
            if self.config.run_haplotype:
                self.logger.info("\n" + "="*60)
                self.logger.info("步骤7: 单倍型分析 | Step 7: Haplotype analysis")
                self.logger.info("="*60)
                self.step7_haplotype_analysis(kmer_matrix)
            
            # 步骤8: 生成报告 | Step 8: Generate report
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤8: 生成报告 | Step 8: Generate report")
            self.logger.info("="*60)
            self.summary_generator.generate_report(start_time, gene_count, kmer_count)
            
            # 检查输出完整性 | Check output completeness
            if not self.file_manager.check_output_completeness():
                self.logger.warning("部分输出文件缺失 | Some output files are missing")
            
            # 清理临时文件 | Clean up temporary files
            self.file_manager.cleanup_temp_files()
            
            # 计算运行时间 | Calculate runtime
            end_time = time.time()
            runtime = int(end_time - start_time)
            runtime_str = f"{runtime // 3600}小时{(runtime % 3600) // 60}分{runtime % 60}秒"
            
            self.logger.info("\n" + "="*60)
            self.logger.info("✓ 高性能k-mer分析流水线完成！| High-performance k-mer analysis pipeline completed!")
            self.logger.info(f"总用时 | Total runtime: {runtime_str}")
            self.logger.info(f"主要结果文件 | Main result file: {self.config.kmer_matrix_final}")
            if not self.config.skip_build:
                self.logger.info("数据库已构建，后续查询将非常快速 | Database built, future queries will be very fast")
                self.logger.info("使用 --skip-build 参数可快速查询新基因 | Use --skip-build flag for fast queries of new genes")
            self.logger.info("="*60)
            
            return True
            
        except Exception as e:
            self.logger.error(f"流水线执行失败 | Pipeline execution failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="""
高性能k-mer数据库查询流水线 | High-Performance K-mer Database Query Pipeline
基于kmtricks + RocksDB的大规模k-mer分析系统 | Large-scale k-mer analysis system based on kmtricks + RocksDB

专为大规模数据集设计 (支持数千个样本) | Designed for large-scale datasets (supporting thousands of samples):
1. 一次构建全基因组k-mer数据库 | Build genome-wide k-mer database once
2. 支持快速查询任意基因的k-mer模式 | Support fast queries of k-mer patterns for any genes
3. 查询速度: 秒级到分钟级 (vs 小时级的暴力搜索) | Query speed: seconds to minutes (vs hours of brute force)

适用场景 | Use Cases:
- 大规模群体基因组学研究 | Large-scale population genomics studies
- 需要重复查询不同基因的场景 | Scenarios requiring repeated queries of different genes
- 对查询速度有高要求的项目 | Projects with high query speed requirements

性能对比 (5000个样本) | Performance Comparison (5000 samples):
- 传统方法 | Traditional method: 每次查询数周 | weeks per query
- 本方法 | This method: 构建一次(1-2天) + 查询(1-5分钟) | build once (1-2 days) + query (1-5 minutes)
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required Arguments')
    required.add_argument(
        '-g', '--gene-fasta',
        required=True,
        help='目标基因FASTA文件路径 | Target gene FASTA file path'
    )
    required.add_argument(
        '-f', '--fastq-dir',
        required=True,
        help='FASTQ文件目录路径 | FASTQ file directory path'
    )
    required.add_argument(
        '-o', '--output-dir',
        required=True,
        help='输出目录路径 | Output directory path'
    )
    
    # k-mer参数 | K-mer parameters
    kmer_group = parser.add_argument_group('k-mer分析参数 | K-mer Analysis Parameters')
    kmer_group.add_argument(
        '-k', '--kmer-size',
        type=int,
        default=51,
        help='k-mer大小 (默认: 51) | k-mer size (default: 51)'
    )
    kmer_group.add_argument(
        '-t', '--threads',
        type=int,
        default=32,
        help='线程数 (默认: 32，建议使用较多线程) | Thread count (default: 32, recommend using more threads)'
    )
    kmer_group.add_argument(
        '-m', '--hard-min',
        type=int,
        default=2,
        help='最小k-mer频次阈值 (默认: 2) | Minimum k-mer frequency threshold (default: 2)'
    )
    
    # 流程控制参数 | Process control parameters
    control_group = parser.add_argument_group('流程控制参数 | Process Control Parameters')
    control_group.add_argument(
        '-p', '--project-name',
        help='项目名称 (默认: 从输出目录名获取) | Project name (default: derived from output directory name)'
    )
    control_group.add_argument(
        '--skip-build',
        action='store_true',
        help='跳过数据库构建步骤 (用于已有数据库的查询) | Skip database build step (for querying existing database)'
    )
    control_group.add_argument(
        '--run-haplotype',
        action='store_true',
        help='运行单倍型聚类分析 | Run haplotype clustering analysis'
    )
    control_group.add_argument(
        '-n', '--n-components',
        type=int,
        default=5,
        help='BGMM最大聚类数 (默认: 5) | Maximum number of BGMM clusters (default: 5)'
    )
    
    args = parser.parse_args()
    
    try:
        analyzer = KmerAnalyzer(**vars(args))
        success = analyzer.run_analysis()
        sys.exit(0 if success else 1)
        
    except KeyboardInterrupt:
        print("\n用户中断操作 | User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行异常 | Program execution exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
