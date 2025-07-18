"""
数据库操作模块 | Database Operations Module
"""

import os
import gzip
from pathlib import Path
from .utils import CommandRunner, FileValidator

class KmtricksRunner:
    """Kmtricks运行器 | Kmtricks Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner, file_validator: FileValidator):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = file_validator
    
    def run_kmtricks_pipeline(self) -> bool:
        """运行kmtricks k-mer构建流水线 | Run kmtricks k-mer construction pipeline"""
        if self.config.skip_build and self.config.kmtricks_run_dir.exists():
            self.logger.info("跳过kmtricks构建步骤 (已存在) | Skipping kmtricks build step (already exists)")
            return True
            
        self.logger.info("开始运行kmtricks k-mer构建流水线 | Starting kmtricks k-mer construction pipeline...")
        self.logger.info("警告: 这一步可能需要数小时到一天时间，但只需要执行一次 | Warning: This step may take hours to a day, but only needs to be run once")
        
        # 检查kmtricks是否可用
        if not self.file_validator.check_tool_available('kmtricks'):
            self.logger.error("请安装kmtricks: conda install -c bioconda kmtricks | Please install kmtricks: conda install -c bioconda kmtricks")
            return False
        
        # 检查FOF文件是否存在
        if not self.config.fof_file.exists():
            self.logger.error(f"FOF文件不存在 | FOF file does not exist: {self.config.fof_file}")
            return False
        
        # 使用绝对路径
        fof_file_abs = self.config.fof_file.resolve()
        kmtricks_run_dir_abs = self.config.kmtricks_run_dir.resolve()
        
        self.logger.info(f"FOF文件路径 | FOF file path: {fof_file_abs}")
        self.logger.info(f"kmtricks运行目录 | kmtricks run directory: {kmtricks_run_dir_abs}")
        
        # 构建k-mer索引
        cmd = [
            'kmtricks', 'pipeline',
            '-t', str(self.config.threads),
            '--file', str(fof_file_abs),
            '--run-dir', str(kmtricks_run_dir_abs),
            '--mode', 'kmer:pa:bin',
            '--hard-min', str(self.config.hard_min),
            '--kmer-size', str(self.config.kmer_size),
            '--cpr'
        ]
        
        try:
            self.logger.info("这可能需要很长时间，请耐心等待 | This may take a long time, please be patient...")
            
            success = self.cmd_runner.run_with_realtime_output(
                cmd, "kmtricks k-mer构建 | kmtricks k-mer construction"
            )
            
            if success:
                self.logger.info("✓ kmtricks k-mer构建完成 | kmtricks k-mer construction completed")
            else:
                self.logger.error("kmtricks pipeline执行失败 | kmtricks pipeline execution failed")
            
            return success
            
        except Exception as e:
            self.logger.error(f"kmtricks执行异常 | kmtricks execution exception: {e}")
            return False
    
    def run_kmtricks_aggregate(self) -> bool:
        """运行kmtricks聚合生成矩阵 | Run kmtricks aggregation to generate matrix"""
        if self.config.skip_build and self.config.kmer_matrix_file.exists():
            self.logger.info("跳过kmtricks聚合步骤 (已存在) | Skipping kmtricks aggregation step (already exists)")
            return True
            
        self.logger.info("运行kmtricks聚合生成k-mer矩阵 | Running kmtricks aggregation to generate k-mer matrix...")
        
        # 使用绝对路径
        kmtricks_run_dir_abs = self.config.kmtricks_run_dir.resolve()
        kmer_matrix_file_abs = self.config.kmer_matrix_file.resolve()
        
        self.logger.info(f"kmtricks运行目录 | kmtricks run directory: {kmtricks_run_dir_abs}")
        self.logger.info(f"输出矩阵文件 | Output matrix file: {kmer_matrix_file_abs}")
        
        cmd = (f"kmtricks aggregate -t {self.config.threads} "
               f"--run-dir {kmtricks_run_dir_abs} "
               f"--pa-matrix kmer --format text --cpr-in | "
               f"gzip > {kmer_matrix_file_abs}")
        
        try:
            success = self.cmd_runner.run(
                cmd, "kmtricks聚合 | kmtricks aggregation", shell=True
            )
            
            if success:
                # 检查输出文件
                if self.config.kmer_matrix_file.exists():
                    file_size = self.config.kmer_matrix_file.stat().st_size / (1024**3)  # GB
                    self.logger.info(f"✓ k-mer矩阵已生成 | k-mer matrix generated: {self.config.kmer_matrix_file}")
                    self.logger.info(f"文件大小 | File size: {file_size:.2f} GB")
                else:
                    self.logger.error("聚合完成但未找到输出文件 | Aggregation completed but output file not found")
                    return False
            
            return success
            
        except Exception as e:
            self.logger.error(f"kmtricks聚合失败 | kmtricks aggregation failed: {e}")
            return False

class RocksDBManager:
    """RocksDB管理器 | RocksDB Manager"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner, file_validator: FileValidator):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = file_validator
    
    def create_rocksdb_database(self) -> bool:
        """创建RocksDB数据库 | Create RocksDB database"""
        if self.config.skip_build and self.config.rocksdb_dir.exists():
            self.logger.info("跳过RocksDB创建步骤 (已存在) | Skipping RocksDB creation step (already exists)")
            return True
            
        self.logger.info("创建RocksDB数据库 | Creating RocksDB database...")
        
        # 生成header文件
        try:
            with gzip.open(self.config.kmer_matrix_file, 'rt') as f:
                header = f.readline().strip()
            with open(self.config.header_file, 'w') as f:
                f.write(header + '\n')
        except Exception as e:
            self.logger.error(f"生成header文件失败 | Failed to generate header file: {e}")
            return False
        
        # 创建导入脚本
        import_script = self.create_rocksdb_import_script()
        
        # 使用绝对路径
        header_file_abs = self.config.header_file.resolve()
        kmer_matrix_file_abs = self.config.kmer_matrix_file.resolve()
        rocksdb_dir_abs = self.config.rocksdb_dir.resolve()
        
        # 运行导入
        cmd = [
            'python3', str(import_script),
            '--input_delimiter', ' ',
            '--header_file', str(header_file_abs),
            '--header_db_key', self.config.project_name,
            str(kmer_matrix_file_abs),
            str(rocksdb_dir_abs)
        ]
        
        try:
            self.logger.info("开始导入k-mer矩阵到RocksDB | Starting import of k-mer matrix to RocksDB...")
            self.logger.info("这可能需要较长时间 | This may take a long time...")
            
            success = self.cmd_runner.run(cmd, "RocksDB导入 | RocksDB import")
            
            if success:
                self.logger.info("✓ RocksDB数据库创建完成 | RocksDB database creation completed")
            
            return success
            
        except Exception as e:
            self.logger.error(f"RocksDB导入失败 | RocksDB import failed: {e}")
            return False
    
    def create_rocksdb_import_script(self) -> Path:
        """创建RocksDB导入脚本 | Create RocksDB import script"""
        script_path = self.config.output_path / "import_kmer_rocksdb.py"
        
        script_content = '''#!/usr/bin/env python3
import argparse
import gzip
import rocksdb
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_delimiter', default=' ')
    parser.add_argument('--header_file', required=True)
    parser.add_argument('--header_db_key', required=True)
    parser.add_argument('matrix_file', help='输入矩阵文件')
    parser.add_argument('db_path', help='RocksDB数据库路径')
    args = parser.parse_args()
    
    # 读取header
    with open(args.header_file, 'r') as f:
        header = f.readline().strip().split(args.input_delimiter)
    
    print(f"Header包含 {len(header)} 个样本")
    
    # 创建RocksDB
    opts = rocksdb.Options()
    opts.create_if_missing = True
    opts.max_open_files = -1
    db = rocksdb.DB(args.db_path, opts)
    
    # 存储header
    db.put(args.header_db_key.encode(), args.input_delimiter.join(header).encode())
    
    # 导入矩阵数据
    print("开始导入k-mer数据...")
    count = 0
    
    with gzip.open(args.matrix_file, 'rt') as f:
        next(f)  # 跳过header
        for line in f:
            parts = line.strip().split(args.input_delimiter)
            if len(parts) >= 2:
                kmer = parts[0]
                values = args.input_delimiter.join(parts[1:])
                db.put(kmer.encode(), values.encode())
                count += 1
                
                if count % 100000 == 0:
                    print(f"已导入 {count} 个k-mer")
    
    print(f"导入完成! 总计 {count} 个k-mer")
    db.close()

if __name__ == "__main__":
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        script_path.chmod(0o755)
        return script_path
    
    def query_kmers_from_database(self) -> bool:
        """从RocksDB数据库查询k-mer | Query k-mers from RocksDB database"""
        self.logger.info("从RocksDB数据库查询k-mer | Querying k-mers from RocksDB database...")
        
        # 创建查询脚本
        search_script = self.create_rocksdb_search_script()
        
        # 使用绝对路径
        rocksdb_dir_abs = self.config.rocksdb_dir.resolve()
        gene_kmer_file_abs = self.config.gene_kmer_file.resolve()
        query_result_file_abs = self.config.query_result_file.resolve()
        
        cmd = [
            'python3', str(search_script),
            '--header_db_key', self.config.project_name,
            '--output', str(query_result_file_abs),
            str(rocksdb_dir_abs),
            str(gene_kmer_file_abs)
        ]
        
        try:
            self.logger.info("开始查询k-mer | Starting k-mer query...")
            success = self.cmd_runner.run(cmd, "k-mer查询 | k-mer query")
            
            if success and self.config.query_result_file.exists():
                self.logger.info(f"✓ 查询完成 | Query completed: {self.config.query_result_file}")
                return True
            else:
                self.logger.error("查询完成但未找到结果文件 | Query completed but result file not found")
                return False
                
        except Exception as e:
            self.logger.error(f"k-mer查询失败 | k-mer query failed: {e}")
            return False
    
    def create_rocksdb_search_script(self) -> Path:
        """创建RocksDB查询脚本 | Create RocksDB search script"""
        script_path = self.config.output_path / "search_kmer_rocksdb.py"
        
        script_content = '''#!/usr/bin/env python3
import argparse
import rocksdb
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--header_db_key', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('db_path', help='RocksDB数据库路径')
    parser.add_argument('kmer_file', help='要查询的k-mer文件')
    args = parser.parse_args()
    
    # 打开数据库
    opts = rocksdb.Options()
    db = rocksdb.DB(args.db_path, opts, read_only=True)
    
    # 获取header
    header_data = db.get(args.header_db_key.encode())
    if header_data is None:
        print(f"错误: 未找到header key {args.header_db_key}")
        sys.exit(1)
    
    header = header_data.decode().split(' ')
    
    # 读取要查询的k-mer
    kmers_to_query = []
    with open(args.kmer_file, 'r') as f:
        next(f)  # 跳过header
        for line in f:
            kmer = line.strip()
            if kmer:
                kmers_to_query.append(kmer)
    
    print(f"需要查询 {len(kmers_to_query)} 个k-mer")
    
    # 查询并输出结果
    with open(args.output, 'w') as f:
        # 写入header
        f.write(' '.join(header) + '\\n')
        
        found_count = 0
        for kmer in kmers_to_query:
            result = db.get(kmer.encode())
            if result is not None:
                f.write(f"{kmer} {result.decode()}\\n")
                found_count += 1
            else:
                # 如果未找到，写入全0行
                zeros = ' '.join(['0'] * (len(header) - 1))
                f.write(f"{kmer} {zeros}\\n")
    
    print(f"查询完成! 找到 {found_count}/{len(kmers_to_query)} 个k-mer")
    db.close()

if __name__ == "__main__":
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        script_path.chmod(0o755)
        return script_path
