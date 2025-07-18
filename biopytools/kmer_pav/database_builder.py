"""
K-mer数据库构建模块 | K-mer Database Builder Module
"""

import os
from pathlib import Path
from typing import List, Dict
from .utils import KMCRunner, FileDetector, find_sequence_files, cleanup_files

class DatabaseBuilder:
    """K-mer数据库构建器 | K-mer Database Builder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.kmc_runner = KMCRunner(config, logger)
        self.file_detector = FileDetector(logger)
    
    def get_database_files(self) -> List[Dict[str, str]]:
        """获取数据库文件列表 | Get database file list"""
        files = find_sequence_files(self.config.database_input, self.config.database_pattern)
        
        if not files:
            raise ValueError(f"未找到数据库文件 | No database files found in: {self.config.database_input}")
        
        db_files = []
        for file_path in files:
            file_format = self.file_detector.detect_format(file_path)
            
            db_files.append({
                'name': Path(file_path).stem,
                'file': file_path,
                'format': file_format
            })
        
        self.logger.info(f"找到 {len(db_files)} 个数据库文件 | Found {len(db_files)} database files")
        return db_files
    
    def build_individual_databases(self, db_files: List[Dict[str, str]]) -> List[str]:
        """为每个数据库文件构建k-mer数据库 | Build k-mer databases for each database file"""
        self.logger.info("开始构建各个文件的k-mer数据库 | Starting to build individual k-mer databases")
        
        db_names = []
        
        for i, db_file in enumerate(db_files):
            db_name = f"db_{i:04d}"
            file_path = db_file['file']
            file_format = db_file['format']
            
            self.logger.info(f"处理数据库文件 {i+1}/{len(db_files)}: {db_file['name']}")
            
            # 确定格式标志 | Determine format flag
            format_flag = "-fm" if file_format == "fasta" else "-fq"
            
            # 使用绝对路径构建输出文件名 | Use absolute paths for output file names
            output_db_path = str(self.config.output_path / db_name)
            
            # 构建KMC命令 | Build KMC command
            cmd = [
                self.config.kmc_path,
                f"-k{self.config.kmer_size}",
                f"-ci{self.config.min_count}",
                f"-cx{self.config.max_count}",
                f"-t{self.config.threads}",
                f"-m{self.config.kmc_memory_gb}",
                format_flag,  # 添加格式标志 | Add format flag
                file_path,
                output_db_path,  # 使用绝对路径 | Use absolute path
                str(self.config.tmp_path)  # 临时目录绝对路径 | Temp directory absolute path
            ]
            
            if self.kmc_runner.run_command(cmd, f"构建数据库 | Build database: {db_file['name']}"):
                # 验证生成的文件是否存在 | Verify generated files exist
                db_pre_file = self.config.output_path / f"{db_name}.kmc_pre"
                db_suf_file = self.config.output_path / f"{db_name}.kmc_suf"
                
                if db_pre_file.exists() and db_suf_file.exists():
                    db_names.append(db_name)
                    self.logger.info(f"数据库文件验证成功 | Database files verified: {db_name}")
                    self.logger.info(f"  - {db_pre_file} (大小: {db_pre_file.stat().st_size} bytes)")
                    self.logger.info(f"  - {db_suf_file} (大小: {db_suf_file.stat().st_size} bytes)")
                else:
                    self.logger.error(f"数据库文件未找到 | Database files not found: {db_name}")
                    self.logger.error(f"查找位置 | Searching at: {self.config.output_path}")
                    self.logger.error(f"期望文件 | Expected files: {db_name}.kmc_pre, {db_name}.kmc_suf")
            else:
                self.logger.error(f"构建数据库失败 | Failed to build database: {db_file['name']}")
                # 继续处理其他文件 | Continue with other files
        
        self.logger.info(f"成功构建 {len(db_names)} 个数据库 | Successfully built {len(db_names)} databases")
        return db_names
    
    def merge_databases(self, db_names: List[str]) -> str:
        """合并所有数据库为统一的k-mer数据库 | Merge all databases into unified k-mer database"""
        if not db_names:
            raise ValueError("没有可合并的数据库 | No databases to merge")
        
        # 调试：列出当前工作目录中的文件 | Debug: list files in current working directory
        self.logger.info(f"当前工作目录内容 | Current working directory contents:")
        try:
            import os
            for file in os.listdir(self.config.output_path):
                if file.startswith('db_') and ('.kmc_pre' in file or '.kmc_suf' in file):
                    self.logger.info(f"  找到KMC文件 | Found KMC file: {file}")
        except Exception as e:
            self.logger.warning(f"无法列出目录内容 | Cannot list directory contents: {e}")
        
        if len(db_names) == 1:
            # 只有一个数据库，直接重命名 | Only one database, rename directly
            final_db_name = "merged_kmer_database"
            self.logger.info(f"只有一个数据库，重命名为最终数据库 | Only one database, renaming to final database")
            
            # 复制数据库文件 | Copy database files
            import shutil
            for suffix in ['.kmc_pre', '.kmc_suf']:
                src = self.config.output_path / f"{db_names[0]}{suffix}"
                dst = self.config.output_path / f"{final_db_name}{suffix}"
                if src.exists():
                    shutil.copy2(src, dst)
                    self.logger.info(f"复制文件 | Copied file: {src} -> {dst}")
                else:
                    self.logger.error(f"源文件不存在 | Source file does not exist: {src}")
            
            return final_db_name
        
        self.logger.info(f"开始合并 {len(db_names)} 个数据库 | Starting to merge {len(db_names)} databases")
        
        # 验证所有数据库文件都存在 | Verify all database files exist
        for db_name in db_names:
            for suffix in ['.kmc_pre', '.kmc_suf']:
                file_path = self.config.output_path / f"{db_name}{suffix}"
                if not file_path.exists():
                    raise RuntimeError(f"数据库文件不存在 | Database file not found: {file_path}")
        
        # 逐步合并数据库 | Merge databases step by step
        current_db = db_names[0]
        
        for i, next_db in enumerate(db_names[1:], 1):
            merged_db = f"merged_temp_{i:04d}"
            
            # 验证输入文件存在 | Verify input files exist
            for db in [current_db, next_db]:
                for suffix in ['.kmc_pre', '.kmc_suf']:
                    file_path = self.config.output_path / f"{db}{suffix}"
                    if not file_path.exists():
                        raise RuntimeError(f"合并输入文件不存在 | Merge input file not found: {file_path}")
            
            # 使用绝对路径构建数据库名称 | Use absolute paths for database names
            current_db_path = str(self.config.output_path / current_db)
            next_db_path = str(self.config.output_path / next_db)
            merged_db_path = str(self.config.output_path / merged_db)
            
            cmd = [
                self.config.kmc_tools_path,
                "simple",
                current_db_path,  # 使用绝对路径 | Use absolute path
                next_db_path,     # 使用绝对路径 | Use absolute path
                "union",
                merged_db_path    # 使用绝对路径 | Use absolute path
            ]
            
            if self.kmc_runner.run_command(cmd, f"合并数据库 {i}/{len(db_names)-1} | Merge database {i}/{len(db_names)-1}"):
                # 清理上一个临时文件 | Clean up previous temporary files
                if current_db.startswith("merged_temp_"):
                    cleanup_files([f"{current_db}.kmc*"], self.config.output_path, self.logger)
                
                current_db = merged_db
            else:
                raise RuntimeError(f"合并数据库失败 | Failed to merge database at step {i}")
        
        # 重命名最终数据库 | Rename final database
        final_db_name = "merged_kmer_database"
        import shutil
        for suffix in ['.kmc_pre', '.kmc_suf']:
            src = self.config.output_path / f"{current_db}{suffix}"
            dst = self.config.output_path / f"{final_db_name}{suffix}"
            if src.exists():
                shutil.move(src, dst)
                self.logger.info(f"重命名文件 | Renamed file: {src} -> {dst}")
            else:
                self.logger.error(f"最终数据库文件不存在 | Final database file not found: {src}")
        
        self.logger.info(f"数据库合并完成 | Database merging completed: {final_db_name}")
        
        # 清理各个单独数据库 | Clean up individual databases
        if not self.config.keep_intermediate:
            for db_name in db_names:
                cleanup_files([f"{db_name}.kmc*"], self.config.output_path, self.logger)
        
        return final_db_name
    
    def build_unified_database(self) -> str:
        """构建统一的k-mer数据库 | Build unified k-mer database"""
        self.logger.info("=" * 60)
        self.logger.info("阶段1: 构建k-mer数据库 | Phase 1: Building k-mer database")
        self.logger.info("=" * 60)
        
        # 1. 获取数据库文件 | Get database files
        db_files = self.get_database_files()
        
        # 2. 构建各个数据库 | Build individual databases
        db_names = self.build_individual_databases(db_files)
        
        # 3. 合并数据库 | Merge databases
        final_db_name = self.merge_databases(db_names)
        
        # 4. 导出k-mer列表（可选，用于统计） | Export k-mer list (optional, for statistics)
        kmer_list_file = f"{final_db_name}_kmers.txt"
        final_db_path = str(self.config.output_path / final_db_name)
        kmer_list_path = str(self.config.output_path / kmer_list_file)
        
        cmd = [
            self.config.kmc_tools_path,
            "transform",
            final_db_path,    # 使用绝对路径 | Use absolute path
            "dump",
            "-s",
            kmer_list_path    # 使用绝对路径 | Use absolute path
        ]
        
        if self.kmc_runner.run_command(cmd, "导出k-mer列表 | Export k-mer list"):
            # 统计k-mer数量 | Count k-mers
            kmer_count = self._count_kmers_in_file(kmer_list_file)
            self.logger.info(f"数据库包含 {kmer_count:,} 个唯一k-mer | Database contains {kmer_count:,} unique k-mers")
            
            # 删除临时文件 | Remove temporary file
            if not self.config.keep_intermediate:
                kmer_list_file_path = self.config.output_path / kmer_list_file
                if kmer_list_file_path.exists():
                    kmer_list_file_path.unlink()
        
        self.logger.info(f"k-mer数据库构建完成 | K-mer database construction completed: {final_db_name}")
        return final_db_name
    
    def _count_kmers_in_file(self, kmer_file: str) -> int:
        """统计文件中的k-mer数量 | Count k-mers in file"""
        try:
            with open(self.config.output_path / kmer_file, 'r') as f:
                count = sum(1 for line in f if line.strip())
            return count
        except Exception as e:
            self.logger.warning(f"无法统计k-mer数量 | Cannot count k-mers: {e}")
            return 0
