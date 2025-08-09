"""
主分析模块 | Main Analysis Module 🧬🔬
"""

import os
import time
from pathlib import Path
from typing import Dict, List, Optional
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from .config import KmerConfig, FileRole, AssignmentStrategy
from .file_manager import FileManager, FileInfo, FileFormat
from .kmc_interface import KMCInterface
from .position_tracker import PositionTracker, KmerInfo
from .output_manager import OutputManager

class KmerAnalyzer:
    """K-mer分析器主类 🤖"""
    
    def __init__(self, config: Optional[KmerConfig] = None):
        """初始化分析器 🚀"""
        self.config = config or KmerConfig()
        self.logger = self._setup_logging()
        
        # 初始化各组件
        self.file_manager = FileManager(self.config)
        self.kmc_interface = KMCInterface(self.config)
        self.position_tracker = PositionTracker(self.config)
        self.output_manager = OutputManager(self.config)
        
        self.logger.info("KmerAnalyzer initialized successfully ✅")
    
    def _setup_logging(self) -> logging.Logger:
        """设置日志 📜"""
        logger = logging.getLogger(__name__)
        
        if not logger.handlers:
            level = logging.DEBUG if self.config.verbose else logging.INFO
            logger.setLevel(level)
            
            # 控制台处理器 💻
            console_handler = logging.StreamHandler()
            console_handler.setLevel(level)
            
            # 文件处理器 📁
            log_file = os.path.join(self.config.output_dir, 'kmer_analysis.log')
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(level)
            
            # 格式化器 🎨
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            console_handler.setFormatter(formatter)
            file_handler.setFormatter(formatter)
            
            logger.addHandler(console_handler)
            logger.addHandler(file_handler)
        
        return logger
    
    def auto_analyze(self, input_paths: List[str], output_dir: Optional[str] = None, 
                    kmer_size: Optional[int] = None, 
                    assignment_strategy: AssignmentStrategy = AssignmentStrategy.INTELLIGENT) -> Dict:
        """
        自动分析模式：智能识别文件角色并执行分析 ✨
        
        Args:
            input_paths: 输入路径列表（支持文件、目录、通配符） 📁
            output_dir: 输出目录 📂
            kmer_size: k-mer长度 📏
            assignment_strategy: 角色分配策略 🤔
        
        Returns:
            分析结果字典 📊
        """
        start_time = time.time()
        
        # 更新配置
        if output_dir:
            self.config.output_dir = output_dir
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        if kmer_size:
            self.config.kmer_size = kmer_size
        
        self.logger.info("🚀 Starting auto analysis mode")
        self.logger.info(f"Input paths: {input_paths}")
        self.logger.info(f"Assignment strategy: {assignment_strategy.value}")
        
        try:
            # 步骤1: 扫描和识别文件 🔍
            all_files = self.file_manager.scan_files(input_paths)
            if not all_files:
                raise ValueError("No valid sequence files found")
            
            # 步骤2: 分配文件角色 🏷️
            assigned_files = self.file_manager.assign_roles(all_files, assignment_strategy)
            kmer_sources = assigned_files['kmer_sources']
            query_targets = assigned_files['query_targets']
            
            if not kmer_sources or not query_targets:
                raise ValueError("Failed to assign file roles properly")
            
            self.logger.info(f"K-mer sources: {len(kmer_sources)} files")
            self.logger.info(f"Query targets: {len(query_targets)} files")
            
            # 步骤3: 执行分析 ⚙️
            results = self._execute_analysis(kmer_sources, query_targets)
            
            # 计算运行时间 ⏳
            runtime = time.time() - start_time
            results['runtime_seconds'] = runtime
            results['config'] = self.config
            
            self.logger.info(f"🎉 Auto analysis completed in {runtime:.1f} seconds")
            return results
            
        except Exception as e:
            self.logger.error(f"❌ Auto analysis failed: {e}")
            raise
    
    def explicit_analyze(self, kmer_sources: List[str], query_targets: List[str], 
                        output_dir: Optional[str] = None, 
                        kmer_size: Optional[int] = None) -> Dict:
        """
        明确角色分析模式：用户指定文件角色 👉
        
        Args:
            kmer_sources: k-mer库来源文件路径列表 🧬
            query_targets: 查询目标文件路径列表 🎯
            output_dir: 输出目录 📂
            kmer_size: k-mer长度 📏
        
        Returns:
            分析结果字典 📊
        """
        # 调试：打印配置信息
        self.logger.info(f"🔧 Debug - fastq_pattern: {self.config.fastq_pattern}")
        self.logger.info(f"🔧 Debug - config: {vars(self.config)}")
        
        start_time = time.time()
        
        # 更新配置
        if output_dir:
            self.config.output_dir = output_dir
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        if kmer_size:
            self.config.kmer_size = kmer_size
        
        self.logger.info("🚀 Starting explicit analysis mode")
        
        try:
            # 扫描k-mer源文件 🔍
            self.logger.info(f"🔧 Debug - scanning kmer_sources: {kmer_sources}")
            source_files = self.file_manager.scan_files(kmer_sources, FileRole.KMER_SOURCE)
            self.logger.info(f"🔧 Debug - source_files count: {len(source_files)}")

            self.logger.info(f"🔧 Debug - scanning query_targets: {len(query_targets)} paths")
            target_files = self.file_manager.scan_files(query_targets, FileRole.QUERY_TARGET)
            self.logger.info(f"🔧 Debug - target_files count before grouping: {len(target_files)}")

            if not source_files:
                self.logger.error("❌ No source files found!")
            if not target_files:
                self.logger.error("❌ No target files found!")
                
            if not source_files or not target_files:
                raise ValueError("No valid source or target files found")

            # 只对目标文件应用pattern grouping
            if self.config.fastq_pattern:
                original_target_count = len(target_files)
                target_files = self.file_manager.group_paired_files_by_pattern(target_files)
                self.logger.info(f"Applied fastq pattern '{self.config.fastq_pattern}': {original_target_count} -> {len(target_files)} files")
                
                if len(target_files) == 0:
                    self.logger.error("❌ No target files after pattern grouping!")

            self.logger.info(f"K-mer sources: {len(source_files)} files")
            self.logger.info(f"Query targets: {len(target_files)} files")
                        
            # 执行分析 ⚙️
            results = self._execute_analysis(source_files, target_files)
            
            # 计算运行时间 ⏳
            runtime = time.time() - start_time
            results['runtime_seconds'] = runtime
            results['config'] = self.config
            
            self.logger.info(f"🎉 Explicit analysis completed in {runtime:.1f} seconds")
            return results
            
        except Exception as e:
            self.logger.error(f"❌ Explicit analysis failed: {e}")
            raise
    
    def _execute_analysis(self, source_files: List[FileInfo], 
                         target_files: List[FileInfo]) -> Dict:
        """执行核心分析流程 ❤️‍🔥"""
        
        # 步骤1: 构建k-mer库 📚
        self.logger.info("Step 1: Building k-mer library")
        kmer_library = self._build_kmer_library(source_files)
        
        # 步骤2: 查询目标文件中的k-mer丰度 📊
        self.logger.info("Step 2: Querying k-mer abundances in target files")
        target_abundances = self._query_target_abundances(target_files, kmer_library)
        
        # 步骤3: 生成输出文件 📄
        self.logger.info("Step 3: Generating output files")
        output_files = self._generate_outputs(kmer_library, target_abundances)
        
        return {
            'kmer_library_size': len(kmer_library),
            'target_samples': len(target_abundances),
            'source_files': [f.path for f in source_files],
            'target_files': [f.path for f in target_files],
            'output_files': output_files
        }
    
    def _build_kmer_library(self, source_files: List[FileInfo]) -> Dict[str, KmerInfo]:
        """构建k-mer库 🏗️"""
        all_kmers = {}
        
        # 并行处理源文件 ⚡
        with ThreadPoolExecutor(max_workers=min(len(source_files), self.config.threads)) as executor:
            future_to_file = {}
            
            for file_info in source_files:
                if file_info.format == FileFormat.FASTA:
                    future = executor.submit(
                        self.position_tracker.extract_fasta_kmers,
                        file_info.path, file_info.sample_name
                    )
                elif file_info.format == FileFormat.FASTQ:
                    future = executor.submit(
                        self.position_tracker.extract_fastq_kmers,
                        file_info.path, file_info.sample_name
                    )
                else:
                    continue
                
                future_to_file[future] = file_info
            
            # 收集结果 📥
            for future in as_completed(future_to_file):
                file_info = future_to_file[future]
                try:
                    file_kmers = future.result()
                    # 合并k-mer 🔄
                    for kmer, info in file_kmers.items():
                        if kmer in all_kmers:
                            all_kmers[kmer].positions.extend(info.positions)
                            all_kmers[kmer].total_count += info.total_count
                        else:
                            all_kmers[kmer] = info
                    
                    self.logger.info(f"Processed {file_info.path}: {len(file_kmers)} k-mers ✅")
                    
                except Exception as e:
                    self.logger.error(f"❌ Failed to process {file_info.path}: {e}")
                    raise
        
        self.logger.info(f"🏆 Built k-mer library with {len(all_kmers)} unique k-mers")
        return all_kmers
    
    def _query_target_abundances(self, target_files: List[FileInfo], 
                               kmer_library: Dict[str, KmerInfo]) -> Dict[str, Dict[str, int]]:
        """查询目标文件中的k-mer丰度 🔎"""
        target_abundances = {}
        kmer_list = list(kmer_library.keys())
        
        for file_info in target_files:
            self.logger.info(f"Querying abundances in {file_info.primary_path}")
            
            try:
                # 使用所有文件路径进行KMC计数
                file_format = "fa" if file_info.format == FileFormat.FASTA else "fq"
                kmc_db = self.kmc_interface.count_kmers(
                    file_info.file_paths,  # 使用所有文件路径
                    f"{self.config.temp_dir}/target_{file_info.sample_name}",
                    file_format
                )
                
                # 查询k-mer丰度 🧐
                abundances = self.kmc_interface.query_kmers(kmc_db, kmer_list)
                target_abundances[file_info.sample_name] = abundances
                
                present_count = sum(1 for count in abundances.values() if count > 0)
                self.logger.info(f"Sample {file_info.sample_name}: {present_count}/{len(kmer_list)} k-mers present ✔️")
                
            except Exception as e:
                self.logger.error(f"❌ Failed to query {file_info.path}: {e}")
                raise
        
        return target_abundances
    
    def _generate_outputs(self, kmer_library: Dict[str, KmerInfo], 
                         target_abundances: Dict[str, Dict[str, int]]) -> Dict[str, str]:
        """生成输出文件 🖨️"""
        output_files = {}
        
        # 生成文件名 🏷️
        timestamp = int(time.time())
        base_name = f"kmer_k{self.config.kmer_size}_{timestamp}"
        
        # 1. K-mer库文件 (FASTA格式) 🧬📄
        if "fasta" in self.config.output_formats:
            kmer_lib_file = os.path.join(self.config.output_dir, f"{base_name}_library.fasta")
            self.output_manager.write_kmer_library(kmer_library, kmer_lib_file)
            output_files['kmer_library'] = kmer_lib_file
        
        # 2. 丰度矩阵 (CSV格式) 📊
        if "csv" in self.config.output_formats:
            abundance_file = os.path.join(self.config.output_dir, f"{base_name}_abundance.csv")
            self.output_manager.write_abundance_matrix(kmer_library, target_abundances, abundance_file)
            output_files['abundance_matrix'] = abundance_file
            
            # 3. 存在/缺失矩阵 ✅❌
            presence_file = os.path.join(self.config.output_dir, f"{base_name}_presence.csv")
            self.output_manager.write_presence_matrix(abundance_file, presence_file)
            output_files['presence_matrix'] = presence_file
        
        # 4. 滑窗分析 🪟
        if self.config.window_sizes and any(info.source_type == 'fasta' for info in kmer_library.values()):
            window_file = os.path.join(self.config.output_dir, f"{base_name}_sliding_window.csv")
            self.output_manager.write_sliding_window_analysis(kmer_library, target_abundances, window_file)
            output_files['sliding_window'] = window_file
        
        # 5. 摘要报告 📝
        if "txt" in self.config.output_formats:
            summary_file = os.path.join(self.config.output_dir, f"{base_name}_summary.txt")
            self.output_manager.write_summary_report(kmer_library, target_abundances, summary_file)
            output_files['summary'] = summary_file
        
        return output_files