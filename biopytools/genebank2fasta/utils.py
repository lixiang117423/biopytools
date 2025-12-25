"""
GenBank序列提取工具函数模块 🛠️ | GenBank Sequence Extraction Utility Functions Module
"""

import logging
import subprocess
import sys
import os
import glob
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

class ExtractorLogger:
    """序列提取日志管理器 📝 | Sequence Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "genbank_extraction.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 🔊 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 📊 | Get logger"""
        return self.logger

class FileManager:
    """文件管理器 📁 | File Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def find_genbank_files(self):
        """查找GenBank文件 🔍 | Find GenBank files"""
        genbank_files = glob.glob(f"{self.config.gb_dir}/*.gb") + glob.glob(f"{self.config.gb_dir}/*.gbk")
        
        if not genbank_files:
            self.logger.error(f"在 {self.config.gb_dir} 中未找到GenBank文件(.gb/.gbk) 🚫")
            return []
        
        self.logger.info(f"找到 {len(genbank_files)} 个GenBank文件 ✅")
        return genbank_files
    
    def create_output_directories(self):
        """创建输出目录 📂 | Create output directories"""
        directories = [
            self.config.cds_dir,
            self.config.pep_dir,
        ]
        
        if self.config.separate_by_sample:
            directories.extend([
                f"{self.config.cds_dir}/by_sample",
                f"{self.config.pep_dir}/by_sample"
            ])
        
        if self.config.separate_by_gene:
            directories.extend([
                f"{self.config.cds_dir}/by_gene",
                f"{self.config.pep_dir}/by_gene"
            ])
        
        for directory in directories:
            os.makedirs(directory, exist_ok=True)
        
        self.logger.info(f"输出目录结构已创建 🏗️")

def parallel_process_files(file_list, process_func, threads, logger):
    """并行处理文件 ⚡ | Parallel file processing"""
    results = []
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # 提交任务 | Submit tasks
        future_to_file = {executor.submit(process_func, file): file for file in file_list}
        
        # 处理结果 | Process results
        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                result = future.result()
                results.append(result)
                logger.info(f"✅ 处理完成: {os.path.basename(file)}")
            except Exception as e:
                logger.error(f"❌ 处理文件 {file} 时出错: {e}")
    
    return results
