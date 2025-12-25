"""
🛠️ 覆盖度分析工具函数模块 | Depth Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import List

class DepthLogger:
    """📝 覆盖度分析日志管理器 | Depth Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "depth_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """🔧 设置日志 | Setup logging"""
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
        """📋 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """⚡ 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()  # 使用绝对路径 | Use absolute path
    
    def run(self, cmd: str, description: str = "") -> bool:
        """🚀 执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        self.logger.info(f"工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False
    
    def _check_and_create_index(self, bam_file: str) -> bool:
        """🔍 检查并创建BAM索引文件 | Check and create BAM index file"""
        # 检查可能的索引文件位置 | Check possible index file locations
        possible_index_files = [
            f"{bam_file}.bai",      # 标准位置 | Standard location
            f"{bam_file}.csi",      # CSI索引 | CSI index
            bam_file.replace('.bam', '.bai'),  # 同名但不同扩展名 | Same name different extension
        ]
        
        # 检查是否已存在索引 | Check if index already exists
        for index_file in possible_index_files:
            if os.path.exists(index_file):
                self.logger.info(f"✅ 找到索引文件: {index_file} | Found index file: {index_file}")
                return True
        
        # 如果没有索引，创建索引 | If no index exists, create one
        self.logger.info(f"📋 未找到索引文件，正在创建... | Index not found, creating...")
        
        index_cmd = f"{self.config.samtools_path} index {bam_file}"
        
        self.logger.info(f"🔧 创建索引命令: {index_cmd} | Creating index command: {index_cmd}")
        
        try:
            result = subprocess.run(
                index_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.cmd_runner.working_dir
            )
            
            self.logger.info(f"✅ 索引创建成功: {bam_file}.bai | Index created successfully: {bam_file}.bai")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 索引创建失败: {e.stderr} | Index creation failed: {e.stderr}")
            self.logger.warning(f"⚠️ 尝试不使用索引进行分析... | Trying analysis without index...")
            return False  # 返回False但不终止，尝试不使用区间筛选 | Return False but don't terminate, try without region filtering
    
    def _get_window_output_file(self) -> str:
        """📊 获取窗口分析输出文件名 | Get window analysis output filename"""
        base_path = Path(self.config.output_file)
        window_file = base_path.with_name(f"{base_path.stem}_windows_{self.config.window_size}bp{base_path.suffix}")
        return str(window_file)
    
    def _collect_window_data(self, window_data: dict, sample_name: str, chrom: str, pos: int, depth: float):
        """📊 收集窗口数据 | Collect window data"""
        # 计算窗口编号 | Calculate window number
        window_start = ((pos - 1) // self.config.window_step) * self.config.window_step + 1
        window_end = window_start + self.config.window_size - 1
        
        # 创建窗口键 | Create window key
        window_key = (sample_name, chrom, window_start, window_end)
        
        if window_key not in window_data:
            window_data[window_key] = {'depths': [], 'positions': []}
        
        window_data[window_key]['depths'].append(depth)
        window_data[window_key]['positions'].append(pos)
    
    def _write_window_data(self, window_data: dict):
        """📝 写入窗口分析数据 | Write window analysis data"""
        window_output_file = self._get_window_output_file()
        
        with open(window_output_file, 'a') as f:
            for (sample_name, chrom, window_start, window_end), data in window_data.items():
                if data['depths']:
                    avg_depth = sum(data['depths']) / len(data['depths'])
                    window_center = (window_start + window_end) / 2
                    data_points = len(data['depths'])
                    
                    f.write(f"{sample_name}\t{chrom}\t{window_start}\t{window_end}\t{window_center:.1f}\t{avg_depth:.3f}\t{data_points}\n")

def check_dependencies(config, logger):
    """🔍 检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.samtools_path, "samtools")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✅ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True

def get_sample_name(file_path: str) -> str:
    """🏷️ 从文件路径提取样品名称 | Extract sample name from file path"""
    basename = os.path.basename(file_path)
    # 移除常见的文件扩展名 | Remove common file extensions
    for ext in ['.bam', '.sam', '.sorted', '.dedup']:
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
    return basename

def format_region(chromosome: str, start_pos: int = None, end_pos: int = None) -> str:
    """📍 格式化samtools区间参数 | Format samtools region parameter"""
    if start_pos is not None and end_pos is not None:
        return f"{chromosome}:{start_pos}-{end_pos}"
    elif chromosome != 'all':
        return chromosome
    else:
        return ""
