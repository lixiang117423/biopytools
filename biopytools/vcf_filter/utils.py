"""
VCF筛选工具函数模块 | VCF Filtering Utility Functions Module
"""

import logging
import sys
import subprocess
import shutil
import time
import gzip
import os
from pathlib import Path
from typing import Optional

class FilterLogger:
    """筛选日志管理器 | Filtering Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "vcf_filter.log", verbose: bool = False):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.verbose = verbose
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 清除已有的handlers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        
        # 创建日志器
        self.logger = logging.getLogger('vcf_filter')
        self.logger.setLevel(logging.INFO)
        
        # 文件handler
        if self.log_file.exists():
            self.log_file.unlink()
        
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.INFO)
        
        # 控制台handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        
        # 格式器
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        # 添加handlers
        self.logger.addHandler(file_handler)
        if self.verbose:
            self.logger.addHandler(console_handler)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class PerformanceLogger:
    """性能优化的日志器 | Performance Optimized Logger"""
    
    def __init__(self, verbose=False, output_dir=None):
        self.verbose = verbose
        self.start_time = time.time()
        
        # 如果提供了输出目录，使用标准日志格式
        if output_dir and verbose:
            self.standard_logger = FilterLogger(Path(output_dir), verbose=True).get_logger()
        else:
            self.standard_logger = None
    
    def info(self, msg):
        if self.standard_logger:
            self.standard_logger.info(msg)
        elif self.verbose:
            elapsed = time.time() - self.start_time
            print(f"[{elapsed:.1f}s] INFO: {msg}")
    
    def warning(self, msg):
        if self.standard_logger:
            self.standard_logger.warning(msg)
        else:
            elapsed = time.time() - self.start_time
            print(f"[{elapsed:.1f}s] WARNING: {msg}")
    
    def error(self, msg):
        if self.standard_logger:
            self.standard_logger.error(msg)
        else:
            elapsed = time.time() - self.start_time
            print(f"[{elapsed:.1f}s] ERROR: {msg}")
    
    def debug(self, msg):
        if self.standard_logger:
            self.standard_logger.debug(msg)
        elif self.verbose:
            elapsed = time.time() - self.start_time
            print(f"[{elapsed:.1f}s] DEBUG: {msg}")

def check_dependencies(logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    missing_deps = []
    
    # 检查Python包 | Check Python packages
    optional_packages = [
        ('pandas', 'Pandas'),
        ('numpy', 'NumPy')
    ]
    
    for package, name in optional_packages:
        try:
            __import__(package)
            logger.info(f"✓ {name} 可用 | available")
        except ImportError:
            logger.warning(f"! {name} 不可用 | not available (可选 | optional)")
    
    # 检查外部软件 | Check external software
    optional_tools = [
        ('bcftools', 'BCFtools'),
        ('tabix', 'Tabix'),
        ('bgzip', 'BGzip')
    ]
    
    for tool, name in optional_tools:
        if shutil.which(tool):
            logger.info(f"✓ {name} 可用 | available")
        else:
            logger.warning(f"! {name} 不可用 | not available (推荐安装 | recommended)")
    
    return True

class FastVCFStats:
    """快速VCF统计器 | Fast VCF Statistics"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def quick_count_variants_system(self, vcf_file: str) -> Optional[int]:
        """使用系统命令快速统计变异数量 | Quick count variants using system commands"""
        try:
            if vcf_file.endswith('.gz'):
                cmd = f"zcat {vcf_file} | grep -v '^#' | wc -l"
            else:
                cmd = f"grep -v '^#' {vcf_file} | wc -l"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                count = int(result.stdout.strip())
                self.logger.info(f"系统命令快速统计 | System command quick count: {count} variants")
                return count
        except Exception as e:
            self.logger.debug(f"系统命令统计失败 | System command count failed: {e}")
        return None
    
    def quick_count_variants_bcftools(self, vcf_file: str) -> Optional[int]:
        """使用bcftools快速统计变异数量 | Quick count variants using bcftools"""
        try:
            if shutil.which('bcftools'):
                cmd = f"bcftools view -H {vcf_file} | wc -l"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
                if result.returncode == 0:
                    count = int(result.stdout.strip())
                    self.logger.info(f"bcftools快速统计 | bcftools quick count: {count} variants")
                    return count
        except Exception as e:
            self.logger.debug(f"bcftools统计失败 | bcftools count failed: {e}")
        return None
    
    def sample_based_estimate(self, vcf_file: str) -> Optional[int]:
        """基于采样的变异数量估算 | Sample-based variant count estimation"""
        try:
            opener = gzip.open if vcf_file.endswith('.gz') else open
            sample_count = 0
            header_lines = 0
            
            with opener(vcf_file, 'rt') as f:
                # 统计头部行数 | Count header lines
                for line in f:
                    if line.startswith('#'):
                        header_lines += 1
                    else:
                        sample_count += 1
                        if sample_count >= 1000:  # 采样1000行
                            break
            
            if sample_count > 0:
                # 基于文件大小估算
                file_size = os.path.getsize(vcf_file)
                # 估算每行大小
                sample_size = file_size // (header_lines + sample_count)
                estimated_total_lines = file_size // sample_size
                estimated_variants = estimated_total_lines - header_lines
                
                self.logger.info(f"采样估算 | Sampling estimate: ~{estimated_variants} variants")
                return max(0, estimated_variants)
                
        except Exception as e:
            self.logger.debug(f"采样估算失败 | Sampling estimate failed: {e}")
        return None
    
    def get_variant_count(self, vcf_file: str, skip_count: bool = True) -> Optional[int]:
        """获取变异数量（优先使用快速方法）| Get variant count (prefer fast methods)"""
        
        if skip_count:
            self.logger.info("跳过变异数量统计以提高性能 | Skipping variant count for performance")
            return None
        
        self.logger.info("快速统计变异数量 | Quick variant counting")
        
        # 方法1: 尝试系统命令
        count = self.quick_count_variants_system(vcf_file)
        if count is not None:
            return count
        
        # 方法2: 尝试bcftools
        count = self.quick_count_variants_bcftools(vcf_file)
        if count is not None:
            return count
        
        # 方法3: 采样估算
        count = self.sample_based_estimate(vcf_file)
        if count is not None:
            return count
        
        self.logger.warning("无法快速统计变异数量，将跳过此步骤 | Cannot quickly count variants, skipping this step")
        return None

def check_plink_availability(plink_path: str, logger):
    """检查plink可用性 | Check plink availability"""
    try:
        result = subprocess.run([plink_path, '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info(f"✓ PLINK 可用 | available: {plink_path}")
            return True
        else:
            logger.warning(f"! PLINK 不可用 | not available: {plink_path}")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        logger.warning(f"! PLINK 不可用 | not available: {plink_path}")
        return False
