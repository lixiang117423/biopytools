"""
🔧 K-mer分析工具函数模块 | K-mer Analysis Utility Functions Module
"""

import os
import gzip
import subprocess
from pathlib import Path
from typing import List

class FileFormatDetector:
    """文件格式检测器 | File Format Detector"""
    
    @staticmethod
    def detect_format(file_path: str) -> str:
        """
        检测文件格式 | Detect file format
        
        Returns:
            'fasta', 'fastq', 'unknown'
        """
        file_path = Path(file_path)
        
        # 检查是否压缩 | Check if compressed
        is_gzipped = file_path.suffix == '.gz'
        
        # 打开文件 | Open file
        try:
            if is_gzipped:
                with gzip.open(file_path, 'rt') as f:
                    first_char = f.read(1)
            else:
                with open(file_path, 'r') as f:
                    first_char = f.read(1)
            
            # 根据第一个字符判断 | Determine based on first character
            if first_char == '>':
                return 'fasta'
            elif first_char == '@':
                return 'fastq'
            else:
                return 'unknown'
        except Exception:
            return 'unknown'
    
    @staticmethod
    def get_all_files(path: str, extensions: List[str] = None) -> List[str]:
        """
        获取目录下所有指定格式的文件 | Get all files with specified extensions
        
        Args:
            path: 文件或目录路径 | File or directory path
            extensions: 扩展名列表 | Extension list
        
        Returns:
            文件路径列表 | List of file paths
        """
        path = Path(path)
        
        if path.is_file():
            return [str(path)]
        
        if extensions is None:
            extensions = ['.fasta', '.fa', '.fna', '.fastq', '.fq', '.gz']
        
        files = []
        for ext in extensions:
            files.extend([str(f) for f in path.glob(f'*{ext}')])
            files.extend([str(f) for f in path.glob(f'**/*{ext}')])
        
        return list(set(files))


class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🔄 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False


def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.metagraph_path, "MetaGraph", "--version"),
        (config.kmc_path, "KMC", ""),
        (config.kmc_tools_path, "KMC-tools", ""),
    ]
    
    missing_deps = []
    
    for cmd, name, version_flag in dependencies:
        try:
            if version_flag:
                result = subprocess.run(
                    [cmd, version_flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            else:
                result = subprocess.run(
                    [cmd],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            logger.info(f"✅ {name} 可用 | available")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            logger.warning(f"⚠️  {name} 未找到 | not found: {cmd}")
            missing_deps.append(name)
    
    # 检查Python依赖 | Check Python dependencies
    try:
        import pandas
        logger.info(f"✅ pandas 可用 | available")
    except ImportError:
        missing_deps.append("pandas")
        logger.warning(f"⚠️  pandas 未安装 | not installed")
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        logger.info("💡 提示 | Hint: conda install -c bioconda metagraph kmc pandas")
        raise RuntimeError(error_msg)
    
    return True
