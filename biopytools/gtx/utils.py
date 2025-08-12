"""
🛠️ GTX WGS分析工具函数模块 | GTX WGS Analysis Utility Functions Module 🛠️
"""

import logging
import subprocess
import sys
import os
import re
from pathlib import Path
from datetime import datetime

class GTXLogger:
    """📝 GTX分析日志管理器 | GTX Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "gtx_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 ✍️ | Setup logging"""
        if self.log_file.exists():
            # 备份现有日志文件 💾 | Backup existing log file
            backup_name = f"gtx_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log.bak"
            self.log_file.rename(self.output_dir / backup_name)
        
        logging.basicConfig(
            level=logging.INFO,
            format='[%(asctime)s] %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 📢 | Get logger"""
        return self.logger

class CommandRunner:
    """🚀 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 ▶️ | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📂 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📜 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📜 标准输出 | Stdout: {e.stdout}")
            return False

class PatternParser:
    """📋 文件模式解析器 | File Pattern Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def extract_sample_name_from_pattern(self, filename: str, pattern: str) -> str:
        """根据模式从文件名提取样品名 | Extract sample name from filename based on pattern"""
        import re
        
        try:
            # 先转义所有特殊字符，但保留*作为占位符
            escaped_pattern = re.escape(pattern)
            # 现在把转义后的\*替换为捕获组
            regex_pattern = escaped_pattern.replace(r'\*', '(.+?)')
            # 添加行首行尾锚点
            regex_pattern = f"^{regex_pattern}$"
            
            self.logger.debug(f"模式转换 | Pattern conversion: {pattern} -> {regex_pattern}")
            
            match = re.match(regex_pattern, filename)
            if match:
                sample_name = match.group(1)
                self.logger.debug(f"成功提取样品名 | Successfully extracted sample name: {filename} -> {sample_name}")
                return sample_name
            else:
                self.logger.debug(f"文件名 {filename} 不匹配模式 {pattern}")
                return None
        except Exception as e:
            self.logger.error(f"模式匹配错误 | Pattern matching error: {e}")
            return None
    
    def build_paired_filename(self, sample_name: str, pattern: str) -> str:
        """根据样品名和模式构建配对文件名 | Build paired filename based on sample name and pattern"""
        return pattern.replace("*", sample_name)

class FileProcessor:
    """📁 文件处理器 | File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.pattern_parser = PatternParser(logger)
    
    def find_fastq_files(self):
        """查找FASTQ文件对 🔍 | Find FASTQ file pairs"""
        self.logger.info("🔍 搜索输入文件 | Searching input files")
        
        input_path = Path(self.config.input_dir)
        r1_files = list(input_path.glob(self.config.read1_pattern))
        
        if not r1_files:
            raise FileNotFoundError(f"在 {self.config.input_dir} 中未找到 {self.config.read1_pattern} 文件")
        
        # 按文件名排序 🔠 | Sort by filename
        r1_files.sort()
        
        self.logger.info(f"找到 {len(r1_files)} 个R1文件 | Found {len(r1_files)} R1 files")
        for r1_file in r1_files[:5]:  # 显示前5个文件作为示例
            self.logger.info(f"  R1文件示例 | R1 file example: {r1_file.name}")
        if len(r1_files) > 5:
            self.logger.info(f"  ... 还有 {len(r1_files) - 5} 个文件 | ... and {len(r1_files) - 5} more files")
        
        # 验证R2文件存在 ✅ | Validate R2 files exist
        file_pairs = []
        for r1_file in r1_files:
            # 从R1文件名中提取样品名 🏷️ | Extract sample name from R1 filename
            sample_name = self.pattern_parser.extract_sample_name_from_pattern(
                r1_file.name, self.config.read1_pattern
            )
            
            if not sample_name:
                self.logger.warning(f"⚠️ 无法从文件名提取样品名 | Cannot extract sample name from filename: {r1_file.name}")
                continue
            
            # 构建R2文件路径 🏗️ | Build R2 file path
            r2_filename = self.pattern_parser.build_paired_filename(sample_name, self.config.read2_pattern)
            r2_file = input_path / r2_filename
            
            if not r2_file.exists():
                self.logger.warning(f"⚠️ 找不到对应的R2文件 | Cannot find corresponding R2 file: {r2_file}")
                continue
            
            file_pairs.append((sample_name, str(r1_file), str(r2_file)))
            self.logger.debug(f"✅ 找到配对文件 | Found paired files: {sample_name}")
        
        total_samples = len(file_pairs)
        self.logger.info(f"✅ 找到 {total_samples} 个样品需要处理 | Found {total_samples} samples to process")
        
        # 显示样品名示例 | Show sample name examples
        if total_samples > 0:
            self.logger.info("样品列表示例 | Sample list examples:")
            for i, (sample_name, r1_file, r2_file) in enumerate(file_pairs[:3]):
                self.logger.info(f"  {i+1}. {sample_name}")
                self.logger.info(f"     R1: {Path(r1_file).name}")
                self.logger.info(f"     R2: {Path(r2_file).name}")
            if total_samples > 3:
                self.logger.info(f"  ... 还有 {total_samples - 3} 个样品 | ... and {total_samples - 3} more samples")
        
        return file_pairs
    
    def check_output_exists(self, sample_name: str) -> bool:
        """检查输出文件是否已存在 🧐 | Check if output files already exist"""
        vcf_file = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
        bam_file = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        return vcf_file.exists() and bam_file.exists()
    
    def get_file_size(self, file_path: str) -> str:
        """获取文件大小 📏 | Get file size"""
        try:
            size_bytes = os.path.getsize(file_path)
            # 转换为人类可读格式 🧑‍💻 | Convert to human readable format
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size_bytes < 1024.0:
                    return f"{size_bytes:.1f} {unit}"
                size_bytes /= 1024.0
            return f"{size_bytes:.1f} TB"
        except:
            return "Unknown"

def check_dependencies(config, logger):
    """检查依赖软件 🧩 | Check dependencies"""
    logger.info("🧩 检查依赖软件 | Checking dependencies")
    
    # 检查GTX程序 💻 | Check GTX program
    if not os.path.exists(config.gtx_path):
        error_msg = f"❌ GTX程序不存在 | GTX program does not exist: {config.gtx_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    # 检查GTX程序是否可执行 🏃 | Check if GTX program is executable
    if not os.access(config.gtx_path, os.X_OK):
        error_msg = f"❌ GTX程序不可执行 | GTX program is not executable: {config.gtx_path}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    logger.info("✅ ✓ GTX程序检查通过 | GTX program check passed")
    return True
