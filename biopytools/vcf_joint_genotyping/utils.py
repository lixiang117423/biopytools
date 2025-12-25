"""
VCF联合分型工具模块 | VCF Joint Genotyping Utilities Module
"""

import subprocess
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Any

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger: logging.Logger, working_dir: str = None):
        self.logger = logger
        self.working_dir = working_dir or str(Path.cwd())
    
    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """⚡ 执行命令 | Execute command"""
        if description:
            self.logger.info(f"⚡ 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"🖥️ 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                timeout=timeout
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"⚠️ 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"⏰ 命令执行超时 | Command execution timeout: {description}")
            return False

class DependencyChecker:
    """依赖检查器 | Dependency Checker"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
    
    def check_software(self, software_path: str, software_name: str) -> bool:
        """检查单个软件 | Check single software"""
        try:
            # 尝试版本检查 | Try version check
            result = subprocess.run([software_path, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.logger.info(f"✓ {software_name} 可用 | available")
                return True
            
            # 尝试帮助信息 | Try help information
            result = subprocess.run([software_path, "-h"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.logger.info(f"✓ {software_name} 可用 | available")
                return True
            
            self.logger.warning(f"⚠ {software_name} 可能不可用 | may not be available")
            return False
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            self.logger.error(f"✗ {software_name} 不可用 | not available")
            return False
    
    def check_all_dependencies(self, config) -> Tuple[bool, List[str]]:
        """🔧 检查所有依赖 | Check all dependencies"""
        self.logger.info("🔧 检查依赖软件 | Checking dependencies")
        
        dependencies = [
            (config.gtx_path, "GTX"),
            (config.gatk_path, "GATK"),
            (config.vcftools_path, "VCFtools"),
            (config.bgzip_path, "bgzip"),
            (config.tabix_path, "tabix"),
            (config.bcftools_path, "BCFtools")
        ]
        
        missing_deps = []
        
        for software_path, software_name in dependencies:
            if not self.check_software(software_path, software_name):
                missing_deps.append(software_name)
        
        if missing_deps:
            self.logger.error(f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}")
            return False, missing_deps
        
        return True, []

class VCFUtils:
    """VCF工具类 | VCF Utilities"""
    
    @staticmethod
    def count_vcf_files(vcf_dir: str) -> int:
        """统计VCF文件数量 | Count VCF files"""
        vcf_files = list(Path(vcf_dir).glob("*.vcf.gz"))
        return len(vcf_files)
    
    @staticmethod
    def get_vcf_files(vcf_dir: str) -> List[Path]:
        """获取VCF文件列表 | Get VCF file list"""
        return list(Path(vcf_dir).glob("*.vcf.gz"))
    
    @staticmethod
    def count_variants_in_vcf(vcf_file: Path) -> int:
        """统计VCF文件中的变异数量 | Count variants in VCF file"""
        try:
            if vcf_file.suffix == '.gz':
                cmd = f"zcat {vcf_file} | grep -v '^#' | wc -l"
            else:
                cmd = f"grep -v '^#' {vcf_file} | wc -l"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return int(result.stdout.strip())
        except:
            return 0
    
    @staticmethod
    def get_file_size(file_path: Path) -> str:
        """获取文件大小 | Get file size"""
        if file_path.exists():
            try:
                result = subprocess.run(f"du -h {file_path}", shell=True, capture_output=True, text=True)
                return result.stdout.split()[0]
            except:
                return "未知 | Unknown"
        return "不存在 | Not exists"

def setup_logging(verbose: bool = False, log_file: str = None) -> logging.Logger:
    """🔧 设置日志记录 | Setup logging"""
    level = logging.DEBUG if verbose else logging.INFO
    
    # 创建logger
    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    
    # 清除可能已存在的handler，避免重复日志
    if logger.handlers:
        logger.handlers.clear()
    
    # 创建格式化器
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # 创建控制台处理器，明确指定输出到stdout
    import sys
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # 如果指定了日志文件，同时输出到文件
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        # 输出日志文件路径信息
        logger.info(f"📁 日志文件 | Log file: {log_file}")
    
    # 防止日志传播到root logger（避免重复输出）
    logger.propagate = False
    
    return logger
