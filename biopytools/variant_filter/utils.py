"""
🔧 变异筛选工具函数模块 | Variant Filtering Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import Optional

class FilterLogger:
    """📝 筛选日志管理器 | Filtering Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "variant_filter.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """⚙️ 设置日志 | Setup logging"""
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
        """📝 获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """💻 命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path, temp_dir: Optional[str] = None):
        self.logger = logger
        self.working_dir = working_dir.resolve()
        self.temp_dir = temp_dir
    
    def run(self, cmd: str, description: str = "") -> bool:
        """🚀 执行命令 | Execute command"""
        if description:
            self.logger.info(f"🔄 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        # ⚙️ 设置环境变量 | Set environment variables
        env = os.environ.copy()
        if self.temp_dir:
            env['TMPDIR'] = self.temp_dir
            env['TMP'] = self.temp_dir
            env['TEMP'] = self.temp_dir
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                env=env
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📄 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            self.logger.error(f"📄 标准输出 | Stdout: {e.stdout}")
            return False

def check_dependencies(config, logger):
    """🔗 检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    # 🔧 基础依赖 | Basic dependencies
    dependencies = [
        (config.vcftools_path, "VCFtools"),
        (config.bcftools_path, "BCFtools"),
        (config.bgzip_path, "bgzip"),
        (config.tabix_path, "tabix")
    ]
    
    # 🧬 如果使用GATK选择变异，检查GATK | Check GATK if using it for variant selection
    if config.use_gatk_select:
        dependencies.append((config.gatk_path, "GATK"))
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            if name == "GATK":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "VCFtools":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            elif name == "BCFtools":
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            else:
                result = subprocess.run([cmd, "--version"], 
                                      capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                logger.info(f"✅ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"❌ 缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True

def get_vcf_stats(vcf_file: str) -> dict:
    """📊 获取VCF文件统计信息 | Get VCF file statistics"""
    stats = {'samples': 0, 'variants': 0}
    
    try:
        if vcf_file.endswith('.gz'):
            with subprocess.Popen(['zcat', vcf_file], stdout=subprocess.PIPE, text=True) as proc:
                for line in proc.stdout:
                    if line.startswith('#CHROM'):
                        stats['samples'] = len(line.strip().split('\t')) - 9
                    elif not line.startswith('#'):
                        stats['variants'] += 1
        else:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#CHROM'):
                        stats['samples'] = len(line.strip().split('\t')) - 9
                    elif not line.startswith('#'):
                        stats['variants'] += 1
    except Exception as e:
        print(f"❌ 获取VCF统计信息失败 | Failed to get VCF statistics: {e}")
    
    return stats
