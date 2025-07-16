"""
PLINK GWAS分析工具函数模块 | PLINK GWAS Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path

class PlinkGWASLogger:
    """PLINK GWAS分析日志管理器 | PLINK GWAS Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "plink_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
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
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd, description: str = "", check: bool = True):
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        self.logger.info(f"命令 | Command: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=isinstance(cmd, str),
                capture_output=True, 
                text=True, 
                check=check,
                cwd=self.working_dir
            )
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"标准错误 | Stderr: {result.stderr}")
            
            if description:
                self.logger.info(f"完成 | Completed: {description}")
            
            return result
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command failed: {cmd_str}")
            self.logger.error(f"错误信息 | Error: {e}")
            if check:
                raise
            return e

class FileManager:
    """文件管理器 | File Manager"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def copy_input_files(self, vcf_file: str, phenotype_file: str):
        """复制输入文件到工作目录 | Copy input files to working directory"""
        self.logger.info("复制输入文件到工作目录 | Copying input files to working directory...")
        
        # 复制VCF文件 | Copy VCF file
        vcf_dest = self.working_dir / "input.vcf.gz"
        shutil.copy2(vcf_file, vcf_dest)
        self.logger.info(f"VCF文件已复制 | VCF file copied: {vcf_dest}")
        
        # 复制表型文件 | Copy phenotype file
        pheno_dest = self.working_dir / "phenotype.txt"
        shutil.copy2(phenotype_file, pheno_dest)
        self.logger.info(f"表型文件已复制 | Phenotype file copied: {pheno_dest}")
        
        return vcf_dest, pheno_dest
    
    def check_file_exists(self, file_path: Path, description: str = "") -> bool:
        """检查文件是否存在 | Check if file exists"""
        if file_path.exists():
            if description:
                self.logger.info(f"✓ {description}存在 | exists: {file_path}")
            return True
        else:
            if description:
                self.logger.error(f"✗ {description}不存在 | does not exist: {file_path}")
            return False
