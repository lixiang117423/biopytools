"""
🛠️ 工具函数模块 | Utility Functions Module
"""

import subprocess
import logging
import os
from pathlib import Path
from typing import Optional

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, dry_run: bool = False):
        self.logger = logger
        self.dry_run = dry_run
    
    def run(self, cmd: str, description: str = "", check_output: bool = False) -> Optional[str]:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"📌 {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        if self.dry_run:
            self.logger.warning("🔍 干运行模式 - 跳过实际执行 | Dry-run mode - Skipping actual execution")
            return None
        
        try:
            if check_output:
                result = subprocess.run(
                    cmd, shell=True, capture_output=True, 
                    text=True, check=True
                )
                return result.stdout.strip()
            else:
                result = subprocess.run(
                    cmd, shell=True, check=True
                )
                self.logger.info(f"✅ 执行成功 | Execution successful")
                return None
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if hasattr(e, 'stderr') and e.stderr:
                self.logger.error(f"错误信息 | Error message: {e.stderr}")
            raise
    
    def run_with_output(self, cmd: str, description: str = "") -> str:
        """执行命令并返回输出 | Execute command and return output"""
        return self.run(cmd, description, check_output=True)

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")
    
    dependencies = [
        # (config.bwa_path, "BWA", "--version"),
        # (config.bwa_path, "BWA"),
        (config.samtools_path, "SAMtools", "--version"),
        (config.gatk_path, "GATK", "--version")
    ]
    
    missing_deps = []
    
    for cmd, name, version_flag in dependencies:
        try:
            result = subprocess.run(
                [cmd, version_flag], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            if result.returncode == 0:
                version = result.stdout.strip().split('\n')[0]
                logger.info(f"  ✅ {name}: {version}")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
            logger.error(f"  ❌ {name}: 未找到 | Not found")
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    logger.info("✅ 所有依赖软件检查完成 | All dependencies checked")
    return True

def check_file_exists(filepath: Path, logger) -> bool:
    """检查文件是否存在 | Check if file exists"""
    if filepath.exists():
        logger.info(f"  ✅ 文件已存在 | File exists: {filepath.name}")
        return True
    return False
