"""
PacBio HiFi结构变异检测工具模块 | PacBio HiFi Structural Variant Detection Utilities Module
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional


class SVCommandRunner:
    """SV分析专用命令执行器 | SV Analysis Command Runner"""
    
    def __init__(self, working_dir: str, logger: logging.Logger):
        self.working_dir = Path(working_dir)
        self.logger = logger
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
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


def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    dependencies = [
        (config.pbsv_path, "pbsv"),
        (config.sniffles_path, "sniffles"),
        (config.cutesv_path, "cuteSV"),
        (config.samtools_path, "samtools"),
        (config.bcftools_path, "bcftools"),
        (config.bgzip_path, "bgzip"),
        (config.tabix_path, "tabix"),
        (config.survivor_path, "SURVIVOR")
    ]
    
    missing_deps = []
    
    for cmd, name in dependencies:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    if missing_deps:
        error_msg = f"缺少依赖软件 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    return True


def get_vcf_stats(vcf_file: Path, logger: logging.Logger) -> Dict[str, int]:
    """获取VCF文件统计信息 | Get VCF file statistics"""
    stats = {}
    
    try:
        # 总SV数量 | Total SV count
        result = subprocess.run(
            f"bcftools view -H {vcf_file} | wc -l",
            shell=True, capture_output=True, text=True, check=True
        )
        stats['total_sv'] = int(result.stdout.strip()) if result.stdout.strip() else 0
        
        # 各类型SV统计 | SV type statistics
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'TRA']
        for sv_type in sv_types:
            result = subprocess.run(
                f"bcftools view -i 'INFO/SVTYPE==\"{sv_type}\"' -H {vcf_file} | wc -l",
                shell=True, capture_output=True, text=True, check=True
            )
            stats[sv_type.lower()] = int(result.stdout.strip()) if result.stdout.strip() else 0
        
    except Exception as e:
        logger.warning(f"获取VCF统计信息失败 | Failed to get VCF stats: {e}")
        stats = {key: 0 for key in ['total_sv', 'del', 'ins', 'dup', 'inv', 'bnd', 'tra']}
    
    return stats
