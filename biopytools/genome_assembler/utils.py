"""
基因组组装工具函数模块 🔧 | Genome Assembly Utility Functions Module
"""

import subprocess
import os
from pathlib import Path

def check_dependencies(logger):
    """检查依赖软件 🔍 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    try:
        result = subprocess.run(['hifiasm', '-V'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info(f"✅ hifiasm 可用 | available: {result.stdout.strip()}")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        logger.error("❌ hifiasm 未安装或不在PATH中 | not installed or not in PATH")
        raise RuntimeError("❌ 缺少必需软件: hifiasm | Missing required software: hifiasm")

def run_command(cmd: str, logger, work_dir: str = None, capture_output: bool = False) -> bool:
    """执行命令 ⚡ | Execute command"""
    try:
        logger.info(f"🔧 执行命令 | Executing command: {cmd}")
        
        if work_dir:
            result = subprocess.run(
                cmd, 
                shell=True, 
                cwd=work_dir,
                capture_output=capture_output,
                text=True,
                check=True
            )
        else:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=capture_output,
                text=True,
                check=True
            )
        
        if capture_output and result.stdout:
            logger.info(f"输出 | Output: {result.stdout}")

        logger.info("✅ 命令执行成功 | Command executed successfully")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"❌ 命令执行失败 | Command execution failed: {e}")
        if e.stderr:
            logger.error(f"错误信息 | Error message: {e.stderr}")
        return False

def get_fasta_stats(fasta_file: str) -> dict:
    """获取FASTA文件统计信息 | Get FASTA file statistics"""
    if not os.path.exists(fasta_file):
        return {}
    
    stats = {}
    
    # 计算序列数和总长度
    with open(fasta_file) as f:
        seq_lengths = []
        current_seq = []
        
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    seq_lengths.append(len(''.join(current_seq)))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_seq:
            seq_lengths.append(len(''.join(current_seq)))
    
    stats['num_seqs'] = len(seq_lengths)
    stats['total_len'] = sum(seq_lengths)
    
    # 计算最长序列
    stats['max_len'] = max(seq_lengths) if seq_lengths else 0
    
    # 计算N50
    if seq_lengths:
        sorted_lengths = sorted(seq_lengths, reverse=True)
        total_len = sum(sorted_lengths)
        half_len = total_len / 2
        
        cum_len = 0
        for length in sorted_lengths:
            cum_len += length
            if cum_len >= half_len:
                stats['n50'] = length
                break
    
    return stats

def format_time(seconds: int) -> str:
    """格式化时间 | Format time"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours}小时 {minutes}分钟 {seconds}秒"
