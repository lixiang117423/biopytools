"""
ALLHiC流水线工具函数模块 🔧 | ALLHiC Pipeline Utility Functions Module
"""

import os
import subprocess
import sys
from pathlib import Path

def check_dependencies(config, logger):
    """检查依赖软件 🔍 | Check dependencies"""
    logger.log_info("检查依赖软件 | Checking dependencies")

    missing_software = []

    # 检查 BWA (使用完整路径，BWA不支持--version参数)
    bwa_cmd = os.path.join(config.bwa_path, "bwa")
    if os.path.exists(bwa_cmd) and os.access(bwa_cmd, os.X_OK):
        try:
            # BWA不支持--version，使用无参数调用来检测
            result = subprocess.run([bwa_cmd],
                                  capture_output=True, text=True, timeout=10)
            # BWA在无参数时会显示帮助信息，返回码不为0但说明程序可运行
            if "Program: bwa" in result.stderr or "Usage: bwa" in result.stderr:
                logger.log_success(f"BWA 可用 | available")
            else:
                missing_software.append("BWA")
        except subprocess.TimeoutExpired:
            missing_software.append("BWA")
    else:
        missing_software.append("BWA")

    # 检查 samtools (假设在PATH中)
    try:
        result = subprocess.run(["samtools", '--version'],
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.log_success(f"Samtools 可用 | available")
        else:
            missing_software.append("Samtools")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        missing_software.append("Samtools")

    # 检查 ALLHiC 主要工具
    allhic_cmd = os.path.join(config.allhic_software_path, "allhic")
    if os.path.exists(allhic_cmd) and os.access(allhic_cmd, os.X_OK):
        logger.log_success(f"ALLHiC 可用 | available")
    else:
        missing_software.append("ALLHiC")

    # 检查 ALLHiC 其他工具 (在bin目录中)
    allhic_tools = [
        ("ALLHiC_prune", "ALLHiC_prune"),
        ("ALLHiC_partition", "ALLHiC_partition"),
        ("ALLHiC_rescue", "ALLHiC_rescue"),
        ("ALLHiC_build", "ALLHiC_build"),
        ("ALLHiC_plot", "ALLHiC_plot")
    ]

    for tool, name in allhic_tools:
        tool_cmd = os.path.join(config.allhic_bin_path, tool)
        if os.path.exists(tool_cmd) and os.access(tool_cmd, os.X_OK):
            logger.log_success(f"{name} 可用 | available")
        else:
            missing_software.append(name)

    # 检查 minimap2 (假设在PATH中)
    try:
        result = subprocess.run(["minimap2", '--version'],
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.log_success(f"Minimap2 可用 | available")
        else:
            missing_software.append("Minimap2")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        missing_software.append("Minimap2")

    # 检查 asmkit
    if os.path.exists(config.asmkit_path) and os.access(config.asmkit_path, os.X_OK):
        logger.log_success("asmkit 可用 | available")
    else:
        missing_software.append("asmkit")

    if missing_software:
        error_msg = f"❌ 缺少必需软件 | Missing required software: {', '.join(missing_software)}"
        logger.log_error(error_msg)

        # 总是继续执行，让用户在运行时看到实际错误
        logger.log_warning("⚠️ 继续执行 - 如果软件实际存在但未检测到，可能在运行时正常工作")
        logger.log_warning("⚠️ Continuing execution - If software exists but wasn't detected, it may work at runtime")

    return True

def link_file(src: str, dest: str, logger) -> bool:
    """链接文件 | Link file"""
    if os.path.exists(src):
        # 如果目标文件已存在，先删除
        if os.path.exists(dest):
            try:
                os.remove(dest)
                logger.log_info(f"删除已存在的目标文件 | Removing existing target file: {dest}")
            except Exception as e:
                logger.log_error(f"删除目标文件失败 | Failed to remove target file: {e}")
                return False

        try:
            os.symlink(os.path.realpath(src), dest)
            logger.log_info(f"创建符号链接 | Created symlink: {dest} -> {src}")
            return True
        except OSError as e:
            logger.log_warning(f"链接失败，尝试复制 | Link failed, trying copy: {e}")
            try:
                # 如果是同一文件，跳过复制
                if os.path.samefile(src, dest):
                    logger.log_info(f"源文件和目标文件相同，跳过 | Source and target are the same file, skipping")
                    return True
                import shutil
                shutil.copy2(src, dest)
                logger.log_info(f"复制文件成功 | File copied: {src} -> {dest}")
                return True
            except Exception as e2:
                logger.log_error(f"复制失败 | Copy failed: {e2}")
                return False
    else:
        logger.log_warning(f"源文件不存在 | Source file not found: {src}")
        return False

def check_file_exists(file_path: str) -> bool:
    """检查文件是否存在 | Check if file exists"""
    return os.path.isfile(file_path)

def setup_conda_environment(logger):
    """设置conda环境 | Setup conda environment"""
    conda_script = "/share/org/YZWL/yzwl_lixg/miniforge3/etc/profile.d/conda.sh"
    if os.path.exists(conda_script):
        try:
            with open(conda_script) as f:
                exec(f.read(), globals())
            
            import subprocess
            result = subprocess.run(["conda", "activate", "biopytools"], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logger.log_info("Conda环境已激活 | Conda environment activated")
            else:
                logger.log_warning("Conda环境激活失败 | Failed to activate conda environment")
        except Exception as e:
            logger.log_warning(f"Conda设置失败 | Conda setup failed: {e}")
    else:
        logger.log_warning("Conda脚本未找到 | Conda script not found")

def clean_work_directory(work_dir: str, logger):
    """清理工作目录中的链接文件 | Clean symlink files in work directory"""
    if os.path.exists(work_dir):
        logger.log_info(f"清理工作目录 | Cleaning work directory: {work_dir}")
        removed_count = 0

        # 遍历工作目录，删除符号链接和普通文件（保留子目录）
        for item in os.listdir(work_dir):
            item_path = os.path.join(work_dir, item)
            if os.path.isfile(item_path) or os.path.islink(item_path):
                try:
                    os.remove(item_path)
                    removed_count += 1
                except Exception as e:
                    logger.log_warning(f"删除文件失败 | Failed to remove file {item}: {e}")

        if removed_count > 0:
            logger.log_info(f"清理完成，删除了 {removed_count} 个文件 | Cleaned {removed_count} files")
        else:
            logger.log_info("目录已经是空的 | Directory is already empty")

def format_file_size(file_path: str) -> str:
    """格式化文件大小 | Format file size"""
    if os.path.exists(file_path):
        size_bytes = os.path.getsize(file_path)
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f}{unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f}TB"
    return "0B"
