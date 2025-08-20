"""
📁 文件处理模块 | File Processing Module
"""

import os
import re
import shutil
import glob
import gzip
import subprocess
from pathlib import Path
from typing import List

class FileProcessor:
    """📂 文件处理器 | File Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    

    # ln -s软链接识别
    def find_fastq_files(self) -> List[tuple]:
        """🔍 查找FASTQ文件 | Find FASTQ files"""
        self.logger.info(f"🔍 查找FASTQ文件 | Finding FASTQ files")
        self.logger.info(f"📁 输入目录: {self.config.input_dir}")
        self.logger.info(f"🔍 文件模式: {self.config.pattern}")
        
        # 构建完整搜索路径
        input_dir = Path(self.config.input_dir).resolve()
        
        self.logger.info(f"📂 输入目录是否存在: {input_dir.exists()}")
        
        # 列出目录内容用于调试
        if input_dir.exists():
            files_in_dir = list(input_dir.glob("*"))
            self.logger.info(f"📄 目录中的文件: {[f.name for f in files_in_dir[:10]]}")  # 只显示前10个
            
            # 🔥 检查软链接情况 | Check symlink status
            symlinks = [f for f in files_in_dir if f.is_symlink()]
            if symlinks:
                self.logger.info(f"🔗 发现{len(symlinks)}个软链接文件 | Found {len(symlinks)} symlink files")
                for link in symlinks[:5]:  # 显示前5个
                    target = link.readlink() if link.is_symlink() else "N/A"
                    target_exists = link.exists()
                    self.logger.info(f"  🔗 {link.name} -> {target} (目标存在: {target_exists})")
        
        # 🔥 使用Path.glob替代glob.glob，更好地处理软链接 | Use Path.glob instead of glob.glob for better symlink handling
        r1_files = list(input_dir.glob(self.config.pattern))
        
        # 🔥 过滤和验证文件 | Filter and validate files
        valid_r1_files = []
        for file_path in r1_files:
            if file_path.is_symlink():
                # 检查软链接目标是否存在 | Check if symlink target exists
                if file_path.exists():
                    self.logger.info(f"✅ 软链接文件有效 | Valid symlink: {file_path.name}")
                    valid_r1_files.append(str(file_path))
                else:
                    self.logger.warning(f"⚠️ 软链接目标不存在 | Symlink target missing: {file_path.name} -> {file_path.readlink()}")
            elif file_path.is_file():
                self.logger.info(f"✅ 常规文件 | Regular file: {file_path.name}")
                valid_r1_files.append(str(file_path))
            else:
                self.logger.warning(f"⚠️ 无效文件 | Invalid file: {file_path.name}")
        
        self.logger.info(f"🔍 找到{len(valid_r1_files)}个有效R1文件")
        
        if not valid_r1_files:
            raise FileNotFoundError(f"📂 未找到匹配的文件 | No files found matching pattern: {self.config.pattern}")
        
        samples = []
        for r1_file in valid_r1_files:
            # 2. 从R1文件提取样本名（*的部分）
            r1_filename = os.path.basename(r1_file)
            
            # 提取样本名：去掉模式中*后面的部分
            pattern_suffix = self.config.pattern.split('*')[1]  # 例如："_1.clean.fq.gz"
            sample_name = r1_filename.replace(pattern_suffix, "")  # 例如："OV8-105"
            
            self.logger.info(f"📝 提取样本名: {sample_name}")
            
            # 3. 查找对应的R2文件
            r2_pattern_suffix = pattern_suffix.replace('_1.', '_2.').replace('_R1.', '_R2.')
            r2_filename = sample_name + r2_pattern_suffix
            r2_file = input_dir / r2_filename
            
            # 🔥 改进R2文件检查，支持软链接 | Improved R2 file check with symlink support
            if r2_file.exists():
                if r2_file.is_symlink():
                    self.logger.info(f"✅ R2软链接文件有效 | Valid R2 symlink: {r2_file.name}")
                samples.append((sample_name, r1_file, str(r2_file)))
                self.logger.info(f"👥 找到双端样本: {sample_name}")
                self.logger.info(f"  📄 R1: {r1_file}")
                self.logger.info(f"  📄 R2: {r2_file}")
            else:
                samples.append((sample_name, r1_file, None))
                self.logger.info(f"📄 找到单端样本: {sample_name}")
                self.logger.info(f"  📄 R1: {r1_file}")
                self.logger.info(f"  ❌ R2不存在: {r2_file}")
        
        if not samples:
            raise FileNotFoundError("❌ 未找到有效的FASTQ文件 | No valid FASTQ files found")
        
        self.config.samples = samples
        return samples
        
    #     return decompressed_files
    def decompress_files(self, sample_files: List[str]) -> List[str]:
        """📦 准备文件列表（不再解压缩）| Prepare file list (no more decompression)"""
        self.logger.info("📦 准备文件列表 | Preparing file list")
        
        # 🔥 直接返回原始文件路径，让Jellyfish通过generator处理压缩文件
        # Return original file paths directly, let Jellyfish handle compressed files via generator
        processed_files = []
        for file_path in sample_files:
            if file_path is None:
                processed_files.append(None)
            else:
                if file_path.endswith('.gz'):
                    self.logger.info(f"📦 压缩文件将通过generator处理 | Compressed file will be processed via generator: {os.path.basename(file_path)}")
                else:
                    self.logger.info(f"📋 普通文件 | Regular file: {os.path.basename(file_path)}")
                processed_files.append(file_path)
        
        return processed_files
    
    def prepare_kmer_library(self) -> Path:
        """🧬 准备k-mer库文件 | Prepare k-mer library file"""
        self.logger.info("🧬 准备k-mer库文件 | Preparing k-mer library file")
        
        # 读取原始文件并转换为大写 | Read original file and convert to uppercase
        upper_kmer_file = self.config.temp_dir / "kmers_upper.fasta"
        
        with open(self.config.kmer_lib, 'r') as f_in:
            with open(upper_kmer_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        f_out.write(line)
                    else:
                        f_out.write(line.upper())
        
        self.logger.info(f"✅ K-mer库文件已准备完成 | K-mer library file prepared: {upper_kmer_file}")
        return upper_kmer_file
    
    
    def check_file_integrity(self, file_path: str) -> bool:
        """🔍 检查文件完整性 | Check file integrity"""
        if not os.path.exists(file_path):
            return False
        
        if file_path.endswith('.gz'):
            try:
                # 使用gzip -t进行快速完整性检查
                result = subprocess.run(['gzip', '-t', file_path], 
                                      capture_output=True, 
                                      text=True, 
                                      timeout=30)
                return result.returncode == 0
            except (subprocess.TimeoutExpired, FileNotFoundError):
                return False
            except Exception:
                return False
        else:
            try:
                # 普通文件快速读取测试
                with open(file_path, 'r') as f:
                    f.read(1024)  # 只读取前1KB
                return True
            except Exception:
                return False
    
    def recheck_files_integrity(self) -> List[tuple]:
        """🔄 重新检查文件完整性并返回有效样本 | Recheck file integrity and return valid samples"""
        self.logger.warning("⚠️ 检测到文件完整性问题，开始重新验证所有文件...")
        
        if not hasattr(self.config, 'samples') or not self.config.samples:
            self.logger.error("❌ 没有样本信息可供重新检查")
            return []
        
        valid_samples = []
        corrupted_files = []
        
        for sample_info in self.config.samples:
            sample_name, r1_file, r2_file = sample_info
            
            # 检查R1文件
            self.logger.info(f"🔍 检查R1文件: {os.path.basename(r1_file)}")
            r1_valid = self.check_file_integrity(r1_file)
            if not r1_valid:
                corrupted_files.append(os.path.basename(r1_file))
                self.logger.error(f"💥 R1文件损坏，跳过样本: {sample_name} ({os.path.basename(r1_file)})")
                continue
            
            # 检查R2文件（如果存在）
            if r2_file:
                self.logger.info(f"🔍 检查R2文件: {os.path.basename(r2_file)}")
                r2_valid = self.check_file_integrity(r2_file)
                if not r2_valid:
                    corrupted_files.append(os.path.basename(r2_file))
                    self.logger.error(f"💥 R2文件损坏，跳过样本: {sample_name} ({os.path.basename(r2_file)})")
                    continue
                
                valid_samples.append((sample_name, r1_file, r2_file))
                self.logger.info(f"✅ 双端样本验证通过: {sample_name}")
            else:
                valid_samples.append((sample_name, r1_file, None))
                self.logger.info(f"✅ 单端样本验证通过: {sample_name}")
        
        if corrupted_files:
            self.logger.warning(f"⚠️ 发现{len(corrupted_files)}个损坏的文件:")
            for corrupted_file in corrupted_files:
                self.logger.warning(f"  💔 {corrupted_file}")
        
        if valid_samples:
            self.logger.info(f"🎉 完整性检查完成，{len(valid_samples)}个样本可用")
            self.config.samples = valid_samples
        else:
            self.logger.error("❌ 完整性检查后没有可用的样本")
        
        return valid_samples
