"""
文件管理和识别模块 | File Management and Detection Module 📂🔍
"""

import os
import gzip
import magic
import glob
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import logging

from .config import FileRole, AssignmentStrategy

class FileFormat(Enum):
    """文件格式枚举 📄"""
    FASTA = "fasta"
    FASTQ = "fastq"
    UNKNOWN = "unknown"

# @dataclass
# class FileInfo:
#     """文件信息类 ℹ️"""
#     path: str
#     format: FileFormat
#     role: FileRole
#     size_bytes: int
#     sample_name: str
#     is_compressed: bool = False
#     estimated_sequences: Optional[int] = None
#     estimated_kmers: Optional[int] = None
@dataclass
class FileInfo:
    """文件信息类"""
    path: Union[str, List[str]]  # 修改为支持文件路径列表
    format: FileFormat
    role: FileRole
    size_bytes: int
    sample_name: str
    is_compressed: bool = False
    estimated_sequences: Optional[int] = None
    estimated_kmers: Optional[int] = None
    
    @property
    def file_paths(self) -> List[str]:
        """获取所有文件路径"""
        if isinstance(self.path, list):
            return self.path
        return [self.path]
    
    @property
    def primary_path(self) -> str:
        """获取主要文件路径（用于显示）"""
        if isinstance(self.path, list):
            return self.path[0]
        return self.path

class FileManager:
    """文件管理器 🗂️"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # 文件扩展名映射 🗺️
        self.fasta_extensions = {'.fa', '.fasta', '.fas', '.fna', '.ffn', '.faa', '.frn'}
        self.fastq_extensions = {'.fq', '.fastq'}
        self.compressed_extensions = {'.gz', '.bz2', '.xz'}
    
    def scan_files(self, input_paths: List[str], default_role: FileRole = FileRole.AUTO_DETECT) -> List[FileInfo]:
        """扫描和识别文件"""
        self.logger.info(f"🔧 Debug - scan_files called with {len(input_paths)} paths, role: {default_role}")
        
        all_files = []
        
        for path_pattern in input_paths:
            self.logger.info(f"🔧 Debug - processing path: {path_pattern}")
            
            # 支持通配符和目录
            if os.path.isdir(path_pattern):
                self.logger.info(f"🔧 Debug - {path_pattern} is a directory")
                # 目录：扫描所有可能的序列文件
                pattern_files = self._scan_directory(path_pattern)
            elif '*' in path_pattern or '?' in path_pattern:
                self.logger.info(f"🔧 Debug - {path_pattern} contains wildcards")
                # 通配符模式
                pattern_files = glob.glob(path_pattern)
                self.logger.info(f"🔧 Debug - glob found {len(pattern_files)} files")
            elif os.path.isfile(path_pattern):
                self.logger.info(f"🔧 Debug - {path_pattern} is a single file")
                # 单个文件
                pattern_files = [path_pattern]
            else:
                self.logger.warning(f"Path not found or invalid: {path_pattern}")
                continue
            
            self.logger.info(f"🔧 Debug - pattern_files for {path_pattern}: {len(pattern_files)} files")
            
            # 处理每个文件
            for file_path in pattern_files:
                self.logger.info(f"🔧 Debug - analyzing file: {file_path}")
                try:
                    file_info = self._analyze_file(file_path, default_role)
                    if file_info:
                        all_files.append(file_info)
                        self.logger.info(f"🔧 Debug - file added: {file_path}")
                    else:
                        self.logger.warning(f"🔧 Debug - file rejected: {file_path}")
                except Exception as e:
                    self.logger.warning(f"Failed to analyze file {file_path}: {e}")
        
        self.logger.info(f"Found {len(all_files)} valid sequence files")
        
        # # 如果指定了pattern，进行文件配对
        # if self.config.fastq_pattern:
        #     all_files = self.group_paired_files_by_pattern(all_files)
        #     self.logger.info(f"After pattern grouping: {len(all_files)} file groups")

        # return all_files

        self.logger.info(f"Found {len(all_files)} valid sequence files")
        return all_files
    
    def _scan_directory(self, directory: str) -> List[str]:
        """扫描目录中的序列文件 🗄️"""
        files = []
        
        for root, _, filenames in os.walk(directory):
            for filename in filenames:
                file_path = os.path.join(root, filename)
                if self._is_sequence_file(filename):
                    files.append(file_path)
        
        return files
    
    def _is_sequence_file(self, filename: str) -> bool:
        """判断是否为序列文件 🧬❓"""
        name_lower = filename.lower()
        
        # 检查压缩文件 📦
        for comp_ext in self.compressed_extensions:
            if name_lower.endswith(comp_ext):
                name_lower = name_lower[:-len(comp_ext)]
                break
        
        # 检查序列文件扩展名 🏷️
        for ext in self.fasta_extensions | self.fastq_extensions:
            if name_lower.endswith(ext):
                return True
        
        return False
    
    def _analyze_file(self, file_path: str, default_role: FileRole) -> Optional[FileInfo]:
        """分析单个文件"""
        self.logger.info(f"🔧 Debug - _analyze_file called for: {file_path}")
        
        if not os.path.exists(file_path):
            self.logger.warning(f"🔧 Debug - file does not exist: {file_path}")
            return None
        
        # 获取文件信息
        file_size = os.path.getsize(file_path)
        sample_name = self._extract_sample_name(file_path)
        is_compressed = self._is_compressed(file_path)
        
        self.logger.info(f"🔧 Debug - file_size: {file_size}, sample_name: {sample_name}")
        
        # 检测文件格式
        file_format = self._detect_format(file_path)
        self.logger.info(f"🔧 Debug - detected format: {file_format}")
        
        if file_format == FileFormat.UNKNOWN:
            self.logger.warning(f"Unknown file format: {file_path}")
            return None
    
        # 估算序列和k-mer数量 📈
        estimated_sequences, estimated_kmers = self._estimate_content(
            file_path, file_format, file_size
        )
        
        return FileInfo(
            path=file_path,
            format=file_format,
            role=default_role,
            size_bytes=file_size,
            sample_name=sample_name,
            is_compressed=is_compressed,
            estimated_sequences=estimated_sequences,
            estimated_kmers=estimated_kmers
        )
    
    def _detect_format(self, file_path: str) -> FileFormat:
        """检测文件格式 🕵️"""
        # 首先基于扩展名判断 🏷️
        format_by_ext = self._detect_format_by_extension(file_path)
        if format_by_ext != FileFormat.UNKNOWN:
            # 验证内容 👀
            if self._verify_format_by_content(file_path, format_by_ext):
                return format_by_ext
        
        # 基于内容检测 📖
        return self._detect_format_by_content(file_path)
    
    def _detect_format_by_extension(self, file_path: str) -> FileFormat:
        """基于扩展名检测格式 🏷️"""
        path_lower = file_path.lower()
        
        # 处理压缩文件 📦
        for comp_ext in self.compressed_extensions:
            if path_lower.endswith(comp_ext):
                path_lower = path_lower[:-len(comp_ext)]
                break
        
        # 检查FASTA扩展名
        for ext in self.fasta_extensions:
            if path_lower.endswith(ext):
                return FileFormat.FASTA
        
        # 检查FASTQ扩展名
        for ext in self.fastq_extensions:
            if path_lower.endswith(ext):
                return FileFormat.FASTQ
        
        return FileFormat.UNKNOWN
    
    def _detect_format_by_content(self, file_path: str) -> FileFormat:
        """基于内容检测格式 📖"""
        try:
            opener = gzip.open if self._is_compressed(file_path) else open
            
            with opener(file_path, 'rt', encoding='utf-8') as f:
                # 读取前几行 📑
                lines = []
                for i, line in enumerate(f):
                    if i >= 10:  # 只读前10行
                        break
                    lines.append(line.strip())
                
                if not lines:
                    return FileFormat.UNKNOWN
                
                # FASTQ格式检测 🧬
                if lines[0].startswith('@') and len(lines) >= 4:
                    # 检查FASTQ的4行结构
                    if (lines[2].startswith('+') and 
                        len(lines[1]) == len(lines[3])):  # 序列和质量长度相等
                        return FileFormat.FASTQ
                
                # FASTA格式检测 🧬
                if lines[0].startswith('>'):
                    return FileFormat.FASTA
                
        except Exception as e:
            self.logger.warning(f"Error detecting format for {file_path}: {e} 😟")
        
        return FileFormat.UNKNOWN
    
    def _verify_format_by_content(self, file_path: str, expected_format: FileFormat) -> bool:
        """验证文件格式 ✅"""
        detected_format = self._detect_format_by_content(file_path)
        return detected_format == expected_format
    
    def _is_compressed(self, file_path: str) -> bool:
        """判断文件是否压缩 📦❓"""
        return any(file_path.lower().endswith(ext) for ext in self.compressed_extensions)
    
    def _extract_sample_name(self, file_path: str) -> str:
        """提取样本名称 🏷️"""
        filename = os.path.basename(file_path)
        
        # 移除扩展名 ✂️
        name = filename
        for ext in self.compressed_extensions:
            if name.lower().endswith(ext):
                name = name[:-len(ext)]
                break
        
        for ext in self.fasta_extensions | self.fastq_extensions:
            if name.lower().endswith(ext):
                name = name[:-len(ext)]
                break
        
        # 处理双端测序文件命名 🧬➡️
        # 例如: sample_R1.fastq.gz -> sample
        name = re.sub(r'[._-][Rr]?[12]$', '', name)
        name = re.sub(r'[._-](forward|reverse)$', '', name, flags=re.IGNORECASE)
        name = re.sub(r'[._-](F|R)$', '', name, flags=re.IGNORECASE)
        
        return name or "unknown"
    
    def _estimate_content(self, file_path: str, file_format: FileFormat, file_size: int) -> Tuple[Optional[int], Optional[int]]:
        """估算文件内容 🧮"""
        try:
            # 简单估算：基于文件大小和格式 📏
            if file_format == FileFormat.FASTA:
                # FASTA: 假设平均序列长度5000bp，50%的非序列内容
                avg_seq_length = 5000
                content_ratio = 0.5
                estimated_sequences = int(file_size * content_ratio / avg_seq_length)
                estimated_kmers = int(file_size * content_ratio / self.config.kmer_size * 0.8)
                
            elif file_format == FileFormat.FASTQ:
                # FASTQ: 假设平均read长度150bp，每个read占4行，质量值占50%
                avg_read_length = 150
                lines_per_read = 4
                content_ratio = 0.25  # 序列内容占比
                estimated_sequences = int(file_size * content_ratio / avg_read_length)
                estimated_kmers = int(estimated_sequences * (avg_read_length - self.config.kmer_size + 1))
            
            else:
                return None, None
            
            return estimated_sequences, estimated_kmers
            
        except Exception:
            return None, None
    
    def assign_roles(self, files: List[FileInfo], strategy: AssignmentStrategy) -> Dict[str, List[FileInfo]]:
        """分配文件角色 🏷️🎭"""
        if strategy == AssignmentStrategy.EXPLICIT:
            return self._assign_explicit(files)
        elif strategy == AssignmentStrategy.SIZE_BASED:
            return self._assign_by_size(files)
        elif strategy == AssignmentStrategy.TYPE_BASED:
            return self._assign_by_type(files)
        elif strategy == AssignmentStrategy.INTELLIGENT:
            return self._assign_intelligent(files)
        elif strategy == AssignmentStrategy.INTERACTIVE:
            return self._assign_interactive(files)
        else:
            raise ValueError(f"Unknown assignment strategy: {strategy}")
    
    def _assign_explicit(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """明确分配（基于已设置的角色） 👉"""
        result = {
            'kmer_sources': [],
            'query_targets': []
        }
        
        for file_info in files:
            if file_info.role == FileRole.KMER_SOURCE:
                result['kmer_sources'].append(file_info)
            elif file_info.role == FileRole.QUERY_TARGET:
                result['query_targets'].append(file_info)
        
        return result
    
    def _assign_by_size(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """基于文件大小分配 ⚖️"""
        threshold_bytes = self.config.size_threshold_gb * 1024**3
        
        result = {
            'kmer_sources': [],
            'query_targets': []
        }
        
        for file_info in files:
            if file_info.size_bytes > threshold_bytes:
                result['query_targets'].append(file_info)
            else:
                result['kmer_sources'].append(file_info)
        
        return result
    
    def _assign_by_type(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """基于文件类型分配 📄"""
        result = {
            'kmer_sources': [],
            'query_targets': []
        }
        
        for file_info in files:
            if file_info.format == FileFormat.FASTA:
                result['kmer_sources'].append(file_info)
            elif file_info.format == FileFormat.FASTQ:
                result['query_targets'].append(file_info)
        
        return result
    
    def _assign_intelligent(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """智能分配策略 🧠"""
        # 综合考虑文件类型、大小、数量等因素 🤔
        fasta_files = [f for f in files if f.format == FileFormat.FASTA]
        fastq_files = [f for f in files if f.format == FileFormat.FASTQ]
        
        result = {
            'kmer_sources': [],
            'query_targets': []
        }
        
        # 策略1: 如果FASTA文件较少且较小，优先作为k-mer源 📚
        if len(fasta_files) <= 10 and all(f.size_bytes < self.config.size_threshold_gb * 1024**3 for f in fasta_files):
            result['kmer_sources'].extend(fasta_files)
            result['query_targets'].extend(fastq_files)
        
        # 策略2: 如果FASTQ文件较少，可能是参考样本 🎯
        elif len(fastq_files) <= 5:
            result['kmer_sources'].extend(fastq_files)
            result['query_targets'].extend(fasta_files)
        
        # 策略3: 默认基于文件类型 🗂️
        else:
            result['kmer_sources'].extend(fasta_files)
            result['query_targets'].extend(fastq_files)
        
        # 如果一方为空，重新分配 🔄
        if not result['kmer_sources'] or not result['query_targets']:
            return self._assign_by_size(files)
        
        return result
    
    def _assign_interactive(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """交互式分配 💬"""
        print("\n=== 文件角色分配 | File Role Assignment ===")
        print(f"Found {len(files)} files:")
        
        for i, file_info in enumerate(files):
            size_mb = file_info.size_bytes / 1024**2
            print(f"{(i+1):2d}. {file_info.path}")
            print(f"     Format: {file_info.format.value}, Size: {size_mb:.1f}MB, Sample: {file_info.sample_name}")
        
        print("\nRole assignment options: 🤔")
        print("1. 🤖 Auto-assign by file type (FASTA->source, FASTQ->target)")
        print("2. ⚖️ Auto-assign by file size (small->source, large->target)")
        print("3. ✍️ Manual assignment")
        
        while True:
            try:
                choice = input("\nSelect option (1-3): ").strip()
                if choice == '1':
                    return self._assign_by_type(files)
                elif choice == '2':
                    return self._assign_by_size(files)
                elif choice == '3':
                    return self._manual_assignment(files)
                else:
                    print("❌ Invalid choice, please enter 1, 2, or 3")
            except KeyboardInterrupt:
                print("\nUser cancelled, using intelligent assignment... 🤖")
                return self._assign_intelligent(files)
    
    def _manual_assignment(self, files: List[FileInfo]) -> Dict[str, List[FileInfo]]:
        """手动分配文件角色 ✍️"""
        result = {
            'kmer_sources': [],
            'query_targets': []
        }
        
        print("\nManual assignment:")
        print("Enter 's' for k-mer source (📚), 't' for query target (🎯), 'skip' to skip file (⏭️)")
        
        for i, file_info in enumerate(files):
            while True:
                role = input(f"File {i+1} ({file_info.sample_name}): ").strip().lower()
                if role in ['s', 'source']:
                    result['kmer_sources'].append(file_info)
                    break
                elif role in ['t', 'target']:
                    result['query_targets'].append(file_info)
                    break
                elif role == 'skip':
                    break
                else:
                    print("❌ Invalid input, enter 's', 't', or 'skip'")
        
        return result
    
    # def group_paired_files_by_pattern(self, files: List[FileInfo]) -> List[FileInfo]:
    #     """基于pattern模式合并配对文件"""
    #     if not self.config.fastq_pattern:
    #         return files
        
    #     pattern = self.config.fastq_pattern
    #     if '*' not in pattern:
    #         return files
        
    #     # 只处理FASTQ文件
    #     fastq_files = [f for f in files if f.format == FileFormat.FASTQ]
    #     if not fastq_files:
    #         return files
        
    #     # 找到匹配pattern的文件
    #     matched_files = []
    #     for file_info in fastq_files:
    #         filename = os.path.basename(file_info.path)
    #         if filename.endswith('_1.fq.gz'):  # 直接匹配_1.fq.gz
    #             matched_files.append(file_info)
        
    #     self.logger.info(f"Found {len(matched_files)} files matching _1.fq.gz pattern")
        
    #     # 为每个_1文件找配对的_2文件
    #     merged_files = []
    #     all_files_by_name = {os.path.basename(f.path): f for f in fastq_files}
        
    #     for file_info in matched_files:
    #         filename = os.path.basename(file_info.path)
    #         sample_name = filename.replace('_1.fq.gz', '')  # 提取样本名
            
    #         # 查找配对文件
    #         pair_filename = f"{sample_name}_2.fq.gz"
    #         sample_files = [file_info]  # 包含_1文件
            
    #         if pair_filename in all_files_by_name:
    #             sample_files.append(all_files_by_name[pair_filename])
            
    #         # 创建合并的FileInfo
    #         if len(sample_files) > 1:
    #             all_paths = [f.path for f in sample_files]
    #             total_size = sum(f.size_bytes for f in sample_files)
                
    #             merged_info = FileInfo(
    #                 path=all_paths,
    #                 format=sample_files[0].format,
    #                 role=sample_files[0].role,
    #                 size_bytes=total_size,
    #                 sample_name=sample_name,
    #                 is_compressed=sample_files[0].is_compressed
    #             )
    #             merged_files.append(merged_info)
    #         else:
    #             merged_files.extend(sample_files)
        
    #     self.logger.info(f"Created {len(merged_files)} merged samples")
    #     return merged_files

    def group_paired_files_by_pattern(self, files: List[FileInfo]) -> List[FileInfo]:
        """基于pattern模式合并配对文件"""
        if not self.config.fastq_pattern:
            return files
        
        pattern = self.config.fastq_pattern
        if '*' not in pattern:
            return files
        
        # 只处理FASTQ文件
        fastq_files = [f for f in files if f.format == FileFormat.FASTQ]
        if not fastq_files:
            return files
        
        self.logger.info(f"Processing {len(fastq_files)} FASTQ files with pattern: {pattern}")
        
        # 支持的命名模式
        pattern_mappings = {
            "*_1.fq.gz": "_1.fq.gz",
            "*_2.fq.gz": "_2.fq.gz", 
            "*_1.clean.fq.gz": "_1.clean.fq.gz",
            "*_2.clean.fq.gz": "_2.clean.fq.gz",
            "*_R1.fq.gz": "_R1.fq.gz",
            "*_R2.fq.gz": "_R2.fq.gz",
            "*_R1.fastq.gz": "_R1.fastq.gz",
            "*_R2.fastq.gz": "_R2.fastq.gz",
            "*_1.fastq.gz": "_1.fastq.gz",
            "*_2.fastq.gz": "_2.fastq.gz",
            "*.forward.fq.gz": ".forward.fq.gz",
            "*.reverse.fq.gz": ".reverse.fq.gz",
            "*_forward.fq.gz": "_forward.fq.gz",
            "*_reverse.fq.gz": "_reverse.fq.gz",
            "*_F.fq.gz": "_F.fq.gz",
            "*_R.fq.gz": "_R.fq.gz"
        }
        
        # 获取当前pattern对应的后缀
        suffix = pattern_mappings.get(pattern)
        if not suffix:
            # 如果不在预定义列表中，通用处理
            suffix = pattern.replace("*", "")
        
        # 提取匹配pattern的样本
        sample_groups = {}
        
        for file_info in fastq_files:
            filename = os.path.basename(file_info.path)
            
            if filename.endswith(suffix):
                sample_name = filename.replace(suffix, "")
                if sample_name not in sample_groups:
                    sample_groups[sample_name] = []
                sample_groups[sample_name].append(file_info)
        
        self.logger.info(f"Found {len(sample_groups)} sample groups")
        
        # 为每个样本收集所有文件
        merged_files = []
        all_files_by_name = {os.path.basename(f.path): f for f in fastq_files}
        
        for sample_name, matched_files in sample_groups.items():
            # 查找该样本的所有文件
            sample_files = []
            for filename, file_info in all_files_by_name.items():
                # 检查文件名是否属于该样本
                if (filename.startswith(sample_name + "_") or 
                    filename.startswith(sample_name + ".") or
                    filename == sample_name + suffix):
                    sample_files.append(file_info)
            
            self.logger.info(f"Sample {sample_name}: found {len(sample_files)} files")
            
            if len(sample_files) > 1:
                # 多文件：创建合并的FileInfo
                all_paths = [f.path for f in sample_files]
                total_size = sum(f.size_bytes for f in sample_files)
                
                merged_info = FileInfo(
                    path=all_paths,
                    format=sample_files[0].format,
                    role=sample_files[0].role,
                    size_bytes=total_size,
                    sample_name=sample_name,
                    is_compressed=sample_files[0].is_compressed
                )
                merged_files.append(merged_info)
            else:
                merged_files.extend(sample_files)
        
        self.logger.info(f"Final result: {len(merged_files)} merged samples")
        return merged_files

    def _extract_sample_from_pattern(self, filepath: str, patterns: List[str]) -> Optional[str]:
        """从文件路径中根据pattern提取样本名"""
        filename = os.path.basename(filepath)
        
        for pattern in patterns:
            # 将pattern转换为正则表达式
            # *_1.fq.gz -> (.+)_1\.fq\.gz
            regex_pattern = pattern.replace('*', '(.+)')
            regex_pattern = regex_pattern.replace('.', r'\.')
            regex_pattern = f'^{regex_pattern}$'
            
            self.logger.debug(f"🔧 Testing pattern '{regex_pattern}' against '{filename}'")
            
            try:
                match = re.match(regex_pattern, filename)
                if match:
                    sample_name = match.group(1)
                    self.logger.info(f"🔧 Pattern matched: '{filename}' -> sample: '{sample_name}'")
                    return sample_name
            except re.error as e:
                self.logger.warning(f"🔧 Regex error for pattern '{pattern}': {e}")
                continue
        
        return None