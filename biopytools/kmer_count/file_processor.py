# """
# 📁 文件处理模块 | File Processing Module
# """

# import os
# import re
# import shutil
# import glob
# import gzip
# from pathlib import Path
# from typing import List

# class FileProcessor:
#     """📂 文件处理器 | File Processor"""
    
#     def __init__(self, config, logger, cmd_runner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     # def find_fastq_files(self) -> List[tuple]:
#     #     """🔍 查找FASTQ文件 | Find FASTQ files"""
#     #     self.logger.info(f"🔍 查找FASTQ文件 | Finding FASTQ files")
#     #     self.logger.info(f"📁 输入目录: {self.config.input_dir}")
#     #     self.logger.info(f"🔍 文件模式: {self.config.pattern}")
        
#     #     # 构建完整搜索路径
#     #     input_dir = Path(self.config.input_dir).resolve()
#     #     search_pattern = str(input_dir / self.config.pattern)
        
#     #     self.logger.info(f"🔍 完整搜索路径: {search_pattern}")
#     #     self.logger.info(f"📂 输入目录是否存在: {input_dir.exists()}")
        
#     #     # 列出目录内容用于调试
#     #     if input_dir.exists():
#     #         files_in_dir = list(input_dir.glob("*"))
#     #         self.logger.info(f"📄 目录中的文件: {[f.name for f in files_in_dir[:10]]}")  # 只显示前10个
        
#     #     # 1. 先找到所有匹配模式的R1文件
#     #     r1_files = glob.glob(search_pattern)
#     #     self.logger.info(f"🔍 找到{len(r1_files)}个R1文件")
        
#     #     if not r1_files:
#     #         raise FileNotFoundError(f"📂 未找到匹配的文件 | No files found matching pattern: {search_pattern}")
        
#     #     samples = []
#     #     for r1_file in r1_files:
#     #         # 2. 从R1文件提取样本名（*的部分）
#     #         r1_filename = os.path.basename(r1_file)
            
#     #         # 提取样本名：去掉模式中*后面的部分
#     #         pattern_suffix = self.config.pattern.split('*')[1]  # 例如："_1.clean.fq.gz"
#     #         sample_name = r1_filename.replace(pattern_suffix, "")  # 例如："OV8-105"
            
#     #         self.logger.info(f"📝 提取样本名: {sample_name}")
            
#     #         # 3. 查找对应的R2文件
#     #         r2_pattern_suffix = pattern_suffix.replace('_1.', '_2.').replace('_R1.', '_R2.')
#     #         r2_filename = sample_name + r2_pattern_suffix
#     #         r2_file = input_dir / r2_filename
            
#     #         if r2_file.exists():
#     #             samples.append((sample_name, r1_file, str(r2_file)))
#     #             self.logger.info(f"👥 找到双端样本: {sample_name}")
#     #             self.logger.info(f"  📄 R1: {r1_file}")
#     #             self.logger.info(f"  📄 R2: {r2_file}")
#     #         else:
#     #             samples.append((sample_name, r1_file, None))
#     #             self.logger.info(f"📄 找到单端样本: {sample_name}")
#     #             self.logger.info(f"  📄 R1: {r1_file}")
#     #             self.logger.info(f"  ❌ R2不存在: {r2_file}")
        
#     #     if not samples:
#     #         raise FileNotFoundError("❌ 未找到有效的FASTQ文件 | No valid FASTQ files found")
        
#     #     self.config.samples = samples
#     #     return samples

#     # ln -s软链接识别
#     def find_fastq_files(self) -> List[tuple]:
#         """🔍 查找FASTQ文件 | Find FASTQ files"""
#         self.logger.info(f"🔍 查找FASTQ文件 | Finding FASTQ files")
#         self.logger.info(f"📁 输入目录: {self.config.input_dir}")
#         self.logger.info(f"🔍 文件模式: {self.config.pattern}")
        
#         # 构建完整搜索路径
#         input_dir = Path(self.config.input_dir).resolve()
        
#         self.logger.info(f"📂 输入目录是否存在: {input_dir.exists()}")
        
#         # 列出目录内容用于调试
#         if input_dir.exists():
#             files_in_dir = list(input_dir.glob("*"))
#             self.logger.info(f"📄 目录中的文件: {[f.name for f in files_in_dir[:10]]}")  # 只显示前10个
            
#             # 🔥 检查软链接情况 | Check symlink status
#             symlinks = [f for f in files_in_dir if f.is_symlink()]
#             if symlinks:
#                 self.logger.info(f"🔗 发现{len(symlinks)}个软链接文件 | Found {len(symlinks)} symlink files")
#                 for link in symlinks[:5]:  # 显示前5个
#                     target = link.readlink() if link.is_symlink() else "N/A"
#                     target_exists = link.exists()
#                     self.logger.info(f"  🔗 {link.name} -> {target} (目标存在: {target_exists})")
        
#         # 🔥 使用Path.glob替代glob.glob，更好地处理软链接 | Use Path.glob instead of glob.glob for better symlink handling
#         r1_files = list(input_dir.glob(self.config.pattern))
        
#         # 🔥 过滤和验证文件 | Filter and validate files
#         valid_r1_files = []
#         for file_path in r1_files:
#             if file_path.is_symlink():
#                 # 检查软链接目标是否存在 | Check if symlink target exists
#                 if file_path.exists():
#                     self.logger.info(f"✅ 软链接文件有效 | Valid symlink: {file_path.name}")
#                     valid_r1_files.append(str(file_path))
#                 else:
#                     self.logger.warning(f"⚠️ 软链接目标不存在 | Symlink target missing: {file_path.name} -> {file_path.readlink()}")
#             elif file_path.is_file():
#                 self.logger.info(f"✅ 常规文件 | Regular file: {file_path.name}")
#                 valid_r1_files.append(str(file_path))
#             else:
#                 self.logger.warning(f"⚠️ 无效文件 | Invalid file: {file_path.name}")
        
#         self.logger.info(f"🔍 找到{len(valid_r1_files)}个有效R1文件")
        
#         if not valid_r1_files:
#             raise FileNotFoundError(f"📂 未找到匹配的文件 | No files found matching pattern: {self.config.pattern}")
        
#         samples = []
#         for r1_file in valid_r1_files:
#             # 2. 从R1文件提取样本名（*的部分）
#             r1_filename = os.path.basename(r1_file)
            
#             # 提取样本名：去掉模式中*后面的部分
#             pattern_suffix = self.config.pattern.split('*')[1]  # 例如："_1.clean.fq.gz"
#             sample_name = r1_filename.replace(pattern_suffix, "")  # 例如："OV8-105"
            
#             self.logger.info(f"📝 提取样本名: {sample_name}")
            
#             # 3. 查找对应的R2文件
#             r2_pattern_suffix = pattern_suffix.replace('_1.', '_2.').replace('_R1.', '_R2.')
#             r2_filename = sample_name + r2_pattern_suffix
#             r2_file = input_dir / r2_filename
            
#             # 🔥 改进R2文件检查，支持软链接 | Improved R2 file check with symlink support
#             if r2_file.exists():
#                 if r2_file.is_symlink():
#                     self.logger.info(f"✅ R2软链接文件有效 | Valid R2 symlink: {r2_file.name}")
#                 samples.append((sample_name, r1_file, str(r2_file)))
#                 self.logger.info(f"👥 找到双端样本: {sample_name}")
#                 self.logger.info(f"  📄 R1: {r1_file}")
#                 self.logger.info(f"  📄 R2: {r2_file}")
#             else:
#                 samples.append((sample_name, r1_file, None))
#                 self.logger.info(f"📄 找到单端样本: {sample_name}")
#                 self.logger.info(f"  📄 R1: {r1_file}")
#                 self.logger.info(f"  ❌ R2不存在: {r2_file}")
        
#         if not samples:
#             raise FileNotFoundError("❌ 未找到有效的FASTQ文件 | No valid FASTQ files found")
        
#         self.config.samples = samples
#         return samples
    
#     # def decompress_files(self, sample_files: List[str]) -> List[str]:
#     #     """📦 解压缩文件到临时目录 | Decompress files to temporary directory"""
#     #     decompressed_files = []
        
#     #     for file_path in sample_files:
#     #         if file_path is None:
#     #             decompressed_files.append(None)
#     #             continue
                
#     #         if file_path.endswith('.gz'):
#     #             # 解压缩文件 | Decompress file
#     #             base_name = os.path.basename(file_path).replace('.gz', '')
#     #             temp_file = self.config.temp_dir / base_name
                
#     #             self.logger.info(f"📦 解压缩文件 | Decompressing file: {file_path}")
#     #             with gzip.open(file_path, 'rt') as f_in:
#     #                 with open(temp_file, 'w') as f_out:
#     #                     shutil.copyfileobj(f_in, f_out)
                
#     #             decompressed_files.append(str(temp_file))
#     #         else:
#     #             # 复制文件到临时目录 | Copy file to temporary directory
#     #             base_name = os.path.basename(file_path)
#     #             temp_file = self.config.temp_dir / base_name
#     #             shutil.copy2(file_path, temp_file)
#     #             decompressed_files.append(str(temp_file))
        
#     #     return decompressed_files
#     def decompress_files(self, sample_files: List[str]) -> List[str]:
#         """📦 准备文件列表（不再解压缩）| Prepare file list (no more decompression)"""
#         self.logger.info("📦 准备文件列表 | Preparing file list")
        
#         # 🔥 直接返回原始文件路径，让Jellyfish通过generator处理压缩文件
#         # Return original file paths directly, let Jellyfish handle compressed files via generator
#         processed_files = []
#         for file_path in sample_files:
#             if file_path is None:
#                 processed_files.append(None)
#             else:
#                 if file_path.endswith('.gz'):
#                     self.logger.info(f"📦 压缩文件将通过generator处理 | Compressed file will be processed via generator: {os.path.basename(file_path)}")
#                 else:
#                     self.logger.info(f"📋 普通文件 | Regular file: {os.path.basename(file_path)}")
#                 processed_files.append(file_path)
        
#         return processed_files
    
#     def prepare_kmer_library(self) -> Path:
#         """🧬 准备k-mer库文件 | Prepare k-mer library file"""
#         self.logger.info("🧬 准备k-mer库文件 | Preparing k-mer library file")
        
#         # 读取原始文件并转换为大写 | Read original file and convert to uppercase
#         upper_kmer_file = self.config.temp_dir / "kmers_upper.fasta"
        
#         with open(self.config.kmer_lib, 'r') as f_in:
#             with open(upper_kmer_file, 'w') as f_out:
#                 for line in f_in:
#                     if line.startswith('>'):
#                         f_out.write(line)
#                     else:
#                         f_out.write(line.upper())
        
#         self.logger.info(f"✅ K-mer库文件已准备完成 | K-mer library file prepared: {upper_kmer_file}")
#         return upper_kmer_file
    
#     def find_input_files(self) -> List[tuple]:
#         """查找输入文件（支持FASTQ和FASTA格式）| Find input files (support FASTQ and FASTA formats)"""
#         self.logger.info(f"查找输入文件 | Finding input files")
#         self.logger.info(f"输入目录: {self.config.input_dir}")
#         self.logger.info(f"文件模式: {self.config.pattern}")
        
#         input_dir = Path(self.config.input_dir).resolve()
        
#         # 查找匹配的文件
#         matched_files = list(input_dir.glob(self.config.pattern))
#         valid_files = []
        
#         for file_path in matched_files:
#             if file_path.exists():
#                 valid_files.append(str(file_path))
        
#         self.logger.info(f"找到{len(valid_files)}个输入文件")
        
#         if not valid_files:
#             raise FileNotFoundError(f"未找到匹配的文件 | No files found matching pattern: {self.config.pattern}")
        
#         samples = []
#         for file_path in valid_files:
#             filename = os.path.basename(file_path)
            
#             # 检测文件类型
#             file_type = self._detect_file_type(file_path)
            
#             if file_type == 'fastq':
#                 # FASTQ文件：查找配对的R2文件
#                 sample_name = self._extract_sample_name(filename)
#                 r2_file = self._find_paired_file(file_path, input_dir)
#                 samples.append((sample_name, file_path, r2_file, 'fastq'))
                
#             elif file_type == 'fasta':
#                 # FASTA文件：单个文件作为一个样本
#                 sample_name = self._extract_sample_name_fasta(filename)
#                 samples.append((sample_name, file_path, None, 'fasta'))
                
#             else:
#                 self.logger.warning(f"无法识别文件类型: {filename}")
        
#         self.config.samples = samples
#         return samples

#     def _detect_file_type(self, file_path: str) -> str:
#         """检测文件类型 | Detect file type"""
#         file_path_lower = file_path.lower()
        
#         # FASTA扩展名
#         fasta_extensions = ['.fa', '.fasta', '.fas', '.fna', '.ffn', '.faa']
#         # FASTQ扩展名  
#         fastq_extensions = ['.fq', '.fastq']
        
#         # 检查扩展名（处理.gz压缩文件）
#         for ext in fasta_extensions:
#             if file_path_lower.endswith(ext) or file_path_lower.endswith(ext + '.gz'):
#                 return 'fasta'
        
#         for ext in fastq_extensions:
#             if file_path_lower.endswith(ext) or file_path_lower.endswith(ext + '.gz'):
#                 return 'fastq'
        
#         # 如果扩展名不明确，尝试读取文件内容判断
#         return self._detect_by_content(file_path)

#     def _detect_by_content(self, file_path: str) -> str:
#         """通过文件内容检测格式 | Detect format by file content"""
#         try:
#             if file_path.endswith('.gz'):
#                 import gzip
#                 with gzip.open(file_path, 'rt') as f:
#                     first_line = f.readline().strip()
#             else:
#                 with open(file_path, 'r') as f:
#                     first_line = f.readline().strip()
            
#             if first_line.startswith('>'):
#                 return 'fasta'
#             elif first_line.startswith('@'):
#                 return 'fastq'
#             else:
#                 return 'unknown'
#         except:
#             return 'unknown'

#     def _extract_sample_name_fasta(self, filename: str) -> str:
#         """从FASTA文件名提取样本名 | Extract sample name from FASTA filename"""
#         # 移除路径和扩展名
#         base_name = os.path.splitext(filename)[0]
#         if base_name.endswith('.gz'):
#             base_name = os.path.splitext(base_name)[0]
        
#         return base_name

"""
📁 文件处理模块 | File Processing Module
修改版本：支持FASTA和FASTQ格式 | Modified version: Support FASTA and FASTQ formats
"""

import os
import re
import shutil
import glob
import gzip
from pathlib import Path
from typing import List

class FileProcessor:
    """📂 文件处理器 | File Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def find_input_files(self) -> List[tuple]:
        """🔍 查找输入文件（支持FASTQ和FASTA格式）| Find input files (support FASTQ and FASTA formats)"""
        self.logger.info(f"🔍 查找输入文件 | Finding input files")
        self.logger.info(f"📁 输入目录: {self.config.input_dir}")
        self.logger.info(f"🔍 文件模式: {self.config.pattern}")
        
        # 构建完整搜索路径
        input_dir = Path(self.config.input_dir).resolve()
        
        self.logger.info(f"📂 输入目录是否存在: {input_dir.exists()}")
        
        # 列出目录内容用于调试
        if input_dir.exists():
            files_in_dir = list(input_dir.glob("*"))
            self.logger.info(f"📄 目录中的文件: {[f.name for f in files_in_dir[:10]]}")
        
        # 使用Path.glob查找文件
        matched_files = list(input_dir.glob(self.config.pattern))
        
        # 过滤和验证文件
        valid_files = []
        for file_path in matched_files:
            if file_path.is_symlink():
                if file_path.exists():
                    self.logger.info(f"✅ 软链接文件有效 | Valid symlink: {file_path.name}")
                    valid_files.append(str(file_path))
                else:
                    self.logger.warning(f"⚠️ 软链接目标不存在 | Symlink target missing: {file_path.name}")
            elif file_path.is_file():
                self.logger.info(f"✅ 常规文件 | Regular file: {file_path.name}")
                valid_files.append(str(file_path))
            else:
                self.logger.warning(f"⚠️ 无效文件 | Invalid file: {file_path.name}")
        
        self.logger.info(f"🔍 找到{len(valid_files)}个有效文件")
        
        if not valid_files:
            raise FileNotFoundError(f"📂 未找到匹配的文件 | No files found matching pattern: {self.config.pattern}")
        
        samples = []
        for file_path in valid_files:
            filename = os.path.basename(file_path)
            
            # 🔥 关键修改：检测文件类型
            file_type = self._detect_file_type(file_path)
            self.logger.info(f"📋 检测到文件类型 | Detected file type: {filename} -> {file_type}")
            
            if file_type == 'fastq':
                # FASTQ文件：使用原有的双端逻辑
                sample_name = self._extract_sample_name_fastq(filename)
                r2_file = self._find_paired_fastq_file(file_path, input_dir)
                samples.append((sample_name, file_path, r2_file, 'fastq'))
                
                if r2_file:
                    self.logger.info(f"👥 找到双端FASTQ样本: {sample_name}")
                    self.logger.info(f"  📄 R1: {file_path}")
                    self.logger.info(f"  📄 R2: {r2_file}")
                else:
                    self.logger.info(f"📄 找到单端FASTQ样本: {sample_name}")
                    self.logger.info(f"  📄 文件: {file_path}")
                    
            elif file_type == 'fasta':
                # FASTA文件：单个文件作为一个样本
                sample_name = self._extract_sample_name_fasta(filename)
                samples.append((sample_name, file_path, None, 'fasta'))
                self.logger.info(f"🧬 找到FASTA样本: {sample_name}")
                self.logger.info(f"  📄 文件: {file_path}")
                
            else:
                self.logger.warning(f"⚠️ 无法识别文件类型，跳过 | Unknown file type, skipping: {filename}")
        
        if not samples:
            raise FileNotFoundError("❌ 未找到有效的输入文件 | No valid input files found")
        
        self.config.samples = samples
        return samples
    
    def _detect_file_type(self, file_path: str) -> str:
        """🔍 检测文件类型 | Detect file type"""
        file_path_lower = file_path.lower()
        
        # FASTA扩展名
        fasta_extensions = ['.fa', '.fasta', '.fas', '.fna', '.ffn', '.faa']
        # FASTQ扩展名  
        fastq_extensions = ['.fq', '.fastq']
        
        # 检查扩展名（处理.gz压缩文件）
        for ext in fasta_extensions:
            if file_path_lower.endswith(ext) or file_path_lower.endswith(ext + '.gz'):
                return 'fasta'
        
        for ext in fastq_extensions:
            if file_path_lower.endswith(ext) or file_path_lower.endswith(ext + '.gz'):
                return 'fastq'
        
        # 如果扩展名不明确，尝试读取文件内容判断
        return self._detect_by_content(file_path)
    
    def _detect_by_content(self, file_path: str) -> str:
        """📖 通过文件内容检测格式 | Detect format by file content"""
        try:
            if file_path.endswith('.gz'):
                import gzip
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline().strip()
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline().strip()
            
            if first_line.startswith('>'):
                return 'fasta'
            elif first_line.startswith('@'):
                return 'fastq'
            else:
                return 'unknown'
        except Exception as e:
            self.logger.warning(f"⚠️ 无法读取文件内容 | Cannot read file content: {e}")
            return 'unknown'
    
    def _extract_sample_name_fastq(self, filename: str) -> str:
        """📝 从FASTQ文件名提取样本名 | Extract sample name from FASTQ filename"""
        # 使用原有的FASTQ逻辑
        pattern_suffix = self.config.pattern.split('*')[1] if '*' in self.config.pattern else ''
        sample_name = filename.replace(pattern_suffix, "") if pattern_suffix else filename
        
        # 移除扩展名
        if sample_name.endswith('.gz'):
            sample_name = sample_name[:-3]
        base_name = os.path.splitext(sample_name)[0]
        
        return base_name
    
    def _extract_sample_name_fasta(self, filename: str) -> str:
        """📝 从FASTA文件名提取样本名 | Extract sample name from FASTA filename"""
        # 移除路径和扩展名
        base_name = os.path.splitext(filename)[0]
        if base_name.endswith('.gz'):
            base_name = os.path.splitext(base_name)[0]
        
        return base_name
    
    def _find_paired_fastq_file(self, r1_file: str, input_dir: Path) -> str:
        """🔍 查找配对的FASTQ文件 | Find paired FASTQ file"""
        r1_filename = os.path.basename(r1_file)
        
        # 提取样本名
        pattern_suffix = self.config.pattern.split('*')[1] if '*' in self.config.pattern else ''
        sample_name = r1_filename.replace(pattern_suffix, "") if pattern_suffix else r1_filename
        
        # 查找对应的R2文件
        r2_pattern_suffix = pattern_suffix.replace('_1.', '_2.').replace('_R1.', '_R2.')
        r2_filename = sample_name + r2_pattern_suffix
        r2_file = input_dir / r2_filename
        
        if r2_file.exists():
            return str(r2_file)
        else:
            return None
    
    # 保持原有的 decompress_files 和 prepare_kmer_library 方法不变
    def decompress_files(self, sample_files: List[str]) -> List[str]:
        """📦 准备文件列表（不再解压缩）| Prepare file list (no more decompression)"""
        self.logger.info("📦 准备文件列表 | Preparing file list")
        
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