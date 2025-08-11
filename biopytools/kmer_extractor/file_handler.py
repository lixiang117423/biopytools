"""
📁 文件处理模块 | File Processing Module
"""

import os
import re
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

class FileTypeDetector:
    """🔎 文件类型检测器 | File Type Detector"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def detect_file_type(self, file_path: str) -> str:
        """检测文件类型 | Detect file type"""
        file_path = Path(file_path)
        
        # 基于扩展名检测 | Detect based on extension
        if file_path.suffix.lower() in ['.fq', '.fastq'] or \
           (file_path.suffix.lower() == '.gz' and file_path.stem.endswith(('.fq', '.fastq'))):
            return 'fastq'
        elif file_path.suffix.lower() in ['.fa', '.fasta', '.fas'] or \
             (file_path.suffix.lower() == '.gz' and file_path.stem.endswith(('.fa', '.fasta', '.fas'))):
            return 'fasta'
        
        # 基于文件内容检测 | Detect based on file content
        try:
            if file_path.suffix.lower() == '.gz':
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline().strip()
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline().strip()
            
            if first_line.startswith('@'):
                return 'fastq'
            elif first_line.startswith('>'):
                return 'fasta'
        except Exception as e:
            self.logger.warning(f"⚠️ 无法检测文件内容 | Cannot detect file content: {e}")
        
        raise ValueError(f"❌ 无法确定文件类型 | Cannot determine file type: {file_path}")

class FastqFileMatcher:
    """🔗 FASTQ文件匹配器 | FASTQ File Matcher"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_pattern(self, pattern: str) -> Tuple[str, str, str]:
        """解析FASTQ文件匹配模式 | Parse FASTQ file matching pattern"""
        self.logger.info(f"🔍 解析匹配模式 | Parsing pattern: {pattern}")
        
        if '*' not in pattern:
            raise ValueError(f"❌ 模式必须包含*作为样品名称占位符 | Pattern must contain * as sample name placeholder: {pattern}")
        
        parts = pattern.split('*')
        if len(parts) != 2:
            raise ValueError(f"❌ 模式只能包含一个* | Pattern can only contain one *: {pattern}")
        
        prefix, suffix = parts
        self.logger.info(f"📝 前缀 | Prefix: '{prefix}', 后缀 | Suffix: '{suffix}'")
        
        # 检测配对标识 | Detect pair indicator
        pair_indicators = ['_1', '_2', '_R1', '_R2', '.1', '.2']
        pair_indicator = None
        
        for indicator in pair_indicators:
            if indicator in suffix:
                pair_indicator = indicator
                self.logger.info(f"🔗 检测到配对标识 | Detected pair indicator: {pair_indicator}")
                break
        
        if not pair_indicator:
            self.logger.info(f"📄 未检测到配对标识，将作为单端模式 | No pair indicator detected, treating as single-end mode")
        
        return prefix, suffix, pair_indicator
    
    def match_paired_files(self, file_list: List[str], pattern: str) -> Dict[str, Tuple[str, str]]:
        """匹配配对的FASTQ文件 | Match paired FASTQ files"""
        self.logger.info(f"🔗 使用模式匹配FASTQ文件 | Matching FASTQ files with pattern: {pattern}")
        self.logger.info(f"📄 总文件数 | Total files: {len(file_list)}")
        
        prefix, suffix, pair_indicator = self.parse_pattern(pattern)
        
        if not pair_indicator:
            # 单端测序 | Single-end sequencing
            self.logger.info("📄 检测到单端测序模式 | Detected single-end sequencing mode")
            samples = {}
            for file_path in file_list:
                file_name = os.path.basename(file_path)
                if file_name.startswith(prefix) and file_name.endswith(suffix):
                    sample_name = file_name[len(prefix):-len(suffix)]
                    samples[sample_name] = (file_path, None)
                    self.logger.info(f"✅ 单端样品 | Single-end sample: {sample_name}")
            return samples
        
        # 双端测序 | Paired-end sequencing
        self.logger.info(f"🔗 检测到双端测序模式 | Detected paired-end sequencing mode: {pair_indicator}")
        
        # 确定配对标识 | Determine pair indicators
        if pair_indicator in ['_1', '_R1', '.1']:
            indicator1, indicator2 = pair_indicator, pair_indicator.replace('1', '2')
        else:
            indicator1, indicator2 = pair_indicator.replace('2', '1'), pair_indicator
        
        self.logger.info(f"🔍 配对标识 | Pair indicators: {indicator1} <-> {indicator2}")
        
        # 匹配文件 | Match files
        samples = {}
        files_dict = {}
        
        for file_path in file_list:
            file_name = os.path.basename(file_path)
            
            # 构建匹配模式 | Build matching pattern  
            pattern1 = pattern.replace('*', '(.+)')
            pattern2 = pattern.replace('*', '(.+)').replace(indicator1, indicator2)
            
            # 匹配文件名
            match1 = re.match('^' + pattern1 + '$', file_name)
            match2 = re.match('^' + pattern2 + '$', file_name)
            
            if match1:
                sample_name = match1.group(1)
                files_dict[f"{sample_name}_1"] = file_path
                self.logger.info(f"📍 找到R1文件 | Found R1 file: {sample_name}")
            elif match2:
                sample_name = match2.group(1) 
                files_dict[f"{sample_name}_2"] = file_path
                self.logger.info(f"📍 找到R2文件 | Found R2 file: {sample_name}")
        
        # 组合配对文件 | Combine paired files
        sample_names = set()
        for key in files_dict.keys():
            sample_names.add(key.rsplit('_', 1)[0])
        
        self.logger.info(f"🧬 识别到的样品名称 | Identified sample names: {list(sample_names)}")
        
        for sample_name in sample_names:
            file1_key = f"{sample_name}_1"
            file2_key = f"{sample_name}_2"
            
            if file1_key in files_dict and file2_key in files_dict:
                samples[sample_name] = (files_dict[file1_key], files_dict[file2_key])
                self.logger.info(f"✅ 找到配对文件 | Found paired files for {sample_name}")
            elif file1_key in files_dict:
                samples[sample_name] = (files_dict[file1_key], None)
                self.logger.info(f"⚠️ 只找到R1文件，作为单端处理 | Only found R1 file, treating as single-end: {sample_name}")
        
        self.logger.info(f"📊 最终匹配结果 | Final matching result: {len(samples)} 个样品 | samples")
        return samples

class FastqFileMerger:
    """📄 FASTQ文件合并器 | FASTQ File Merger"""
    
    def __init__(self, logger, output_dir: Path):
        self.logger = logger
        self.output_dir = output_dir
    
    def merge_fastq_files(self, fastq_samples: Dict[str, Tuple[str, str]]) -> Tuple[str, str]:
        """合并FASTQ文件 | Merge FASTQ files"""
        self.logger.info("🔗 开始合并FASTQ文件 | Starting FASTQ file merging")
        
        # 收集所有R1和R2文件
        r1_files = []
        r2_files = []
        
        for sample_name, (file1, file2) in fastq_samples.items():
            if file1:
                r1_files.append(file1)
                self.logger.info(f"📄 R1文件: {sample_name} -> {os.path.basename(file1)}")
            if file2:
                r2_files.append(file2)
                self.logger.info(f"📄 R2文件: {sample_name} -> {os.path.basename(file2)}")
        
        # 合并R1文件
        merged_r1 = str(self.output_dir / "merged_R1.fq")
        self._merge_files(r1_files, merged_r1, "R1")
        
        # 合并R2文件（如果存在）
        merged_r2 = None
        if r2_files:
            merged_r2 = str(self.output_dir / "merged_R2.fq")
            self._merge_files(r2_files, merged_r2, "R2")
        
        return merged_r1, merged_r2
    
    def _merge_files(self, input_files: List[str], output_file: str, read_type: str):
        """合并文件列表 | Merge file list"""
        self.logger.info(f"🔗 合并{read_type}文件: {len(input_files)}个文件 -> {os.path.basename(output_file)}")
        
        with open(output_file, 'w') as outf:
            for i, input_file in enumerate(input_files, 1):
                self.logger.info(f"  📄 处理文件 {i}/{len(input_files)}: {os.path.basename(input_file)}")
                
                # 处理压缩和非压缩文件
                if input_file.endswith('.gz'):
                    with gzip.open(input_file, 'rt') as inf:
                        shutil.copyfileobj(inf, outf)
                else:
                    with open(input_file, 'r') as inf:
                        shutil.copyfileobj(inf, outf)
        
        self.logger.info(f"✅ {read_type}文件合并完成 | {read_type} file merging completed: {output_file}")
