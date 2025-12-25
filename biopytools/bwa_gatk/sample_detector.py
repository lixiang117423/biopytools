# """
# 🔬 样本检测和Read Group管理模块 | Sample Detection and Read Group Management Module
# """

# import os
# import re
# from pathlib import Path
# from typing import List, Dict, Tuple

# class SampleDetector:
#     """样本检测器 | Sample Detector"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
#         self.input_path = Path(config.input_path)
    
#     def detect_samples(self) -> List[Dict]:
#         """检测样本 | Detect samples"""
#         self.logger.info("=" * 80)
#         self.logger.info("🔬 检测样本信息 | Detecting sample information")
#         self.logger.info("=" * 80)
        
#         if self.input_path.is_file():
#             # 单个文件 | Single file
#             samples = self._detect_from_file()
#         else:
#             # 目录 | Directory
#             samples = self._detect_from_directory()
        
#         if not samples:
#             raise ValueError("❌ 未检测到任何样本 | No samples detected")
        
#         self.logger.info(f"✅ 检测到 {len(samples)} 个样本 | Detected {len(samples)} samples")
        
#         # 检测测序类型 | Detect sequencing type
#         seq_type = self._detect_sequencing_type(samples[0])
#         self.logger.info(f"📊 测序类型 | Sequencing type: {seq_type}")
        
#         for sample in samples:
#             sample['seq_type'] = seq_type
#             self.logger.info(f"  📌 {sample['name']}: {sample['r1']}" + 
#                            (f" + {sample['r2']}" if sample['r2'] else ""))
        
#         return samples
    
#     def _detect_from_file(self) -> List[Dict]:
#         """从单个文件检测 | Detect from single file"""
#         filename = self.input_path.name
#         sample_name = self._extract_sample_name(filename)
        
#         # 判断是R1还是R2 | Determine if R1 or R2
#         if '_1.' in filename or '_R1' in filename:
#             # 尝试找到对应的R2 | Try to find corresponding R2
#             r2_file = self._find_paired_file(self.input_path)
#             return [{
#                 'name': sample_name,
#                 'r1': str(self.input_path),
#                 'r2': str(r2_file) if r2_file else None
#             }]
#         else:
#             # 单端或未知 | Single-end or unknown
#             return [{
#                 'name': sample_name,
#                 'r1': str(self.input_path),
#                 'r2': None
#             }]
    
#     def _detect_from_directory(self) -> List[Dict]:
#         """从目录检测样本 | Detect samples from directory"""
#         # 查找所有R1文件 | Find all R1 files
#         r1_pattern = re.compile(r'.*_1\.clean\.fq\.gz$|.*_R1.*\.f(ast)?q\.gz$')
#         r1_files = [f for f in self.input_path.glob('*.fq.gz') if r1_pattern.match(f.name)]
#         r1_files.extend([f for f in self.input_path.glob('*.fastq.gz') if r1_pattern.match(f.name)])
        
#         samples = []
#         for r1_file in sorted(r1_files):
#             sample_name = self._extract_sample_name(r1_file.name)
#             r2_file = self._find_paired_file(r1_file)
            
#             samples.append({
#                 'name': sample_name,
#                 'r1': str(r1_file),
#                 'r2': str(r2_file) if r2_file else None
#             })
        
#         # 如果没有找到R1文件，尝试所有fastq文件 | If no R1 files found, try all fastq files
#         if not samples:
#             all_fastq = list(self.input_path.glob('*.fq.gz')) + list(self.input_path.glob('*.fastq.gz'))
#             for fq_file in sorted(all_fastq):
#                 sample_name = self._extract_sample_name(fq_file.name)
#                 samples.append({
#                     'name': sample_name,
#                     'r1': str(fq_file),
#                     'r2': None
#                 })
        
#         return samples
    
#     def _extract_sample_name(self, filename: str) -> str:
#         """提取样本名 | Extract sample name"""
#         # 移除常见的测序文件后缀 | Remove common sequencing file suffixes
#         patterns = [
#             r'_1\.clean\.fq\.gz$',
#             r'_2\.clean\.fq\.gz$',
#             r'_R[12].*\.f(ast)?q\.gz$',
#             r'\.f(ast)?q\.gz$'
#         ]
        
#         sample_name = filename
#         for pattern in patterns:
#             sample_name = re.sub(pattern, '', sample_name)
        
#         return sample_name
    
#     def _find_paired_file(self, r1_file: Path) -> Path:
#         """查找配对的R2文件 | Find paired R2 file"""
#         r1_name = r1_file.name
        
#         # 尝试不同的R2命名模式 | Try different R2 naming patterns
#         r2_patterns = [
#             r1_name.replace('_1.clean.fq.gz', '_2.clean.fq.gz'),
#             r1_name.replace('_R1', '_R2'),
#             r1_name.replace('_1.', '_2.'),
#         ]
        
#         for r2_name in r2_patterns:
#             r2_file = r1_file.parent / r2_name
#             if r2_file.exists():
#                 return r2_file
        
#         return None
    
#     def _detect_sequencing_type(self, sample: Dict) -> str:
#         """检测测序类型 | Detect sequencing type"""
#         if sample['r2']:
#             return "双端测序 | Paired-end"
#         else:
#             return "单端测序 | Single-end"
    
#     def get_read_group(self, sample_name: str) -> str:
#         """生成Read Group信息 | Generate Read Group information"""
#         # Read Group格式 | Read Group format
#         # @RG\tID:sample_name\tSM:sample_name\tPL:ILLUMINA\tLB:lib1
#         rg = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:lib1"
#         return rg

"""
🔬 样本检测和Read Group管理模块 | Sample Detection and Read Group Management Module
"""

import os
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple # 建议添加 Optional

class SampleDetector:
    """样本检测器 | Sample Detector"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.input_path = Path(config.input_path)
    
    def detect_samples(self) -> List[Dict]:
        """检测样本 | Detect samples"""
        self.logger.info("=" * 80)
        self.logger.info("🔬 检测样本信息 | Detecting sample information")
        self.logger.info("=" * 80)
        
        if self.input_path.is_file():
            samples = self._detect_from_file()
        else:
            samples = self._detect_from_directory()
        
        if not samples:
            raise ValueError("❌ 未检测到任何样本 | No samples detected")
        
        self.logger.info(f"✅ 检测到 {len(samples)} 个样本 | Detected {len(samples)} samples")
        
        seq_type = self._detect_sequencing_type(samples[0])
        self.logger.info(f"📊 测序类型 | Sequencing type: {seq_type}")
        
        for sample in samples:
            sample['seq_type'] = seq_type
            self.logger.info(f"  📌 {sample['name']}: {sample['r1']}" + 
                           (f" + {sample['r2']}" if sample['r2'] else ""))
        
        return samples
    
    def _detect_from_file(self) -> List[Dict]:
        """从单个文件检测 | Detect from single file"""
        filename = self.input_path.name
        
        # 使用更稳健的R1模式判断
        if re.search(r'[_.]R1[_.]|[_.]1[_.]', filename, re.IGNORECASE):
            sample_name = self._extract_sample_name(filename)
            r2_file = self._find_paired_file(self.input_path)
            return [{
                'name': sample_name,
                'r1': str(self.input_path),
                'r2': str(r2_file) if r2_file else None
            }]
        else:
            # 否则，默认为单端或R2文件（不单独处理R2）
            sample_name = self._extract_sample_name(filename)
            return [{
                'name': sample_name,
                'r1': str(self.input_path),
                'r2': None
            }]
            
    def _detect_from_directory(self) -> List[Dict]:
        """从目录检测样本 | Detect samples from directory"""
        # 1. 查找所有可能的fastq文件，无论是否压缩
        all_fastq_files = list(self.input_path.glob('*.fq')) + \
                          list(self.input_path.glob('*.fastq')) + \
                          list(self.input_path.glob('*.fq.gz')) + \
                          list(self.input_path.glob('*.fastq.gz'))
        
        if not all_fastq_files:
            return []

        # 2. 定义一个更通用的R1文件正则表达式，不区分大小写
        r1_pattern = re.compile(r'.*(_R1|_1)[\._].*', re.IGNORECASE)
        
        r1_files = sorted([f for f in all_fastq_files if r1_pattern.match(f.name)])
        
        samples_dict = {} # 使用字典来防止样本重复添加
        
        # 优先处理配对文件
        if r1_files:
            for r1_file in r1_files:
                sample_name = self._extract_sample_name(r1_file.name)
                # 如果样本已存在，则跳过，避免将R2误判为新的R1样本
                if sample_name in samples_dict:
                    continue
                
                r2_file = self._find_paired_file(r1_file)
                samples_dict[sample_name] = {
                    'name': sample_name,
                    'r1': str(r1_file),
                    'r2': str(r2_file) if r2_file else None
                }
        
        # 如果没有找到R1文件，或处理完R1后仍有剩余文件，将它们作为单端处理
        processed_files = set()
        for sample in samples_dict.values():
            processed_files.add(Path(sample['r1']).name)
            if sample['r2']:
                processed_files.add(Path(sample['r2']).name)

        remaining_files = [f for f in all_fastq_files if f.name not in processed_files]
        
        for fq_file in sorted(remaining_files):
            sample_name = self._extract_sample_name(fq_file.name)
            if sample_name not in samples_dict:
                samples_dict[sample_name] = {
                    'name': sample_name,
                    'r1': str(fq_file),
                    'r2': None
                }
                
        return list(samples_dict.values())

    def _extract_sample_name(self, filename: str) -> str:
        """提取样本名 (更通用版本) | Extract sample name (more generic version)"""
        # 移除常见的测序文件后缀，包括 .gz 和非 .gz
        # 使用 re.sub 的一个技巧，逐步去除后缀
        name = filename
        # 首先去除可选的 .gz
        if name.endswith('.gz'):
            name = name[:-3]
        # 然后去除 .fq 或 .fastq
        if name.endswith('.fq'):
            name = name[:-3]
        elif name.endswith('.fastq'):
            name = name[:-6]
            
        # 移除 R1/R2 或 1/2 标识符及其前面的分隔符
        # 例如: sample_R1.clean -> sample
        #       sample.1 -> sample
        name = re.sub(r'[_.]R[12]$', '', name)
        name = re.sub(r'[_.]R[12][_.]', '_', name) # 处理中间的R1/R2
        name = re.sub(r'[_.][12]$', '', name)
        
        # 移除特定的clean等标识
        name = name.replace('.clean', '')
        
        return name

    def _find_paired_file(self, r1_file: Path) -> Optional[Path]:
        """查找配对的R2文件 (更通用版本) | Find paired R2 file (more generic version)"""
        r1_name = r1_file.name
        
        # 生成可能的R2文件名，考虑大小写和不同分隔符
        possible_r2_names = set()
        
        # 规则 1: _1 -> _2
        if '_1' in r1_name:
            possible_r2_names.add(r1_name.replace('_1', '_2', 1))
        # 规则 2: .1 -> .2
        if '.1' in r1_name:
            possible_r2_names.add(r1_name.replace('.1', '.2', 1))
        # 规则 3: _R1 -> _R2 (不区分大小写)
        if re.search('_r1', r1_name, re.IGNORECASE):
            possible_r2_names.add(re.sub('_r1', '_r2', r1_name, count=1, flags=re.IGNORECASE))
        
        for r2_name in possible_r2_names:
            if r2_name == r1_name: continue # 避免替换失败导致死循环
            
            r2_file = r1_file.with_name(r2_name)
            if r2_file.exists():
                return r2_file
        
        return None

    def _detect_sequencing_type(self, sample: Dict) -> str:
        """检测测序类型 | Detect sequencing type"""
        if sample.get('r2'): # 使用 .get 更安全
            return "双端测序 | Paired-end"
        else:
            return "单端测序 | Single-end"
    
    def get_read_group(self, sample_name: str) -> str:
        """生成Read Group信息 | Generate Read Group information"""
        rg = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:lib1"
        return rg