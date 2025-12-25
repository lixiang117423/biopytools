"""
📊 覆盖度分析配置管理模块 | Depth Analysis Configuration Management Module
"""

import os
import re
import glob
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

@dataclass
class DepthConfig:
    """📊 覆盖度分析配置类 | Depth Analysis Configuration Class"""
    
    # 输入文件 | Input files
    input_files: Union[str, List[str]]
    output_file: str = './depth_results.txt'
    
    # 分析参数 | Analysis parameters
    chromosome: str = 'all'  # 染色体名称 | Chromosome name
    region: str = 'all'      # 染色体区间 | Chromosome region (e.g., 100:1235)
    threads: int = 88        # 线程数 | Number of threads
    
    # samtools参数 | samtools parameters
    samtools_path: str = 'samtools'
    quality_threshold: int = 0   # 最小质量值 | Minimum quality score
    mapping_quality: int = 0     # 最小比对质量 | Minimum mapping quality
    
    # 输出格式 | Output format
    output_format: str = 'txt'   # 输出格式 txt | Output format
    compress_output: bool = False # 压缩输出 | Compress output
    
    # 滑窗分析参数 | Sliding window analysis parameters
    enable_window_analysis: bool = False  # 是否启用窗口分析 | Enable window analysis
    window_size: int = 1000      # 窗口大小(bp) | Window size (bp)
    window_step: int = 0         # 步长(bp)，0表示无重叠 | Step size (bp), 0 means no overlap
    
    # 内部属性 | Internal attributes
    base_name: str = 'depth_analysis'
    
    def __post_init__(self):
        """🔧 初始化后处理 | Post-initialization processing"""
        # 处理输入文件列表 | Process input file list
        self.input_files = self._process_input_files(self.input_files)
        
        # 处理输出文件路径 | Process output file path
        self.output_file = self._process_output_file(self.output_file)
        
        # 创建输出目录 | Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 解析区间信息 | Parse region information
        self.start_pos = None
        self.end_pos = None
        if self.region != 'all':
            self._parse_region()
    
    def _process_input_files(self, input_files: Union[str, List[str]]) -> List[str]:
        """🔍 处理输入文件，支持文件和文件夹 | Process input files, support files and directories"""
        if isinstance(input_files, str):
            input_files = [input_files]
        
        processed_files = []
        
        for item in input_files:
            item = os.path.normpath(os.path.abspath(item))
            
            if os.path.isfile(item):
                # 单个文件 | Single file
                if item.lower().endswith(('.bam', '.sam')):
                    processed_files.append(item)
                else:
                    raise ValueError(f"不支持的文件格式，仅支持BAM/SAM文件 | Unsupported file format, only BAM/SAM files are supported: {item}")
            
            elif os.path.isdir(item):
                # 文件夹，自动识别BAM/SAM文件 | Directory, auto-detect BAM/SAM files
                bam_files = glob.glob(os.path.join(item, "*.bam"))
                sam_files = glob.glob(os.path.join(item, "*.sam"))
                
                found_files = bam_files + sam_files
                found_files = [os.path.normpath(os.path.abspath(f)) for f in found_files]
                
                if not found_files:
                    raise ValueError(f"文件夹中未找到BAM/SAM文件 | No BAM/SAM files found in directory: {item}")
                
                processed_files.extend(found_files)
                print(f"📁 在目录 {item} 中找到 {len(found_files)} 个BAM/SAM文件")
            
            else:
                # 可能是通配符模式 | Might be wildcard pattern
                glob_files = glob.glob(item)
                if glob_files:
                    for glob_file in glob_files:
                        if glob_file.lower().endswith(('.bam', '.sam')):
                            processed_files.append(os.path.normpath(os.path.abspath(glob_file)))
                else:
                    raise ValueError(f"文件或目录不存在 | File or directory does not exist: {item}")
        
        if not processed_files:
            raise ValueError("未找到有效的BAM/SAM文件 | No valid BAM/SAM files found")
        
        # 去重并排序 | Remove duplicates and sort
        processed_files = sorted(list(set(processed_files)))
        
        return processed_files
    
    def _process_output_file(self, output_path: str) -> str:
        """📄 处理输出文件路径，支持目录输入 | Process output file path, support directory input"""
        output_path = os.path.normpath(os.path.abspath(output_path))
        
        # 检查是否是目录 | Check if it's a directory
        if os.path.isdir(output_path):
            # 生成默认文件名 | Generate default filename
            if self.chromosome != 'all' and self.region != 'all':
                filename = f"depth_{self.chromosome}_{self.region.replace(':', '-')}.txt"
            elif self.chromosome != 'all':
                filename = f"depth_{self.chromosome}.txt"
            else:
                filename = "depth_results.txt"
            
            output_path = os.path.join(output_path, filename)
            print(f"📄 输出目录检测，自动生成文件名: {filename}")
        
        return output_path
    
    def _parse_region(self):
        """🔍 解析染色体区间 | Parse chromosome region"""
        if ':' in self.region:
            try:
                start, end = self.region.split(':')
                self.start_pos = int(start)
                self.end_pos = int(end)
                if self.start_pos >= self.end_pos:
                    raise ValueError(f"起始位置必须小于结束位置 | Start position must be less than end position: {self.region}")
            except ValueError as e:
                raise ValueError(f"区间格式错误，应为 start:end | Invalid region format, should be start:end: {self.region}")
        else:
            raise ValueError(f"区间格式错误，应为 start:end | Invalid region format, should be start:end: {self.region}")
    
    def validate(self):
        """✅ 验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input files
        for input_file in self.input_files:
            if not os.path.exists(input_file):
                errors.append(f"输入文件不存在 | Input file does not exist: {input_file}")
            elif not input_file.lower().endswith(('.bam', '.sam')):
                errors.append(f"输入文件必须是BAM或SAM格式 | Input file must be BAM or SAM format: {input_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.quality_threshold < 0:
            errors.append(f"质量阈值不能为负数 | Quality threshold cannot be negative: {self.quality_threshold}")
        
        if self.mapping_quality < 0:
            errors.append(f"比对质量不能为负数 | Mapping quality cannot be negative: {self.mapping_quality}")
        
        # 检查窗口分析参数 | Check window analysis parameters
        if self.enable_window_analysis:
            if self.window_size <= 0:
                errors.append(f"窗口大小必须为正整数 | Window size must be positive: {self.window_size}")
            
            if self.window_step < 0:
                errors.append(f"步长不能为负数 | Step size cannot be negative: {self.window_step}")
            
            if self.window_step > 0 and self.window_step > self.window_size:
                errors.append(f"步长不应大于窗口大小 | Step size should not be larger than window size: step={self.window_step}, window={self.window_size}")
        
        # 检查输出格式 | Check output format
        if self.output_format not in ['txt']:
            errors.append(f"不支持的输出格式 | Unsupported output format: {self.output_format}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
