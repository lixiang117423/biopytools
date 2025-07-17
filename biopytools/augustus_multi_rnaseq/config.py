"""
Augustus多转录组预测配置管理模块 | Augustus Multiple RNA-seq Prediction Configuration Module
"""

import os
import glob
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional

@dataclass
class AugustusConfig:
    """Augustus多转录组预测配置类 | Augustus Multiple RNA-seq Prediction Configuration Class"""
    
    # 必需文件 | Required files
    genome_file: str
    species_model: str  # Augustus训练的模型名 | Augustus trained model name
    
    # 输入方式 (二选一) | Input method (choose one)
    input_dir: Optional[str] = None  # FASTQ文件目录 | FASTQ files directory
    config_file: Optional[str] = None  # 样本配置文件 | Sample configuration file
    
    # 文件模式参数 | File pattern parameters
    pattern: str = "*.R1.fastq.gz"  # R1文件匹配模式 | R1 file pattern
    
    output_dir: str = "augustus_multi_rnaseq"
    
    # 处理参数 | Processing parameters
    threads: int = 8
    hisat2_index: Optional[str] = None  # 如果不提供，会自动生成 | Auto-generate if not provided
    
    # Augustus参数 | Augustus parameters
    alternatives_from_evidence: bool = True
    allow_hinted_splicesites: str = "atac"
    gff3_output: bool = True
    
    # 过滤参数 | Filtering parameters
    min_intron_support: int = 2
    filter_bam: bool = True
    
    # 内部属性 | Internal attributes
    samples: List[Dict] = None
    
    # def __post_init__(self):
    #     """初始化后处理 | Post-initialization processing"""
    #     self.output_path = Path(self.output_dir)
    #     self.output_path.mkdir(parents=True, exist_ok=True)
        
    #     # 标准化路径 | Normalize paths
    #     self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
    #     self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
    #     # 只在路径不为None时处理 | Only process paths when not None
    #     if self.input_dir is not None:
    #         self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        
    #     if self.config_file is not None:
    #         self.config_file = os.path.normpath(os.path.abspath(self.config_file))
        
    #     # 设置默认的HISAT2索引路径 | Set default HISAT2 index path
    #     if self.hisat2_index is None:
    #         genome_name = os.path.splitext(os.path.basename(self.genome_file))[0]
    #         self.hisat2_index = os.path.join(os.path.dirname(self.genome_file), f"{genome_name}_hisat2_index")

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths (先标准化再创建Path对象)
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 使用绝对路径创建Path对象 | Create Path object with absolute path
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 只在路径不为None时处理 | Only process paths when not None
        if self.input_dir is not None:
            self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        
        if self.config_file is not None:
            self.config_file = os.path.normpath(os.path.abspath(self.config_file))
        
        # 设置默认的HISAT2索引路径 | Set default HISAT2 index path
        if self.hisat2_index is None:
            genome_name = os.path.splitext(os.path.basename(self.genome_file))[0]
            self.hisat2_index = os.path.join(os.path.dirname(self.genome_file), f"{genome_name}_hisat2_index")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查输入方式 | Check input method
        if not self.input_dir and not self.config_file:
            errors.append("必须指定输入目录(-i/--input-dir)或配置文件(-c/--config)之一 | Must specify either input directory (-i/--input-dir) or config file (-c/--config)")
        
        if self.input_dir and self.config_file:
            errors.append("不能同时指定输入目录和配置文件，请选择其中一种方式 | Cannot specify both input directory and config file, please choose one method")
        
        # 检查输入目录或配置文件是否存在 | Check if input directory or config file exists
        if self.input_dir and not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
        if self.config_file and not os.path.exists(self.config_file):
            errors.append(f"配置文件不存在 | Config file does not exist: {self.config_file}")
        
        # 检查线程数 | Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")
        
        # 检查模型名称 | Check model name
        if not self.species_model:
            errors.append("必须指定Augustus物种模型 | Must specify Augustus species model")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def parse_sample_config(self):
        """解析样本配置 | Parse sample configuration"""
        if self.input_dir:
            # 自动发现样本 | Auto-discover samples
            return self._auto_discover_samples()
        elif self.config_file:
            # 从配置文件读取 | Read from config file
            return self._parse_config_file()
        else:
            raise ValueError("必须指定输入目录或配置文件 | Must specify input directory or config file")
    
    def _auto_discover_samples(self):
        """自动发现样本 | Auto-discover samples"""
        samples = []
        
        try:
            # 获取所有R1文件 | Get all R1 files
            r1_pattern = os.path.join(self.input_dir, self.pattern)
            r1_files = glob.glob(r1_pattern)
            r1_files.sort()
            
            if not r1_files:
                raise ValueError(f"在目录 {self.input_dir} 中未找到匹配模式 {self.pattern} 的文件 | No files matching pattern {self.pattern} found in directory {self.input_dir}")
            
            # 检测R2文件模式 | Detect R2 file pattern
            r2_pattern = self._get_r2_pattern(self.pattern)
            
            for r1_file in r1_files:
                # 提取样本名称 | Extract sample name
                sample_name = self._extract_sample_name(r1_file, self.pattern)
                
                # 构建R2文件路径 | Build R2 file path
                r2_file = self._get_r2_file_path(r1_file, r2_pattern)
                
                # 检查R2文件是否存在 | Check if R2 file exists
                if not os.path.exists(r2_file):
                    raise ValueError(f"找不到样本 {sample_name} 对应的R2文件 | Cannot find R2 file for sample {sample_name}: {r2_file}")
                
                samples.append({
                    'name': sample_name,
                    'r1_file': os.path.abspath(r1_file),
                    'r2_file': os.path.abspath(r2_file)
                })
        
        except Exception as e:
            raise ValueError(f"自动发现样本失败 | Failed to auto-discover samples: {e}")
        
        if not samples:
            raise ValueError("未找到有效的样本 | No valid samples found")
        
        self.samples = samples
        return samples
    
    def _parse_config_file(self):
        """从配置文件解析样本 | Parse samples from config file"""
        samples = []
        
        try:
            with open(self.config_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    # 跳过空行和注释行 | Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 3:
                        sample_name = parts[0]
                        r1_file = parts[1]
                        r2_file = parts[2]
                        
                        # 检查文件是否存在 | Check if files exist
                        if not os.path.exists(r1_file):
                            raise ValueError(f"R1文件不存在 | R1 file does not exist: {r1_file}")
                        if not os.path.exists(r2_file):
                            raise ValueError(f"R2文件不存在 | R2 file does not exist: {r2_file}")
                        
                        samples.append({
                            'name': sample_name,
                            'r1_file': os.path.abspath(r1_file),
                            'r2_file': os.path.abspath(r2_file)
                        })
                    else:
                        raise ValueError(f"第{line_num}行格式错误 | Line {line_num} format error: {line}")
        
        except Exception as e:
            raise ValueError(f"解析样本配置文件失败 | Failed to parse sample config file: {e}")
        
        if not samples:
            raise ValueError("未找到有效的样本配置 | No valid sample configuration found")
        
        self.samples = samples
        return samples
    
    def _get_r2_pattern(self, r1_pattern):
        """根据R1模式推断R2模式 | Infer R2 pattern from R1 pattern"""
        # 常见的R1/R2模式映射 | Common R1/R2 pattern mappings
        r1_r2_mappings = [
            ('R1', 'R2'),
            ('_1', '_2'),
            ('.1', '.2'),
            ('read1', 'read2'),
            ('_r1', '_r2'),
        ]
        
        for r1_marker, r2_marker in r1_r2_mappings:
            if r1_marker in r1_pattern:
                return r1_pattern.replace(r1_marker, r2_marker)
        
        # 如果没有找到明确的R1标记，尝试智能推断
        raise ValueError(f"无法从R1模式推断R2模式 | Cannot infer R2 pattern from R1 pattern: {r1_pattern}")
    
    def _extract_sample_name(self, r1_file, pattern):
        """从R1文件名提取样本名称 | Extract sample name from R1 filename"""
        basename = os.path.basename(r1_file)
        
        # 移除路径中的通配符部分 | Remove wildcard parts from pattern
        pattern_basename = os.path.basename(pattern)
        
        # 将模式转换为正则表达式 | Convert pattern to regex
        regex_pattern = pattern_basename.replace('*', '(.+)')
        regex_pattern = regex_pattern.replace('.', r'\.')
        
        match = re.match(regex_pattern, basename)
        if match:
            return match.group(1)
        
        # 如果正则匹配失败，使用简单的字符串替换 | Fallback to simple string replacement
        return basename.replace(pattern_basename.replace('*', ''), '')
    
    def _get_r2_file_path(self, r1_file, r2_pattern):
        """根据R1文件路径和R2模式生成R2文件路径 | Generate R2 file path from R1 file and R2 pattern"""
        r1_basename = os.path.basename(r1_file)
        r1_pattern_basename = os.path.basename(self.pattern)
        r2_pattern_basename = os.path.basename(r2_pattern)
        
        # 提取样本名称 | Extract sample name
        sample_name = self._extract_sample_name(r1_file, self.pattern)
        
        # 构建R2文件名 | Build R2 filename
        r2_basename = r2_pattern_basename.replace('*', sample_name)
        
        # 构建完整的R2路径 | Build full R2 path
        r2_file = os.path.join(os.path.dirname(r1_file), r2_basename)
        
        return r2_file
