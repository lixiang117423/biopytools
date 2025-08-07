"""
基因组共线性分析配置模块 | Genome Collinearity Analysis Configuration Module
"""

import os
from pathlib import Path
from typing import List, Optional

class CollinearityConfig:
    """共线性分析配置类 | Collinearity Analysis Configuration Class"""
    
    def __init__(self, **kwargs):
        # 基本参数 | Basic parameters
        self.sample_order_file: str = kwargs.get('sample_order_file', '')
        self.output_dir: str = kwargs.get('output_dir', './collinearity_output')
        self.chromosome: Optional[str] = kwargs.get('chromosome', None)
        self.threads: int = kwargs.get('threads', 4)
        
        # 工具路径 | Tool paths
        self.minimap2_path: str = kwargs.get('minimap2_path', 'minimap2')
        self.samtools_path: str = kwargs.get('samtools_path', 'samtools')
        self.syri_path: str = kwargs.get('syri_path', 'syri')
        self.plotsr_path: str = kwargs.get('plotsr_path', 'plotsr')
        
        # 分析参数 | Analysis parameters
        self.minimap2_preset: str = kwargs.get('minimap2_preset', 'asm5')
        self.plotsr_format: str = kwargs.get('plotsr_format', 'png')
        self.figure_width: int = kwargs.get('figure_width', 12)
        self.figure_height: int = kwargs.get('figure_height', 8)
        self.line_width: float = kwargs.get('line_width', 1.5)
        
        # 过滤参数 | Filtering parameters
        self.skip_synteny: bool = kwargs.get('skip_synteny', False)
        self.min_alignment_length: int = kwargs.get('min_alignment_length', 1000)
        
        # 处理路径 | Process paths
        self._process_paths()
        
        # 初始化样本列表 | Initialize sample list
        self.sample_list: List[str] = []
        self.genome_paths: dict = {}
        
        # 加载样本顺序 | Load sample order
        if self.sample_order_file:
            self._load_sample_order()
    
    def _process_paths(self):
        """处理路径为绝对路径 | Process paths to absolute paths"""
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        if self.sample_order_file:
            self.sample_order_file = os.path.normpath(os.path.abspath(self.sample_order_file))
        
        # 创建输出目录 | Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
    
    def _load_sample_order(self):
        """加载样本顺序文件 | Load sample order file"""
        if not os.path.exists(self.sample_order_file):
            raise FileNotFoundError(f"样本顺序文件不存在 | Sample order file not found: {self.sample_order_file}")
        
        print(f"📁 正在读取样本文件 | Reading sample file: {self.sample_order_file}")
        
        with open(self.sample_order_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                original_line = line
                line = line.strip()
                
                print(f"🔍 处理第 {line_num} 行 | Processing line {line_num}: '{original_line.rstrip()}'")
                
                if line and not line.startswith('#'):
                    # 支持制表符或空格分隔 | Support tab or space separation
                    if '\t' in line:
                        parts = line.split('\t')
                        print(f"  📋 使用制表符分隔，得到 {len(parts)} 部分 | Tab-separated, got {len(parts)} parts")
                    else:
                        # 使用空格分隔，但只分割成2部分（路径和名称）
                        # Split by space, but only into 2 parts (path and name)
                        parts = line.split(None, 1)
                        print(f"  📋 使用空格分隔，得到 {len(parts)} 部分 | Space-separated, got {len(parts)} parts")
                    
                    if len(parts) >= 2:
                        genome_path = parts[0].strip()
                        sample_name = parts[1].strip()
                        
                        print(f"  📍 基因组路径 | Genome path: {genome_path}")
                        print(f"  🏷️ 样本名称 | Sample name: {sample_name}")
                        
                        # 检查基因组文件是否存在 | Check if genome file exists
                        if not os.path.exists(genome_path):
                            raise FileNotFoundError(f"❌ 基因组文件不存在 | Genome file not found: {genome_path}")
                        
                        self.sample_list.append(sample_name)
                        self.genome_paths[sample_name] = os.path.abspath(genome_path)
                        print(f"  ✅ 成功添加样本 | Successfully added sample: {sample_name}")
                    else:
                        print(f"  ❌ 跳过格式不正确的行 {line_num} | Skipping incorrectly formatted line {line_num}: {line}")
                else:
                    if line.startswith('#'):
                        print(f"  💬 跳过注释行 | Skipping comment line")
                    else:
                        print(f"  ⏭️ 跳过空行 | Skipping empty line")
        
        print(f"🎯 总共加载了 {len(self.sample_list)} 个样本 | Total loaded {len(self.sample_list)} samples: {self.sample_list}")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查样本顺序文件 | Check sample order file
        if not self.sample_order_file:
            errors.append("❌ 必须提供样本顺序文件 | Sample order file is required")
        elif not os.path.exists(self.sample_order_file):
            errors.append(f"❌ 样本顺序文件不存在 | Sample order file does not exist: {self.sample_order_file}")
        
        # 检查样本数量 | Check sample count
        if len(self.sample_list) < 2:
            errors.append(f"❌ 至少需要2个样本进行比较 | At least 2 samples required for comparison, got: {len(self.sample_list)}")
        
        # 检查线程数 | Check thread count
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        # 检查图像尺寸 | Check figure dimensions
        if self.figure_width <= 0 or self.figure_height <= 0:
            errors.append(f"❌ 图像尺寸必须为正数 | Figure dimensions must be positive: {self.figure_width}x{self.figure_height}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
