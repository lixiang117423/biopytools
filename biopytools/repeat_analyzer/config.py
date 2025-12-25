"""
重复序列分析配置管理模块 | Repeat Sequence Analysis Configuration Management Module
"""

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class RepeatConfig:
    """重复序列分析配置类 | Repeat Sequence Analysis Configuration Class"""
    
    # 🧬 输入文件 | Input files
    genome_file: str
    output_dir: str = './repeat_output'
    
    # ⚙️ 分析参数 | Analysis parameters
    threads: int = 88
    skip_modeler: bool = False  # 跳过RepeatModeler步骤 | Skip RepeatModeler step
    skip_ltr: bool = False      # 跳过LTR分析步骤 | Skip LTR analysis step
    
    # 🔧 工具路径 | Tool paths - 自动检测或手动指定
    repeatmodeler_path: str = 'RepeatModeler'
    ltr_finder_path: str = 'ltr_finder'
    ltrharvest_path: str = 'gt ltrharvest'
    ltr_retriever_path: str = 'LTR_retriever'
    repeatmasker_path: str = 'RepeatMasker'
    tesorter_path: str = 'TEsorter'
    
    # 📁 内部属性 | Internal attributes
    base_name: str = 'repeat_analysis'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 从基因组文件名生成base_name | Generate base_name from genome filename
        genome_name = Path(self.genome_file).stem
        if genome_name.endswith('.fasta') or genome_name.endswith('.fa'):
            genome_name = Path(genome_name).stem
        self.base_name = f"{genome_name}_repeat"
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genome file
        if not os.path.exists(self.genome_file):
            errors.append(f"❌ 基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 检查线程数 | Check thread count
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    # def detect_tool_paths(self):
    #     """自动检测工具路径 | Auto-detect tool paths"""
    #     tools = {
    #         'repeatmodeler_path': ['RepeatModeler', '~/.local/bin/RepeatModeler'],
    #         'ltr_finder_path': ['ltr_finder', '~/.local/bin/ltr_finder'],
    #         'ltrharvest_path': ['gt ltrharvest', '~/.local/bin/gt ltrharvest'],
    #         'ltr_retriever_path': ['LTR_retriever', '~/.local/bin/LTR_retriever'],
    #         'repeatmasker_path': ['RepeatMasker', '~/.local/bin/RepeatMasker'],
    #         'tesorter_path': ['TEsorter', '~/.local/bin/TEsorter']
    #     }
        
    #     detected_tools = {}
        
    #     for attr_name, possible_paths in tools.items():
    #         for tool_path in possible_paths:
    #             expanded_path = os.path.expanduser(tool_path)
    #             if shutil.which(expanded_path.split()[0]):  # 检查第一个命令
    #                 detected_tools[attr_name] = expanded_path
    #                 break
    #         else:
    #             # 如果没找到，保持原有值 | Keep original value if not found
    #             detected_tools[attr_name] = getattr(self, attr_name)
        
    #     # 更新配置 | Update configuration
    #     for attr_name, tool_path in detected_tools.items():
    #         setattr(self, attr_name, tool_path)
        
    #     return detected_tools

    def detect_tool_paths(self):
        """自动检测工具路径 | Auto-detect tool paths"""
        # 通用路径列表，而不是硬编码特定用户路径
        common_paths = [
            '',  # PATH中的命令
            '~/.local/bin/',
            '/usr/local/bin/',
            '/opt/local/bin/',
            '~/bin/',
        ]
        
        tools = {
            'repeatmodeler_path': ['RepeatModeler'],
            'ltr_finder_path': ['ltr_finder'],
            'ltrharvest_path': ['gt ltrharvest'],  # 复合命令
            'ltr_retriever_path': ['LTR_retriever'],
            'repeatmasker_path': ['RepeatMasker'],
            'tesorter_path': ['TEsorter']
        }
        
        detected_tools = {}
        
        for attr_name, tool_names in tools.items():
            found = False
            for tool_name in tool_names:
                # 先检查PATH
                if shutil.which(tool_name.split()[0]):
                    detected_tools[attr_name] = tool_name
                    found = True
                    break
                
                # 然后检查常见路径
                for path_prefix in common_paths[1:]:  # 跳过空字符串
                    full_path = os.path.expanduser(path_prefix + tool_name)
                    if os.path.exists(full_path) or shutil.which(full_path.split()[0]):
                        detected_tools[attr_name] = full_path
                        found = True
                        break
                
                if found:
                    break
            
            if not found:
                detected_tools[attr_name] = getattr(self, attr_name)  # 保持原值
