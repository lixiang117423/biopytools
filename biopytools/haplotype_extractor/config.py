"""
VCF单体型提取配置管理模块 | VCF Haplotype Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple

@dataclass
class HaplotypeConfig:
    """单体型提取配置类 | Haplotype Extraction Configuration Class"""
    
    # 输入文件 | Input files
    vcf_file: str
    output_file: str
    position_file: Optional[str] = None
    
    # 单个位点参数 (当不使用position_file时) | Single position parameters
    chromosome: Optional[str] = None
    position: Optional[int] = None
    
    # 工具路径 | Tool paths
    bcftools_path: str = 'bcftools'
    
    # 输出控制 | Output control
    verbose: bool = True
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))
        
        if self.position_file:
            self.position_file = os.path.normpath(os.path.abspath(self.position_file))
        
        # 创建输出目录 | Create output directory
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查VCF文件 | Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        # 检查位置参数 | Check position parameters
        if self.position_file:
            if not os.path.exists(self.position_file):
                errors.append(f"位置文件不存在 | Position file does not exist: {self.position_file}")
        else:
            # 如果没有位置文件，必须提供单个位点信息
            if not self.chromosome or self.position is None:
                errors.append("必须提供位置文件或染色体+位置信息 | Must provide position file or chromosome+position")
            elif self.position <= 0:
                errors.append(f"位点位置必须为正整数 | Position must be positive: {self.position}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def get_positions(self) -> List[Tuple[str, int]]:
        """获取要查询的位置列表 | Get list of positions to query"""
        positions = []
        
        if self.position_file:
            with open(self.position_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
                
                if not lines:
                    raise ValueError(f"位置文件为空 | Position file is empty: {self.position_file}")
                
                # 检测是否有表头：如果第一行第二列不是数字，则认为是表头
                start_idx = 0
                if lines:
                    first_parts = lines[0].split('\t')
                    if len(first_parts) >= 2:
                        try:
                            int(first_parts[1])  # 尝试转换第二列为整数
                            # 如果成功，说明第一行是数据行
                            start_idx = 0
                        except ValueError:
                            # 如果失败，说明第一行是表头
                            start_idx = 1
                
                # 解析数据行
                for line in lines[start_idx:]:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        try:
                            chrom = parts[0]
                            pos = int(parts[1])
                            positions.append((chrom, pos))
                        except ValueError as e:
                            raise ValueError(f"位置文件格式错误 | Invalid position file format: {line}, 错误: {e}")
        else:
            positions.append((self.chromosome, self.position))
        
        return positions
