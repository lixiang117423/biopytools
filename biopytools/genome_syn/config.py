"""
配置管理模块 | Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any

@dataclass 
class GenomeSynConfig:
    """基因组共线性分析配置类 | Genome Synteny Analysis Configuration Class"""
    
    # 必需参数 | Required parameters
    sample_map: Optional[str] = None
    config_file: Optional[str] = None
    output_dir: str = "./genome_syn_output"
    
    # 比对参数 | Alignment parameters
    aligner: str = "minimap2"
    alignment_mode: str = "chain"  # chain, star, all_vs_all
    threads: int = 16
    min_length: int = 5000
    
    # 可视化参数 | Visualization parameters  
    canvas_width: Optional[int] = None
    canvas_height: Optional[int] = None
    output_formats: List[str] = None
    
    # 高级参数 | Advanced parameters
    generate_config: bool = False
    aligner_params: Dict[str, Any] = None
    
    # 染色体过滤参数 | Chromosome filtering parameters
    chromosomes: Optional[List[int]] = None
    _chromosome_str: Optional[str] = None
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 设置默认输出格式 | Set default output formats
        if self.output_formats is None:
            self.output_formats = ["svg", "png"]
        
        # 设置默认比对器参数 | Set default aligner parameters
        if self.aligner_params is None:
            self.aligner_params = self._get_default_aligner_params()
        
        # 解析染色体参数 | Parse chromosome parameters
        if self._chromosome_str:
            self.chromosomes = self._parse_chromosome_string(self._chromosome_str)
    
    def _parse_chromosome_string(self, chr_str: str) -> List[int]:
        """解析染色体字符串 | Parse chromosome string
        
        支持格式:
        - "1,2,3" : 指定第1、2、3号染色体
        - "1-5" : 指定第1到5号染色体
        - "1,3-5,7" : 混合格式
        """
        chromosomes = []
        for part in chr_str.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    chromosomes.extend(range(start, end + 1))
                except ValueError:
                    raise ValueError(f"无效的染色体范围格式: {part}")
            else:
                try:
                    chromosomes.append(int(part))
                except ValueError:
                    raise ValueError(f"无效的染色体编号: {part}")
        
        # 去重、排序并验证
        chromosomes = sorted(set(chromosomes))
        if any(c <= 0 for c in chromosomes):
            raise ValueError("染色体编号必须为正整数")
        
        return chromosomes
    
    def _get_default_aligner_params(self) -> Dict[str, Any]:
        """获取默认比对器参数 | Get default aligner parameters"""
        default_params = {
            "minimap2": {
                "preset": "asm5",
                "min_length": self.min_length
            },
            "mummer": {
                "match_type": "mumreference",  # mum, mumreference, maxmatch
                "min_match": 20,              # minimum match length
                "min_cluster": 65,            # minimum cluster length
                "max_gap": 90,                # maximum gap between matches
                "prefix": "mummer_out"        # output prefix
            },
            "mcscanx": {
                "evalue": "1e-5",
                "min_match": 5
            },
            "syri": {
                "min_sv_size": 50,
                "min_alignment_length": 1000
            }
        }
        return default_params.get(self.aligner, {})
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input files
        if not self.sample_map and not self.config_file:
            errors.append("必须提供sample_map或config_file")
        
        if self.sample_map and not os.path.exists(self.sample_map):
            errors.append(f"Sample map文件不存在: {self.sample_map}")
        
        if self.config_file and not os.path.exists(self.config_file):
            errors.append(f"配置文件不存在: {self.config_file}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append("线程数必须为正整数")
        
        if self.min_length <= 0:
            errors.append("最小长度必须为正整数")
        
        # 检查比对器支持 | Check aligner support
        supported_aligners = ["minimap2", "mummer", "mcscanx", "syri"]
        if self.aligner not in supported_aligners:
            errors.append(f"不支持的比对器: {self.aligner}，支持的比对器: {', '.join(supported_aligners)}")
        
        # 检查染色体参数 | Check chromosome parameters
        if self.chromosomes and len(self.chromosomes) == 0:
            errors.append("染色体列表不能为空")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def get_chromosome_suffix(self) -> str:
        """获取染色体后缀用于文件命名 | Get chromosome suffix for file naming"""
        if self.chromosomes:
            return "_chr" + "_".join(map(str, self.chromosomes))
        return ""