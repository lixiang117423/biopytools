"""
ANNOVAR注释配置管理模块|ANNOVAR Annotation Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class ANNOVARConfig:
    """ANNOVAR注释配置类|ANNOVAR Annotation Configuration Class"""

    # 必需文件|Required files (按新顺序|in new order)
    vcf_file: str
    gff3_file: str
    genome_file: str
    build_ver: str
    
    # 路径配置|Path configuration
    annovar_path: str = '/share/org/YZWL/yzwl_lixg/software/annovar/annovar'
    database_path: str = './database'
    output_dir: str = './annovar_output'
    
    # 处理参数|Processing parameters
    qual_threshold: int = 20
    skip_gff_fix: bool = False
    skip_gff_cleaning: bool = False  # 新增：跳过GFF3清理步骤|New: Skip GFF3 cleaning step
    skip_vcf_filter: bool = True
    
    # 步骤控制|Step control
    step: Optional[int] = None  # 1-4, None表示运行全部步骤|None means run all steps
    
    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径|Normalize paths
        self.gff3_file = os.path.normpath(os.path.abspath(self.gff3_file))
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.annovar_path = os.path.normpath(os.path.abspath(self.annovar_path))
        self.database_path = os.path.normpath(os.path.abspath(self.database_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 清理build_ver参数|Clean build_ver parameter
        if '/' in self.build_ver or '\\' in self.build_ver:
            self.build_ver = os.path.basename(self.build_ver)
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        
        # 检查必需文件|Check required files
        required_files = [
            ('GFF3文件|GFF3 file', self.gff3_file),
            ('基因组文件|Genome file', self.genome_file),
            ('VCF文件|VCF file', self.vcf_file),
        ]
        
        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")
        
        # 检查ANNOVAR路径|Check ANNOVAR path
        if not os.path.exists(self.annovar_path):
            errors.append(f"ANNOVAR路径不存在|ANNOVAR path does not exist: {self.annovar_path}")
        
        # 检查步骤参数|Check step parameter
        if self.step is not None and self.step not in [1, 2, 3, 4]:
            errors.append(f"无效的步骤编号|Invalid step number: {self.step} (应为1-4|should be 1-4)")
        
        # 数据库路径警告|Database path warning
        if not os.path.isabs(self.database_path) and not os.path.exists(self.database_path):
            # 这只是警告，不是错误|This is just a warning, not an error
            pass
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True