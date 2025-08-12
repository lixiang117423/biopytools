"""
⚙️ GTX WGS分析配置管理模块 | GTX WGS Analysis Configuration Management Module ⚙️
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class GTXConfig:
    """📝 GTX分析配置类 | GTX Analysis Configuration Class"""
    
    # 📂 输入输出路径 | Input/Output paths
    input_dir: str
    output_dir: str
    reference: str
    
    # ⚙️ GTX参数 | GTX parameters
    gtx_path: str = "/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx"
    threads: int = 88
    
    # 🔬 质量控制参数 | Quality control parameters
    min_confidence_threshold: int = 30
    min_base_quality: int = 20
    ploidy: int = 2
    pcr_indel_model: str = "CONSERVATIVE"
    
    # 📁 文件模式参数 | File pattern parameters
    read1_pattern: str = "*_1.fq.gz"
    read2_pattern: str = "*_2.fq.gz"
    
    # 🗑️ 临时目录 | Temporary directory
    tmp_dir: Optional[str] = None
    
    # 🔒 内部属性 | Internal attributes
    bam_output_dir: Optional[Path] = None
    vcf_output_dir: Optional[Path] = None
    
    def __post_init__(self):
        """初始化后处理 🚀 | Post-initialization processing"""
        # 标准化路径 🛤️ | Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        
        # 创建输出目录结构 🏗️ | Create output directory structure
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 设置子目录 🗂️ | Set subdirectories
        self.bam_output_dir = self.output_path / "bam"
        self.vcf_output_dir = self.output_path / "vcf"
        self.bam_output_dir.mkdir(parents=True, exist_ok=True)
        self.vcf_output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置临时目录 🗑️ | Set temporary directory
        if self.tmp_dir is None:
            self.tmp_dir = str(self.output_path / "tmp")
        self.tmp_path = Path(self.tmp_dir)
        self.tmp_path.mkdir(parents=True, exist_ok=True)
        
    def validate(self):
        """验证配置参数 ✅ | Validate configuration parameters"""
        errors = []
        
        # 检查输入目录 📂 | Check input directory
        if not os.path.exists(self.input_dir):
            errors.append(f"❌ 输入目录不存在 | Input directory does not exist: {self.input_dir}")
        
        # 检查参考基因组文件 🧬 | Check reference genome file
        if not os.path.exists(self.reference):
            errors.append(f"❌ 参考基因组文件不存在 | Reference genome file does not exist: {self.reference}")
        
        # 检查GTX程序 💻 | Check GTX program
        if not os.path.exists(self.gtx_path):
            errors.append(f"❌ GTX程序不存在 | GTX program does not exist: {self.gtx_path}")
        
        # 检查参数范围 📏 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if not 0 <= self.min_confidence_threshold <= 100:
            errors.append(f"❌ 置信度阈值必须在0-100之间 | Confidence threshold must be between 0-100: {self.min_confidence_threshold}")
        
        if not 0 <= self.min_base_quality <= 100:
            errors.append(f"❌ 碱基质量阈值必须在0-100之间 | Base quality threshold must be between 0-100: {self.min_base_quality}")
        
        if self.ploidy <= 0:
            errors.append(f"❌ 倍性必须为正整数 | Ploidy must be positive: {self.ploidy}")
        
        # 检查文件模式 📁 | Check file patterns
        if "*" not in self.read1_pattern:
            errors.append(f"❌ read1模式必须包含*通配符 | read1 pattern must contain * wildcard: {self.read1_pattern}")
        
        if "*" not in self.read2_pattern:
            errors.append(f"❌ read2模式必须包含*通配符 | read2 pattern must contain * wildcard: {self.read2_pattern}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
