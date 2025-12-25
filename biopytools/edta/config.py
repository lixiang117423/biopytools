"""
🌾 EDTA配置管理模块 | EDTA Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class EDTAConfig:
    """EDTA分析配置类 | EDTA Analysis Configuration Class"""
    
    # 🔹 必需参数 | Required parameters
    genome: str
    output_dir: str = './edta_output'
    
    # 🔹 基本参数 | Basic parameters (与EDTA默认值一致)
    species: str = "others"  # Rice|Maize|others
    step: str = "all"  # all|filter|final|anno
    overwrite: int = 0  # 0|1
    cds: Optional[str] = None
    curatedlib: Optional[str] = None
    rmlib: Optional[str] = None
    sensitive: int = 0  # 0|1
    anno: int = 1  # 0|1 (默认开启注释)
    rmout: Optional[str] = None
    maxdiv: int = 40  # 0-100
    evaluate: int = 1  # 0|1 (默认开启评估)
    exclude: Optional[str] = None
    force: int = 1  # 0|1
    u: float = 1.3e-8  # neutral mutation rate
    
    # 🔹 依赖路径 | Dependency paths
    repeatmodeler: Optional[str] = None
    repeatmasker: Optional[str] = None
    annosine: Optional[str] = None
    ltrretriever: Optional[str] = None
    
    # 🔹 系统参数 | System parameters
    threads: int = 88  # 默认线程数88
    debug: int = 0  # 0|1
    check_dependencies: bool = False
    
    # 🔹 批量处理参数 | Batch processing parameters
    genome_list: Optional[List[str]] = None
    batch_mode: bool = False
    
    # 🔹 断点续跑参数 | Resume parameters
    resume: bool = False
    resume_step: Optional[str] = None
    
    # 🔹 后处理参数 | Post-processing parameters
    generate_stats: bool = True
    generate_plots: bool = True
    compare_results: bool = False
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        if isinstance(self.genome, str):
            self.genome = os.path.normpath(os.path.abspath(self.genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 处理可选文件路径
        optional_files = [self.cds, self.curatedlib, self.rmlib, self.rmout, self.exclude]
        for i, file_path in enumerate(optional_files):
            if file_path:
                optional_files[i] = os.path.normpath(os.path.abspath(file_path))
        self.cds, self.curatedlib, self.rmlib, self.rmout, self.exclude = optional_files
        
        # 处理批量模式 | Handle batch mode
        if self.genome_list:
            self.batch_mode = True
            self.genome_list = [os.path.normpath(os.path.abspath(g)) for g in self.genome_list]

        # if self.force == 0:
        #     self.force = 1
        #     print("检测到EDTA文件检查问题，自动启用--force 1参数")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查基因组文件 | Check genome file
        if self.batch_mode:
            for genome in self.genome_list:
                if not os.path.exists(genome):
                    errors.append(f"基因组文件不存在 | Genome file does not exist: {genome}")
        else:
            if not os.path.exists(self.genome):
                errors.append(f"基因组文件不存在 | Genome file does not exist: {self.genome}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if not 0 <= self.maxdiv <= 100:
            errors.append(f"最大分化度必须在0-100之间 | Maximum divergence must be between 0-100: {self.maxdiv}")
        
        if self.species not in ["Rice", "Maize", "others"]:
            errors.append(f"物种参数必须是Rice、Maize或others | Species must be Rice, Maize or others: {self.species}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
    
    def to_edta_command(self, genome_file: str = None) -> List[str]:
        """转换为EDTA命令参数 | Convert to EDTA command parameters"""
        cmd = ["EDTA.pl"]
        
        # 使用指定的基因组文件或默认基因组文件
        target_genome = genome_file if genome_file else self.genome
        cmd.extend(["--genome", target_genome])
        
        # 添加基本参数 | Add basic parameters
        cmd.extend(["--species", self.species])
        cmd.extend(["--step", self.step])
        cmd.extend(["--overwrite", str(self.overwrite)])
        cmd.extend(["--sensitive", str(self.sensitive)])
        cmd.extend(["--anno", str(self.anno)])
        cmd.extend(["--maxdiv", str(self.maxdiv)])
        cmd.extend(["--evaluate", str(self.evaluate)])
        cmd.extend(["--force", str(self.force)])
        cmd.extend(["--u", str(self.u)])
        cmd.extend(["--threads", str(self.threads)])
        cmd.extend(["--debug", str(self.debug)])
        
        # 添加可选文件参数 | Add optional file parameters
        optional_params = [
            ("--cds", self.cds),
            ("--curatedlib", self.curatedlib),
            ("--rmlib", self.rmlib),
            ("--rmout", self.rmout),
            ("--exclude", self.exclude)
        ]
        
        for param, value in optional_params:
            if value:
                cmd.extend([param, value])
        
        # 添加依赖路径参数 | Add dependency path parameters
        dependency_params = [
            ("--repeatmodeler", self.repeatmodeler),
            ("--repeatmasker", self.repeatmasker),
            ("--annosine", self.annosine),
            ("--ltrretriever", self.ltrretriever)
        ]
        
        for param, value in dependency_params:
            if value:
                cmd.extend([param, value])
        
        return cmd
