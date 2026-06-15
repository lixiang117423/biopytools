"""
EDTA转座子注释配置模块|EDTA TE Annotation Configuration Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class EDTAConfig:
    """EDTA配置类|EDTA Configuration Class"""

    # 必需参数|Required parameters
    genome: str

    # 可选参数|Optional parameters
    species: str = "others"
    step: str = "all"
    overwrite: int = 0
    cds: Optional[str] = None
    curatedlib: Optional[str] = None
    rmlib: Optional[str] = None
    sensitive: int = 0
    anno: int = 0
    rmout: Optional[str] = None
    maxdiv: int = 40
    evaluate: int = 0
    exclude: Optional[str] = None
    force: int = 0
    u: float = 1.3e-8
    threads: int = 12
    debug: int = 0
    output_dir: str = "./edta_output"

    # EDTA路径配置|EDTA path configuration
    edta_path: Optional[str] = None
    repeatmasker_path: Optional[str] = None
    repeatmodeler_path: Optional[str] = None
    ltrretriever_path: Optional[str] = None
    annosine_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认EDTA路径|Set default EDTA path
        if self.edta_path is None:
            conda_env = os.environ.get('CONDA_DEFAULT_ENV', '')
            if conda_env == 'EDTA_v.2.2.2' or 'EDTA' in conda_env:
                # 尝试从conda环境查找|Try to find from conda environment
                conda_prefix = os.environ.get('CONDA_PREFIX', '')
                if conda_prefix:
                    self.edta_path = os.path.join(conda_prefix, 'share', 'EDTA')
                else:
                    # 使用已知的安装路径|Use known installation path
                    self.edta_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA'
            else:
                self.edta_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA'

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        self.genome = os.path.abspath(self.genome) if self.genome else None
        if self.cds:
            self.cds = os.path.abspath(self.cds)
        if self.curatedlib:
            self.curatedlib = os.path.abspath(self.curatedlib)
        if self.rmlib:
            self.rmlib = os.path.abspath(self.rmlib)
        if self.rmout:
            self.rmout = os.path.abspath(self.rmout)
        if self.exclude:
            self.exclude = os.path.abspath(self.exclude)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.genome:
            errors.append("基因组文件未指定|Genome file not specified")

        if self.genome and not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 验证可选文件参数|Validate optional file parameters
        if self.cds and not os.path.exists(self.cds):
            errors.append(f"CDS文件不存在|CDS file not found: {self.cds}")

        if self.curatedlib and not os.path.exists(self.curatedlib):
            errors.append(f"筛选库文件不存在|Curated library file not found: {self.curatedlib}")

        if self.rmlib and not os.path.exists(self.rmlib):
            errors.append(f"RepeatModeler库文件不存在|RepeatModeler library file not found: {self.rmlib}")

        if self.rmout and not os.path.exists(self.rmout):
            errors.append(f"RepeatMasker输出文件不存在|RepeatMasker output file not found: {self.rmout}")

        if self.exclude and not os.path.exists(self.exclude):
            errors.append(f"排除区域文件不存在|Exclude regions file not found: {self.exclude}")

        # 验证EDTA路径|Validate EDTA path
        if self.edta_path and not os.path.exists(self.edta_path):
            errors.append(f"EDTA路径不存在|EDTA path not found: {self.edta_path}")

        # 验证参数范围|Validate parameter ranges
        if self.step not in ['all', 'filter', 'final', 'anno', 'ALL', 'FILTER', 'FINAL', 'ANNO']:
            errors.append(f"无效的步骤参数|Invalid step parameter: {self.step}")

        if self.overwrite not in [0, 1]:
            errors.append(f"overwrite参数必须为0或1|overwrite must be 0 or 1: {self.overwrite}")

        if self.sensitive not in [0, 1]:
            errors.append(f"sensitive参数必须为0或1|sensitive must be 0 or 1: {self.sensitive}")

        if self.anno not in [0, 1]:
            errors.append(f"anno参数必须为0或1|anno must be 0 or 1: {self.anno}")

        if self.evaluate not in [0, 1]:
            errors.append(f"evaluate参数必须为0或1|evaluate must be 0 or 1: {self.evaluate}")

        if self.force not in [0, 1]:
            errors.append(f"force参数必须为0或1|force must be 0 or 1: {self.force}")

        if self.maxdiv < 0 or self.maxdiv > 100:
            errors.append(f"maxdiv参数必须在0-100之间|maxdiv must be between 0 and 100: {self.maxdiv}")

        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Thread count must be positive: {self.threads}")

        if self.u <= 0:
            errors.append(f"突变率必须大于0|Mutation rate must be positive: {self.u}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class PanEDTAConfig:
    """PanEDTA配置类|PanEDTA Configuration Class"""

    # 必需参数|Required parameters
    genome_list: str

    # 可选参数|Optional parameters
    cds: Optional[str] = None
    curatedlib: Optional[str] = None
    fl_copy: int = 3
    anno: int = 1
    overwrite: int = 0
    threads: int = 12
    output_dir: str = "./panedta_output"

    # EDTA路径配置|EDTA path configuration
    edta_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置默认EDTA路径|Set default EDTA path
        if self.edta_path is None:
            conda_env = os.environ.get('CONDA_DEFAULT_ENV', '')
            if 'EDTA' in conda_env:
                conda_prefix = os.environ.get('CONDA_PREFIX', '')
                if conda_prefix:
                    self.edta_path = os.path.join(conda_prefix, 'share', 'EDTA')
                else:
                    self.edta_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA'
            else:
                self.edta_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA'

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径参数|Normalize path parameters
        self.genome_list = os.path.abspath(self.genome_list) if self.genome_list else None
        if self.cds:
            self.cds = os.path.abspath(self.cds)
        if self.curatedlib:
            self.curatedlib = os.path.abspath(self.curatedlib)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证必需参数|Validate required parameters
        if not self.genome_list:
            errors.append("基因组列表文件未指定|Genome list file not specified")

        if self.genome_list and not os.path.exists(self.genome_list):
            errors.append(f"基因组列表文件不存在|Genome list file not found: {self.genome_list}")

        # 验证可选文件参数|Validate optional file parameters
        if self.cds and not os.path.exists(self.cds):
            errors.append(f"CDS文件不存在|CDS file not found: {self.cds}")

        if self.curatedlib and not os.path.exists(self.curatedlib):
            errors.append(f"筛选库文件不存在|Curated library file not found: {self.curatedlib}")

        # 验证EDTA路径|Validate EDTA path
        if self.edta_path and not os.path.exists(self.edta_path):
            errors.append(f"EDTA路径不存在|EDTA path not found: {self.edta_path}")

        # 验证参数范围|Validate parameter ranges
        if self.anno not in [0, 1]:
            errors.append(f"anno参数必须为0或1|anno must be 0 or 1: {self.anno}")

        if self.overwrite not in [0, 1]:
            errors.append(f"overwrite参数必须为0或1|overwrite must be 0 or 1: {self.overwrite}")

        if self.fl_copy < 0:
            errors.append(f"fl_copy参数必须大于等于0|fl_copy must be >= 0: {self.fl_copy}")

        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
