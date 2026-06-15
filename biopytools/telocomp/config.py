"""
TeloComp配置管理模块|TeloComp Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class TeloCompConfig:
    """TeloComp配置类|TeloComp Configuration Class"""

    # 必需参数|Required parameters
    genome: str
    output_dir: str

    # 可选参数|Optional parameters
    ont: Optional[str] = None  # ONT数据文件|ONT data file
    hifi: Optional[str] = None  # HiFi数据文件|HiFi data file

    # Filter_1 参数|Filter_1 parameters
    threads: int = 12  # 线程数|Number of threads
    motifs: Optional[List[str]] = None  # 端粒重复序列|Telomeric repeat motifs
    max_break: int = 50  # 最大可容忍断裂长度|Maximum tolerable fracture length
    min_clip: int = 1  # 最小切割长度|Minimum cutting length

    # Filter_2 参数|Filter_2 parameters
    coverage: int = 100  # 覆盖度参数(0-100)|Coverage parameter (0-100)
    parallels: int = 5  # 并行处理参数|Parallel processing parameter
    min_ratio: float = 0.2  # 原始基因组序列长度与reads长度的比例|Ratio of genome to reads

    # Complement 参数|Complement parameters
    motif: str = "CCCTAAA"  # 端粒重复序列(植物默认)|Telomeric repeat sequence (plant default)
    motif_num: int = 7  # 端粒重复序列碱基数|Number of bases in telomere motif

    # 软件路径|Software paths
    conda_env: str = "~/miniforge3/envs/telocomp"
    telocomp_bin: str = "~/software/telocomp/TeloComp-1.0.0/bin"
    genomesyn_bin: str = "~/software/GenomeSyn/GenomeSyn-main/GenomeSyn-1.2.7/bin"

    # 流程控制|Pipeline control
    skip_filter: bool = False  # 跳过Filter步骤|Skip filter steps
    skip_assembly: bool = True  # 跳过Assembly步骤(默认跳过)|Skip assembly step (skip by default)
    skip_complement: bool = False  # 跳过Complement步骤|Skip complement step
    run_visualization: bool = True  # 运行可视化|Run visualization

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        if self.genome:
            self.genome = os.path.normpath(os.path.abspath(self.genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.conda_env = expand_path(self.conda_env)
        self.telocomp_bin = os.path.normpath(os.path.abspath(expand_path(self.telocomp_bin)))
        self.genomesyn_bin = os.path.normpath(os.path.abspath(expand_path(self.genomesyn_bin)))

        # 设置默认motifs|Set default motifs
        # 如果用户没有指定motifs，使用motif参数构建motifs列表
        # 如果用户指定了motif参数但没有指定motifs，自动添加反向互补序列
        if self.motifs is None:
            # 根据motif参数构建motifs列表（包含正向和反向互补）
            self.motifs = [self.motif, self._get_reverse_complement(self.motif)]

    def _get_reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列|Get reverse complement sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, base) for base in reversed(sequence.upper()))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome}")

        # 检查FAI文件|Check FAI file (如果不存在会自动创建|Will be created automatically if not exists)
        # 注释：FAI索引检查移至主程序，将自动创建|FAI index check moved to main, will be auto-created

        # 检查输入数据|Check input data
        if not self.skip_filter:
            if not self.ont and not self.hifi:
                errors.append(f"必须提供ONT或HiFi数据|Must provide ONT or HiFi data")

            if self.ont and not os.path.exists(self.ont):
                errors.append(f"ONT数据文件不存在|ONT data file does not exist: {self.ont}")

            if self.hifi and not os.path.exists(self.hifi):
                errors.append(f"HiFi数据文件不存在|HiFi data file does not exist: {self.hifi}")

        # 检查软件路径|Check software paths
        if not os.path.exists(self.telocomp_bin):
            errors.append(f"TeloComp bin目录不存在|TeloComp bin directory does not exist: {self.telocomp_bin}")

        if not os.path.exists(self.genomesyn_bin):
            errors.append(f"GenomeSyn bin目录不存在|GenomeSyn bin directory does not exist: {self.genomesyn_bin}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if not (0 <= self.coverage <= 100):
            errors.append(f"覆盖度参数必须在0-100之间|Coverage parameter must be between 0 and 100: {self.coverage}")

        if self.motif_num <= 0:
            errors.append(f"端粒重复序列碱基数必须为正数|Motif number must be positive: {self.motif_num}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
