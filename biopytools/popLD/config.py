"""
PopLDdecay配置管理模块 | PopLDdecay Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class PopLDdecayConfig:
    """PopLDdecay配置类 | PopLDdecay Configuration Class"""

    # 必需参数 | Required parameters
    input_vcf: str
    output_stat: str

    # 软件路径配置 | Software path configuration
    poplddecay_path: str = '/share/org/YZWL/yzwl_lixg/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay'

    # 可选处理参数 | Optional processing parameters
    sub_pop: Optional[str] = None  # 亚群样本列表文件 | Subgroup sample list file
    max_dist: int = 300  # 最大距离(kb) | Max distance in kb
    maf: float = 0.005  # 最小次等位基因频率 | Min minor allele frequency
    het: float = 0.88  # 最大杂合率 | Max ratio of het allele
    miss: float = 0.25  # 最大缺失率 | Max ratio of miss allele
    ehh: Optional[str] = None  # EHH起始位点 | EHH start site

    # 输出控制参数 | Output control parameters
    out_filter_snp: bool = False  # 输出最终SNP | Output final SNP
    out_type: int = 1  # 输出类型(1-8) | Output type
    method: int = 1  # 算法方法(1或2) | Algorithm method (1 or 2)

    # 其他配置 | Other configuration
    threads: int = 1  # 线程数 | Number of threads
    keep_intermediate: bool = False  # 保留中间文件 | Keep intermediate files

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.input_vcf = os.path.normpath(os.path.abspath(self.input_vcf))
        self.poplddecay_path = os.path.normpath(os.path.abspath(self.poplddecay_path))

        # 输出文件路径处理
        output_dir = os.path.dirname(os.path.abspath(self.output_stat))
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        self.output_stat = os.path.normpath(os.path.abspath(self.output_stat))

        # 亚群样本列表文件路径处理
        if self.sub_pop:
            self.sub_pop = os.path.normpath(os.path.abspath(self.sub_pop))

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查输入VCF文件 | Check input VCF file
        if not os.path.exists(self.input_vcf):
            errors.append(f"输入VCF文件不存在 | Input VCF file does not exist: {self.input_vcf}")

        # 检查PopLDdecay软件 | Check PopLDdecay software
        if not os.path.exists(self.poplddecay_path):
            errors.append(f"PopLDdecay软件不存在 | PopLDdecay software does not exist: {self.poplddecay_path}")

        # 检查亚群样本列表文件 | Check subgroup sample list file
        if self.sub_pop and not os.path.exists(self.sub_pop):
            errors.append(f"亚群样本列表文件不存在 | Subgroup sample list file does not exist: {self.sub_pop}")

        # 检查参数范围 | Check parameter ranges
        if self.max_dist <= 0:
            errors.append(f"MaxDist必须大于0 | MaxDist must be > 0: {self.max_dist}")

        if not (0 <= self.maf <= 0.5):
            errors.append(f"MAF必须在0-0.5之间 | MAF must be between 0-0.5: {self.maf}")

        if not (0 <= self.het <= 1):
            errors.append(f"Het必须在0-1之间 | Het must be between 0-1: {self.het}")

        if not (0 <= self.miss <= 1):
            errors.append(f"Miss必须在0-1之间 | Miss must be between 0-1: {self.miss}")

        if self.out_type not in [1, 2, 3, 4, 5, 6, 7, 8]:
            errors.append(f"OutType必须在1-8之间 | OutType must be between 1-8: {self.out_type}")

        if self.method not in [1, 2]:
            errors.append(f"Method必须是1或2 | Method must be 1 or 2: {self.method}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
