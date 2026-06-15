"""
Swave配置管理模块|Swave Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class SwaveConfig:
    """Swave配置类|Swave Configuration Class"""

    # 必需参数|Required parameters
    assemblies_tsv: str  # 样本组装文件TSV
    ref_fasta: str  # 参考基因组FASTA
    gfa_file: str  # 泛基因组图GFA文件
    gfa_source: str  # GFA来源: minigraph, cactus, pggb

    # 路径配置|Path configuration
    swave_path: str = '~/software/swave/Swave-main'
    output_dir: str = './swave_output'

    # 可选参数|Optional parameters
    decomposed_vcf: Optional[str] = None  # cactus/pggb需要
    output_mode: str = 'auto'  # auto, population, single
    spec_samples: List[str] = None  # 指定样本

    # SV检测参数|SV detection parameters
    min_sv_size: int = 50
    max_sv_size: int = 1000000
    max_sv_comps: int = 5

    # 处理选项|Processing options
    dup_to_ins: bool = False  # 将duplication报告为insertion
    remove_small: bool = False  # 移除小于min_sv_size的节点
    force_reverse: bool = False  # 强制调用反向映射snarls

    # 性能参数|Performance parameters
    threads: int = 12

    # 外部工具路径|External tool paths
    minigraph_path: str = 'minigraph'
    gfatools_path: str = 'gfatools'

    # 高级选项|Advanced options
    spec_snarl: Optional[str] = None  # 只调用特定snarl
    spec_path: Optional[str] = None  # 只调用特定path

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径并转为绝对路径（subprocess的cwd会改变工作目录，相对路径必须提前解析）
        # Expand paths and convert to absolute paths (subprocess cwd changes working dir)
        self.swave_path = os.path.abspath(expand_path(self.swave_path))
        self.assemblies_tsv = os.path.abspath(expand_path(self.assemblies_tsv))
        self.ref_fasta = os.path.abspath(expand_path(self.ref_fasta))
        self.gfa_file = os.path.abspath(expand_path(self.gfa_file))
        self.output_dir = os.path.abspath(expand_path(self.output_dir))

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 规范化spec_samples|Normalize spec_samples
        if self.spec_samples is None:
            self.spec_samples = ["all"]

        # 验证gfa_source|Validate gfa_source
        valid_sources = ['minigraph', 'cactus', 'pggb']
        if self.gfa_source not in valid_sources:
            raise ValueError(f"无效的gfa_source|Invalid gfa_source: {self.gfa_source}. "
                           f"必须为|Must be one of {valid_sources}")

        # cactus/pggb需要decomposed_vcf|cactus/pggb require decomposed_vcf
        if self.gfa_source in ['cactus', 'pggb'] and self.decomposed_vcf:
            self.decomposed_vcf = os.path.abspath(expand_path(self.decomposed_vcf))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        required_files = [
            ('样本组装TSV文件|Assemblies TSV file', self.assemblies_tsv),
            ('参考基因组FASTA文件|Reference FASTA file', self.ref_fasta),
            ('GFA文件|GFA file', self.gfa_file),
        ]

        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查Swave路径|Check Swave path
        swave_py = os.path.join(self.swave_path, 'Swave.py')
        if not os.path.exists(swave_py):
            errors.append(f"Swave主程序不存在|Swave main script not found: {swave_py}")

        # 检查cactus/pggb的decomposed_vcf|Check decomposed_vcf for cactus/pggb
        if self.gfa_source in ['cactus', 'pggb']:
            if not self.decomposed_vcf:
                errors.append(f"{self.gfa_source}图需要decomposed_vcf参数|"
                           f"{self.gfa_source} graph requires decomposed_vcf parameter")
            elif not os.path.exists(self.decomposed_vcf):
                errors.append(f"Decomposed VCF文件不存在|Decomposed VCF file not found: {self.decomposed_vcf}")

        # 验证参数范围|Validate parameter ranges
        if self.min_sv_size < 0:
            errors.append(f"min_sv_size必须>=0|min_sv_size must be >= 0")

        if self.max_sv_size < self.min_sv_size:
            errors.append(f"max_sv_size必须>=min_sv_size|max_sv_size must be >= min_sv_size")

        if self.threads < 1:
            errors.append(f"线程数必须>=1|Thread count must be >= 1")

        if self.output_mode not in ['auto', 'population', 'single']:
            errors.append(f"无效的output_mode|Invalid output_mode: {self.output_mode}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
