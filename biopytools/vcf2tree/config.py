"""
VCF转树工具配置管理模块|VCF to Tree Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def _expand_path(path: str) -> str:
    """展开路径中的~和环境变量|Expand ~ and environment variables in path

    本地实现避免导入common.paths（其中import yaml可能挂起）
    |Local implementation to avoid importing common.paths (yaml import may hang)
    """
    expanded = os.path.expandvars(os.path.expanduser(path))
    if '/' not in expanded and '.' not in expanded:
        return expanded
    return os.path.abspath(expanded)


@dataclass
class Vcf2TreeConfig:
    """VCF转树配置类|VCF to Tree Configuration Class"""

    # 必需参数|Required parameters
    input_file: str

    # 建树方法|Tree method（默认IQ-TREE|default IQ-TREE）
    method: str = 'iqtree'  # 'fasttree' or 'iqtree'

    # 输出参数|Output parameters
    output_dir: str = './vcf2tree_output'

    # 运行参数|Run parameters
    threads: int = 12

    # VCF过滤参数|VCF filtering parameters
    min_samples_locus: int = 4
    outgroup: str = ""

    # FastTree参数|FastTree parameters
    fasttree_path: str = 'fasttree'
    fasttree_params: str = ''

    # IQ-TREE参数|IQ-TREE parameters
    # 默认使用conda环境iqtree_v.3.0.1中的IQ-TREE
    # |Default uses IQ-TREE from conda env iqtree_v.3.0.1
    iqtree_path: str = '~/miniforge3/envs/iqtree_v.3.0.1/bin/iqtree'
    iqtree_bootstrap: int = 1000
    iqtree_model: Optional[str] = None  # None = ModelFinder auto

    # 内部属性|Internal attributes
    base_name: str = 'variants'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.fasttree_path = _expand_path(self.fasttree_path)
        self.iqtree_path = _expand_path(self.iqtree_path)

        # 创建输出目录|Create output directories
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 提取基础名称|Extract base name from input
        input_basename = os.path.basename(self.input_file)
        for ext in ['.vcf.gz', '.vcf']:
            if input_basename.endswith(ext):
                input_basename = input_basename[:-len(ext)]
                break
        self.base_name = input_basename

        # 输出目录|Output directories
        self.step1_dir = self.output_path / '01_vcf2fasta'
        self.step2_dir = self.output_path / '02_tree'
        self.log_dir = self.output_path / '99_logs'
        self.info_dir = self.output_path / '00_pipeline_info'

        for d in [self.step1_dir, self.step2_dir, self.log_dir, self.info_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # 关键输出文件|Key output files
        self.snps_fa = self.step1_dir / f"{self.base_name}.snps.fa"
        self.tree_nwk = self.step2_dir / f"{self.base_name}.{self.method}.nwk"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.input_file):
            errors.append(f"输入VCF文件不存在|Input VCF file does not exist: {self.input_file}")

        if self.method not in ('fasttree', 'iqtree'):
            errors.append(f"建树方法必须为fasttree或iqtree|Method must be fasttree or iqtree: {self.method}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if self.min_samples_locus <= 0:
            errors.append(f"最小样本数必须为正整数|Min samples per locus must be positive: {self.min_samples_locus}")

        if self.iqtree_bootstrap < 0:
            errors.append(f"Bootstrap次数不能为负数|Bootstrap replicates cannot be negative: {self.iqtree_bootstrap}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
