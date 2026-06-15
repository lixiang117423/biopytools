"""TreeMix配置管理模块|TreeMix Configuration Management Module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class TreemixConfig:
    """TreeMix配置类|TreeMix Configuration Class"""

    # ============ 输入文件 ============
    vcf_file: str = ""                # VCF输入文件 (prepare/all用)
    treemix_input: str = ""           # 预处理的TreeMix输入文件 (run用)
    cluster_file: str = ""            # 样本分组文件; 不提供则自动推断
    pop_delimiter: str = "_"          # 自动推断群体时的分隔符

    # ============ 输出 ============
    output_dir: str = "./treemix_output"

    # ============ 软件路径 (conda环境) ============
    treemix_path: str = field(
        default_factory=lambda: get_tool_path(
            'treemix',
            '~/miniforge3/envs/treemix_v.1.13/bin/treemix',
            'TREEMIX_PATH',
        )
    )
    plink_path: str = field(
        default_factory=lambda: get_tool_path(
            'plink',
            '~/miniforge3/envs/Population_genetics/bin/plink',
            'PLINK_PATH',
        )
    )
    r_path: str = field(
        default_factory=lambda: get_tool_path(
            'R',
            '~/miniforge3/envs/treemix_v.1.13/bin/R',
            'R_PATH',
        )
    )
    bcftools_path: str = field(
        default_factory=lambda: get_tool_path(
            'bcftools',
            '~/miniforge3/envs/bcftools_v.1.22/bin/bcftools',
            'BCFTOOLS_PATH',
        )
    )

    # ============ Prepare参数 (LD过滤) ============
    ld_window: int = 50               # LD窗口大小 (SNP数)
    ld_step: int = 10                 # LD步长
    ld_r2: float = 0.2                # LD r2阈值

    # ============ TreeMix参数 ============
    m_max: int = 10                   # 测试m=0..m_max
    m: Optional[int] = None           # 指定m值 (用于单独运行或绘图)
    root: str = ""                    # 外群群体名
    bootstrap: int = 1000             # bootstrap重复次数
    noss: bool = False                # 关闭样本量校正 (单样本群体时使用)
    global_opt: bool = True           # 全局拓扑重排
    se: bool = False                  # 计算迁移权重标准误 (计算量大)
    k: int = 500                      # SNP block大小
    threads: int = 12                 # 并行bootstrap线程数
    seed: Optional[int] = None        # 随机种子
    replicates: int = 10              # 每个m值的重复次数 (scan阶段)
    plotting_funcs_r: str = ""        # plotting_funcs.R路径 (为空则自动查找)

    # ============ 内部路径 (自动生成) ============
    prepare_dir: Optional[str] = field(default=None, init=False)
    treemix_dir: Optional[str] = field(default=None, init=False)
    logs_dir: Optional[str] = field(default=None, init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.vcf_file = expand_path(self.vcf_file) if self.vcf_file else ""
        self.treemix_input = expand_path(self.treemix_input) if self.treemix_input else ""
        self.cluster_file = expand_path(self.cluster_file) if self.cluster_file else ""
        self.output_dir = expand_path(self.output_dir) if self.output_dir else ""
        self.treemix_path = expand_path(self.treemix_path)
        self.plink_path = expand_path(self.plink_path)
        self.r_path = expand_path(self.r_path)
        self.bcftools_path = expand_path(self.bcftools_path)
        if self.plotting_funcs_r:
            self.plotting_funcs_r = expand_path(self.plotting_funcs_r)

        # 创建输出目录结构|Create output directory structure
        out = Path(self.output_dir)
        self.prepare_dir = str(out / "01_prepare")
        self.treemix_dir = str(out / "02_treemix")
        self.logs_dir = str(out / "99_logs")

        for d in [self.prepare_dir, self.treemix_dir, self.logs_dir]:
            Path(d).mkdir(parents=True, exist_ok=True)

    def validate_prepare(self):
        """验证prepare步骤参数|Validate prepare step parameters"""
        errors = []
        if not self.vcf_file:
            errors.append("需要VCF输入文件 -i/--input|VCF input file required")
        elif not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file not found: {self.vcf_file}")
        if not os.path.exists(self.treemix_path):
            errors.append(f"treemix不存在|treemix not found: {self.treemix_path}")
        if not os.path.exists(self.plink_path):
            errors.append(f"plink不存在|plink not found: {self.plink_path}")
        if self.cluster_file and not os.path.exists(self.cluster_file):
            errors.append(f"分组文件不存在|Cluster file not found: {self.cluster_file}")
        if errors:
            raise ValueError("\n".join(errors))
        return True

    def validate_run(self):
        """验证run步骤参数|Validate run step parameters"""
        errors = []
        if not self.treemix_input:
            errors.append("需要TreeMix输入文件 -i/--input|TreeMix input file required")
        elif not os.path.exists(self.treemix_input):
            errors.append(f"TreeMix输入文件不存在|Input file not found: {self.treemix_input}")
        if not os.path.exists(self.treemix_path):
            errors.append(f"treemix不存在|treemix not found: {self.treemix_path}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
