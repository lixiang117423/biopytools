"""
选择性扫荡检测配置管理模块|Selective Sweep Detection Configuration Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path


@dataclass
class SweepMergeConfig:
    """选择性扫荡检测配置类|Selective Sweep Detection Config Class"""

    # 必需参数|Required parameters
    input_vcf: str          # 输入VCF|Input VCF (支持.gz|.gz supported)
    pop_info: str           # 群体分组信息(样品ID<TAB>分组,无表头)|Pop info (sample<TAB>group, no header)

    # 输出参数|Output parameters
    output_dir: str = './selective_sweep_output'

    # 运行参数|Run parameters
    threads: int = 12

    # 窗口与候选区域参数|Window and candidate region parameters
    win: int = 50000            # 窗口大小|Window size
    step: int = 50000           # 窗口步长|Window step (vcftools --*-step)
    top_quantile: float = 0.01  # 候选阈值分位数|Candidate threshold quantile
    merge_gap: int = 100000     # 候选窗口合并最大间隔|Max gap for merging candidate windows

    # 过滤参数|Filtering parameters
    min_maf: float = 0.05       # 最小等位基因频率|Min MAF
    max_missing: float = 0.10   # 最大缺失率|Max missing rate

    # RAiSD参数|RAiSD parameters
    raisd_window: int = 50      # RAiSD SNP窗口|RAiSD SNP window (default 50)
    raisd_min_samples: int = 15  # 低样本量阈值(<此值MU分量默认排除)|Low sample threshold (MU excluded below)
    include_mu_low_n: bool = False  # 低样本群体是否仍加入MU分量|Include MU for low-n pops

    # SweeD参数|SweeD parameters
    sweed_grid: int = 10000     # CLR计算网格点数|CLR grid points
    sweed_min_samples: int = 15  # 低样本量阈值(<此值CLR分量默认排除)|Low sample threshold (CLR excluded below)
    include_sweed_low_n: bool = False  # 低样本群体是否仍加入CLR分量|Include CLR for low-n pops
    sweed_folded: bool = True   # 折叠SFS(无祖先状态时正确)|Folded SFS (correct without ancestral state)

    # 工具路径|Tool paths (conda env selective_sweep)
    bcftools_path: str = '~/miniforge3/envs/selective_sweep/bin/bcftools'
    vcftools_path: str = '~/miniforge3/envs/selective_sweep/bin/vcftools'
    raisd_path: str = '~/miniforge3/envs/selective_sweep/bin/RAiSD'
    sweed_path: str = '~/software/sweed/SweeD-P'

    # 日志|Logging
    log_level: str = 'INFO'

    # 内部目录|Internal directories (filled in __post_init__)
    info_dir: Optional[Path] = None
    filter_dir: Optional[Path] = None
    stats_dir: Optional[Path] = None
    sweep_dir: Optional[Path] = None
    log_dir: Optional[Path] = None
    tmp_dir: Optional[Path] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径(~必须展开)|Expand tool paths (~ must expand)
        self.bcftools_path = expand_path(self.bcftools_path)
        self.vcftools_path = expand_path(self.vcftools_path)
        self.raisd_path = expand_path(self.raisd_path)
        self.sweed_path = expand_path(self.sweed_path)

        # 标准化输入输出路径(展开~并转绝对路径)|Normalize input/output paths
        # 注意|Note: expand_path 对裸名(无/无.)原样返回相对路径,必须无条件 abspath,
        # 否则 tmp_dir 等子目录为相对路径,CommandRunner 以 output_dir 为 cwd 时错位
        self.input_vcf = os.path.abspath(expand_path(self.input_vcf))
        self.pop_info = os.path.abspath(expand_path(self.pop_info))
        self.output_dir = os.path.abspath(expand_path(self.output_dir))

        # 创建输出子目录|Create output subdirectories
        self.info_dir = Path(self.output_dir) / '00_pipeline_info'
        self.filter_dir = Path(self.output_dir) / '01_filter'
        self.stats_dir = Path(self.output_dir) / '02_stats'
        self.sweep_dir = Path(self.output_dir) / '03_sweep'
        self.log_dir = Path(self.output_dir) / '99_logs'
        self.tmp_dir = Path(self.output_dir) / 'tmp'
        for d in [self.info_dir, self.filter_dir, self.stats_dir,
                  self.sweep_dir, self.log_dir, self.tmp_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def validate(self) -> bool:
        """验证配置参数(收集所有错误一次抛出)|Validate params (collect all errors, raise once)"""
        errors = []
        for attr, desc in [('input_vcf', '输入VCF文件|Input VCF file'),
                           ('pop_info', '群体分组文件|Population info file')]:
            if not os.path.exists(getattr(self, attr)):
                errors.append(f"{desc}不存在|Not found: {getattr(self, attr)}")
        if self.win <= 0:
            errors.append(f"窗口大小必须为正数|Window size must be positive: {self.win}")
        if self.step <= 0:
            errors.append(f"窗口步长必须为正数|Window step must be positive: {self.step}")
        if not 0 < self.top_quantile < 1:
            errors.append(f"候选阈值分位数必须在(0,1)之间|Top quantile must be in (0,1): {self.top_quantile}")
        if self.merge_gap < 0:
            errors.append(f"合并间隔不能为负|Merge gap must be >= 0: {self.merge_gap}")
        if not 0 <= self.min_maf < 1:
            errors.append(f"MAF阈值必须在[0,1)之间|Min MAF must be in [0,1): {self.min_maf}")
        if not 0 <= self.max_missing <= 1:
            errors.append(f"缺失率阈值必须在[0,1]之间|Max missing must be in [0,1]: {self.max_missing}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")
        if self.sweed_grid <= 0:
            errors.append(f"SweeD网格点数必须为正数|SweeD grid must be positive: {self.sweed_grid}")
        if self.sweed_min_samples <= 0:
            errors.append(f"SweeD低样本阈值必须为正数|SweeD min samples must be positive: {self.sweed_min_samples}")
        if self.log_level not in {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}:
            errors.append(f"日志级别无效|Invalid log level: {self.log_level}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
