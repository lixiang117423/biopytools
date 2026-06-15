"""TreeMix工具函数模块|TreeMix Utility Functions Module"""

import gzip
import logging
import os
import sys
from collections import OrderedDict
from pathlib import Path
from typing import List, Optional, Tuple

from .config import TreemixConfig


# ============================================================
# 日志管理器|Logger Manager
# ============================================================

class TreemixLogger:
    """TreeMix日志管理器|TreeMix Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self._setup_logging(log_level)

    def _setup_logging(self, log_level: str):
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger("TreeMix")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file handler - 全部级别
        if self.log_file:
            Path(self.log_file).parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self) -> logging.Logger:
        return self.logger


# ============================================================
# VCF样本解析|VCF Sample Parsing
# ============================================================

def parse_vcf_samples(vcf_file: str) -> List[str]:
    """从VCF文件中提取样本名称|Extract sample names from VCF file

    Args:
        vcf_file: VCF文件路径 (支持 .vcf 和 .vcf.gz)

    Returns:
        样本名称列表
    """
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                if len(fields) > 9:
                    return fields[9:]

    raise ValueError(f"未找到样本信息|No samples found in VCF: {vcf_file}")


def infer_populations(
    samples: List[str],
    delimiter: str = "_",
) -> List[Tuple[str, str]]:
    """从样本名自动推断群体信息|Auto-infer populations from sample names

    Args:
        samples: 样本名称列表
        delimiter: 分隔符, 第一部分作为群体名

    Returns:
        [(sample_id, population), ...] 列表
    """
    result = []
    for sample in samples:
        parts = sample.split(delimiter)
        if len(parts) < 2:
            raise ValueError(
                f"样本名 '{sample}' 无法按分隔符 '{delimiter}' 拆分, "
                f"请提供 --cluster 文件|Sample '{sample}' cannot be split "
                f"by delimiter '{delimiter}', please provide --cluster file"
            )
        pop = parts[0]
        result.append((sample, pop))
    return result


# ============================================================
# PLINK文件操作|PLINK File Operations
# ============================================================

def modify_fam_file(fam_file: str, cluster_data: List[Tuple[str, str]]) -> int:
    """修改plink .fam文件的FID列为群体名|Modify FID column in plink .fam file

    plink从VCF导入后fam中FID=IID=样本名, 需要将FID改为群体名,
    这样plink --within才能按群体正确分组计算频率

    Args:
        fam_file: .fam文件路径
        cluster_data: [(sample_id, population), ...]

    Returns:
        修改的样本数
    """
    sample_to_pop = {sid: pop for sid, pop in cluster_data}

    with open(fam_file, 'r') as f:
        lines = f.readlines()

    modified = 0
    new_lines = []
    for line in lines:
        fields = line.rstrip('\n').split()
        if len(fields) >= 6:
            iid = fields[1]
            if iid in sample_to_pop and fields[0] != sample_to_pop[iid]:
                fields[0] = sample_to_pop[iid]
                modified += 1
            new_lines.append(' '.join(fields) + '\n')
        else:
            new_lines.append(line)

    with open(fam_file, 'w') as f:
        f.writelines(new_lines)

    return modified


def write_pop_cov(
    cluster_data: List[Tuple[str, str]],
    out_file: str,
):
    """写plink --within的pop.cov文件|Write plink --within pop.cov file

    两列格式: FID(population) IID(sample_id)
    FID列须与修改后的fam文件FID一致

    Args:
        cluster_data: [(sample_id, population), ...]
        out_file: 输出文件路径
    """
    with open(out_file, 'w') as f:
        for sample_id, population in cluster_data:
            f.write(f"{population} {sample_id} {population}\n")


def generate_pop_order(
    cluster_data: List[Tuple[str, str]],
    out_file: str,
):
    """生成群体顺序文件|Generate population order file

    从分组数据中提取排序后的唯一群体名列表

    Args:
        cluster_data: [(sample_id, population), ...]
        out_file: 输出文件路径
    """
    pops = sorted(set(pop for _, pop in cluster_data))
    with open(out_file, 'w') as f:
        for pop in pops:
            f.write(f"{pop}\n")
    return pops


# ============================================================
# 频率格式转换 (替代 plink2treemix.py)
# ============================================================

def convert_plink_to_treemix(
    frq_strat_file: str,
    out_file: str,
) -> List[str]:
    """将plink .frq.strat转换为TreeMix输入格式|Convert plink .frq.strat to TreeMix input

    基于plink2treemix.py (https://github.com/taotaoyuan/treemix) 改写为Python 3

    plink .frq.strat格式 (v1.9):
        CHR  SNP  CLST  A1  A2  MAF  NCHROBS

    TreeMix输入格式 (等位基因计数):
        pop1 pop2 pop3 ...
        c1_a,c1_b  c2_a,c2_b  c3_a,c3_b  ...

    Args:
        frq_strat_file: plink --freq --within 生成的 .frq.strat 或 .frq.strat.gz
        out_file: TreeMix格式输出文件 (.gz)

    Returns:
        排序后的群体名列表
    """
    # snp_id -> {pop: (c_minor, c_major)}
    snps = OrderedDict()

    open_func = gzip.open if frq_strat_file.endswith('.gz') else open
    mode = 'rt' if frq_strat_file.endswith('.gz') else 'r'

    with open_func(frq_strat_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('CHR'):
                continue
            fields = line.split()
            if len(fields) < 7:
                continue

            snp_id = fields[1]
            pop = fields[2]
            maf = float(fields[5])
            nchr = int(fields[6])

            c_minor = round(maf * nchr)
            c_major = nchr - c_minor

            if snp_id not in snps:
                snps[snp_id] = {}
            snps[snp_id][pop] = (c_minor, c_major)

    # 收集群体名 (保持出现顺序, 然后排序)
    all_pops = []
    seen = set()
    for pop_dict in snps.values():
        for pop in pop_dict:
            if pop not in seen:
                all_pops.append(pop)
                seen.add(pop)
    all_pops.sort()

    # 写入TreeMix格式
    out_open = gzip.open if out_file.endswith('.gz') else open
    out_mode = 'wt' if out_file.endswith('.gz') else 'w'

    with out_open(out_file, out_mode) as f:
        # 第一行: 群体名
        f.write(' '.join(all_pops) + '\n')
        # 后续行: 每个SNP的等位基因计数
        for snp_id, pop_dict in snps.items():
            counts = []
            for pop in all_pops:
                c1, c2 = pop_dict.get(pop, (0, 0))
                counts.append(f"{c1},{c2}")
            f.write(' '.join(counts) + '\n')

    return all_pops


# ============================================================
# 结果汇总|Result Summary
# ============================================================

def parse_llik(llik_file: str) -> Optional[float]:
    """解析.llik文件获取log likelihood值|Parse .llik file for log likelihood value"""
    if not os.path.exists(llik_file):
        return None
    with open(llik_file, 'r') as f:
        content = f.read()
    for line in content.strip().split('\n'):
        if 'likelihood' in line.lower():
            parts = line.split(':')
            if len(parts) >= 2:
                try:
                    return float(parts[-1].strip())
                except ValueError:
                    pass
    return None
