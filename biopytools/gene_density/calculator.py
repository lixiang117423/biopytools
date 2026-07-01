"""基因密度计算核心逻辑|Gene density calculator core logic

纯函数实现分箱/密度计算,便于单测;I/O与编排放在 GeneDensityCalculator 类中
|Pure functions for binning/density (easily testable); I/O and orchestration in GeneDensityCalculator
"""

import math
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import format_number

# 模块版本|Module version
_MODULE_VERSION = "1.0.0"


def midpoint(start: int, end: int) -> int:
    """计算基因中点(向下取整)|Compute gene midpoint (floored)

    Args:
        start: 起始位置(1-based)|Start position (1-based)
        end: 终止位置(1-based)|End position (1-based)

    Returns:
        int: 中点坐标|Midpoint coordinate
    """
    return (start + end) // 2


def window_index(pos: int, window_size: int) -> int:
    """计算位置所属窗口索引(1-based闭区间)|Window index for a position (1-based closed)

    窗口0覆盖 [1, window_size],窗口1覆盖 [window_size+1, 2*window_size],以此类推
    |Window 0 covers [1, window_size], window 1 covers [window_size+1, 2*window_size], etc.

    Args:
        pos: 位置(1-based)|Position (1-based)
        window_size: 窗口大小|Window size

    Returns:
        int: 窗口索引(从0开始)|Window index (0-based)
    """
    return (pos - 1) // window_size


def bin_genes(
    genes_by_chrom: Dict[str, List[Tuple[int, int]]],
    window_size: int
) -> Dict[Tuple[str, int], int]:
    """按基因中点分箱,每个基因只计一次|Bin genes by midpoint, each gene counted once

    Args:
        genes_by_chrom: {chrom: [(start, end), ...]}
        window_size: 窗口大小|Window size

    Returns:
        dict: {(chrom, window_index): gene_count}
    """
    binned: Dict[Tuple[str, int], int] = {}
    for chrom, coords in genes_by_chrom.items():
        for start, end in coords:
            idx = window_index(midpoint(start, end), window_size)
            key = (chrom, idx)
            binned[key] = binned.get(key, 0) + 1
    return binned


def compute_chrom_lengths(
    genes_by_chrom: Dict[str, List[Tuple[int, int]]],
    genome_lengths: Dict[str, int]
) -> Dict[str, int]:
    """合并染色体长度:优先genome长度,否则回退到该chrom最大基因end
    |Merge chromosome lengths: prefer genome lengths, else fallback to max gene end

    Args:
        genes_by_chrom: {chrom: [(start, end), ...]}
        genome_lengths: 从.fai/FASTA解析的长度|Lengths parsed from .fai/FASTA

    Returns:
        dict: {chrom: length}
    """
    if genome_lengths:
        result = dict(genome_lengths)
        # genome未覆盖但GFF有基因的染色体:回退到最大基因end
        # |Chroms with genes but absent from genome: fallback to max gene end
        for chrom, coords in genes_by_chrom.items():
            if chrom not in result and coords:
                result[chrom] = max(e for _, e in coords)
        return result

    # 无genome文件:全部回退到最大基因end|No genome file: all fallback to max gene end
    return {chrom: max(e for _, e in coords) for chrom, coords in genes_by_chrom.items() if coords}


def build_window_rows(
    chrom_lengths: Dict[str, int],
    binned: Dict[Tuple[str, int], int],
    window_size: int
) -> List[Dict]:
    """构建窗口行(含空窗口与末尾截断窗口)|Build window rows (incl. empty & truncated last windows)

    Args:
        chrom_lengths: {chrom: length}
        binned: {(chrom, window_index): gene_count} 来自bin_genes()|from bin_genes()
        window_size: 窗口大小|Window size

    Returns:
        list of dict: 每行 {chrom, start, end, gene_count, genes_per_Mb}
    """
    rows: List[Dict] = []
    for chrom, length in chrom_lengths.items():
        n_windows = math.ceil(length / window_size)
        for idx in range(n_windows):
            start = idx * window_size + 1
            end = min((idx + 1) * window_size, length)
            win_len = end - start + 1
            count = binned.get((chrom, idx), 0)
            density = count / (win_len / 1_000_000)
            rows.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'gene_count': count,
                'genes_per_Mb': round(density, 4),
            })
    return rows


def parse_genome_lengths(path) -> Dict[str, int]:
    """解析染色体长度文件(.fai或FASTA)|Parse chromosome length file (.fai or FASTA)

    Args:
        path: 文件路径(.fai优先按扩展名识别,否则按FASTA解析);None返回空
              |File path (.fai detected by extension, else FASTA); None -> empty

    Returns:
        dict: {chrom: length}
    """
    if not path:
        return {}

    try:
        with open(path, 'r') as f:
            if str(path).endswith('.fai'):
                return _parse_fai(f)
            return _parse_fasta(f)
    except Exception:
        return {}


def _parse_fai(f) -> Dict[str, int]:
    """解析.fai文件|Parse .fai file (col0=name, col1=length)"""
    lengths: Dict[str, int] = {}
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            try:
                lengths[parts[0]] = int(parts[1])
            except ValueError:
                continue
    return lengths


def _parse_fasta(f) -> Dict[str, int]:
    """解析FASTA文件统计序列长度|Parse FASTA and count sequence lengths"""
    lengths: Dict[str, int] = {}
    current = None
    n = 0
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if current is not None:
                lengths[current] = n
            current = line[1:].split()[0]
            n = 0
        else:
            n += len(line.strip())
    if current is not None:
        lengths[current] = n
    return lengths


def parse_gff3_file(
    path: str,
    feature_type: str
) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
    """解析GFF3文件,收集目标feature的坐标|Parse GFF3, collect coords of target feature

    Args:
        path: GFF3文件路径|GFF3 file path
        feature_type: 目标feature类型(GFF第三列)|Target feature type (GFF column 3)

    Returns:
        tuple: ({chrom: [(start, end), ...]}, total_count)
    """
    genes_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    total = 0
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != feature_type:
                continue
            chrom = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            genes_by_chrom.setdefault(chrom, []).append((start, end))
            total += 1
    return genes_by_chrom, total


class GeneDensityCalculator:
    """基因密度计算器|Gene density calculator

    组合纯函数完成解析→分箱→密度→TSV→绘图→版本信息的编排
    |Orchestrates pure functions: parse -> bin -> density -> TSV -> plot -> versions
    """

    def __init__(self, config, logger):
        """
        初始化计算器|Initialize calculator

        Args:
            config: GeneDensityConfig|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def run(self) -> List[Dict]:
        """运行完整流水线|Run full pipeline

        Returns:
            list of dict: 窗口行|Window rows [{chrom, start, end, gene_count, genes_per_Mb}]
        """
        start_time = datetime.now()
        self.logger.info("=" * 60)
        self.logger.info("基因密度计算|Gene density calculation")
        self.logger.info("=" * 60)
        self.logger.info(f"GFF文件|GFF file: {self.config.gff_file}")
        self.logger.info(f"窗口大小|Window size: {format_number(self.config.window_size)} bp")
        self.logger.info(f"Feature类型|Feature type: {self.config.feature_type}")
        if self.config.genome_file:
            self.logger.info(f"基因组长度文件|Genome length file: {self.config.genome_file}")

        # 步骤1: 解析GFF|Step 1: parse GFF
        self.logger.info("步骤1|Step 1: 解析GFF文件|Parsing GFF file")
        genes_by_chrom, total = parse_gff3_file(self.config.gff_file, self.config.feature_type)
        self.logger.info(
            f"目标feature总数|Total target features ({self.config.feature_type}): "
            f"{format_number(total)}"
        )
        self.logger.info(f"覆盖染色体|Covered chromosomes: {len(genes_by_chrom)}")

        if total == 0:
            self.logger.warning(
                f"未找到任何 {self.config.feature_type} 记录,将输出空表|"
                f"No {self.config.feature_type} records found, will output empty table"
            )

        # 染色体长度(步骤1子步骤)|Chromosome lengths (sub-step of step 1)
        genome_lengths = parse_genome_lengths(self.config.genome_file)
        if self.config.genome_file and not genome_lengths:
            self.logger.warning(
                f"基因组长度文件读取失败,回退到GFF最大坐标|Genome length file unreadable, "
                f"fallback to max GFF coordinate: {self.config.genome_file}"
            )
        chrom_lengths = compute_chrom_lengths(genes_by_chrom, genome_lengths)

        # 分箱与窗口(步骤1子步骤)|Bin and build windows (sub-step of step 1)
        binned = bin_genes(genes_by_chrom, self.config.window_size)
        rows = build_window_rows(chrom_lengths, binned, self.config.window_size)

        # 步骤2: 写TSV|Step 2: write TSV
        self.logger.info("步骤2|Step 2: 写入密度表|Writing density table")
        self._write_tsv(rows)
        self.logger.info(f"密度表已保存|Density table saved: {self.config.density_tsv}")

        # 步骤3: 绘图|Step 3: plot
        if self.config.generate_plot:
            self.logger.info("步骤3|Step 3: 绘制密度图|Plotting density")
            if self._plot(rows):
                self.logger.info(f"密度图已保存|Density plot saved: {self.config.plot_png}")
            else:
                self.logger.warning(
                    "绘图已跳过(无可绘制数据或matplotlib不可用)|"
                    "Plot skipped (no data or matplotlib unavailable)"
                )

        # 收尾: 版本信息|Finalize: versions (§12.5)
        self._write_versions(start_time)

        self.logger.info("分析完成|Analysis completed")
        return rows

    def _write_tsv(self, rows: List[Dict]):
        """写入TSV密度表|Write TSV density table"""
        Path(self.config.density_tsv).parent.mkdir(parents=True, exist_ok=True)

        header = ['chrom', 'start', 'end', 'gene_count', 'genes_per_Mb']
        lines = ['\t'.join(header)]
        for r in rows:
            lines.append('\t'.join([
                r['chrom'],
                str(r['start']),
                str(r['end']),
                str(r['gene_count']),
                f"{r['genes_per_Mb']:.4f}",
            ]))

        with open(self.config.density_tsv, 'w') as f:
            f.write('\n'.join(lines) + '\n')

    def _plot(self, rows: List[Dict]) -> bool:
        """绘制各染色体基因密度图|Plot per-chromosome gene density

        Returns:
            bool: 是否成功绘图|Whether plot was produced
        """
        try:
            import matplotlib
            matplotlib.use('Agg')  # 无显示环境|headless backend
            import matplotlib.pyplot as plt
        except ImportError:
            return False

        if not rows:
            return False

        # 按染色体分组(保留首次出现顺序)|Group by chrom (preserve first-seen order)
        by_chrom = OrderedDict()
        for r in rows:
            by_chrom.setdefault(r['chrom'], []).append(r)

        n = len(by_chrom)
        fig, axes = plt.subplots(n, 1, figsize=(12, max(2.0, 1.2 * n)), squeeze=False)
        for ax, (chrom, chrom_rows) in zip(axes[:, 0], by_chrom.items()):
            xs = [(r['start'] + r['end']) / 2 for r in chrom_rows]
            ys = [r['genes_per_Mb'] for r in chrom_rows]
            ax.fill_between(xs, ys, step='mid', alpha=0.6, color='#4C72B0')
            ax.set_title(chrom, fontsize=9)
            ax.set_ylabel('genes/Mb', fontsize=8)
            ax.tick_params(labelsize=7)

        axes[-1, 0].set_xlabel('position (bp)', fontsize=8)
        fig.tight_layout()

        Path(self.config.plot_png).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(self.config.plot_png, dpi=150)
        plt.close(fig)
        return True

    def _write_versions(self, start_time: datetime):
        """写入软件版本与参数元数据|Write software versions and parameters (§12.5)"""
        import sys as _sys
        import yaml

        try:
            import matplotlib as _mpl
            mpl_version = _mpl.__version__
        except ImportError:
            mpl_version = "unavailable"

        end_time = datetime.now()
        runtime_seconds = int((end_time - start_time).total_seconds())

        info = {
            'pipeline': {
                'name': 'biopytools gene_density',
                'version': _MODULE_VERSION,
            },
            'tools': {
                'matplotlib': {'version': mpl_version},
            },
            'parameters': {
                'window_size': self.config.window_size,
                'feature_type': self.config.feature_type,
                'gff_file': self.config.gff_file,
                'genome_file': self.config.genome_file,
                'generate_plot': self.config.generate_plot,
            },
            'execution': {
                'python': _sys.version.split()[0],
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime_seconds,
            },
        }

        Path(self.config.versions_yml).parent.mkdir(parents=True, exist_ok=True)
        with open(self.config.versions_yml, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
