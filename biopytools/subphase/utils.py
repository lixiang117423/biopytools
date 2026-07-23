"""SubPhaser工具函数|SubPhaser utility functions"""

import logging
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import SUBGENOME_LABELS


class SubPhaserLogger:
    """SubPhaser日志管理器|SubPhaser Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger("subphaser")
        logger.setLevel(level)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在conda环境|Detect conda environment from command path"""
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name
    return None


def build_conda_command(conda_env: str, command: str, args: list) -> list:
    """构建conda run命令|Build conda run command"""
    return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args


# ===== 自动模式工具函数|Auto mode utility functions =====

def parse_fasta_info(fasta_path: str) -> List[Tuple[str, int]]:
    """解析FASTA文件获取染色体名称和长度|Parse FASTA for chromosome names and lengths

    Returns:
        [(chrom_name, length), ...] 按长度降序排列
    """
    chroms = []
    current_name = None
    current_len = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    chroms.append((current_name, current_len))
                current_name = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_name is not None:
            chroms.append((current_name, current_len))

    chroms.sort(key=lambda x: x[1], reverse=True)
    return chroms


def generate_auto_config(
    chroms: List[Tuple[str, int]],
    nsg: int,
    output_path: str,
    min_chrom_size: int = 0,
) -> str:
    """自动生成SubPhaser配置文件|Auto-generate SubPhaser config file

    按染色体大小排序后轮询分配到nsg列，大染色体优先配对（同源染色体大小相近）。
    小于min_chrom_size的序列会被过滤掉。
    Distribute chromosomes sorted by size across nsg columns (large chromosomes
    are paired first, since homologous chromosomes tend to be similar in size).
    Sequences smaller than min_chrom_size are excluded.

    Returns:
        生成的配置文件路径
    """
    if min_chrom_size > 0:
        chroms = [(name, length) for name, length in chroms if length >= min_chrom_size]

    n_chroms = len(chroms)
    n_rows = (n_chroms + nsg - 1) // nsg

    rows = []
    for row_idx in range(n_rows):
        columns = []
        for col_idx in range(nsg):
            chrom_idx = row_idx * nsg + col_idx
            if chrom_idx < n_chroms:
                name, _ = chroms[chrom_idx]
                new_id = str(chrom_idx + 1)
                columns.append(f"{new_id}|{name}")
            else:
                columns.append("")
        rows.append("\t".join(columns))

    config_content = "\n".join(rows) + "\n"

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(config_content)

    return output_path


def find_chrom_subgenome_file(output_dir: str, prefix: Optional[str] = None) -> Optional[str]:
    """查找chrom-subgenome.tsv输出文件|Find chrom-subgenome.tsv output file

    SubPhaser可能直接输出到-o目录，也可能在prefix+子目录中。
    """
    candidates = [output_dir]
    if prefix:
        candidates.append(os.path.join(output_dir, prefix + "subphaser_output"))

    for search_dir in candidates:
        if not os.path.isdir(search_dir):
            continue
        for f in os.listdir(search_dir):
            if f.endswith(".chrom-subgenome.tsv"):
                return os.path.join(search_dir, f)

    for entry in os.listdir(output_dir):
        candidate = os.path.join(output_dir, entry)
        if os.path.isdir(candidate):
            for f in os.listdir(candidate):
                if f.endswith(".chrom-subgenome.tsv"):
                    return os.path.join(candidate, f)

    return None


def parse_chrom_subgenome(tsv_path: str) -> Dict[str, str]:
    """解析chrom-subgenome.tsv获取分组结果|Parse chrom-subgenome.tsv for subgenome assignments

    Returns:
        {chromosome_name: subgenome_label} e.g. {"1": "SG1", "6": "SG2"}
    """
    assignments = {}
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                sg = parts[1]
                assignments[chrom] = sg
    return assignments


def build_subgenome_groups(
    assignments: Dict[str, str],
    chroms: List[Tuple[str, int]],
) -> Dict[str, List[Tuple[str, int]]]:
    """按亚基因组分组染色体|Group chromosomes by subgenome

    Args:
        assignments: {chrom_id: subgenome} from parse_chrom_subgenome
        chroms: [(original_name, length), ...] from parse_fasta_info

    Returns:
        {subgenome: [(original_name, length), ...]} 每组内按长度降序
    """
    name_to_length = {name: length for name, length in chroms}

    groups = defaultdict(list)
    for chrom_id, sg in assignments.items():
        original_name = _resolve_original_name(chrom_id, chroms)
        length = name_to_length.get(original_name, 0)
        groups[sg].append((original_name, length))

    for sg in groups:
        groups[sg].sort(key=lambda x: x[1], reverse=True)

    return dict(groups)


def _resolve_original_name(chrom_id: str, chroms: List[Tuple[str, int]]) -> str:
    """将SubPhaser分配的染色体编号还原为原始名称|Resolve SubPhaser chrom ID to original name"""
    for idx, (name, _) in enumerate(chroms):
        if str(idx + 1) == chrom_id:
            return name
    return chrom_id


def assign_chr_names(
    groups: Dict[str, List[Tuple[str, int]]],
    nsg: int,
) -> Dict[str, str]:
    """分配Chr1A/B/C命名|Assign Chr1A/B/C naming

    每个亚基因组内按长度降序排列，同一排名位置的染色体互为同源组。
    Sort within each subgenome by size descending; same rank = homologous group.

    Args:
        groups: {subgenome: [(original_name, length), ...]}
        nsg: number of subgenomes

    Returns:
        {original_name: new_name} e.g. {"CM032900.1": "Chr1A"}
    """
    sg_keys = sorted(groups.keys())
    max_homologs = max(len(v) for v in groups.values()) if groups else 0

    name_map = {}
    for homolog_idx in range(max_homologs):
        for sg_idx, sg_key in enumerate(sg_keys):
            sg_label = SUBGENOME_LABELS[sg_idx] if sg_idx < len(SUBGENOME_LABELS) else f"SG{sg_idx + 1}"
            if homolog_idx < len(groups[sg_key]):
                original_name, _ = groups[sg_key][homolog_idx]
                new_name = f"Chr{homolog_idx + 1}{sg_label}"
                name_map[original_name] = new_name

    return name_map


def write_renamed_fasta(
    fasta_path: str,
    name_map: Dict[str, str],
    output_dir: str,
) -> str:
    """输出重命名后的FASTA文件|Output renamed FASTA file

    所有染色体合并为一个文件，染色体名称替换为Chr1A/B/C格式。
    All chromosomes merged into one file with Chr1A/B/C names.

    Returns:
        输出文件路径
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_path = os.path.join(output_dir, "phased_genome.fa")

    with open(fasta_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                original_name = line[1:].split()[0]
                new_name = name_map.get(original_name, original_name)
                rest = line[1:].split(None, 1)
                description = f" {rest[1]}" if len(rest) > 1 else ""
                fout.write(f">{new_name}{description}")
            else:
                fout.write(line)

    return output_path


def generate_phased_config(
    groups: Dict[str, List[Tuple[str, int]]],
    name_map: Dict[str, str],
    nsg: int,
    output_path: str,
) -> str:
    """生成最终的phasing配置文件|Generate final phased config file

    基于phasing结果，按同源组组织为SubPhaser格式的配置文件。
    Generate SubPhaser-format config based on phasing results.
    """
    sg_keys = sorted(groups.keys())
    max_homologs = max(len(v) for v in groups.values()) if groups else 0

    lines = []
    for homolog_idx in range(max_homologs):
        columns = []
        for sg_idx, sg_key in enumerate(sg_keys):
            sg_label = SUBGENOME_LABELS[sg_idx] if sg_idx < len(SUBGENOME_LABELS) else f"SG{sg_idx + 1}"
            if homolog_idx < len(groups[sg_key]):
                original_name, _ = groups[sg_key][homolog_idx]
                new_name = name_map.get(original_name, original_name)
                columns.append(f"{new_name}|{original_name}")
            else:
                columns.append("")
        lines.append("\t".join(columns))

    with open(output_path, 'w') as f:
        f.write("\n".join(lines) + "\n")

    return output_path


def write_phasing_result(
    groups: Dict[str, List[Tuple[str, int]]],
    name_map: Dict[str, str],
    assignments: Dict[str, str],
    chroms: List[Tuple[str, int]],
    output_path: str,
):
    """输出phasing结果汇总TSV|Write phasing result summary TSV"""
    with open(output_path, 'w') as f:
        f.write("#original_name\tchrom_id\tsubgenome\tphased_name\tlength\n")
        name_to_length = {name: length for name, length in chroms}
        for chrom_id, sg in sorted(assignments.items(), key=lambda x: (x[1], x[0])):
            original_name = _resolve_original_name(chrom_id, chroms)
            new_name = name_map.get(original_name, original_name)
            length = name_to_length.get(original_name, 0)
            f.write(f"{original_name}\t{chrom_id}\t{sg}\t{new_name}\t{length}\n")


def validate_with_parents(
    groups: Dict[str, List[Tuple[str, int]]],
    parental_genomes: List[str],
    genome_fasta: str,
    threads: int = 8,
) -> Dict[str, str]:
    """用父本基因组验证亚基因组分配|Validate subgenome assignments with parental genomes

    通过minimap2比对计算每个亚基因组与各父本的相似度，确定A/B标签。
    Align each subgenome to parental genomes using minimap2 to determine A/B labels.

    Args:
        groups: {subgenome: [(original_name, length), ...]}
        parental_genomes: [parent_A_path, parent_B_path]
        genome_fasta: original genome FASTA path
        threads: number of threads for minimap2

    Returns:
        {original_sg_label: new_sg_label} e.g. {"SG1": "A", "SG2": "B"}
    """
    parent_a, parent_b = parental_genomes[0], parental_genomes[1]
    sg_keys = sorted(groups.keys())

    scores = {}
    for sg_key in sg_keys:
        sg_chroms = [name for name, _ in groups[sg_key]]
        sg_fasta = _extract_chroms_fasta(genome_fasta, sg_chroms)

        score_a = _alignment_coverage(sg_fasta, parent_a, threads)
        score_b = _alignment_coverage(sg_fasta, parent_b, threads)
        scores[sg_key] = (score_a, score_b)

        os.unlink(sg_fasta)

    label_map = {}
    used_labels = set()
    parent_names = [
        os.path.basename(parent_a).replace('.fa', '').replace('.fasta', ''),
        os.path.basename(parent_b).replace('.fa', '').replace('.fasta', ''),
    ]

    for sg_key in sg_keys:
        score_a, score_b = scores[sg_key]
        if score_a >= score_b:
            preferred, alt = 0, 1
        else:
            preferred, alt = 1, 0

        label_a = SUBGENOME_LABELS[preferred] if preferred < len(SUBGENOME_LABELS) else f"SG{preferred + 1}"
        label_b = SUBGENOME_LABELS[alt] if alt < len(SUBGENOME_LABELS) else f"SG{alt + 1}"

        if label_a not in used_labels:
            label_map[sg_key] = label_a
            used_labels.add(label_a)
        else:
            label_map[sg_key] = label_b
            used_labels.add(label_b)

    return label_map


def _extract_chroms_fasta(fasta_path: str, chrom_names: List[str]) -> str:
    """从基因组FASTA中提取指定染色体|Extract specific chromosomes from genome FASTA"""
    import tempfile

    chrom_set = set(chrom_names)
    tmp = tempfile.NamedTemporaryFile(suffix='.fa', delete=False, mode='w', dir=os.path.dirname(fasta_path))

    include = False
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                include = name in chrom_set
            if include:
                tmp.write(line)

    tmp.close()
    return tmp.name


def _alignment_coverage(query_fasta: str, target_fasta: str, threads: int) -> float:
    """计算比对覆盖度|Calculate alignment coverage using minimap2"""
    try:
        paf = subprocess.run(
            ['minimap2', '-x', 'asm5', '-t', str(threads),
             '--secondary=no', target_fasta, query_fasta],
            capture_output=True, text=True, timeout=3600,
        )
        if paf.returncode != 0:
            return 0.0

        total_aligned = 0
        query_len = 0

        with open(query_fasta) as f:
            for line in f:
                if line.startswith('>'):
                    continue
                query_len += len(line.strip())

        for line in paf.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 11:
                total_aligned += int(parts[10])

        return total_aligned / query_len if query_len > 0 else 0.0

    except (subprocess.TimeoutExpired, FileNotFoundError):
        return 0.0
