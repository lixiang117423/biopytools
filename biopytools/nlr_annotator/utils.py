"""
NLR-Annotator工具函数模块|NLR-Annotator Utility Functions Module
"""

import glob
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# 结果TSV表头|Result TSV header
TSV_HEADER = 'gene_id\tnlr_id\ttype\tstart\tend\tstrand\tmotifs'


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path
    """
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

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


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class NLRLogger:
    """NLR-Annotator日志管理器|NLR-Annotator Logger Manager"""

    def __init__(self, output_dir: str = "./output", log_name: str = "nlr_annotator.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / "99_logs" / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.log_name = log_name  # 用于独立logger命名,避免批量模式劫持共享logger|Unique logger name avoids batch hijack
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        logger = logging.getLogger(self.log_name)  # 每样本独立logger,避免批量模式劫持顶层logger|Per-sample logger, avoids hijacking shared logger
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        class _MaxLevelFilter(logging.Filter):
            def __init__(self, max_level):
                super().__init__()
                self.max_level = max_level

            def filter(self, record):
                return record.levelno <= self.max_level

        # stdout handler - 仅INFO|stdout handler - INFO only
        # INFO→stdout→.out, WARNING+→stderr→.err
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(_MaxLevelFilter(logging.INFO))
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def clean_output(output_file: str, logger: logging.Logger):
    """
    清洗NLR-Annotator输出：加表头、去重并排序motif列表|Clean NLR-Annotator output: add header, deduplicate and sort motifs
    """
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        logger.warning(f"输出文件为空，跳过清洗|Output file is empty, skipping clean")
        return

    with open(output_file, 'r') as f:
        lines = f.read().strip().split('\n')

    cleaned = [TSV_HEADER]
    for line in lines:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 7:
            cleaned.append(line.strip())
            continue
        # motif列：去重、排序、去掉motif_前缀|Deduplicate, sort, remove motif_ prefix
        motifs = sorted(set(m.replace('motif_', '') for m in fields[6].split(',')))
        fields[6] = ','.join(motifs)
        cleaned.append('\t'.join(fields))

    with open(output_file, 'w') as f:
        f.write('\n'.join(cleaned) + '\n')

    logger.info(f"输出已清洗(含表头)|Output cleaned (with header): {output_file}")


def removed_tsv_path(output_file: str) -> str:
    """
    被包含冗余调用留档文件路径|Audit file path for contained-call removals

    Args:
        output_file: 结果TSV路径|Result TSV path
    """
    if output_file.endswith('.tsv'):
        return output_file[:-len('.tsv')] + '.removed.tsv'
    return output_file + '.removed.tsv'


def gff_path_for(output_file: str) -> str:
    """
    GFF3输出文件路径|GFF3 output file path

    Args:
        output_file: 结果TSV路径|Result TSV path
    """
    if output_file.endswith('.tsv'):
        return output_file[:-len('.tsv')] + '.gff'
    return output_file + '.gff'


# GFF3属性值保留字符百分号编码(,仅在motifs多值分隔时按需还原,单值一律编码)
# |Percent-encode GFF3 reserved chars (comma encoded too; motifs joins encoded values itself)
_GFF3_ESCAPES = str.maketrans({
    '%': '%25', ';': '%3B', '=': '%3D', '&': '%26', ',': '%2C',
    '\t': '%09', '\n': '%0A', '\r': '%0D',
})


def _escape_gff3_attr(value: str) -> str:
    """GFF3属性值百分号编码|Percent-encode a GFF3 attribute value"""
    return value.translate(_GFF3_ESCAPES)


def tsv_to_gff(tsv_path: str, gff_path: str, logger: logging.Logger) -> Tuple[int, int]:
    """
    结果TSV转GFF3,与结果表逐行一致|Convert result TSV to GFF3, row-for-row consistent

    一行数据记录输出一条gene记录:seqid=gene_id,坐标沿用TSV(同为1-based),
    NLR类型放nlr_type属性(GFF3第3列必须是SO术语)。TSV是唯一数据源,
    经清洗+冗余过滤后生成,GFF天然与结果表/汇总口径一致
    |Each data row emits one gene record: seqid=gene_id, coordinates verbatim
    from the TSV (both 1-based), NLR type in the nlr_type attribute (GFF3 col3
    must be an SO term). The TSV is the single source of truth, so the GFF
    matches the cleaned+filtered result table by construction.

    幂等:每次运行覆盖重写;断点续传跳过java时也重新生成
    |Idempotent: overwritten on every run; regenerated even when checkpoint skips java

    Args:
        tsv_path: 结果TSV路径|Result TSV path
        gff_path: 输出GFF3路径|Output GFF3 path
        logger: 日志器|Logger

    Returns:
        (写入记录数, 跳过畸形行数)|(records written, malformed rows skipped)
    """
    path = Path(tsv_path)
    if not path.exists():
        logger.warning(f"结果文件不存在，跳过GFF生成|Result file missing, skip GFF: {tsv_path}")
        return (0, 0)

    records = []
    skipped = 0
    for line in path.read_text().strip().split('\n'):
        if not line.strip() or line.startswith('gene_id\t'):
            continue
        fields = line.split('\t')
        # gene_id/nlr_id/type/start/end/strand;第7列motifs可缺(6列行)|col7 motifs optional
        if len(fields) >= 6:
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                start = None
            if start is not None and start <= end:
                seq_id, nlr_id, nlr_type, strand = fields[0], fields[1], fields[2], fields[5]
                motifs = fields[6] if len(fields) >= 7 else ''
                attributes = (f"ID={_escape_gff3_attr(nlr_id)}"
                              f";Name={_escape_gff3_attr(nlr_id)}"
                              f";nlr_type={_escape_gff3_attr(nlr_type)}")
                if motifs:
                    # 逐值转义后用裸逗号join,符合GFF3多值属性惯例|escape per value, join with bare commas per GFF3 multi-value convention
                    escaped_motifs = ','.join(_escape_gff3_attr(m)
                                              for m in motifs.split(',') if m)
                    if escaped_motifs:
                        attributes += f";motifs={escaped_motifs}"
                records.append('\t'.join([
                    seq_id, 'nlr_annotator', 'gene', str(start), str(end),
                    '.', strand, '.', attributes,
                ]))
                continue
        skipped += 1
        logger.warning(f"坐标或列数异常，该行跳过不出GFF|Bad coordinates or columns, row skipped in GFF: {line[:80]}")

    with open(gff_path, 'w') as f:
        f.write('##gff-version 3\n')
        for record in records:
            f.write(record + '\n')

    logger.info(f"GFF已生成(与结果表逐行一致)|GFF generated (row-for-row consistent with TSV): "
                f"{gff_path} ({len(records)} 条记录|records"
                + (f", 跳过|skipped {skipped} 畸形行|malformed rows" if skipped else "") + ")")
    return (len(records), skipped)


def filter_contained_calls(output_file: str, logger: logging.Logger) -> Tuple[int, int]:
    """
    过滤被完整基因包含的冗余NLR调用|Filter redundant NLR calls fully contained in another call

    NLR-Annotator的motif链接算法在密集/串联NLR位点会把同一基因内部的motif子集
    (如TIR-only短片段)单独打包成一条记录;本函数按"同序列上坐标完全包含"关系
    剔除这些冗余调用(链式包含只留最外层,同坐标重复留靠前一条),被剔除记录
    留档到 *.removed.tsv(含contained_by列)。幂等:重复运行对已过滤文件无副作用
    |The motif-chaining step of NLR-Annotator can emit sub-calls of the same gene
    (e.g. TIR-only fragments) at dense/clustered loci; this drops calls whose
    interval is fully contained in another call on the same sequence (chained
    containment keeps the outermost; identical intervals keep the first),
    archiving removals to *.removed.tsv (with a contained_by column). Idempotent.

    Args:
        output_file: 结果TSV路径|Result TSV path
        logger: 日志器|Logger

    Returns:
        (保留行数, 剔除行数)|(kept rows, removed rows)
    """
    path = Path(output_file)
    if not path.exists() or path.stat().st_size == 0:
        logger.warning(f"结果文件不存在或为空，跳过冗余过滤|Result file missing or empty, skip filtering: {output_file}")
        return (0, 0)

    lines = path.read_text().strip().split('\n')

    # 解析数据行;畸形行(坐标非整数/倒置/列数不足)保留不过滤|Parse data rows; malformed rows (non-int/inverted/short coords) kept unfiltered
    records = []  # (line_idx, gene_id, start, end)
    malformed = 0
    for idx, line in enumerate(lines):
        if not line.strip() or line.startswith('gene_id\t'):
            continue
        fields = line.split('\t')
        if len(fields) >= 5:
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                start = end = None
            if start is not None and start <= end:
                records.append((idx, fields[0], start, end))
                continue
        malformed += 1
        logger.warning(f"坐标无法解析，该行保留不过滤|Unparseable coordinates, line kept unfiltered: {line[:80]}")

    # 按序列分组,组内按(start升序,end降序)扫描,维护最大end:命中即被完全包含
    # |Group by sequence; sort (start asc, end desc) and sweep max end: a hit means full containment
    groups = {}  # gene_id -> [record, ...]
    for record in records:
        groups.setdefault(record[1], []).append(record)

    removed_by_container = {}  # line_idx -> container nlr_id
    for seq_records in groups.values():
        ordered = sorted(seq_records, key=lambda r: (r[2], -r[3]))
        max_end = -1
        container_id = None
        for idx, _seq, start, end in ordered:
            if end <= max_end:
                removed_by_container[idx] = container_id
            else:
                max_end = end
                container_id = lines[idx].split('\t')[1]

    removed_count = len(removed_by_container)
    if removed_count:
        kept_lines = [line for idx, line in enumerate(lines) if idx not in removed_by_container]
        path.write_text('\n'.join(kept_lines) + '\n')

        audit_file = removed_tsv_path(output_file)
        audit_lines = [TSV_HEADER + '\tcontained_by']
        for idx, container in sorted(removed_by_container.items()):
            audit_lines.append(lines[idx] + '\t' + container)
            removed_id, removed_seq = lines[idx].split('\t')[1], lines[idx].split('\t')[0]
            logger.info(f"剔除被包含调用|Removed contained call: {removed_id} "
                        f"contained in {container} ({removed_seq})")
        with open(audit_file, 'w') as f:
            f.write('\n'.join(audit_lines) + '\n')
        logger.info(f"冗余留档已写入|Audit file written: {audit_file}")

    kept_count = len(records) + malformed - removed_count
    logger.info(f"冗余包含调用过滤|Contained-call filtering: 剔除|removed {removed_count} "
                f"条|records, 保留|kept {kept_count} 条|records")
    return (kept_count, removed_count)


def extract_sample_name(filename: str, sample_suffix: str) -> str:
    """
    从文件名提取样本名|Extract sample name from filename

    Args:
        filename: 文件名（含扩展名）|Filename with extension
        sample_suffix: 匹配后缀，如"*.cds.fa"|Match suffix, e.g. "*.cds.fa"
    """
    escaped = re.escape(sample_suffix).replace(r'\*', '(.*)')
    match = re.match(escaped, filename)
    if match and match.group(1):
        return match.group(1)
    return Path(filename).stem


def collect_input_files(input_path: str, sample_suffix: str, logger: logging.Logger) -> List[Tuple[str, str]]:
    """
    收集输入文件列表|Collect input file list

    支持单文件或目录批处理|Supports single file or directory batch processing

    Args:
        input_path: 输入文件或目录路径|Input file or directory path
        sample_suffix: 目录模式下文件匹配后缀|File match suffix in directory mode
        logger: 日志器|Logger

    Returns:
        [(file_path, sample_name), ...]
    """
    path = Path(input_path)

    if path.is_file():
        sample_name = path.stem
        logger.info(f"发现单个输入文件|Found single input file: {path.name}")
        return [(str(path), sample_name)]

    if path.is_dir():
        pattern = str(path / sample_suffix)
        files = sorted(glob.glob(pattern))
        if not files:
            raise ValueError(f"目录中未找到匹配'{sample_suffix}'的文件|No files matching '{sample_suffix}' found in directory")

        results = []
        for f in files:
            name = extract_sample_name(Path(f).name, sample_suffix)
            results.append((f, name))

        logger.info(f"发现批处理文件|Found batch files: {len(results)} 个文件|files")
        return results

    raise ValueError(f"输入路径无效|Invalid input path: {input_path}")


# 结果文件后缀(merge-only 模式按此收集已有结果)|Result file suffix for merge-only collection
RESULT_SUFFIX = '.nlr_annotator.tsv'


def _sample_from_result_name(filename: str) -> str:
    """
    从 {sample}.nlr_annotator.tsv 提取样本名|Extract sample name from result filename

    Args:
        filename: 结果文件名(如 S01.nlr_annotator.tsv)|Result filename
    """
    if filename.endswith(RESULT_SUFFIX):
        return filename[:-len(RESULT_SUFFIX)]
    return Path(filename).stem


def collect_result_files(input_path: str, logger: logging.Logger) -> List[Tuple[str, str]]:
    """
    收集已有NLR结果TSV(merge-only用)|Collect existing NLR result TSVs for merge-only

    自动识别两种目录形态|Auto-detects two directory shapes:
      - by-sample(批处理输出): dir/{sample}/{sample}.nlr_annotator.tsv
      - 平铺|flat: dir/{sample}.nlr_annotator.tsv

    Args:
        input_path: 结果文件或目录|Result file or directory
        logger: 日志器|Logger

    Returns:
        [(sample_name, tsv_path), ...] 按样本名排序|sorted by sample name

    Raises:
        ValueError: 目录中无任何结果文件|No result files in directory
    """
    path = Path(input_path)

    if path.is_file():
        sample = _sample_from_result_name(path.name)
        logger.info(f"发现单个结果文件|Found single result file: {path.name}")
        return [(sample, str(path))]

    if path.is_dir():
        # 同时收集平铺与子目录形态,按绝对路径去重|collect flat + subdir shapes, dedup by abspath
        candidates = set()
        for pattern in ('*.nlr_annotator.tsv', '*/*.nlr_annotator.tsv'):
            for f in glob.glob(str(path / pattern)):
                candidates.add(os.path.abspath(f))

        if not candidates:
            raise ValueError(
                f"目录中未找到结果文件 '*.nlr_annotator.tsv'|"
                f"No '*.nlr_annotator.tsv' result files found in: {input_path}"
            )

        results = [(_sample_from_result_name(os.path.basename(f)), f) for f in candidates]
        results.sort(key=lambda x: x[0])  # 按样本名稳定排序|stable sort by sample name

        logger.info(f"发现结果文件|Found result files: {len(results)} 个文件|files")
        return results

    raise ValueError(f"输入路径无效|Invalid input path: {input_path}")


def generate_summary(sample_results: List[Tuple[str, str]], output_path: Path, logger: logging.Logger):
    """
    生成多样本汇总文件|Generate multi-sample summary file

    Args:
        sample_results: [(sample_name, tsv_path), ...]
        output_path: 输出目录|Output directory
        logger: 日志器|Logger
    """
    summary_file = output_path / "nlr_annotator_summary.tsv"

    all_rows = ['gene_id\tnlr_id\tsample\ttype\tstart\tend\tstrand\tmotifs']
    for sample_name, tsv_path in sample_results:
        if not os.path.exists(tsv_path):
            continue
        with open(tsv_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('gene_id'):
                    continue
                fields = line.split('\t')
                if len(fields) >= 6:
                    row = f"{fields[0]}\t{fields[1]}\t{sample_name}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}"
                    if len(fields) >= 7:
                        row += f"\t{fields[6]}"
                    all_rows.append(row)

    with open(summary_file, 'w') as f:
        f.write('\n'.join(all_rows) + '\n')

    logger.info(f"汇总文件已生成|Summary file generated: {summary_file}")
    logger.info(f"汇总NLR记录数|Total NLR records in summary: {len(all_rows) - 1}")
