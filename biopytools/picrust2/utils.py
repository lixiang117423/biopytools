"""
PICRUSt2工具函数模块|PICRUSt2 Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|conda environment name or None
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
            command_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command_name)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令|Build conda run command

    Args:
        command: 命令路径|Command path
        args: 参数列表|Argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        command_name = os.path.basename(command)
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command_name] + args
    else:
        return [command] + args


class Picrust2Logger:
    """PICRUSt2日志管理器|PICRUSt2 Logger Manager"""

    def __init__(self, log_file_path: str, log_level: str = "INFO"):
        self.log_file = Path(log_file_path)
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)
        formatter = logging.Formatter(log_format, datefmt=date_format)

        logger = logging.getLogger("Picrust2Pipeline")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        # stdout: INFO | stdout: INFO
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr: WARNING+ | stderr: WARNING+
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file: ALL | file: ALL
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = None):
        self.logger = logger
        self.working_dir = working_dir or "."

    def run_command(self, cmd: list, description: str = "") -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description

        Returns:
            (success, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        cmd_str = ' '.join(str(c) for c in cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )
            self.logger.info(f"命令执行成功|Command executed successfully: {description or 'Command'}")
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout[:500]}")
            return True, result.stdout, result.stderr

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description or 'Command'}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                stderr_file = os.path.join(self.working_dir, f"stderr_{int(time.time())}.log")
                try:
                    with open(stderr_file, 'w', encoding='utf-8') as f:
                        f.write(e.stderr)
                    self.logger.error(f"完整错误信息已写入|Full stderr written to: {stderr_file}")
                    tail = e.stderr[-2000:] if len(e.stderr) > 2000 else e.stderr
                    self.logger.error(f"错误信息(尾部)|Error message (tail): {tail}")
                except Exception:
                    self.logger.error(f"错误信息|Error message: {e.stderr[:5000]}")
            return False, e.stdout or "", e.stderr or ""

        except Exception as e:
            self.logger.error(f"执行异常|Execution exception: {e}")
            return False, "", str(e)


def format_number(num: int) -> str:
    """格式化数字为大单位|Format number to large units"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def detect_input_format(filepath: str) -> str:
    """
    检测输入特征表文件格式|Detect input feature table file format

    支持格式|Supported formats:
        - BIOM (.biom)
        - TSV (.tsv, .txt, .csv)
        - Excel (.xlsx, .xls)

    Args:
        filepath: 输入文件路径|Input file path

    Returns:
        str: 格式名称 ('biom', 'tsv', 'excel')
    """
    ext = os.path.splitext(filepath)[1].lower()

    # .gz压缩文件，检查内层扩展名
    if ext == '.gz':
        inner = os.path.splitext(os.path.basename(filepath).rstrip('.gz'))[1].lower()
        if inner == '.tsv' or inner == '.txt':
            return 'tsv'
        return 'tsv'  # gz压缩的默认按tsv处理

    if ext == '.biom':
        return 'biom'
    elif ext in ('.xlsx', '.xls'):
        return 'excel'
    else:
        # .tsv, .txt, .csv 或无扩展名 → 按tsv处理
        # PICRUSt2内部会自动区分tsv和mothur shared格式
        return 'tsv'


def excel_to_tsv(excel_path: str, output_dir: str, logger=None) -> str:
    """
    将Excel特征表转换为TSV格式|Convert Excel feature table to TSV format

    要求|Requirements:
        - 第一列为序列ID (行名)|First column is sequence ID (row name)
        - 第一行为样本名 (列名)|First row is sample names (column names)
        - 其余单元格为数值 (read counts)|Remaining cells are numeric (read counts)

    Args:
        excel_path: Excel文件路径|Excel file path
        output_dir: 输出目录|Output directory
        logger: 日志器|Logger

    Returns:
        str: 转换后的TSV文件路径|Path to converted TSV file
    """
    import pandas as pd

    if logger:
        logger.info(f"检测到Excel输入，转换为TSV格式|Detected Excel input, converting to TSV: {excel_path}")

    ext = os.path.splitext(excel_path)[1].lower()

    # 根据扩展名选择引擎，xlrd不可用时给出明确提示
    # Select engine by extension, give clear message if xlrd unavailable
    if ext == '.xls':
        try:
            import xlrd  # noqa: F401
            engine = 'xlrd'
        except ImportError:
            raise RuntimeError(
                f"读取.xls格式需要xlrd库，请将文件另存为.xlsx格式，或安装xlrd: pip install xlrd|"
                f"Reading .xls requires xlrd. Please save as .xlsx or install: pip install xlrd"
            )
    else:
        engine = 'openpyxl'

    df = pd.read_excel(excel_path, dtype=str, engine=engine)

    # 确保第一列作为行索引|Ensure first column is set as index
    first_col = str(df.columns[0])
    df.set_index(first_col, drop=True, inplace=True)

    # 检查是否有Unnamed列（Excel空列导致的问题）
    # Check for Unnamed columns (caused by empty columns in Excel)
    unnamed_cols = [c for c in df.columns if str(c).startswith('Unnamed')]
    if unnamed_cols:
        if logger:
            logger.warning(f"移除 {len(unnamed_cols)} 个Unnamed列|Removing {len(unnamed_cols)} Unnamed columns")
        df.drop(columns=unnamed_cols, inplace=True)

    # 检查数据是否为数值|Check if data is numeric
    for col in df.columns:
        if col == df.index.name:
            continue
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

    # 生成输出文件名|Generate output file name
    base_name = os.path.splitext(os.path.basename(excel_path))[0].replace('.xls', '')
    tsv_path = os.path.join(output_dir, f"{base_name}.tsv")

    df.to_csv(tsv_path, sep='\t', index=True)

    if logger:
        logger.info(f"Excel转TSV完成|Excel to TSV conversion completed: {tsv_path}")
        logger.info(f"序列数|Sequences: {len(df)}, 样本数|Samples: {len(df.columns)}")

    return tsv_path


def prepare_input_table(input_path: str, output_dir: str, logger=None) -> str:
    """
    自动识别输入特征表格式并准备PICRUSt2可用文件|
    Auto-detect input table format and prepare PICRUSt2-compatible file

    Args:
        input_path: 输入文件路径|Input file path
        output_dir: 输出目录(用于存放转换后的文件)|Output dir for converted files
        logger: 日志器|Logger

    Returns:
        str: PICRUSt2可用的文件路径|PICRUSt2-compatible file path
    """
    fmt = detect_input_format(input_path)

    if logger:
        logger.info(f"输入特征表格式|Input table format: {fmt}")

    if fmt == 'excel':
        return excel_to_tsv(input_path, output_dir, logger)
    else:
        return input_path


def _extract_fasta_ids(fasta_path: str) -> set:
    """
    从FASTA文件中提取序列ID|Extract sequence IDs from FASTA file

    Args:
        fasta_path: FASTA文件路径|FASTA file path

    Returns:
        set: 序列ID集合|Set of sequence IDs
    """
    ids = set()
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].split()[0].strip())
    return ids


def clean_fasta_headers(fasta_path: str, output_dir: str, logger=None) -> str:
    """
    清理FASTA文件的header，只保留ID部分（空格前的内容）|
    Clean FASTA headers, keep only the ID portion (before first whitespace)

    PICRUSt2要求FASTA header中的ID与BIOM表observation ID完全一致，
    如果header包含描述信息（如'>OTU1 T1R1_10381'），需要去除描述部分。
    PICRUSt2 requires FASTA header IDs to exactly match BIOM observation IDs.
    If headers contain descriptions (e.g. '>OTU1 T1R1_10381'), strip them.

    Args:
        fasta_path: FASTA文件路径|FASTA file path
        output_dir: 输出目录|Output directory
        logger: 日志器|Logger

    Returns:
        str: 原文件路径（无需清理时）或清理后文件路径|Original or cleaned FASTA path
    """
    needs_clean = False
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].rstrip('\n\r')
                if ' ' in header or '\t' in header:
                    needs_clean = True
                    break

    if not needs_clean:
        return fasta_path

    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    out_path = os.path.join(output_dir, f"{base_name}_cleaned.fasta")

    with open(fasta_path, 'r') as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                fout.write(f">{seq_id}\n")
            else:
                fout.write(line)

    if logger:
        logger.info(f"FASTA header已清理|FASTA headers cleaned: {out_path}")
    return out_path


def preprocess_biom_table(biom_path: str, fasta_path: str, output_dir: str, logger=None) -> str:
    """
    预处理BIOM特征表，自动检测并修复轴方向|
    Preprocess BIOM table, auto-detect and fix axis orientation

    检查BIOM表的行名(observation IDs)是否与FASTA序列ID匹配，
    如果不匹配则尝试转置表格（行样本名列OTU ID）。
    Check if BIOM row IDs match FASTA sequence IDs,
    if not try transposing (rows=samples, columns=OTU IDs).

    Args:
        biom_path: BIOM文件路径|BIOM file path
        fasta_path: FASTA文件路径|FASTA file path
        output_dir: 输出目录|Output directory
        logger: 日志器|Logger

    Returns:
        str: 处理后的BIOM文件路径|Path to processed BIOM file
    """
    import biom

    if logger:
        logger.info(f"预处理BIOM特征表|Preprocessing BIOM table: {biom_path}")

    # 读取BIOM表|Read BIOM table
    table = biom.load_table(biom_path)
    obs_ids = set(str(id_) for id_ in table.ids(axis='observation'))

    # 读取FASTA ID|Read FASTA IDs
    fasta_ids = _extract_fasta_ids(fasta_path)
    if logger:
        logger.info(f"FASTA序列数|FASTA sequences: {len(fasta_ids)}")

    # 检查行名(observation IDs)与FASTA ID的交集
    # Check overlap between row IDs and FASTA IDs
    overlap = obs_ids & fasta_ids
    if logger:
        logger.info(f"BIOM行名(observation IDs)与FASTA ID交集|BIOM row-FASTA ID overlap: {len(overlap)}")

    if len(overlap) > len(fasta_ids) * 0.5:
        # 超过一半匹配，不需要转置|More than half match, no transpose needed
        if logger:
            logger.info(f"BIOM表方向正确，无需转置|BIOM table orientation correct, no transpose needed")
        return biom_path

    # 行名不匹配，尝试转置|Row IDs don't match, try transpose
    if logger:
        logger.info(f"BIOM行名与FASTA ID不匹配，尝试转置表格|BIOM rows don't match FASTA IDs, attempting transpose")
        logger.info(f"BIOM原始维度|BIOM original shape: {table.shape[0]} observations x {table.shape[1]} samples")

    table_t = table.transpose()
    obs_ids_t = set(str(id_) for id_ in table_t.ids(axis='observation'))
    overlap_t = obs_ids_t & fasta_ids

    if logger:
        logger.info(f"转置后维度|Transposed shape: {table_t.shape[0]} observations x {table_t.shape[1]} samples")
        logger.info(f"转置后行名与FASTA ID交集|Transposed row-FASTA ID overlap: {len(overlap_t)}")

    if len(overlap_t) > len(fasta_ids) * 0.5:
        # 转置后匹配，保存转置后的表|Match after transpose, save transposed table
        base_name = os.path.splitext(os.path.basename(biom_path))[0]
        out_path = os.path.join(output_dir, f"{base_name}_transposed.biom")
        with open(out_path, 'w') as f:
            f.write(table_t.to_json('biopytools_picrust2_wrapper'))
        if logger:
            logger.info(f"BIOM表已转置并保存|BIOM table transposed and saved: {out_path}")
        return out_path

    # 转置后仍不匹配，报错|Still no match after transpose, error
    raise ValueError(
        f"BIOM特征表与FASTA序列ID不匹配|BIOM table and FASTA sequence IDs don't match\n"
        f"BIOM行名示例|BIOM row IDs sample: {list(obs_ids)[:5]}\n"
        f"BIOM列名示例|BIOM column IDs sample: {list(obs_ids_t)[:5]}\n"
        f"FASTA ID示例|FASTA IDs sample: {list(fasta_ids)[:5]}\n"
        f"请检查特征表和序列文件的ID是否一致|Please verify that table and sequence IDs are consistent"
    )


def format_pathway_table(pathway_dir: str, logger=None) -> str:
    """
    格式化通路丰度表：添加通路名称、取整、均值，原地整合为一个TSV文件|
    Format pathway abundance table: add pathway names, round values, mean column.
    Raw PICRUSt2 output is preserved as _raw.tsv.

    Args:
        pathway_dir: 通路输出目录|Pathway output directory
        logger: 日志器|Logger

    Returns:
        str: 格式化后的TSV文件路径|Path to formatted TSV file
    """
    import pandas as pd
    import gzip

    src = os.path.join(pathway_dir, 'path_abun_unstrat.tsv')
    if not os.path.exists(src):
        if logger:
            logger.warning(f"通路丰度表不存在|Pathway abundance table not found: {src}")
        return None

    # 读取MetaCyc通路名称映射
    map_path = os.path.expanduser(
        '~/miniforge3/envs/picrust_v.2.6.3/lib/python3.12/site-packages/'
        'picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz'
    )
    name_map = {}
    if os.path.exists(map_path):
        with gzip.open(map_path, 'rt') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t', 1)
                if len(parts) == 2:
                    name_map[parts[0]] = parts[1]

    # 保留原始文件|Preserve raw file
    raw_path = os.path.join(pathway_dir, 'path_abun_unstrat_raw.tsv')
    os.rename(src, raw_path)

    # 读取通路丰度表
    df = pd.read_csv(raw_path, sep='\t')
    sample_cols = [c for c in df.columns if c != 'pathway']

    # 添加通路名称|Add pathway readable name
    df.insert(1, 'Pathway Name',
              df['pathway'].map(name_map).fillna(df['pathway']))

    # 数值取整|Round values
    for col in sample_cols:
        df[col] = df[col].round(2)

    # 添加样本均值|Add sample mean
    df['Mean'] = df[sample_cols].mean(axis=1).round(2)

    # 按均值降序排列|Sort by mean descending
    df = df.sort_values('Mean', ascending=False).reset_index(drop=True)

    # 写入格式化后的TSV|Write formatted TSV
    df.to_csv(src, sep='\t', index=False)

    if logger:
        logger.info(f"通路丰度表已格式化|Pathway table formatted: {src} (原始数据保留为 {raw_path})")
    return src


def reorganize_outputs(work_dir: str, output_config: dict, logger=None) -> bool:
    """
    将PICRUSt2原始输出重组到标准编号目录|Reorganize PICRUSt2 outputs into standard numbered directories

    Args:
        work_dir: PICRUSt2 pipeline的原始输出目录|PICRUSt2 pipeline raw output directory
        output_config: 目标目录映射|Target directory mapping
            {
                'placement_dir': '01_placement',
                'hsp_dir': '02_hsp',
                'metagenome_dir': '03_metagenome',
                'pathway_dir': '04_pathway',
            }
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Whether successful
    """
    if logger:
        logger.info("开始重组输出文件|Starting output reorganization")

    if not os.path.exists(work_dir):
        if logger:
            logger.warning(f"PICRUSt2工作目录不存在|Work directory not found: {work_dir}")
        return False

    # 定义文件到目标目录的映射规则
    # Define file-to-directory mapping rules
    file_mapping = {}

    # 序列放置文件|Sequence placement files
    placement_files = []
    for f in os.listdir(work_dir):
        if f.endswith('.tre'):
            placement_files.append(f)

    # intermediate/ 目录中的放置相关文件
    interm_dir = os.path.join(work_dir, 'intermediate')
    placement_subdirs = ['place_seqs_bac', 'place_seqs_arc', 'place_seqs']
    for subdir in placement_subdirs:
        src = os.path.join(interm_dir, subdir)
        if os.path.exists(src):
            placement_files.append(os.path.join('intermediate', subdir))

    # HSP相关文件|HSP-related files
    hsp_files = []
    for f in os.listdir(work_dir):
        if 'marker_predicted' in f or 'nsti' in f:
            hsp_files.append(f)
        elif f.endswith('_predicted.tsv.gz') and 'metagenome' not in f:
            hsp_files.append(f)

    # 宏基因组预测文件|Metagenome prediction files
    metagenome_files = []
    metagenome_subdir = os.path.join(work_dir, 'metagenome_out')
    if os.path.exists(metagenome_subdir):
        for f in os.listdir(metagenome_subdir):
            metagenome_files.append(os.path.join('metagenome_out', f))
    for f in os.listdir(work_dir):
        if 'pred_metagenome' in f or 'seqtab_norm' in f or 'weighted_nsti' in f:
            metagenome_files.append(f)

    # 通路推断文件|Pathway inference files
    pathway_files = []
    pathway_subdir = os.path.join(work_dir, 'pathways_out')
    if os.path.exists(pathway_subdir):
        for f in os.listdir(pathway_subdir):
            pathway_files.append(os.path.join('pathways_out', f))
    for f in os.listdir(work_dir):
        if 'path_abun' in f or 'path_cov' in f:
            pathway_files.append(f)

    # 执行文件移动|Move files
    moved_count = 0
    all_mappings = [
        (placement_files, 'placement_dir'),
        (hsp_files, 'hsp_dir'),
        (metagenome_files, 'metagenome_dir'),
        (pathway_files, 'pathway_dir'),
    ]

    for files, target_key in all_mappings:
        target_dir = output_config.get(target_key)
        if not target_dir:
            continue
        Path(target_dir).mkdir(parents=True, exist_ok=True)

        for f in files:
            src = os.path.join(work_dir, f)
            if os.path.exists(src):
                # 如果是目录，移动整个目录
                if os.path.isdir(src):
                    dst = os.path.join(target_dir, os.path.basename(f))
                    if os.path.exists(dst):
                        shutil.rmtree(dst)
                    shutil.move(src, dst)
                else:
                    dst = os.path.join(target_dir, os.path.basename(f))
                    if os.path.exists(dst):
                        os.remove(dst)
                    shutil.move(src, dst)
                moved_count += 1

    if logger:
        logger.info(f"输出文件重组完成|Output reorganization completed: 移动 {moved_count} 个文件/目录|moved {moved_count} files/directories")

    # 检查work_dir是否为空，清理空目录
    # Check if work_dir is empty, clean up empty dirs
    try:
        remaining = os.listdir(work_dir)
        # 移除空子目录
        for item in remaining:
            item_path = os.path.join(work_dir, item)
            if os.path.isdir(item_path):
                try:
                    shutil.rmtree(item_path)
                except OSError:
                    pass
        # 如果work_dir空了，删除它
        remaining = os.listdir(work_dir)
        if not remaining:
            shutil.rmtree(work_dir)
            if logger:
                logger.info("临时工作目录已清理|Temporary work directory cleaned up")
    except OSError:
        pass

    return True


def generate_software_versions_yml(output_dir: str, config, start_time) -> str:
    """
    生成software_versions.yml|Generate software_versions.yml

    Args:
        output_dir: 输出目录|Output directory
        config: Picrust2Config实例|Picrust2Config instance
        start_time: 开始时间|Start time

    Returns:
        str: 版本文件路径|Version file path
    """
    import yaml
    from datetime import datetime

    end_time = datetime.now()
    runtime_seconds = int((end_time - start_time).total_seconds())

    # 获取PICRUSt2版本|Get PICRUSt2 version
    picrust2_version = "unknown"
    try:
        result = subprocess.run(
            ['conda', 'run', '-n', 'picrust_v.2.6.3', '--no-capture-output',
             'picrust2_pipeline.py', '--version'],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode == 0:
            picrust2_version = result.stdout.strip().split('\n')[0]
    except Exception:
        pass

    info = {
        'pipeline': {
            'name': 'biopytools picrust2',
            'version': '1.0.0'
        },
        'tools': {
            'picrust2': {
                'version': picrust2_version,
                'path': config.picrust2_path
            }
        },
        'parameters': {
            'study_fasta': config.study_fasta,
            'input_table': config.input_table,
            'threads': config.threads,
            'max_nsti': config.max_nsti,
            'in_traits': config.in_traits,
            'placement_tool': config.placement_tool,
            'hsp_method': config.hsp_method,
            'stratified': config.stratified,
            'skip_pathways': config.skip_pathways,
        },
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': runtime_seconds
        }
    }

    info_dir = Path(output_dir) / '00_pipeline_info'
    info_dir.mkdir(parents=True, exist_ok=True)
    info_file = info_dir / 'software_versions.yml'

    with open(info_file, 'w', encoding='utf-8') as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True)

    return str(info_file)
