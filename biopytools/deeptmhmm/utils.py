"""
DeepTMHMM 1.0跨膜螺旋/信号肽预测工具函数模块|DeepTMHMM 1.0 Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Optional, List


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect conda env for a command, return env name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 优先从命令完整路径检测|First detect from full command path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 兜底: 搜索所有conda环境|Fallback: search all conda envs
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', os.path.basename(command))
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令|Build conda run command

    必须传完整路径(含/envs/), 不能用basename, 否则无法识别环境|Must pass full path
    (containing /envs/), never a basename, or env detection fails

    Args:
        command: 命令完整路径|Full command path
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class DeeptmhmmLogger:
    """DeepTMHMM日志管理器|DeepTMHMM Logger Manager"""

    def __init__(self, output_dir, log_name: str = "deeptmhmm.log"):
        self.output_dir = output_dir
        self.log_file = os.path.join(str(output_dir), log_name)
        self._setup_logging()

    def _setup_logging(self):
        """设置日志(stdout INFO + stderr WARNING+ + 文件 DEBUG)|Setup logging"""
        self.logger = logging.getLogger("Deeptmhmm")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        return self.logger


# predict.py输出文件 -> 最终目录规范命名映射|predict.py output -> final dir naming map
_OUTPUT_NAME_MAP = {
    'predicted_topologies.3line': '{prefix}_deeptmhmm_topologies.3line',
    'TMRs.gff3': '{prefix}_deeptmhmm_tmr.gff3',
    'deeptmhmm_results.md': '{prefix}_deeptmhmm_results.md',
}


def run_deeptmhmm(logger: logging.Logger, config) -> bool:
    """
    运行DeepTMHMM预测|Run DeepTMHMM prediction

    predict.py要求--output-dir不能预先存在, 因此先在输出目录下建唯一父目录,
    传一个尚不存在的子目录给predict.py由其创建, 跑完再把结果搬到最终目录|predict.py
    requires --output-dir not pre-exist, so create a unique parent under output dir,
    pass a non-existent subdir for predict.py to create, then move results to final dir.

    Args:
        logger: 日志对象|Logger
        config: DeeptmhmmConfig实例|DeeptmhmmConfig instance

    Returns:
        是否成功|Whether successful
    """
    # 唯一父目录(建于输出目录所在文件系统, 避免大蛋白组撑爆/tmp)
    # Unique parent on the same filesystem as output (avoids filling /tmp for large proteomes)
    parent_dir = tempfile.mkdtemp(prefix='_deeptmhmm_run_', dir=config.output_dir)
    work_dir = os.path.join(parent_dir, 'predict_out')  # 尚不存在, 交给predict.py创建|not created yet
    logger.info(f"临时工作目录|Temp work dir: {work_dir}")

    # 构建命令: conda run -n <env> --no-capture-output <python_bin> <predict_py> --fasta .. --output-dir ..
    args = [
        config.predict_py,
        '--fasta', config.input_file,
        '--output-dir', work_dir,
    ]
    cmd = build_conda_command(config.python_bin, args)

    logger.info("执行|Executing: DeepTMHMM预测|DeepTMHMM prediction")
    logger.info(f"命令|Command: {' '.join(cmd)}")

    success = False
    try:
        # 不捕获输出, 让predict.py进度直接流式写到作业的.out (长任务实时可见)
        # Do not capture; let predict.py progress stream to the job's .out (visible in real time)
        result = subprocess.run(cmd, check=False, cwd=config.deeptmhmm_dir)

        if result.returncode != 0:
            logger.error(
                f"DeepTMHMM执行失败|DeepTMHMM execution failed "
                f"(exit code: {result.returncode})"
            )
            return False

        success = True

        # 搬运结果到最终目录|Move results to final dir
        moved = _move_outputs(work_dir, config)
        for dst in moved:
            logger.info(f"已生成|Generated: {os.path.basename(dst)}")

        missing = [
            src_name for src_name in _OUTPUT_NAME_MAP
            if not os.path.exists(os.path.join(work_dir, src_name))
        ]
        if missing:
            logger.warning(f"部分输出未生成|Some outputs not generated: {', '.join(missing)}")

        return True

    except FileNotFoundError as e:
        logger.error(f"命令未找到|Command not found: {e}")
        logger.error("通常是conda不在PATH, 确认计算节点conda已初始化|Usually conda not in PATH, ensure conda is initialized on the compute node")
        return False
    finally:
        # 清理临时父目录(含predict.py产生的probabilities/embeddings等中间文件)
        # Clean temp parent (including predict.py intermediates like probabilities/embeddings)
        try:
            shutil.rmtree(parent_dir)
        except Exception as e:
            logger.warning(f"清理临时目录失败|Failed to clean temp dir: {e}")
            logger.warning(f"可手动删除|You can remove it manually: {parent_dir}")


def _move_outputs(work_dir: str, config) -> List[str]:
    """把predict.py输出搬到最终目录并按规范改名|Move predict.py outputs to final dir with spec naming"""
    moved = []
    for src_name, dst_template in _OUTPUT_NAME_MAP.items():
        src = os.path.join(work_dir, src_name)
        dst = os.path.join(config.output_dir, dst_template.format(prefix=config.output_prefix))
        if os.path.exists(src):
            shutil.move(src, dst)
            moved.append(dst)
    return moved


def parse_deeptmhmm_output(
    topologies_3line: str,
    tmr_gff3: str,
    logger: logging.Logger = None,
) -> list:
    """
    解析DeepTMHMM输出|Parse DeepTMHMM outputs

    Args:
        topologies_3line: predicted_topologies.3line 文件路径|3line file path
        tmr_gff3: TMRs.gff3 文件路径|gff3 file path
        logger: 日志对象|Logger

    Returns:
        记录列表, 每条含 id/length/protein_type/pred_tmhs/signal_peptide/tm_regions|Record list
    """
    # 解析3line: >id | TYPE / SEQUENCE / LABELS
    records_order = []  # 保序|preserve order
    three_line = {}
    try:
        with open(topologies_3line, 'r') as fh:
            lines = [ln.rstrip('\n') for ln in fh]
    except FileNotFoundError:
        if logger:
            logger.error(f"三行拓扑文件不存在|3line file not found: {topologies_3line}")
        return []

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('>'):
            # >id | TYPE
            header = line[1:]
            if ' | ' in header:
                seq_id, ptype = header.rsplit(' | ', 1)
            else:
                seq_id, ptype = header.strip(), ''
            seq_id = seq_id.strip()
            # 下一行是序列|next line is sequence
            length = 0
            if i + 1 < len(lines) and not lines[i + 1].startswith('>'):
                length = len(lines[i + 1].strip())
            three_line[seq_id] = {
                'length': length,
                'protein_type': ptype.strip(),
            }
            records_order.append(seq_id)
            i += 3
            continue
        i += 1

    # 解析gff3: 统计TMhelix / SignalPeptide
    gff = {}
    try:
        with open(tmr_gff3, 'r') as fh:
            for raw in fh:
                line = raw.rstrip('\n')
                if not line or line.startswith('#') or line.strip() == '//':
                    continue
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                seqid, ftype, start, end = fields[0], fields[2], fields[3], fields[4]
                rec = gff.setdefault(seqid, {'tm_regions': [], 'signal': None})
                if ftype == 'TMhelix':
                    rec['tm_regions'].append(f"{start}-{end}")
                elif ftype == 'SignalPeptide':
                    rec['signal'] = (start, end)
    except FileNotFoundError:
        if logger:
            logger.warning(f"GFF3文件不存在|GFF3 file not found: {tmr_gff3}")

    # 合并|Merge
    records = []
    for seq_id in records_order:
        tl = three_line[seq_id]
        g = gff.get(seq_id, {'tm_regions': [], 'signal': None})
        tm_regions = g['tm_regions']
        if g['signal']:
            signal_peptide = f"yes ({g['signal'][0]}-{g['signal'][1]})"
        else:
            signal_peptide = 'no'
        records.append({
            'id': seq_id,
            'length': tl['length'],
            'protein_type': tl['protein_type'],
            'pred_tmhs': len(tm_regions),
            'signal_peptide': signal_peptide,
            'tm_regions': ';'.join(tm_regions) if tm_regions else '-',
        })

    return records


def write_clean_tsv(records: list, output_file: str):
    """
    写入整理后的TSV|Write cleaned TSV

    Args:
        records: parse_deeptmhmm_output 返回的记录|records from parse_deeptmhmm_output
        output_file: 输出文件路径|output file path
    """
    header = [
        'ID', 'Length', 'Protein_Type',
        'Pred_TMHs', 'Signal_Peptide', 'TM_Regions',
    ]
    with open(output_file, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in records:
            row = [
                r['id'],
                str(r['length']),
                r['protein_type'],
                str(r['pred_tmhs']),
                r['signal_peptide'],
                r['tm_regions'],
            ]
            fh.write('\t'.join(row) + '\n')
