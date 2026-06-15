"""Phytophthora效应子鉴定工具函数模块|Phytophthora Effector Identification Utility Functions Module"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple


class PhytoEffectorLogger:
    """Phytophthora效应子鉴定日志管理器|Phytophthora Effector Identification Logger Manager"""

    def __init__(self, log_file=None, verbose=False):
        self.log_file = log_file
        self.verbose = verbose
        self.logger = None
        self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器|Setup logger"""
        self.logger = logging.getLogger("phyto_effector")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO if not self.verbose else logging.DEBUG)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        if self.log_file:
            try:
                log_dir = os.path.dirname(os.path.abspath(self.log_file))
                Path(log_dir).mkdir(parents=True, exist_ok=True)
                file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
                file_handler.setLevel(logging.DEBUG)
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)
            except Exception as e:
                self.logger.warning(f"无法创建日志文件|Cannot create log file: {e}")

    def get_logger(self):
        """获取日志记录器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中|Detect if command is in a conda environment"""
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
    """构建conda run命令|Build conda run command"""
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def run_command(cmd: List[str], logger: logging.Logger, description: str = "",
                cwd: str = None, timeout: int = 86400) -> Tuple[bool, str, str]:
    """执行命令|Execute command"""
    if description:
        logger.info(f"执行|Executing: {description}")
    logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd, shell=False, capture_output=True, text=True,
            check=True, cwd=cwd, timeout=timeout
        )
        return True, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {description}")
        logger.error(f"错误输出|Error output: {e.stderr[:500]}")
        return False, e.stdout, e.stderr
    except subprocess.TimeoutExpired:
        logger.error(f"命令超时|Command timeout: {description}")
        return False, '', 'timeout'
    except FileNotFoundError:
        logger.error(f"命令未找到|Command not found: {cmd[0]}")
        return False, '', f'Command not found: {cmd[0]}'


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    if num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def parse_fasta(file_path: str) -> Iterator[Tuple[str, str]]:
    """解析FASTA文件|Parse FASTA file"""
    current_id = None
    current_seq = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id is not None:
                        yield current_id, ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id is not None:
                yield current_id, ''.join(current_seq)
    except Exception as e:
        raise IOError(f"读取FASTA文件失败|Failed to read FASTA file: {e}")


def parse_fasta_to_dict(file_path: str) -> Dict[str, str]:
    """解析FASTA文件为字典|Parse FASTA file to dictionary"""
    return {sid: seq for sid, seq in parse_fasta(file_path)}


def extract_sequences_by_id(fasta_path: str, ids: Set[str], output_path: str) -> int:
    """按ID从FASTA中提取序列|Extract sequences by ID from FASTA"""
    count = 0
    with open(output_path, 'w', encoding='utf-8') as out_f:
        for sid, seq in parse_fasta(fasta_path):
            if sid in ids:
                out_f.write(f">{sid}\n{seq}\n")
                count += 1
    return count


def merge_fasta_files(input_files: List[str], output_path: str) -> Tuple[int, Dict[str, str]]:
    """合并多个FASTA文件并记录来源|Merge multiple FASTA files and track sources"""
    source_map = {}
    total = 0
    with open(output_path, 'w', encoding='utf-8') as out_f:
        for fpath in input_files:
            stem = Path(fpath).stem
            for sid, seq in parse_fasta(fpath):
                out_f.write(f">{sid}\n{seq}\n")
                source_map[sid] = stem
                total += 1
    return total, source_map


def parse_domtblout_hits(domtblout_file: str) -> Set[str]:
    """解析domtblout文件，返回命中目标ID集合|Parse domtblout file, return set of hit target IDs"""
    hits = set()
    try:
        with open(domtblout_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split()
                if len(fields) >= 1:
                    hits.add(fields[0])
    except Exception as e:
        raise IOError(f"解析domtblout文件失败|Failed to parse domtblout file: {e}")
    return hits


def parse_domtblout_details(domtblout_file: str) -> List[Dict]:
    """解析domtblout文件，返回详细信息|Parse domtblout file, return detailed info"""
    results = []
    try:
        with open(domtblout_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split()
                if len(fields) >= 23:
                    results.append({
                        'target': fields[0],
                        'target_accession': fields[1],
                        'query': fields[3],
                        'full_evalue': float(fields[6]),
                        'full_score': float(fields[7]),
                        'full_bias': float(fields[8]),
                        'domain_evalue': float(fields[11]),
                        'domain_score': float(fields[12]),
                        'hmm_from': int(fields[15]),
                        'hmm_to': int(fields[16]),
                        'ali_from': int(fields[17]),
                        'ali_to': int(fields[18]),
                    })
    except Exception as e:
        raise IOError(f"解析domtblout文件失败|Failed to parse domtblout file: {e}")
    return results


def parse_signalp_output(signalp_dir: str, prefix: str = 'combined_input') -> Dict[str, Dict]:
    """解析SignalP输出结果|Parse SignalP output results"""
    result = {}

    # 尝试多种可能的文件名|Try multiple possible filenames
    candidates = [
        os.path.join(signalp_dir, f"{prefix}_summary.txt"),
        os.path.join(signalp_dir, "prediction_results.txt"),
    ]
    # 也搜索目录下任何 *summary* 文件|Also search for any *summary* file
    if not any(os.path.exists(f) for f in candidates):
        for fname in os.listdir(signalp_dir):
            if 'summary' in fname.lower() or 'prediction_results' in fname.lower():
                candidates.append(os.path.join(signalp_dir, fname))

    summary_file = None
    for f in candidates:
        if os.path.exists(f):
            summary_file = f
            break

    if not summary_file:
        return result

    try:
        with open(summary_file, 'r', encoding='utf-8') as f:
            header_found = False
            header_fields = []
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if not header_found:
                    header_found = True
                    header_fields = fields
                    continue
                if len(fields) < 2:
                    continue

                protein_id = fields[0].strip()
                prediction = fields[1].strip() if len(fields) > 1 else 'OTHER'

                # 兼容多种SignalP输出格式|Compatible with multiple SignalP output formats
                # 格式1: "SP(Sec/SPI)" / "OTHER"
                # 格式2: "SP" / "OTHER"
                has_sp = prediction.upper().startswith('SP')

                entry = {
                    'has_signal_peptide': has_sp,
                    'prediction': prediction,
                    'cs_position': '-',
                    'y_score': '-',
                }

                # 尝试解析各列|Try to parse columns
                for i, hf in enumerate(header_fields):
                    if i >= len(fields):
                        break
                    key = hf.strip().lower()
                    if 'cs' in key and 'pos' in key:
                        entry['cs_position'] = fields[i].strip()
                    elif 'sp(sec/spi)' in key:
                        entry['sp_score'] = fields[i].strip()

                result[protein_id] = entry
    except Exception as e:
        raise IOError(f"解析SignalP结果失败|Failed to parse SignalP results: {e}")

    return result


class RxLRMotifScanner:
    """RxLR基序扫描器|RxLR Motif Scanner"""

    RXLR_PATTERNS = {
        'RxLR': re.compile(r'R.LR'),
        'QxLR': re.compile(r'Q.LR'),
        'GxLR': re.compile(r'G.LR'),
    }
    EER_PATTERN = re.compile(r'EER')

    def __init__(self, window_start=20, window_end=120):
        self.window_start = window_start
        self.window_end = window_end

    def scan(self, seq_id: str, sequence: str) -> Dict:
        """扫描单条序列|Scan a single sequence"""
        window = sequence[self.window_start:self.window_end] if len(sequence) > self.window_start else ''

        rxlr_matches = []
        eer_matches = []

        if window:
            for name, pattern in self.RXLR_PATTERNS.items():
                for m in pattern.finditer(window):
                    rxlr_matches.append({
                        'type': name,
                        'position': self.window_start + m.start() + 1,
                        'context': window[max(0, m.start() - 5):m.end() + 5]
                    })
            for m in self.EER_PATTERN.finditer(window):
                eer_matches.append({
                    'type': 'EER',
                    'position': self.window_start + m.start() + 1,
                    'context': window[max(0, m.start() - 5):m.end() + 5]
                })

        has_rxlr = len(rxlr_matches) > 0
        has_eer = len(eer_matches) > 0
        rxlr_types = sorted(set(m['type'] for m in rxlr_matches))

        return {
            'seq_id': seq_id,
            'length': len(sequence),
            'has_rxlr': has_rxlr,
            'has_eer': has_eer,
            'is_candidate': has_rxlr or has_eer,
            'rxlr_types': ';'.join(rxlr_types) if rxlr_types else 'None',
            'rxlr_positions': ';'.join(str(m['position']) for m in rxlr_matches) if rxlr_matches else 'None',
            'eer_positions': ';'.join(str(m['position']) for m in eer_matches) if eer_matches else 'None',
        }


def split_fasta_chunks(fasta_path: str, chunk_size: int = 300) -> List[str]:
    """按序列数分割FASTA文件(SignalP 3.0限制4000序列)|Split FASTA into chunks by sequence count (SP3 limit 4000)"""
    import tempfile
    tmp_dir = tempfile.mkdtemp(prefix='sp3_chunks_')
    chunks = []
    current = []
    seq_count = 0
    chunk_num = 1

    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>') and current:
                if seq_count >= chunk_size:
                    out_path = os.path.join(tmp_dir, f'chunk_{chunk_num:03d}.fa')
                    with open(out_path, 'w') as out:
                        out.writelines(current)
                    chunks.append(out_path)
                    chunk_num += 1
                    current = []
                    seq_count = 0
            current.append(line)
            if line.startswith('>'):
                seq_count += 1
        if current:
            out_path = os.path.join(tmp_dir, f'chunk_{chunk_num:03d}.fa')
            with open(out_path, 'w') as out:
                out.writelines(current)
            chunks.append(out_path)

    return chunks


def parse_signalp3_short(text: str, sprob_threshold: float = 0.9) -> Dict[str, Dict]:
    """解析SignalP 3.0 short格式输出|Parse SignalP 3.0 short format output"""
    results = {}
    for line in text.split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        hmm_fields = parts[-1].split()
        if len(hmm_fields) < 6:
            continue
        pid = hmm_fields[0]
        prediction = hmm_fields[1]
        try:
            sprob = float(hmm_fields[5])
        except (ValueError, IndexError):
            sprob = 0.0

        has_sp = (prediction == 'S' and sprob >= sprob_threshold)
        results[pid] = {
            'has_signal_peptide': has_sp,
            'prediction': f'SP3_HMM_{prediction}',
            'cs_position': '-',
            'y_score': str(sprob),
            'sp_score': str(sprob),
        }
    return results


def _run_signalp3_single(signalp3_path: str, fasta_path: str, dest: str,
                          logger: logging.Logger, label: str):
    """运行单次SignalP 3.0|Run a single SignalP 3.0 invocation"""
    cmd = build_conda_command(signalp3_path, [
        '-t', 'euk', '-f', 'short',
        '-destination', dest, fasta_path,
    ])
    success, stdout, stderr = run_command(cmd, logger, label)
    return success, stdout, stderr


def _retry_signalp3_chunk(signalp3_path: str, chunk_path: str, logger: logging.Logger,
                           sprob_threshold: float, chunk_label: str,
                           tmp_base: str, retry_idx: int) -> Dict[str, Dict]:
    """失败chunk递归拆半重试|Recursively retry failed chunk by splitting in half"""
    import tempfile
    import uuid

    retry_id = uuid.uuid4().hex[:8]

    seqs = list(parse_fasta(chunk_path))
    if len(seqs) <= 1:
        logger.debug(f"SignalP 3.0 {chunk_label}: 单序列仍失败|single seq still failed, giving up")
        return {}

    mid = len(seqs) // 2
    sub_results = {}

    for part_idx, (start, end) in enumerate([(0, mid), (mid, len(seqs))]):
        sub_seqs = seqs[start:end]
        sub_path = os.path.join(tmp_base, f'{chunk_label}_{retry_id}_{retry_idx}_{part_idx}.fa')
        with open(sub_path, 'w') as f:
            for sid, seq in sub_seqs:
                f.write(f'>{sid}\n{seq}\n')

        sub_dest = os.path.join(tmp_base, f'{chunk_label}_{retry_id}_{retry_idx}_{part_idx}_dest')
        os.makedirs(sub_dest, exist_ok=True)

        success, stdout, stderr = _run_signalp3_single(
            signalp3_path, sub_path, sub_dest, logger,
            f"SignalP 3.0 {chunk_label} 重试|retry {part_idx+1}/2 ({len(sub_seqs)}条|seqs)"
        )

        if success and 'error running HOW' not in (stdout + stderr):
            sub_results.update(parse_signalp3_short(stdout, sprob_threshold))
        else:
            sub_results.update(_retry_signalp3_chunk(
                signalp3_path, sub_path, logger, sprob_threshold,
                f"{chunk_label}.{part_idx}", tmp_base, retry_idx + 1
            ))

    return sub_results


def run_signalp3(signalp3_path: str, fasta_path: str, logger: logging.Logger,
                 sprob_threshold: float = 0.9) -> Dict[str, Dict]:
    """运行SignalP 3.0(自动分chunk+失败递归重试)|Run SignalP 3.0 with chunking and recursive retry"""
    import tempfile

    seq_count = sum(1 for _ in parse_fasta(fasta_path))
    if seq_count <= 3800:
        chunks = [fasta_path]
        logger.info(f"SignalP 3.0: {seq_count}条序列(无需分chunk)|{seq_count} seqs, no chunking needed")
    else:
        chunks = split_fasta_chunks(fasta_path, chunk_size=300)
        logger.info(f"SignalP 3.0: {seq_count}条序列分割为{len(chunks)}个chunk(300条/chunk)|{seq_count} seqs split into {len(chunks)} chunks")

    all_results = {}
    chunk_failed = 0
    tmp_base = tempfile.mkdtemp(prefix='sp3_run_')

    for i, chunk in enumerate(chunks):
        dest = os.path.join(tmp_base, f'dest_{i:03d}')
        os.makedirs(dest, exist_ok=True)

        success, stdout, stderr = _run_signalp3_single(
            signalp3_path, chunk, dest, logger,
            f"SignalP 3.0 chunk {i+1}/{len(chunks)}"
        )

        if success and 'error running HOW' not in (stdout + stderr):
            all_results.update(parse_signalp3_short(stdout, sprob_threshold))
        else:
            chunk_failed += 1
            logger.warning(f"SignalP 3.0 chunk {i+1}/{len(chunks)} 失败，开始递归重试|chunk failed, starting recursive retry")
            retry_results = _retry_signalp3_chunk(
                signalp3_path, chunk, logger, sprob_threshold,
                f"chunk{i+1}", tmp_base, 0
            )
            recovered = len(retry_results)
            all_results.update(retry_results)
            logger.info(f"SignalP 3.0 chunk {i+1}/{len(chunks)} 重试完成|retry done: 恢复{recovered}条|recovered {recovered} seqs")

    shutil.rmtree(tmp_base, ignore_errors=True)

    if chunk_failed > 0:
        logger.info(f"SignalP 3.0: {chunk_failed}/{len(chunks)}个chunk失败后重试|chunks failed then retried")

    sp_count = sum(1 for v in all_results.values() if v['has_signal_peptide'])
    logger.info(
        f"SignalP 3.0完成|completed: {len(all_results)}/{seq_count}条预测|predictions, "
        f"{sp_count}条含信号肽|with SP (Sprob>={sprob_threshold})"
    )
    return all_results


def save_signalp3_compatible(results: Dict[str, Dict], output_path: str):
    """将SP3结果保存为SP6兼容格式|Save SP3 results in SP6-compatible format"""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# SignalP-3.0\tConverted to SP6-compatible format\n")
        f.write("# ID\tPrediction\tSP_score\n")
        for pid in sorted(results.keys()):
            info = results[pid]
            pred = "SP" if info['has_signal_peptide'] else "OTHER"
            score = info.get('sp_score', '0')
            f.write(f"{pid}\t{pred}\t{score}\n")


def merge_signalp_results(*result_dicts) -> Dict[str, Dict]:
    """合并多个SignalP结果(并集: 任一版本SP+即为SP+)|Merge SignalP results (union)"""
    merged = {}
    for results in result_dicts:
        for pid, info in results.items():
            if pid in merged:
                if info['has_signal_peptide'] and not merged[pid]['has_signal_peptide']:
                    merged[pid] = info
            else:
                merged[pid] = info
    return merged


def is_step_completed(output_file: str) -> bool:
    """检查步骤是否已完成|Check if step is completed"""
    return os.path.exists(output_file)


def parse_blastp_tabular(tabular_file: str) -> Set[str]:
    """解析BLASTP -outfmt 6输出，返回命中目标ID集合|Parse BLASTP -outfmt 6, return set of subject IDs"""
    hits = set()
    try:
        with open(tabular_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) >= 2:
                    hits.add(fields[1].split('/')[0])
    except Exception as e:
        raise IOError(f"解析BLASTP结果失败|Failed to parse BLASTP results: {e}")
    return hits


def parse_tmhmm_output(tmhmm_file: str) -> Dict[str, List[Tuple[int, int]]]:
    """解析TMHMM输出，返回蛋白ID到TM螺旋区间列表的映射|Parse TMHMM output, return protein ID -> list of TM helix (start, end)"""
    tm_helices = {}
    try:
        with open(tmhmm_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4 and fields[2] == 'TMhelix':
                    protein_id = fields[0].split('/')[0]
                    pos_fields = fields[3].split()
                    start = int(pos_fields[0])
                    end = int(pos_fields[1])
                    if protein_id not in tm_helices:
                        tm_helices[protein_id] = []
                    tm_helices[protein_id].append((start, end))
    except Exception as e:
        raise IOError(f"解析TMHMM结果失败|Failed to parse TMHMM results: {e}")
    return tm_helices


def has_tm_outside_sp(tm_helices: List[Tuple[int, int]], sp_cleavage_site: int) -> bool:
    """检查是否存在信号肽区外的跨膜域|Check if any TM helix exists outside signal peptide region"""
    if not tm_helices:
        return False
    for start, end in tm_helices:
        if start > sp_cleavage_site:
            return True
    return False


def generate_software_versions(output_dir: str, tools: Dict[str, str],
                               params: Dict, start_time=None):
    """生成software_versions.yml|Generate software_versions.yml"""
    try:
        import yaml
        from datetime import datetime
    except ImportError:
        return

    if start_time is None:
        start_time = datetime.now()

    versions = {}
    for name, path in tools.items():
        try:
            cmd = build_conda_command(path, ['--version'])
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            version = r.stdout.strip().split('\n')[0] if r.returncode == 0 else 'unknown'
        except Exception:
            version = 'unknown'
        versions[name] = {'version': version, 'path': path}

    end_time = datetime.now()
    info = {
        'pipeline': {
            'name': 'biopytools phyto_effector',
            'version': '1.0.0'
        },
        'tools': versions,
        'parameters': params,
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': int((end_time - start_time).total_seconds())
        }
    }

    info_dir = Path(output_dir) / '00_pipeline_info'
    info_dir.mkdir(parents=True, exist_ok=True)
    info_file = info_dir / 'software_versions.yml'
    with open(info_file, 'w', encoding='utf-8') as f:
        yaml.dump(info, f, default_flow_style=False, allow_unicode=True)


def merge_candidate_files(sample_dirs: List[str], effector_type: str,
                          merged_output_dir: str, logger=None) -> Optional[str]:
    """合并多个样本的候选效应子TSV文件|Merge candidate TSVs from multiple samples

    Args:
        sample_dirs: 各样本输出目录列表|List of per-sample output directories
        effector_type: 效应子类型|Effector type (rxlr/crn/nlp/protease/scp/elicitin/yxsl)
        merged_output_dir: 合并文件输出目录|Directory for merged output file
        logger: 日志记录器|Logger instance

    Returns:
        合并文件路径，无数据时返回None|Merged file path, or None if no data
    """
    import pandas as pd

    if effector_type == 'rxlr':
        candidates_rel = os.path.join('06_candidates', 'rxlr_candidates.tsv')
    elif effector_type == 'crn':
        candidates_rel = os.path.join('03_candidates', 'crn_candidates.tsv')
    else:
        from .generic_finder import EFFECTOR_TYPE_CONFIG
        tc = EFFECTOR_TYPE_CONFIG.get(effector_type, {})
        candidates_rel = os.path.join(tc.get('output_dir', ''), tc.get('candidates_file', ''))

    all_dfs = []
    for sample_dir in sample_dirs:
        tsv_path = os.path.join(sample_dir, effector_type, candidates_rel)
        if os.path.exists(tsv_path):
            try:
                df = pd.read_csv(tsv_path, sep='\t')
                if not df.empty:
                    all_dfs.append(df)
            except Exception as e:
                if logger:
                    logger.warning(f"读取候选文件失败|Failed to read candidate file: {tsv_path}: {e}")
        elif logger:
            logger.warning(f"未找到候选文件|Candidate file not found: {tsv_path}")

    if not all_dfs:
        if logger:
            logger.warning(f"无{effector_type}候选结果可合并|No {effector_type} candidates to merge")
        return None

    merged_df = pd.concat(all_dfs, ignore_index=True)
    os.makedirs(merged_output_dir, exist_ok=True)
    merged_file = os.path.join(merged_output_dir, f'{effector_type}_candidates_all_samples.tsv')
    merged_df.to_csv(merged_file, sep='\t', index=False)

    if logger:
        logger.info(
            f"合并{effector_type}候选|Merged {effector_type} candidates: "
            f"{len(all_dfs)}个样本|samples, {len(merged_df)}条记录|records -> {merged_file}"
        )
    return merged_file
