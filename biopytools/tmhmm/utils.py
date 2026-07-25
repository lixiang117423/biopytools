"""
TMHMM工具类|TMHMM Utilities
"""

import logging
import sys
import os
import subprocess


class TmhmmLogger:
    """TMHMM日志管理器|TMHMM Logger Manager"""

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        # 日志统一放99_logs子目录(§12.2.3)|Logs live under 99_logs/
        self.log_dir = os.path.join(output_dir, "99_logs")
        os.makedirs(self.log_dir, exist_ok=True)

        log_file = os.path.join(self.log_dir, "tmhmm.log")
        self._setup_logging(log_file)

    def _setup_logging(self, log_file: str):
        self.logger = logging.getLogger("Tmhmm")
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

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        return self.logger


def run_tmhmm(
    logger: logging.Logger,
    tmhmm_path: str,
    input_file: str,
    output_file: str,
    noplot: bool = True,
) -> bool:
    """执行TMHMM预测|Execute TMHMM prediction

    tmhmm直接调用(不用conda run包装):其Perl驱动仅用核心模块Getopt::Long,
    内部decodeanhmm为静态链接二进制,不依赖conda环境的动态库。
    |tmhmm is called directly (no conda run): its Perl driver uses only the
    core Getopt::Long module and the internal decodeanhmm is a statically
    linked binary, so it does not depend on the conda env's shared libs.

    Args:
        logger: 日志对象|Logger object
        tmhmm_path: tmhmm可执行文件路径|tmhmm binary path
        input_file: 输入蛋白质FASTA文件|Input protein FASTA file
        output_file: 输出文件路径|Output file path
        noplot: 不生成图形|No plot generation

    Returns:
        是否成功|Whether successful
    """
    cmd = [tmhmm_path]
    if noplot:
        cmd.append('-noplot')
    cmd.append(input_file)

    # 以output_dir为cwd: tmhmm的plot文件写到cwd($opt_workdir="."),否则散落到调用方CWD
    # |cwd=output_dir: tmhmm writes plot files to cwd, else they scatter in the caller's CWD
    cwd = os.path.dirname(os.path.abspath(output_file))

    cmd_display = ' '.join(cmd) + f' > {output_file}'
    logger.info(f"命令|Command: {cmd_display}")

    # 原子写: 先写.tmp,成功才replace,失败清理半成品(避免断点续传解析截断文件,§10.2)
    # |Atomic write: write .tmp, replace on success, clean up on failure
    tmp_file = output_file + '.tmp'
    try:
        with open(tmp_file, 'w') as fh:
            result = subprocess.run(
                cmd,
                stdout=fh,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                cwd=cwd,
            )
        if result.stderr:
            logger.debug(f"tmhmm stderr|tmhmm stderr: {result.stderr.strip()}")
        os.replace(tmp_file, output_file)
        return True
    except subprocess.CalledProcessError as e:
        _safe_remove(tmp_file)
        logger.error(f"tmhmm执行失败|tmhmm execution failed: {(e.stderr or '').strip()}")
        return False
    except FileNotFoundError:
        _safe_remove(tmp_file)
        logger.error(f"tmhmm不存在|tmhmm not found: {tmhmm_path}")
        return False


def _safe_remove(path: str):
    """安全删除文件(忽略不存在)|Safely remove a file (ignore missing)"""
    try:
        os.remove(path)
    except OSError:
        pass


def generate_software_version_yml(config, tool_versions=None) -> str:
    """生成software_versions.yml到00_pipeline_info(§12.5)|Generate software_versions.yml"""
    import yaml
    from datetime import datetime

    info_dir = os.path.join(config.output_dir, "00_pipeline_info")
    os.makedirs(info_dir, exist_ok=True)
    output_file = os.path.join(info_dir, "software_versions.yml")

    if tool_versions is None:
        tool_versions = {}
        try:
            result = subprocess.run(
                [config.tmhmm_path], capture_output=True, text=True, timeout=30)
            ver = (result.stdout or result.stderr or "").strip() or "unknown"
            tool_versions["tmhmm"] = {"version": ver, "path": config.tmhmm_path}
        except Exception:
            tool_versions["tmhmm"] = {"version": "unknown", "path": config.tmhmm_path}

    info = {
        "pipeline": {"name": "biopytools tmhmm", "version": "1.0.0"},
        "tools": tool_versions,
        "parameters": {
            "noplot": config.noplot,
            "output_prefix": config.output_prefix,
        },
        "execution": {
            "start_time": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        },
    }

    with open(output_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(info, f, default_flow_style=False, sort_keys=False)
    return output_file


def parse_tmhmm_output(raw_file: str) -> list:
    """解析tmhmm原始输出|Parse tmhmm raw output

    Args:
        raw_file: tmhmm原始输出文件路径|tmhmm raw output file path

    Returns:
        解析后的记录列表|Parsed record list
    """
    import re

    records = []
    current = {}

    with open(raw_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('# '):
                m_id = re.match(r'^# (\S+) Length:\s*(\d+)', line)
                if m_id:
                    if current:
                        records.append(current)
                    current = {
                        'id': m_id.group(1),
                        'length': int(m_id.group(2)),
                        'pred_tmhs': 0,
                        'exp_aas_in_tmhs': 0.0,
                        'exp_aas_first60': 0.0,
                        'prob_n_in': 0.0,
                        'signal_seq': '',
                    }
                    continue

                m = re.match(r'^# \S+ Number of predicted TMHs:\s*(\d+)', line)
                if m and current:
                    current['pred_tmhs'] = int(m.group(1))
                    continue

                m = re.match(r'^# \S+ Exp number of AAs in TMHs:\s*([\d.]+)', line)
                if m and current:
                    current['exp_aas_in_tmhs'] = float(m.group(1))
                    continue

                m = re.match(r'^# \S+ Exp number, first 60 AAs:\s*([\d.]+)', line)
                if m and current:
                    current['exp_aas_first60'] = float(m.group(1))
                    continue

                m = re.match(r'^# \S+ Total prob of N-in:\s*([\d.]+)', line)
                if m and current:
                    current['prob_n_in'] = float(m.group(1))
                    continue

                if 'POSSIBLE N-term signal sequence' in line and current:
                    current['signal_seq'] = 'yes'

    if current:
        records.append(current)

    return records


def write_clean_tsv(records: list, output_file: str):
    """写入整理后的TSV文件|Write cleaned TSV file

    Args:
        records: 解析后的记录列表|Parsed record list
        output_file: 输出文件路径|Output file path
    """
    header = ['ID', 'Length', 'Pred_TMHs', 'Exp_AAs_in_TMHs', 'Exp_AAs_first60', 'Prob_N_in', 'Signal_Seq']
    with open(output_file, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in records:
            row = [
                r['id'],
                str(r['length']),
                str(r['pred_tmhs']),
                f"{r['exp_aas_in_tmhs']:.2f}",
                f"{r['exp_aas_first60']:.2f}",
                f"{r['prob_n_in']:.5f}",
                r['signal_seq'],
            ]
            fh.write('\t'.join(row) + '\n')
