"""
Phobius工具函数模块|Phobius utility functions module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令所在conda环境,返回环境名|Detect conda env for a command, return env name

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名或None|conda env name or None
    """
    # 优先从完整路径检测|First detect from full path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r"/envs/([^/]+)", cmd_path)
        if match:
            return match.group(1)

    # 兜底: 搜索所有conda环境|Fallback: search all conda envs
    conda_base = os.environ.get("CONDA_EXE")
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, "envs")
        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, "bin", os.path.basename(command))
                if os.path.exists(env_bin):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令(列表形式)|Build conda run command (list form)

    必须传完整路径(含/envs/), 不能用basename|Must pass full path (with /envs/), not basename

    Args:
        command: 命令完整路径|Full command path
        args: 命令参数|Command arguments

    Returns:
        完整命令列表(配合subprocess.run(shell=False))|Full command list (for shell=False)
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ["conda", "run", "-n", conda_env, "--no-capture-output", command] + args
    return [command] + args


class PhobiusLogger:
    """Phobius日志管理器|Phobius logger manager"""

    def __init__(self, output_path, log_name: str = "phobius.log"):
        self.output_path = Path(str(output_path))
        logs_dir = self.output_path / "99_logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = logs_dir / log_name
        self._setup_logging()

    def _setup_logging(self):
        """设置日志(stdout INFO + stderr WARNING+ + 文件 DEBUG)|Setup logging"""
        self.logger = logging.getLogger("Phobius")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        file_handler = logging.FileHandler(self.log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger|Get logger"""
        return self.logger


def run_phobius(logger, phobius_path: str, mode: str,
                input_file: str, output_file: str) -> bool:
    """
    执行phobius预测|Run phobius prediction

    Args:
        logger: logger对象|logger object
        phobius_path: phobius.pl完整路径|full path to phobius.pl
        mode: '-short' 或 '-long'|'-short' or '-long'
        input_file: 输入蛋白质FASTA|input protein FASTA
        output_file: 输出文件路径|output file path

    Returns:
        是否成功|whether successful
    """
    cmd = build_conda_command(phobius_path, [mode, input_file])
    logger.info(f"执行|Executing: phobius {mode}")
    logger.info(f"命令|Command: {' '.join(cmd)} > {output_file}")

    try:
        with open(output_file, "w") as fh:
            result = subprocess.run(
                cmd,
                stdout=fh,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        if result.stderr:
            logger.debug(f"phobius stderr: {result.stderr.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"phobius执行失败|phobius failed: {(e.stderr or '').strip()}")
        return False
    except FileNotFoundError:
        logger.error(f"phobius不存在|phobius not found: {phobius_path}")
        return False


def parse_short_output(short_file: str) -> Dict[str, dict]:
    """
    解析phobius -short输出|Parse phobius -short output

    跳过banner(前3行)与列头行,逐数据行解析|Skip banner + header, parse data lines

    Returns:
        {id: {"tm_short": int, "sp_raw": str, "prediction": str}}
    """
    result: Dict[str, dict] = {}
    with open(short_file, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # 跳过banner与列头|Skip banner and header
            if stripped.startswith("Phobius ver") or stripped.startswith("(c)"):
                continue
            if stripped.startswith("SEQENCE"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            try:
                tm = int(parts[1])
            except ValueError:
                continue
            pid = parts[0]
            sp_raw = parts[2]
            prediction = " ".join(parts[3:]) if len(parts) > 3 else ""
            result[pid] = {"tm_short": tm, "sp_raw": sp_raw, "prediction": prediction}
    return result


def parse_long_output(long_file: str) -> Dict[str, dict]:
    """
    解析phobius -long FT特征表|Parse phobius -long FT feature table

    从FT SIGNAL/TRANSMEM提取坐标|Extract coords from FT SIGNAL/TRANSMEM

    Returns:
        {id: {"tm": int, "sp": "Y"/"N", "sp_region": str, "tm_regions": List[str]}}
    """
    result: Dict[str, dict] = {}
    current_id = None
    tm_regions: List[str] = []
    sp_region = None

    def _flush():
        # 结束当前记录,写入result|Finish current record, store in result
        nonlocal current_id, tm_regions, sp_region
        if current_id is not None:
            result[current_id] = {
                "tm": len(tm_regions),
                "sp": "Y" if sp_region else "N",
                "sp_region": sp_region if sp_region else "-",
                "tm_regions": list(tm_regions),
            }
        current_id = None
        tm_regions = []
        sp_region = None

    with open(long_file, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("ID "):
                _flush()
                current_id = stripped.split()[1]
            elif stripped.startswith("FT "):
                fields = stripped.split()
                # FT  SIGNAL  start end  /  FT  TRANSMEM  start end
                if len(fields) >= 4 and fields[1] == "SIGNAL":
                    sp_region = f"{fields[2]}-{fields[3]}"
                elif len(fields) >= 4 and fields[1] == "TRANSMEM":
                    tm_regions.append(f"{fields[2]}-{fields[3]}")
            elif stripped == "//":
                _flush()
    _flush()  # 兜底处理文件末尾|flush trailing record
    return result


def build_tsv_records(short_map: Dict[str, dict], long_map: Dict[str, dict]) -> List[dict]:
    """
    合并short/long两路为TSV记录|Merge short/long into TSV records

    TM/SP/坐标以long为准,Prediction取short;ID顺序跟随short_map|TM/SP/coords from long,
    Prediction from short; order follows short_map

    Returns:
        记录列表|list of record dicts
    """
    records: List[dict] = []
    seen = set()
    # 先按short顺序|First in short order
    for pid, short_rec in short_map.items():
        seen.add(pid)
        long_rec = long_map.get(pid, {})
        tm = long_rec.get("tm", short_rec.get("tm_short", 0))
        tm_regions_list = long_rec.get("tm_regions", [])
        records.append({
            "id": pid,
            "tm": tm,
            "sp": long_rec.get("sp", "N"),
            "sp_region": long_rec.get("sp_region", "-"),
            "tm_regions": ";".join(tm_regions_list) if tm_regions_list else "-",
            "prediction": short_rec.get("prediction", ""),
        })
    # 补long-only的ID(缺short)|Add long-only ids (missing short)
    for pid, long_rec in long_map.items():
        if pid in seen:
            continue
        tm_regions_list = long_rec.get("tm_regions", [])
        records.append({
            "id": pid,
            "tm": long_rec.get("tm", 0),
            "sp": long_rec.get("sp", "N"),
            "sp_region": long_rec.get("sp_region", "-"),
            "tm_regions": ";".join(tm_regions_list) if tm_regions_list else "-",
            "prediction": "",
        })
    return records


def write_clean_tsv(records: List[dict], output_file: str) -> None:
    """
    写出解析后的TSV|Write cleaned TSV

    表头: ID  TM  SP  SP_region  TM_regions  Prediction
    """
    header = ["ID", "TM", "SP", "SP_region", "TM_regions", "Prediction"]
    with open(output_file, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in records:
            fh.write("\t".join([
                r["id"],
                str(r["tm"]),
                r["sp"],
                r["sp_region"],
                r["tm_regions"],
                r["prediction"],
            ]) + "\n")
