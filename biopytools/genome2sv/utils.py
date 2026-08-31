"""genome2sv 工具函数|genome2sv utilities

日志管理器、conda 环境包装、SV 统计纯函数。
|Logger, conda wrapping, and pure SV-stat helpers.
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

# SURVIVOR 统计的 SV 类型|SV types tracked for SURVIVOR stats
_VALID_SVTYPES = {"INS", "DEL", "INV", "DUP", "BND"}
_SUMMARY_COLS = ["sample", "total", "INS", "DEL", "INV", "DUP", "BND", "OTHER"]


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("genome2sv")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)   # WARNING+ → 超算 .err
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:                              # 全级别 → 文件
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """检测命令所在 conda 环境|Detect conda env of a command.

    先从 which 路径的 /envs/<name>/ 解析,否则搜索所有 conda 环境兜底。
    |Parse /envs/<name>/ from which path, else search all envs.
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        m = re.search(r"/envs/([^/]+)/", cmd_path)
        if m:
            return m.group(1)
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        base = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(base, "envs")
        if os.path.isdir(envs_dir):
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, "bin", command)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """构建 conda run 命令(必带 --no-capture-output)|Build conda run command."""
    env = get_conda_env(command)
    if env:
        return ["conda", "run", "-n", env, "--no-capture-output", command] + args
    return [command] + args


def check_dependencies(config, logger: logging.Logger) -> bool:
    """检查关键工具可用性(版本/可执行性探测)|Check key tools via version/exec probe.

    SURVIVOR 为子命令式 CLI 且无 --version:空参必非0退出(打印用法),故按"能执行"判定
    而非 rc==0,避免误报缺失终止流程。|SURVIVOR is a subcommand CLI with no --version;
    no-args exits non-zero (prints usage), so judge by "runs" not rc==0 to avoid
    false-missing aborts.
    """
    logger.info("检查依赖工具|Checking dependencies")
    # (path, args, expect_rc0): True=需 rc==0(有 --version);False=能执行即可(子命令式)|
    # (path, args, expect_rc0): True=need rc==0 (has --version); False=runs is enough
    tools = {
        "minimap2": (config.minimap2_path, ["--version"], True),
        "samtools": (config.samtools_path, ["--version"], True),
        "svim-asm": (config.svim_asm_path, ["--version"], True),
        "survivor": (config.survivor_path, [], False),
    }
    missing = []
    for name, (path, args, expect_rc0) in tools.items():
        try:
            cmd = build_conda_command(path, args)
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            available = (r.returncode == 0) if expect_rc0 else True
            if available:
                detail = (r.stdout or r.stderr).strip().splitlines()
                ver = detail[0][:60] if detail else "(子命令式CLI|subcommand CLI)"
                logger.info(f"{name} 可用|available: {ver}")
            else:
                missing.append(name)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            missing.append(name)
    if missing:
        logger.error(
            f"缺失工具|Missing tools: {', '.join(missing)} "
            f"(align 域环境或环境变量覆盖|align domain env or env-var overrides)")
        return False
    return True


def parse_svtype_stats(vcf_path: str) -> dict:
    """解析 VCF 统计 SVTYPE 分布|Parse VCF for SVTYPE distribution.

    Returns:
        dict: keys total/INS/DEL/INV/DUP/BND/OTHER;文件不存在返回全零。
        |dict with keys total/INS/DEL/INV/DUP/BND/OTHER; zeros if file missing.
    """
    counts = {k: 0 for k in _SUMMARY_COLS[1:]}
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                svtype = None
                for entry in fields[7].split(";"):
                    if entry.startswith("SVTYPE="):
                        svtype = entry.split("=", 1)[1]
                        break
                counts["total"] += 1
                if svtype in _VALID_SVTYPES:
                    counts[svtype] += 1
                else:
                    counts["OTHER"] += 1
    except FileNotFoundError:
        pass
    return counts


def build_survivor_input(sample_vcf_map: dict, output_file: str) -> List[str]:
    """写 SURVIVOR 输入文件(每行一个 VCF 绝对路径)|Write SURVIVOR input file.

    Args:
        sample_vcf_map: {sample_name: vcf_path} 仅含成功样本|successful samples only
        output_file: 输出文件路径|output path
    Returns:
        写入的绝对路径列表(按样本名排序)|list of absolute paths written (sorted)
    """
    paths = []
    with open(output_file, "w") as fh:
        for sample in sorted(sample_vcf_map):
            vcf = os.path.abspath(sample_vcf_map[sample])
            fh.write(vcf + "\n")
            paths.append(vcf)
    return paths


def format_sv_summary_tsv(rows: List[Tuple[str, dict]]) -> str:
    """格式化 SV 统计为 TSV|Format SV stats as TSV.

    Args:
        rows: [(name, counts_dict), ...]
    Returns:
        TSV 文本(表头 + 数据行)|TSV text with header and rows
    """
    lines = ["\t".join(_SUMMARY_COLS)]
    for name, counts in rows:
        cells = [name if c == "sample" else str(counts.get(c, 0))
                 for c in _SUMMARY_COLS]
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


# ---------- SV 序列提取 + PAV 矩阵|SV sequence extraction & PAV matrix ----------

_REVCOMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    """反转互补|Reverse complement."""
    return seq.translate(_REVCOMP_TABLE)[::-1]


class FaidxReader:
    """基于 .fai 的 FASTA 随机读取器|.fai-backed random-access FASTA reader.

    步骤 0 已为参考建好 .fai;此处纯 Python seek,避免逐条调 samtools faidx。
    |The reference is faidx-ed in step 0; plain Python seek here avoids one
    samtools faidx invocation per SV.
    """

    def __init__(self, fasta_path: str):
        fai_path = fasta_path + ".fai"
        if not os.path.exists(fai_path):
            raise FileNotFoundError(f"fai 索引不存在|fai missing: {fai_path}")
        self._fh = open(fasta_path, "rb")
        self._index = {}  # name -> (length, offset, linebases, linewidth)
        with open(fai_path) as fai:
            for line in fai:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                name, length, offset, linebases, linewidth = parts[:5]
                self._index[name] = (int(length), int(offset),
                                     int(linebases), int(linewidth))

    def fetch(self, chrom: str, start1: int, end1: int) -> str:
        """取 1-based 闭区间序列(end 超长自动截断)|Fetch 1-based closed interval (clamped)."""
        if chrom not in self._index:
            raise KeyError(f"染色体不在 fai|chrom missing from fai: {chrom}")
        length, offset, linebases, linewidth = self._index[chrom]
        start0 = max(0, start1 - 1)
        end0 = min(end1, length)
        if end0 <= start0:
            return ""
        # 字节区间:起点行内偏移 + 起始整行偏移;终点同理(含)|byte range via line math
        byte_start = offset + (start0 // linebases) * linewidth + (start0 % linebases)
        last0 = end0 - 1
        byte_end = offset + (last0 // linebases) * linewidth + (last0 % linebases) + 1
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start)
        return raw.decode().replace("\n", "").replace("\r", "").upper()

    def close(self) -> None:
        """关闭句柄|Close handle."""
        self._fh.close()


def gt_present(gt_field: str) -> int:
    """GT 字段转 PAV(含 allele 1 → 1,缺失/纯参考 → 0)|GT to PAV 0/1."""
    gt = gt_field.split(":", 1)[0]
    alleles = gt.replace("|", "/").split("/")
    return 1 if any(a == "1" for a in alleles) else 0


def stable_sv_id(svtype: str, number: int) -> str:
    """自增稳定 SV id(DUP 子类型归一)|Stable auto-increment SV id (DUP normalized).

    SURVIVOR 合并后原 VCF ID 可能重复,PAV 矩阵与序列 FASTA 共用此 id 互相对应。
    |SURVIVOR-merged IDs may collide; PAV matrix and sequence FASTA share this id.
    """
    base = svtype.split(":", 1)[0]
    return f"pan_sv.{base}.{number:05d}"


def parse_info_str(info: str) -> dict:
    """INFO 字符串转 dict|INFO string to dict."""
    out = {}
    for entry in info.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            out[key] = value
    return out


def _is_symbolic(allele: str) -> bool:
    """判断等位基因是否符号化/占位|Allele symbolic or placeholder."""
    return allele.startswith("<") or allele == "N" or allele == ""


def extract_sv_sequence(svtype: str, chrom: str, pos: int, ref: str, alt: str,
                        info: dict, reader: Optional["FaidxReader"]
                        ) -> Optional[Tuple[str, str]]:
    """按 SV 类型提取序列|Extract sequence for one SV by type.

    策略(与 svim-asm sequence-alleles 输出对齐)|Strategy:
      INS → ALT 剥 anchor;DEL → REF 字段;INV → ALT(已是 revcomp);
      DUP* → 参考 [POS,END] 重复单元;符号化时 DEL/INV/DUP 回退坐标提取。
      |INS→ALT minus anchor; DEL→REF field; INV→ALT (already revcomp);
      DUP*→reference [POS,END] unit; symbolic alleles fall back to region.

    Returns:
        (sequence, source) 或 None(无法提取:BND/未知/符号化 INS)|None if not extractable
    """
    try:
        end = int(info.get("END", pos))
    except (TypeError, ValueError):
        end = pos   # 畸形 END(逗号列表/非数字)回退起点,序列提取按单碱基跳过
    base_type = svtype.split(":", 1)[0]

    def from_region(invert: bool = False) -> Optional[Tuple[str, str]]:
        if reader is None:
            return None
        try:
            seq = reader.fetch(chrom, pos, end)
        except (KeyError, OSError, ValueError):
            return None
        if not seq:
            return None
        return (revcomp(seq), "region_revcomp") if invert else (seq, "region")

    if base_type == "INS":
        if _is_symbolic(alt) or len(alt) <= len(ref):
            return None  # 参考中没有插入序列,无法回退|no reference source for INS
        return alt[len(ref):].upper(), "alt"
    if base_type == "DEL":
        if not _is_symbolic(ref) and len(ref) > 1:
            return ref.upper(), "ref"
        return from_region()
    if base_type == "INV":
        if not _is_symbolic(alt):
            return alt.upper(), "alt"
        return from_region(invert=True)
    if base_type == "DUP":
        return from_region()
    return None  # BND/未知类型无区间|BND/unknown have no interval


def format_pav_matrix(rows: List[Tuple], sample_names: List[str]) -> str:
    """PAV 主矩阵 TSV(带元数据列)|PAV matrix TSV with metadata columns.

    Args:
        rows: [(sv_id, chrom, pos, end, svtype, svlen, [0/1,...]), ...]
    """
    header = ["sv_id", "chrom", "pos", "end", "svtype", "svlen"] + sample_names
    lines = ["\t".join(header)]
    for sv_id, chrom, pos, end, svtype, svlen, pav in rows:
        cells = [sv_id, chrom, str(pos), str(end), svtype, str(svlen)]
        cells += [str(p) for p in pav]
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


def format_pav_binary(rows: List[Tuple], sample_names: List[str]) -> str:
    """纯 0/1 PAV 矩阵 TSV(R 可直接 as.matrix)|Pure 0/1 PAV matrix TSV."""
    lines = ["\t".join(["sv_id"] + sample_names)]
    for row in rows:
        lines.append("\t".join([row[0]] + [str(p) for p in row[-1]]))
    return "\n".join(lines) + "\n"


def write_software_versions(config, logger: logging.Logger, output_path: str,
                            start_time=None) -> None:
    """生成 software_versions.yml|Generate software_versions.yml.

    探测 6 个工具版本 + 记录参数与运行时间,写入 00_pipeline_info/software_versions.yml。
    |Probe 6 tool versions, record parameters and runtime; write yml.
    """
    from datetime import datetime
    import yaml
    tools = {
        "minimap2": (config.minimap2_path, ["--version"]),
        "samtools": (config.samtools_path, ["--version"]),
        "svim-asm": (config.svim_asm_path, ["--version"]),
        "bcftools": (config.bcftools_path, ["--version"]),
        "bedtools": (config.bedtools_path, ["--version"]),
        "survivor": (config.survivor_path, []),
    }
    versions = {}
    for name, (path, args) in tools.items():
        try:
            cmd = build_conda_command(path, args)
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            raw = (r.stdout or r.stderr).strip()
            ver = raw.splitlines()[0] if raw else "unknown"
            versions[name] = {"version": ver, "path": path}
        except Exception as e:
            logger.warning(f"版本探测失败|Version probe failed [{name}]: {e}")
            versions[name] = {"version": "unknown", "path": path}
    # 关键参数(getattr 容错,适配测试桩)|key params (getattr-tolerant for stubs)
    param_keys = ["ref_sample", "preset", "svim_mode", "threads", "max_dist",
                  "min_support", "survivor_type", "survivor_strand", "est_dist",
                  "min_sv_length", "svim_min_sv_size"]
    info = {
        "pipeline": {"name": "biopytools genome2sv", "version": "1.1.0"},
        "tools": versions,
        "parameters": {k: getattr(config, k, None) for k in param_keys},
    }
    end_time = datetime.now()
    if start_time is not None:
        info["execution"] = {
            "start_time": start_time.strftime("%Y-%m-%d %H:%M:%S"),
            "end_time": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "runtime_seconds": int((end_time - start_time).total_seconds()),
        }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
