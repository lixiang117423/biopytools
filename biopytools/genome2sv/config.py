"""genome2sv 配置与 fof 解析|genome2sv configuration & fof parsing"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

from ..common.paths import get_tool_path, expand_path


def parse_fof(fof_path: str) -> List[Tuple[str, str]]:
    """解析样本清单(name<TAB>path)|Parse sample manifest.

    每行: name<TAB>path;`#` 注释行与空行忽略;非注释行缺 TAB 视为格式错误抛 ValueError。
    |Each line: name<TAB>path; `#` comments and blank lines ignored; a
    non-comment line without TAB raises ValueError.

    Returns:
        [(name, path), ...] 保持文件顺序|entries in file order
    """
    entries: List[Tuple[str, str]] = []
    with open(fof_path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" not in line:
                raise ValueError(
                    f"fof 第 {lineno} 行格式错误(需 name<TAB>path)|"
                    f"malformed fof line {lineno} (need name<TAB>path): {line!r}")
            parts = line.split("\t")
            # 展开 fof 内路径的 ~ 与环境变量,否则 validate() 的存在性检查会误报
            # |Expand ~ and env vars in fof paths or validate() existence check misfires
            name, path = parts[0].strip(), expand_path(parts[1].strip())
            if name and path:
                entries.append((name, path))
    return entries


@dataclass
class Genome2SVConfig:
    """assembly-to-assembly SV calling 配置|Assembly-to-assembly SV calling config"""

    # 必需|Required
    input_fof: str       # 样本清单(name<TAB>path)|sample manifest
    ref_sample: str      # 参考样本名(须在 fof 中)|reference sample name

    # 输出|Output
    output_dir: str = "./genome2sv_output"

    # 运行|Runtime
    threads: int = 12
    preset: str = "asm5"          # minimap2 预设|minimap2 preset (asm5/asm10/asm20)
    svim_mode: str = "haploid"    # svim-asm 模式|svim-asm mode (haploid/diploid)
    log_level: str = "INFO"

    # SURVIVOR merge 参数(v1.0.7, 7 参)|SURVIVOR merge params (v1.0.7)
    max_dist: int = 1000          # 断点最大距离(bp)|max breakpoint distance
    min_support: int = 1          # 最小支持调用数|min supporting callers
    survivor_type: int = 1        # SV 类型一致(1)/任意(0)|same type (1) or any (0)
    survivor_strand: int = 1      # 链方向一致(1)/任意(0)|same strand (1) or any (0)
    est_dist: int = 1             # 按长度估距(1=yes)|estimate dist by size
    min_sv_length: int = 50       # 最小 SV 长度(bp)|min SV length

    # svim-asm 额外参数|svim-asm extra param
    svim_min_sv_size: int = 40    # svim-asm --min_sv_size

    # 侧翼序列(步骤6)|flank sequences (step 6)
    flank: int = 300              # SV 上下游侧翼长度(bp)|up/downstream flank bp

    # 工具路径:固定 align 域环境全路径,环境变量/用户配置可覆盖(§13.2.3 传全路径)
    # 裸命令名会被 get_conda_env 按 PATH/listdir 顺序随机解析(曾漂移到 mga/
    # Augustus),管道混 env 后只能靠父进程 PATH 侥幸命中 minimap2。
    # 六工具已全部并入 align 域环境(svim-asm/SURVIVOR 2026-08-31 自 sv_calling
    # 补齐,见 envs/align.yml 与 envs/README.md 踩坑实录)。
    # |Tool paths: pinned to the align domain env (env var/user config may
    # override). Bare names resolved nondeterministically (drifted to mga/
    # Augustus envs). All six tools merged into align (svim-asm/SURVIVOR
    # added from sv_calling on 2026-08-31).
    minimap2_path: str = field(default_factory=lambda: get_tool_path(
        "minimap2", "~/miniforge3/envs/align/bin/minimap2", "MINIMAP2_PATH"))
    samtools_path: str = field(default_factory=lambda: get_tool_path(
        "samtools", "~/miniforge3/envs/align/bin/samtools", "SAMTOOLS_PATH"))
    svim_asm_path: str = field(default_factory=lambda: get_tool_path(
        "svim-asm", "~/miniforge3/envs/align/bin/svim-asm", "SVIM_ASM_PATH"))
    bcftools_path: str = field(default_factory=lambda: get_tool_path(
        "bcftools", "~/miniforge3/envs/align/bin/bcftools", "BCFTOOLS_PATH"))
    bedtools_path: str = field(default_factory=lambda: get_tool_path(
        "bedtools", "~/miniforge3/envs/align/bin/bedtools", "BEDTOOLS_PATH"))
    survivor_path: str = field(default_factory=lambda: get_tool_path(
        "survivor", "~/miniforge3/envs/align/bin/SURVIVOR", "SURVIVOR_PATH"))

    def __post_init__(self):
        """展开路径、建目录、解析 fof|Expand paths, make dirs, parse fof"""
        self.input_fof = expand_path(self.input_fof)
        self.output_dir = expand_path(self.output_dir)
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 步骤目录|step directories
        self.reference_dir = self.output_path / "reference"
        self.alignment_dir = self.output_path / "01_alignment"
        self.svim_dir = self.output_path / "02_svim"
        self.merged_dir = self.output_path / "03_merged"
        self.stats_dir = self.output_path / "04_stats"
        self.sv_seq_dir = self.output_path / "05_sv_sequences"
        self.flank_dir = self.output_path / "06_sv_flanks"
        self.logs_dir = self.output_path / "99_logs"
        self.info_dir = self.output_path / "00_pipeline_info"
        for d in (self.reference_dir, self.alignment_dir, self.svim_dir,
                  self.merged_dir, self.stats_dir, self.sv_seq_dir,
                  self.flank_dir, self.logs_dir, self.info_dir):
            d.mkdir(parents=True, exist_ok=True)

        # 解析 fof(文件不存在则留空,交 validate 报错)|parse fof
        self._fof_entries: List[Tuple[str, str]] = []
        self.samples: Dict[str, str] = {}
        if Path(self.input_fof).is_file():
            self._fof_entries = parse_fof(self.input_fof)
            self.samples = dict(self._fof_entries)
        self.reference_path: str = self.samples.get(self.ref_sample, "")

    @property
    def reference_fasta(self) -> str:
        """reference/ 下参考 fasta 路径|reference fasta path under reference/"""
        if not self.reference_path:
            return ""
        return str(self.reference_dir / Path(self.reference_path).name)

    @property
    def reference_fai(self) -> str:
        """参考 fai 路径|reference fai path"""
        return self.reference_fasta + ".fai" if self.reference_path else ""

    def validate(self) -> bool:
        """校验配置|Validate; raise ValueError on invalid"""
        errors = []
        if not Path(self.input_fof).is_file():
            errors.append(f"fof 文件不存在|fof not found: {self.input_fof}")
        if not self.ref_sample:
            errors.append("参考样本名为空|ref sample name empty")
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Threads must be positive: {self.threads}")
        if self.preset not in ("asm5", "asm10", "asm20"):
            errors.append(f"无效 preset|Invalid preset: {self.preset} (asm5/asm10/asm20)")
        if self.svim_mode not in ("haploid", "diploid"):
            errors.append(f"无效 svim 模式|Invalid svim mode: {self.svim_mode} (haploid/diploid)")
        if self.survivor_type not in (0, 1):
            errors.append(f"survivor_type 必须为 0 或 1|must be 0 or 1: {self.survivor_type}")
        if self.survivor_strand not in (0, 1):
            errors.append(f"survivor_strand 必须为 0 或 1|must be 0 or 1: {self.survivor_strand}")
        if self.max_dist < 0 or self.min_sv_length < 0:
            errors.append("max_dist/min_sv_length 不能为负|cannot be negative")
        if self.flank < 0:
            errors.append(f"flank 不能为负|flank cannot be negative: {self.flank}")
        # 仅当 fof 可解析且参考名非空时做样本校验|sample checks only if parseable
        if Path(self.input_fof).is_file() and self.ref_sample:
            names = [n for n, _ in self._fof_entries]
            if len(names) != len(set(names)):
                errors.append("fof 存在重名样本|duplicate sample names in fof")
            if self.ref_sample not in self.samples:
                errors.append(f"参考样本不在 fof 中|ref sample not in fof: {self.ref_sample}")
            elif not any(n != self.ref_sample for n in self.samples):
                errors.append("fof 中除参考外无查询样本|no query samples besides reference")
            for name, path in self._fof_entries:
                if not Path(path).exists():
                    errors.append(f"样本路径不存在|sample path missing [{name}]: {path}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
