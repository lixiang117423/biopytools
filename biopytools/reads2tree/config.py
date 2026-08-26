"""reads2tree 配置|reads2tree configuration"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

from ..common.paths import expand_path, get_tool_path
from .utils import scan_fastq_dir


@dataclass
class Reads2TreeConfig:
    """fastq 目录→物种树配置|FASTQ-dir to species-tree config

    输入目录自动识别双端 fastq(_R1/_R2、_1/_2、read1/read2、_R1_001 等),
    WASTER 直接从 reads 建树(免组装免比对,1.5X 覆盖度即可)
    |Auto-detect paired-end fastq in the input dir; WASTER builds the species
    tree directly from reads (no assembly/alignment; 1.5X coverage suffices)
    """

    # 必需|Required
    input_dir: str

    # 输出|Output
    output_dir: str = "./reads2tree_output"

    # 运行|Runtime
    threads: int = 12
    root: str = ""               # 外群物种名|outgroup species name
    branch_length: bool = False  # 03 步 waster_branchlength|optional branch lengths
    samples_map: str = ""        # 个体→物种映射文件|individual-to-species map
    merge: bool = False          # 重叠双端用 BBMerge 合并(默认 cat)|BBMerge overlapping reads
    log_level: str = "INFO"

    # 工具路径(独立二进制/脚本,绝对路径直调不走 conda run)
    # |Tool paths (standalone binaries; direct call, no conda run)
    waster_path: str = field(
        default_factory=lambda: get_tool_path(
            "waster", "~/software/ASTER/bin/waster", "WASTER_PATH"))
    bbmerge_path: str = field(
        default_factory=lambda: get_tool_path(
            "bbmerge.sh", "~/miniforge3/envs/bbmap_v.39.81/bin/bbmerge.sh",
            "BBMERGE_PATH"))

    def __post_init__(self):
        """展开路径、建目录、扫描输入|Expand paths, make dirs, scan input"""
        self.input_dir = expand_path(self.input_dir)
        self.output_dir = expand_path(self.output_dir)
        if self.samples_map:
            self.samples_map = expand_path(self.samples_map)

        out = Path(self.output_dir)
        self.step_input_dir = out / "01_input"
        self.uncompressed_dir = self.step_input_dir / "uncompressed"
        self.waster_dir = out / "02_waster"
        self.branch_dir = out / "03_branch_length"
        self.info_dir = out / "00_pipeline_info"
        self.logs_dir = out / "99_logs"
        dirs = [self.step_input_dir, self.uncompressed_dir, self.waster_dir,
                self.info_dir, self.logs_dir]
        if self.branch_length:
            dirs.append(self.branch_dir)
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)

        # 关键产物路径|key artifact paths
        self.input_tsv = str(self.step_input_dir / "input.tsv")
        self.samples_map_tsv = str(self.step_input_dir / "samples_map.tsv")
        self.species_tree = str(self.waster_dir / "waster.species_tree.nw")
        self.waster_log = str(self.waster_dir / "waster.log")
        self.bl_tree = str(self.branch_dir / "waster_branchlength.species_tree.nw")
        self.bl_log = str(self.branch_dir / "waster_branchlength.log")

        # 扫描输入(目录不存在则留空,交 validate 报错)|scan input
        self.raw_fastq: List[str] = []
        self.paired: List[Tuple[str, List[str], List[str]]] = []
        self.singles: List[str] = []
        self.ignored_files: List[str] = []
        if Path(self.input_dir).is_dir():
            (self.raw_fastq, self.ignored_files,
             self.paired, self.singles) = scan_fastq_dir(self.input_dir)

    @property
    def waster_branchlength_path(self) -> str:
        """由 waster 同目录推导|Derived from waster's directory"""
        return str(Path(self.waster_path).parent / "waster_branchlength")

    def load_samples_map(self) -> Dict[str, str]:
        """加载映射(缺失/格式错返回空,错误由 validate 报)|Load map (empty if absent)"""
        try:
            from .utils import parse_samples_map
            return parse_samples_map(self.samples_map) if self.samples_map else {}
        except Exception:
            return {}

    def validate(self) -> bool:
        """校验配置(错误收集后一次抛)|Validate; raise ValueError listing all errors"""
        errors = []
        if not Path(self.input_dir).is_dir():
            errors.append(f"输入目录不存在|Input dir not found: {self.input_dir}")
        if not self.raw_fastq:
            errors.append(f"输入目录中未找到 fastq 文件|No fastq files found in "
                          f"{self.input_dir}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Threads must be positive: {self.threads}")
        if self.merge and not Path(self.bbmerge_path).exists():
            errors.append(f"BBMerge 工具不存在(--merge 需要)|bbmerge.sh not found "
                          f"(--merge requires it): {self.bbmerge_path}")
        if self.samples_map and not Path(self.samples_map).is_file():
            errors.append(f"样本映射文件不存在|Sample map not found: {self.samples_map}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
