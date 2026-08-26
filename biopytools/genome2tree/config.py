"""genome2tree 配置|genome2tree configuration"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

from ..common.paths import expand_path, get_tool_path
from .utils import SampleFile, parse_samples_map, scan_input_dir


@dataclass
class Genome2TreeConfig:
    """基因组目录→物种树配置|Genome-dir to species-tree config"""

    # 必需|Required
    input_dir: str

    # 输出|Output
    output_dir: str = "./genome2tree_output"

    # 运行|Runtime
    threads: int = 12
    root: str = ""               # 外群物种名|outgroup species name
    branch_length: bool = False  # 03 步 waster_branchlength|optional branch lengths
    samples_map: str = ""        # 个体→物种映射文件|individual-to-species map
    log_level: str = "INFO"

    # 工具路径(独立 C++ 二进制,绝对路径直调不走 conda run)
    # |Tool path (standalone C++ binary; direct call, no conda run)
    waster_path: str = field(
        default_factory=lambda: get_tool_path(
            "waster", "~/software/ASTER/bin/waster", "WASTER_PATH"))

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
        self.raw_samples: List[SampleFile] = []
        self.ignored_files: List[str] = []
        if Path(self.input_dir).is_dir():
            self.raw_samples, self.ignored_files = scan_input_dir(self.input_dir)

    @property
    def waster_branchlength_path(self) -> str:
        """由 waster 同目录推导|Derived from waster's directory"""
        return str(Path(self.waster_path).parent / "waster_branchlength")

    def load_samples_map(self) -> Dict[str, str]:
        """加载映射(缺失/格式错返回空,错误由 validate 报)|Load map (empty if absent)"""
        if self.samples_map and Path(self.samples_map).is_file():
            try:
                return parse_samples_map(self.samples_map)
            except ValueError:
                return {}
        return {}

    def validate(self) -> bool:
        """校验配置(错误收集后一次抛)|Validate; raise ValueError listing all errors"""
        errors = []
        if not Path(self.input_dir).is_dir():
            errors.append(f"输入目录不存在|Input dir not found: {self.input_dir}")
        else:
            stems = [s.stem for s in self.raw_samples]
            if not stems:
                errors.append(
                    "输入目录无序列文件(.fa/.fasta/.fna/.fas,可.gz)|"
                    "No sequence files found in input dir")
            dups = sorted({x for x in stems if stems.count(x) > 1})
            if dups:
                errors.append(f"样本重名(去后缀去.gz后冲突)|Duplicate stems: {dups}")
            if len(set(stems)) < 3:
                errors.append(
                    f"样本数不足(至少3个不同样本,当前{len(set(stems))})|"
                    f"Need >=3 distinct samples, got {len(set(stems))}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Threads must be positive: {self.threads}")
        if not Path(self.waster_path).exists():
            errors.append(f"waster 不存在|waster not found: {self.waster_path}")
        if self.branch_length and not Path(self.waster_branchlength_path).exists():
            errors.append(
                f"waster_branchlength 不存在|waster_branchlength not found: "
                f"{self.waster_branchlength_path}")
        if self.samples_map:
            if not Path(self.samples_map).is_file():
                errors.append(f"映射文件不存在|samples map not found: {self.samples_map}")
            else:
                try:
                    mapping = parse_samples_map(self.samples_map)
                except ValueError as e:
                    errors.append(str(e))
                else:
                    stems = {s.stem for s in self.raw_samples}
                    unknown = sorted(set(mapping) - stems)
                    if unknown:
                        errors.append(
                            f"映射含未知个体(输入目录无此文件)|"
                            f"map entries not in input dir: {unknown}")
                    if self.raw_samples:
                        species = {mapping.get(s, s) for s in stems}
                        if len(species) < 3:
                            errors.append(
                                f"映射后物种数不足(至少3,当前{len(species)})|"
                                f"Need >=3 distinct species after mapping, got {len(species)}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
