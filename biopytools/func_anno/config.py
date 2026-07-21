"""func_anno 配置|func_anno configuration."""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path


@dataclass
class FuncAnnoConfig:
    """功能注释流水线配置|Functional annotation pipeline config.

    输入蛋白序列 → interproscan(结构域) + eggnog-mapper(GO/KEGG) → 标准 GO/KEGG 表.
    |Protein input → interproscan (domains) + eggnog-mapper (GO/KEGG) → standard tables.
    """

    # 必需|required
    input_file: str
    output_dir: str
    # 通用|general
    threads: int = 12
    sample_name: Optional[str] = None
    # 复用已有结果(避免重跑)|reuse existing results
    ips_result: Optional[str] = None       # 已有 IPS 结果目录|existing IPS dir
    eggnog_result: Optional[str] = None    # 已有 .emapper.annotations|existing annotations
    skip_ips: bool = False                 # 跳过 IPS(只要 GO/KEGG)|skip IPS
    skip_eggnog: bool = False              # 跳过 eggnog(只整理已有)|skip eggnog
    # KEGG|KEGG
    kegg_map: Optional[str] = None         # 外部 KEGG 映射 TSV(补 category)|external KEGG map
    # eggnog 透传|eggnog passthrough
    data_dir: Optional[str] = None
    mode: str = "mmseqs"
    emapper_path: Optional[str] = None

    def __post_init__(self):
        """初始化后: 展开~路径, 建目录结构|Post-init: expand ~ paths, build dirs."""
        # ⚠️ 关键: 展开所有含~的路径|CRITICAL: expand all ~ paths
        self.input_file = expand_path(self.input_file)
        self.output_dir = expand_path(self.output_dir)
        if self.ips_result:
            self.ips_result = expand_path(self.ips_result)
        if self.eggnog_result:
            self.eggnog_result = expand_path(self.eggnog_result)
        if self.kegg_map:
            self.kegg_map = expand_path(self.kegg_map)
        if self.data_dir:
            self.data_dir = expand_path(self.data_dir)
        if self.emapper_path:
            self.emapper_path = expand_path(self.emapper_path)

        # sample_name 默认从输入文件名|default from input stem
        if not self.sample_name:
            self.sample_name = Path(self.input_file).stem

        # by-sample 目录结构(遵循 §12)|by-sample dir layout
        self.output_path = Path(self.output_dir)
        self.sample_dir = self.output_path / self.sample_name
        self.ips_dir = self.sample_dir / "01_interproscan"
        self.eggnog_dir = self.sample_dir / "02_eggnog"
        self.tables_dir = self.sample_dir / "03_tables"
        self.logs_dir = self.sample_dir / "99_logs"

        self.sample_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置|Validate config."""
        errors = []
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")
        if self.threads <= 0:
            errors.append("线程数必须为正|Threads must be positive")
        if self.mode not in ("mmseqs", "diamond", "hmmer"):
            errors.append(f"不支持的搜索模式|Unsupported mode: {self.mode}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
