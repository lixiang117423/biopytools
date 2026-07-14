"""
序列提取配置|Sequence extraction configuration
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import get_tool_path, expand_path


@dataclass
class SeqExtractConfig:
    """序列提取配置|Sequence extraction config"""

    # 必需参数|Required parameters
    input_query: str           # 查询：单个ID、ID文件或BED文件|Query: single ID, ID file, or BED file
    sequence_file: str         # 目标序列FASTA文件|Target sequence FASTA file
    output_file: Optional[str] = None  # 输出文件(可选,自动推导)|Output file (optional, auto-derived)

    # 可选参数|Optional parameters
    force_bed: bool = False    # 强制BED模式(跳过自动检测)|Force BED mode (skip auto-detection)
    seqkit_path: str = field(
        default_factory=lambda: get_tool_path(
            "seqkit",
            "seqkit",
            "SEQKIT_PATH",
        )
    )

    # 内部状态|Internal state
    query_type: str = ""       # "single_id" | "id_file" | "bed_file"
    effective_output: str = "" # 实际使用的输出路径|Actual output path used

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.seqkit_path = expand_path(self.seqkit_path)
        self._detect_query_type()
        self._derive_output()

    def _detect_query_type(self):
        """自动检测查询类型|Auto-detect query type"""
        if self.force_bed:
            if not os.path.isfile(self.input_query):
                raise ValueError(
                    f"强制BED模式但文件不存在|Force BED mode but file not found: {self.input_query}"
                )
            self.query_type = "bed_file"
            return

        # 不是文件 -> 单ID模式|Not a file -> single ID mode
        if not os.path.isfile(self.input_query):
            self.query_type = "single_id"
            return

        # 是文件，读取第一行判断|Is a file, read first line to determine
        try:
            with open(self.input_query, 'r') as f:
                first_line = f.readline().strip()
        except Exception as e:
            raise ValueError(
                f"无法读取查询文件|Cannot read query file: {self.input_query}\n{e}"
            )

        if not first_line:
            raise ValueError(f"查询文件为空|Query file is empty: {self.input_query}")

        # Tab分隔列数 >= 2 -> BED模式|Tab-separated columns >= 2 -> BED mode
        columns = first_line.split('\t')
        if len(columns) >= 2:
            self.query_type = "bed_file"
        else:
            self.query_type = "id_file"

    def _derive_output(self):
        """自动推导输出文件名|Auto-derive output filename"""
        if self.output_file:
            self.effective_output = self.output_file
            return

        query_stem = Path(self.input_query).stem if os.path.isfile(self.input_query) else self.input_query
        subject_stem = Path(self.sequence_file).stem
        self.effective_output = f"{query_stem}.{subject_stem}.fa"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 序列文件必须存在|Sequence file must exist
        if not os.path.isfile(self.sequence_file):
            errors.append(f"序列文件不存在|Sequence file not found: {self.sequence_file}")

        # BED/ID文件模式时文件必须存在|BED/ID file mode requires file to exist
        if self.query_type in ("bed_file", "id_file"):
            if not os.path.isfile(self.input_query):
                errors.append(f"查询文件不存在|Query file not found: {self.input_query}")

        # 输出目录必须可创建|Output directory must be creatable
        output_dir = os.path.dirname(os.path.abspath(self.effective_output))
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
            except Exception as e:
                errors.append(f"无法创建输出目录|Cannot create output directory: {output_dir}\n{e}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
