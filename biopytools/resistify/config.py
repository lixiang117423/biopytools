"""
Resistify配置管理模块|Resistify Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Tuple
from ..common.paths import expand_path, get_tool_path


FASTA_EXTENSIONS = {'.fa', '.fasta', '.faa', '.pep', '.fna'}


@dataclass
class ResistifyConfig:
    """Resistify配置类|Resistify Configuration Class"""

    # 必需参数|Required parameters
    input_file: str  # 支持单个文件或目录|Supports single file or directory

    # 输出配置|Output configuration
    output_dir: str = './resistify_output'
    output_prefix: str = 'resistify_results'

    # Resistify工具参数|Resistify tool parameters
    resistify_path: str = field(
        default_factory=lambda: get_tool_path('resistify', '~/miniforge3/envs/resistify_v.1.3.0/bin/resistify', 'RESISTIFY_PATH')
    )
    skip_resistify: bool = False
    threads: int = 12

    # 序列提取选项|Sequence extraction options
    extract_nlr_sequences: bool = False
    extract_nbarc_sequences: bool = False

    # 筛选选项|Filtering options
    filter_classification: Optional[str] = None
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    min_lrr_length: Optional[int] = None

    # 输出格式|Output format
    output_tsv: bool = True
    output_csv: bool = True
    output_excel: bool = True

    # 是否包含motifs详情|Whether to include motifs details
    include_motifs: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_file = os.path.normpath(os.path.abspath(expand_path(self.input_file)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        self.resistify_path = os.path.normpath(os.path.abspath(expand_path(self.resistify_path)))

        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    @property
    def is_directory(self) -> bool:
        """是否为目录输入|Whether input is a directory"""
        return os.path.isdir(self.input_file)

    @property
    def input_files(self) -> List[str]:
        """获取输入FASTA文件列表|Get list of input FASTA files"""
        if not self.is_directory:
            return [self.input_file]
        files = []
        for f in sorted(Path(self.input_file).iterdir()):
            if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS:
                files.append(str(f))
        return files

    @property
    def sample_names(self) -> List[str]:
        """获取样本名称列表|Get list of sample names"""
        return [Path(f).stem for f in self.input_files]

    def sample_output_dir(self, sample_name: str) -> str:
        """获取样本输出目录|Get sample output directory"""
        return os.path.join(self.output_dir, sample_name)

    def get_batch_tasks(self) -> List[Tuple[str, str, str]]:
        """
        获取批量处理任务列表|Get list of batch processing tasks

        Returns:
            List of (sample_name, input_path, resistify_output_dir) tuples
        """
        if not self.is_directory:
            return [(Path(self.input_file).stem, self.input_file, self.resistify_output_dir)]

        if self.skip_resistify:
            # 跳过模式：查找包含results.tsv的子目录|Skip mode: find subdirectories with results.tsv
            tasks = []
            for d in sorted(Path(self.input_file).iterdir()):
                if d.is_dir() and (d / 'results.tsv').exists():
                    tasks.append((d.name, str(d), str(d)))
            return tasks
        else:
            return [
                (Path(f).stem, f, self.sample_output_dir(Path(f).stem))
                for f in self.input_files
            ]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径存在性|Check input path existence
        if not os.path.exists(self.input_file):
            errors.append(f"输入路径不存在|Input path not found: {self.input_file}")

        if self.is_directory:
            if not self.skip_resistify:
                if len(self.input_files) == 0:
                    errors.append(
                        f"输入目录中未找到FASTA文件(支持: {', '.join(sorted(FASTA_EXTENSIONS))})|"
                        f"No FASTA files found in input directory (supported: {', '.join(sorted(FASTA_EXTENSIONS))}): {self.input_file}"
                    )
            else:
                if len(self.get_batch_tasks()) == 0:
                    errors.append(
                        f"输入目录中未找到包含results.tsv的子目录|"
                        f"No subdirectories with results.tsv found: {self.input_file}"
                    )
        else:
            if not self.skip_resistify:
                if not os.path.isfile(self.input_file):
                    errors.append(f"输入文件不存在|Input file not found: {self.input_file}")
            else:
                # 仅解析模式：input_file作为resistify输出目录|Parse-only mode: input_file as resistify output dir
                if not os.path.isdir(self.input_file):
                    errors.append(f"Resistify输出目录不存在|Resistify output directory not found: {self.input_file}")

                required_files = ['results.tsv', 'domains.tsv']
                for fname in required_files:
                    file_path = os.path.join(self.input_file, fname)
                    if not os.path.exists(file_path):
                        errors.append(f"必需文件不存在|Required file not found: {fname}")

        # 单文件模式下检查序列提取文件|Check sequence extraction files in single-file mode
        if not self.is_directory:
            resistify_dir = self.output_dir if not self.skip_resistify else self.input_file
            if self.extract_nlr_sequences:
                nlr_fasta = os.path.join(resistify_dir, 'nlr.fasta')
                if not os.path.exists(nlr_fasta):
                    errors.append(f"NLR序列文件不存在|NLR fasta file not found: {nlr_fasta}")

            if self.extract_nbarc_sequences:
                nbarc_fasta = os.path.join(resistify_dir, 'nbarc.fasta')
                if not os.path.exists(nbarc_fasta):
                    errors.append(f"NB-ARC序列文件不存在|NB-ARC fasta file not found: {nbarc_fasta}")

        # 检查长度参数|Check length parameters
        if self.min_length is not None and self.min_length <= 0:
            errors.append(f"最小长度必须为正数|Min length must be positive: {self.min_length}")

        if self.max_length is not None and self.max_length <= 0:
            errors.append(f"最大长度必须为正数|Max length must be positive: {self.max_length}")

        if self.min_length and self.max_length and self.min_length > self.max_length:
            errors.append(f"最小长度不能大于最大长度|Min length cannot be greater than max length")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    @property
    def resistify_output_dir(self) -> str:
        """Resistify工具的输出目录|Resistify tool output directory"""
        if self.skip_resistify:
            return self.input_file
        return self.output_dir
