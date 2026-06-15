"""
MEME Parser配置管理模块|MEME Parser Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from ..common.paths import expand_path


@dataclass
class MemeParserConfig:
    """MEME Parser配置类|MEME Parser Configuration Class"""

    # 必需参数|Required parameters
    input_file: str  # 可以是FASTA文件或MEME输出文件(xml/txt)

    # 输出配置|Output configuration
    output_prefix: str = 'meme_results'
    output_dir: str = '.'

    # 运行模式|Run mode
    parse_only: bool = False  # True=只解析现有输出，False=运行MEME

    # MEME软件路径|MEME software path
    meme_path: str = '~/miniforge3/envs/meme_v.5.5.9/bin/meme'

    # MEME参数|MEME parameters
    protein: bool = True  # -protein 输入序列为蛋白质
    dna: bool = False  # -dna 输入序列为DNA
    mod: str = 'zoops'  # -mod motif分布模式
    nmotifs: int = 10  # -nmotifs motif数量
    minw: int = 6  # -minw 最小motif宽度
    maxw: int = 50  # -maxw 最大motif宽度
    objfun: str = 'classic'  # -objfun 目标函数
    markov_order: int = 0  # -markov_order Markov链阶数

    # 输出格式|Output format
    output_tsv: bool = True
    output_csv: bool = True
    output_excel: bool = True

    # 序列提取选项|Sequence extraction options
    extract_motif_seqs: bool = True  # 提取motif序列|Extract motif sequences (default: True)
    input_fasta: str = None  # 原始FASTA文件路径（解析模式时需要）|Original FASTA file path (needed for parse-only mode)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 展开meme_path中的~|Expand ~ in meme_path
        self.meme_path = expand_path(self.meme_path)

        # 如果是运行模式，输出目录就是MEME的输出目录
        if not self.parse_only:
            # 设置MEME输出目录为output_dir
            pass

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        # 检查MEME软件（仅在运行模式）|Check MEME software (run mode only)
        if not self.parse_only:
            if not self._check_meme_installed():
                errors.append(f"MEME软件未找到或不可执行|MEME software not found or not executable: {self.meme_path}")

        # 检查mod参数|Check mod parameter
        valid_mods = ['zoops', 'anr', 'oor']
        if self.mod not in valid_mods:
            errors.append(f"无效的mod参数|Invalid mod parameter: {self.mod} (valid: {valid_mods})")

        # 检查序列类型参数互斥|Check sequence type parameters are mutually exclusive
        if self.protein and self.dna:
            errors.append(f"-protein和-dna参数不能同时使用|-protein and -dna parameters cannot be used together")

        # 检查objfun参数|Check objfun parameter
        valid_objfuns = ['classic', 'de', 'ce', 'cd']
        if self.objfun not in valid_objfuns:
            errors.append(f"无效的objfun参数|Invalid objfun parameter: {self.objfun} (valid: {valid_objfuns})")

        # 检查宽度参数|Check width parameters
        if self.minw <= 0:
            errors.append(f"最小宽度必须为正数|Min width must be positive: {self.minw}")

        if self.maxw <= 0:
            errors.append(f"最大宽度必须为正数|Max width must be positive: {self.maxw}")

        if self.minw > self.maxw:
            errors.append(f"最小宽度不能大于最大宽度|Min width cannot be greater than max width")

        if self.nmotifs <= 0:
            errors.append(f"motif数量必须为正数|Number of motifs must be positive: {self.nmotifs}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def _check_meme_installed(self):
        """检查MEME是否已安装|Check if MEME is installed"""
        import shutil
        return shutil.which(self.meme_path) is not None

    def get_meme_output_dir(self):
        """获取MEME输出目录路径|Get MEME output directory path"""
        if self.parse_only:
            # 解析模式：输入文件所在目录
            return str(Path(self.input_file).parent)
        else:
            # 运行模式：使用配置的输出目录
            return str(Path(self.output_dir) / f"{self.output_prefix}_meme_out")

    def get_meme_xml_path(self):
        """获取MEME XML文件路径|Get MEME XML file path"""
        if self.parse_only:
            # 解析模式：输入文件就是XML/TXT
            return self.input_file
        else:
            # 运行模式：XML文件在输出目录中
            return str(Path(self.get_meme_output_dir()) / "meme.xml")

    def get_meme_txt_path(self):
        """获取MEME TXT文件路径|Get MEME TXT file path"""
        if self.parse_only:
            return self.input_file
        else:
            return str(Path(self.get_meme_output_dir()) / "meme.txt")
