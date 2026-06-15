"""
三代转录组比对配置管理模块|Long RNA-seq Alignment Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class LongRNASeqConfig:
    """三代转录组比对配置类|Long RNA-seq Alignment Configuration Class"""

    # 必需参数|Required parameters
    input_file: str  # 输入文件（BAM或FASTQ）或文件夹|Input file (BAM or FASTQ) or directory
    ref_genome: str  # 参考基因组文件|Reference genome file
    output_dir: str  # 输出目录|Output directory
    sample_name: str = None  # 样本名称|Sample name (可选，默认从输入文件名提取)

    # 比对参数|Alignment parameters
    threads: int = 64  # 线程数|Number of threads
    max_intron: int = 100000  # 最大intron长度|Maximum intron length
    min_mapq: int = 20  # 最小mapping quality|Minimum mapping quality
    secondary: bool = True  # 是否输出次优比对|Whether to output secondary alignments

    # minimap2额外参数|Additional minimap2 parameters
    minimap2_path: str = "minimap2"  # minimap2路径|minimap2 path
    samtools_path: str = "samtools"  # samtools路径|samtools path

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 规范化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.ref_genome = os.path.normpath(os.path.abspath(self.ref_genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 检查是否为文件夹|Check if input is directory
        self.is_directory = os.path.isdir(self.input_file)

        # 如果是文件夹，查找所有符合条件的文件|If directory, find all matching files
        if self.is_directory:
            self.input_files = self._find_input_files()
            if not self.input_files:
                raise ValueError(f"文件夹内未找到BAM或FASTQ文件|No BAM or FASTQ files found in directory: {self.input_file}")
        else:
            self.input_files = [self.input_file]

        # 检测输入文件类型|Detect input file type (for the first file if directory)
        self.input_type = self._detect_input_type()

        # 创建输出目录结构|Create output directory structure
        self.output_path = Path(self.output_dir)
        self.result_dir = self.output_path
        self.log_dir = self.result_dir / "logs"
        self.stats_dir = self.result_dir / "stats"
        self.tmp_dir = self.result_dir / "tmp"

        for dir_path in [self.result_dir, self.log_dir, self.stats_dir, self.tmp_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        # 如果未指定样本名，从输入文件提取|Extract sample name from input file if not specified
        if self.sample_name is None:
            if self.is_directory:
                # 文件夹模式：使用文件夹名称作为样本名|Directory mode: use folder name as sample name
                self.sample_name = Path(self.input_file).name
            else:
                # 单文件模式：从文件名提取|Single file mode: extract from filename
                self.sample_name = self._extract_sample_name(self.input_file)

        # 日志文件|Log file
        self.log_file = self.log_dir / f"alignment_{self.sample_name}.log"

    def _extract_sample_name(self, file_path: str) -> str:
        """
        从输入文件路径提取样本名称|Extract sample name from input file path

        Args:
            file_path: 文件路径|File path

        Returns:
            str: 样本名称|Sample name
        """
        path = Path(file_path)
        name = path.name  # 获取文件名（含扩展名）|Get filename with extension

        # 去除扩展名|Remove extensions
        # 处理 .fastq.gz, .fq.gz 等双重扩展名|Handle double extensions like .fastq.gz, .fq.gz
        while name.endswith('.gz'):
            name = name[:-3]  # 去掉 .gz
        # 去掉最后的扩展名|Remove final extension
        name = str(Path(name).stem)  # 使用stem去掉扩展名，例如 sample.bam -> sample

        # 去除常见的中间标记|Remove common intermediate markers
        # 按优先级排序，先匹配长的|Sort by priority, match longer ones first
        intermediate_patterns = [
            '.clean', '.trimmed', '.filtered',
            '_clean', '_trimmed', '_filtered',
            '.fastq', '.fq',  # 有时候文件名里还会包含扩展名
        ]

        for pattern in intermediate_patterns:
            if name.endswith(pattern):
                name = name[:-len(pattern)]
                # 去掉后重新检查是否还有其他标记|After removal, recheck for other markers
                # 例如 sample_1.clean.fastq -> sample_1 -> sample
                break

        # 去除常见的测序标记|Remove common sequencing markers
        # R1/R2, 1/2 等标记|Markers like R1/R2, 1/2
        patterns_to_remove = ['_R1', '_R2', '_1', '_2', '-R1', '-R2', '-1', '-2']

        for pattern in patterns_to_remove:
            if name.endswith(pattern):
                name = name[:-len(pattern)]
                break

        # 如果样本名为空，使用默认值|Use default if sample name is empty
        if not name:
            name = "sample"

        return name

    def _find_input_files(self) -> list:
        """
        在文件夹内查找所有BAM和FASTQ文件|Find all BAM and FASTQ files in directory

        Returns:
            list: 文件路径列表|List of file paths
        """
        import glob

        input_dir = self.input_file
        files = []

        # 支持的扩展名|Supported extensions
        bam_patterns = ['*.bam', '*.BAM']
        fastq_patterns = ['*.fq', '*.fastq', '*.fq.gz', '*.fastq.gz', '*.FQ.gz', '*.FASTQ.GZ', '*.clean.fastq.gz']

        # 查找BAM文件|Find BAM files
        for pattern in bam_patterns:
            files.extend(glob.glob(os.path.join(input_dir, pattern)))

        # 查找FASTQ文件|Find FASTQ files
        for pattern in fastq_patterns:
            found = glob.glob(os.path.join(input_dir, pattern))
            for f in found:
                # 排除R2文件（只保留R1）|Exclude R2 files (keep only R1)
                if not self._is_r2_file(f):
                    files.append(f)

        return sorted(files)

    def _is_r2_file(self, file_path: str) -> bool:
        """
        判断是否为R2文件|Check if file is R2

        Args:
            file_path: 文件路径|File path

        Returns:
            bool: 是否为R2文件|Whether is R2 file
        """
        name = Path(file_path).name.upper()
        r2_patterns = ['_R2.', '_2.', '-R2.', '-2.']

        for pattern in r2_patterns:
            if pattern in name:
                return True

        return False

    def _detect_input_type(self) -> str:
        """
        检测输入文件类型|Detect input file type

        Returns:
            str: 'bam', 'fastq_single', 'fastq_paired', 'directory'
        """
        # 如果是文件夹，直接返回|If directory, return directly
        if self.is_directory:
            return 'directory'

        path = Path(self.input_file)
        file_lower = str(path).lower()

        # BAM文件|BAM file
        if file_lower.endswith('.bam'):
            return 'bam'

        # FASTQ文件 (支持压缩格式)|FASTQ file (supports compressed formats)
        if file_lower.endswith('.fq.gz') or file_lower.endswith('.fastq.gz') or file_lower.endswith('.fq') or file_lower.endswith('.fastq'):
            # 检查是否为双端测序
            # 文件名模式检查|Filename pattern check
            stem = path.stem  # 去掉最后的扩展名
            if stem.endswith('.gz'):
                stem = Path(stem).stem

            # 常见的双端测序命名模式|Common paired-end naming patterns
            paired_patterns = ['_R1', '_R2', '_1', '_2', '.r1', '.r2', '-R1', '-R2', '-1', '-2']

            for pattern in paired_patterns:
                if pattern.upper() in stem.upper():
                    # 找到R1或R2，说明是双端测序|Found R1 or R2, indicating paired-end
                    return 'fastq_paired'

            return 'fastq_single'

        raise ValueError(f"无法识别的输入文件格式|Unrecognized input file format: {self.input_file}. 支持的格式|Supported formats: .bam, .fq, .fastq, .fq.gz, .fastq.gz")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查参考基因组|Check reference genome
        if not os.path.exists(self.ref_genome):
            errors.append(f"参考基因组文件不存在|Reference genome not found: {self.ref_genome}")

        # 检查输入路径|Check input path
        if not os.path.exists(self.input_file):
            errors.append(f"输入路径不存在|Input path not found: {self.input_file}")

        # 如果是文件夹，检查是否找到文件|If directory, check if files found
        if self.is_directory:
            if not self.input_files:
                errors.append(f"文件夹内未找到BAM或FASTQ文件|No BAM or FASTQ files found in directory: {self.input_file}")
        else:
            # 如果是单个文件且是双端FASTQ，检查R2文件是否存在|If single file and paired-end FASTQ, check R2 file exists
            if self.input_type == 'fastq_paired':
                r2_file = self._find_r2_file()
                if r2_file and not os.path.exists(r2_file):
                    errors.append(f"未找到配对的R2文件|R2 file not found: {r2_file}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Threads must be positive: {self.threads}")

        if self.max_intron <= 0:
            errors.append(f"最大intron长度必须为正整数|Max intron must be positive: {self.max_intron}")

        if self.min_mapq < 0 or self.min_mapq > 60:
            errors.append(f"最小MAPQ必须在0-60之间|Min MAPQ must be between 0-60: {self.min_mapq}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def _find_r2_file(self) -> Optional[str]:
        """
        查找配对的R2文件|Find paired R2 file

        Returns:
            str: R2文件路径|R2 file path，如果不存在返回None|return None if not found
        """
        path = Path(self.input_file)
        name = path.name
        stem = path.stem  # 例如：sample_R1.fastq.gz -> sample_R1

        # 尝试各种替换模式|Try various replacement patterns
        r2_patterns = [
            stem.replace('_R1', '_R2'),
            stem.replace('_R2', '_R1'),
            stem.replace('_1', '_2'),
            stem.replace('_2', '_1'),
            stem.replace('-R1', '-R2'),
            stem.replace('-R2', '-R1'),
            stem.replace('-1', '-2'),
            stem.replace('-2', '-1'),
            stem.replace('.r1', '.r2'),
            stem.replace('.r2', '.r1'),
        ]

        # 检查压缩扩展名|Check compression extension
        if str(path).endswith('.gz'):
            ext = '.fastq.gz'
        else:
            ext = '.fastq'

        for r2_stem in r2_patterns:
            if r2_stem != stem:  # 确保确实做了替换|Ensure replacement was made
                r2_path = path.parent / (r2_stem + ext)
                if r2_path.exists():
                    return str(r2_path)

        return None

    def get_minimap2_params(self) -> List[str]:
        """获取minimap2参数列表|Get minimap2 parameters list"""
        params = [
            "-a",  # 输出SAM格式(包含header)|Output SAM format (with header)
            "-x", "splice",  # 剪接比对模式(长读段)|Splice alignment mode (long reads)
            "--eqx",  # 使用=和X表示匹配/错配|Use = and X for match/mismatch
            "--MD",  # 生成MD标签|Generate MD tag
            "--cs",  # 生成CIGAR字符串|Generate CIGAR string
            "-Y",  # 使用软剪切|Use soft clipping
            "--sam-hit-only",  # 只输出有比对的reads|Only output aligned reads
            "-G", str(self.max_intron),  # 最大intron长度|Maximum intron length
            "-M", "0.9",  # 最小比对分数|Min alignment score fraction
            "-p", "0.1",  # 次要比对分数阈值|Secondary alignment score threshold
            "-t", str(self.threads),  # 线程数|Number of threads
        ]

        if not self.secondary:
            params.append("--secondary=no")

        return params
