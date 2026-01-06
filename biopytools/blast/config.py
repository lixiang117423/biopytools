"""
BLAST分析配置模块|BLAST Analysis Configuration Module
标准化BLAST配置类，遵循参数命名规范
"""

import os
import re
from typing import Optional
from ..core.config import BaseConfig


class BLASTConfig(BaseConfig):
    """BLAST分析配置类|BLAST Analysis Configuration Class"""

    def __init__(
        self,
        input: Optional[str] = None,
        output: Optional[str] = None,
        reference: Optional[str] = None,
        prefix: str = "blast_output",
        threads: int = 12,
        quality: float = 1e-5,
        memory: str = "8G",
        sample_id: Optional[str] = None,
        sample_name: Optional[str] = None,
        read_group: Optional[str] = None,
        min_quality: float = 20,
        min_length: int = 50,
        min_depth: int = 10,
        max_depth: int = 1000,
        mapping_quality: int = 20,
        log_level: str = "INFO",
        log_file: Optional[str] = None,
        force: bool = False,
        dry_run: bool = False,
        keep_intermediate: bool = False,
        tmp_dir: str = "/tmp",
        timeout: Optional[int] = None,
        # BLAST特定参数
        blast_type: str = "blastn",
        evalue: float = 1e-5,
        max_target_seqs: int = 10,
        word_size: Optional[int] = None,
        target_db_type: str = "nucl",
        min_identity: float = 70.0,
        min_coverage: float = 50.0,
        high_quality_evalue: float = 1e-10,
        input_suffix: str = "*.fa",
        auto_detect_samples: bool = True,
        sample_name_pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
        sample_map_file: Optional[str] = None,
        # BLAST工具路径
        makeblastdb_path: str = "makeblastdb",
        blastn_path: str = "blastn",
        blastp_path: str = "blastp",
        blastx_path: str = "blastx",
        tblastn_path: str = "tblastn",
        tblastx_path: str = "tblastx",
        # 比对可视化参数|Alignment visualization parameters
        alignment_output: str = "both",
        alignment_width: int = 80,
        alignment_min_identity: float = 0.0,
        alignment_min_coverage: float = 0.0,
        alignment_max_per_sample: int = 100,
        html_theme: str = "modern"
    ):
        """
        初始化BLAST配置|Initialize BLAST configuration

        Args:
            input: 输入文件路径|Input file path
            output: 输出目录路径|Output directory path
            reference: 目标数据库文件路径|Target database file path
            prefix: 输出文件前缀|Output prefix
            threads: 线程数|Number of threads
            quality: 质量阈值|Quality threshold
            memory: 内存限制|Memory limit
            sample_id: 样本ID|Sample ID
            sample_name: 样本名称|Sample name
            read_group: Read Group信息|Read Group information
            min_quality: 最小质量值|Minimum quality value
            min_length: 最小序列长度|Minimum sequence length
            min_depth: 最小测序深度|Minimum sequencing depth
            max_depth: 最大测序深度|Maximum sequencing depth
            mapping_quality: 最小mapping质量|Minimum mapping quality
            log_level: 日志级别|Log level
            log_file: 日志文件路径|Log file path
            force: 强制覆盖已存在文件|Force overwrite existing files
            dry_run: 模拟运行不执行|Dry run without execution
            keep_intermediate: 保留中间文件|Keep intermediate files
            tmp_dir: 临时文件目录|Temporary directory
            timeout: 超时时间(秒)|Timeout in seconds
            blast_type: BLAST程序类型|BLAST program type
            evalue: E-value阈值|E-value threshold
            max_target_seqs: 最大目标序列数|Maximum target sequences
            word_size: 词大小|Word size
            target_db_type: 目标数据库类型|Target database type
            min_identity: 最小序列相似度|Minimum sequence identity
            min_coverage: 最小覆盖度|Minimum coverage
            high_quality_evalue: 高质量比对E-value阈值|High quality alignment E-value
            input_suffix: 输入文件后缀模式|Input file suffix pattern
            auto_detect_samples: 自动检测样品名称|Auto-detect sample names
            sample_name_pattern: 样品名称提取正则表达式|Sample name extraction regex
            sample_map_file: 样品映射文件|Sample mapping file
            makeblastdb_path: makeblastdb程序路径|makeblastdb program path
            blastn_path: blastn程序路径|blastn program path
            blastp_path: blastp程序路径|blastp program path
            blastx_path: blastx程序路径|blastx program path
            tblastn_path: tblastn程序路径|tblastn program path
            tblastx_path: tblastx程序路径|tblastx program path
            alignment_output: 比对可视化输出格式|Alignment visualization output format
            alignment_width: 比对每行显示的字符数|Characters per line in alignment display
            alignment_min_identity: 比对可视化最小相似度过滤|Minimum identity for alignment visualization
            alignment_min_coverage: 比对可视化最小覆盖度过滤|Minimum coverage for alignment visualization
            alignment_max_per_sample: 每个样品最多显示的比对数|Maximum alignments to display per sample
            html_theme: HTML主题样式|HTML theme style
        """
        super().__init__()

        # 基本参数|Basic parameters
        self.input = input
        self.output = output
        self.reference = reference
        self.prefix = prefix
        self.threads = self.validate_threads(threads)
        self.quality = self.validate_quality(quality)
        self.memory = memory
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.read_group = read_group
        self.min_quality = min_quality
        self.min_length = min_length
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.mapping_quality = mapping_quality

        # 日志参数|Logging parameters
        self.log_level = log_level
        self.log_file = log_file
        self.force = force
        self.dry_run = dry_run
        self.keep_intermediate = keep_intermediate

        # 执行控制参数|Execution control parameters
        self.tmp_dir = tmp_dir
        self.timeout = timeout

        # BLAST特定参数|BLAST-specific parameters
        self.blast_type = blast_type
        self.evalue = self.validate_quality(evalue, 0, 1)
        self.max_target_seqs = max_target_seqs
        self.word_size = word_size
        self.target_db_type = target_db_type
        self.min_identity = self.validate_quality(min_identity, 0, 100)
        self.min_coverage = self.validate_quality(min_coverage, 0, 100)
        self.high_quality_evalue = self.validate_quality(high_quality_evalue, 0, 1)
        self.input_suffix = input_suffix
        self.auto_detect_samples = auto_detect_samples
        self.sample_name_pattern = sample_name_pattern
        self.sample_map_file = sample_map_file

        # 工具路径参数|Tool path parameters
        self.makeblastdb_path = makeblastdb_path
        self.blastn_path = blastn_path
        self.blastp_path = blastp_path
        self.blastx_path = blastx_path
        self.tblastn_path = tblastn_path
        self.tblastx_path = tblastx_path

        # 比对可视化参数|Alignment visualization parameters
        self.alignment_output = self._validate_alignment_output(alignment_output)
        self.alignment_width = alignment_width
        self.alignment_min_identity = self.validate_quality(alignment_min_identity, 0, 100)
        self.alignment_min_coverage = self.validate_quality(alignment_min_coverage, 0, 100)
        self.alignment_max_per_sample = alignment_max_per_sample
        self.html_theme = self._validate_html_theme(html_theme)

        # 设置默认输出路径|Set default output path
        if self.output is None:
            self.output = self.prefix

        # 验证BLAST类型|Validate BLAST type
        self._validate_blast_type()

        # 根据BLAST类型自动设置数据库类型|Auto-set database type based on BLAST type
        self._auto_set_db_type()

        # 设置默认word_size|Set default word_size
        self._set_default_word_size()

        # 标记是否自动生成映射文件|Mark if auto-generated map
        self.auto_generated_map = False

    @property
    def output_path(self):
        """输出目录的Path对象|Path object for output directory"""
        from pathlib import Path
        return Path(self.output)

    @property
    def alignment_output_dir(self):
        """比对可视化输出子目录|Alignment visualization output subdirectory"""
        return "alignments"

    def _validate_blast_type(self):
        """验证BLAST类型|Validate BLAST type"""
        valid_blast_types = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
        if self.blast_type not in valid_blast_types:
            raise ValueError(f"无效的BLAST类型|Invalid BLAST type: {self.blast_type}")

    def _set_default_word_size(self):
        """设置默认word_size|Set default word_size"""
        if self.word_size is None:
            # 根据BLAST类型设置默认word_size
            default_word_sizes = {
                'blastn': 11,
                'tblastx': 11,
                'blastp': 3,
                'blastx': 3,
                'tblastn': 3
            }
            self.word_size = default_word_sizes.get(self.blast_type, 11)

    def _auto_set_db_type(self):
        """根据BLAST类型自动设置数据库类型|Auto-set database type based on BLAST type"""
        # 如果用户没有明确指定target_db_type为非默认值，则自动推断
        # 如果用户明确指定了，则使用用户指定的值
        # 但为了兼容性，我们根据blast_type给出建议的默认值

        # blast_type到数据库类型的映射
        # blastn: 查询nucl vs 数据库nucl
        # blastp: 查询prot vs 数据库prot
        # blastx: 查询nucl vs 数据库prot
        # tblastn: 查询prot vs 数据库nucl
        # tblastx: 查询prot vs 数据库nucl
        blast_to_db = {
            'blastn': 'nucl',
            'blastp': 'prot',
            'blastx': 'prot',
            'tblastn': 'nucl',
            'tblastx': 'nucl'
        }

        # 自动设置正确的数据库类型
        self.target_db_type = blast_to_db.get(self.blast_type, 'nucl')

    def _validate_alignment_output(self, alignment_output: str) -> str:
        """验证比对可视化输出格式|Validate alignment visualization output format"""
        valid_choices = ['none', 'text', 'html', 'both']
        if alignment_output not in valid_choices:
            raise ValueError(f"无效的比对输出格式|Invalid alignment output format: {alignment_output}")
        return alignment_output

    def _validate_html_theme(self, html_theme: str) -> str:
        """验证HTML主题|Validate HTML theme"""
        valid_themes = ['modern', 'classic', 'dark']
        if html_theme not in valid_themes:
            raise ValueError(f"无效的HTML主题|Invalid HTML theme: {html_theme}")
        return html_theme

    def get_output_path(self, filename: str) -> str:
        """
        获取输出文件的完整路径|Get full path of output file

        Args:
            filename: 文件名|File name

        Returns:
            str: 完整路径|Full path
        """
        return os.path.join(self.output, filename)

    def get_target_db_path(self) -> str:
        """
        获取目标数据库路径|Get target database path

        Returns:
            str: 数据库路径|Database path
        """
        base_name = os.path.splitext(os.path.basename(self.reference))[0]
        return os.path.join(self.output, f"{base_name}.db")

    def get_blast_output_path(self) -> str:
        """
        获取BLAST输出路径|Get BLAST output path

        Returns:
            str: 输出路径|Output path
        """
        base_name = os.path.splitext(os.path.basename(self.reference))[0]
        return os.path.join(self.output, f"{base_name}_{self.blast_type}_results.tsv")

    def get_summary_output_path(self) -> str:
        """
        获取汇总输出路径|Get summary output path

        Returns:
            str: 输出路径|Output path
        """
        return os.path.join(self.output, "blast_summary_results.tsv")

    def get_alignment_output_path(self, output_type: str = "html") -> str:
        """
        获取比对可视化输出路径|Get alignment visualization output path

        Args:
            output_type: 输出类型 (html, text)|Output type (html, text)

        Returns:
            str: 输出路径|Output path
        """
        if output_type == "html":
            return os.path.join(self.output, "alignments", "html", "index.html")
        elif output_type == "text":
            return os.path.join(self.output, "alignments", "text", "all_samples_alignments.txt")
        else:
            raise ValueError(f"不支持的输出类型|Unsupported output type: {output_type}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        # 调用父类验证
        super().validate()

        errors = []

        # 验证输入文件|Validate input files
        if not self.input and not self.sample_map_file:
            errors.append("必须指定输入文件(-i/--input)或样品映射文件|Must specify input file or sample mapping file")

        if self.input:
            try:
                self.validate_input_file(self.input, "输入文件|Input file")
            except FileNotFoundError as e:
                errors.append(str(e))

        # 验证目标数据库文件|Validate target database file
        if self.reference:
            try:
                self.validate_input_file(self.reference, "目标数据库文件|Target database file")
            except FileNotFoundError as e:
                errors.append(str(e))

        # 验证样品映射文件|Validate sample mapping file
        if self.sample_map_file and not self.auto_generated_map:
            try:
                self.validate_input_file(self.sample_map_file, "样品映射文件|Sample mapping file")
            except FileNotFoundError as e:
                errors.append(str(e))

        # 验证参数范围|Validate parameter ranges
        if self.evalue <= 0 or self.evalue > 1:
            errors.append(f"E-value阈值必须在0-1之间|E-value must be between 0 and 1: {self.evalue}")

        if self.high_quality_evalue <= 0 or self.high_quality_evalue > 1:
            errors.append(f"高质量E-value阈值必须在0-1之间|High quality E-value must be between 0 and 1: {self.high_quality_evalue}")

        if not 0 <= self.min_identity <= 100:
            errors.append(f"最小相似度必须在0-100之间|Minimum identity must be between 0 and 100: {self.min_identity}")

        if not 0 <= self.min_coverage <= 100:
            errors.append(f"最小覆盖度必须在0-100之间|Minimum coverage must be between 0 and 100: {self.min_coverage}")

        if self.max_target_seqs < 1:
            errors.append(f"最大目标序列数必须大于0|Max target sequences must be greater than 0: {self.max_target_seqs}")

        if self.word_size < 1:
            errors.append(f"词大小必须大于0|Word size must be greater than 0: {self.word_size}")

        # 验证比对可视化参数|Validate alignment visualization parameters
        if self.alignment_width < 1:
            errors.append(f"比对宽度必须大于0|Alignment width must be greater than 0: {self.alignment_width}")

        if self.alignment_max_per_sample < 1:
            errors.append(f"每个样品最大比对数必须大于0|Maximum alignments per sample must be greater than 0: {self.alignment_max_per_sample}")

        # 验证正则表达式|Validate regex pattern
        try:
            re.compile(self.sample_name_pattern)
        except re.error as e:
            errors.append(f"样品名称正则表达式无效|Invalid sample name regex: {e}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def __str__(self) -> str:
        """配置的字符串表示|String representation of configuration"""
        return f"""BLASTConfig:
    输入文件|Input: {self.input}
    目标数据库|Target Database: {self.reference}
    输出目录|Output Directory: {self.output}
    BLAST类型|BLAST Type: {self.blast_type}
    线程数|Threads: {self.threads}
    E-value阈值|E-value: {self.evalue}
    最大目标序列数|Max Target Seqs: {self.max_target_seqs}
    最小相似度|Min Identity: {self.min_identity}%
    最小覆盖度|Min Coverage: {self.min_coverage}%
    自动检测样品|Auto Detect Samples: {self.auto_detect_samples}
    词大小|Word Size: {self.word_size}
    样品映射文件|Sample Map File: {self.sample_map_file}
    保持中间文件|Keep Intermediate: {self.keep_intermediate}
"""