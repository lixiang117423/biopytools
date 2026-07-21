"""
BLAST分析配置模块|BLAST Analysis Configuration Module
标准化BLAST配置类,遵循参数命名规范|Standardized BLAST configuration class following naming conventions
"""

import os
import re
import logging
from typing import Optional, Tuple
from Bio import SeqIO
from ..core.config import BaseConfig
from ..common.paths import get_tool_path, expand_path


class BLASTConfig(BaseConfig):
    """BLAST分析配置类|BLAST Analysis Configuration Class"""

    def __init__(
        self,
        input: Optional[str] = None,
        output_dir: Optional[str] = None,
        reference: Optional[str] = None,
        sample_name: Optional[str] = None,
        sample_map_file: Optional[str] = None,
        threads: int = 12,
        log_level: str = "INFO",
        log_file: Optional[str] = None,
        force: bool = False,
        dry_run: bool = False,
        # BLAST特定参数|BLAST-specific parameters
        blast_type: str = None,
        evalue: float = 1e-5,
        max_target_seqs: int = 10,
        word_size: Optional[int] = None,
        target_db_type: str = None,
        min_identity: float = 70.0,
        min_coverage: float = 50.0,
        high_quality_evalue: float = 1e-10,
        input_suffix: str = "*.fa",
        auto_detect_samples: bool = True,
        sample_name_pattern: str = r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
        # BLAST工具路径|BLAST tool paths
        makeblastdb_path: str = None,
        blastn_path: str = None,
        blastp_path: str = None,
        blastx_path: str = None,
        tblastn_path: str = None,
        tblastx_path: str = None,
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
            output_dir: 输出目录路径|Output directory path
            reference: 目标数据库文件路径|Target database file path
            sample_name: 单文件输入时的样品名称|Sample name for single-file input
            sample_map_file: 样品映射文件|Sample mapping file
            threads: 线程数|Number of threads
            log_level: 日志级别|Log level
            log_file: 日志文件路径|Log file path
            force: 强制覆盖已存在文件|Force overwrite existing files
            dry_run: 模拟运行不执行|Dry run without execution
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
        self.output = output_dir
        self.reference = reference
        self.sample_name = sample_name
        self.sample_map_file = sample_map_file
        self.threads = self.validate_threads(threads)

        # 日志参数|Logging parameters
        self.log_level = log_level
        self.log_file = log_file
        self.force = force
        self.dry_run = dry_run

        # BLAST特定参数|BLAST-specific parameters
        self.blast_type = blast_type
        self.evalue = self._validate_evalue(evalue)
        self.max_target_seqs = max_target_seqs
        self.word_size = word_size
        self.target_db_type = target_db_type
        self.min_identity = self.validate_quality(min_identity, 0, 100)
        self.min_coverage = self.validate_quality(min_coverage, 0, 100)
        self.high_quality_evalue = self._validate_evalue(high_quality_evalue)
        self.input_suffix = input_suffix
        self.auto_detect_samples = auto_detect_samples
        self.sample_name_pattern = sample_name_pattern

        # 工具路径参数|Tool path parameters
        blast_env = '~/miniforge3/envs/Blast_v.2.16.0/bin'
        self.makeblastdb_path = makeblastdb_path or get_tool_path('makeblastdb', f'{blast_env}/makeblastdb', 'MAKEBLASTDB_PATH')
        self.blastn_path = blastn_path or get_tool_path('blastn', f'{blast_env}/blastn', 'BLASTN_PATH')
        self.blastp_path = blastp_path or get_tool_path('blastp', f'{blast_env}/blastp', 'BLASTP_PATH')
        self.blastx_path = blastx_path or get_tool_path('blastx', f'{blast_env}/blastx', 'BLASTX_PATH')
        self.tblastn_path = tblastn_path or get_tool_path('tblastn', f'{blast_env}/tblastn', 'TBLASTN_PATH')
        self.tblastx_path = tblastx_path or get_tool_path('tblastx', f'{blast_env}/tblastx', 'TBLASTX_PATH')

        # 比对可视化参数|Alignment visualization parameters
        self.alignment_output = self._validate_alignment_output(alignment_output)
        self.alignment_width = alignment_width
        self.alignment_min_identity = self.validate_quality(alignment_min_identity, 0, 100)
        self.alignment_min_coverage = self.validate_quality(alignment_min_coverage, 0, 100)
        self.alignment_max_per_sample = alignment_max_per_sample
        self.html_theme = self._validate_html_theme(html_theme)

        # 设置默认输出路径|Set default output path
        if self.output is None:
            self.output = "./blast_output"

        # 分层输出目录(§12.2)|Layered output directories
        self.pipeline_info_dir = os.path.join(self.output, "00_pipeline_info")
        self.db_dir = os.path.join(self.output, "01_database")
        self.blast_dir = os.path.join(self.output, "02_blast")
        self.alignments_dir = os.path.join(self.output, "03_alignments")
        self.logs_dir = os.path.join(self.output, "99_logs")

        # 自动检测BLAST类型(仅当用户未指定时)|Auto-detect BLAST type (only when user didn't specify)
        self._auto_detect_blast_type()

        # 验证BLAST类型|Validate BLAST type
        self._validate_blast_type()

        # 根据BLAST类型自动设置数据库类型(仅在用户未指定时)|Auto-set db type (only when user did not specify)
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
        """比对可视化输出子目录名(供text/html generator拼接)|Alignment output subdir name"""
        return "03_alignments"

    def _validate_evalue(self, evalue: float) -> float:
        """验证E-value(允许0表示完美匹配,不设上限)|Validate E-value (0 allowed for exact match, no upper limit)"""
        if evalue < 0:
            raise ValueError(f"E-value阈值不能为负数|E-value must be non-negative: {evalue}")
        return evalue

    def _detect_sequence_type(self, fasta_file: str, check_sequences: int = 5) -> Tuple[str, float]:
        """
        自动检测FASTA文件的序列类型|Auto-detect sequence type in FASTA file

        Args:
            fasta_file: FASTA文件路径|FASTA file path
            check_sequences: 检查的序列数量|Number of sequences to check

        Returns:
            Tuple[str, float]: ('dna'/'protein'/'unknown', 置信度|confidence)
        """
        logger = logging.getLogger(__name__)

        try:
            sequences = []
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(str(record.seq).upper())
                if len(sequences) >= check_sequences:
                    break

            if not sequences:
                logger.error(f"FASTA文件中没有序列|No sequences found in FASTA file: {fasta_file}")
                return "unknown", 0.0

            total_chars = 0
            dna_chars = 0
            protein_chars = 0

            for seq in sequences:
                for char in seq:
                    total_chars += 1
                    if char in "ATCGN":
                        dna_chars += 1
                    elif char in "DEFGHIKLMNPQRSVWY*":
                        protein_chars += 1
                    else:
                        # 未知字符按蛋白处理(蛋白字母表更宽)|Treat unknown as protein (larger alphabet)
                        protein_chars += 1

            if total_chars == 0:
                return "unknown", 0.0

            dna_ratio = dna_chars / total_chars
            protein_ratio = protein_chars / total_chars

            if dna_ratio > 0.90:
                return "dna", dna_ratio
            elif protein_ratio > 0.30:
                return "protein", protein_ratio
            else:
                logger.warning(f"无法确定序列类型|Cannot determine sequence type (DNA: {dna_ratio:.2f}, Protein: {protein_ratio:.2f}): {fasta_file}")
                return "unknown", 0.0

        except Exception as e:
            logger.error(f"序列类型检测失败|Sequence type detection failed: {e}")
            return "unknown", 0.0

    def _collect_fasta_files(self, path: str) -> list:
        """收集FASTA文件列表|Collect FASTA file list

        当path是文件时返回单元素列表,目录时返回所有匹配文件
        """
        import glob as glob_mod
        if os.path.isfile(path):
            return [path]
        if os.path.isdir(path):
            files = []
            for ext in ('*.fa', '*.fasta', '*.fna', '*.pep', '*.faa'):
                files.extend(glob_mod.glob(os.path.join(path, ext)))
            return files
        return [path]

    def _detect_type_from_files(self, files: list, label: str) -> str:
        """对多个文件逐一检测序列类型,用多数投票决定|Detect sequence type from multiple files by majority vote"""
        logger = logging.getLogger(__name__)
        from collections import Counter
        counts = Counter()
        for f in files:
            seq_type, _ = self._detect_sequence_type(f)
            if seq_type != "unknown":
                counts[seq_type] += 1

        if not counts:
            logger.warning(f"无法确定{label}序列类型|Cannot determine {label} sequence type")
            return "unknown"

        most_common_type, most_common_count = counts.most_common(1)[0]
        if len(counts) > 1:
            logger.warning(f"{label}文件中存在混合序列类型|Mixed sequence types in {label} files: "
                           f"{dict(counts)}, 使用多数类型|using majority type: {most_common_type}")
        return most_common_type

    def _auto_detect_blast_type(self):
        """根据输入文件序列类型自动推断BLAST类型|Auto-infer BLAST type from input file sequence types"""
        logger = logging.getLogger(__name__)

        if self.blast_type is not None:
            return

        if not self.input or not self.reference:
            return

        query_files = self._collect_fasta_files(self.input)
        ref_files = self._collect_fasta_files(self.reference)
        query_type = self._detect_type_from_files(query_files, "查询|query")
        ref_type = self._detect_type_from_files(ref_files, "参考|reference")

        type_map = {
            ("dna", "dna"): "blastn",
            ("protein", "protein"): "blastp",
            ("dna", "protein"): "blastx",
            ("protein", "dna"): "tblastn",
        }

        detected = type_map.get((query_type, ref_type))

        if detected:
            self.blast_type = detected
            logger.info(f"自动检测BLAST类型|Auto-detected BLAST type: {detected} "
                        f"(查询|query: {query_type}, 参考|reference: {ref_type})")
        else:
            self.blast_type = "blastn"
            logger.warning(f"无法自动检测BLAST类型,使用默认值blastn|Cannot auto-detect BLAST type, using default blastn "
                           f"(查询|query: {query_type}, 参考|reference: {ref_type})")

    def _validate_blast_type(self):
        """验证BLAST类型|Validate BLAST type"""
        valid_blast_types = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
        if self.blast_type not in valid_blast_types:
            raise ValueError(f"无效的BLAST类型|Invalid BLAST type: {self.blast_type}")

    def _set_default_word_size(self):
        """设置默认word_size|Set default word_size"""
        if self.word_size is None:
            # 根据BLAST类型设置默认word_size|Set default word_size by BLAST type
            default_word_sizes = {
                'blastn': 11,
                'tblastx': 11,
                'blastp': 3,
                'blastx': 3,
                'tblastn': 3
            }
            self.word_size = default_word_sizes.get(self.blast_type, 11)

    def _auto_set_db_type(self):
        """根据BLAST类型自动设置数据库类型(仅在用户未指定时)|Auto-set db type only when user did not specify"""
        # 用户显式指定target_db_type时尊重其选择,不再覆盖|Respect user-specified target_db_type
        if self.target_db_type is not None:
            return

        # blast_type到数据库类型的映射|blast_type to dbtype mapping
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

        self.target_db_type = blast_to_db.get(self.blast_type, 'nucl')

    def get_log_level(self) -> int:
        """获取日志级别(显式log_level优先,否则按verbose/quiet推导)|Get log level (explicit log_level wins, else from verbose/quiet)"""
        # CLI的--log-level或main.py由-v推导的log_level优先于基类的verbose/quiet逻辑
        # |Explicit log_level (from --log-level or -v derivation) takes precedence over base verbose/quiet
        if self.log_level:
            return getattr(logging, str(self.log_level).upper(), logging.INFO)
        return super().get_log_level()

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

    def get_target_db_path(self) -> str:
        """
        获取目标数据库路径|Get target database path

        Returns:
            str: 数据库路径|Database path
        """
        base_name = os.path.splitext(os.path.basename(self.reference))[0]
        return os.path.join(self.db_dir, f"{base_name}.db")

    def get_summary_output_path(self) -> str:
        """
        获取汇总输出路径|Get summary output path

        Returns:
            str: 输出路径|Output path
        """
        return os.path.join(self.blast_dir, "blast_summary_results.tsv")

    def get_alignment_output_path(self, output_type: str = "html") -> str:
        """
        获取比对可视化输出路径|Get alignment visualization output path

        Args:
            output_type: 输出类型 (html, text)|Output type (html, text)

        Returns:
            str: 输出路径|Output path
        """
        if output_type == "html":
            return os.path.join(self.alignments_dir, "html", "index.html")
        elif output_type == "text":
            return os.path.join(self.alignments_dir, "text", "all_samples_alignments.txt")
        else:
            raise ValueError(f"不支持的输出类型|Unsupported output type: {output_type}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        # 调用父类验证|Call parent validation
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
        if self.evalue < 0:
            errors.append(f"E-value阈值不能为负数|E-value must be non-negative: {self.evalue}")

        if self.high_quality_evalue < 0:
            errors.append(f"高质量E-value阈值不能为负数|High quality E-value must be non-negative: {self.high_quality_evalue}")

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
"""
