"""
 序列格式检测与转换模块|Sequence Format Detection and Conversion Module

 输入格式自动检测(phylip/fasta/vcf); VCF输入时复用vcf2phylip模块转换为
 PHYLIP矩阵, 其余格式直通(RAxML原生支持)
 |Auto-detect input format (phylip/fasta/vcf); VCF input is converted to a
 PHYLIP matrix via the vcf2phylip module, other formats pass through
 (RAxML-native)
"""

import gzip
from pathlib import Path

from ..vcf2phylip import VCFConverter


def detect_input_format(sequence_file: str) -> str:
    """检测输入文件格式|Detect input file format

    先看后缀(.vcf/.vcf.gz), 再看文件首行(##fileformat=VCF, >, "int int")
    |Check suffix (.vcf/.vcf.gz) first, then the first line (##fileformat=VCF, >, "int int")

    Returns:
        'vcf' | 'fasta' | 'phylip'

    Raises:
        ValueError: 无法识别格式|Unrecognized format
    """
    name_lower = Path(sequence_file).name.lower()
    if name_lower.endswith('.vcf') or name_lower.endswith('.vcf.gz'):
        return 'vcf'

    opener = gzip.open if name_lower.endswith('.gz') else open
    with opener(sequence_file, 'rt') as f:
        first_line = f.readline().strip()

    if first_line.startswith('##fileformat=VCF') or first_line.startswith('#CHROM'):
        return 'vcf'
    if first_line.startswith('>'):
        return 'fasta'

    parts = first_line.split()
    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
        return 'phylip'

    raise ValueError(
        f"无法识别输入文件格式|Unrecognized input file format: {sequence_file}"
        " (支持|supported: PHYLIP/FASTA/VCF)"
    )


def count_vcf_samples(vcf_file: str) -> int:
    """从VCF头部统计样本数|Count samples from the VCF header"""
    opener = gzip.open if vcf_file.lower().endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return max(0, len(line.rstrip('\n').split('\t')) - 9)
    return 0


class FormatConverter:
    """输入格式转换器|Input Format Converter"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def prepare_alignment(self) -> str:
        """检测格式并按需转换, 返回可直接喂给RAxML的对齐文件路径
        |Detect format and convert when needed; return the RAxML-ready alignment path"""
        input_format = self.config.input_format
        if input_format in (None, 'auto'):
            input_format = detect_input_format(self.config.sequence_file)
        self.logger.info(f" 检测到输入格式|Detected input format: {input_format}")

        if input_format == 'vcf':
            return self._convert_vcf_to_phylip()
        if input_format in ('fasta', 'phylip'):
            if input_format == 'fasta':
                self.logger.info(" FASTA直通RAxML原生输入|FASTA passes through as RAxML-native input")
            return self.config.sequence_file

        raise ValueError(f" 不支持的输入格式|Unsupported input format: {input_format}")

    def _convert_vcf_to_phylip(self) -> str:
        """调用vcf2phylip模块完成VCF到PHYLIP转换|Convert VCF to PHYLIP via the vcf2phylip module"""
        target = self._converted_phylip_path()

        # 断点续传: 已转换过则跳过|Checkpoint resume: skip if already converted
        if Path(target).exists():
            self.logger.info(f" 跳过已完成转换|Skipping completed conversion: {target}")
            return target

        self.logger.info(
            " 开始VCF到PHYLIP转换(复用vcf2phylip模块)"
            "|Starting VCF to PHYLIP conversion (via vcf2phylip module)"
        )

        converter = VCFConverter(
            input_file=self.config.sequence_file,
            output_dir=str(self.config.output_path),
            output_prefix=f"{self.config.output_name}_converted",
            min_samples_locus=self._effective_min_samples(),
            outgroup=self.config.outgroup or "",
            resolve_IUPAC=self.config.resolve_iupac,
            threads=self.config.threads,
        )
        try:
            converter.run_conversion()
        except SystemExit as e:
            # vcf2phylip失败路径走sys.exit(1), 转为正常异常上抛
            # |vcf2phylip failures call sys.exit(1); re-raise as a regular exception
            raise RuntimeError(
                f"VCF到PHYLIP转换失败|VCF to PHYLIP conversion failed (exit={e.code})"
            ) from None

        if not Path(target).exists():
            raise RuntimeError(
                f"VCF到PHYLIP转换失败: 未生成目标文件|VCF to PHYLIP conversion failed,"
                f" output missing: {target}"
            )

        self.logger.info(f" VCF转换成功|VCF conversion successful: {target}")
        return target

    def _effective_min_samples(self) -> int:
        """与vcf2phylip同规则收缩min: 超过样本数时降到样本数(影响输出文件名)
        |Shrink min like vcf2phylip: capped at sample count (affects output filename)"""
        num_samples = count_vcf_samples(self.config.sequence_file)
        effective = min(self.config.min_samples_locus, num_samples)
        if effective != self.config.min_samples_locus:
            self.logger.info(
                f" 最小样本数已调整为|Min samples adjusted to: {effective}"
            )
        return effective

    def _converted_phylip_path(self) -> str:
        """转换产物路径(vcf2phylip固定追加_min{{m}}后缀)
        |Converted output path (vcf2phylip always appends _min{{m}})"""
        num_samples = count_vcf_samples(self.config.sequence_file)
        effective_min = min(self.config.min_samples_locus, num_samples)
        return str(
            self.config.output_path
            / f"{self.config.output_name}_converted_min{effective_min}.phy"
        )
