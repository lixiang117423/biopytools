"""
转录本从头组装核心逻辑模块|Transcript De Novo Assembly Core Logic Module

流程|Pipeline:
    Step 1: hisat2-build 构建基因组索引|Build genome index
    Step 2: hisat2 --dta 比对|Alignment with --dta
    Step 3: samtools sort 格式转换|Format conversion and sorting
    Step 4: stringtie 逐样本组装|Per-sample de novo assembly
    Step 5: stringtie --merge 合并GTF|Merge GTFs
    Step 6: gffread 提取转录本序列|Extract transcript sequences
"""

import os
import glob
import re
from typing import List, Dict, Optional

from .utils import CommandRunner, build_conda_command


class HISAT2Indexer:
    """HISAT2基因组索引构建器|HISAT2 Genome Index Builder"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_index(self) -> str:
        """
        构建HISAT2基因组索引|Build HISAT2 genome index

        Returns:
            索引前缀路径|Index prefix path
        """
        genome_path = self.config.genome_file
        threads = self.config.threads

        genome_dir = os.path.dirname(genome_path)
        genome_name = os.path.splitext(os.path.basename(genome_path))[0]
        index_prefix = os.path.join(self.config.output_dir, "01_hisat2_index", f"{genome_name}")

        # 断点续传：检查索引是否已存在|Checkpoint: check if index exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2索引已存在，跳过|HISAT2 index already exists, skipping: {index_prefix}")
            return index_prefix

        # 确保输出目录存在|Ensure output directory exists
        os.makedirs(os.path.dirname(index_prefix), exist_ok=True)

        # 构建命令|Build command
        hisat2_args = [
            'hisat2-build',
            '-p', str(threads),
            genome_path,
            index_prefix
        ]
        cmd_list = build_conda_command('hisat2-build', hisat2_args[1:])
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(cmd, "构建HISAT2基因组索引|Building HISAT2 genome index")

        if not success:
            raise RuntimeError("HISAT2索引构建失败|HISAT2 index building failed")

        self.logger.info(f"HISAT2索引构建完成|HISAT2 index building completed: {index_prefix}")
        return index_prefix


class HISAT2Aligner:
    """HISAT2序列比对器|HISAT2 Aligner"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def align(self, index_prefix: str, fastq1: str, fastq2: str,
              output_sam: str) -> bool:
        """
        运行HISAT2比对（开启--dta）|Run HISAT2 alignment with --dta

        Args:
            index_prefix: HISAT2索引前缀|HISAT2 index prefix
            fastq1: Read1文件路径|Read1 file path
            fastq2: Read2文件路径|Read2 file path
            output_sam: 输出SAM文件路径|Output SAM file path

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads
        timeout = self.config.sample_timeout

        os.makedirs(os.path.dirname(output_sam), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_sam) and os.path.getsize(output_sam) > 0:
            if not self.config.force:
                self.logger.info(f"比对结果已存在，跳过|Alignment result already exists, skipping: {output_sam}")
                return True

        # hisat2 --dta 优化StringTie组装|hisat2 --dta optimizes StringTie assembly
        hisat2_args = [
            'hisat2',
            '-x', index_prefix,
            '-1', fastq1,
            '-2', fastq2,
            '-p', str(threads),
            '--dta',
            '-S', output_sam
        ]
        cmd_list = build_conda_command('hisat2', hisat2_args[1:])
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(
            cmd,
            f"HISAT2比对|HISAT2 alignment -> {output_sam}",
            timeout=timeout
        )

        return success


class SamtoolsProcessor:
    """Samtools格式转换器|Samtools Format Converter"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def sort_and_index(self, input_sam: str, output_bam: str) -> bool:
        """
        SAM转排序BAM并建索引|Convert SAM to sorted BAM and index

        Args:
            input_sam: 输入SAM文件|Input SAM file
            output_bam: 输出BAM文件|Output BAM file

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads
        timeout = self.config.sample_timeout

        os.makedirs(os.path.dirname(output_bam), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_bam) and os.path.exists(f"{output_bam}.bai"):
            if not self.config.force:
                self.logger.info(f"排序BAM已存在，跳过|Sorted BAM already exists, skipping: {output_bam}")
                return True

        # samtools sort
        sort_args = [
            'samtools', 'sort',
            '-@', str(threads),
            '-o', output_bam,
            input_sam
        ]
        sort_cmd_list = build_conda_command('samtools', sort_args[1:])
        sort_cmd = ' '.join(sort_cmd_list)

        success = self.cmd_runner.run(
            sort_cmd,
            f"samtools排序|samtools sort -> {output_bam}",
            timeout=timeout
        )
        if not success:
            return False

        # samtools index
        index_args = [
            'samtools', 'index',
            '-@', str(threads),
            output_bam
        ]
        index_cmd_list = build_conda_command('samtools', index_args[1:])
        index_cmd = ' '.join(index_cmd_list)

        success = self.cmd_runner.run(
            index_cmd,
            f"samtools索引|samtools index -> {output_bam}.bai",
            timeout=timeout
        )

        return success


class StringTieAssembler:
    """StringTie逐样本组装器|StringTie Per-Sample Assembler"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def assemble(self, bam_file: str, output_gtf: str) -> bool:
        """
        逐样本从头组装转录本|De novo transcript assembly per sample (no GTF reference)

        Args:
            bam_file: 排序后的BAM文件|Sorted BAM file
            output_gtf: 输出GTF文件|Output GTF file

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads
        timeout = self.config.sample_timeout

        os.makedirs(os.path.dirname(output_gtf), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_gtf) and os.path.getsize(output_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"样本GTF已存在，跳过|Sample GTF already exists, skipping: {output_gtf}")
                return True

        # 不使用-G参数，进行从头组装|No -G parameter, de novo assembly
        stringtie_args = [
            'stringtie',
            bam_file,
            '-p', str(threads),
            '-o', output_gtf,
            '-l', os.path.splitext(os.path.basename(output_gtf))[0]
        ]
        cmd_list = build_conda_command('stringtie', stringtie_args[1:])
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(
            cmd,
            f"StringTie从头组装|StringTie de novo assembly -> {output_gtf}",
            timeout=timeout
        )

        return success


class StringTieMerger:
    """StringTie GTF合并器|StringTie GTF Merger"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def merge(self, gtf_list_file: str, genome_file: str, output_merged_gtf: str) -> bool:
        """
        合并所有样本的GTF|Merge all sample GTFs

        Args:
            gtf_list_file: GTF列表文件（每行一个GTF路径）|GTF list file (one path per line)
            genome_file: 基因组FASTA文件|Genome FASTA file
            output_merged_gtf: 合并后的GTF文件|Merged GTF file

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads

        os.makedirs(os.path.dirname(output_merged_gtf), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_merged_gtf) and os.path.getsize(output_merged_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"合并GTF已存在，跳过|Merged GTF already exists, skipping: {output_merged_gtf}")
                return True

        # stringtie --merge 需要基因组序列作为参考|stringtie --merge needs genome as reference
        stringtie_args = [
            'stringtie',
            '--merge',
            '-p', str(threads),
            '-G', genome_file,
            '-o', output_merged_gtf,
            gtf_list_file
        ]
        cmd_list = build_conda_command('stringtie', stringtie_args[1:])
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(
            cmd,
            f"StringTie合并GTF|StringTie merge GTF -> {output_merged_gtf}"
        )

        return success


class GFFReadExtractor:
    """GFFRead转录本序列提取器|GFFRead Transcript Sequence Extractor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_transcripts(self, merged_gtf: str, genome_file: str,
                            output_fa: str) -> bool:
        """
        从merged GTF和基因组提取转录本cDNA序列|Extract transcript cDNA sequences from merged GTF and genome

        Args:
            merged_gtf: 合并后的GTF文件|Merged GTF file
            genome_file: 基因组FASTA文件|Genome FASTA file
            output_fa: 输出转录本FASTA文件|Output transcript FASTA file

        Returns:
            是否成功|Whether successful
        """
        os.makedirs(os.path.dirname(output_fa), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_fa) and os.path.getsize(output_fa) > 0:
            if not self.config.force:
                self.logger.info(f"转录本序列已存在，跳过|Transcript sequences already exists, skipping: {output_fa}")
                return True

        gffread_args = [
            'gffread',
            merged_gtf,
            '-g', genome_file,
            '-w', output_fa
        ]
        cmd_list = build_conda_command('gffread', gffread_args[1:])
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(
            cmd,
            f"gffread提取转录本序列|gffread extract transcript sequences -> {output_fa}"
        )

        return success


class SampleParser:
    """样本解析器|Sample Parser"""

    def __init__(self, logger):
        self.logger = logger

    def parse_input_samples(self, input_dir: str, fastq_pattern: str = None) -> List[Dict]:
        """
        从输入目录解析样本信息|Parse sample information from input directory

        Args:
            input_dir: 输入目录|Input directory
            fastq_pattern: FASTQ文件命名模式|FASTQ file naming pattern

        Returns:
            样本信息列表|Sample information list
        """
        samples = []

        if fastq_pattern:
            samples = self._parse_with_pattern(input_dir, fastq_pattern)
        else:
            samples = self._parse_with_default_patterns(input_dir)

        return samples

    def _parse_with_pattern(self, input_dir: str, fastq_pattern: str) -> List[Dict]:
        """使用指定模式解析样本|Parse samples with specified pattern"""
        samples = []

        if '*' not in fastq_pattern:
            self.logger.error(f"文件模式必须包含*作为样本名占位符|File pattern must contain * as sample name placeholder: {fastq_pattern}")
            return samples

        parts = fastq_pattern.split('*')
        if len(parts) != 2:
            self.logger.error(f"文件模式只能包含一个*占位符|File pattern can only contain one * placeholder: {fastq_pattern}")
            return samples

        prefix = parts[0]
        suffix = parts[1]

        # 检测read标识符|Detect read indicators
        read_pairs = [
            ("_1", "_2"), ("R1", "R2"), (".1", ".2"),
            ("_f1", "_r2"), ("_f1", "_f2"),
            ("_F1", "_F2"),
        ]
        read1_indicator = None
        read2_indicator = None

        for r1, r2 in read_pairs:
            if r1 in suffix:
                read1_indicator = r1
                read2_indicator = r2
                break

        if not read1_indicator:
            self.logger.error(
                f"无法识别read标识符|Cannot identify read indicator in pattern '{fastq_pattern}'")
            return samples

        search_pattern = os.path.join(input_dir, fastq_pattern)
        fastq_files = sorted(glob.glob(search_pattern))

        self.logger.info(f"搜索模式|Search pattern: {search_pattern}")
        self.logger.info(f"找到 {len(fastq_files)} 个read1文件|Found {len(fastq_files)} read1 files")

        for fq1 in fastq_files:
            basename = os.path.basename(fq1)

            # 提取样本名|Extract sample name
            sample_name = basename
            if prefix:
                sample_name = sample_name.replace(prefix, "", 1)
            if suffix:
                sample_name = sample_name.replace(suffix, "", 1)

            # 构建read2文件名|Build read2 filename
            read2_suffix = suffix.replace(read1_indicator, read2_indicator)
            read2_filename = prefix + sample_name + read2_suffix
            fq2 = os.path.join(input_dir, read2_filename)

            if os.path.exists(fq2):
                samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
            else:
                self.logger.warning(f"找不到配对的read2文件|Cannot find paired read2 file: {fq2}")

        return samples

    def _parse_with_default_patterns(self, input_dir: str) -> List[Dict]:
        """使用默认模式解析样本|Parse samples with default patterns"""
        patterns = [
            "*_1.fq.gz", "*_R1.fq.gz", "*.R1.fastq.gz",
            "*_1.fastq.gz", "*.1.fq.gz", "*_R1.fq",
            "*_1.clean.fq.gz", "*_f1.fq.gz",
        ]

        for pattern in patterns:
            samples = self._parse_with_pattern(input_dir, pattern)
            if samples:
                self.logger.info(f"使用模式|Using pattern: {pattern}")
                return samples

        return []
