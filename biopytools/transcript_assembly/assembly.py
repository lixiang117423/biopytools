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
import statistics
import subprocess
from typing import List, Dict, Optional

from .utils import CommandRunner, build_conda_command


READ_LENGTH_THRESHOLD = 500  # 短/长读分界 bp|short/long threshold bp
READ_SAMPLE_COUNT = 1000     # 采样 read 数|sample read count


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
        hisat2_args = ['-p', str(threads), genome_path, index_prefix]
        cmd_list = build_conda_command(self.config.hisat2_build_bin, hisat2_args)
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
        hisat2_args = ['-x', index_prefix, '-1', fastq1, '-2', fastq2,
                       '-p', str(threads), '--dta', '-S', output_sam]
        cmd_list = build_conda_command(self.config.hisat2_bin, hisat2_args)
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
        sort_args = ['sort', '-@', str(threads), '-o', output_bam, input_sam]
        sort_cmd_list = build_conda_command(self.config.samtools_bin, sort_args)
        sort_cmd = ' '.join(sort_cmd_list)

        success = self.cmd_runner.run(
            sort_cmd,
            f"samtools排序|samtools sort -> {output_bam}",
            timeout=timeout
        )
        if not success:
            return False

        # samtools index
        index_args = ['index', '-@', str(threads), output_bam]
        index_cmd_list = build_conda_command(self.config.samtools_bin, index_args)
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

    def assemble(self, bam_file: str, output_gtf: str, read_type: str,
                 guide_gff: str = None) -> bool:
        """
        逐样本组装|Per-sample assembly

        Args:
            bam_file: 排序后的BAM|sorted BAM
            output_gtf: 输出GTF|output GTF
            read_type: short|long (long → -L)
            guide_gff: 参考注释(可选 → -G)|reference annotation (optional → -G)

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads
        timeout = self.config.sample_timeout

        os.makedirs(os.path.dirname(output_gtf), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_gtf) and os.path.getsize(output_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"样本GTF已存在，跳过|Sample GTF exists, skipping: {output_gtf}")
                return True

        stringtie_args = [bam_file, '-p', str(threads), '-o', output_gtf,
                          '-l', os.path.splitext(os.path.basename(output_gtf))[0]]
        if read_type == 'long':
            stringtie_args.append('-L')
        if guide_gff:
            stringtie_args.extend(['-G', guide_gff])

        cmd_list = build_conda_command(self.config.stringtie_bin, stringtie_args)
        cmd = ' '.join(cmd_list)

        success = self.cmd_runner.run(
            cmd,
            f"StringTie组装|StringTie assembly ({read_type}"
            f"{' guided' if guide_gff else ''}) -> {output_gtf}",
            timeout=timeout
        )

        return success


class StringTieMerger:
    """StringTie GTF合并器|StringTie GTF Merger"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def merge(self, gtf_list_file: str, guide_gff: str,
              output_merged_gtf: str) -> bool:
        """
        合并GTF|Merge GTFs

        Args:
            gtf_list_file: GTF列表(每行一个路径)|GTF list (one path per line)
            guide_gff: 参考注释(可选 → -G)|reference annotation (optional → -G)
            output_merged_gtf: 合并后GTF|merged GTF

        Returns:
            是否成功|Whether successful
        """
        threads = self.config.threads
        os.makedirs(os.path.dirname(output_merged_gtf), exist_ok=True)

        # 断点续传|Checkpoint
        if os.path.exists(output_merged_gtf) and os.path.getsize(output_merged_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"合并GTF已存在，跳过|Merged GTF exists, skipping: {output_merged_gtf}")
                return True

        # --merge: -G 只接注释,不接 genome.fa (修复原 bug)|-G takes annotation not genome (bugfix)
        stringtie_args = ['--merge', '-p', str(threads), '-o', output_merged_gtf]
        if guide_gff:
            stringtie_args.extend(['-G', guide_gff])
        stringtie_args.append(gtf_list_file)

        cmd_list = build_conda_command(self.config.stringtie_bin, stringtie_args)
        cmd = ' '.join(cmd_list)

        return self.cmd_runner.run(
            cmd, f"StringTie合并GTF|StringTie merge -> {output_merged_gtf}")


class GFF3Converter:
    """GTF → GFF3 转换器(主输出)|GTF to GFF3 converter (primary output)"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def convert(self, gtf_file: str, output_gff3: str) -> bool:
        """GTF转GFF3|Convert GTF to GFF3 (no genome needed)"""
        os.makedirs(os.path.dirname(output_gff3), exist_ok=True)

        if os.path.exists(output_gff3) and os.path.getsize(output_gff3) > 0:
            if not self.config.force:
                self.logger.info(f"GFF3已存在，跳过|GFF3 exists, skipping: {output_gff3}")
                return True

        gffread_args = [gtf_file, '-o', output_gff3, '-F']
        cmd_list = build_conda_command(self.config.gffread_bin, gffread_args)
        cmd = ' '.join(cmd_list)

        return self.cmd_runner.run(
            cmd, f"gffread转GFF3|gffread GTF→GFF3 -> {output_gff3}")


class TranscriptExtractor:
    """cDNA 提取器(可选)|cDNA extractor (optional)"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract(self, gtf_file: str, genome_file: str, output_fa: str) -> bool:
        """提取转录本cDNA|Extract transcript cDNA sequences"""
        os.makedirs(os.path.dirname(output_fa), exist_ok=True)

        if os.path.exists(output_fa) and os.path.getsize(output_fa) > 0:
            if not self.config.force:
                self.logger.info(f"转录本序列已存在，跳过|Transcripts exist, skipping: {output_fa}")
                return True

        gffread_args = [gtf_file, '-g', genome_file, '-w', output_fa]
        cmd_list = build_conda_command(self.config.gffread_bin, gffread_args)
        cmd = ' '.join(cmd_list)

        return self.cmd_runner.run(
            cmd, f"gffread提取cDNA|gffread extract cDNA -> {output_fa}")


class TransDecoderRunner:
    """TransDecoder CDS 预测器|TransDecoder CDS predictor (transcripts→gene/mRNA/CDS)"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def predict(self, transcripts_fa: str, merged_gtf: str,
                output_dir: str) -> bool:
        """
        预测 CDS 并输出基因组坐标 GFF3|Predict CDS, output genome-coord GFF3

        Flow: merged.gtf→alignment.gff3 + LongOrfs + Predict + ORF→genome
        (基础版,无 blastp/hmmscan 同源 DB|basic, no homology DBs)

        Args:
            transcripts_fa: cDNA FASTA(gffread -w 产物)|cDNA FASTA
            merged_gtf: 合并GTF(生成 alignment.gff3)|merged GTF
            output_dir: TransDecoder 输出目录|output dir
        """
        os.makedirs(output_dir, exist_ok=True)
        base = os.path.basename(transcripts_fa)  # transcripts.fa
        genome_gff3 = os.path.join(output_dir, f"{base}.transdecoder.genome.gff3")

        # 断点续传:最终基因组坐标 GFF3 存在则跳过|Checkpoint
        if os.path.exists(genome_gff3) and os.path.getsize(genome_gff3) > 0:
            if not self.config.force:
                self.logger.info(f"TransDecoder已完成，跳过|TransDecoder done, skipping: {genome_gff3}")
                return True

        timeout = self.config.sample_timeout

        # 1. merged.gtf → alignment.gff3(转录本→基因组比对 GFF3)|gtf→alignment.gff3
        align_gff3 = os.path.join(output_dir, "alignment.gff3")
        cmd = ' '.join(build_conda_command(
            self.config.transdecoder_gtf2gff_bin, [merged_gtf])) + f" > {align_gff3}"
        if not self.cmd_runner.run(cmd, f"TransDecoder: GTF→alignment GFF3 -> {align_gff3}"):
            return False

        # 2. LongOrfs(找长 ORF)|find long ORFs
        cmd = ' '.join(build_conda_command(
            self.config.transdecoder_longorfs_bin, ['-t', transcripts_fa, '-O', output_dir]))
        if not self.cmd_runner.run(cmd, f"TransDecoder.LongOrfs -> {output_dir}", timeout=timeout):
            return False

        # 3. Predict(基础版,无同源 DB)|Predict ORFs (basic, no homology DBs)
        cmd = ' '.join(build_conda_command(
            self.config.transdecoder_predict_bin, ['-t', transcripts_fa, '-O', output_dir]))
        if not self.cmd_runner.run(cmd, f"TransDecoder.Predict -> {output_dir}", timeout=timeout):
            return False

        # 4. ORF → 基因组坐标(出 gene/mRNA/CDS/exon)|ORF→genome coords
        td_gff3 = os.path.join(output_dir, f"{base}.transdecoder.gff3")
        cmd = ' '.join(build_conda_command(
            self.config.transdecoder_orf2genome_bin, [td_gff3, align_gff3, transcripts_fa])) + f" > {genome_gff3}"
        if not self.cmd_runner.run(cmd, f"TransDecoder: ORF→基因组坐标 -> {genome_gff3}"):
            return False

        return True


class BamInputParser:
    """BAM 输入解析与校验器|BAM input parser and validator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def validate_bam(self, bam_file: str) -> bool:
        """校验 BAM 完整性 + coordinate-sorted|Validate BAM is coordinate-sorted"""
        # samtools quickcheck
        qc_cmd_list = build_conda_command(self.config.samtools_bin, ['quickcheck', bam_file])
        qc = subprocess.run(' '.join(qc_cmd_list), shell=True, capture_output=True, text=True)
        if qc.returncode != 0:
            raise ValueError(f"BAM 损坏或不可读|BAM corrupt or unreadable: {bam_file}")

        # 检查 SO:coordinate|check sort order
        hdr_cmd_list = build_conda_command(self.config.samtools_bin, ['view', '-H', bam_file])
        hdr = subprocess.run(' '.join(hdr_cmd_list), shell=True, capture_output=True, text=True)
        if qc.returncode == 0 and 'SO:coordinate' not in hdr.stdout:
            raise ValueError(
                f"BAM 非坐标排序|BAM not coordinate-sorted (需先 samtools sort): {bam_file}")
        return True

    def ensure_index(self, bam_file: str) -> bool:
        """无 .bai 则建索引|Index if .bai missing"""
        if os.path.exists(f"{bam_file}.bai"):
            return True
        idx_args = ['index', '-@', str(self.config.threads), bam_file]
        cmd_list = build_conda_command(self.config.samtools_bin, idx_args)
        cmd = ' '.join(cmd_list)
        self.logger.info(f"命令|Command: {cmd}")
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"samtools index 失败|index failed: {bam_file} ({r.stderr})")
        return True

    def detect_read_type(self, bam_file: str) -> str:
        """采样读长中位数判断短/长读|Detect short/long by median read length"""
        view_cmd = ' '.join(build_conda_command(self.config.samtools_bin, ['view', bam_file]))
        # 单个 conda run + 系统 head/awk,无 conda|conda(§13.2.1)|single conda run + sys head/awk
        pipe = f"{view_cmd} | head -{READ_SAMPLE_COUNT} | awk '{{print length($10)}}'"
        self.logger.info(f"命令|Command: {pipe}")
        r = subprocess.run(pipe, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"读长检测失败|read-length detection failed: {bam_file} ({r.stderr})")
        lengths = [int(x) for x in r.stdout.split() if x.strip().isdigit()]
        if not lengths:
            raise ValueError(f"BAM 无 read 或为空|BAM has no reads (empty?): {bam_file}")
        median = statistics.median(lengths)
        rtype = 'long' if median >= READ_LENGTH_THRESHOLD else 'short'
        self.logger.info(f"读长检测|Read-length: {bam_file} 中位数|median={median} → {rtype}")
        return rtype

    def parse_bam_samples(self, bam_files: List[str], read_type: str) -> List[Dict]:
        """每个 BAM → {name, bam, read_type}|Build sample list from BAMs"""
        samples = []
        for bam in bam_files:
            self.validate_bam(bam)
            self.ensure_index(bam)
            name = os.path.splitext(os.path.basename(bam))[0]
            rtype = read_type if read_type != 'auto' else self.detect_read_type(bam)
            samples.append({'name': name, 'bam': bam, 'read_type': rtype})
        return samples


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
