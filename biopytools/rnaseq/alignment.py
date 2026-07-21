"""
RNA-seq|RNA-seq Alignment Module
"""

import os
from .utils import CommandRunner, FileValidator

class HISAT2Indexer:
    """HISAT2|HISAT2 Index Builder"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def extract_splice_sites(self, gtf_file: str, output_file: str) -> bool:
        """从GTF文件提取剪接位点|Extract splice sites from GTF file"""
        try:
            self.logger.info(f"提取剪接位点|Extracting splice sites: {gtf_file} -> {output_file}")

            # 检查输出文件|Check output file
            if os.path.exists(output_file):
                self.logger.info(f"剪接位点文件已存在，跳过|Splice sites file already exists, skipping: {output_file}")
                return True

            # HISAT2脚本|HISAT2 script
            cmd = f"extract_splice_sites.py {gtf_file} > {output_file}"
            success = self.cmd_runner.run(cmd, "提取剪接位点|Extract splice sites")

            if success:
                self.logger.info(f"剪接位点提取完成|Splice sites extraction completed: {output_file}")

            return success

        except Exception as e:
            self.logger.error(f"剪接位点提取出错|Error extracting splice sites: {e}")
            return False
    
    def extract_exons(self, gtf_file: str, output_file: str) -> bool:
        """从GTF文件提取外显子|Extract exons from GTF file"""
        try:
            self.logger.info(f"提取外显子|Extracting exons: {gtf_file} -> {output_file}")

            # 检查输出文件|Check output file
            if os.path.exists(output_file):
                self.logger.info(f"外显子文件已存在，跳过|Exons file already exists, skipping: {output_file}")
                return True

            # HISAT2脚本|HISAT2 script
            cmd = f"extract_exons.py {gtf_file} > {output_file}"
            success = self.cmd_runner.run(cmd, "提取外显子|Extract exons")

            if success:
                self.logger.info(f"外显子提取完成|Exons extraction completed: {output_file}")

            return success

        except Exception as e:
            self.logger.error(f"外显子提取出错|Error extracting exons: {e}")
            return False
    
    def build_hisat2_index(self) -> str:
        """构建HISAT2基因组索引|Build HISAT2 genome index"""
        genome_path = self.config.genome_file
        gtf_file = self.config.gtf_file
        threads = self.config.threads

        genome_dir = os.path.dirname(genome_path)
        genome_name = os.path.splitext(os.path.basename(genome_path))[0]
        index_prefix = os.path.join(genome_dir, f"{genome_name}.hisat2.index")

        # 检查索引是否已存在|Check if index already exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2索引已存在|HISAT2 index already exists: {index_prefix}")
            return index_prefix

        # 准备剪接位点和外显子文件路径|Prepare splice sites and exons file paths
        splice_sites_file = os.path.join(genome_dir, f"{genome_name}.ss")
        exons_file = os.path.join(genome_dir, f"{genome_name}.exon")

        # 提取剪接位点|Extract splice sites
        if not self.extract_splice_sites(gtf_file, splice_sites_file):
            self.logger.warning("剪接位点提取失败，将不使用剪接位点构建索引|Splice sites extraction failed, building index without splice sites")
            splice_sites_file = None

        # 提取外显子|Extract exons
        if not self.extract_exons(gtf_file, exons_file):
            self.logger.warning("外显子提取失败，将不使用外显子构建索引|Exons extraction failed, building index without exons")
            exons_file = None

        # 构建HISAT2索引命令|Build HISAT2 index command
        cmd_parts = [f"hisat2-build"]

        # 添加剪接位点参数|Add splice sites parameter
        if splice_sites_file and os.path.exists(splice_sites_file):
            cmd_parts.append(f"--ss {splice_sites_file}")
            self.logger.info(f"使用剪接位点文件|Using splice sites file: {splice_sites_file}")

        # 添加外显子参数|Add exons parameter
        if exons_file and os.path.exists(exons_file):
            cmd_parts.append(f"--exon {exons_file}")
            self.logger.info(f"使用外显子文件|Using exons file: {exons_file}")

        # 添加线程和基本参数|Add threads and basic parameters
        cmd_parts.extend([f"-p {threads}", genome_path, index_prefix])

        # 组合命令|Combine command
        cmd = " ".join(cmd_parts)

        self.logger.info(f"构建HISAT2索引|Building HISAT2 index with command: {cmd}")
        self.cmd_runner.run(cmd, "构建HISAT2索引|Building HISAT2 index")

        self.logger.info(f"HISAT2索引构建完成|HISAT2 index building completed: {index_prefix}")
        return index_prefix

class HISAT2Aligner:
    """HISAT2序列比对器|HISAT2 Aligner"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)

    def run_hisat2_mapping(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str) -> bool:
        """运行HISAT2比对和排序|Run HISAT2 alignment and sorting"""
        threads = self.config.threads
        timeout = getattr(self.config, 'sample_timeout', None)  # 获取超时配置|Get timeout config

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(output_bam)
        os.makedirs(output_dir, exist_ok=True)

        # 检查BAM文件是否已存在|Check if BAM file already exists
        if self.file_validator.check_file_exists(output_bam, "BAM文件|BAM file"):
            return True

        # 根据FASTQ文件大小决定管道策略|Decide pipe strategy based on FASTQ file size
        max_fastq_size = max(
            os.path.getsize(fastq1) if os.path.exists(fastq1) else 0,
            os.path.getsize(fastq2) if os.path.exists(fastq2) else 0
        )
        pipe_split_threshold = getattr(self.config, 'pipe_split_threshold', 10_000_000_000)

        if max_fastq_size > pipe_split_threshold:
            return self._run_hisat2_split_pipeline(index_prefix, fastq1, fastq2, output_bam, threads, timeout)
        else:
            return self._run_hisat2_pipe(index_prefix, fastq1, fastq2, output_bam, threads, timeout)

    def _run_hisat2_split_pipeline(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str, threads: int, timeout: int) -> bool:
        """逐步执行，无管道，适用于大数据量|Step-by-step execution without pipe for large data"""
        sam_file = output_bam.replace(".sorted.bam", ".sam")

        # 第一步：HISAT2比对输出SAM文件|Step 1: HISAT2 alignment to SAM file
        hisat2_cmd = (
            f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads} -S {sam_file}"
        )

        success = self.cmd_runner.run(
            hisat2_cmd,
            f"HISAT2比对|HISAT2 alignment -> {sam_file}",
            timeout=timeout
        )

        if not success:
            return False

        # 第二步：SAM转排序BAM|Step 2: Convert SAM to sorted BAM
        sort_cmd = (
            f"samtools sort -@ {threads} -m 2G -o {output_bam} {sam_file}"
        )

        success = self.cmd_runner.run(
            sort_cmd,
            f"BAM排序|BAM sorting -> {output_bam}",
            timeout=timeout
        )

        # 清理中间SAM文件|Clean up intermediate SAM file
        if os.path.exists(sam_file):
            os.remove(sam_file)

        if success:
            self.logger.info(f"比对完成|Alignment completed: {output_bam}")

        return success

    def _run_hisat2_pipe(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str, threads: int, timeout: int) -> bool:
        """管道方式执行，适用于小基因组|Pipe mode for small genomes"""
        cmd = (
            f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads}|"
            f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
        )

        success = self.cmd_runner.run(
            cmd,
            f"HISAT2比对和排序|HISAT2 alignment and sorting -> {output_bam}",
            timeout=timeout
        )

        if success:
            self.logger.info(f"比对完成|Alignment completed: {output_bam}")

        return success