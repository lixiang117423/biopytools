"""
BAM转FASTQ模块|BAM to FASTQ Conversion Module
从分类后的BAM文件中提取FASTQ文件|Extract FASTQ files from classified BAM files
"""

import os
from typing import List, Tuple
from .utils import CommandRunner


class BamToFastqExtractor:
    """BAM转FASTQ提取器|BAM to FASTQ Extractor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_fastq_from_bam(self, bam_file: str, output_prefix: str,
                               paired: bool = True,
                               compress: bool = True,
                               min_mapq: int = 0) -> bool:
        """从BAM文件提取FASTQ|Extract FASTQ from BAM file

        Args:
            bam_file: 输入BAM文件路径|Input BAM file path
            output_prefix: 输出文件前缀|Output file prefix
            paired: 是否为双端测序|Is paired-end sequencing (default: True)
            compress: 是否压缩输出|Compress output (default: True)
            min_mapq: 最小mapping quality|Minimum mapping quality (default: 0)

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            self.logger.info(f"从BAM提取FASTQ|Extracting FASTQ from BAM: {bam_file}")
            self.logger.info(f"输出前缀|Output prefix: {output_prefix}")
            self.logger.info(f"模式|Mode: {'双端|paired-end' if paired else '单端|single-end'}")

            # 检查BAM文件是否存在|Check if BAM file exists
            if not os.path.exists(bam_file):
                self.logger.error(f"BAM文件不存在|BAM file not found: {bam_file}")
                return False

            # 构建samtools命令|Build samtools command
            cmd = self._build_samtools_command(
                bam_file, output_prefix, paired, compress, min_mapq
            )

            # 执行命令|Execute command
            timeout = 7200  # 2小时超时|2 hours timeout
            success = self.cmd_runner.run(cmd, "提取FASTQ|Extracting FASTQ", timeout=timeout)

            if success:
                self.logger.info(f"FASTQ提取完成|FASTQ extraction completed: {output_prefix}")
                # 验证输出文件|Verify output files
                return self._verify_output(output_prefix, paired, compress)
            else:
                self.logger.error(f"FASTQ提取失败|FASTQ extraction failed")
                return False

        except Exception as e:
            self.logger.error(f"提取FASTQ时出错|Error extracting FASTQ: {e}")
            return False

    def _build_samtools_command(self, bam_file: str, output_prefix: str,
                                paired: bool, compress: bool,
                                min_mapq: int) -> str:
        """构建samtools命令|Build samtools command

        Args:
            bam_file: 输入BAM文件|Input BAM file
            output_prefix: 输出前缀|Output prefix
            paired: 是否双端|Is paired-end
            compress: 是否压缩|Is compressed
            min_mapq: 最小质量|Minimum quality

        Returns:
            str: samtools命令|samtools command
        """
        threads = self.config.threads

        # 确定输出文件扩展名|Determine output file extension
        ext = "fastq.gz" if compress else "fastq"

        # 构建命令|Build command
        if paired:
            # 双端模式|Paired-end mode
            r1_file = f"{output_prefix}_R1.{ext}"
            r2_file = f"{output_prefix}_R2.{ext}"
            singleton_file = f"{output_prefix}_singleton.{ext}"

            # 使用bash执行，因为命令使用了进程替换|Use bash for process substitution
            cmd = (
                f"bash -c \"samtools collate -@ {threads} -O -u "
                f"<(samtools view -@ {threads} -F 2816 -b {bam_file}) "
                f"/tmp/bam2fq_collate_{os.path.basename(output_prefix)} | "
                f"samtools fastq -@ {threads} "
            )

            if compress:
                cmd += f"-1 >(gzip -c > {r1_file}) "
                cmd += f"-2 >(gzip -c > {r2_file}) "
                cmd += f"-s >(gzip -c > {singleton_file}) "
            else:
                cmd += f"-1 {r1_file} "
                cmd += f"-2 {r2_file} "
                cmd += f"-s {singleton_file} "

            cmd += "-0 /dev/null -n -\""

        else:
            # 单端模式|Single-end mode
            se_file = f"{output_prefix}.{ext}"

            cmd = (
                f"samtools view -@ {threads} -F 2048 -b {bam_file} | "
                f"samtools fastq -@ {threads} -n -"
            )

            if compress:
                cmd += f" | gzip -c > {se_file}"
            else:
                cmd += f" > {se_file}"

        # 添加质量过滤|Add quality filter
        if min_mapq > 0:
            # 在view命令中添加-q参数，注意参数顺序：-@ threads -q min_mapq
            # Add -q parameter to view command, note param order: -@ threads -q min_mapq
            if paired:
                # 双端模式使用2816标志|Paired-end mode uses 2816 flag
                cmd = cmd.replace(
                    f"-@ {threads} -F 2816",
                    f"-@ {threads} -q {min_mapq} -F 2816"
                )
            else:
                # 单端模式使用2048标志|Single-end mode uses 2048 flag
                cmd = cmd.replace(
                    f"-@ {threads} -F 2048",
                    f"-@ {threads} -q {min_mapq} -F 2048"
                )

        return cmd

    def _verify_output(self, output_prefix: str, paired: bool,
                      compress: bool) -> bool:
        """验证输出文件|Verify output files

        Args:
            output_prefix: 输出前缀|Output prefix
            paired: 是否双端|Is paired-end
            compress: 是否压缩|Is compressed

        Returns:
            bool: 文件是否都存在|Whether all files exist
        """
        ext = "fastq.gz" if compress else "fastq"

        if paired:
            files_to_check = [
                f"{output_prefix}_R1.{ext}",
                f"{output_prefix}_R2.{ext}",
                f"{output_prefix}_singleton.{ext}"
            ]
        else:
            files_to_check = [f"{output_prefix}.{ext}"]

        all_exist = True
        for file_path in files_to_check:
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"输出文件已生成|Output file generated: {file_path} ({size} bytes)")
            else:
                self.logger.warning(f"输出文件未生成|Output file not generated: {file_path}")
                all_exist = False

        return all_exist

    def extract_sample_fastqs(self, sample_name: str,
                              classification_dir: str,
                              output_dir: str) -> bool:
        """提取单个样本的所有FASTQ文件|Extract all FASTQ files for a sample

        从分类后的BAM文件中提取FASTQ|Extract FASTQ from classified BAM files:
        - {sample}.{species1_name}.bam → {sample}.{species1_name}_R?.fastq.gz
        - {sample}.{species2_name}.bam → {sample}.{species2_name}_R?.fastq.gz
        - {sample}.ambiguous.bam → {sample}.ambiguous_R?.fastq.gz
        - {sample}.unassigned.bam → {sample}.unassigned_R?.fastq.gz

        Args:
            sample_name: 样本名|Sample name
            classification_dir: 分类目录|Classification directory
            output_dir: 输出目录|Output directory

        Returns:
            bool: 是否全部成功|Whether all succeeded
        """
        try:
            # 创建输出目录|Create output directory
            os.makedirs(output_dir, exist_ok=True)

            # 定义BAM文件类型|Define BAM file types
            bam_types = [
                (self.config.species1_name, f"{sample_name}.{self.config.species1_name}.bam"),
                (self.config.species2_name, f"{sample_name}.{self.config.species2_name}.bam"),
                ("ambiguous", f"{sample_name}.ambiguous.bam"),
                ("unassigned", f"{sample_name}.unassigned.bam")
            ]

            all_success = True

            for bam_type, bam_filename in bam_types:
                bam_path = os.path.join(classification_dir, bam_filename)

                # 检查BAM文件是否存在|Check if BAM file exists
                if not os.path.exists(bam_path):
                    self.logger.warning(f"BAM文件不存在，跳过|BAM file not found, skipping: {bam_path}")
                    continue

                # 定义输出前缀|Define output prefix
                output_prefix = os.path.join(output_dir, f"{sample_name}.{bam_type}")

                # 提取FASTQ|Extract FASTQ
                self.logger.info(f"提取{bam_type} FASTQ|Extracting {bam_type} FASTQ: {sample_name}")
                success = self.extract_fastq_from_bam(
                    bam_file=bam_path,
                    output_prefix=output_prefix,
                    paired=True,  # 双端测序|Paired-end
                    compress=True,
                    min_mapq=self.config.min_mapq
                )

                if not success:
                    self.logger.error(f"{bam_type} FASTQ提取失败|{bam_type} FASTQ extraction failed")
                    all_success = False

            return all_success

        except Exception as e:
            self.logger.error(f"提取样本FASTQ时出错|Error extracting sample FASTQ: {e}")
            return False

    def extract_all_samples(self, samples: List[str]) -> bool:
        """提取所有样本的FASTQ文件|Extract FASTQ files for all samples

        Args:
            samples: 样本列表|Sample list

        Returns:
            bool: 是否全部成功|Whether all succeeded
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始提取FASTQ文件|Starting FASTQ extraction")
            self.logger.info("=" * 60)

            # 创建输出目录|Create output directory
            fastq_dir = os.path.join(self.config.output_dir, "06.extracted_fastq")
            os.makedirs(fastq_dir, exist_ok=True)

            all_success = True

            for sample_name in samples:
                self.logger.info(f"处理样本|Processing sample: {sample_name}")

                classification_dir = os.path.join(
                    self.config.output_dir,
                    "02.classification",
                    sample_name
                )

                # 检查分类目录是否存在|Check if classification directory exists
                if not os.path.exists(classification_dir):
                    self.logger.warning(f"分类目录不存在，跳过|Classification directory not found, skipping: {classification_dir}")
                    continue

                # 提取样本FASTQ|Extract sample FASTQ
                success = self.extract_sample_fastqs(
                    sample_name=sample_name,
                    classification_dir=classification_dir,
                    output_dir=fastq_dir
                )

                if not success:
                    self.logger.error(f"样本{sample_name} FASTQ提取失败|Sample {sample_name} FASTQ extraction failed")
                    all_success = False

            if all_success:
                self.logger.info("=" * 60)
                self.logger.info("FASTQ提取完成|FASTQ extraction completed")
                self.logger.info(f"输出目录|Output directory: {fastq_dir}")
                self.logger.info("=" * 60)
            else:
                self.logger.warning("=" * 60)
                self.logger.warning("部分样本FASTQ提取失败|Some samples FASTQ extraction failed")
                self.logger.warning("=" * 60)

            return all_success

        except Exception as e:
            self.logger.error(f"提取FASTQ时出错|Error extracting FASTQ: {e}")
            return False
