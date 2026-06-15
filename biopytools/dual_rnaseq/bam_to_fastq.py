"""
BAM转FASTQ模块|BAM to FASTQ Conversion Module
从分类后的BAM文件中提取FASTQ文件|Extract FASTQ files from classified BAM files
"""

import os
from typing import List
from .utils import CommandRunner


class BamToFastqExtractor:
    """BAM转FASTQ提取器|BAM to FASTQ Extractor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_fastq_from_bam(self, bam_file: str, output_prefix: str,
                               paired: bool = True,
                               compress: bool = True) -> bool:
        """从BAM文件提取FASTQ|Extract FASTQ from BAM file"""
        try:
            if not os.path.exists(bam_file):
                self.logger.error(f"BAM文件不存在|BAM file not found: {bam_file}")
                return False

            ext = "fastq.gz" if compress else "fastq"
            threads = self.config.threads

            if paired:
                r1_file = f"{output_prefix}_1.clean.{ext}"
                r2_file = f"{output_prefix}_2.clean.{ext}"

                if compress:
                    r1_tmp = f"{output_prefix}_1.clean.fq"
                    r2_tmp = f"{output_prefix}_2.clean.fq"

                    self.logger.info(f"从BAM提取FASTQ|Extracting FASTQ from BAM: {bam_file}")
                    cmd = f"samtools fastq -@ {threads} -1 {r1_tmp} -2 {r2_tmp} {bam_file}"
                    self.logger.info(f"命令|Command: {cmd}")
                    success = self.cmd_runner.run(cmd, "提取FASTQ|Extracting FASTQ", timeout=7200)
                    if not success:
                        return False

                    self.logger.info("压缩FASTQ文件|Compressing FASTQ files")
                    for fq in [r1_tmp, r2_tmp]:
                        cmd = f"gzip -f {fq}"
                        self.logger.info(f"命令|Command: {cmd}")
                        if not self.cmd_runner.run(cmd, f"压缩|Compressing {os.path.basename(fq)}", timeout=7200):
                            return False

                    return self._verify_output(output_prefix, paired, compress)
                else:
                    cmd = f"samtools fastq -@ {threads} -1 {r1_file} -2 {r2_file} {bam_file}"
                    self.logger.info(f"从BAM提取FASTQ|Extracting FASTQ from BAM: {bam_file}")
                    self.logger.info(f"命令|Command: {cmd}")
                    success = self.cmd_runner.run(cmd, "提取FASTQ|Extracting FASTQ", timeout=7200)
                    if not success:
                        return False
                    return self._verify_output(output_prefix, paired, compress)

            else:
                se_file = f"{output_prefix}.{ext}"
                if compress:
                    se_tmp = f"{output_prefix}.fastq"
                    cmd = f"samtools fastq -@ {threads} -o {se_tmp} {bam_file}"
                    self.logger.info(f"从BAM提取FASTQ|Extracting FASTQ from BAM: {bam_file}")
                    self.logger.info(f"命令|Command: {cmd}")
                    success = self.cmd_runner.run(cmd, "提取FASTQ|Extracting FASTQ", timeout=7200)
                    if not success:
                        return False

                    cmd = f"gzip -f {se_tmp}"
                    self.logger.info("压缩FASTQ文件|Compressing FASTQ files")
                    self.logger.info(f"命令|Command: {cmd}")
                    if not self.cmd_runner.run(cmd, f"压缩|Compressing {os.path.basename(se_tmp)}", timeout=7200):
                        return False

                    return self._verify_output(output_prefix, paired, compress)
                else:
                    cmd = f"samtools fastq -@ {threads} -o {se_file} {bam_file}"
                    self.logger.info(f"从BAM提取FASTQ|Extracting FASTQ from BAM: {bam_file}")
                    self.logger.info(f"命令|Command: {cmd}")
                    success = self.cmd_runner.run(cmd, "提取FASTQ|Extracting FASTQ", timeout=7200)
                    if not success:
                        return False
                    return self._verify_output(output_prefix, paired, compress)

        except Exception as e:
            self.logger.error(f"提取FASTQ时出错|Error extracting FASTQ: {e}")
            return False

    def _verify_output(self, output_prefix: str, paired: bool,
                      compress: bool) -> bool:
        """验证输出文件|Verify output files"""
        ext = "fastq.gz" if compress else "fastq"

        if paired:
            files_to_check = [
                f"{output_prefix}_1.clean.{ext}",
                f"{output_prefix}_2.clean.{ext}"
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
        """提取单个样本的所有FASTQ文件|Extract all FASTQ files for a sample"""
        try:
            os.makedirs(output_dir, exist_ok=True)

            bam_types = [
                (self.config.species1_name, f"{sample_name}.{self.config.species1_name}.bam"),
                (self.config.species2_name, f"{sample_name}.{self.config.species2_name}.bam"),
                ("ambiguous", f"{sample_name}.ambiguous.bam"),
                ("unassigned", f"{sample_name}.unassigned.bam")
            ]

            all_success = True

            for bam_type, bam_filename in bam_types:
                bam_path = os.path.join(classification_dir, bam_filename)

                if not os.path.exists(bam_path):
                    self.logger.warning(f"BAM文件不存在，跳过|BAM file not found, skipping: {bam_path}")
                    continue

                output_prefix = os.path.join(output_dir, f"{sample_name}.{bam_type}")
                output_r1 = f"{output_prefix}_1.clean.fq.gz"
                output_r2 = f"{output_prefix}_2.clean.fq.gz"

                if os.path.exists(output_r1) and os.path.exists(output_r2):
                    self.logger.info(f"跳过已完成|Skipping completed: {bam_type} ({sample_name})")
                    continue

                self.logger.info(f"提取{bam_type} FASTQ|Extracting {bam_type} FASTQ: {sample_name}")
                success = self.extract_fastq_from_bam(
                    bam_file=bam_path,
                    output_prefix=output_prefix,
                    paired=True,
                    compress=True
                )

                if not success:
                    self.logger.error(f"{bam_type} FASTQ提取失败|{bam_type} FASTQ extraction failed")
                    all_success = False

            return all_success

        except Exception as e:
            self.logger.error(f"提取样本FASTQ时出错|Error extracting sample FASTQ: {e}")
            return False

    def extract_all_samples(self, samples: List[str]) -> bool:
        """提取所有样本的FASTQ文件|Extract FASTQ files for all samples"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始提取FASTQ文件|Starting FASTQ extraction")
            self.logger.info("=" * 60)

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

                if not os.path.exists(classification_dir):
                    self.logger.warning(f"分类目录不存在，跳过|Classification directory not found, skipping: {classification_dir}")
                    continue

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
