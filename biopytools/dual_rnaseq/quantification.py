"""
互作转录组定量分析模块|Dual RNA-seq Quantification Module
"""

import os
import re
from .utils import CommandRunner, FileValidator


class DualQuantifier:
    """双物种定量分析器|Dual-species Quantifier"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)

    def run_stringtie(self, bam_file: str, gtf_file: str, output_gtf: str, species_name: str) -> bool:
        """运行StringTie定量|Run StringTie quantification"""
        threads = self.config.threads

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(output_gtf)
        os.makedirs(output_dir, exist_ok=True)

        # 检查GTF文件是否已存在|Check if GTF file already exists
        if self.file_validator.check_file_exists(output_gtf, f"GTF文件|GTF file"):
            return True

        cmd = f"stringtie -p {threads} -G {gtf_file} -o {output_gtf} -e {bam_file}"
        # 添加24小时超时限制|Add 24-hour timeout limit
        success = self.cmd_runner.run(cmd, f"StringTie定量{species_name}|StringTie quantification for {species_name}", timeout=86400)

        if success:
            self.logger.info(f"{species_name}定量完成|{species_name} quantification completed: {output_gtf}")

        return success

    def extract_gtf_values(self, gtf_file: str, sample_name: str, output_file: str) -> bool:
        """从GTF文件提取基因ID、FPKM和TPM值|Extract gene ID, FPKM and TPM values from GTF file"""
        try:
            # 确保输出目录存在|Ensure output directory exists
            os.makedirs(os.path.dirname(output_file), exist_ok=True)

            self.logger.info(f"提取FPKM和TPM值|Extracting FPKM and TPM values: {gtf_file}")

            # 正则表达式匹配键值对|Regular expression to match key-value pairs
            pattern = re.compile(r'(\w+)\s+"([^"]+)"')

            extracted_count = 0
            with open(gtf_file, "r") as infile, open(output_file, "w") as outfile:
                # 首先写入表头|First write the header
                outfile.write("gene_id\ttranscript_id\tcov\tFPKM\tTPM\tsample\n")

                for line in infile:
                    # 跳过注释行和不包含FPKM的行|Skip comment lines and lines without FPKM
                    if line.startswith("#") or "FPKM" not in line:
                        continue

                    # 分割行（保留第9列的完整属性）|Split line (keep complete attributes in column 9)
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue

                    # 解析属性列|Parse attributes column
                    attributes = {}
                    for match in pattern.finditer(cols[8]):
                        key, value = match.groups()
                        attributes[key] = value

                    # 获取所需的值（确保所有键都存在）|Get required values (ensure all keys exist)
                    required_keys = ["gene_id", "transcript_id", "cov", "FPKM", "TPM"]
                    if all(key in attributes for key in required_keys):
                        gene_id = attributes["gene_id"]
                        transcript_id = attributes["transcript_id"]
                        cov = attributes["cov"]
                        fpkm = attributes["FPKM"]
                        tpm = attributes["TPM"]

                        # 写入输出文件（制表符分隔）|Write to output file (tab-separated)
                        outfile.write(f"{gene_id}\t{transcript_id}\t{cov}\t{fpkm}\t{tpm}\t{sample_name}\n")
                        extracted_count += 1

            self.logger.info(f"GTF值提取完成|GTF value extraction completed: {output_file}")
            self.logger.info(f"提取了|Extracted {extracted_count}个转录本|transcripts")
            return True

        except Exception as e:
            self.logger.error(f"GTF值提取出错|Error extracting GTF values: {e}")
            return False

    def quantify_sample(self, sample_name: str) -> bool:
        """对单个样本进行双物种定量|Perform dual-species quantification for a single sample"""
        self.logger.info(f"开始样本定量|Starting sample quantification: {sample_name}")

        # 获取分类结果BAM文件|Get classified BAM files
        classification_dir = os.path.join(self.config.output_dir, "02.classification", sample_name)
        species1_bam = os.path.join(classification_dir, f"{sample_name}.{self.config.species1_name}.bam")
        species2_bam = os.path.join(classification_dir, f"{sample_name}.{self.config.species2_name}.bam")

        # 检查BAM文件是否存在|Check if BAM files exist
        if not os.path.exists(species1_bam):
            self.logger.error(f"找不到{self.config.species1_name} BAM文件|Cannot find {self.config.species1_name} BAM file: {species1_bam}")
            return False

        if not os.path.exists(species2_bam):
            self.logger.error(f"找不到{self.config.species2_name} BAM文件|Cannot find {self.config.species2_name} BAM file: {species2_bam}")
            return False

        # 创建定量输出目录|Create quantification output directory
        quant_dir = os.path.join(self.config.output_dir, "04.quantification")
        species1_quant_dir = os.path.join(quant_dir, self.config.species1_name)
        species2_quant_dir = os.path.join(quant_dir, self.config.species2_name)
        os.makedirs(species1_quant_dir, exist_ok=True)
        os.makedirs(species2_quant_dir, exist_ok=True)

        # 物种1定量|Species 1 quantification
        self.logger.info(f"定量{self.config.species1_name}|Quantifying {self.config.species1_name}")
        species1_gtf = os.path.join(species1_quant_dir, f"{sample_name}.gtf")
        species1_fpkm = os.path.join(species1_quant_dir, f"{sample_name}.fpkm.txt")

        if not self.run_stringtie(species1_bam, self.config.species1_gtf, species1_gtf, self.config.species1_name):
            self.logger.error(f"{self.config.species1_name} StringTie定量失败|{self.config.species1_name} StringTie quantification failed")
            return False

        if not self.extract_gtf_values(species1_gtf, sample_name, species1_fpkm):
            self.logger.error(f"{self.config.species1_name} FPKM提取失败|{self.config.species1_name} FPKM extraction failed")
            return False

        # 物种2定量|Species 2 quantification
        self.logger.info(f"定量{self.config.species2_name}|Quantifying {self.config.species2_name}")
        species2_gtf = os.path.join(species2_quant_dir, f"{sample_name}.gtf")
        species2_fpkm = os.path.join(species2_quant_dir, f"{sample_name}.fpkm.txt")

        if not self.run_stringtie(species2_bam, self.config.species2_gtf, species2_gtf, self.config.species2_name):
            self.logger.error(f"{self.config.species2_name} StringTie定量失败|{self.config.species2_name} StringTie quantification failed")
            return False

        if not self.extract_gtf_values(species2_gtf, sample_name, species2_fpkm):
            self.logger.error(f"{self.config.species2_name} FPKM提取失败|{self.config.species2_name} FPKM extraction failed")
            return False

        self.logger.info(f"样本定量完成|Sample quantification completed: {sample_name}")

        return True

    def quantify_all_samples(self, samples: list) -> bool:
        """对所有样本进行定量|Perform quantification for all samples"""
        self.logger.info(f"开始对所有样本进行定量|Starting quantification for all samples: {len(samples)}个样本|samples")

        success_count = 0

        for sample in samples:
            sample_name = sample["name"]
            if self.quantify_sample(sample_name):
                success_count += 1
            else:
                self.logger.warning(f"样本定量失败，跳过|Sample quantification failed, skipping: {sample_name}")

        self.logger.info(f"定量完成|Quantification completed: {success_count}/{len(samples)}个样本成功|samples successful")

        return success_count > 0
