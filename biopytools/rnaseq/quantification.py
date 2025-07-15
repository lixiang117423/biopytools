"""
RNA-seq定量分析模块 | RNA-seq Quantification Module
"""

import os
import re
from .utils import CommandRunner, FileValidator

class StringTieQuantifier:
    """StringTie定量器 | StringTie Quantifier"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def run_stringtie(self, bam_file: str, output_gtf: str) -> bool:
        """运行StringTie定量 | Run StringTie quantification"""
        gtf_file = self.config.gtf_file
        threads = self.config.threads
        
        # 确保输出目录存在 | Ensure output directory exists
        output_dir = os.path.dirname(output_gtf)
        os.makedirs(output_dir, exist_ok=True)

        # 检查GTF文件是否已存在 | Check if GTF file already exists
        if self.file_validator.check_file_exists(output_gtf, "GTF文件 | GTF file"):
            return True

        cmd = f"stringtie -p {threads} -G {gtf_file} -o {output_gtf} -e {bam_file}"
        success = self.cmd_runner.run(cmd, f"StringTie定量 | StringTie quantification -> {output_gtf}")
        
        if success:
            self.logger.info(f"定量完成 | Quantification completed: {output_gtf}")
        
        return success

class GTFValueExtractor:
    """GTF值提取器 | GTF Value Extractor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def extract_gtf_values(self, gtf_file: str, sample_name: str, output_file: str) -> bool:
        """从GTF文件中提取基因ID、FPKM和TPM值 | Extract gene ID, FPKM and TPM values from GTF file"""
        try:
            # 确保输出目录存在 | Ensure output directory exists
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            # 用于匹配键值对的正则表达式 | Regular expression to match key-value pairs
            pattern = re.compile(r'(\w+)\s+"([^"]+)"')

            with open(gtf_file, "r") as infile, open(output_file, "w") as outfile:
                # 1. 首先写入头部 | First write the header
                outfile.write("gene_id\ttranscript_id\tcov\tFPKM\tTPM\tsample\n")

                for line in infile:
                    # 跳过注释行和没有FPKM的行 | Skip comment lines and lines without FPKM
                    if line.startswith("#") or "FPKM" not in line:
                        continue

                    # 分割行（保持第9列的完整属性） | Split line (keep complete attributes in column 9)
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue

                    # 解析属性列 | Parse attributes column
                    attributes = {}
                    for match in pattern.finditer(cols[8]):
                        key, value = match.groups()
                        attributes[key] = value

                    # 获取所需值（确保所有键都存在） | Get required values (ensure all keys exist)
                    required_keys = ["gene_id", "transcript_id", "cov", "FPKM", "TPM"]
                    if all(key in attributes for key in required_keys):
                        gene_id = attributes["gene_id"]
                        transcript_id = attributes["transcript_id"]
                        cov = attributes["cov"]
                        fpkm = attributes["FPKM"]
                        tpm = attributes["TPM"]

                        # 写入输出文件（制表符分隔） | Write to output file (tab-separated)
                        outfile.write(f"{gene_id}\t{transcript_id}\t{cov}\t{fpkm}\t{tpm}\t{sample_name}\n")

            self.logger.info(f"GTF值提取完成 | GTF value extraction completed: {output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"提取GTF值时出错 | Error extracting GTF values: {e}")
            return False
