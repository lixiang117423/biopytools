"""
RNA-seq比对模块 | RNA-seq Alignment Module
"""

# import os
# from .utils import CommandRunner, FileValidator

# class HISAT2Indexer:
#     """HISAT2索引构建器 | HISAT2 Index Builder"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.file_validator = FileValidator(logger)
    
#     def build_hisat2_index(self) -> str:
#         """构建HISAT2基因组索引 | Build HISAT2 genome index"""
#         genome_path = self.config.genome_file
#         threads = self.config.threads
        
#         genome_dir = os.path.dirname(genome_path)
#         genome_name = os.path.splitext(os.path.basename(genome_path))[0]
#         index_prefix = os.path.join(genome_dir, f"{genome_name}.hisat2.index")

#         # 检查索引是否已存在 | Check if index already exists
#         if os.path.exists(f"{index_prefix}.1.ht2"):
#             self.logger.info(f"HISAT2索引已存在 | HISAT2 index already exists: {index_prefix}")
#             return index_prefix

#         cmd = f"hisat2-build -p {threads} {genome_path} {index_prefix}"
#         self.cmd_runner.run(cmd, "构建HISAT2索引 | Building HISAT2 index")
        
#         self.logger.info(f"HISAT2索引构建完成 | HISAT2 index building completed: {index_prefix}")
#         return index_prefix

# class HISAT2Aligner:
#     """HISAT2比对器 | HISAT2 Aligner"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.file_validator = FileValidator(logger)
    
#     def run_hisat2_mapping(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str) -> bool:
#         """运行HISAT2比对和排序 | Run HISAT2 alignment and sorting"""
#         threads = self.config.threads
        
#         # 确保输出目录存在 | Ensure output directory exists
#         output_dir = os.path.dirname(output_bam)
#         os.makedirs(output_dir, exist_ok=True)

#         # 检查BAM文件是否已存在 | Check if BAM file already exists
#         if self.file_validator.check_file_exists(output_bam, "BAM文件 | BAM file"):
#             return True

#         cmd = (
#             f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads} | "
#             f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
#         )

#         success = self.cmd_runner.run(cmd, f"HISAT2比对和排序 | HISAT2 alignment and sorting -> {output_bam}")
        
#         if success:
#             self.logger.info(f"比对完成 | Alignment completed: {output_bam}")
        
#         return success

"""
RNA-seq比对模块 | RNA-seq Alignment Module
"""

import os
from .utils import CommandRunner, FileValidator

class HISAT2Indexer:
    """HISAT2索引构建器 | HISAT2 Index Builder"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def extract_splice_sites(self, gtf_file: str, output_file: str) -> bool:
        """从GTF文件提取剪接位点信息 | Extract splice sites from GTF file"""
        try:
            self.logger.info(f"提取剪接位点信息 | Extracting splice sites: {gtf_file} -> {output_file}")
            
            # 检查output_file是否已存在
            if os.path.exists(output_file):
                self.logger.info(f"剪接位点文件已存在，跳过 | Splice sites file already exists, skipping: {output_file}")
                return True
            
            # 使用HISAT2自带的提取脚本
            cmd = f"extract_splice_sites.py {gtf_file} > {output_file}"
            success = self.cmd_runner.run(cmd, "提取剪接位点信息 | Extract splice sites")
            
            if success:
                self.logger.info(f"剪接位点信息提取完成 | Splice sites extraction completed: {output_file}")
            
            return success
            
        except Exception as e:
            self.logger.error(f"提取剪接位点信息时出错 | Error extracting splice sites: {e}")
            return False
    
    def extract_exons(self, gtf_file: str, output_file: str) -> bool:
        """从GTF文件提取外显子信息 | Extract exons from GTF file"""
        try:
            self.logger.info(f"提取外显子信息 | Extracting exons: {gtf_file} -> {output_file}")
            
            # 检查output_file是否已存在
            if os.path.exists(output_file):
                self.logger.info(f"外显子文件已存在，跳过 | Exons file already exists, skipping: {output_file}")
                return True
            
            # 使用HISAT2自带的提取脚本
            cmd = f"extract_exons.py {gtf_file} > {output_file}"
            success = self.cmd_runner.run(cmd, "提取外显子信息 | Extract exons")
            
            if success:
                self.logger.info(f"外显子信息提取完成 | Exons extraction completed: {output_file}")
            
            return success
            
        except Exception as e:
            self.logger.error(f"提取外显子信息时出错 | Error extracting exons: {e}")
            return False
    
    def build_hisat2_index(self) -> str:
        """构建HISAT2基因组索引 | Build HISAT2 genome index"""
        genome_path = self.config.genome_file
        gtf_file = self.config.gtf_file
        threads = self.config.threads
        
        genome_dir = os.path.dirname(genome_path)
        genome_name = os.path.splitext(os.path.basename(genome_path))[0]
        index_prefix = os.path.join(genome_dir, f"{genome_name}.hisat2.index")

        # 检查索引是否已存在 | Check if index already exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2索引已存在 | HISAT2 index already exists: {index_prefix}")
            return index_prefix

        # 准备剪接位点和外显子文件路径 | Prepare splice sites and exons file paths
        splice_sites_file = os.path.join(genome_dir, f"{genome_name}.ss")
        exons_file = os.path.join(genome_dir, f"{genome_name}.exon")
        
        # 提取剪接位点信息 | Extract splice sites
        if not self.extract_splice_sites(gtf_file, splice_sites_file):
            self.logger.warning("剪接位点提取失败，将使用基本模式构建索引 | Splice sites extraction failed, building index without splice sites")
            splice_sites_file = None
        
        # 提取外显子信息 | Extract exons
        if not self.extract_exons(gtf_file, exons_file):
            self.logger.warning("外显子信息提取失败，将使用基本模式构建索引 | Exons extraction failed, building index without exons")
            exons_file = None
        
        # 构建HISAT2索引命令 | Build HISAT2 index command
        cmd_parts = [f"hisat2-build"]
        
        # 添加剪接位点参数 | Add splice sites parameter
        if splice_sites_file and os.path.exists(splice_sites_file):
            cmd_parts.append(f"--ss {splice_sites_file}")
            self.logger.info(f"使用剪接位点文件 | Using splice sites file: {splice_sites_file}")
        
        # 添加外显子参数 | Add exons parameter
        if exons_file and os.path.exists(exons_file):
            cmd_parts.append(f"--exon {exons_file}")
            self.logger.info(f"使用外显子文件 | Using exons file: {exons_file}")
        
        # 添加线程数和基本参数 | Add threads and basic parameters
        cmd_parts.extend([f"-p {threads}", genome_path, index_prefix])
        
        # 组合命令 | Combine command
        cmd = " ".join(cmd_parts)
        
        self.logger.info(f"构建HISAT2索引 | Building HISAT2 index with command: {cmd}")
        self.cmd_runner.run(cmd, "构建HISAT2索引 | Building HISAT2 index")
        
        self.logger.info(f"HISAT2索引构建完成 | HISAT2 index building completed: {index_prefix}")
        return index_prefix

class HISAT2Aligner:
    """HISAT2比对器 | HISAT2 Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def run_hisat2_mapping(self, index_prefix: str, fastq1: str, fastq2: str, output_bam: str) -> bool:
        """运行HISAT2比对和排序 | Run HISAT2 alignment and sorting"""
        threads = self.config.threads
        
        # 确保输出目录存在 | Ensure output directory exists
        output_dir = os.path.dirname(output_bam)
        os.makedirs(output_dir, exist_ok=True)

        # 检查BAM文件是否已存在 | Check if BAM file already exists
        if self.file_validator.check_file_exists(output_bam, "BAM文件 | BAM file"):
            return True

        cmd = (
            f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads} | "
            f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
        )

        success = self.cmd_runner.run(cmd, f"HISAT2比对和排序 | HISAT2 alignment and sorting -> {output_bam}")
        
        if success:
            self.logger.info(f"比对完成 | Alignment completed: {output_bam}")
        
        return success