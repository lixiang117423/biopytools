"""
互作转录组索引构建模块|Dual RNA-seq Index Building Module
"""

import os
from .utils import CommandRunner, FileValidator


class DualIndexBuilder:
    """双物种HISAT2索引构建器|Dual-species HISAT2 Index Builder"""

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

    def normalize_genome_fasta(self, genome_file: str, output_file: str) -> bool:
        """
        标准化基因组FASTA文件（将小写字母转为大写）|Normalize genome FASTA file (convert lowercase to uppercase)

        Args:
            genome_file: 输入基因组FASTA文件|Input genome FASTA file
            output_file: 输出标准化文件|Output normalized file

        Returns:
            是否成功|Whether successful
        """
        try:
            with open(genome_file, 'r') as f_in, open(output_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        # 保留header不变|Keep header unchanged
                        f_out.write(line)
                    else:
                        # 将序列转为大写|Convert sequence to uppercase
                        f_out.write(line.upper())
            self.logger.info(f"基因组文件标准化完成|Genome file normalized: {output_file}")
            return True
        except Exception as e:
            self.logger.error(f"基因组文件标准化失败|Failed to normalize genome file: {e}")
            return False

    def build_hisat2_index(self, genome_file: str, gtf_file: str, species_name: str) -> str:
        """构建HISAT2基因组索引|Build HISAT2 genome index"""
        threads = self.config.threads

        # 使用01.index子目录|Use 01.index subdirectory
        index_dir = os.path.join(self.config.output_dir, "01.index")
        os.makedirs(index_dir, exist_ok=True)

        genome_basename = os.path.splitext(os.path.basename(genome_file))[0]
        index_prefix = os.path.join(index_dir, f"{species_name}.hisat2.index")

        # 检查索引是否已存在|Check if index already exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2索引已存在|HISAT2 index already exists: {index_prefix}")
            return index_prefix

        # 标准化基因组文件（将小写转为大写）|Normalize genome file (convert lowercase to uppercase)
        normalized_genome = os.path.join(index_dir, f"{species_name}.normalized.fa")
        self.logger.info(f"标准化{species_name}基因组文件|Normalizing {species_name} genome file")
        if not self.normalize_genome_fasta(genome_file, normalized_genome):
            self.logger.error(f"{species_name}基因组文件标准化失败|{species_name} genome file normalization failed")
            return None

        # 使用标准化的基因组文件|Use normalized genome file
        genome_file = normalized_genome

        # 准备剪接位点和外显子文件路径|Prepare splice sites and exons file paths
        splice_sites_file = os.path.join(index_dir, f"{species_name}.ss")
        exons_file = os.path.join(index_dir, f"{species_name}.exon")

        # 提取剪接位点|Extract splice sites
        if not self.extract_splice_sites(gtf_file, splice_sites_file):
            self.logger.warning("剪接位点提取失败，将不使用剪接位点构建索引|Splice sites extraction failed, building index without splice sites")
            splice_sites_file = None

        # 提取外显子|Extract exons
        if not self.extract_exons(gtf_file, exons_file):
            self.logger.warning("外显子提取失败，将不使用外显子构建索引|Exons extraction failed, building index without exons")
            exons_file = None

        # 检查文件是否为空|Check if files are empty
        if splice_sites_file and os.path.exists(splice_sites_file):
            if os.path.getsize(splice_sites_file) == 0:
                self.logger.warning(f"剪接位点文件为空（{splice_sites_file}），将不使用剪接位点|Splice sites file is empty, building index without splice sites")
                splice_sites_file = None

        if exons_file and os.path.exists(exons_file):
            if os.path.getsize(exons_file) == 0:
                self.logger.warning(f"外显子文件为空（{exons_file}），将不使用外显子|Exons file is empty, building index without exons")
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
        cmd_parts.extend([f"-p {threads}", genome_file, index_prefix])

        # 组合命令|Combine command
        cmd = " ".join(cmd_parts)

        self.logger.info(f"构建{species_name} HISAT2索引|Building {species_name} HISAT2 index with command: {cmd}")
        success = self.cmd_runner.run(cmd, f"构建{species_name} HISAT2索引|Building {species_name} HISAT2 index")

        if success:
            self.logger.info(f"{species_name} HISAT2索引构建完成|{species_name} HISAT2 index building completed: {index_prefix}")
            return index_prefix
        else:
            self.logger.error(f"{species_name} HISAT2索引构建失败|{species_name} HISAT2 index building failed")
            return None

    def build_dual_indexes(self):
        """构建双物种HISAT2索引|Build dual-species HISAT2 indexes"""
        self.logger.info("开始构建双物种HISAT2索引|Starting to build dual-species HISAT2 indexes")

        # 构建物种1索引|Build species 1 index
        self.logger.info(f"构建{self.config.species1_name}索引|Building {self.config.species1_name} index")
        species1_index = self.build_hisat2_index(
            self.config.species1_genome,
            self.config.species1_gtf,
            self.config.species1_name
        )

        if not species1_index:
            self.logger.error(f"{self.config.species1_name}索引构建失败|{self.config.species1_name} index building failed")
            return None, None

        # 构建物种2索引|Build species 2 index
        self.logger.info(f"构建{self.config.species2_name}索引|Building {self.config.species2_name} index")
        species2_index = self.build_hisat2_index(
            self.config.species2_genome,
            self.config.species2_gtf,
            self.config.species2_name
        )

        if not species2_index:
            self.logger.error(f"{self.config.species2_name}索引构建失败|{self.config.species2_name} index building failed")
            return None, None

        self.logger.info("双物种HISAT2索引构建完成|Dual-species HISAT2 indexes building completed")

        return species1_index, species2_index
