"""
ANNOVAR数据处理模块|ANNOVAR Data Processing Module
"""

import os
import shutil
from pathlib import Path
from .utils import CommandRunner, GFF3Validator, build_conda_command


class GFF3Processor:
    """GFF3处理器|GFF3 Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.gff3_validator = GFF3Validator(logger)

    def _get_gff_chromosomes(self, gff3_file: str) -> set:
        """提取GFF3中所有包含gene特征的染色体|Extract chromosomes with gene features from GFF3"""
        chromosomes = set()
        with open(gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) >= 3 and fields[2] == 'gene':
                    chromosomes.add(fields[0])
        return chromosomes

    def _get_genepred_chromosomes(self, genepred_file: str) -> set:
        """提取GenePred文件中的染色体|Extract chromosomes from GenePred file"""
        chromosomes = set()
        with open(genepred_file, 'r') as f:
            for line in f:
                fields = line.split('\t')
                if len(fields) >= 2:
                    chromosomes.add(fields[1])
        return chromosomes

    def _validate_genepred_coverage(self, gff3_file: str, genepred_file: str) -> bool:
        """校验GenePred是否覆盖GFF3中的所有染色体|Validate GenePred covers all GFF3 chromosomes"""
        gff_chroms = self._get_gff_chromosomes(gff3_file)
        genepred_chroms = self._get_genepred_chromosomes(genepred_file)

        if not gff_chroms:
            self.logger.warning("GFF3中未找到gene特征|No gene features found in GFF3")
            return True

        genepred_count = sum(1 for _ in open(genepred_file))
        self.logger.info(
            f"GFF基因染色体数|GFF gene chromosomes: {len(gff_chroms)}, "
            f"GenePred染色体数|GenePred chromosomes: {len(genepred_chroms)}, "
            f"GenePred条目数|GenePred entries: {genepred_count}"
        )

        missing = gff_chroms - genepred_chroms
        if missing:
            self.logger.error(
                f"GenePred缺失{len(missing)}条染色体的基因模型|"
                f"GenePred missing gene models for {len(missing)} chromosomes: "
                f"{sorted(missing)}"
            )
            self.logger.error(
                "gff3ToGenePred转换可能不完整，建议检查GFF3文件格式或重试|"
                "gff3ToGenePred conversion may be incomplete, check GFF3 format or retry"
            )
            return False

        self.logger.info(
            f"染色体覆盖校验通过|Chromosome coverage check passed: "
            f"{len(genepred_chroms)}条染色体, {genepred_count}个转录本|"
            f"{len(genepred_chroms)} chromosomes, {genepred_count} transcripts"
        )
        return True

    def gff3_to_genepred(self):
        """
        GFF3文件转GenPred文件|Convert GFF3 file to GenePred format

        关键: 所有清洗/修复都在 output_dir 内的工作副本上进行,绝不修改用户输入的GFF3|
        Key: all cleaning/fixing happens on a working copy inside output_dir; the user's
        input GFF3 is never modified.
        """
        gff3_file = self.config.gff3_file
        output_dir = self.config.output_dir
        build_ver = self.config.build_ver
        output_file = os.path.join(output_dir, f"{build_ver}_refGene.txt")

        # 工作副本(始终在output_dir,输入文件保持不动)|working copy in output_dir, input untouched
        working_gff = os.path.join(output_dir, f"{build_ver}.cleaned.gff3")

        self.logger.info(f"GFF3文件(输入)|GFF3 file (input): {gff3_file}")
        self.logger.info(f"工作副本|Working copy: {working_gff}")
        self.logger.info(f"输出文件|Output file: {output_file}")

        os.makedirs(output_dir, exist_ok=True)

        if not self.config.skip_gff_cleaning:
            # 清洗+修复坐标: 读输入,写出工作副本(不碰输入)|clean+fix coords: read input, write working copy
            if not self.gff3_validator.clean_and_fix_gff3(gff3_file, working_gff):
                self.logger.error("GFF3文件清理失败|GFF3 file cleaning failed")
                return False
        else:
            self.logger.info("跳过GFF3文件清洗(用户指定),复制输入为工作副本|"
                             "Skipping GFF3 cleaning (user specified), copying input as working copy")
            shutil.copy2(gff3_file, working_gff)

        # header与CDS phase修复都在工作副本上进行(输出目录内,安全)|fix header/phase on the working copy (in output_dir, safe)
        self.gff3_validator.check_gff3_header(working_gff)

        if not self.config.skip_gff_fix:
            self.gff3_validator.fix_gff3_cds_phase(working_gff)
        else:
            self.logger.info("跳过GFF3文件修复(用户指定)|Skipping GFF3 file fix (user specified)")

        # 重定向后续步骤使用工作副本(gff3ToGenePred + step2的gffread)|redirect downstream to working copy
        self.config.gff3_file = working_gff

        command = ' '.join(build_conda_command('gff3ToGenePred', [
            '-warnAndContinue', '-maxParseErrors=-1', '-maxConvertErrors=-1',
            working_gff, output_file
        ]))

        success = self.cmd_runner.run(command, "GFF3转GenPred格式|GFF3 to GenPred conversion")
        if success:
            self.config.genepred_file = output_file
            self.logger.info(f"GenPred文件已生成|GenPred file generated: {output_file}")
            # 覆盖率校验基于工作副本|coverage check against the working copy
            self._validate_genepred_coverage(working_gff, output_file)
        return success


class SequenceExtractor:
    """序列提取器|Sequence Extractor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_transcript_sequences(self):
        """提取转录本序列、蛋白序列和CDS序列|Extract transcript, protein and CDS sequences"""
        genome_file = self.config.genome_file
        genepred_file = self.config.genepred_file
        gff3_file = self.config.gff3_file
        output_dir = self.config.output_dir
        build_ver = self.config.build_ver
        output_file = os.path.join(output_dir, f"{build_ver}_refGeneMrna.fa")
        pep_file = os.path.join(output_dir, f"{build_ver}_refGenePep.fa")
        cds_file = os.path.join(output_dir, f"{build_ver}_refGeneCds.fa")
        annovar_path = self.config.annovar_path

        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"基因组文件|Genome file: {genome_file}")
        self.logger.info(f"GenPred文件|Genepred file: {genepred_file}")
        self.logger.info(f"输出序列文件|Output sequence file: {output_file}")
        self.logger.info(f"输出蛋白文件|Output protein file: {pep_file}")
        self.logger.info(f"输出CDS文件|Output CDS file: {cds_file}")

        # 提取mRNA序列(ANNOVAR脚本,用perl直接调用)|Extract mRNA sequences via ANNOVAR perl script
        command = (f"perl {annovar_path}/retrieve_seq_from_fasta.pl "
                   f"--format refGene --seqfile {genome_file} "
                   f"{genepred_file} -outfile {output_file}")

        success = self.cmd_runner.run(command, "提取转录本序列|Extract transcript sequences")
        if success:
            self.config.mrna_file = output_file
            self.logger.info(f"转录本序列文件已生成|Transcript sequence file generated: {output_file}")

        if not success:
            return False

        # 提取蛋白序列(gffread,经build_conda_command解析env)|Extract proteins via gffread (conda-wrapped)
        gffread_cmd = ' '.join(build_conda_command(
            self.config.gffread_path,
            ['-g', genome_file, '-y', pep_file, gff3_file]
        ))
        success = self.cmd_runner.run(gffread_cmd, "提取蛋白序列|Extract protein sequences")
        if success:
            self.config.pep_file = pep_file
            self.logger.info(f"蛋白序列文件已生成|Protein sequence file generated: {pep_file}")

        if not success:
            return False

        # 提取CDS序列|Extract CDS sequences
        cds_cmd = ' '.join(build_conda_command(
            self.config.gffread_path,
            ['-g', genome_file, '-x', cds_file, gff3_file]
        ))
        success = self.cmd_runner.run(cds_cmd, "提取CDS序列|Extract CDS sequences")
        if success:
            self.config.cds_file = cds_file
            self.logger.info(f"CDS序列文件已生成|CDS sequence file generated: {cds_file}")

        return success


class VCFProcessor:
    """VCF处理器|VCF Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def filter_and_convert_vcf(self):
        """过滤并转换VCF格式|Filter and convert VCF format"""
        vcf_file = self.config.vcf_file
        qual_threshold = self.config.qual_threshold
        annovar_path = self.config.annovar_path
        output_dir = self.config.output_dir
        skip_filter = self.config.skip_vcf_filter

        os.makedirs(output_dir, exist_ok=True)

        # 获取VCF文件的基础名称（不含扩展名）|Get VCF file base name (without extension)
        vcf_basename = os.path.splitext(os.path.basename(vcf_file))[0]
        if vcf_basename.endswith('.vcf'):
            vcf_basename = vcf_basename[:-4]  # 移除.vcf后缀|Remove .vcf suffix

        if skip_filter:
            # 跳过过滤，直接使用原始VCF文件|Skip filtering, use original VCF file directly
            self.logger.info("跳过VCF过滤步骤，直接使用输入的VCF文件|Skipping VCF filtering, using input VCF file directly")
            filtered_vcf = vcf_file
        else:
            # 步骤3a: 过滤VCF文件|Step 3a: Filter VCF file
            filtered_vcf = os.path.join(output_dir, f"{vcf_basename}.filtered.gz")
            filter_command = ' '.join(build_conda_command(
                'bcftools',
                ['filter', '-i', f'QUAL>={qual_threshold}', vcf_file, '-O', 'z', '-o', filtered_vcf]
            ))

            self.logger.info(f"过滤后VCF文件|Filtered VCF file: {filtered_vcf}")

            if not self.cmd_runner.run(filter_command, f"过滤VCF文件|Filter VCF file (QUAL>={qual_threshold})"):
                return False

        # 步骤3b: 转换为ANNOVAR格式|Step 3b: Convert to ANNOVAR format
        annovar_vcf = os.path.join(output_dir, f"{vcf_basename}.annovar.vcf")
        convert_command = (f"perl {annovar_path}/convert2annovar.pl "
                           f"-format vcf4 -allsample -withfreq {filtered_vcf} > {annovar_vcf}")

        self.logger.info(f"输入VCF文件|Input VCF file: {filtered_vcf}")
        self.logger.info(f"ANNOVAR格式文件|ANNOVAR format file: {annovar_vcf}")

        success = self.cmd_runner.run(convert_command, "转换VCF为ANNOVAR格式|Convert VCF to ANNOVAR format")
        if success:
            self.config.annovar_vcf = annovar_vcf
            self.config.vcf_basename = vcf_basename
            self.logger.info(f"ANNOVAR格式文件已生成|ANNOVAR format file generated: {annovar_vcf}")
        return success
