"""
ANNOVAR数据处理模块|ANNOVAR Data Processing Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, GFF3Validator, build_conda_command
from ..common.paths import get_tool_path

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
		"""GFF3文件转GenPred文件|Convert GFF3 file to GenePred format"""
		gff3_file = self.config.gff3_file
		output_dir = self.config.output_dir
		build_ver = self.config.build_ver
		output_file = os.path.join(output_dir, f"{build_ver}_refGene.txt")

		self.logger.info(f"GFF3文件|GFF3 file: {gff3_file}")
		self.logger.info(f"输出目录|Output directory: {output_dir}")
		self.logger.info(f"输出文件|Output file: {output_file}")

		# 确保输出目录存在|Ensure output directory exists
		os.makedirs(output_dir, exist_ok=True)

		# 清理和修复GFF3文件（除非用户指定跳过）|Clean and fix GFF3 file (unless user specifies to skip)
		if not self.config.skip_gff_cleaning:
			if not self.gff3_validator.clean_and_fix_gff3(gff3_file):
				self.logger.error("GFF3文件清理失败|GFF3 file cleaning failed")
				return False
		else:
			self.logger.info("跳过GFF3文件清理（用户指定）|Skipping GFF3 file cleaning (user specified)")

		# 检查GFF3文件格式|Check GFF3 file format
		self.gff3_validator.check_gff3_header(gff3_file)

		# 修复CDS phase问题（除非用户指定跳过）|Fix CDS phase issues (unless user specifies to skip)
		if not self.config.skip_gff_fix:
			self.gff3_validator.fix_gff3_cds_phase(gff3_file)
		else:
			self.logger.info("跳过GFF3文件修复（用户指定）|Skipping GFF3 file fix (user specified)")

		command = ' '.join(build_conda_command('gff3ToGenePred', [
			'-warnAndContinue', '-maxParseErrors=-1', '-maxConvertErrors=-1',
			gff3_file, output_file
		]))

		success = self.cmd_runner.run(command, "GFF3转GenPred格式|GFF3 to GenPred conversion")
		if success:
			self.config.genepred_file = output_file
			self.logger.info(f"GenPred文件已生成|GenPred file generated: {output_file}")
			self._validate_genepred_coverage(gff3_file, output_file)
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

		# 确保输出目录存在|Ensure output directory exists
		os.makedirs(output_dir, exist_ok=True)

		self.logger.info(f"基因组文件|Genome file: {genome_file}")
		self.logger.info(f"GenPred文件|GenPred file: {genepred_file}")
		self.logger.info(f"输出序列文件|Output sequence file: {output_file}")
		self.logger.info(f"输出蛋白文件|Output protein file: {pep_file}")
		self.logger.info(f"输出CDS文件|Output CDS file: {cds_file}")

		# 提取mRNA序列|Extract mRNA sequences
		command = (f"perl {annovar_path}/retrieve_seq_from_fasta.pl "
				  f"--format refGene --seqfile {genome_file} "
				  f"{genepred_file} -outfile {output_file}")

		success = self.cmd_runner.run(command, "提取转录本序列|Extract transcript sequences")
		if success:
			self.config.mrna_file = output_file
			self.logger.info(f"转录本序列文件已生成|Transcript sequence file generated: {output_file}")

		if not success:
			return False

		# 提取蛋白序列|Extract protein sequences
		gffread_path = get_tool_path('gffread', '~/miniforge3/envs/RNA_Seq/bin/gffread', 'GFFREAD_PATH')
		pep_command = f"{gffread_path} -g {genome_file} -y {pep_file} {gff3_file}"

		success = self.cmd_runner.run(pep_command, "提取蛋白序列|Extract protein sequences")
		if success:
			self.config.pep_file = pep_file
			self.logger.info(f"蛋白序列文件已生成|Protein sequence file generated: {pep_file}")

		if not success:
			return False

		# 提取CDS序列|Extract CDS sequences
		cds_command = f"{gffread_path} -g {genome_file} -x {cds_file} {gff3_file}"

		success = self.cmd_runner.run(cds_command, "提取CDS序列|Extract CDS sequences")
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

		# 确保输出目录存在|Ensure output directory exists
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
