"""
ANNOVAR数据处理模块 | ANNOVAR Data Processing Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, GFF3Validator

class GFF3Processor:
    """GFF3处理器 | GFF3 Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.gff3_validator = GFF3Validator(logger)
    
    def gff3_to_genepred(self):
        """GFF3文件转GenPred文件 | Convert GFF3 file to GenPred format"""
        gff3_file = self.config.gff3_file
        output_dir = self.config.output_dir
        build_ver = self.config.build_ver
        output_file = os.path.join(output_dir, f"{build_ver}_refGene.txt")
        
        self.logger.info(f"📂 GFF3文件 | GFF3 file: {gff3_file}")
        self.logger.info(f"📁 输出目录 | Output directory: {output_dir}")
        self.logger.info(f"📄 输出文件 | Output file: {output_file}")
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 清理和修复GFF3文件（除非用户指定跳过） | Clean and fix GFF3 file (unless user specifies to skip)
        if not self.config.skip_gff_cleaning:
            if not self.gff3_validator.clean_and_fix_gff3(gff3_file):
                self.logger.error("❌ GFF3文件清理失败 | GFF3 file cleaning failed")
                return False
        else:
            self.logger.info("⏭️ 跳过GFF3文件清理（用户指定） | Skipping GFF3 file cleaning (user specified)")
        
        # 检查GFF3文件格式 | Check GFF3 file format
        self.gff3_validator.check_gff3_header(gff3_file)
        
        # 修复CDS phase问题（除非用户指定跳过） | Fix CDS phase issues (unless user specifies to skip)
        if not self.config.skip_gff_fix:
            self.gff3_validator.fix_gff3_cds_phase(gff3_file)
        else:
            self.logger.info("⏭️ 跳过GFF3文件修复（用户指定） | Skipping GFF3 file fix (user specified)")
        
        command = f"gff3ToGenePred -warnAndContinue {gff3_file} {output_file}"
        
        success = self.cmd_runner.run(command, "🔄 GFF3转GenPred格式 | GFF3 to GenPred conversion")
        if success:
            self.config.genepred_file = output_file
            self.logger.info(f"✅ GenPred文件已生成 | GenPred file generated: {output_file}")
        return success

class SequenceExtractor:
    """序列提取器 | Sequence Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def extract_transcript_sequences(self):
        """提取转录本序列 | Extract transcript sequences"""
        if not hasattr(self.config, 'genepred_file'):
            self.logger.error("🚫 GenPred文件不存在，请先执行步骤1 | GenPred file does not exist, please run step 1 first")
            return False
        
        genome_file = self.config.genome_file
        genepred_file = self.config.genepred_file
        output_dir = self.config.output_dir
        build_ver = self.config.build_ver
        output_file = os.path.join(output_dir, f"{build_ver}_refGeneMrna.fa")
        annovar_path = self.config.annovar_path
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"🧬 基因组文件 | Genome file: {genome_file}")
        self.logger.info(f"📄 GenPred文件 | GenPred file: {genepred_file}")
        self.logger.info(f"🧬 输出序列文件 | Output sequence file: {output_file}")
        
        command = (f"perl {annovar_path}/retrieve_seq_from_fasta.pl "
                  f"--format refGene --seqfile {genome_file} "
                  f"{genepred_file} -outfile {output_file}")
        
        success = self.cmd_runner.run(command, "🧬 提取转录本序列 | Extract transcript sequences")
        if success:
            self.config.mrna_file = output_file
            self.logger.info(f"✅ 转录本序列文件已生成 | Transcript sequence file generated: {output_file}")
        return success

class VCFProcessor:
    """VCF处理器 | VCF Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def filter_and_convert_vcf(self):
        """过滤并转换VCF格式 | Filter and convert VCF format"""
        vcf_file = self.config.vcf_file
        qual_threshold = self.config.qual_threshold
        annovar_path = self.config.annovar_path
        output_dir = self.config.output_dir
        skip_filter = self.config.skip_vcf_filter
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 获取VCF文件的基础名称（不含扩展名） | Get VCF file base name (without extension)
        vcf_basename = os.path.splitext(os.path.basename(vcf_file))[0]
        if vcf_basename.endswith('.vcf'):
            vcf_basename = vcf_basename[:-4]  # 移除.vcf后缀 | Remove .vcf suffix
        
        if skip_filter:
            # 跳过过滤，直接使用原始VCF文件 | Skip filtering, use original VCF file directly
            self.logger.info("⏭️ 跳过VCF过滤步骤，直接使用输入的VCF文件 | Skipping VCF filtering, using input VCF file directly")
            filtered_vcf = vcf_file
        else:
            # 步骤3a: 过滤VCF文件 | Step 3a: Filter VCF file
            filtered_vcf = os.path.join(output_dir, f"{vcf_basename}.filtered.gz")
            filter_command = (f"bcftools filter -i 'QUAL>={qual_threshold}' "
                             f"{vcf_file} -O z -o {filtered_vcf}")
            
            self.logger.info(f"🔍 过滤后VCF文件 | Filtered VCF file: {filtered_vcf}")
            
            if not self.cmd_runner.run(filter_command, f"🔍 过滤VCF文件 | Filter VCF file (QUAL>={qual_threshold})"):
                return False
        
        # 步骤3b: 转换为ANNOVAR格式 | Step 3b: Convert to ANNOVAR format
        annovar_vcf = os.path.join(output_dir, f"{vcf_basename}.annovar.vcf")
        convert_command = (f"perl {annovar_path}/convert2annovar.pl "
                          f"-format vcf4 -allsample -withfreq {filtered_vcf} > {annovar_vcf}")
        
        self.logger.info(f"📥 输入VCF文件 | Input VCF file: {filtered_vcf}")
        self.logger.info(f"📤 ANNOVAR格式文件 | ANNOVAR format file: {annovar_vcf}")
        
        success = self.cmd_runner.run(convert_command, "🔄 转换VCF为ANNOVAR格式 | Convert VCF to ANNOVAR format")
        if success:
            self.config.annovar_vcf = annovar_vcf
            self.config.vcf_basename = vcf_basename
            self.logger.info(f"✅ ANNOVAR格式文件已生成 | ANNOVAR format file generated: {annovar_vcf}")
        return success