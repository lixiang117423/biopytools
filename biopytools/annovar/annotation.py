"""
ANNOVAR注释核心模块 | ANNOVAR Annotation Core Module
"""

import os
from .utils import CommandRunner

class VariantAnnotator:
    """变异注释器 | Variant Annotator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def annotate_variants(self):
        """注释变异 | Annotate variants"""
        if not hasattr(self.config, 'annovar_vcf'):
            self.logger.error("ANNOVAR格式VCF文件不存在，请先执行步骤3 | ANNOVAR format VCF file does not exist, please run step 3 first")
            return False
        
        annovar_vcf = self.config.annovar_vcf
        annovar_path = self.config.annovar_path
        database_path = self.config.database_path
        build_ver = self.config.build_ver
        output_dir = self.config.output_dir
        vcf_basename = self.config.vcf_basename
        output_prefix = os.path.join(output_dir, vcf_basename)
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"注释输出前缀 | Annotation output prefix: {output_prefix}")
        
        command = (f"perl {annovar_path}/annotate_variation.pl "
                  f"--geneanno -dbtype refGene --buildver {build_ver} "
                  f"{annovar_vcf} {database_path} -out {output_prefix}")
        
        success = self.cmd_runner.run(command, "ANNOVAR变异注释 | ANNOVAR variant annotation")
        
        if success:
            # 记录输出文件 | Record output files
            output_files = [
                f"{output_prefix}.exonic_variant_function",
                f"{output_prefix}.variant_function",
                f"{output_prefix}.log"
            ]
            
            existing_files = [f for f in output_files if os.path.exists(f)]
            self.config.output_files = existing_files
            
            self.logger.info("注释完成，输出文件 | Annotation completed, output files:")
            for file in existing_files:
                self.logger.info(f"  - {file}")
        
        return success
