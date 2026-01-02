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
            self.logger.error("🚫 ANNOVAR格式VCF文件不存在，请先执行步骤3 | ANNOVAR format VCF file does not exist, please run step 3 first")
            return False
        
        annovar_vcf = self.config.annovar_vcf
        annovar_path = self.config.annovar_path
        build_ver = self.config.build_ver
        output_dir = self.config.output_dir
        vcf_basename = self.config.vcf_basename
        output_prefix = os.path.join(output_dir, vcf_basename)
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"📝 注释输出前缀 | Annotation output prefix: {output_prefix}")
        
        # 检查必需的refGene文件是否存在 | Check if required refGene file exists
        refgene_file = os.path.join(output_dir, f"{build_ver}_refGene.txt")
        if not os.path.exists(refgene_file):
            self.logger.error(f"🚫 必需的基因注释文件不存在 | Required gene annotation file does not exist: {refgene_file}")
            self.logger.error("请先执行步骤1生成GenePred文件 | Please run step 1 first to generate GenePred file")
            return False
        
        # 使用输出目录作为数据库路径，因为refGene文件就在那里 | Use output directory as database path since refGene file is there
        command = (f"perl {annovar_path}/annotate_variation.pl "
                  f"--geneanno -dbtype refGene --buildver {build_ver} "
                  f"{annovar_vcf} {output_dir} -out {output_prefix}")
        
        success = self.cmd_runner.run(command, "🧬 ANNOVAR变异注释 | ANNOVAR variant annotation")
        
        if success:
            # 记录输出文件 | Record output files
            output_files = [
                f"{output_prefix}.exonic_variant_function",
                f"{output_prefix}.variant_function",
                f"{output_prefix}.log"
            ]
            
            existing_files = [f for f in output_files if os.path.exists(f)]
            self.config.output_files = existing_files
            
            self.logger.info("✅ 注释完成，输出文件 | Annotation completed, output files:")
            for file in existing_files:
                self.logger.info(f"  📄 {file}")
        
        return success