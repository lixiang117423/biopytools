"""
VCF分离器模块 | VCF Separator Module
"""

from pathlib import Path

class VCFSeparator:
    """VCF文件分离器 - 分离SNP和INDEL | VCF File Separator - Separate SNPs and INDELs"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 定义输出文件名 | Define output file names
        # self.raw_snp_file = self.config.output_path / f"{self.config.base_name}.raw.snp.vcf.gz"
        # self.raw_indel_file = self.config.output_path / f"{self.config.base_name}.raw.indel.vcf.gz"
        # After
        self.raw_snp_file = Path(f"{self.config.base_name}.raw.snp.vcf.gz")
        self.raw_indel_file = Path(f"{self.config.base_name}.raw.indel.vcf.gz")
    
    def separate_variants(self) -> bool:
        """分离SNP和INDEL | Separate SNPs and INDELs"""
        self.logger.info("=" * 60)
        self.logger.info("🔬 步骤1: 分离SNP和INDEL | Step 1: Separating SNPs and INDELs")
        self.logger.info("=" * 60)
        
        # 提取SNP | Extract SNPs
        snp_cmd = (
            f"{self.config.bcftools_path} view "
            f"--types snps "
            f"--threads {self.config.threads} "
            f"-Oz -o {self.raw_snp_file} "
            f"{self.config.vcf_file}"
        )
        
        if not self.cmd_runner.run(snp_cmd, "提取SNP变异 | Extracting SNP variants"):
            return False
        
        # 提取INDEL | Extract INDELs
        indel_cmd = (
            f"{self.config.bcftools_path} view "
            f"--types indels "
            f"--threads {self.config.threads} "
            f"-Oz -o {self.raw_indel_file} "
            f"{self.config.vcf_file}"
        )
        
        if not self.cmd_runner.run(indel_cmd, "提取INDEL变异 | Extracting INDEL variants"):
            return False
        
        # 创建索引 | Create index
        self.logger.info("📇 创建VCF索引 | Creating VCF index")
        
        for vcf_file in [self.raw_snp_file, self.raw_indel_file]:
            index_cmd = f"{self.config.bcftools_path} index -f {vcf_file}"
            if not self.cmd_runner.run(index_cmd, f"索引文件 | Indexing {vcf_file.name}"):
                return False
        
        self.logger.info("✅ SNP和INDEL分离完成 | SNP and INDEL separation completed")
        return True
