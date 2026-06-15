"""
VCF分离器模块|VCF Separator Module
"""

from pathlib import Path

class VCFSeparator:
    """VCF文件分离器 - 分离SNP和INDEL|VCF File Separator - Separate SNPs and INDELs"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 定义输出文件名|Define output file names
        # self.raw_snp_file = self.config.output_path / f"{self.config.base_name}.raw.snp.vcf.gz"
        # self.raw_indel_file = self.config.output_path / f"{self.config.base_name}.raw.indel.vcf.gz"
        # After
        self.raw_snp_file = Path(f"{self.config.base_name}.raw.snp.vcf.gz")
        self.raw_indel_file = Path(f"{self.config.base_name}.raw.indel.vcf.gz")
    
    def separate_variants(self) -> bool:
        """分离SNP和INDEL|Separate SNPs and INDELs"""
        self.logger.info("=" * 60)
        self.logger.info("步骤1: 分离SNP和INDEL|Step 1: Separating SNPs and INDELs")
        self.logger.info("=" * 60)

        # 根据variant_type决定提取哪些变异|Extract variants based on variant_type
        files_to_index = []

        # 提取SNP|Extract SNPs
        if self.config.variant_type in ['both', 'snp_only']:
            snp_cmd = (
                f"{self.config.bcftools_path} view "
                f"--types snps "
                f"--threads {self.config.threads} "
                f"-Oz -o {self.raw_snp_file} "
                f"{self.config.vcf_file}"
            )

            if not self.cmd_runner.run(snp_cmd, "提取SNP变异|Extracting SNP variants"):
                return False
            files_to_index.append(self.raw_snp_file)
        else:
            self.logger.info("跳过SNP提取|Skipping SNP extraction")

        # 提取INDEL|Extract INDELs
        if self.config.variant_type in ['both', 'indel_only']:
            indel_cmd = (
                f"{self.config.bcftools_path} view "
                f"--types indels "
                f"--threads {self.config.threads} "
                f"-Oz -o {self.raw_indel_file} "
                f"{self.config.vcf_file}"
            )

            if not self.cmd_runner.run(indel_cmd, "提取INDEL变异|Extracting INDEL variants"):
                return False
            files_to_index.append(self.raw_indel_file)
        else:
            self.logger.info("跳过INDEL提取|Skipping INDEL extraction")

        # 创建索引|Create index
        if files_to_index:
            self.logger.info("创建VCF索引|Creating VCF index")
            for vcf_file in files_to_index:
                index_cmd = f"{self.config.bcftools_path} index -f {vcf_file}"
                if not self.cmd_runner.run(index_cmd, f"索引文件|Indexing {vcf_file.name}"):
                    return False

        self.logger.info("变异分离完成|Variant separation completed")
        return True
