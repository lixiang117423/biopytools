"""
🧬 变异筛选配置模块 | Variant Filtering Configuration Module
"""

import os
from pathlib import Path
from typing import Optional, Union

class VariantFilterConfig:
    """🧬 变异筛选配置类 | Variant Filtering Configuration Class"""
    
    def __init__(
        self,
        input_vcf: str,
        output_dir: str = "./variant_filter_output",
        output_prefix: str = "filtered_variants",
        variant_type: str = "SNP",  # SNP, INDEL, or BOTH
        # 🔧 工具路径 | Tool paths
        gatk_path: str = "gatk",
        vcftools_path: str = "vcftools", 
        bcftools_path: str = "bcftools",
        bgzip_path: str = "bgzip",
        tabix_path: str = "tabix",
        # 🗂️ 变异选择参数 | Variant selection parameters
        java_options: str = "-Xmx128g",
        use_gatk_select: bool = False,  # 默认不使用GATK选择变异
        exclude_filtered: bool = False,
        exclude_non_variants: bool = False,
        remove_unused_alternates: bool = False,
        restrict_alleles_to: str = "ALL",
        select_random_fraction: Optional[float] = None,
        max_indel_size: Optional[int] = None,
        min_indel_size: Optional[int] = None,
        max_filtered_genotypes: Optional[int] = None,
        min_filtered_genotypes: Optional[int] = None,
        max_fraction_filtered_genotypes: Optional[float] = None,
        min_fraction_filtered_genotypes: Optional[float] = None,
        max_nocall_number: Optional[int] = None,
        max_nocall_fraction: Optional[float] = None,
        # 📉 vcftools 过滤参数 | vcftools filtering parameters
        maf: Optional[float] = None,  # 将根据变异类型自动设置
        max_missing: Optional[float] = None,  # 将根据变异类型自动设置
        hwe: float = 0.001,
        min_meanDP: int = 5,
        max_meanDP: int = 50,
        minQ: int = 30,
        minGQ: int = 20,
        max_alleles: int = 2,
        min_alleles: int = 2,
        max_maf: Optional[float] = None,
        hardy_weinberg_p: Optional[float] = None,
        min_DP: Optional[int] = None,
        max_DP: Optional[int] = None,
        indv_min_depth: Optional[int] = None,
        indv_max_depth: Optional[int] = None,
        thin: Optional[int] = None,
        # 🧩 染色体过滤参数 | Chromosome filtering parameters
        chromosome_filter: bool = False,
        chromosome_pattern: str = "^OV",  # 过滤以OV开头的染色体
        # 📦 压缩和索引参数 | Compression and indexing parameters
        compress_output: bool = True,
        create_index: bool = True,
        # ⚙️ 其他参数 | Other parameters
        threads: int = 1,
        memory: str = "128g",
        temp_dir: Optional[str] = None,
        keep_intermediate: bool = False,
        verbose: bool = True
    ):
        # 📁 文件路径设置 | File path settings
        self.input_vcf = os.path.normpath(os.path.abspath(input_vcf))
        self.output_dir = Path(output_dir).resolve()
        self.output_prefix = output_prefix
        self.variant_type = variant_type.upper()
        
        # 🎯 根据变异类型设置默认参数 | Set default parameters based on variant type
        if maf is None:
            self.maf = 0.05 if self.variant_type == "SNP" else 0.05
        else:
            self.maf = maf
            
        if max_missing is None:
            self.max_missing = 0.9 if self.variant_type == "SNP" else 1.0
        else:
            self.max_missing = max_missing
        
        # 🔧 工具路径 | Tool paths
        self.gatk_path = gatk_path
        self.vcftools_path = vcftools_path
        self.bcftools_path = bcftools_path
        self.bgzip_path = bgzip_path
        self.tabix_path = tabix_path
        
        # 🗂️ 变异选择参数 | Variant selection parameters
        self.java_options = java_options
        self.use_gatk_select = use_gatk_select
        self.exclude_filtered = exclude_filtered
        self.exclude_non_variants = exclude_non_variants
        self.remove_unused_alternates = remove_unused_alternates
        self.restrict_alleles_to = restrict_alleles_to
        self.select_random_fraction = select_random_fraction
        self.max_indel_size = max_indel_size
        self.min_indel_size = min_indel_size
        self.max_filtered_genotypes = max_filtered_genotypes
        self.min_filtered_genotypes = min_filtered_genotypes
        self.max_fraction_filtered_genotypes = max_fraction_filtered_genotypes
        self.min_fraction_filtered_genotypes = min_fraction_filtered_genotypes
        self.max_nocall_number = max_nocall_number
        self.max_nocall_fraction = max_nocall_fraction
        
        # 📉 vcftools参数 | vcftools parameters
        self.hwe = hwe
        self.min_meanDP = min_meanDP
        self.max_meanDP = max_meanDP
        self.minQ = minQ
        self.minGQ = minGQ
        self.max_alleles = max_alleles
        self.min_alleles = min_alleles
        self.max_maf = max_maf
        self.hardy_weinberg_p = hardy_weinberg_p
        self.min_DP = min_DP
        self.max_DP = max_DP
        self.indv_min_depth = indv_min_depth
        self.indv_max_depth = indv_max_depth
        self.thin = thin
        
        # 🧩 染色体过滤参数 | Chromosome filtering parameters
        self.chromosome_filter = chromosome_filter
        self.chromosome_pattern = chromosome_pattern
        
        # 📦 压缩工具参数 | Compression tool parameters
        self.compress_output = compress_output
        self.create_index = create_index
        
        # ⚙️ 其他参数 | Other parameters
        self.threads = threads
        self.memory = memory
        self.temp_dir = temp_dir
        self.keep_intermediate = keep_intermediate
        self.verbose = verbose
        
        # 💾 生成输出文件路径 | Generate output file paths
        self.output_dir.mkdir(parents=True, exist_ok=True)
        variant_suffix = self.variant_type.lower() if self.variant_type != "BOTH" else "variants"
        self.selected_vcf = self.output_dir / f"{self.output_prefix}.{variant_suffix}.vcf.gz"
        # VCFtools会在输出文件名后自动添加.recode.vcf | VCFtools automatically adds .recode.vcf to output filename
        self.filtered_vcf = self.output_dir / f"{self.output_prefix}.filtered.{variant_suffix}.recode.vcf"
        self.final_vcf = self.output_dir / f"{self.output_prefix}.filtered.{variant_suffix}.recode.vcf.gz"
        if self.chromosome_filter:
            self.chr_filtered_vcf = self.output_dir / f"{self.output_prefix}.filtered.{variant_suffix}.recode.chr.vcf.gz"
        
    def validate(self):
        """✅ 验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 🔍 检查输入VCF文件 | Check input VCF file
        if not os.path.exists(self.input_vcf):
            errors.append(f"❌ 输入VCF文件不存在 | Input VCF file does not exist: {self.input_vcf}")
        
        # 🔍 检查变异类型 | Check variant type
        if self.variant_type not in ["SNP", "INDEL", "BOTH"]:
            errors.append(f"❌ 不支持的变异类型 | Unsupported variant type: {self.variant_type}")
        
        # 📊 检查参数范围 | Check parameter ranges
        if not 0 <= self.maf <= 0.5:
            errors.append(f"❌ MAF值必须在0-0.5之间 | MAF must be between 0-0.5: {self.maf}")
            
        if not 0 <= self.max_missing <= 1:
            errors.append(f"❌ max_missing必须在0-1之间 | max_missing must be between 0-1: {self.max_missing}")
            
        if not 0 <= self.hwe <= 1:
            errors.append(f"❌ HWE p值必须在0-1之间 | HWE p-value must be between 0-1: {self.hwe}")
            
        if self.min_meanDP < 0:
            errors.append(f"❌ 最小平均深度必须为非负数 | Minimum mean depth must be non-negative: {self.min_meanDP}")
            
        if self.max_meanDP <= self.min_meanDP:
            errors.append(f"❌ 最大平均深度必须大于最小平均深度 | Max mean depth must be greater than min mean depth")
            
        if self.minQ < 0:
            errors.append(f"❌ 最小质量值必须为非负数 | Minimum quality must be non-negative: {self.minQ}")
            
        if self.minGQ < 0:
            errors.append(f"❌ 最小基因型质量值必须为非负数 | Minimum genotype quality must be non-negative: {self.minGQ}")
            
        if self.threads < 1:
            errors.append(f"❌ 线程数必须大于0 | Number of threads must be greater than 0: {self.threads}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
