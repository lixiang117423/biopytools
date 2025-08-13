"""
📉 VCFtools变异过滤模块 | VCFtools Variant Filtering Module
"""

from .utils import CommandRunner

class VariantFilter:
    """🔍 VCFtools变异过滤器 | VCFtools Variant Filter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def filter_variants(self) -> bool:
        """🔍 使用VCFtools过滤变异 | Filter variants using VCFtools"""
        self.logger.info("🚀 开始使用VCFtools过滤变异 | Starting variant filtering with VCFtools")
        
        # 🔧 构建基本命令 | Build basic command
        cmd_parts = [
            self.config.vcftools_path,
            f"--gzvcf {self.config.selected_vcf}",
            f"--maf {self.config.maf}",
            f"--max-missing {self.config.max_missing}",
            f"--hwe {self.config.hwe}",
            f"--min-meanDP {self.config.min_meanDP}",
            f"--max-meanDP {self.config.max_meanDP}",
            f"--minQ {self.config.minQ}",
            f"--minGQ {self.config.minGQ}"
        ]
        
        # 🎯 根据变异类型添加特定过滤 | Add variant-type specific filtering
        if self.config.variant_type == "SNP":
            cmd_parts.append("--remove-indels")
            self.logger.info("🧬 添加SNP特定过滤：移除INDEL | Adding SNP-specific filter: remove INDELs")
        elif self.config.variant_type == "INDEL":
            cmd_parts.append("--keep-only-indels")
            self.logger.info("🧬 添加INDEL特定过滤：只保留INDEL | Adding INDEL-specific filter: keep only INDELs")
        # BOTH类型不添加额外过滤 | No additional filtering for BOTH type
        
        # ⚙️ 添加可选参数 | Add optional parameters
        if self.config.max_alleles != 2:
            cmd_parts.append(f"--max-alleles {self.config.max_alleles}")
            
        if self.config.min_alleles != 2:
            cmd_parts.append(f"--min-alleles {self.config.min_alleles}")
            
        if self.config.max_maf is not None:
            cmd_parts.append(f"--max-maf {self.config.max_maf}")
            
        if self.config.hardy_weinberg_p is not None:
            cmd_parts.append(f"--hardy-weinberg-p {self.config.hardy_weinberg_p}")
            
        if self.config.min_DP is not None:
            cmd_parts.append(f"--minDP {self.config.min_DP}")
            
        if self.config.max_DP is not None:
            cmd_parts.append(f"--maxDP {self.config.max_DP}")
            
        if self.config.indv_min_depth is not None:
            cmd_parts.append(f"--minDP {self.config.indv_min_depth}")
            
        if self.config.indv_max_depth is not None:
            cmd_parts.append(f"--maxDP {self.config.indv_max_depth}")
            
        if self.config.thin is not None:
            cmd_parts.append(f"--thin {self.config.thin}")
        
        # 💾 添加输出参数 | Add output parameters
        # VCFtools会自动添加.recode.vcf后缀 | VCFtools automatically adds .recode.vcf suffix
        output_prefix = f"{self.config.output_prefix}.filtered.{self.config.variant_type.lower()}"
        cmd_parts.extend([
            "--recode --recode-INFO-all",
            f"--out {self.config.output_dir / output_prefix}"
        ])
        
        # 🔗 组合命令 | Combine command
        cmd = " \\\n    ".join(cmd_parts)
        
        success = self.cmd_runner.run(cmd, f"📉 VCFtools变异过滤 | VCFtools variant filtering")
        
        if success:
            self.logger.info(f"✅ 变异过滤完成 | Variant filtering completed: {self.config.filtered_vcf}")
        
        return success
