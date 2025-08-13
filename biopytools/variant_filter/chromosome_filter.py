"""
🧩 染色体过滤模块 | Chromosome Filtering Module
"""

from .utils import CommandRunner

class ChromosomeFilter:
    """🧩 染色体过滤器 | Chromosome Filter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def filter_chromosomes(self) -> bool:
        """🧩 过滤指定模式的染色体 | Filter chromosomes with specified pattern"""
        if not self.config.chromosome_filter:
            self.logger.info("⏭️ 跳过染色体过滤 | Skipping chromosome filtering")
            return True
            
        self.logger.info(f"🧩 过滤染色体，模式: {self.config.chromosome_pattern} | Filtering chromosomes with pattern: {self.config.chromosome_pattern}")
        
        # 🔍 使用zcat和grep过滤指定染色体 | Use zcat and grep to filter specific chromosomes
        cmd = (
            f"zcat {self.config.final_vcf} | "
            f"grep -E \"^#|{self.config.chromosome_pattern}\" | "
            f"bgzip > {self.config.chr_filtered_vcf}"
        )
        
        success = self.cmd_runner.run(cmd, f"🧩 染色体过滤 | Chromosome filtering")
        
        if success:
            self.logger.info(f"✅ 染色体过滤完成 | Chromosome filtering completed: {self.config.chr_filtered_vcf}")
            
            # 🏷️ 为染色体过滤后的文件创建索引 | Create index for chromosome-filtered file
            if self.config.create_index:
                index_cmd = f"{self.config.tabix_path} -p vcf {self.config.chr_filtered_vcf}"
                self.cmd_runner.run(index_cmd, "🏷️ 为染色体过滤文件创建索引 | Create index for chromosome-filtered file")
        
        return success
