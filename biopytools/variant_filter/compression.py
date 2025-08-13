"""
📦 压缩和索引模块 | Compression and Indexing Module
"""

from .utils import CommandRunner

class CompressionHandler:
    """📦 压缩处理器 | Compression Handler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def compress_and_index(self) -> bool:
        """📦 压缩并索引VCF文件 | Compress and index VCF file"""
        if not self.config.compress_output:
            self.logger.info("⏭️ 跳过压缩步骤 | Skipping compression step")
            return True
            
        self.logger.info("🚀 开始压缩和索引VCF文件 | Starting VCF compression and indexing")
        
        # 📦 压缩文件 | Compress file
        compress_cmd = f"{self.config.bgzip_path} {self.config.filtered_vcf}"
        
        success = self.cmd_runner.run(compress_cmd, "📦 bgzip压缩VCF文件 | bgzip compress VCF file")
        
        if not success:
            return False
        
        # 🏷️ 创建索引 | Create index
        if self.config.create_index:
            index_cmd = f"{self.config.tabix_path} -p vcf {self.config.final_vcf}"
            
            success = self.cmd_runner.run(index_cmd, "🏷️ tabix索引VCF文件 | tabix index VCF file")
            
            if success:
                self.logger.info(f"✅ 压缩和索引完成 | Compression and indexing completed: {self.config.final_vcf}")
        else:
            self.logger.info(f"✅ 压缩完成 | Compression completed: {self.config.final_vcf}")
        
        return success
