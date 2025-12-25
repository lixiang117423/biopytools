"""
🔧 BAM预处理模块 | BAM Preprocessing Module
"""

from pathlib import Path
from .utils import CommandRunner, check_file_exists

class BAMPreprocessor:
    """BAM预处理器 | BAM Preprocessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def process_sample(self, sample_name: str, sam_file: Path) -> Path:
        """处理单个样本的BAM | Process single sample BAM"""
        self.logger.info("=" * 80)
        self.logger.info(f"🔧 BAM预处理 | BAM preprocessing: {sample_name}")
        self.logger.info("=" * 80)
        
        # 转换为BAM | Convert to BAM
        bam_file = self._sam_to_bam(sample_name, sam_file)
        
        # 排序 | Sort
        sorted_bam = self._sort_bam(sample_name, bam_file)
        
        # 标记重复 | Mark duplicates
        dedup_bam = self._mark_duplicates(sample_name, sorted_bam)
        
        # 建立索引 | Build index
        self._index_bam(dedup_bam)
        
        # 清理中间文件 | Clean up intermediate files
        self._cleanup(sam_file, bam_file)
        
        return dedup_bam
    
    def _sam_to_bam(self, sample_name: str, sam_file: Path) -> Path:
        """SAM转BAM | Convert SAM to BAM"""
        self.logger.info("📦 SAM转BAM | Converting SAM to BAM")
        
        bam_file = self.config.temp_dir / f"{sample_name}.bam"
        
        if not self.config.force_restart and check_file_exists(bam_file, self.logger):
            return bam_file
        
        cmd = (f"{self.config.samtools_path} view "
               f"-@ {self.config.threads} "
               f"-b {sam_file} "
               f"-o {bam_file}")
        
        self.cmd_runner.run(cmd, "SAM转BAM | SAM to BAM")
        return bam_file
    
    def _sort_bam(self, sample_name: str, bam_file: Path) -> Path:
        """排序BAM | Sort BAM"""
        self.logger.info("🔄 排序BAM | Sorting BAM")
        
        sorted_bam = self.config.temp_dir / f"{sample_name}.sorted.bam"
        
        if not self.config.force_restart and check_file_exists(sorted_bam, self.logger):
            return sorted_bam
        
        cmd = (f"{self.config.samtools_path} sort "
               f"-@ {self.config.threads} "
               f"-o {sorted_bam} "
               f"{bam_file}")
        
        self.cmd_runner.run(cmd, "排序BAM | Sort BAM")
        return sorted_bam
    
    def _mark_duplicates(self, sample_name: str, sorted_bam: Path) -> Path:
        """标记重复 | Mark duplicates"""
        self.logger.info("🏷️  标记重复序列 | Marking duplicates")
        
        dedup_bam = self.config.bam_dir / f"{sample_name}.dedup.bam"
        metrics_file = self.config.stats_dir / f"{sample_name}.dedup_metrics.txt"
        
        if not self.config.force_restart and check_file_exists(dedup_bam, self.logger):
            return dedup_bam
        
        cmd = (f"{self.config.gatk_path} MarkDuplicates "
               f"-I {sorted_bam} "
               f"-O {dedup_bam} "
               f"-M {metrics_file} "
               f"--TMP_DIR {self.config.temp_dir}")
        
        self.cmd_runner.run(cmd, "标记重复 | Mark duplicates")
        return dedup_bam
    
    def _index_bam(self, bam_file: Path):
        """建立BAM索引 | Index BAM"""
        self.logger.info("📑 建立BAM索引 | Indexing BAM")
        
        bai_file = Path(str(bam_file) + '.bai')
        
        if not self.config.force_restart and check_file_exists(bai_file, self.logger):
            return
        
        cmd = f"{self.config.samtools_path} index {bam_file}"
        self.cmd_runner.run(cmd, "建立BAM索引 | Index BAM")
    
    def _cleanup(self, sam_file: Path, bam_file: Path):
        """清理中间文件 | Clean up intermediate files"""
        if not self.config.dry_run:
            self.logger.info("🧹 清理中间文件 | Cleaning up intermediate files")
            
            # 删除SAM文件 | Delete SAM file
            if sam_file.exists():
                sam_file.unlink()
                self.logger.info(f"  🗑️  已删除 | Deleted: {sam_file.name}")
            
            # 删除未排序的BAM | Delete unsorted BAM
            if bam_file.exists():
                bam_file.unlink()
                self.logger.info(f"  🗑️  已删除 | Deleted: {bam_file.name}")
