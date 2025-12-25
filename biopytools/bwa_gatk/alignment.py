"""
🧬 BWA比对模块 | BWA Alignment Module
"""

from pathlib import Path
from .utils import CommandRunner, check_file_exists

class BWAAligner:
    """BWA比对器 | BWA Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def align_sample(self, sample: dict, read_group: str) -> Path:
        """比对单个样本 | Align single sample"""
        sample_name = sample['name']
        
        self.logger.info("=" * 80)
        self.logger.info(f"🧬 比对样本 | Aligning sample: {sample_name}")
        self.logger.info("=" * 80)
        
        sam_file = self.config.temp_dir / f"{sample_name}.sam"
        
        # 检查断点续跑 | Check for resume
        if not self.config.force_restart and check_file_exists(sam_file, self.logger):
            self.logger.info("⏩ 跳过比对步骤（文件已存在）| Skipping alignment (file exists)")
            return sam_file
        
        # 构建BWA命令 | Build BWA command
        if sample['seq_type'] == "双端测序 | Paired-end":
            cmd = (f"{self.config.bwa_path} mem "
                   f"-t {self.config.threads} "
                   f"-R '{read_group}' "
                   f"{self.config.reference} "
                   f"{sample['r1']} {sample['r2']} "
                   f"> {sam_file}")
        else:
            cmd = (f"{self.config.bwa_path} mem "
                   f"-t {self.config.threads} "
                   f"-R '{read_group}' "
                   f"{self.config.reference} "
                   f"{sample['r1']} "
                   f"> {sam_file}")
        
        self.cmd_runner.run(cmd, f"BWA比对 | BWA alignment: {sample_name}")
        
        return sam_file
