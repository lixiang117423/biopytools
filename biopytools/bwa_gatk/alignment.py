"""
 BWA比对模块|BWA Alignment Module
"""

from pathlib import Path
from .utils import CommandRunner, check_file_exists
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command

class BWAAligner:
    """BWA比对器|BWA Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def align_sample(self, sample: dict, read_group: str) -> Path:
        """比对单个样本|Align single sample"""
        sample_name = sample['name']
        
        self.logger.info("=" * 80)
        self.logger.info(f" 比对样本|Aligning sample: {sample_name}")
        self.logger.info("=" * 80)
        
        sam_file = self.config.temp_dir / f"{sample_name}.sam"
        
        # 检查断点续跑|Check for resume
        if not self.config.force_restart and check_file_exists(sam_file, self.logger):
            self.logger.info("跳过比对步骤（文件已存在）| Skipping alignment (file exists)")
            return sam_file
        
        # 构建BWA命令(conda环境自动包装, SAM输出重定向到文件)
        # |Build BWA command (auto conda wrap, redirect SAM to file)
        args = ["mem", "-t", str(self.config.threads),
                "-R", read_group, self.config.reference]
        if sample['seq_type'] == "双端测序|Paired-end":
            args += [sample['r1'], sample['r2']]
        else:
            args += [sample['r1']]
        cmd = build_conda_command(self.config.bwa_path, args)
        
        # 用公共执行器重定向输出|Redirect output via common runner
        success, _, _ = self.cmd_runner.run(
            cmd, f"BWA比对|BWA alignment: {sample_name}", output_file=str(sam_file)
        )
        
        return sam_file
