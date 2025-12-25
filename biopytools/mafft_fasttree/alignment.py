"""
🧩 多序列比对模块 | Multiple Sequence Alignment Module
"""

from pathlib import Path

class MAFFTAligner:
    """MAFFT多序列比对器 | MAFFT Multiple Sequence Aligner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_alignment(self, input_file: Path, output_file: Path) -> bool:
        """运行MAFFT比对 | Run MAFFT alignment"""
        self.logger.info("=" * 60)
        self.logger.info("🧩 开始多序列比对 | Starting multiple sequence alignment")
        self.logger.info("=" * 60)
        
        # 构建MAFFT命令
        # mafft_cmd = f"{self.config.mafft_path} {self.config.mafft_params} "
        # mafft_cmd += f"--thread {self.config.threads} "
        # mafft_cmd += f"{input_file} > {output_file}"
        mafft_cmd = f"{self.config.mafft_path} {self.config.mafft_params} "
        mafft_cmd += f"--thread {self.config.threads} "
        mafft_cmd += f"{input_file.resolve()} > {output_file.resolve()}"
        
        success = self.cmd_runner.run(
            mafft_cmd,
            description="MAFFT多序列比对 | MAFFT multiple sequence alignment"
        )
        
        if success:
            self.logger.info(f"✅ MAFFT比对完成 | MAFFT alignment completed: {output_file}")
        else:
            self.logger.error(f"❌ MAFFT比对失败 | MAFFT alignment failed")
        
        return success
