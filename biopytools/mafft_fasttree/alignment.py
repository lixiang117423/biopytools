"""
多序列比对模块|Multiple Sequence Alignment Module
"""

import subprocess
from pathlib import Path

from ..common.conda_runner import build_conda_command


class MAFFTAligner:
    """MAFFT多序列比对器|MAFFT Multiple Sequence Aligner"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_alignment(self, input_file: Path, output_file: Path) -> bool:
        """运行MAFFT比对|Run MAFFT alignment"""
        self.logger.info("=" * 60)
        self.logger.info("开始多序列比对|Starting multiple sequence alignment")
        self.logger.info("=" * 60)
        
        # 构建MAFFT命令|Build MAFFT command
        mafft_args = self.config.mafft_params.split() + [
            "--thread", str(self.config.threads),
            str(input_file.resolve()),
        ]
        cmd = build_conda_command(self.config.mafft_path, mafft_args)
        self.logger.info(f"命令|Command: {' '.join(cmd)} > {output_file.resolve()}")

        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            success = (result.returncode == 0)
            if not success and result.stderr:
                self.logger.error(f" MAFFT比对失败|MAFFT alignment failed: {result.stderr}")
        except Exception as e:
            self.logger.error(f" MAFFT比对失败|MAFFT alignment failed: {e}")
            return False
        
        if success:
            self.logger.info(f" MAFFT比对完成|MAFFT alignment completed: {output_file}")
        else:
            self.logger.error(f" MAFFT比对失败|MAFFT alignment failed")
        
        return success
