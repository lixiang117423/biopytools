"""
Minimap2比对模块|Minimap2 Alignment Module
"""

import os
from .utils import CommandRunner
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command
from ..common.conda_runner import build_conda_command

class Minimap2Aligner:
    """Minimap2比对器|Minimap2 Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_alignment(self):
        """运行Minimap2比对|Run Minimap2 alignment"""
        target_genome = self.config.target_genome
        query_genome = self.config.query_genome
        output_paf = self.config.paf_file
        preset = self.config.preset
        threads = self.config.threads
        minimap2_path = self.config.minimap2_path
        
        self.logger.info(f"目标基因组|Target genome: {target_genome}")
        self.logger.info(f"查询基因组|Query genome: {query_genome}")
        self.logger.info(f"输出PAF文件|Output PAF file: {output_paf}")
        self.logger.info(f"预设参数|Preset: {preset}")
        self.logger.info(f"线程数|Threads: {threads}")
        
        # 检查PAF文件是否已存在|Check if PAF file already exists
        if os.path.exists(output_paf):
            self.logger.info(f"PAF文件已存在，跳过比对|PAF file already exists, skipping alignment: {output_paf}")
            return True
        
        # 构建minimap2命令(conda环境自动包装, 输出重定向到文件)
        # |Build minimap2 command (auto conda wrap, redirect stdout to file)
        args = [
            "-x", preset,
            "-t", str(threads),
            target_genome,
            query_genome,
        ]
        command = build_conda_command(minimap2_path, args)

        success, _, _ = self.cmd_runner.run(
            command,
            "Minimap2全基因组比对|Minimap2 whole genome alignment",
            output_file=output_paf,
        )
        
        if success:
            self.logger.info(f"比对完成，PAF文件已生成|Alignment completed, PAF file generated: {output_paf}")
            
            # 检查PAF文件是否为空|Check if PAF file is empty
            if os.path.getsize(output_paf) == 0:
                self.logger.warning("警告：PAF文件为空，可能没有找到比对结果|Warning: PAF file is empty, no alignments found")
                return False
        
        return success
