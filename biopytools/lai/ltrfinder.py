"""
LTR_FINDER运行器|LTR_FINDER Runner
"""

import subprocess
import os
from pathlib import Path
from .utils import build_conda_command


class LTRFinderRunner:
    """LTR_FINDER运行器类|LTR_FINDER Runner Class"""

    def __init__(self, config, logger, finder_config):
        """
        初始化LTR_FINDER运行器|Initialize LTR_FINDER runner

        Args:
            config: LAI配置对象|LAI configuration object
            logger: 日志器|Logger
            finder_config: LTR_FINDER配置对象|LTR_FINDER configuration object
        """
        self.config = config
        self.logger = logger
        self.finder_config = finder_config

    def run(self) -> Path:
        """
        运行LTR_FINDER识别LTR候选|Run LTR_FINDER to identify LTR candidates

        Returns:
            Path: 输出文件路径|Output file path
        """
        self.logger.info("=" * 60)
        self.logger.info("开始LTR_FINDER流程|Starting LTR_FINDER pipeline")
        self.logger.info("=" * 60)

        output_prefix = self.config.genome_path.name
        output_file = self.config.output_path / f"{output_prefix}.finder.combine.scn"

        # 使用绝对路径
        # Use absolute path
        genome_file = str(self.config.genome_path)

        cmd = [
            self.config.ltr_finder_path,
            "-seq", genome_file,
            "-threads", str(self.config.threads),
            "-size", str(self.finder_config.size),
            "-time", str(self.finder_config.time)
        ]

        if self.finder_config.harvest_out:
            cmd.append("-harvest_out")

        try:
            # 检查并删除可能存在的符号链接|Check and remove potential symlink
            genome_link = self.config.output_path / self.config.genome_path.name
            if genome_link.exists() or genome_link.is_symlink():
                self.logger.info(f"删除已存在的符号链接|Removing existing symlink: {genome_link}")
                try:
                    os.remove(genome_link)
                except Exception as e:
                    self.logger.warning(f"删除符号链接失败|Failed to remove symlink: {e}")

            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")
            self.logger.info(f"工作目录|Working directory: {self.config.output_path}")

            # 在输出目录中运行|Run in output directory
            # 自动包装conda环境的命令|Auto-wrap conda environment commands
            cmd_name = os.path.basename(cmd[0])
            wrapped_cmd = build_conda_command(cmd_name, cmd[1:])

            result = subprocess.run(
                wrapped_cmd,
                cwd=str(self.config.output_path),
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                check=True
            )

            # LTR_FINDER_parallel会生成.combine.scn文件
            # LTR_FINDER_parallel generates .combine.scn file
            if not output_file.exists():
                self.logger.error(f"输出文件未生成|Output file not generated: {output_file}")
                return None

            self.logger.info(f"LTR_FINDER完成|LTR_FINDER completed: {output_file}")

            # 统计候选数量|Count candidates
            candidate_count = 0
            with open(output_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        candidate_count += 1

            self.logger.info(f"找到|Found {candidate_count} 个LTR候选|LTR candidates")

            self.logger.info("=" * 60)
            self.logger.info("LTR_FINDER流程完成|LTR_FINDER pipeline completed")
            self.logger.info("=" * 60)

            return output_file

        except subprocess.CalledProcessError as e:
            self.logger.error(f"LTR_FINDER失败|LTR_FINDER failed: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"运行LTR_FINDER时出错|Error running LTR_FINDER: {e}")
            return None
