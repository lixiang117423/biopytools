"""
LTR_retriever运行器|LTR_retriever Runner
"""

import subprocess
import os
from pathlib import Path
from .utils import build_conda_command


class LTRRetrieverRunner:
    """LTR_retriever运行器类|LTR_retriever Runner Class"""

    def __init__(self, config, logger, retriever_config):
        """
        初始化LTR_retriever运行器|Initialize LTR_retriever runner

        Args:
            config: LAI配置对象|LAI configuration object
            logger: 日志器|Logger
            retriever_config: LTR_retriever配置对象|LTR_retriever configuration object
        """
        self.config = config
        self.logger = logger
        self.retriever_config = retriever_config

    def run(self, raw_ltr_file: Path) -> dict:
        """
        运行LTR_retriever筛选LTR|Run LTR_retriever to filter LTRs

        Args:
            raw_ltr_file: 原始LTR候选文件|Raw LTR candidate file

        Returns:
            dict: 输出文件路径字典|Output file paths dictionary
        """
        self.logger.info("=" * 60)
        self.logger.info("开始LTR_retriever流程|Starting LTR_retriever pipeline")
        self.logger.info("=" * 60)

        cmd = [
            "perl",
            self.config.ltr_retriever_path,
            "-genome", str(self.config.genome_path),
            "-inharvest", str(raw_ltr_file),
            "-threads", str(self.config.threads),
            "-outdir", str(self.config.output_path)
        ]

        try:
            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            # 自动包装conda环境的命令|Auto-wrap conda environment commands
            cmd_name = os.path.basename(cmd[0])  # perl
            wrapped_cmd = build_conda_command(cmd_name, cmd[1:])

            result = subprocess.run(
                wrapped_cmd,
                cwd=str(self.config.output_path),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )

            self.logger.info("LTR_retriever完成|LTR_retriever completed")

            # 检查输出文件|Check output files
            prefix = self.config.genome_path.name
            pass_list = self.config.output_path / f"{prefix}.pass.list"
            out_gff = self.config.output_path / f"{prefix}.out.gff"
            out_file = self.config.output_path / f"{prefix}.out"
            lai_file = self.config.output_path / f"{prefix}.out.LAI"

            outputs = {
                'pass_list': pass_list if pass_list.exists() else None,
                'out_gff': out_gff if out_gff.exists() else None,
                'out': out_file if out_file.exists() else None,
                'lai': lai_file if lai_file.exists() else None
            }

            # 统计完整LTR数量|Count intact LTRs
            if pass_list.exists():
                intact_count = 0
                with open(pass_list, 'r') as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            intact_count += 1
                self.logger.info(f"找到|Found {intact_count} 个完整LTR-RT|intact LTR-RTs")

            self.logger.info("=" * 60)
            self.logger.info("LTR_retriever流程完成|LTR_retriever pipeline completed")
            self.logger.info("=" * 60)

            return outputs

        except subprocess.CalledProcessError as e:
            self.logger.error(f"LTR_retriever失败|LTR_retriever failed: {e.stderr}")
            return {}
        except Exception as e:
            self.logger.error(f"运行LTR_retriever时出错|Error running LTR_retriever: {e}")
            return {}
