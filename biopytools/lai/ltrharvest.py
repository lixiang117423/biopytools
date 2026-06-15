"""
LTRharvest运行器|LTRharvest Runner
"""

import subprocess
import os
from pathlib import Path
from .utils import run_command, build_conda_command


class LTRHarvestRunner:
    """LTRharvest运行器类|LTRharvest Runner Class"""

    def __init__(self, config, logger, harvest_config):
        """
        初始化LTRharvest运行器|Initialize LTRharvest runner

        Args:
            config: LAI配置对象|LAI configuration object
            logger: 日志器|Logger
            harvest_config: LTRharvest配置对象|LTRharvest configuration object
        """
        self.config = config
        self.logger = logger
        self.harvest_config = harvest_config

    def run_index(self) -> bool:
        """
        运行suffixerator建立索引|Run suffixerator to build index

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤1/2: 建立基因组索引|Step 1/2: Building genome index")

        index_name = str(self.config.genome_path)

        cmd = [
            self.config.gt_path,
            "suffixerator",
            "-db", str(self.config.genome_path),
            "-indexname", index_name,
            "-tis",
            "-suf",
            "-lcp",
            "-des",
            "-ssp",
            "-sds",
            "-dna",
            "-memlimit", f"{self.harvest_config.index_memory}MB"
        ]

        try:
            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")
            result = run_command(cmd, self.logger, check=True)

            self.logger.info("索引建立完成|Index building completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"索引建立失败|Index building failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"建立索引时出错|Error building index: {e}")
            return False

    def run_harvest(self) -> Path:
        """
        运行LTRharvest识别LTR候选|Run LTRharvest to identify LTR candidates

        Returns:
            Path: 输出文件路径|Output file path
        """
        self.logger.info("步骤2/2: 运行LTRharvest识别LTR候选|Step 2/2: Running LTRharvest to identify LTR candidates")

        output_file = self.config.output_path / f"{self.config.genome_path.name}.harvest.scn"

        index_name = str(self.config.genome_path)

        cmd = [
            self.config.gt_path,
            "ltrharvest",
            "-index", index_name,
            "-minlenltr", str(self.harvest_config.minlenltr),
            "-maxlenltr", str(self.harvest_config.maxlenltr),
            "-mintsd", str(self.harvest_config.mintsd),
            "-maxtsd", str(self.harvest_config.maxtsd),
            "-motif", self.harvest_config.motif,
            "-motifmis", str(self.harvest_config.motifmis),
            "-similar", str(self.harvest_config.similar),
            "-vic", str(self.harvest_config.vic),
            "-seed", str(self.harvest_config.seed),
            "-seqids", "yes"
        ]

        try:
            self.logger.info(f"运行命令|Running command: {' '.join(cmd)} > {output_file}")

            # 使用重定向输出到文件|Use redirection to output to file
            # 自动包装conda环境的命令|Auto-wrap conda environment commands
            cmd_name = os.path.basename(cmd[0])
            wrapped_cmd = build_conda_command(cmd_name, cmd[1:])

            with open(output_file, 'w') as f:
                result = subprocess.run(
                    wrapped_cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True
                )

            if not output_file.exists():
                self.logger.error(f"输出文件未生成|Output file not generated: {output_file}")
                return None

            self.logger.info(f"LTRharvest完成|LTRharvest completed: {output_file}")

            # 统计候选数量|Count candidates
            candidate_count = 0
            with open(output_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        candidate_count += 1

            self.logger.info(f"找到|Found {candidate_count} 个LTR候选|LTR candidates")

            return output_file

        except subprocess.CalledProcessError as e:
            self.logger.error(f"LTRharvest失败|LTRharvest failed: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"运行LTRharvest时出错|Error running LTRharvest: {e}")
            return None

    def run(self) -> Path:
        """
        运行完整的LTRharvest流程|Run complete LTRharvest pipeline

        Returns:
            Path: 输出文件路径|Output file path
        """
        self.logger.info("=" * 60)
        self.logger.info("开始LTRharvest流程|Starting LTRharvest pipeline")
        self.logger.info("=" * 60)

        # 步骤1: 建立索引|Step 1: Build index
        if not self.run_index():
            return None

        # 步骤2: 运行LTRharvest|Step 2: Run LTRharvest
        output_file = self.run_harvest()
        if not output_file:
            return None

        self.logger.info("=" * 60)
        self.logger.info("LTRharvest流程完成|LTRharvest pipeline completed")
        self.logger.info("=" * 60)

        return output_file
