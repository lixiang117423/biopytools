"""
LAI计算器|LAI Calculator
"""

import subprocess
import os
from pathlib import Path
from .utils import build_conda_command


class LAICalculatorRunner:
    """LAI计算器类|LAI Calculator Runner Class"""

    def __init__(self, config, logger, lai_config):
        """
        初始化LAI计算器|Initialize LAI calculator

        Args:
            config: LAI配置对象|LAI configuration object
            logger: 日志器|Logger
            lai_config: LAI计算配置对象|LAI calculation configuration object
        """
        self.config = config
        self.logger = logger
        self.lai_config = lai_config

    def run(self, pass_list: Path, out_file: Path) -> Path:
        """
        运行LAI计算|Run LAI calculation

        Args:
            pass_list: 完整LTR-RT列表文件|Intact LTR-RT list file
            out_file: RepeatMasker注释文件|RepeatMasker annotation file

        Returns:
            Path: LAI输出文件路径|LAI output file path
        """
        self.logger.info("=" * 60)
        self.logger.info("开始LAI计算|Starting LAI calculation")
        self.logger.info("=" * 60)

        cmd = [
            "perl",
            self.config.lai_path,
            "-genome", str(self.config.genome_path),
            "-intact", str(pass_list),
            "-all", str(out_file),
            "-window", str(self.lai_config.window),
            "-step", str(self.lai_config.step),
            "-t", str(self.lai_config.threads)
        ]

        # 添加可选参数|Add optional parameters
        if self.lai_config.quick:
            cmd.append("-q")
        if self.lai_config.qquick:
            cmd.append("-qq")
        if self.lai_config.mono:
            cmd.extend(["-mono", str(self.lai_config.mono)])
        if self.lai_config.iden is not None:
            cmd.extend(["-iden", str(self.lai_config.iden)])
        if self.lai_config.totltr is not None:
            cmd.extend(["-totLTR", str(self.lai_config.totltr)])
        if self.lai_config.genome_size is not None:
            cmd.extend(["-genome_size", str(self.lai_config.genome_size)])

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

            self.logger.info("LAI计算完成|LAI calculation completed")

            # 检查输出文件|Check output file
            lai_file = self.config.output_path / f"{out_file.name}.LAI"

            if not lai_file.exists():
                self.logger.error(f"LAI文件未生成|LAI file not generated: {lai_file}")
                return None

            # 解析LAI结果|Parse LAI result
            self._parse_lai_result(lai_file)

            self.logger.info("=" * 60)
            self.logger.info("LAI计算流程完成|LAI calculation pipeline completed")
            self.logger.info("=" * 60)

            return lai_file

        except subprocess.CalledProcessError as e:
            self.logger.error(f"LAI计算失败|LAI calculation failed: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"运行LAI计算时出错|Error running LAI calculation: {e}")
            return None

    def _parse_lai_result(self, lai_file: Path):
        """
        解析LAI结果文件|Parse LAI result file

        Args:
            lai_file: LAI结果文件|LAI result file
        """
        self.logger.info("LAI结果|LAI results:")
        self.logger.info("-" * 60)

        try:
            with open(lai_file, 'r') as f:
                for line in f:
                    if line.strip():
                        self.logger.info(line.strip())

            # 提取全基因组LAI值|Extract whole-genome LAI value
            with open(lai_file, 'r') as f:
                for line in f:
                    if "whole_genome" in line:
                        parts = line.split()
                        if len(parts) >= 7:
                            raw_lai = parts[5]
                            lai = parts[6]
                            self.logger.info(f"全基因组Raw LAI|Whole-genome Raw LAI: {raw_lai}")
                            self.logger.info(f"全基因组LAI|Whole-genome LAI: {lai}")
                        break

        except Exception as e:
            self.logger.error(f"解析LAI结果时出错|Error parsing LAI result: {e}")
