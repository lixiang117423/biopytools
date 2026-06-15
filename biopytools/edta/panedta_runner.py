"""
PanEDTA泛基因组转座子注释运行器|PanEDTA Pan-genome TE Annotation Runner
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Optional

from .config import PanEDTAConfig
from .utils import EDTALogger, CommandRunner


class PanEDTARunner:
    """PanEDTA泛基因组转座子注释运行器|PanEDTA Pan-genome TE Annotation Runner"""

    def __init__(self, config: PanEDTAConfig):
        """初始化PanEDTA运行器|Initialize PanEDTA Runner

        Args:
            config: PanEDTA配置对象|PanEDTA configuration object
        """
        self.config = config
        self.logger_manager = EDTALogger(config.output_path, log_name="panedta_analysis.log")
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger, config.output_path)

    def build_command(self) -> list:
        """构建panEDTA命令|Build panEDTA command

        Returns:
            panEDTA命令列表|panEDTA command list
        """
        cmd = []

        # 查找panEDTA.sh|Find panEDTA.sh
        panedta_sh = self._find_panedta_sh()

        if not panedta_sh or not os.path.exists(panedta_sh):
            raise RuntimeError(f"未找到panEDTA.sh|panEDTA.sh not found at: {panedta_sh}")

        # 使用bash执行panEDTA.sh|Use bash to execute panEDTA.sh
        cmd = ["bash", panedta_sh]

        # 必需参数|Required parameters
        cmd.extend(["-g", self.config.genome_list])

        # 可选参数|Optional parameters
        if self.config.cds:
            cmd.extend(["-c", self.config.cds])
        if self.config.curatedlib:
            cmd.extend(["-l", self.config.curatedlib])
        cmd.extend(["-f", str(self.config.fl_copy)])
        cmd.extend(["-a", str(self.config.anno)])
        cmd.extend(["-o", str(self.config.overwrite)])
        cmd.extend(["-t", str(self.config.threads)])

        return cmd

    def _find_panedta_sh(self) -> Optional[str]:
        """查找panEDTA.sh脚本|Find panEDTA.sh script

        Returns:
            panEDTA.sh路径或None|panEDTA.sh path or None
        """
        # 1. 检查配置中的edta_path|Check edta_path in config
        if self.config.edta_path:
            panedta_sh = os.path.join(self.config.edta_path, "panEDTA.sh")
            if os.path.exists(panedta_sh):
                return panedta_sh

        # 2. 检查conda环境|Check conda environment
        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if conda_prefix:
            panedta_sh = os.path.join(conda_prefix, 'share', 'EDTA', 'panEDTA.sh')
            if os.path.exists(panedta_sh):
                return panedta_sh

        # 3. 使用默认路径|Use default path
        default_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA/panEDTA.sh'
        if os.path.exists(default_path):
            return default_path

        return None

    def run(self) -> bool:
        """运行PanEDTA分析|Run PanEDTA analysis

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始PanEDTA泛基因组转座子注释分析|Starting PanEDTA pan-genome TE annotation analysis")
        self.logger.info("=" * 60)

        try:
            # 验证基因组列表|Validate genome list
            self._validate_genome_list()

            # 构建命令|Build command
            cmd = self.build_command()

            # 打印参数信息|Print parameter information
            self._log_parameters()

            # 执行PanEDTA|Execute PanEDTA
            self.logger.info("开始执行panEDTA.sh|Starting panEDTA.sh execution")
            result = self.cmd_runner.run(
                cmd,
                description="PanEDTA泛基因组转座子注释|PanEDTA pan-genome TE annotation",
                check=True,
                show_progress=True
            )

            # 检查结果|Check results
            if result.returncode == 0:
                self.logger.info("PanEDTA分析完成|PanEDTA analysis completed")
                self._check_outputs()
                return True
            else:
                self.logger.error(f"PanEDTA分析失败|PanEDTA analysis failed with return code: {result.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"PanEDTA运行出错|PanEDTA runtime error: {str(e)}")
            return False

    def _validate_genome_list(self):
        """验证基因组列表文件|Validate genome list file"""
        self.logger.info("验证基因组列表文件|Validating genome list file")

        if not os.path.exists(self.config.genome_list):
            raise FileNotFoundError(f"基因组列表文件不存在|Genome list file not found: {self.config.genome_list}")

        with open(self.config.genome_list, 'r') as f:
            lines = f.readlines()

        genome_count = 0
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 1:
                continue

            genome_file = parts[0]
            if not os.path.exists(genome_file):
                self.logger.warning(f"基因组文件不存在|Genome file not found: {genome_file}")
            else:
                genome_count += 1

        self.logger.info(f"共找到{genome_count}个有效基因组文件|Found {genome_count} valid genome files")

        if genome_count == 0:
            raise ValueError("未找到有效的基因组文件|No valid genome files found")

    def _log_parameters(self):
        """记录参数信息|Log parameter information"""
        self.logger.info("PanEDTA运行参数|PanEDTA running parameters:")
        self.logger.info(f"  基因组列表|Genome list: {self.config.genome_list}")
        self.logger.info(f"  CDS文件|CDS: {self.config.cds if self.config.cds else 'None'}")
        self.logger.info(f"  筛选库|Curated library: {self.config.curatedlib if self.config.curatedlib else 'None'}")
        self.logger.info(f"  拷贝数阈值|Copy number cutoff: {self.config.fl_copy}")
        self.logger.info(f"  注释|Annotation: {self.config.anno}")
        self.logger.info(f"  覆盖|Overwrite: {self.config.overwrite}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")

    def _check_outputs(self):
        """检查输出文件|Check output files"""
        self.logger.info("检查PanEDTA输出文件|Checking PanEDTA output files")

        # PanEDTA的输出文件检查较为复杂，这里只检查目录结构
        output_path = Path(self.config.output_dir)

        if output_path.exists():
            self.logger.info(f"输出目录存在|Output directory exists: {output_path}")

            # 列出输出目录中的文件|List files in output directory
            files = list(output_path.glob("*"))
            self.logger.info(f"输出目录中共有{len(files)}个文件/目录|Found {len(files)} files/directories in output directory")
        else:
            self.logger.warning(f"输出目录不存在|Output directory not found: {output_path}")
