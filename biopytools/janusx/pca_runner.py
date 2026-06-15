"""
JanusX PCA运行器|JanusX PCA Runner

执行主成分分析|Execute principal component analysis
"""

import os
from pathlib import Path
from typing import List

from .config import JanusXPCAConfig
from .utils import JanusXLogger, CommandRunner, check_janusx_dependencies


class JanusXPCARunner:
    """JanusX PCA运行器|JanusX PCA Runner"""

    def __init__(self, config: JanusXPCAConfig):
        """初始化PCA运行器|Initialize PCA Runner

        Args:
            config: PCA配置对象|PCA configuration object
        """
        self.config = config
        self.logger_manager = JanusXLogger(
            config.output_path,
            log_file="janusx_pca.log"
        )
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger, config.output_path)

        # 检查JanusX依赖|Check JanusX dependencies
        try:
            check_janusx_dependencies(config.janusx_path, self.logger)
        except RuntimeError as e:
            self.logger.error(f"依赖检查失败|Dependency check failed: {e}")
            raise

    def build_command(self) -> List[str]:
        """构建PCA命令|Build PCA command

        Returns:
            命令列表|Command list
        """
        cmd = [self.config.janusx_path, "pca"]

        # 输入源优先级|Input source priority: grm > pcfile > genotype
        if self.config.grm:
            cmd.extend(["-grm", self.config.grm])
        elif self.config.pcfile:
            cmd.extend(["-pcfile", self.config.pcfile])
        else:
            # 基因型输入|Genotype input
            if self.config.genotype_type == "vcf":
                cmd.extend(["-vcf", self.config.genotype])
            else:  # bfile
                cmd.extend(["-bfile", self.config.genotype])

        # 可选参数|Optional parameters
        if self.config.dim != 3:
            cmd.extend(["-dim", str(self.config.dim)])

        if self.config.plot:
            cmd.append("-plot")

        if self.config.plot3d:
            cmd.append("-plot3D")

        if self.config.group:
            cmd.extend(["-group", self.config.group])

        if self.config.color != 1:
            cmd.extend(["-color", str(self.config.color)])

        # 输出参数|Output parameters
        cmd.extend(["-o", self.config.output_dir])

        if self.config.prefix:
            cmd.extend(["-prefix", self.config.prefix])

        return cmd

    def run(self) -> bool:
        """运行PCA分析|Run PCA analysis

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始JanusX PCA分析|Starting JanusX PCA analysis")
        self.logger.info("=" * 60)

        # 输出配置信息|Output configuration info
        if self.config.grm:
            self.logger.info(f"输入来源|Input source: GRM文件 ({self.config.grm})")
        elif self.config.pcfile:
            self.logger.info(f"输入来源|Input source: PCA文件 ({self.config.pcfile})")
        else:
            self.logger.info(f"输入来源|Input source: 基因型文件 ({self.config.genotype})")

        self.logger.info(f"主成分数量|Number of PCs: {self.config.dim}")
        self.logger.info(f"生成2D图|Generate 2D plot: {self.config.plot}")
        self.logger.info(f"生成3D图|Generate 3D plot: {self.config.plot3d}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        try:
            # 构建命令|Build command
            cmd = self.build_command()
            self.logger.info(f"构建命令|Command built: {' '.join(cmd)}")

            # 执行命令|Execute command
            success = self.cmd_runner.run_command(
                cmd,
                description="运行PCA分析|Running PCA analysis"
            )

            if success:
                self.logger.info("=" * 60)
                self.logger.info("PCA分析完成|PCA analysis completed")
                self.logger.info("=" * 60)

                # 列出输出文件|List output files
                self._list_output_files()

                return True
            else:
                self.logger.error("PCA分析失败|PCA analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"PCA分析异常|PCA analysis error: {e}")
            return False

    def _list_output_files(self):
        """列出输出文件|List output files"""
        self.logger.info("输出文件|Output files:")

        output_path = Path(self.config.output_dir)
        prefix = self.config.prefix

        # 查找所有输出文件|Find all output files
        output_files = list(output_path.glob(f"{prefix}.*"))

        if output_files:
            for file in sorted(output_files):
                size = file.stat().st_size if file.exists() else 0
                size_kb = size / 1024
                self.logger.info(f"  {file.name} ({size_kb:.2f} KB)")
        else:
            self.logger.warning("未找到输出文件|No output files found")
