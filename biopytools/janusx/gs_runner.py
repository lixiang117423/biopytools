"""
JanusX基因组选择运行器|JanusX Genomic Selection Runner

执行基因组选择分析|Execute genomic selection analysis
"""

import os
from pathlib import Path
from typing import List

from .config import JanusXGSConfig
from .utils import JanusXLogger, CommandRunner, check_janusx_dependencies


class JanusXGSRunner:
    """JanusX基因组选择运行器|JanusX Genomic Selection Runner"""

    def __init__(self, config: JanusXGSConfig):
        """初始化GS运行器|Initialize GS Runner

        Args:
            config: GS配置对象|GS configuration object
        """
        self.config = config
        self.logger_manager = JanusXLogger(
            config.output_path,
            log_file="janusx_gs.log"
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
        """构建GS命令|Build GS command

        Returns:
            命令列表|Command list
        """
        cmd = [self.config.janusx_path, "gs"]

        # 基因型输入|Genotype input
        if self.config.genotype_type == "vcf":
            cmd.extend(["-vcf", self.config.genotype])
        else:  # bfile
            cmd.extend(["-bfile", self.config.genotype])

        # 表型文件|Phenotype file
        cmd.extend(["-p", self.config.pheno])

        # GS模型|GS models
        for model in self.config.models:
            if model == "GBLUP":
                cmd.append("-GBLUP")
            elif model == "rrBLUP":
                cmd.append("-rrBLUP")
            elif model == "BayesA":
                cmd.append("-BayesA")
            elif model == "BayesB":
                cmd.append("-BayesB")
            elif model == "BayesCpi":
                cmd.append("-BayesCpi")

        # 可选参数|Optional parameters
        if self.config.maf != 0.02:
            cmd.extend(["-maf", str(self.config.maf)])

        if self.config.geno != 0.05:
            cmd.extend(["-geno", str(self.config.geno)])

        if self.config.pcd:
            cmd.append("-pcd")

        if self.config.ncol is not None:
            cmd.extend(["-n", str(self.config.ncol)])

        if self.config.cv is not None:
            cmd.extend(["-cv", str(self.config.cv)])

        if self.config.plot:
            cmd.append("-plot")

        # 输出参数|Output parameters
        cmd.extend(["-o", self.config.output_dir])

        if self.config.prefix:
            cmd.extend(["-prefix", self.config.prefix])

        return cmd

    def run(self) -> bool:
        """运行GS分析|Run GS analysis

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始JanusX基因组选择分析|Starting JanusX GS analysis")
        self.logger.info("=" * 60)

        # 输出配置信息|Output configuration info
        self.logger.info(f"基因型文件|Genotype file: {self.config.genotype}")
        self.logger.info(f"表型文件|Phenotype file: {self.config.pheno}")
        self.logger.info(f"GS模型|GS models: {', '.join(self.config.models)}")
        self.logger.info(f"MAF阈值|MAF threshold: {self.config.maf}")
        self.logger.info(f"缺失率阈值|Missing rate threshold: {self.config.geno}")
        self.logger.info(f"PCA降维|PCA dimensionality reduction: {self.config.pcd}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        try:
            # 构建命令|Build command
            cmd = self.build_command()
            self.logger.info(f"构建命令|Command built: {' '.join(cmd)}")

            # 执行命令|Execute command
            success = self.cmd_runner.run_command(
                cmd,
                description="运行GS分析|Running GS analysis"
            )

            if success:
                self.logger.info("=" * 60)
                self.logger.info("GS分析完成|GS analysis completed")
                self.logger.info("=" * 60)

                # 列出输出文件|List output files
                self._list_output_files()

                return True
            else:
                self.logger.error("GS分析失败|GS analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"GS分析异常|GS analysis error: {e}")
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
                size_mb = size / (1024 * 1024)
                self.logger.info(f"  {file.name} ({size_mb:.2f} MB)")
        else:
            self.logger.warning("未找到输出文件|No output files found")
