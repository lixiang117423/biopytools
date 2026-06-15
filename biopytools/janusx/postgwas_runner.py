"""
JanusX PostGWAS运行器|JanusX PostGWAS Runner

执行GWAS结果可视化和注释|Execute GWAS result visualization and annotation
"""

import os
from pathlib import Path
from typing import List

from .config import JanusXPostGWASConfig
from .utils import JanusXLogger, CommandRunner, check_janusx_dependencies


class JanusXPostGWASRunner:
    """JanusX PostGWAS运行器|JanusX PostGWAS Runner"""

    def __init__(self, config: JanusXPostGWASConfig):
        """初始化PostGWAS运行器|Initialize PostGWAS Runner

        Args:
            config: PostGWAS配置对象|PostGWAS configuration object
        """
        self.config = config
        self.logger_manager = JanusXLogger(
            config.output_path,
            log_file="janusx_postgwas.log"
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
        """构建PostGWAS命令|Build PostGWAS command

        Returns:
            命令列表|Command list
        """
        cmd = [self.config.janusx_path, "postGWAS"]

        # GWAS结果文件|GWAS result files
        for file in self.config.files:
            cmd.extend(["-f", file])

        # 可选参数|Optional parameters
        if self.config.chr_col != "#CHROM":
            cmd.extend(["-chr", self.config.chr_col])

        if self.config.pos_col != "POS":
            cmd.extend(["-pos", self.config.pos_col])

        if self.config.pvalue_col != "p":
            cmd.extend(["-pvalue", self.config.pvalue_col])

        if self.config.threshold is not None:
            cmd.extend(["-threshold", str(self.config.threshold)])

        if self.config.noplot:
            cmd.append("-noplot")

        if self.config.color != 0:
            cmd.extend(["-color", str(self.config.color)])

        if self.config.highlight:
            cmd.extend(["-hl", self.config.highlight])

        if self.config.format != "png":
            cmd.extend(["-format", self.config.format])

        if self.config.anno:
            cmd.extend(["-a", self.config.anno])

        if self.config.anno_broaden is not None:
            cmd.extend(["-ab", str(self.config.anno_broaden)])

        if self.config.desc_item != "description":
            cmd.extend(["-descItem", self.config.desc_item])

        # 输出参数|Output parameters
        cmd.extend(["-o", self.config.output_dir])
        cmd.extend(["-prefix", self.config.prefix])

        if self.config.threads != -1:
            cmd.extend(["-t", str(self.config.threads)])

        return cmd

    def run(self) -> bool:
        """运行PostGWAS分析|Run PostGWAS analysis

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始JanusX PostGWAS分析|Starting JanusX PostGWAS analysis")
        self.logger.info("=" * 60)

        # 输出配置信息|Output configuration info
        self.logger.info(f"GWAS结果文件|GWAS result files: {len(self.config.files)}")
        for file in self.config.files:
            self.logger.info(f"  - {file}")
        self.logger.info(f"染色体列名|Chromosome column: {self.config.chr_col}")
        self.logger.info(f"位置列名|Position column: {self.config.pos_col}")
        self.logger.info(f"P值列名|P-value column: {self.config.pvalue_col}")
        self.logger.info(f"显著性阈值|Significance threshold: {self.config.threshold}")
        self.logger.info(f"生成图表|Generate plots: {not self.config.noplot}")
        self.logger.info(f"输出格式|Output format: {self.config.format}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        try:
            # 构建命令|Build command
            cmd = self.build_command()
            self.logger.info(f"构建命令|Command built: {' '.join(cmd)}")

            # 执行命令|Execute command
            success = self.cmd_runner.run_command(
                cmd,
                description="运行PostGWAS分析|Running PostGWAS analysis"
            )

            if success:
                self.logger.info("=" * 60)
                self.logger.info("PostGWAS分析完成|PostGWAS analysis completed")
                self.logger.info("=" * 60)

                # 列出输出文件|List output files
                self._list_output_files()

                return True
            else:
                self.logger.error("PostGWAS分析失败|PostGWAS analysis failed")
                return False

        except Exception as e:
            self.logger.error(f"PostGWAS分析异常|PostGWAS analysis error: {e}")
            return False

    def _list_output_files(self):
        """列出输出文件|List output files"""
        self.logger.info("输出文件|Output files:")

        output_path = Path(self.config.output_dir)

        # 查找所有输出文件|Find all output files
        output_files = []

        # 查找图表文件|Find plot files
        if not self.config.noplot:
            for file in self.config.files:
                file_basename = Path(file).stem
                manh_file = output_path / f"{file_basename}.manh.{self.config.format}"
                qq_file = output_path / f"{file_basename}.qq.{self.config.format}"
                if manh_file.exists():
                    output_files.append(manh_file)
                if qq_file.exists():
                    output_files.append(qq_file)

        # 查找注释文件|Find annotation files
        if self.config.threshold is not None:
            anno_files = list(output_path.glob(f"*.{self.config.threshold}.anno.tsv"))
            output_files.extend(anno_files)

        if output_files:
            for file in sorted(output_files):
                size = file.stat().st_size if file.exists() else 0
                size_kb = size / 1024
                self.logger.info(f"  {file.name} ({size_kb:.2f} KB)")
        else:
            self.logger.warning("未找到输出文件|No output files found")
