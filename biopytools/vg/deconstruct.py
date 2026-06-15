"""
VG Deconstruct - 从变异图生成VCF|VG Deconstruct - Export VCF from Graph
"""

from pathlib import Path
from .config import VGDeconstructConfig
from .utils import run_vg_command, get_file_size_mb


class VGDeconstructRunner:
    """VG Deconstruct运行器|VG Deconstruct Runner"""

    def __init__(self, config: VGDeconstructConfig, logger):
        """
        初始化运行器|Initialize runner

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def run(self) -> bool:
        """
        运行VG deconstruct|Run VG deconstruct

        Returns:
            是否成功|Whether successful
        """
        try:
            self.logger.info("")
            self.logger.info("开始导出VCF|Starting VCF export")
            self.logger.info("-" * 60)

            # 检查输出文件是否已存在|Check if output file already exists
            if Path(self.config.output).exists():
                self.logger.info("")
                self.logger.info("输出文件已存在，跳过|Output file already exists, skipping")
                file_size = get_file_size_mb(self.config.output)
                self.logger.info(f"  文件|File: {self.config.output} ({file_size:.2f} MB)")
                return True

            # 构建命令|Build command
            cmd_parts = [
                "vg deconstruct",
                f"-t {self.config.threads}",
                f"-a {self.config.reference_path}",
                f"-e",  # 输出所有等位基因|Output all alleles
            ]

            # 样本选择|Sample selection
            if self.config.samples:
                for sample in self.config.samples:
                    cmd_parts.append(f"-s {sample}")

            if self.config.progress:
                cmd_parts.append("-p")

            # 添加输入和输出文件|Add input and output files
            cmd_parts.append(self.config.input_graph)
            command = " ".join(cmd_parts) + f" > {self.config.output}"

            # 执行命令|Execute command
            self.logger.info("")
            self.logger.info(f"执行命令|Executing command:")
            self.logger.info(f"  {command}")
            self.logger.info("")

            success, stdout, stderr = run_vg_command(
                self.config.vg_env,
                command,
                self.logger
            )

            if success:
                self.logger.info("")
                self.logger.info("VCF导出成功|VCF export completed successfully")

                # 检查输出文件|Check output file
                if Path(self.config.output).exists():
                    file_size = get_file_size_mb(self.config.output)
                    self.logger.info(f"  输出文件|Output file: {self.config.output}")
                    self.logger.info(f"  文件大小|File size: {file_size:.2f} MB")
                else:
                    self.logger.warning("输出文件未找到|Output file not found")

                return True
            else:
                self.logger.error("")
                self.logger.error("VCF导出失败|VCF export failed")
                return False

        except Exception as e:
            self.logger.error(f"导出过程异常|Export exception: {e}")
            return False
