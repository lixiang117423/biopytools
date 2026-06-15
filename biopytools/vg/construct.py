"""
VG Construct - 变异图构建|VG Construct - Graph Construction
"""

from pathlib import Path
from .config import VGConstructConfig
from .utils import run_vg_command, get_file_size_mb


class VGConstructRunner:
    """VG构建运行器|VG Construction Runner"""

    def __init__(self, config: VGConstructConfig, logger):
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
        运行VG construct|Run VG construct

        Returns:
            是否成功|Whether successful
        """
        try:
            self.logger.info("")
            self.logger.info("开始构建变异图|Starting variation graph construction")
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
                "vg construct",
                f"-r {self.config.reference}",
                f"-v {self.config.vcf}",
                f"-t {self.config.threads}",
                f"--node-max {self.config.node_max}"
            ]

            # 可选参数|Optional parameters
            if self.config.region:
                cmd_parts.append(f"-R {self.config.region}")

            if self.config.alt_paths:
                cmd_parts.append("-a")

            if self.config.progress:
                cmd_parts.append("-p")

            # 重定向输出|Redirect output
            command = " ".join(cmd_parts) + f" > {self.config.output}"

            # 执行命令|Execute command
            self.logger.info("")
            self.logger.info(f"执行命令|Executing command:")
            self.logger.info(f"  {command}")
            self.logger.info("")

            success, stdout, stderr = run_vg_command(
                self.config.vg_env,
                command,
                self.logger,
                cwd=str(Path(self.config.output).parent)
            )

            if success:
                self.logger.info("")
                self.logger.info("变异图构建成功|Graph construction completed successfully")

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
                self.logger.error("变异图构建失败|Graph construction failed")
                return False

        except Exception as e:
            self.logger.error(f"构建过程异常|Construction exception: {e}")
            return False
