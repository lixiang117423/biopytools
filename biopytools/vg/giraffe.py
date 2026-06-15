"""
VG Giraffe - 快速序列比对|VG Giraffe - Fast Read Alignment
"""

from pathlib import Path
from .config import VGGiraffeConfig
from .utils import run_vg_command, get_file_size_mb


class VGGiraffeRunner:
    """VG Giraffe运行器|VG Giraffe Runner"""

    def __init__(self, config: VGGiraffeConfig, logger):
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
        运行VG giraffe|Run VG giraffe

        Returns:
            是否成功|Whether successful
        """
        try:
            self.logger.info("")
            self.logger.info("开始序列比对|Starting read alignment")
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
                "vg giraffe",
                f"-x {self.config.graph}.xg",
                f"-H {self.config.graph}.min",
                f"-d {self.config.graph}.dist",
                f"-f {self.config.reads}",
                f"-t {self.config.threads}"
            ]

            # 双端测序|Paired-end reads
            if self.config.reads2:
                cmd_parts.append(f"-f {self.config.reads2}")

            # 片段长度参数|Fragment length parameters
            if self.config.fragment_length > 0:
                cmd_parts.append(f"-l {self.config.fragment_length}")

            if self.config.fragment_std_dev > 0:
                cmd_parts.append(f"-s {self.config.fragment_std_dev}")

            # 最小相似度|Min identity
            if self.config.min_identity > 0:
                cmd_parts.append(f"-m {self.config.min_identity}")

            # 输出格式|Output format
            if self.config.output_format == "GAF":
                cmd_parts.append("-a")
            else:
                cmd_parts.append("-G")  # GAM format

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
                self.logger
            )

            if success:
                self.logger.info("")
                self.logger.info("序列比对成功|Read alignment completed successfully")

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
                self.logger.error("序列比对失败|Read alignment failed")
                return False

        except Exception as e:
            self.logger.error(f"比对过程异常|Alignment exception: {e}")
            return False
