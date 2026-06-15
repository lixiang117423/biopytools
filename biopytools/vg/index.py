"""
VG Index - 图索引创建|VG Index - Graph Index Creation
"""

from pathlib import Path
from .config import VGIndexConfig
from .utils import run_vg_command, get_file_size_mb


class VGIndexRunner:
    """VG索引运行器|VG Index Runner"""

    def __init__(self, config: VGIndexConfig, logger):
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
        运行VG index|Run VG index

        Returns:
            是否成功|Whether successful
        """
        try:
            self.logger.info("")
            self.logger.info("开始创建图索引|Starting graph indexing")
            self.logger.info("-" * 60)

            # 确定要创建的索引类型|Determine index types to create
            index_types = []
            if self.config.xg:
                index_types.append('xg')
            if self.config.gcsa:
                index_types.append('gcsa')
            if self.config.gbwt:
                index_types.append('gbwt')
            if self.config.giraffe:
                index_types.extend(['xg', 'giraffe'])

            if not index_types:
                self.logger.error("必须指定至少一种索引类型|Must specify at least one index type")
                return False

            self.logger.info(f"索引类型|Index types: {', '.join(index_types)}")

            # 构建命令|Build command
            cmd_parts = [
                "vg index",
                f"-t {self.config.threads}",
                f"-i {self.config.input_graph}"
            ]

            # XG索引|XG index
            if self.config.xg or self.config.giraffe:
                xg_file = f"{self.config.output_prefix}.xg"
                if Path(xg_file).exists():
                    self.logger.info(f"XG索引已存在，跳过|XG index exists, skipping: {xg_file}")
                else:
                    cmd_parts.append(f"-x {xg_file}")

            # GCSA索引|GCSA index
            if self.config.gcsa:
                gcsa_file = f"{self.config.output_prefix}.gcsa"
                if Path(gcsa_file).exists():
                    self.logger.info(f"GCSA索引已存在，跳过|GCSA index exists, skipping: {gcsa_file}")
                else:
                    cmd_parts.append(f"-g {gcsa_file}")
                    cmd_parts.append(f"-k {self.config.kmer_size}")
                    cmd_parts.append(f"-E {self.config.edge_max}")

            # GBWT索引|GBWT index
            if self.config.gbwt:
                gbwt_file = f"{self.config.output_prefix}.gbwt"
                if Path(gbwt_file).exists():
                    self.logger.info(f"GBWT索引已存在，跳过|GBWT index exists, skipping: {gbwt_file}")
                else:
                    cmd_parts.append(f"-b {gbwt_file}")

            # GIRAFFE索引|GIRAFFE indexes
            if self.config.giraffe:
                min_file = f"{self.config.output_prefix}.min"
                dist_file = f"{self.config.output_prefix}.dist"

                if Path(min_file).exists() and Path(dist_file).exists():
                    self.logger.info(f"GIRAFFE索引已存在，跳过|GIRAFFE indexes exist, skipping")
                else:
                    cmd_parts.append(f"-m {min_file}")
                    cmd_parts.append(f"-d {dist_file}")

            if self.config.progress:
                cmd_parts.append("-p")

            command = " ".join(cmd_parts)

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
                self.logger.info("图索引创建成功|Graph indexing completed successfully")

                # 列出创建的索引文件|List created index files
                self.logger.info("")
                self.logger.info("创建的索引文件|Created index files:")

                for index_type in index_types:
                    if index_type == 'xg':
                        self._check_and_report(f"{self.config.output_prefix}.xg")
                    elif index_type == 'gcsa':
                        self._check_and_report(f"{self.config.output_prefix}.gcsa")
                    elif index_type == 'gbwt':
                        self._check_and_report(f"{self.config.output_prefix}.gbwt")
                    elif index_type == 'giraffe':
                        self._check_and_report(f"{self.config.output_prefix}.min")
                        self._check_and_report(f"{self.config.output_prefix}.dist")

                return True
            else:
                self.logger.error("")
                self.logger.error("图索引创建失败|Graph indexing failed")
                return False

        except Exception as e:
            self.logger.error(f"索引创建异常|Indexing exception: {e}")
            return False

    def _check_and_report(self, file_path: str):
        """检查并报告文件|Check and report file"""
        if Path(file_path).exists():
            file_size = get_file_size_mb(file_path)
            self.logger.info(f"  ✓ {Path(file_path).name} ({file_size:.2f} MB)")
        else:
            self.logger.warning(f"  ✗ {Path(file_path).name} (未找到|not found)")
