"""
HiTE 单基因组运行器|HiTE Single-genome Runner

通过 singularity 直接挂载调用 HiTE,核心运行逻辑
Calls HiTE via singularity direct-mount; core run logic
"""

import os
from typing import List
from .utils import (
    build_singularity_command,
    run_command,
    validate_singularity,
    validate_hite_sif,
)


class HiteRunner:
    """
    HiTE 运行器|HiTE Runner

    构建 HiTE 参数,通过 singularity exec 直接挂载运行,支持断点续传
    Builds HiTE args, runs via singularity exec direct-mount, supports checkpoint resume
    """

    # HiTE main.py 默认参数(用于判断是否追加高级参数)
    # HiTE main.py defaults (to decide whether to append advanced params)
    _DEFAULT_CHUNK_SIZE = 400
    _DEFAULT_MIU = 1.3e-8
    _DEFAULT_MIN_TE_LEN = 80
    _DEFAULT_TE_TYPE = "all"

    def __init__(self, config, logger):
        """
        初始化运行器|Initialize runner

        Args:
            config: HiteConfig 实例|HiteConfig instance
            logger: 日志器|logger
        """
        self.config = config
        self.logger = logger

    def _build_hite_args(self) -> List[str]:
        """
        构建 HiTE main.py 参数列表|Build HiTE main.py argument list

        将 HiteConfig 转为 HiTE 命令行参数;高级参数仅在非默认值时追加
        Convert HiteConfig to HiTE CLI args; advanced params only when non-default
        """
        c = self.config
        args: List[str] = []

        # 核心参数:genome/out_dir/work_dir 都用主机绝对路径(已 identity 挂载)
        # Core args: host absolute paths (identity-mounted)
        args.extend(["--genome", c.genome])
        args.extend(["--out_dir", c.hite_out_dir])
        args.extend(["--work_dir", c.work_dir])
        args.extend(["--thread", str(c.threads)])

        # 布尔参数|Boolean params (HiTE 用 1/0)
        args.extend(["--plant", "1" if c.plant else "0"])
        args.extend(["--annotate", "1" if c.annotate else "0"])
        args.extend(["--recover", "1" if c.recover else "0"])
        args.extend(["--domain", "1" if c.domain else "0"])
        args.extend(["--remove_nested", "1" if c.remove_nested else "0"])

        # 高级参数:仅非默认值时追加|Advanced params only when non-default
        if c.chunk_size != self._DEFAULT_CHUNK_SIZE:
            args.extend(["--chunk_size", str(c.chunk_size)])
        if c.miu != self._DEFAULT_MIU:
            args.extend(["--miu", str(c.miu)])
        if c.min_te_len != self._DEFAULT_MIN_TE_LEN:
            args.extend(["--min_TE_len", str(c.min_te_len)])
        if c.te_type != self._DEFAULT_TE_TYPE:
            args.extend(["--te_type", c.te_type])

        # 可选 curated 库|Optional curated library
        if c.curated_lib:
            args.extend(["--curated_lib", c.curated_lib])

        # debug 模式|debug mode
        if c.debug:
            args.extend(["--debug", "1"])

        return args

    def _build_command(self) -> List[str]:
        """构建完整 singularity 命令|Build full singularity command"""
        c = self.config
        # identity 挂载:基因组目录 + HiTE 输出目录 + /tmp
        # identity mount: genome dir + HiTE output dir + /tmp
        genome_dir = os.path.dirname(c.genome) or "."
        bind_paths = [genome_dir, c.hite_out_dir, "/tmp"]
        internal_cmd = ["python", "/HiTE/main.py"] + self._build_hite_args()
        return build_singularity_command(
            singularity_path=c.singularity_path,
            sif_file=c.sif_file,
            bind_paths=bind_paths,
            internal_cmd=internal_cmd,
        )

    def _is_completed(self) -> bool:
        """
        检查主结果是否已存在(模块级断点续传)|Check if main output exists (module-level checkpoint)

        通过 01_hite/confident_TE.cons.fa 是否存在判断
        """
        from pathlib import Path
        main_output = os.path.join(self.config.hite_out_dir, "confident_TE.cons.fa")
        return Path(main_output).exists()

    def _validate_environment(self) -> bool:
        """验证 singularity 与 SIF 可用|Validate singularity and SIF availability"""
        if not validate_singularity(self.config.singularity_path, self.logger):
            return False
        if not validate_hite_sif(self.config.sif_file, self.logger):
            return False
        return True

    def run(self) -> bool:
        """
        运行 HiTE 分析|Run HiTE analysis

        模块级断点续传:主结果已存在则跳过;否则 singularity exec 调用 HiTE

        Returns:
            成功返回 True|True on success

        Raises:
            RuntimeError: 环境验证失败|if environment validation fails
            subprocess.CalledProcessError: HiTE 执行失败|if HiTE execution fails
        """
        # 模块级断点续传|Module-level checkpoint
        if self._is_completed():
            self.logger.info(
                f"跳过已完成步骤|Skipping completed step: "
                f"{self.config.hite_out_dir}/confident_TE.cons.fa 已存在|already exists"
            )
            return True

        if not self._validate_environment():
            raise RuntimeError(
                "环境验证失败|Environment validation failed "
                "(singularity 或 SIF 不可用|singularity or SIF unavailable)"
            )

        self.logger.info("=" * 80)
        self.logger.info("开始 HiTE 转座子检测|Starting HiTE TE detection")
        self.logger.info("=" * 80)

        cmd = self._build_command()
        run_command(
            cmd, self.logger,
            description="HiTE 转座子检测与注释|HiTE TE detection and annotation",
        )

        self.logger.info("=" * 80)
        self.logger.info("HiTE 检测完成|HiTE detection completed")
        self.logger.info("=" * 80)
        return True
