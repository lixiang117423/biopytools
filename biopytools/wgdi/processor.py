"""
WGDI核心处理器模块|WGDI Core Processor Module
"""

from pathlib import Path
from .config import WGDIConfig, DotPlotConfig, CollinearityConfig, CalKsConfig
from .utils import WGDILogger, CommandRunner, WGDIConfGenerator


class WGDIProcessor:
    """WGDI核心处理器|WGDI Core Processor"""

    def __init__(self, config: WGDIConfig, logger=None, cmd_runner=None):
        """
        初始化处理器|Initialize processor

        Args:
            config: 配置对象|Configuration object
            logger: 日志器（可选）|Logger (optional)
            cmd_runner: 命令执行器（可选）|Command runner (optional)
        """
        self.config = config

        # 初始化日志|Initialize logging
        if logger is None:
            log_file = Path(config.working_dir) / "wgdi.log"
            self.logger_manager = WGDILogger(str(log_file))
            self.logger = self.logger_manager.get_logger()
        else:
            self.logger = logger

        # 初始化命令执行器|Initialize command runner
        if cmd_runner is None:
            self.cmd_runner = CommandRunner(self.logger, config.working_dir)
        else:
            self.cmd_runner = cmd_runner

    def run_dotplot(self, config: DotPlotConfig) -> bool:
        """
        运行DotPlot分析|Run DotPlot analysis

        Args:
            config: DotPlot配置对象|DotPlot configuration object

        Returns:
            bool: 是否成功|Whether succeeded
        """
        self.logger.info("开始DotPlot分析|Starting DotPlot analysis")

        # 验证配置|Validate configuration
        config.validate()

        # 生成配置文件|Generate configuration file
        conf_generator = WGDIConfGenerator()
        conf_content = conf_generator.generate_dotplot_conf(config)

        conf_file = Path(config.working_dir) / "dotplot.conf"
        with open(conf_file, 'w') as f:
            f.write(conf_content)

        self.logger.info(f"配置文件已生成|Configuration file generated: {conf_file}")

        # 构建WGDI命令|Build WGDI command
        cmd = f"{self.config.wgdi_path} -d {conf_file}"

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, "DotPlot分析|DotPlot analysis")

        if success:
            self.logger.info(f"DotPlot分析完成|DotPlot analysis completed: {config.savefig}")

        return success

    def run_collinearity(self, config: CollinearityConfig) -> bool:
        """
        运行Collinearity分析|Run Collinearity analysis

        Args:
            config: Collinearity配置对象|Collinearity configuration object

        Returns:
            bool: 是否成功|Whether succeeded
        """
        self.logger.info("开始Collinearity分析|Starting Collinearity analysis")

        # 验证配置|Validate configuration
        config.validate()

        # 生成配置文件|Generate configuration file
        conf_generator = WGDIConfGenerator()
        conf_content = conf_generator.generate_collinearity_conf(config)

        conf_file = Path(config.working_dir) / "collinearity.conf"
        with open(conf_file, 'w') as f:
            f.write(conf_content)

        self.logger.info(f"配置文件已生成|Configuration file generated: {conf_file}")

        # 构建WGDI命令|Build WGDI command
        cmd = f"{self.config.wgdi_path} -icl {conf_file}"

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, "Collinearity分析|Collinearity analysis")

        if success:
            self.logger.info(f"Collinearity分析完成|Collinearity analysis completed: {config.savefile}")

        return success

    def run_calks(self, config: CalKsConfig) -> bool:
        """
        运行CalKs分析|Run CalKs analysis

        Args:
            config: CalKs配置对象|CalKs configuration object

        Returns:
            bool: 是否成功|Whether succeeded
        """
        self.logger.info("开始CalKs分析|Starting CalKs analysis")

        # 验证配置|Validate configuration
        config.validate()

        # 生成配置文件|Generate configuration file
        conf_generator = WGDIConfGenerator()
        conf_content = conf_generator.generate_calks_conf(config)

        conf_file = Path(config.working_dir) / "calks.conf"
        with open(conf_file, 'w') as f:
            f.write(conf_content)

        self.logger.info(f"配置文件已生成|Configuration file generated: {conf_file}")

        # 构建WGDI命令|Build WGDI command
        cmd = f"{self.config.wgdi_path} -ks {conf_file}"

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, "CalKs分析|CalKs analysis")

        if success:
            self.logger.info(f"CalKs分析完成|CalKs analysis completed: {config.savefile}")

        return success
