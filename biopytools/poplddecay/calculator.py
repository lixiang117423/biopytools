"""
PopLDdecay计算模块|PopLDdecay Calculator Module

运行PopLDdecay软件进行LD衰减计算|Run PopLDdecay software for LD decay calculation
"""

import subprocess
from pathlib import Path

from .utils import build_conda_command


class PopLDdecayCalculator:
    """PopLDdecay计算器|PopLDdecay Calculator"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def run(self) -> bool:
        """运行PopLDdecay计算|Run PopLDdecay calculation

        Returns:
            bool: 是否成功|Success
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始LD衰减计算|Start LD decay calculation")
            self.logger.info("=" * 60)

            # 构建命令|Build command
            cmd = self._build_command()

            # 打印命令|Print command
            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            # 运行命令|Run command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            # 检查结果|Check result
            if result.returncode != 0:
                self.logger.error(f"PopLDdecay运行失败|PopLDdecay run failed")
                self.logger.error(f"错误信息|Error: {result.stderr}")
                return False

            # 输出结果|Output result
            if result.stdout:
                self.logger.info(f"输出|Output:\n{result.stdout}")

            # 检查输出文件|Check output files
            output_files = self.config.get_output_files()
            stat_file = Path(output_files['stat'])

            if not stat_file.exists():
                self.logger.error(f"输出文件未生成|Output file not generated: {stat_file}")
                return False

            self.logger.info("=" * 60)
            self.logger.info("LD衰减计算完成|LD decay calculation completed")
            self.logger.info("=" * 60)
            self.logger.info(f"统计文件|Statistics file: {stat_file}")

            return True

        except Exception as e:
            self.logger.error(f"运行PopLDdecay时出错|Error running PopLDdecay: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _build_command(self) -> list:
        """构建命令行|Build command line

        Returns:
            list: 命令列表|Command list
        """
        # 构建PopLDdecay参数|Build PopLDdecay arguments
        args = []

        # 输入文件|Input file
        if self.config.input_type == "vcf":
            args.extend(["-InVCF", str(self.config.input_path)])
        else:
            args.extend(["-InGenotype", str(self.config.input_path)])

        # 输出文件|Output file
        output_prefix = str(self.config.output_path)
        args.extend(["-OutStat", output_prefix])

        # 过滤参数|Filter parameters
        args.extend(["-MaxDist", str(self.config.max_dist)])
        args.extend(["-MAF", str(self.config.min_maf)])
        args.extend(["-Het", str(self.config.max_het)])
        args.extend(["-Miss", str(self.config.max_miss)])

        # 子群体|Subpopulation
        if self.config.subpop_file:
            args.extend(["-SubPop", str(self.config.subpop_path)])

        # 输出类型|Output type
        args.extend(["-OutType", str(self.config.out_type)])

        # 其他参数|Other parameters
        if self.config.ehh_site != "NA":
            args.extend(["-EHH", self.config.ehh_site])

        if self.config.output_filtered_snp:
            args.append("-OutFilterSNP")

        # 使用conda run包装命令|Wrap command with conda run
        cmd = build_conda_command(self.config.poplddecay_path, args)

        return cmd
