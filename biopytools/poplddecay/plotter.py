"""
PopLDdecay绘图模块|PopLDdecay Plotter Module

使用Plot_OnePop.pl或Plot_MultiPop.pl绘制LD衰减图|Plot LD decay figure using Plot_OnePop.pl or Plot_MultiPop.pl
"""

import subprocess
from pathlib import Path
from typing import Optional, List


class PopLDdecayPlotter:
    """PopLDdecay绘图器|PopLDdecay Plotter"""

    def __init__(self, config, logger):
        """初始化绘图器|Initialize plotter

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def plot_single_population(self) -> bool:
        """绘制单个群体的LD衰减图|Plot LD decay figure for single population

        Returns:
            bool: 是否成功|Success
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始绘制LD衰减图|Start plotting LD decay figure")
            self.logger.info("=" * 60)

            # 构建命令|Build command
            cmd = self._build_plot_command()

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
                self.logger.error(f"绘图失败|Plotting failed")
                self.logger.error(f"错误信息|Error: {result.stderr}")
                return False

            # 输出结果|Output result
            if result.stdout:
                self.logger.info(f"输出|Output:\n{result.stdout}")

            # 检查输出文件|Check output files
            output_files = self.config.get_output_files()
            fig_file = Path(output_files['figure'])

            if not fig_file.exists():
                self.logger.warning(f"图像文件未生成|Figure file not generated: {fig_file}")
                # 不返回False，因为统计文件已经生成|Don't return False as stat file is generated
            else:
                self.logger.info(f"图像文件已保存|Figure file saved: {fig_file}")

            self.logger.info("=" * 60)
            self.logger.info("LD衰减图绘制完成|LD decay figure plotting completed")
            self.logger.info("=" * 60)

            return True

        except Exception as e:
            self.logger.error(f"绘制LD衰减图时出错|Error plotting LD decay figure: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _build_plot_command(self) -> list:
        """构建绘图命令|Build plotting command

        Returns:
            list: 命令列表|Command list
        """
        cmd = ['perl', self.config.plot_one_pop_path]

        # 输入文件|Input file
        output_files = self.config.get_output_files()
        stat_file = output_files['stat']
        cmd.extend(["-inFile", stat_file])

        # 输出文件|Output file
        output_prefix = str(self.config.output_path)
        cmd.extend(["-output", output_prefix])

        # 绘图参数|Plotting parameters
        cmd.extend(["-bin1", str(self.config.bin1)])
        cmd.extend(["-bin2", str(self.config.bin2)])
        cmd.extend(["-break", str(self.config.break_point)])
        cmd.extend(["-measure", self.config.measure])
        cmd.extend(["-method", self.config.method])
        cmd.extend(["-percent", str(self.config.percentile)])

        # 最大X坐标|Max X coordinate
        if self.config.max_x is not None:
            cmd.extend(["-maxX", str(self.config.max_x)])

        return cmd

    def plot_multiple_populations(self, pop_list: List[tuple]) -> bool:
        """绘制多个群体的LD衰减图|Plot LD decay figure for multiple populations

        Args:
            pop_list: 群体列表|Population list [(stat_file, pop_id), ...]

        Returns:
            bool: 是否成功|Success
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始绘制多群体LD衰减图|Start plotting multi-population LD decay figure")
            self.logger.info("=" * 60)

            # 创建临时列表文件|Create temporary list file
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.list', delete=False) as f:
                for stat_file, pop_id in pop_list:
                    f.write(f"{stat_file}\t{pop_id}\n")
                list_file = f.name

            # 构建命令|Build command
            output_prefix = str(self.config.output_path)
            cmd = [
                'perl', self.config.plot_multi_pop_path,
                '-inList', list_file,
                '-output', output_prefix
            ]

            # 打印命令|Print command
            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            # 运行命令|Run command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            # 删除临时文件|Delete temporary file
            Path(list_file).unlink()

            # 检查结果|Check result
            if result.returncode != 0:
                self.logger.error(f"绘图失败|Plotting failed")
                self.logger.error(f"错误信息|Error: {result.stderr}")
                return False

            # 输出结果|Output result
            if result.stdout:
                self.logger.info(f"输出|Output:\n{result.stdout}")

            self.logger.info("=" * 60)
            self.logger.info("多群体LD衰减图绘制完成|Multi-population LD decay figure plotting completed")
            self.logger.info("=" * 60)

            return True

        except Exception as e:
            self.logger.error(f"绘制多群体LD衰减图时出错|Error plotting multi-population LD decay figure: {e}")
            import traceback
            traceback.print_exc()
            return False
