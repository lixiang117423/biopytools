"""
Pixy计算器模块|Pixy Calculator Module
"""

import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
import time


class PixyCalculator:
    """Pixy计算器类|Pixy Calculator Class"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: Pixy配置对象|Pixy configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def _build_pixy_command(self, stats: str) -> List[str]:
        """构建pixy命令|Build pixy command

        Args:
            stats: 统计量类型（pi/dxy/fst）|Statistic type (pi/dxy/fst)

        Returns:
            List[str]: pixy命令列表|pixy command list
        """
        # 从conda环境路径提取环境名|Extract environment name from conda environment path
        env_name = Path(self.config.conda_env).name

        cmd = [
            "conda", "run", "-n", env_name, "--no-capture-output",
            "pixy",
            "--stats", stats,
            "--vcf", str(self.config.vcf_path),
            "--populations", str(self.config.pop_path),
            "--output_folder", str(self.config.output_path),
            "--n_cores", str(self.config.threads)
        ]

        # 添加窗口参数|Add window parameters
        if self.config.window_size:
            cmd.extend(["--window_size", str(self.config.window_size)])

        if self.config.bed_file:
            cmd.extend(["--bed_file", str(self.config.bed_path)])

        if self.config.sites_file:
            cmd.extend(["--sites_file", str(self.config.sites_path)])

        # 添加质控参数|Add QC parameters
        if self.config.min_samples > 0:
            cmd.extend(["--min_samples", str(self.config.min_samples)])

        if self.config.max_missing < 1.0:
            cmd.extend(["--max_missing", str(self.config.max_missing)])

        if self.config.min_maf > 0.0:
            cmd.extend(["--min_maf", str(self.config.min_maf)])

        if self.config.zscore_window:
            cmd.extend(["--zscore_window", str(self.config.zscore_window)])

        # 添加染色体参数|Add chromosome parameters
        if self.config.chromosome_list:
            cmd.extend(["--chromosomes", ",".join(self.config.chromosome_list)])

        # 添加bypass_invariant_check参数（不需要值）|Add bypass_invariant_check parameter (no value needed)
        if self.config.bypass_invariant_check:
            cmd.append("--bypass_invariant_check")

        return cmd

    def _run_command(self, cmd: List[str], description: str) -> bool:
        """运行命令|Run command

        Args:
            cmd: 命令列表|Command list
            description: 命令描述|Command description

        Returns:
            bool: 是否成功|Whether successful
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 记录完整命令（INFO级别）|Log complete command (INFO level)
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding='utf-8',
                errors='replace'
            )

            elapsed_time = time.time() - start_time

            if result.returncode == 0:
                self.logger.info(f"{description} 完成|completed successfully (耗时|time: {elapsed_time:.2f}秒|s)")

                # 输出pixy的日志|Output pixy logs
                if result.stdout:
                    self.logger.debug(f"Pixy输出|Pixy stdout:\n{result.stdout}")

                return True
            else:
                self.logger.error(f"{description} 失败|failed")
                self.logger.error(f"返回码|Return code: {result.returncode}")

                # 输出错误信息|Output error messages
                if result.stderr:
                    self.logger.error(f"错误信息|Error:\n{result.stderr}")

                if result.stdout:
                    self.logger.error(f"输出信息|Output:\n{result.stdout}")

                return False

        except subprocess.TimeoutExpired:
            self.logger.error(f"{description} 超时|timeout")
            return False
        except Exception as e:
            self.logger.error(f"{description} 异常|error: {e}")
            return False

    def calculate_pi(self) -> Dict:
        """计算pi（核苷酸多样性）|Calculate pi (nucleotide diversity)

        Returns:
            dict: 计算结果|Calculation results
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算pi|Calculating pi (nucleotide diversity)")
        self.logger.info("=" * 60)

        cmd = self._build_pixy_command("pi")
        success = self._run_command(cmd, "pi计算|pi calculation")

        # 查找输出文件|Find output files
        output_file = None
        if success:
            # pixy输出文件命名格式: [prefix]pi.txt
            potential_files = list(self.config.output_path.glob("*pi*"))
            if potential_files:
                output_file = potential_files[0]
                self.logger.info(f"输出文件|Output file: {output_file}")

        return {
            'success': success,
            'statistic': 'pi',
            'output_file': str(output_file) if output_file else None
        }

    def calculate_dxy(self) -> Dict:
        """计算dxy（群体间核苷酸差异）|Calculate dxy (nucleotide divergence)

        Returns:
            dict: 计算结果|Calculation results
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算dxy|Calculating dxy (nucleotide divergence)")
        self.logger.info("=" * 60)

        cmd = self._build_pixy_command("dxy")
        success = self._run_command(cmd, "dxy计算|dxy calculation")

        # 查找输出文件|Find output files
        output_file = None
        if success:
            potential_files = list(self.config.output_path.glob("*dxy*"))
            if potential_files:
                output_file = potential_files[0]
                self.logger.info(f"输出文件|Output file: {output_file}")

        return {
            'success': success,
            'statistic': 'dxy',
            'output_file': str(output_file) if output_file else None
        }

    def calculate_fst(self) -> Dict:
        """计算fst（遗传分化系数）|Calculate fst (genetic differentiation)

        Returns:
            dict: 计算结果|Calculation results
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算fst|Calculating fst (genetic differentiation)")
        self.logger.info("=" * 60)

        cmd = self._build_pixy_command("fst")
        success = self._run_command(cmd, "fst计算|fst calculation")

        # 查找输出文件|Find output files
        output_file = None
        if success:
            potential_files = list(self.config.output_path.glob("*fst*"))
            if potential_files:
                output_file = potential_files[0]
                self.logger.info(f"输出文件|Output file: {output_file}")

        return {
            'success': success,
            'statistic': 'fst',
            'output_file': str(output_file) if output_file else None
        }

    def calculate_watterson_theta(self) -> Dict:
        """计算Watterson's theta|Calculate Watterson's theta

        Returns:
            dict: 计算结果|Calculation results
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算Watterson's theta|Calculating Watterson's theta")
        self.logger.info("=" * 60)

        cmd = self._build_pixy_command("watterson_theta")
        success = self._run_command(cmd, "Watterson's theta计算|Watterson's theta calculation")

        # 查找输出文件|Find output files
        output_file = None
        if success:
            potential_files = list(self.config.output_path.glob("*Watterson*"))
            if potential_files:
                output_file = potential_files[0]
                self.logger.info(f"输出文件|Output file: {output_file}")

        return {
            'success': success,
            'statistic': 'watterson_theta',
            'output_file': str(output_file) if output_file else None
        }

    def calculate_tajima_d(self) -> Dict:
        """计算Tajima's D|Calculate Tajima's D

        Returns:
            dict: 计算结果|Calculation results
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算Tajima's D|Calculating Tajima's D")
        self.logger.info("=" * 60)

        cmd = self._build_pixy_command("tajima_d")
        success = self._run_command(cmd, "Tajima's D计算|Tajima's D calculation")

        # 查找输出文件|Find output files
        output_file = None
        if success:
            potential_files = list(self.config.output_path.glob("*TajimaD*"))
            if potential_files:
                output_file = potential_files[0]
                self.logger.info(f"输出文件|Output file: {output_file}")

        return {
            'success': success,
            'statistic': 'tajima_d',
            'output_file': str(output_file) if output_file else None
        }

    def run_all_calculations(self) -> Dict[str, Dict]:
        """运行所有选中的统计量计算|Run all selected statistics calculations

        Returns:
            dict: 所有计算结果|All calculation results
        """
        results = {}
        overall_success = True

        # pi计算|pi calculation
        if self.config.calc_pi:
            result = self.calculate_pi()
            results['pi'] = result
            if not result['success']:
                overall_success = False

        # dxy计算|dxy calculation
        if self.config.calc_dxy:
            result = self.calculate_dxy()
            results['dxy'] = result
            if not result['success']:
                overall_success = False

        # fst计算|fst calculation
        if self.config.calc_fst:
            result = self.calculate_fst()
            results['fst'] = result
            if not result['success']:
                overall_success = False

        # Watterson's theta计算|Watterson's theta calculation
        if self.config.calc_watterson_theta:
            result = self.calculate_watterson_theta()
            results['watterson_theta'] = result
            if not result['success']:
                overall_success = False

        # Tajima's D计算|Tajima's D calculation
        if self.config.calc_tajima_d:
            result = self.calculate_tajima_d()
            results['tajima_d'] = result
            if not result['success']:
                overall_success = False

        results['overall_success'] = overall_success

        return results
