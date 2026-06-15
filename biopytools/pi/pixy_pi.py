"""
Pixy Pi计算模块|Pixy Pi Calculation Module
使用pixy计算群体内核苷酸多样性(pi)
Calculate within-population nucleotide diversity (pi) using pixy
"""

import os
import subprocess
from pathlib import Path
from typing import List

from .config import PiConfig
from .utils import (
    PiRow, CommandRunner, build_conda_command,
    parse_population_file, is_step_completed, format_number
)


class PixyPiCalculator:
    """Pixy Pi计算器|Pixy Pi Calculator"""

    def __init__(self, config: PiConfig, logger):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger)
        self.results: List[PiRow] = []
        self._bypass_invariant = False

    def calculate(self) -> List[PiRow]:
        """
        执行pixy pi计算|Execute pixy pi calculation

        Returns:
            Pi结果列表|List of Pi results
        """
        self.logger.info("开始pixy pi计算|Starting pixy pi calculation")

        # pixy不支持全基因组模式|pixy doesn't support genome-wide mode
        if self.config.window_size is None:
            self.logger.warning(
                "pixy不支持全基因组模式，跳过pixy计算|"
                "pixy does not support genome-wide mode, skipping pixy calculation"
            )
            return self.results

        # 解析群体文件|Parse population file
        pop_to_samples = parse_population_file(self.config.pop_file)
        populations = sorted(pop_to_samples.keys())

        self.logger.info(
            f"检测到{len(populations)}个群体|Detected {len(populations)} populations: "
            f"{', '.join(populations)}"
        )
        self.logger.info(
            f"计算模式|Calculation mode: 窗口模式 (window_size={self.config.window_size}bp)"
        )

        # pixy不支持步长参数（步长始终等于窗口大小）|pixy doesn't support step (step always equals window_size)
        if self.config.window_step is not None and self.config.window_step != self.config.window_size:
            effective_step = self.config.window_size
            self.logger.warning(
                f"pixy不支持自定义步长，将使用步长={effective_step}bp（等于窗口大小）|"
                f"pixy does not support custom step, using step={effective_step}bp (equals window_size)"
            )

        # 自动检测VCF不变位点|Auto-detect VCF invariant sites
        self._check_invariant_sites()

        # 构建并执行pixy命令|Build and execute pixy command
        output_dir = str(self.config.output_path / '02_pixy')
        output_file = os.path.join(output_dir, 'pi_output.txt')

        # 断点续传|Checkpoint resume
        if is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: pixy pi calculation")
            self._parse_pixy_output(output_file)
            return self.results

        # 构建pixy命令|Build pixy command
        args = [
            '--stats', 'pi',
            '--vcf', self.config.vcf_file,
            '--populations', self.config.pop_file,
            '--output_folder', output_dir,
            '--window_size', str(self.config.window_size),
            '--n_cores', str(self.config.threads)
        ]

        # 如果需要绕过不变位点检查|If need to bypass invariant sites check
        if self._bypass_invariant:
            args.append('--bypass_invariant_check')
            self.logger.info("将绕过不变位点检查|Will bypass invariant sites check")

        cmd = build_conda_command(self.config.pixy_env, args)
        success = self.runner.run(cmd, "pixy pi计算|pixy pi calculation", timeout=7200)

        if success:
            # pixy输出文件名: pi_output.txt 在 output_folder 中
            # pixy output file: pi_output.txt in output_folder
            if os.path.exists(output_file):
                self._parse_pixy_output(output_file)
            else:
                self.logger.warning(
                    f"pixy输出文件未找到|pixy output file not found: {output_file}"
                )
        else:
            self.logger.error("pixy pi计算失败|pixy pi calculation failed")

        return self.results

    def _check_invariant_sites(self):
        """
        检查VCF是否包含不变位点，决定是否绕过检查
        Check if VCF contains invariant sites, decide whether to bypass check

        使用bcftools检查前100个位点是否有ALT="."的位点
        Use bcftools to check if first 100 sites have ALT="."
        """
        self.logger.info("检查VCF是否包含不变位点|Checking if VCF contains invariant sites")

        try:
            # 使用pixy所在conda环境中的bcftools|Use bcftools from pixy's conda environment
            check_cmd = build_conda_command(self.config.pixy_env, [
                'bcftools', 'view', '-H', self.config.vcf_file
            ])
            # 通过管道传递给head和grep|Pipe through head and grep
            check_pipeline = (
                f"{' '.join(check_cmd)} | head -n 100 | grep -c '\\t\\.\\t'"
            )

            result = subprocess.run(
                check_pipeline,
                shell=True,
                capture_output=True,
                text=True,
                timeout=120
            )

            # grep返回0表示找到了不变位点|grep returns 0 if invariant sites found
            if result.returncode == 0:
                count = int(result.stdout.strip())
                self.logger.info(
                    f"  VCF包含不变位点（前100个位点中发现{count}个）|"
                    f"VCF contains invariant sites (found {count} in first 100 sites)"
                )
                self._bypass_invariant = False
            else:
                self.logger.warning(
                    "  VCF不包含不变位点（ALT=\".\"），将自动绕过检查|"
                    "VCF does not contain invariant sites (ALT=\".\"), will bypass check automatically"
                )
                self.logger.warning(
                    "  注意：这可能导致pi统计结果不准确|"
                    "Note: This may cause inaccurate pi statistics"
                )
                self._bypass_invariant = True

        except subprocess.TimeoutExpired:
            self.logger.warning(
                "  VCF不变位点检查超时，假设不包含不变位点|"
                "VCF invariant sites check timeout, assuming no invariant sites"
            )
            self._bypass_invariant = True
        except Exception as e:
            self.logger.warning(
                f"  VCF不变位点检查失败: {e}，假设不包含不变位点|"
                f"VCF invariant sites check failed: {e}, assuming no invariant sites"
            )
            self._bypass_invariant = True

    def _parse_pixy_output(self, output_file: str):
        """
        解析pixy pi输出文件|Parse pixy pi output file

        pixy pi_output.txt格式 (tab-separated):
        pop   chrom   window_pos_1   window_pos_2   avg_pi   n_sites   ...
        """
        n_rows = 0
        with open(output_file, 'r') as f:
            # 读取header确定列位置|Read header to determine column positions
            header = f.readline().strip()
            headers = header.split('\t')

            # 找到关键列的索引|Find indices of key columns
            try:
                idx_pop = headers.index('pop')
                idx_chrom = headers.index('chrom')
                idx_window_start = headers.index('window_pos_1')
                idx_window_end = headers.index('window_pos_2')
                idx_avg_pi = headers.index('avg_pi')
                idx_n_sites = headers.index('n_sites')
            except ValueError as e:
                self.logger.error(f"pixy输出文件格式错误|pixy output file format error: {e}")
                return

            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) <= max(idx_pop, idx_chrom, idx_window_start,
                                     idx_window_end, idx_avg_pi, idx_n_sites):
                    continue

                try:
                    pi_val = parts[idx_avg_pi]
                    if pi_val.upper() == 'NA' or pi_val == '':
                        continue
                    pi_float = float(pi_val)

                    row = PiRow(
                        population=parts[idx_pop],
                        chromosome=parts[idx_chrom],
                        window_start=int(parts[idx_window_start]),
                        window_end=int(parts[idx_window_end]),
                        pi_value=pi_float,
                        n_sites=int(parts[idx_n_sites]) if parts[idx_n_sites].upper() != 'NA' else 0,
                        source="pixy"
                    )
                    self.results.append(row)
                    n_rows += 1
                except (ValueError, IndexError):
                    continue

        self.logger.info(
            f"pixy pi计算完成，解析到{n_rows}条结果|"
            f"pixy pi calculation completed, parsed {n_rows} results"
        )
