"""
Vcftools Pi计算模块|Vcftools Pi Calculation Module
使用vcftools计算群体内核苷酸多样性(pi)
Calculate within-population nucleotide diversity (pi) using vcftools
"""

import os
from pathlib import Path
from typing import Dict, List

from .config import PiConfig
from .utils import (
    PiRow, CommandRunner, build_conda_command,
    parse_population_file, is_step_completed, format_number
)


class VcftoolsPiCalculator:
    """Vcftools Pi计算器|Vcftools Pi Calculator"""

    def __init__(self, config: PiConfig, logger, chrom_lengths: Dict[str, int] = None):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger)
        self.chrom_lengths = chrom_lengths or {}
        self.results: List[PiRow] = []

    def calculate(self) -> List[PiRow]:
        """
        执行vcftools pi计算|Execute vcftools pi calculation

        Returns:
            Pi结果列表|List of Pi results
        """
        self.logger.info("开始vcftools pi计算|Starting vcftools pi calculation")

        # 解析群体文件|Parse population file
        pop_to_samples = parse_population_file(self.config.pop_file)
        populations = sorted(pop_to_samples.keys())

        self.logger.info(
            f"检测到{len(populations)}个群体|Detected {len(populations)} populations: "
            f"{', '.join(populations)}"
        )

        # 判断模式|Determine mode
        is_windowed = self.config.window_size is not None
        mode_name = "窗口模式" if is_windowed else "全基因组模式"
        self.logger.info(f"计算模式|Calculation mode: {mode_name}")

        # 对每个群体分别计算|Calculate for each population separately
        for pop_name in populations:
            samples = pop_to_samples[pop_name]
            self.logger.info(
                f"计算群体|Calculating population: {pop_name} ({len(samples)}个样本|samples)"
            )

            if is_windowed:
                self._calculate_windowed(pop_name, samples)
            else:
                self._calculate_genome_wide(pop_name, samples)

        self.logger.info(
            f"vcftools pi计算完成，共{len(self.results)}条结果|"
            f"vcftools pi calculation completed, {len(self.results)} results"
        )
        return self.results

    def _calculate_windowed(self, pop_name: str, samples: List[str]):
        """
        窗口模式pi计算|Windowed pi calculation

        Args:
            pop_name: 群体名称|Population name
            samples: 样本列表|Sample list
        """
        output_prefix = str(self.config.output_path / '01_vcftools' / pop_name)
        keep_file = f"{output_prefix}.keep.txt"

        # 断点续传：检查输出文件是否已存在
        # Checkpoint resume: check if output file already exists
        output_file = f"{output_prefix}.windowed.pi"
        if is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: {pop_name} windowed pi")
            self._parse_windowed_output(output_file, pop_name)
            return

        # 写入keep文件|Write keep file
        with open(keep_file, 'w') as f:
            for sample in samples:
                f.write(f"{sample}\n")

        # 构建vcftools命令|Build vcftools command
        args = [
            '--gzvcf', self.config.vcf_file,
            '--keep', keep_file,
            '--window-pi', str(self.config.window_size),
        ]

        # 添加步长（仅步长!=窗口大小时需要指定）|Add step (only when step != window_size)
        if self.config.window_step is not None and self.config.window_step != self.config.window_size:
            args.extend(['--window-pi-step', str(self.config.window_step)])

        args.extend(['--out', output_prefix])

        # 添加质控参数|Add QC parameters
        if self.config.maf > 0:
            args.extend(['--maf', str(self.config.maf)])
        if self.config.max_missing < 1.0:
            args.extend(['--max-missing', str(self.config.max_missing)])

        cmd = build_conda_command(self.config.vcftools_path, args)
        success = self.runner.run(cmd, f"vcftools窗口pi计算|vcftools windowed pi for {pop_name}")

        if success and os.path.exists(output_file):
            self._parse_windowed_output(output_file, pop_name)
        else:
            self.logger.warning(
                f"群体{pop_name}的vcftools窗口pi计算失败|"
                f"vcftools windowed pi calculation failed for {pop_name}"
            )

        # 清理keep文件（如果不需要保留中间文件）|Clean up keep file if not keeping intermediates
        if not self.config.keep_intermediate and os.path.exists(keep_file):
            os.remove(keep_file)

    def _calculate_genome_wide(self, pop_name: str, samples: List[str]):
        """
        全基因组模式pi计算|Genome-wide pi calculation

        Args:
            pop_name: 群体名称|Population name
            samples: 样本列表|Sample list
        """
        output_prefix = str(self.config.output_path / '01_vcftools' / pop_name)
        keep_file = f"{output_prefix}.keep.txt"

        # 断点续传|Checkpoint resume
        output_file = f"{output_prefix}.sites.pi"
        if is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: {pop_name} site pi")
            self._parse_site_pi_output(output_file, pop_name)
            return

        # 写入keep文件|Write keep file
        with open(keep_file, 'w') as f:
            for sample in samples:
                f.write(f"{sample}\n")

        # 构建vcftools命令|Build vcftools command
        args = [
            '--gzvcf', self.config.vcf_file,
            '--keep', keep_file,
            '--site-pi',
            '--out', output_prefix
        ]

        # 添加质控参数|Add QC parameters
        if self.config.maf > 0:
            args.extend(['--maf', str(self.config.maf)])
        if self.config.max_missing < 1.0:
            args.extend(['--max-missing', str(self.config.max_missing)])

        cmd = build_conda_command(self.config.vcftools_path, args)
        success = self.runner.run(cmd, f"vcftools site pi计算|vcftools site pi for {pop_name}")

        if success and os.path.exists(output_file):
            self._parse_site_pi_output(output_file, pop_name)
        else:
            self.logger.warning(
                f"群体{pop_name}的vcftools site pi计算失败|"
                f"vcftools site pi calculation failed for {pop_name}"
            )

        # 清理keep文件|Clean up keep file
        if not self.config.keep_intermediate and os.path.exists(keep_file):
            os.remove(keep_file)

    def _parse_windowed_output(self, output_file: str, pop_name: str):
        """
        解析vcftools窗口pi输出|Parse vcftools windowed pi output

        vcftools .windowed.pi格式:
        CHROM    BIN_START    BIN_END    N_VARIANTS    N_MONOMORPHIC    PI
        """
        n_rows = 0
        with open(output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('CHROM'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 6:
                    try:
                        pi_val = float(parts[5])
                        row = PiRow(
                            population=pop_name,
                            chromosome=parts[0],
                            window_start=int(parts[1]),
                            window_end=int(parts[2]),
                            pi_value=pi_val,
                            n_sites=int(parts[3]),
                            source="vcftools"
                        )
                        self.results.append(row)
                        n_rows += 1
                    except (ValueError, IndexError):
                        continue

        self.logger.info(
            f"  {pop_name}: 解析到{n_rows}个窗口的pi值|"
            f"{pop_name}: parsed pi values for {n_rows} windows"
        )

    def _parse_site_pi_output(self, output_file: str, pop_name: str):
        """
        解析vcftools site-pi输出，计算每条染色体和全基因组的pi
        Parse vcftools site-pi output, calculate per-chromosome and genome-wide pi

        pi计算公式: pi = sum(位点pi之和) / 染色体(基因组)碱基总数
        Pi formula: pi = sum(site_pi) / total_bases_of_chromosome(or_genome)

        Args:
            output_file: vcftools .sites.pi输出文件路径|vcftools .sites.pi output file path
            pop_name: 群体名称|Population name
        """
        # 收集每条染色体的pi值之和|Collect sum of pi values per chromosome
        chrom_pi_sum: Dict[str, float] = {}
        chrom_pi_count: Dict[str, int] = {}

        with open(output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('CHROM'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        chrom = parts[0]
                        pi_val = float(parts[2])
                        # 跳过nan值（单态位点/全部缺失）|Skip nan values (monomorphic/missing sites)
                        if pi_val != pi_val:  # nan != nan is True
                            continue
                        if chrom not in chrom_pi_sum:
                            chrom_pi_sum[chrom] = 0.0
                            chrom_pi_count[chrom] = 0
                        chrom_pi_sum[chrom] += pi_val
                        chrom_pi_count[chrom] += 1
                    except (ValueError, IndexError):
                        continue

        # 计算每条染色体的pi = sum(pi) / 染色体碱基数
        # Per-chromosome pi = sum(pi) / chromosome length in bases
        total_pi_sum = 0.0
        total_pi_sites = 0
        total_bases = 0

        for chrom in sorted(chrom_pi_sum.keys()):
            pi_sum = chrom_pi_sum[chrom]
            n_sites = chrom_pi_count[chrom]

            # 获取染色体长度，使用fai文件中的值|Get chromosome length from fai file
            chrom_len = self.chrom_lengths.get(chrom)
            if chrom_len is None:
                self.logger.warning(
                    f"  {chrom}: 在fai文件中未找到染色体长度，跳过|"
                    f"{chrom}: chromosome length not found in fai file, skipping"
                )
                continue

            # pi = sum(所有位点的pi) / 染色体总碱基数
            # pi = sum(pi at all sites) / total bases of chromosome
            pi = pi_sum / chrom_len

            row = PiRow(
                population=pop_name,
                chromosome=chrom,
                window_start=None,
                window_end=None,
                pi_value=pi,
                n_sites=n_sites,
                source="vcftools"
            )
            self.results.append(row)
            total_pi_sum += pi_sum
            total_pi_sites += n_sites
            total_bases += chrom_len

            self.logger.info(
                f"  {pop_name} {chrom}: pi={format_number(pi)} "
                f"({n_sites}个位点|sites, 染色体长度|chrom_length={chrom_len})"
            )

        # 全基因组pi = sum(所有位点pi) / 全基因组总碱基数
        # Genome-wide pi = sum(all site pi) / total genome bases
        if total_bases > 0:
            genome_wide_pi = total_pi_sum / total_bases
            row = PiRow(
                population=pop_name,
                chromosome="genome_wide",
                window_start=None,
                window_end=None,
                pi_value=genome_wide_pi,
                n_sites=total_pi_sites,
                source="vcftools"
            )
            self.results.append(row)

            self.logger.info(
                f"  {pop_name} genome_wide: pi={format_number(genome_wide_pi)} "
                f"({total_pi_sites}个位点|sites, "
                f"基因组大小|genome_size={total_bases})"
            )
