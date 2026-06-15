"""
Fst计算核心逻辑模块|Fst Calculation Core Logic Module
"""

import os
import random
from pathlib import Path
from typing import Dict, List, Tuple
from .utils import (
    CommandRunner, build_conda_command, create_plink_population_file,
    parse_population_file, analyze_populations, decide_bootstrap_strategy,
    decide_thinning_strategy, filter_population_file
)


class FstCalculator:
    """Fst计算器|Fst Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger, config.output_path)

        # 获取样本到群体的映射|Get sample to population mapping
        self.sample_to_pop = parse_population_file(config.pop_file)

        # 分析群体信息，过滤小样本群体|Analyze populations, filter small ones
        self.pop_stats = analyze_populations(
            self.sample_to_pop,
            config.min_samples,
            config.exclude_pops
        )

        self.filtered_populations = sorted(self.pop_stats['pop_counts'].keys())

        # 输出群体统计信息|Output population statistics
        self.logger.info(f"原始群体数|Original populations: {len(set(self.sample_to_pop.values()))}")
        self.logger.info(f"过滤后群体数|Filtered populations: {len(self.filtered_populations)}")
        self.logger.info(f"群体统计|Population statistics:")
        for pop, count in self.pop_stats['pop_counts'].items():
            self.logger.info(f"  {pop}: {count} 个样本|samples")

        # 决定bootstrap策略|Decide bootstrap strategy
        self.should_bootstrap, self.n_iterations = decide_bootstrap_strategy(
            self.pop_stats['ratio'],
            self.config.enable_bootstrap
        )

        self.logger.info(f"样本量比值|Sample size ratio: {self.pop_stats['ratio']:.1f}")
        self.logger.info(f"Bootstrap策略|Bootstrap strategy: {self.should_bootstrap} ({self.n_iterations} 迭代|iterations)")

    def step1_vcf_to_plink(self) -> bool:
        """步骤1: VCF转PLINK格式|Step 1: VCF to PLINK format"""
        self.logger.info("步骤1: VCF转PLINK格式|Step 1: Converting VCF to PLINK format")

        # 设置输出前缀|Set output prefix
        if self.config.keep_intermediate:
            output_prefix = str(self.config.intermediate_dir / 'mydata')
        else:
            output_prefix = str(self.config.output_path / 'mydata')

        # 构建命令|Build command
        args = [
            '--vcf', self.config.vcf_file,
            '--make-bed',
            '--out', output_prefix
        ]

        cmd = build_conda_command(self.config.plink_path, args)

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, "VCF转PLINK格式|VCF to PLINK conversion")

        if success:
            self.plink_prefix = output_prefix
            self.logger.info(f"PLINK文件创建完成|PLINK files created: {output_prefix}")

        return success

    def step2_quality_control(self) -> bool:
        """步骤2: 质量控制|Step 2: Quality control"""
        if not self.config.enable_qc:
            self.logger.info("步骤2: 跳过质量控制（已禁用）|Step 2: Skipping quality control (disabled)")
            # 直接使用步骤1的输出|Use output from step 1 directly
            return True

        self.logger.info("步骤2: 质量控制|Step 2: Quality control")

        # 设置输出前缀|Set output prefix
        if self.config.keep_intermediate:
            output_prefix = str(self.config.intermediate_dir / 'mydata_qc')
        else:
            output_prefix = str(self.config.output_path / 'mydata_qc')

        # 构建质控命令|Build QC command
        args = [
            '--bfile', self.plink_prefix,
            '--maf', str(self.config.maf),
            '--geno', str(self.config.geno),
            '--mind', str(self.config.mind),
            '--hwe', str(self.config.hwe),
            '--make-bed',
            '--out', output_prefix
        ]

        cmd = build_conda_command(self.config.plink_path, args)

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, "质量控制|Quality control")

        if success:
            self.plink_prefix = output_prefix
            self.logger.info(f"质控完成|QC completed: {output_prefix}")

        return success

    def step3_ld_pruning(self) -> bool:
        """步骤3: LD pruning (可选)|Step 3: LD pruning (optional)"""
        if not self.config.enable_ld_prune:
            self.logger.info("跳过LD pruning|Skipping LD pruning")
            self.plink_pruned_prefix = self.plink_prefix
            return True

        self.logger.info("步骤3: LD pruning|Step 3: LD pruning")

        # 设置输出前缀|Set output prefix
        if self.config.keep_intermediate:
            prune_output = str(self.config.intermediate_dir / 'ld_prune')
            pruned_prefix = str(self.config.intermediate_dir / 'mydata_pruned')
        else:
            prune_output = str(self.config.output_path / 'ld_prune')
            pruned_prefix = str(self.config.output_path / 'mydata_pruned')

        # 计算SNP数量|Count SNPs
        snp_count = sum(1 for line in open(f"{self.plink_prefix}.bim"))
        self.logger.info(f"SNP数量|SNP count: {snp_count}")

        # 决定是否需要LD pruning|Decide if LD pruning is needed
        if snp_count < 500000:
            self.logger.info("SNP数量<50万，跳过LD pruning|SNP count < 500k, skipping LD pruning")
            self.plink_pruned_prefix = self.plink_prefix
            return True

        # 执行LD pruning|Perform LD pruning
        args = [
            '--bfile', self.plink_prefix,
            '--indep-pairwise',
            str(self.config.ld_window_size),
            str(self.config.ld_step_size),
            str(self.config.ld_r2_threshold),
            '--out', prune_output
        ]

        cmd = build_conda_command(self.config.plink_path, args)
        success = self.cmd_runner.run(cmd, "LD pruning计算|LD pruning calculation")

        if not success:
            self.logger.warning("LD pruning失败，使用原始数据|LD pruning failed, using original data")
            self.plink_pruned_prefix = self.plink_prefix
            return True

        # 提取pruning后的SNP|Extract pruned SNPs
        args = [
            '--bfile', self.plink_prefix,
            '--extract', f"{prune_output}.prune.in",
            '--make-bed',
            '--out', pruned_prefix
        ]

        cmd = build_conda_command(self.config.plink_path, args)
        success = self.cmd_runner.run(cmd, "提取pruning后SNP|Extract pruned SNPs")

        if success:
            self.plink_pruned_prefix = pruned_prefix
            # 统计pruning后的SNP数量|Count SNPs after pruning
            pruned_snp_count = sum(1 for line in open(f"{pruned_prefix}.bim"))
            self.logger.info(f"LD pruning后SNP数|SNPs after LD pruning: {pruned_snp_count}")

        return success

    def step4_calculate_fst(self, iteration: int = None, keep_samples: List[str] = None) -> bool:
        """
        步骤4: 计算Fst|Step 4: Calculate Fst

        Args:
            iteration: bootstrap迭代编号|Bootstrap iteration number
            keep_samples: 要保留的样本列表|List of samples to keep
        """
        if iteration is not None:
            self.logger.info(f"步骤4: 计算Fst (迭代|iteration {iteration})|Step 4: Calculating Fst (iteration {iteration})")
        else:
            self.logger.info("步骤4: 计算Fst|Step 4: Calculating Fst statistics")

        # 获取群体对列表|Get population pairs list
        pop_pairs = self._get_population_pairs()

        # 如果只有两个群体，使用原来的逻辑|If only two populations, use original logic
        if len(self.filtered_populations) == 2:
            return self._calculate_fst_single(iteration, keep_samples)
        else:
            # 多群体：为每对群体分别计算FST|Multiple populations: calculate FST for each pair separately
            return self._calculate_fst_pairwise(iteration, keep_samples, pop_pairs)

    def _calculate_fst_single(self, iteration: int = None, keep_samples: List[str] = None) -> bool:
        """
        计算单个FST值（两个群体或整体）|Calculate single FST value (two populations or overall)

        Args:
            iteration: bootstrap迭代编号|Bootstrap iteration number
            keep_samples: 要保留的样本列表|List of samples to keep

        Returns:
            bool: 是否成功|Whether successful
        """
        # 创建过滤后的群体文件|Create filtered population file
        if iteration is not None:
            # Bootstrap模式：文件放到bootstrap_iterations文件夹|Bootstrap mode: files in bootstrap_iterations folder
            bootstrap_dir = self.config.output_path / 'bootstrap_iterations'
            bootstrap_dir.mkdir(exist_ok=True)

            plink_pop_file = bootstrap_dir / f'population_iter{iteration}.plink'
            filtered_sample_to_pop = {s: self.sample_to_pop[s] for s in keep_samples}
            create_plink_population_file(filtered_sample_to_pop, str(plink_pop_file))

            # 设置输出文件|Set output file
            fst_output = bootstrap_dir / f'fst_iter{iteration}'
        else:
            # 直接计算模式|Direct calculation mode
            # 过滤群体文件|Filter population file
            filtered_pop_file = self.config.output_path / 'population_filtered.plink'
            filter_population_file(
                self.config.pop_file,
                str(filtered_pop_file),
                self.filtered_populations
            )
            plink_pop_file = filtered_pop_file

            # 设置输出文件|Set output file
            fst_output = self.config.output_path / 'fst_result'

        # 构建Fst计算命令|Build Fst calculation command
        args = [
            '--bfile', self.plink_pruned_prefix,
            '--fst',
            '--within', str(plink_pop_file),
            '--out', str(fst_output)
        ]

        cmd = build_conda_command(self.config.plink_path, args)

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, f"Fst计算|Fst calculation (iter {iteration})" if iteration else "Fst计算|Fst calculation")

        return success

    def _calculate_fst_pairwise(self, iteration: int = None, keep_samples: List[str] = None,
                                 pop_pairs: List[Tuple[str, str]] = None) -> bool:
        """
        为每对群体分别计算FST|Calculate FST for each population pair separately

        Args:
            iteration: bootstrap迭代编号|Bootstrap iteration number
            keep_samples: 要保留的样本列表|List of samples to keep
            pop_pairs: 群体对列表|List of population pairs

        Returns:
            bool: 是否成功|Whether successful
        """
        if not pop_pairs:
            pop_pairs = self._get_population_pairs()

        # Bootstrap模式：文件放到bootstrap_iterations文件夹|Bootstrap mode: files in bootstrap_iterations folder
        if iteration is not None:
            bootstrap_dir = self.config.output_path / 'bootstrap_iterations'
            bootstrap_dir.mkdir(exist_ok=True)

            # 为每对群体创建群体文件并计算FST|Create population file and calculate FST for each pair
            all_success = True
            for pop1, pop2 in pop_pairs:
                # 过滤出这两个群体的样本|Filter samples for these two populations
                if iteration is not None:
                    filtered_samples = [s for s in keep_samples
                                       if self.sample_to_pop[s] in [pop1, pop2]]
                    filtered_sample_to_pop = {s: self.sample_to_pop[s] for s in filtered_samples}
                else:
                    filtered_sample_to_pop = {s: p for s, p in self.sample_to_pop.items()
                                             if p in [pop1, pop2]}

                # 创建临时群体文件|Create temporary population file
                pair_name = f"{pop1}_{pop2}"
                plink_pop_file = bootstrap_dir / f'population_iter{iteration}_{pair_name}.plink'
                create_plink_population_file(filtered_sample_to_pop, str(plink_pop_file))

                # 设置输出文件|Set output file
                fst_output = bootstrap_dir / f'fst_iter{iteration}_{pair_name}'

                # 构建Fst计算命令|Build Fst calculation command
                args = [
                    '--bfile', self.plink_pruned_prefix,
                    '--fst',
                    '--within', str(plink_pop_file),
                    '--out', str(fst_output)
                ]

                cmd = build_conda_command(self.config.plink_path, args)

                # 执行命令|Execute command
                success = self.cmd_runner.run(cmd, f"Fst计算|Fst calculation: {pop1} vs {pop2} (iter {iteration})")
                if not success:
                    all_success = False

            return all_success
        else:
            # 非bootstrap模式，直接使用原来的逻辑|Non-bootstrap mode, use original logic
            return self._calculate_fst_single(iteration, keep_samples)

    def bootstrap_sampling(self, iteration: int) -> List[str]:
        """
        Bootstrap抽样|Bootstrap sampling

        Args:
            iteration: 迭代编号|Iteration number

        Returns:
            抽取的样本ID列表|List of sampled sample IDs
        """
        # 为每个群体抽取min_n个样本|Sample min_n samples from each population
        sampled_samples = []
        min_n = self.pop_stats['min_n']

        for pop in self.filtered_populations:
            # 获取该群体的所有样本|Get all samples from this population
            pop_samples = [s for s, p in self.sample_to_pop.items() if p == pop]
            # 随机抽取min_n个|Randomly sample min_n
            sampled = random.sample(pop_samples, min(min_n, len(pop_samples)))
            sampled_samples.extend(sampled)

        return sampled_samples

    def run_bootstrap(self) -> Dict[str, str]:
        """
        运行bootstrap抽样Fst计算
        Run bootstrap sampling Fst calculation

        Returns:
            输出文件字典|Output files dictionary
        """
        self.logger.info(f"开始Bootstrap Fst计算 ({self.n_iterations} 次迭代)|Starting Bootstrap Fst calculation ({self.n_iterations} iterations)")

        all_results_file = self.config.output_path / 'all_fst_results.txt'
        pairwise_results = {}  # 存储两两群体的Fst值|Store pairwise population Fst values

        for i in range(1, self.n_iterations + 1):
            # Bootstrap抽样|Bootstrap sampling
            keep_samples = self.bootstrap_sampling(i)

            # 计算Fst|Calculate Fst
            if not self.step4_calculate_fst(i, keep_samples):
                self.logger.error(f"迭代{i}失败|Iteration {i} failed")
                continue

            # 提取并保存结果|Extract and save results
            # 检查是否为多群体模式|Check if multiple populations mode
            if len(self.filtered_populations) > 2:
                # 多群体：从每对群体的文件中提取FST值|Multiple populations: extract FST from each pair's file
                pop_pairs = self._get_population_pairs()
                for pop1, pop2 in pop_pairs:
                    pair_name = f"{pop1}_{pop2}"
                    fst_file = self.config.output_path / 'bootstrap_iterations' / f'fst_iter{i}_{pair_name}.fst'
                    if fst_file.exists():
                        iteration_results = self._extract_fst_from_file_pair(fst_file, i, (pop1, pop2))

                        # 汇总结果|Aggregate results
                        for pair, fst_value in iteration_results.items():
                            if pair not in pairwise_results:
                                pairwise_results[pair] = []
                            pairwise_results[pair].append(fst_value)
            else:
                # 两个群体：从单个文件中提取FST值|Two populations: extract FST from single file
                fst_file = self.config.output_path / 'bootstrap_iterations' / f'fst_iter{i}.fst'
                iteration_results = self._extract_fst_from_file(fst_file, i)

                # 汇总结果|Aggregate results
                for pair, fst_value in iteration_results.items():
                    if pair not in pairwise_results:
                        pairwise_results[pair] = []
                    pairwise_results[pair].append(fst_value)

        # 计算统计量并保存|Calculate statistics and save
        self._save_bootstrap_summary(pairwise_results, all_results_file)

        self.logger.info("Bootstrap Fst计算完成|Bootstrap Fst calculation completed")

        return {
            'bootstrap_results': str(all_results_file),
            'bootstrap_summary': str(all_results_file)
        }

    def _extract_fst_from_file(self, fst_file: Path, iteration: int) -> Dict[Tuple[str, str], float]:
        """
        从PLINK .fst文件中提取Fst值
        Extract Fst values from PLINK .fst file

        Args:
            fst_file: .fst文件路径|.fst file path
            iteration: 迭代编号|Iteration number

        Returns:
            群体对到Fst值的映射|Population pair to Fst value mapping
        """
        results = {}

        if not fst_file.exists():
            return results

        try:
            with open(fst_file, 'r') as f:
                # 跳过表头|Skip header
                lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

                # 计算平均Fst|Calculate mean Fst
                fst_values = []
                for line in lines[1:]:  # 跳过第一行表头|Skip first header line
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            # 最后一列是Fst值|Last column is Fst value
                            fst_val_str = parts[-1]
                            # 跳过nan和缺失值|Skip nan and missing values
                            if fst_val_str.lower() == 'nan':
                                continue
                            fst_val = float(fst_val_str)
                            # 检查NMISS列（第4列），如果为0则跳过该位点
                            # Check NMISS column (4th column), skip site if NMISS is 0
                            nmiss = int(parts[3])
                            if nmiss == 0:
                                continue
                            fst_values.append(fst_val)
                        except (ValueError, IndexError):
                            continue

                if fst_values:
                    # 对于bootstrap，我们计算所有位点的平均Fst
                    # For bootstrap, we calculate mean Fst across all sites
                    mean_fst = sum(fst_values) / len(fst_values)

                    # 假设只有两个群体比较，或者取所有群体的平均值
                    # 这里简化处理：使用所有群体的平均Fst作为key
                    # In practice, need to handle multiple population pairs
                    pop_pairs = self._get_population_pairs()
                    for pair in pop_pairs:
                        results[pair] = mean_fst
                    self.logger.debug(f"从|From {fst_file.name} 提取了|extracted {len(fst_values)} 个有效FST值|valid FST values，平均值|mean: {mean_fst:.6f}")
                else:
                    self.logger.warning(f"从|From {fst_file.name} 未能提取到有效FST值|Failed to extract valid FST values")

        except Exception as e:
            self.logger.error(f"提取Fst结果失败|Failed to extract Fst results: {e}")

        return results

    def _get_population_pairs(self) -> List[Tuple[str, str]]:
        """
        获取所有群体对
        Get all population pairs

        Returns:
            群体对列表|List of population pairs
        """
        pairs = []
        pops = self.filtered_populations

        if len(pops) == 2:
            pairs.append((pops[0], pops[1]))
        else:
            # 多群体：生成所有两两组合|Multiple populations: generate all pairs
            for i in range(len(pops)):
                for j in range(i + 1, len(pops)):
                    pairs.append((pops[i], pops[j]))

        return pairs

    def _extract_fst_from_file_pair(self, fst_file: Path, iteration: int,
                                    pop_pair: Tuple[str, str]) -> Dict[Tuple[str, str], float]:
        """
        从单个群体对的PLINK .fst文件中提取Fst值
        Extract Fst value from a single population pair's PLINK .fst file

        Args:
            fst_file: .fst文件路径|.fst file path
            iteration: 迭代编号|Iteration number
            pop_pair: 群体对|Population pair

        Returns:
            群体对到Fst值的映射|Population pair to Fst value mapping
        """
        results = {}

        if not fst_file.exists():
            return results

        try:
            with open(fst_file, 'r') as f:
                # 跳过表头|Skip header
                lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

                # 计算平均Fst|Calculate mean Fst
                fst_values = []
                for line in lines[1:]:  # 跳过第一行表头|Skip first header line
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            # 最后一列是Fst值|Last column is Fst value
                            fst_val_str = parts[-1]
                            # 跳过nan和缺失值|Skip nan and missing values
                            if fst_val_str.lower() == 'nan':
                                continue
                            fst_val = float(fst_val_str)
                            # 检查NMISS列（第4列），如果为0则跳过该位点
                            # Check NMISS column (4th column), skip site if NMISS is 0
                            nmiss = int(parts[3])
                            if nmiss == 0:
                                continue
                            fst_values.append(fst_val)
                        except (ValueError, IndexError):
                            continue

                if fst_values:
                    # 计算该群体对的平均Fst|Calculate mean Fst for this population pair
                    mean_fst = sum(fst_values) / len(fst_values)
                    results[pop_pair] = mean_fst
                    self.logger.debug(f"从|From {fst_file.name} ({pop_pair[0]} vs {pop_pair[1]}) 提取了|extracted {len(fst_values)} 个有效FST值|valid FST values，平均值|mean: {mean_fst:.6f}")
                else:
                    self.logger.warning(f"从|From {fst_file.name} ({pop_pair[0]} vs {pop_pair[1]}) 未能提取到有效FST值|Failed to extract valid FST values")

        except Exception as e:
            self.logger.error(f"提取Fst结果失败|Failed to extract Fst results: {e}")

        return results

    def _save_bootstrap_summary(self, pairwise_results: Dict[Tuple[str, str], List[float]],
                                 output_file: Path):
        """
        保存bootstrap统计结果
        Save bootstrap statistical results

        Args:
            pairwise_results: 群体对到Fst值列表的映射|Population pair to Fst value list mapping
            output_file: 输出文件路径|Output file path
        """
        import numpy as np

        with open(output_file, 'w') as f:
            # 写入表头|Write header
            f.write("Population1\tPopulation2\tMean_Fst\tStd\tMin\tMax\tMedian\n")

            # 计算统计量|Calculate statistics
            for (pop1, pop2), fst_values in sorted(pairwise_results.items()):
                if not fst_values:
                    continue

                fst_array = np.array(fst_values)
                mean_fst = np.mean(fst_array)
                std_fst = np.std(fst_array)
                min_fst = np.min(fst_array)
                max_fst = np.max(fst_array)
                median_fst = np.median(fst_array)

                f.write(f"{pop1}\t{pop2}\t{mean_fst:.6f}\t{std_fst:.6f}\t"
                       f"{min_fst:.6f}\t{max_fst:.6f}\t{median_fst:.6f}\n")

                self.logger.info(f"  {pop1} vs {pop2}: "
                               f"Mean={mean_fst:.6f}, Std={std_fst:.6f}, "
                               f"Min={min_fst:.6f}, Max={max_fst:.6f}")

        self.logger.info(f"Bootstrap统计结果已保存|Bootstrap statistics saved: {output_file}")

    def run_analysis(self) -> Dict[str, str]:
        """
        运行完整的Fst计算流程
        Run complete Fst calculation pipeline

        Returns:
            输出文件字典|Output files dictionary
        """
        self.logger.info("开始Fst计算流程|Starting Fst calculation pipeline")

        output_files = {}

        # 步骤1: VCF转PLINK格式|Step 1: VCF to PLINK
        if not self.step1_vcf_to_plink():
            return output_files

        output_files['plink_bed'] = f"{self.plink_prefix}.bed"
        output_files['plink_bim'] = f"{self.plink_prefix}.bim"
        output_files['plink_fam'] = f"{self.plink_prefix}.fam"

        # 步骤2: 质量控制|Step 2: Quality control
        if not self.step2_quality_control():
            return output_files

        output_files['plink_qc_bed'] = f"{self.plink_prefix}.bed"
        output_files['plink_qc_bim'] = f"{self.plink_prefix}.bim"
        output_files['plink_qc_fam'] = f"{self.plink_prefix}.fam"

        # 步骤3: LD pruning|Step 3: LD pruning
        if not self.step3_ld_pruning():
            return output_files

        if self.plink_pruned_prefix != self.plink_prefix:
            output_files['plink_pruned_bed'] = f"{self.plink_pruned_prefix}.bed"
            output_files['plink_pruned_bim'] = f"{self.plink_pruned_prefix}.bim"
            output_files['plink_pruned_fam'] = f"{self.plink_pruned_prefix}.fam"

        # 步骤4: Fst计算|Step 4: Fst calculation
        if self.should_bootstrap:
            # Bootstrap模式|Bootstrap mode
            bootstrap_files = self.run_bootstrap()
            output_files.update(bootstrap_files)
        else:
            # 直接计算|Direct calculation
            if not self.step4_calculate_fst():
                return output_files

            output_files['fst'] = f"{self.config.output_path / 'fst_result'}.fst"
            output_files['plink_pop_file'] = str(self.config.output_path / 'population_filtered.plink')

        self.logger.info("Fst计算流程完成|Fst calculation pipeline completed")

        return output_files
