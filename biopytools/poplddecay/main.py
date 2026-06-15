"""
PopLDdecay主程序模块|PopLDdecay Main Module
"""

from pathlib import Path

from .config import PopLDdecayConfig
from .utils import PopLDdecayLogger, validate_subpop_file, parse_subpop_file, check_perl_modules
from .calculator import PopLDdecayCalculator
from .plotter import PopLDdecayPlotter
from .ld_threshold import LDThresholdRecommender


class PopLDdecayRunner:
    """PopLDdecay运行器|PopLDdecay Runner"""

    def __init__(self, config: PopLDdecayConfig, logger):
        """初始化运行器|Initialize runner

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 初始化各模块|Initialize modules
        self.calculator = PopLDdecayCalculator(config, logger)
        self.plotter = PopLDdecayPlotter(config, logger)

    def run(self) -> bool:
        """运行完整的LD衰减分析流程|Run complete LD decay analysis pipeline

        Returns:
            bool: 是否成功|Success
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始PopLDdecay分析|Start PopLDdecay analysis")
            self.logger.info("=" * 60)

            # 存储输出文件信息|Store output file information
            output_results = []

            # 如果指定了子群体文件，先计算所有样本的LD衰减|If subpopulation file is specified, calculate LD decay for all samples first
            if self.config.subpop_file:
                self.logger.info("[步骤1/5] 验证子群体文件|[Step 1/5] Validating subpopulation file")
                if not validate_subpop_file(self.config.subpop_path, self.logger):
                    return False

                # 解析子群体文件|Parse subpopulation file
                populations = parse_subpop_file(self.config.subpop_path, self.logger)

                if not populations:
                    self.logger.error("未能识别子群体|Failed to identify subpopulations")
                    return False

                # 步骤2: 计算所有样本的LD衰减|Step 2: Calculate LD decay for all samples
                self.logger.info("[步骤2/5] 计算所有样本的LD衰减|[Step 2/5] Calculating LD decay for all samples")
                all_config = self._create_all_samples_config()
                all_calculator = PopLDdecayCalculator(all_config, self.logger)

                if all_calculator.run():
                    output_files = all_config.get_output_files()
                    output_results.append(('all', '所有样本|All samples', output_files))
                else:
                    self.logger.error("所有样本LD衰减计算失败|Failed to calculate LD decay for all samples")
                    return False

                # 步骤3: 为每个子群体分别计算LD衰减|Step 3: Calculate LD decay for each subpopulation
                self.logger.info("[步骤3/5] 计算各子群体的LD衰减|[Step 3/5] Calculating LD decay for each subpopulation")

                for pop_name, samples in sorted(populations.items()):
                    self.logger.info(f"计算群体|Calculating population: {pop_name} ({len(samples)}个样本|samples)")

                    # 为当前群体创建配置|Create config for current population
                    pop_config = self._create_population_config(pop_name, samples)
                    pop_calculator = PopLDdecayCalculator(pop_config, self.logger)

                    if pop_calculator.run():
                        output_files = pop_config.get_output_files()
                        output_results.append(('population', pop_name, output_files))
                    else:
                        self.logger.warning(f"群体{pop_name}计算失败|Failed to calculate population {pop_name}")
                        # 继续计算其他群体|Continue calculating other populations

                # 步骤4: 绘图（如果启用）|Step 4: Plot (if enabled)
                if self.config.plot:
                    self.logger.info("[步骤4/6] 绘制LD衰减图|[Step 4/6] Plotting LD decay figures")
                    self._plot_multiple_populations(output_results)
                else:
                    self.logger.info("[步骤4/6] 跳过绘图|[Step 4/6] Skipping plotting")

                # 步骤5: 整合所有群体的stat文件|Step 5: Merge all population stat files
                self.logger.info("[步骤5/6] 整合所有群体的LD衰减结果|[Step 5/6] Merging LD decay results from all populations")
                merged_stat_file = self._merge_stat_files(output_results)

            else:
                # 没有子群体文件，直接计算|No subpopulation file, calculate directly
                self.logger.info("[步骤1/2] 计算LD衰减|[Step 1/2] Calculating LD decay")
                if not self.calculator.run():
                    return False

                output_files = self.config.get_output_files()
                output_results.append(('single', '单一群体|Single population', output_files))

                # 步骤2: 绘图（如果启用）|Step 2: Plot (if enabled)
                if self.config.plot:
                    self.logger.info("[步骤2/2] 绘制LD衰减图|[Step 2/2] Plotting LD decay figure")
                    if check_perl_modules(self.logger):
                        if not self.plotter.plot_single_population():
                            self.logger.warning("绘图失败，但统计文件已生成|Plotting failed, but statistics file generated")
                    else:
                        self.logger.warning("Perl模块不可用，跳过绘图|Perl modules not available, skipping plotting")
                else:
                    self.logger.info("[步骤2/2] 跳过绘图|[Step 2/2] Skipping plotting")

                # 单一群体不需要整合|No need to merge for single population
                merged_stat_file = None

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("PopLDdecay分析完成|PopLDdecay analysis completed")
            self.logger.info("=" * 60)

            # 输出文件|Output files
            self.logger.info(f"输出文件|Output files:")
            for pop_type, pop_name, files in output_results:
                self.logger.info(f"  {pop_name}:")
                self.logger.info(f"    统计文件|Statistics: {files['stat']}")
                if self.config.plot:
                    self.logger.info(f"    图像文件|Figure: {files['figure']}")

            # 输出整合文件|Output merged file
            if merged_stat_file:
                self.logger.info(f"  整合文件|Merged file: {merged_stat_file}")

            # 步骤6: 阈值推荐|Step 6: Threshold recommendation
            if self.config.recommend_threshold:
                self.logger.info("")
                self.logger.info("[步骤6/6] LD阈值推荐|[Step 6/6] LD threshold recommendation")
                self._generate_threshold_recommendations(output_results)

            return True

        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _create_all_samples_config(self):
        """创建计算所有样本的配置对象|Create configuration object for all samples calculation

        Returns:
            PopLDdecayConfig: 新的配置对象（不包含子群体）|New configuration object (without subpopulation)
        """
        from copy import deepcopy

        # 深拷贝配置|Deep copy configuration
        all_config = deepcopy(self.config)

        # 修改输出文件前缀，添加_all后缀|Modify output prefix, add _all suffix
        output_prefix_str = str(all_config.output_path)
        all_config.output_path = Path(output_prefix_str + "_all")

        # 清空子群体文件|Clear subpopulation file
        all_config.subpop_file = ""
        all_config.subpop_path = None

        return all_config

    def _create_population_config(self, pop_name: str, samples: list):
        """创建计算单个群体的配置对象|Create configuration object for single population calculation

        Args:
            pop_name: 群体名称|Population name
            samples: 样本列表|Sample list

        Returns:
            PopLDdecayConfig: 新的配置对象|New configuration object
        """
        from copy import deepcopy
        import tempfile

        # 深拷贝配置|Deep copy configuration
        pop_config = deepcopy(self.config)

        # 修改输出文件前缀，添加群体名称|Modify output prefix, add population name
        output_prefix_str = str(pop_config.output_path)
        # 清理群体名称中的特殊字符|Clean special characters in population name
        clean_pop_name = pop_name.replace(' ', '_').replace('/', '_')
        pop_config.output_path = Path(f"{output_prefix_str}_{clean_pop_name}")

        # 创建临时子群体文件|Create temporary subpopulation file
        temp_subpop_file = tempfile.NamedTemporaryFile(
            mode='w',
            suffix=f'_{clean_pop_name}.txt',
            delete=False,
            dir=pop_config.output_path.parent
        )

        # 写入样本列表|Write sample list
        for sample in samples:
            temp_subpop_file.write(f"{sample}\n")

        temp_subpop_file.close()

        # 设置子群体文件|Set subpopulation file
        pop_config.subpop_file = str(temp_subpop_file.name)
        pop_config.subpop_path = Path(temp_subpop_file.name)

        return pop_config

    def _get_sample_count_from_input(self) -> int:
        """从输入文件获取样本数量|Get sample count from input file

        Returns:
            int: 样本数量|Number of samples
        """
        try:
            if self.config.input_type == 'vcf':
                # 从VCF文件头获取样本数|Get sample count from VCF header
                import gzip
                open_func = gzip.open if self.config.input_path.suffix == '.gz' else open
                with open_func(self.config.input_path, 'rt') as f:
                    for line in f:
                        if line.startswith('#CHROM'):
                            parts = line.strip().split('\t')
                            # VCF格式: 从第10列开始是样本|VCF format: samples start from column 10
                            return max(0, len(parts) - 9)
            else:
                # Genotype格式: 计算行数|Genotype format: count lines
                with open(self.config.input_path, 'r') as f:
                    return sum(1 for line in f if not line.startswith('#') and line.strip())
        except Exception as e:
            self.logger.warning(f"无法从输入文件获取样本数量: {e}|Failed to get sample count from input file: {e}")
            return 0

    def _recommend_threshold_for_stat(self, stat_file: Path, pop_name: str, n: int) -> dict:
        """为单个stat文件推荐LD阈值|Recommend LD threshold for single stat file

        Args:
            stat_file: stat文件路径|Stat file path
            pop_name: 群体名称|Population name
            n: 样本数量|Number of samples

        Returns:
            dict: 推荐结果|Recommendation result
        """
        try:
            # 创建recommender|Create recommender
            # PopLDdecay输出格式: #Dist, Mean_r^2|PopLDdecay output format: #Dist, Mean_r^2
            recommender = LDThresholdRecommender(
                n=n,
                dist_col="#Dist",  # PopLDdecay格式|PopLDdecay format
                value_col="Mean_r^2",  # PopLDdecay格式|PopLDdecay format
                dist_unit="bp",  # PopLDdecay输出单位是bp|PopLDdecay output unit is bp
                verbose=False  # 不输出详细日志|No verbose output
            )

            # 读取stat文件|Read stat file
            result = recommender.recommend(stat_file)

            # 添加群体名称和样本数|Add population name and sample count
            result['population'] = pop_name
            result['n'] = n

            return result

        except Exception as e:
            self.logger.warning(f"为{pop_name}推荐阈值失败|Failed to recommend threshold for {pop_name}: {e}")
            return None

    def _generate_threshold_recommendations(self, output_results):
        """生成阈值推荐汇总表|Generate threshold recommendation summary table

        Args:
            output_results: 输出结果列表|Output results list [(pop_type, pop_name, files), ...]
        """
        self.logger.info("=" * 60)
        self.logger.info("开始LD阈值推荐|Start LD threshold recommendation")
        self.logger.info("=" * 60)

        # 存储所有推荐结果|Store all recommendation results
        recommendations = []

        # 为每个群体推荐阈值|Recommend threshold for each population
        for pop_type, pop_name, files in output_results:
            stat_file = Path(files['stat'])

            # 获取样本数量|Get sample count
            if pop_type == 'all':
                # 总群体|All samples
                n = self._get_sample_count_from_input()
            elif pop_type == 'population':
                # 子群体 - 需要从子群体文件获取样本数|Subpopulation - get sample count from subpop file
                populations = parse_subpop_file(self.config.subpop_path, self.logger)
                n = len(populations.get(pop_name, []))
            else:
                # 单一群体|Single population
                n = self._get_sample_count_from_input()

            if n == 0:
                self.logger.warning(f"无法获取{pop_name}的样本数量，跳过阈值推荐|Cannot get sample count for {pop_name}, skipping threshold recommendation")
                continue

            # 推荐阈值|Recommend threshold
            result = self._recommend_threshold_for_stat(stat_file, pop_name, n)

            if result:
                recommendations.append(result)

        # 生成汇总表|Generate summary table
        if recommendations:
            self._output_threshold_summary(recommendations)
        else:
            self.logger.warning("未能生成任何阈值推荐|Failed to generate any threshold recommendations")

    def _output_threshold_summary(self, recommendations: list):
        """输出阈值推荐汇总|Output threshold recommendation summary

        Args:
            recommendations: 推荐结果列表|Recommendation result list
        """
        import pandas as pd

        # 构建汇总表|Build summary table
        summary_data = []
        for rec in recommendations:
            best = rec['recommended']
            summary_data.append({
                'Population': rec.get('population', 'Unknown'),
                'N_Samples': rec.get('n', 'N/A'),
                'Rho': f"{rec['rho']:.4f}",
                'Background_r2': f"{rec['background_r2']:.4f}",
                'Recommended_Threshold': f"{best['threshold'].iloc[0]:.2f}",
                'Decay_Distance_kb': f"{best['decay_kb'].iloc[0]:.1f}",
                'BG_Ratio': f"{best['bg_ratio'].iloc[0]:.2f}",
                'GWAS_Window_kb': f"+/- {rec['gwas_window_kb']}"
            })

        summary_df = pd.DataFrame(summary_data)

        # 输出到日志|Output to log
        self.logger.info("")
        self.logger.info("LD阈值推荐汇总|LD Threshold Recommendation Summary:")
        self.logger.info("=" * 100)
        self.logger.info("\n" + summary_df.to_string(index=False))
        self.logger.info("")
        self.logger.info("=" * 100)

        # 保存到文件|Save to file
        output_files = self.config.get_output_files()
        output_file = Path(output_files['threshold'])

        try:
            summary_df.to_csv(output_file, sep='\t', index=False)
            self.logger.info(f"阈值推荐汇总已保存到|Threshold recommendation summary saved to: {output_file}")
        except Exception as e:
            self.logger.error(f"保存阈值推荐汇总失败|Failed to save threshold recommendation summary: {e}")

    def _plot_multiple_populations(self, output_results):
        """绘制多群体LD衰减图|Plot multi-population LD decay figures

        Args:
            output_results: 输出结果列表|Output results list [(pop_type, pop_name, files), ...]
        """
        # 检查Perl模块|Check Perl modules
        if not check_perl_modules(self.logger):
            self.logger.warning("Perl模块不可用，跳过绘图|Perl modules not available, skipping plotting")
            return

        # 准备群体列表（只包含population类型，不包含all）|Prepare population list (only include population type, not all)
        pop_list = []
        for pop_type, pop_name, files in output_results:
            if pop_type == 'population':  # 只为子群体绘图|Only plot subpopulations
                stat_file = files['stat']
                pop_id = pop_name  # 使用群体名称作为ID
                pop_list.append((stat_file, pop_id))

        # 如果没有子群体，只绘制all|If no subpopulations, only plot all
        if not pop_list:
            for pop_type, pop_name, files in output_results:
                if pop_type == 'all':
                    all_config = self._create_all_samples_config()
                    all_plotter = PopLDdecayPlotter(all_config, self.logger)
                    if not all_plotter.plot_single_population():
                        self.logger.warning("绘图失败，但统计文件已生成|Plotting failed, but statistics file generated")
                    return

        # 使用多群体绘图|Use multi-population plotting
        if len(pop_list) > 1:
            self.logger.info(f"绘制{len(pop_list)}个群体的LD衰减图|Plotting LD decay figures for {len(pop_list)} populations")
            if not self.plotter.plot_multiple_populations(pop_list):
                self.logger.warning("多群体绘图失败，尝试分别绘制各群体|Multi-population plotting failed, trying to plot each population separately")

                # 备用方案：分别绘制每个群体|Fallback: plot each population separately
                for pop_type, pop_name, files in output_results:
                    if pop_type == 'population':
                        # 创建新的plotter配置|Create new plotter config
                        from copy import deepcopy
                        plot_config = deepcopy(self.config)
                        plot_config.output_path = Path(files['stat']).with_suffix('')

                        # 创建临时plotter|Create temporary plotter
                        temp_plotter = PopLDdecayPlotter(plot_config, self.logger)

                        if not temp_plotter.plot_single_population():
                            self.logger.warning(f"绘制{pop_name}失败|Failed to plot {pop_name}")
        else:
            # 单个群体，直接使用单群体绘图|Single population, use single population plotting directly
            if not self.plotter.plot_single_population():
                self.logger.warning("绘图失败，但统计文件已生成|Plotting failed, but statistics file generated")

    def _merge_stat_files(self, output_results: list) -> str:
        """整合所有群体的stat文件为汇总表|Merge all population stat files into summary table

        输出TSV文件，包含Population/Dist/Mean_r2列，可直接用于ggplot2绘图
        Output TSV file with Population/Dist/Mean_r2 columns, ready for ggplot2 plotting

        Args:
            output_results: 输出结果列表|Output results list [(pop_type, pop_name, files), ...]

        Returns:
            str: 汇总表文件路径|Summary table file path, 如果失败返回None|returns None if failed
        """
        try:
            import pandas as pd

            # 存储所有数据|Store all data
            all_data = []

            for pop_type, pop_name, files in output_results:
                stat_file = Path(files['stat'])

                # 确定群体标识|Determine population identifier
                if pop_type == 'all':
                    population_id = 'all'
                elif pop_type == 'population':
                    population_id = pop_name
                else:
                    population_id = 'single'

                try:
                    # 读取stat文件|Read stat file
                    if stat_file.suffix == '.gz':
                        df = pd.read_csv(stat_file, sep='\t', compression='gzip')
                    else:
                        df = pd.read_csv(stat_file, sep='\t')

                    # 重命名列，去掉特殊字符|Rename columns, remove special characters
                    column_mapping = {
                        '#Dist': 'Dist',
                        'Mean_r^2': 'Mean_r2',
                        'Mean_D\'': 'Mean_Dprime',
                        'Sum_r^2': 'Sum_r2',
                        'Sum_D\'': 'Sum_Dprime'
                    }
                    df.rename(columns=column_mapping, inplace=True)

                    # 保留所有列，供后续分析使用|Keep all columns for downstream analysis

                    # 添加群体列|Add population column
                    df.insert(0, 'Population', population_id)

                    all_data.append(df)

                    self.logger.info(f"  已添加|Added: {population_id} ({len(df)} 行|rows)")

                except Exception as e:
                    self.logger.warning(f"读取{stat_file}失败|Failed to read {stat_file}: {e}")
                    continue

            if not all_data:
                self.logger.error("没有成功读取任何stat文件|Failed to read any stat files")
                return None

            # 合并所有数据|Merge all data
            merged_df = pd.concat(all_data, ignore_index=True)

            # 输出文件路径|Output file path
            summary_file = self.config.output_path.parent / f"{self.config.output_path.name}_summary.tsv"

            # 保存汇总表（未压缩TSV，方便R/ggplot2直接读取）
            # Save summary table (uncompressed TSV, easy for R/ggplot2 to read directly)
            merged_df.to_csv(summary_file, sep='\t', index=False)

            self.logger.info(f"汇总表保存成功|Summary table saved successfully")
            self.logger.info(f"  文件路径|File path: {summary_file}")
            self.logger.info(f"  总行数|Total rows: {len(merged_df)}")
            self.logger.info(f"  包含{len(all_data)}个群体|Contains {len(all_data)} populations")
            self.logger.info(f"  列|Columns: {', '.join(merged_df.columns.tolist())}")
            self.logger.info(f"  ggplot2示例|ggplot2 example: ggplot(df, aes(Dist, Mean_r2, color=Population)) + geom_line()")

            return str(summary_file)

        except Exception as e:
            self.logger.error(f"整合stat文件失败|Failed to merge stat files: {e}")
            import traceback
            traceback.print_exc()
            return None


def main(args) -> int:
    """主函数入口|Main function entry point

    Args:
        args: 命令行参数|Command line arguments

    Returns:
        int: 退出码|Exit code (0=success, 1=failure)
    """
    # 创建配置对象|Create configuration object
    config = PopLDdecayConfig(
        input_file=args.input,
        output_prefix=args.output,
        input_type=args.type,
        max_dist=args.max_dist,
        min_maf=args.min_maf,
        max_het=args.max_het,
        max_miss=args.max_miss,
        subpop_file=args.subpop if args.subpop else "",
        out_type=args.out_type,
        plot=not args.no_plot,
        bin1=args.bin1,
        bin2=args.bin2,
        break_point=args.break_point,
        max_x=args.max_x,
        measure=args.measure,
        method=args.method,
        percentile=args.percentile,
        recommend_threshold=args.recommend_threshold
    )

    # 验证配置|Validate configuration
    try:
        config.validate()
    except ValueError as e:
        print(f"配置错误|Configuration error:\n{e}")
        return 1

    # 设置日志|Setup logging
    output_files = config.get_output_files()
    log_file = Path(output_files['log'])

    logger_manager = PopLDdecayLogger(log_file)
    logger = logger_manager.get_logger()

    # 运行分析|Run analysis
    runner = PopLDdecayRunner(config, logger)
    success = runner.run()

    return 0 if success else 1
