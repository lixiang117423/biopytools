"""
端粒识别主程序模块|Telomere Finder Main Module
"""

import os
import sys
import time
from pathlib import Path
from typing import Optional
from .config import TelomereFinderConfig
from .utils import (TelomereLogger, CommandRunner, OutputValidator, CladeDatabase,
                    select_prioritized_telomere_repeat, get_plant_telomere_priority_list,
                    evaluate_telomere_search_results, build_conda_command)


class TelomereFinder:
    """端粒识别主类|Main Telomere Finder Class"""

    def __init__(self, **kwargs):
        """
        初始化端粒识别器|Initialize telomere finder

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = TelomereFinderConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = TelomereLogger(
            self.config.output_dir,
            self.config.verbose,
            self.config.log_file
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

        # 初始化输出验证器|Initialize output validator
        self.output_validator = OutputValidator(self.logger)

        # 输出配置信息|Output configuration information
        self._log_configuration()

    def _log_configuration(self):
        """记录配置信息|Log configuration information"""
        self.logger.info("=" * 60)
        self.logger.info("端粒识别分析|Telomere Finder Analysis")
        self.logger.info("=" * 60)
        self.logger.info(f"模式|Mode: {self.config.mode}")
        self.logger.info(f"基因组文件|Genome file: {self.config.genome_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
        self.logger.info(f"tidk 路径|tidk path: {self.config.tidk_path}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info("=" * 60)

    def _build_explore_command(self) -> list:
        """
        构建 explore 模式命令|Build explore mode command

        Returns:
            list: 命令列表|Command list
        """
        output_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_explore.tsv")

        # 构建命令参数列表|Build command arguments list
        args = [
            'explore',
            self.config.genome_file,
            '--minimum', str(self.config.explore_min_length),
            '--maximum', str(self.config.explore_max_length),
            '--threshold', str(self.config.explore_threshold),
            '--distance', str(self.config.explore_distance)
        ]

        # 使用build_conda_command包装命令|Use build_conda_command to wrap command
        command = build_conda_command(self.config.tidk_path, args)

        return command, output_file

    def _build_find_command(self) -> list:
        """
        构建 find 模式命令|Build find mode command

        Returns:
            list: 命令列表|Command list
        """
        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_telomeric_repeat_windows.tsv"
        )

        # 构建命令参数列表|Build command arguments list
        args = [
            'find',
            self.config.genome_file,
            '--clade', self.config.clade,
            '--window', str(self.config.window_size),
            '--output', self.config.output_prefix,
            '--dir', self.config.output_dir
        ]

        # 使用build_conda_command包装命令|Use build_conda_command to wrap command
        command = build_conda_command(self.config.tidk_path, args)

        return command, output_file

    def _build_search_command(self) -> list:
        """
        构建 search 模式命令|Build search mode command

        Returns:
            list: 命令列表|Command list
        """
        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_search_windows.{self.config.output_format}"
        )

        # 构建命令参数列表|Build command arguments list
        args = [
            'search',
            self.config.genome_file,
            '--string', self.config.search_string,
            '--window', str(self.config.window_size),
            '--output', self.config.output_prefix,
            '--dir', self.config.output_dir,
            '--extension', self.config.output_format
        ]

        # 使用build_conda_command包装命令|Use build_conda_command to wrap command
        command = build_conda_command(self.config.tidk_path, args)

        return command, output_file

    def _build_plot_command(self) -> list:
        """
        构建 plot 模式命令|Build plot mode command

        Returns:
            list: 命令列表|Command list
        """
        # tidk plot 生成的文件名格式|tidk plot output file format
        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.svg"
        )

        # 构建命令参数列表|Build command arguments list
        args = [
            'plot',
            '--tsv', self.config.tsv_file,
            '--height', str(self.config.plot_height),
            '--width', str(self.config.plot_width),
            '--fontsize', str(self.config.plot_fontsize),
            '--strokewidth', str(self.config.plot_strokewidth),
            '--output', os.path.join(self.config.output_dir, self.config.output_prefix)
        ]

        # 使用build_conda_command包装命令|Use build_conda_command to wrap command
        command = build_conda_command(self.config.tidk_path, args)

        return command, output_file

    def run_explore(self) -> bool:
        """
        运行 explore 模式|Run explore mode

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤|Step: 探索端粒重复序列|Explore telomeric repeats")
        self.logger.info("=" * 60)

        try:
            command, output_file = self._build_explore_command()

            self.logger.info(f"参数|Parameters:")
            self.logger.info(f"  最小长度|Min length: {self.config.explore_min_length}")
            self.logger.info(f"  最大长度|Max length: {self.config.explore_max_length}")
            self.logger.info(f"  阈值|Threshold: {self.config.explore_threshold}")
            self.logger.info(f"  距离比例|Distance ratio: {self.config.explore_distance}")

            # 运行命令并捕获输出|Run command and capture output
            self.logger.info(f"执行中|Executing: 探索端粒重复序列|Exploring telomeric repeats")
            self.logger.debug(f"完整命令|Full command: {' '.join(command)}")

            import subprocess
            start_time = time.time()

            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                shell=False  # 使用列表形式时必须使用shell=False|Must use shell=False with list
            )

            elapsed_time = time.time() - start_time

            # 处理 stderr 输出|Handle stderr output
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if line and 'Warning' not in line and 'warning' not in line:
                        self.logger.warning(line)

            # 将 stdout 保存到文件|Save stdout to file
            if result.stdout:
                with open(output_file, 'w') as f:
                    f.write(result.stdout)

                # 同时输出到日志|Also output to log
                for line in result.stdout.strip().split('\n'):
                    if line:
                        self.logger.info(line)

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed with return code {result.returncode}")
                return False

            self.logger.info(f"命令执行完成|Command completed in {elapsed_time:.2f} 秒|seconds")

            # 验证输出|Validate output
            if os.path.exists(output_file):
                stats = self.output_validator.get_file_stats(output_file)
                self.logger.info(f"输出文件|Output file: {output_file}")
                if 'lines' in stats:
                    self.logger.info(f"结果行数|Result lines: {stats['lines']}")
                self.logger.info("探索完成|Exploration completed successfully")
            else:
                self.logger.warning(f"输出文件未生成|Output file not generated: {output_file}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"探索失败|Exploration failed: {str(e)}")
            return False

    def run_find(self) -> bool:
        """
        运行 find 模式|Run find mode

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤|Step: 查找端粒重复序列|Find telomeric repeats")
        self.logger.info("=" * 60)

        try:
            command, output_file = self._build_find_command()

            self.logger.info(f"参数|Parameters:")
            self.logger.info(f"  分类群|Clade: {self.config.clade}")
            self.logger.info(f"  窗口大小|Window size: {self.config.window_size}")

            # 运行命令|Run command
            success = self.cmd_runner.run_command(
                command,
                description="查找端粒重复序列|Finding telomeric repeats"
            )

            if success:
                # 验证输出|Validate output
                expected_files = [
                    os.path.join(self.config.output_dir, f"{self.config.output_prefix}_telomeric_repeat_windows.tsv"),
                    os.path.join(self.config.output_dir, f"{self.config.output_prefix}_telomeric_repeat_counts.tsv")
                ]
                self.output_validator.validate_output_files(expected_files)
                self.logger.info("查找完成|Finding completed successfully")

            return success

        except Exception as e:
            self.logger.error(f"查找失败|Finding failed: {str(e)}")
            return False

    def run_search(self) -> bool:
        """
        运行 search 模式|Run search mode

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤|Step: 搜索端粒重复序列|Search telomeric repeats")
        self.logger.info("=" * 60)

        try:
            command, output_file = self._build_search_command()

            self.logger.info(f"参数|Parameters:")
            self.logger.info(f"  搜索字符串|Search string: {self.config.search_string}")
            self.logger.info(f"  窗口大小|Window size: {self.config.window_size}")
            self.logger.info(f"  输出格式|Output format: {self.config.output_format}")

            # 运行命令|Run command
            success = self.cmd_runner.run_command(
                command,
                description="搜索端粒重复序列|Searching telomeric repeats"
            )

            if success:
                # 验证输出|Validate output
                # tidk search 生成的实际文件名|Actual file names generated by tidk search
                expected_files = [
                    os.path.join(self.config.output_dir, f"{self.config.output_prefix}_telomeric_repeat_windows.{self.config.output_format}"),
                    os.path.join(self.config.output_dir, f"{self.config.output_prefix}_telomeric_repeat_counts.{self.config.output_format}")
                ]
                self.output_validator.validate_output_files(expected_files)
                self.logger.info("搜索完成|Search completed successfully")

            return success

        except Exception as e:
            self.logger.error(f"搜索失败|Search failed: {str(e)}")
            return False

    def run_plot(self) -> bool:
        """
        运行 plot 模式|Run plot mode

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤|Step: 绘制端粒分布图|Plot telomere distribution")
        self.logger.info("=" * 60)

        try:
            command, output_file = self._build_plot_command()

            self.logger.info(f"参数|Parameters:")
            self.logger.info(f"  TSV 文件|TSV file: {self.config.tsv_file}")
            self.logger.info(f"  图像高度|Plot height: {self.config.plot_height}")
            self.logger.info(f"  图像宽度|Plot width: {self.config.plot_width}")
            self.logger.info(f"  字体大小|Font size: {self.config.plot_fontsize}")
            self.logger.info(f"  线条宽度|Stroke width: {self.config.plot_strokewidth}")

            # 运行命令|Run command
            success = self.cmd_runner.run_command(
                command,
                description="绘制端粒分布图|Plotting telomere distribution"
            )

            if success:
                # 验证输出|Validate output
                if os.path.exists(output_file):
                    stats = self.output_validator.get_file_stats(output_file)
                    self.logger.info(f"输出文件|Output file: {output_file}")
                    self.logger.info(f"文件大小|File size: {stats['size_mb']:.2f} MB")
                    self.logger.info("绘图完成|Plotting completed successfully")
                else:
                    self.logger.warning(f"输出文件未生成|Output file not generated: {output_file}")
                    return False

            return success

        except Exception as e:
            self.logger.error(f"绘图失败|Plotting failed: {str(e)}")
            return False

    def run_pipeline(self) -> bool:
        """
        运行完整的自动化流程|Run complete automated pipeline
        1. 探索端粒序列|Explore telomeric repeats
        2. 多轮尝试搜索端粒位置|Multi-round search for telomere positions
        3. 绘制分布图|Plot distribution

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始自动化端粒识别流程|Starting automated telomere identification pipeline")
        self.logger.info("=" * 60)

        # 步骤1: 探索端粒序列|Step 1: Explore telomeric repeats
        self.logger.info("")
        self.logger.info("步骤1/3|Step 1/3: 探索端粒重复序列|Exploring telomeric repeats")
        self.logger.info("")

        explore_success = self.run_explore()
        if not explore_success:
            self.logger.error("探索失败，终止流程|Exploration failed, terminating pipeline")
            return False

        # 读取探索结果|Read exploration results
        explore_output = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_explore.tsv")
        if not os.path.exists(explore_output):
            self.logger.error(f"探索结果文件不存在|Exploration output not found: {explore_output}")
            return False

        # 解析探索结果，获取候选端粒序列|Parse exploration results to get candidate telomeric repeats
        explore_candidates = []
        try:
            with open(explore_output, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    for line in lines[1:]:  # 跳过header|Skip header
                        parts = line.strip().split('\t')
                        if len(parts) >= 1 and parts[0].strip():
                            explore_candidates.append(parts[0].strip())

                    if not explore_candidates:
                        self.logger.error("探索结果为空|Exploration results are empty")
                        return False

                    self.logger.info(f"发现 {len(explore_candidates)} 个探索候选序列|Found {len(explore_candidates)} explored candidate sequences")
                    for i, repeat in enumerate(explore_candidates[:5], 1):
                        self.logger.info(f"  候选 {i}|Candidate {i}: {repeat}")
                else:
                    self.logger.error("探索结果为空|Exploration results are empty")
                    return False
        except Exception as e:
            self.logger.error(f"读取探索结果失败|Failed to read exploration results: {str(e)}")
            return False

        # 步骤2: 多轮尝试搜索|Step 2: Multi-round search attempts
        self.logger.info("")
        self.logger.info("步骤2/3|Step 2/3: 智能搜索端粒位置|Smart telomere position search")
        self.logger.info("")

        # 构建尝试序列列表|Build sequence list to try
        plant_priority_list = get_plant_telomere_priority_list()
        self.logger.info(f"植物端粒优先级列表|Plant telomere priority list: {plant_priority_list[:10]}...")

        # 优先级策略|Priority strategy:
        # 1. 典型植物端粒 TTTAGGG|Typical plant telomere
        # 2. 其他植物端粒（最多5个）| Other plant telomeres (max 5)
        # 3. 探索结果的第一个|First explored result
        sequences_to_try = []

        # 第一优先: TTTAGGG|First priority: TTTAGGG
        if 'TTTAGGG' not in sequences_to_try:
            sequences_to_try.append('TTTAGGG')

        # 第二优先: 其他植物端粒（最多5个）| Second priority: other plant telomeres (max 5)
        for seq in plant_priority_list:
            if seq not in ['TTTAGGG', 'CCCTAAA'] and seq not in sequences_to_try:
                sequences_to_try.append(seq)
                if len(sequences_to_try) >= 6:  # TTTAGGG + 5 others
                    break

        # 最后: 探索结果的第一个|Last: first explored result
        if explore_candidates and explore_candidates[0] not in sequences_to_try:
            sequences_to_try.append(explore_candidates[0])

        self.logger.info(f"将按顺序尝试 {len(sequences_to_try)} 个端粒序列|Will try {len(sequences_to_try)} telomere sequences in order:")
        for i, seq in enumerate(sequences_to_try, 1):
            source = "Plant" if seq in plant_priority_list else "Explored"
            self.logger.info(f"  {i}. {seq} ({source})")

        # 多轮尝试|Multi-round attempts
        successful_repeat = None
        successful_windows_file = None

        original_mode = self.config.mode
        original_search_string = self.config.search_string

        for attempt, repeat_seq in enumerate(sequences_to_try, 1):
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info(f"尝试 {attempt}/{len(sequences_to_try)}|Attempt {attempt}/{len(sequences_to_try)}: {repeat_seq}")
            self.logger.info("=" * 60)

            # 执行搜索|Perform search
            self.config.mode = 'search'
            self.config.search_string = repeat_seq

            search_success = self.run_search()

            if not search_success:
                self.logger.warning(f"使用 {repeat_seq} 搜索失败|Search with {repeat_seq} failed")
                if attempt < len(sequences_to_try):
                    self.logger.info(f"继续尝试下一个序列|Continuing to next sequence...")
                continue

            # 评估搜索结果|Evaluate search results
            windows_file = os.path.join(
                self.config.output_dir,
                f"{self.config.output_prefix}_telomeric_repeat_windows.tsv"
            )

            evaluation = evaluate_telomere_search_results(windows_file, self.logger)

            if evaluation['valid']:
                self.logger.info("")
                self.logger.info(f"成功! 使用端粒序列 {repeat_seq} 找到了足够的端粒信号|Success! Found sufficient telomere signals using {repeat_seq}")
                successful_repeat = repeat_seq
                successful_windows_file = windows_file
                break
            else:
                self.logger.warning(f"使用 {repeat_seq} 的搜索结果不理想|Search results with {repeat_seq} are not ideal")
                self.logger.warning(f"原因|Reason: {evaluation['reason']}")
                if attempt < len(sequences_to_try):
                    self.logger.info(f"尝试下一个序列|Trying next sequence...")

        # 恢复原始配置|Restore original config
        self.config.mode = original_mode
        self.config.search_string = original_search_string

        # 检查是否成功|Check if successful
        if not successful_repeat:
            self.logger.error("")
            self.logger.error("=" * 60)
            self.logger.error("所有尝试均失败|All attempts failed")
            self.logger.error("=" * 60)
            return False

        # 步骤3: 绘制分布图|Step 3: Plot distribution
        self.logger.info("")
        self.logger.info("步骤3/3|Step 3/3: 绘制端粒分布图|Plotting telomere distribution")
        self.logger.info("")

        if not successful_windows_file or not os.path.exists(successful_windows_file):
            self.logger.error(f"TSV 文件不存在|TSV file not found: {successful_windows_file}")
            return False

        # 临时修改配置为 plot 模式|Temporarily modify config to plot mode
        original_tsv = self.config.tsv_file
        self.config.mode = 'plot'
        self.config.tsv_file = successful_windows_file

        plot_success = self.run_plot()

        # 恢复原始配置|Restore original config
        self.config.mode = original_mode
        self.config.tsv_file = original_tsv

        if not plot_success:
            self.logger.warning("绘图失败，但流程已完成|Plotting failed, but pipeline completed")
            return True  # 绘图失败不算整体失败

        # 流程完成|Pipeline completed
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("自动化流程全部完成|Automated pipeline completed successfully")
        self.logger.info("=" * 60)
        self.logger.info("")
        self.logger.info("最终使用的端粒序列|Final telomere repeat used: " + successful_repeat)
        self.logger.info("")
        self.logger.info("结果文件|Output files:")
        self.logger.info(f"  1. 端粒序列探索|Telomere exploration: {explore_output}")
        self.logger.info(f"  2. 端粒位置统计|Position statistics: {successful_windows_file}")
        self.logger.info(f"  3. 端粒分布图|Distribution plot: {os.path.join(self.config.output_dir, f'{self.config.output_prefix}.svg')}")
        self.logger.info("")

        return True

    def run(self) -> bool:
        """
        运行端粒识别分析|Run telomere finder analysis

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            if self.config.mode == 'explore':
                return self.run_explore()
            elif self.config.mode == 'find':
                return self.run_find()
            elif self.config.mode == 'search':
                return self.run_search()
            elif self.config.mode == 'plot':
                return self.run_plot()
            elif self.config.mode == 'pipeline':
                return self.run_pipeline()
            else:
                self.logger.error(f"未知的模式|Unknown mode: {self.config.mode}")
                return False

        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {str(e)}")
            return False


def main():
    """主函数入口|Main function entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='端粒识别分析工具|Telomere Finder Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genome.fa -o results
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-g', '--genome', required=True,
                         help='基因组序列文件|Genome sequence file (FASTA format)')

    # 模式选择|Mode selection
    mode_group = parser.add_argument_group('分析模式|Analysis mode')
    mode_group.add_argument('-m', '--mode',
                           choices=['explore', 'find', 'search', 'plot', 'pipeline'],
                           default='pipeline',
                           help='分析模式|Analysis mode')

    # 输出参数|Output parameters
    output_group = parser.add_argument_group('输出参数|Output parameters')
    output_group.add_argument('-o', '--output-dir', default='./telomere_output',
                             help='输出目录|Output directory')
    output_group.add_argument('-p', '--prefix', default='telomere',
                             help='输出前缀|Output prefix')

    # Explore模式参数|Explore mode parameters
    explore_group = parser.add_argument_group('Explore模式参数|Explore mode parameters')
    explore_group.add_argument('--explore-min', type=int, default=5,
                              help='最小重复长度|Minimum repeat length')
    explore_group.add_argument('--explore-max', type=int, default=12,
                              help='最大重复长度|Maximum repeat length')
    explore_group.add_argument('--explore-threshold', type=int, default=100,
                              help='重复阈值|Repeat threshold')
    explore_group.add_argument('--explore-distance', type=float, default=0.01,
                              help='距离比例|Distance ratio')

    # Find模式参数|Find mode parameters
    find_group = parser.add_argument_group('Find模式参数|Find mode parameters')
    find_group.add_argument('-c', '--clade',
                           help='分类群名称|Clade name')
    find_group.add_argument('-w', '--window', type=int, default=10000,
                           help='窗口大小|Window size')

    # Search模式参数|Search mode parameters
    search_group = parser.add_argument_group('Search模式参数|Search mode parameters')
    search_group.add_argument('-s', '--search-string',
                             help='搜索字符串|Search string')
    search_group.add_argument('--format', choices=['tsv', 'bedgraph'], default='tsv',
                             help='输出格式|Output format')

    # Plot模式参数|Plot mode parameters
    plot_group = parser.add_argument_group('Plot模式参数|Plot mode parameters')
    plot_group.add_argument('-t', '--tsv',
                           help='TSV文件路径|TSV file path')
    plot_group.add_argument('--plot-height', type=int, default=200,
                           help='图像高度|Plot height')
    plot_group.add_argument('--plot-width', type=int, default=1000,
                           help='图像宽度|Plot width')
    plot_group.add_argument('--plot-fontsize', type=int, default=12,
                           help='字体大小|Font size')
    plot_group.add_argument('--plot-strokewidth', type=int, default=2,
                           help='线条宽度|Stroke width')

    # 软件配置|Software configuration
    software_group = parser.add_argument_group('软件配置|Software configuration')
    software_group.add_argument('--tidk-path',
                               default='~/miniforge3/envs/tidk_v.0.2.65/bin/tidk',
                               help='tidk软件路径|tidk software path')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument('-v', '--verbose', action='store_true',
                          help='详细输出模式|Verbose output mode')
    log_group.add_argument('--log-file',
                          help='日志文件路径|Log file path')

    # 其他选项|Other options
    other_group = parser.add_argument_group('其他选项|Other options')
    other_group.add_argument('--print-clades', action='store_true',
                            help='打印支持的分类群列表|Print supported clades list')

    args = parser.parse_args()

    # 打印分类群|Print clades
    if args.print_clades:
        CladeDatabase.print_clade_table()
        sys.exit(0)

    # 构建配置字典|Build configuration dictionary
    config = {
        'genome_file': args.genome,
        'mode': args.mode,
        'output_dir': args.output_dir,
        'output_prefix': args.prefix,
        'tidk_path': args.tidk_path,
        'verbose': args.verbose,
        'log_file': args.log_file,
        'explore_min_length': args.explore_min,
        'explore_max_length': args.explore_max,
        'explore_threshold': args.explore_threshold,
        'explore_distance': args.explore_distance,
        'clade': args.clade,
        'window_size': args.window,
        'search_string': args.search_string,
        'output_format': args.format,
        'tsv_file': args.tsv,
        'plot_height': args.plot_height,
        'plot_width': args.plot_width,
        'plot_fontsize': args.plot_fontsize,
        'plot_strokewidth': args.plot_strokewidth,
    }

    # 创建分析器并运行|Create analyzer and run
    try:
        finder = TelomereFinder(**config)
        success = finder.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"错误|Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
