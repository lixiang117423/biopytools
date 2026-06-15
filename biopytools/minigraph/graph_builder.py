"""
Minigraph核心功能模块|Minigraph Core Functions Module
"""

import os
from typing import List, Optional
from .config import (MinigraphBuildConfig, MinigraphCallConfig,
                     MinigraphBubbleConfig, MinigraphMapConfig)
from .utils import build_conda_command, CommandRunner


class MinigraphGraphBuilder:
    """Minigraph泛基因组图构建器|Minigraph Pangenome Graph Builder"""

    def __init__(self, config: MinigraphBuildConfig, logger, cmd_runner: CommandRunner):
        """
        初始化图构建器|Initialize graph builder

        Args:
            config: 构建配置对象|Build configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_graph(self) -> bool:
        """
        构建泛基因组图|Build pangenome graph

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始构建泛基因组图|Starting pangenome graph construction")

        # 构建minigraph命令|Build minigraph command
        cmd = self._build_command()

        # 执行命令，输出到文件|Execute command with output to file
        success = self.cmd_runner.run_to_file(
            cmd,
            output_file=self.config.output_gfa,
            description="构建泛基因组图|Build pangenome graph"
        )

        if not success:
            self.logger.error("泛基因组图构建失败|Pangenome graph construction failed")
            return False

        # 验证输出文件|Verify output file
        if not os.path.exists(self.config.output_gfa):
            self.logger.error(f"输出GFA文件未生成|Output GFA file not generated: {self.config.output_gfa}")
            return False

        # 统计图信息|Get graph statistics
        self._log_graph_stats()

        self.logger.info("泛基因组图构建成功|Pangenome graph constructed successfully")
        return True

    def _build_command(self) -> List[str]:
        """
        构建minigraph命令|Build minigraph command

        Returns:
            命令列表|Command list
        """
        # 基础命令|Base command
        minigraph_cmd = build_conda_command(self.config.minigraph_path, [])
        # minigraph_cmd 是完整的 conda run 命令或直接命令
        # 例如: ['conda', 'run', '-n', 'swave_v.1.2', 'minigraph']
        # 或: ['minigraph']

        # 添加预设参数|Add preset parameters
        cmd = minigraph_cmd + [f'-{self.config.preset}']

        # 添加其他参数|Add other parameters
        if self.config.min_identity != 0.9:
            cmd.extend(['-i', str(self.config.min_identity)])

        if self.config.min_aln_len != 100000:
            cmd.extend(['-l', str(self.config.min_aln_len)])

        if self.config.max_gap != 1000000:
            cmd.extend(['-g', str(self.config.max_gap)])

        # 性能参数|Performance parameters
        cmd.extend(['-t', str(self.config.threads)])

        if self.config.batch_size:
            cmd.extend(['-K', f'{self.config.batch_size}M'])

        # 如果是追加模式，添加现有GFA作为第一个输入|If append mode, add existing GFA as first input
        input_files = [self.config.ref_fasta] + self.config.sample_fastas

        # 添加输入文件|Add input files
        cmd.extend(input_files)

        return cmd

    def _log_graph_stats(self):
        """记录图统计信息|Log graph statistics"""
        try:
            with open(self.config.output_gfa, 'r') as f:
                lines = f.readlines()

            # 统计不同行类型|Count different line types
            segment_count = sum(1 for line in lines if line.startswith('S'))
            link_count = sum(1 for line in lines if line.startswith('L'))
            path_count = sum(1 for line in lines if line.startswith('P'))

            self.logger.info(f"图统计|Graph statistics:")
            self.logger.info(f"  片段数|Segments: {segment_count}")
            self.logger.info(f"  连接数|Links: {link_count}")
            self.logger.info(f"  路径数|Paths: {path_count}")

        except Exception as e:
            self.logger.warning(f"无法获取图统计信息|Cannot get graph statistics: {e}")


class MinigraphSVCaller:
    """Minigraph SV调用器|Minigraph SV Caller"""

    def __init__(self, config: MinigraphCallConfig, logger, cmd_runner: CommandRunner):
        """
        初始化SV调用器|Initialize SV caller

        Args:
            config: 调用配置对象|Call configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def call_samples(self) -> bool:
        """
        为所有样本调用SV|Call SVs for all samples

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始SV调用|Starting SV calling")

        success_count = 0

        for i, sample_fasta in enumerate(self.config.sample_fastas, 1):
            self.logger.info(f"处理样本 {i}/{len(self.config.sample_fastas)}|Processing sample {i}/{len(self.config.sample_fastas)}: {sample_fasta}")

            if self._call_single_sample(sample_fasta):
                success_count += 1
            else:
                self.logger.warning(f"样本 {i} 调用失败|Sample {i} calling failed: {sample_fasta}")

        total = len(self.config.sample_fastas)
        self.logger.info(f"SV调用完成|SV calling completed: {success_count}/{total} 成功|succeeded")

        return success_count == total

    def _call_single_sample(self, sample_fasta: str) -> bool:
        """
        为单个样本调用SV|Call SV for a single sample

        Args:
            sample_fasta: 样本FASTA文件|Sample FASTA file

        Returns:
            是否成功|Whether successful
        """
        # 获取样本名|Get sample name
        sample_name = os.path.splitext(os.path.basename(sample_fasta))[0]
        output_bed = os.path.join(self.config.output_dir, f"{sample_name}.bed")

        # 构建命令|Build command
        minigraph_cmd = build_conda_command(self.config.minigraph_path, [])
        cmd = minigraph_cmd + [f'-x{self.config.preset}']

        if self.config.call_mode:
            cmd.append('--call')

        cmd.extend(['-t', str(self.config.threads)])
        cmd.extend([self.config.graph_gfa, sample_fasta])

        # 执行命令，输出到文件|Execute command with output to file
        success = self.cmd_runner.run_to_file(
            cmd,
            output_file=output_bed,
            description=f"SV调用|SV calling: {sample_name}"
        )

        if success and os.path.exists(output_bed):
            self.logger.info(f"输出文件|Output: {output_bed}")
            return True

        return False


class MinigraphBubbleExtractor:
    """Minigraph bubble提取器|Minigraph Bubble Extractor"""

    def __init__(self, config: MinigraphBubbleConfig, logger, cmd_runner: CommandRunner):
        """
        初始化bubble提取器|Initialize bubble extractor

        Args:
            config: 提取配置对象|Extraction configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def extract_bubbles(self) -> bool:
        """
        提取SV bubbles|Extract SV bubbles

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始提取SV bubbles|Starting SV bubble extraction")

        # 构建gfatools命令|Build gfatools command
        gfatools_cmd = build_conda_command(self.config.gfatools_path, [])
        cmd = gfatools_cmd + ['bubble', self.config.graph_gfa]

        # 执行命令，输出到文件|Execute command with output to file
        success = self.cmd_runner.run_to_file(
            cmd,
            output_file=self.config.output_bed,
            description="提取SV bubbles|Extract SV bubbles"
        )

        if not success:
            self.logger.error("SV bubble提取失败|SV bubble extraction failed")
            return False

        # 验证输出文件|Verify output file
        if not os.path.exists(self.config.output_bed):
            self.logger.error(f"输出BED文件未生成|Output BED file not generated: {self.config.output_bed}")
            return False

        # 统计bubble数量|Count bubbles
        self._log_bubble_stats()

        self.logger.info("SV bubble提取成功|SV bubble extraction completed successfully")
        return True

    def _log_bubble_stats(self):
        """记录bubble统计信息|Log bubble statistics"""
        try:
            with open(self.config.output_bed, 'r') as f:
                line_count = sum(1 for line in f if line.strip())

            self.logger.info(f"Bubble统计|Bubble statistics:")
            self.logger.info(f"  总数|Total: {line_count}")

        except Exception as e:
            self.logger.warning(f"无法获取bubble统计信息|Cannot get bubble statistics: {e}")


class MinigraphMapper:
    """Minigraph序列映射器|Minigraph Sequence Mapper"""

    def __init__(self, config: MinigraphMapConfig, logger, cmd_runner: CommandRunner):
        """
        初始化映射器|Initialize mapper

        Args:
            config: 映射配置对象|Mapping configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def map_sequences(self) -> bool:
        """
        映射序列到图|Map sequences to graph

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始序列映射|Starting sequence mapping")

        # 构建minigraph命令|Build minigraph command
        minigraph_cmd = build_conda_command(self.config.minigraph_path, [])
        cmd = minigraph_cmd + [f'-x{self.config.preset}']

        # 添加参数|Add parameters
        if self.config.max_intron_len:
            cmd.extend(['-G', str(self.config.max_intron_len)])

        cmd.extend(['-t', str(self.config.threads)])

        if self.config.batch_size:
            cmd.extend(['-K', f'{self.config.batch_size}M'])

        # 输入和输出|Input and output
        cmd.extend([self.config.graph_gfa])
        cmd.extend(self.config.query_fastas)

        # 执行命令，输出到文件|Execute command with output to file
        success = self.cmd_runner.run_to_file(
            cmd,
            output_file=self.config.output_gaf,
            description="序列映射|Sequence mapping"
        )

        if not success:
            self.logger.error("序列映射失败|Sequence mapping failed")
            return False

        # 验证输出文件|Verify output file
        if not os.path.exists(self.config.output_gaf):
            self.logger.error(f"输出GAF文件未生成|Output GAF file not generated: {self.config.output_gaf}")
            return False

        # 统计映射数量|Count mappings
        self._log_mapping_stats()

        self.logger.info("序列映射成功|Sequence mapping completed successfully")
        return True

    def _log_mapping_stats(self):
        """记录映射统计信息|Log mapping statistics"""
        try:
            with open(self.config.output_gaf, 'r') as f:
                line_count = sum(1 for line in f if line.strip())

            self.logger.info(f"映射统计|Mapping statistics:")
            self.logger.info(f"  映射数|Mappings: {line_count}")

        except Exception as e:
            self.logger.warning(f"无法获取映射统计信息|Cannot get mapping statistics: {e}")
