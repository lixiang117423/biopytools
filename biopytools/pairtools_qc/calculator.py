"""
Hi-C数据质量控制评估核心计算模块|Hi-C QC Assessment Core Calculation Module
"""

import os
import subprocess
from typing import Dict, Any, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class PairtoolsQCCalculator:
    """Hi-C数据质量控制评估计算器|Hi-C QC Assessment Calculator"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: PairtoolsQCConfig配置对象|PairtoolsQCConfig object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger
        self.stats_data = {}
        self.bam_converted = False

    def _is_bam_file(self) -> bool:
        """检查是否为BAM/SAM文件|Check if input is BAM/SAM file"""
        file_path = self.config.pairs_file.lower()
        return file_path.endswith('.bam') or file_path.endswith('.sam')

    def _convert_bam_to_pairs(self) -> str:
        """将BAM文件转换为pairs文件|Convert BAM file to pairs file

        Returns:
            str: 转换后的pairs文件路径|Converted pairs file path
        """
        self.logger.info("检测到BAM/SAM输入文件，开始转换为pairs格式|Detected BAM/SAM input, converting to pairs format")

        # 检查chrom.sizes文件|Check chrom.sizes file
        if not hasattr(self.config, 'chroms_path') or not self.config.chroms_path:
            raise ValueError(
                "处理BAM文件需要提供chrom.sizes文件|Processing BAM files requires chrom.sizes file\n"
                "请使用 --chroms-path 参数指定|Please specify with --chroms-path parameter\n"
                "或提供fasta文件自动生成|Or provide a fasta file to auto-generate"
            )

        if not os.path.exists(self.config.chroms_path):
            raise ValueError(f"Chrom.sizes文件不存在|Chrom.sizes file not found: {self.config.chroms_path}")

        # 生成输出文件名|Generate output filename
        base_name = os.path.basename(self.config.pairs_file).replace('.bam', '').replace('.sam', '')
        pairs_file = os.path.join(self.config.output_dir, f"{base_name}.pairs.gz")

        # 构建pairtools parse命令|Build pairtools parse command
        pairtools_cmd = os.path.basename(self.config.pairtools_path)
        cmd_args = [
            'parse',
            '-c', self.config.chroms_path,
            '-o', pairs_file,
            '--drop-sam',  # 去掉SAM信息以减少文件大小|Drop SAM info to reduce file size
            self.config.pairs_file
        ]

        # 自动包装conda环境的命令|Auto-wrap conda environment commands
        wrapped_cmd = build_conda_command(pairtools_cmd, cmd_args)

        self.logger.info(f"运行pairtools parse|Running pairtools parse")
        self.logger.debug(f"命令|Command: {' '.join(wrapped_cmd)}")

        try:
            # 执行转换|Execute conversion
            result = subprocess.run(
                wrapped_cmd,
                check=True,
                capture_output=True,
                text=True
            )

            self.logger.info(f"BAM转换完成|BAM conversion completed: {pairs_file}")
            self.bam_converted = True
            return pairs_file

        except subprocess.CalledProcessError as e:
            self.logger.error(f"pairtools parse执行失败|pairtools parse execution failed: {e.stderr}")
            raise

    def _get_input_file(self) -> str:
        """获取输入文件（如果需要则转换BAM）|Get input file (convert BAM if needed)

        Returns:
            str: 输入文件路径|Input file path
        """
        if self._is_bam_file():
            # 需要先转换BAM|Need to convert BAM first
            return self._convert_bam_to_pairs()
        else:
            # 直接使用pairs文件|Use pairs file directly
            return self.config.pairs_file

    def run_pairtools_stats(self, input_file: str = None) -> str:
        """运行pairtools stats命令|Run pairtools stats command

        Args:
            input_file: 输入文件路径（如果为None则自动获取）|Input file path (auto-get if None)

        Returns:
            str: 统计文件路径|Statistics file path
        """
        if input_file is None:
            input_file = self._get_input_file()

        self.logger.info("运行pairtools stats|Running pairtools stats")

        # 生成统计文件路径|Generate stats file path
        base_name = os.path.basename(input_file).replace('.gz', '').replace('.pairs', '')
        stats_file = os.path.join(self.config.output_dir, f"{base_name}_stats.tsv")

        # 构建命令|Build command
        pairtools_cmd = os.path.basename(self.config.pairtools_path)
        cmd_args = [
            'stats',
            '-o', stats_file,
            input_file
        ]

        # 自动包装conda环境的命令|Auto-wrap conda environment commands
        wrapped_cmd = build_conda_command(pairtools_cmd, cmd_args)

        self.logger.debug(f"命令|Command: {' '.join(wrapped_cmd)}")

        try:
            # 执行命令|Execute command
            result = subprocess.run(
                wrapped_cmd,
                check=True,
                capture_output=True,
                text=True
            )

            self.logger.info(f"统计文件已生成|Statistics file generated: {stats_file}")
            return stats_file

        except subprocess.CalledProcessError as e:
            self.logger.error(f"pairtools stats执行失败|pairtools stats execution failed: {e.stderr}")
            raise

    def parse_stats_file(self, stats_file: str) -> Dict[str, Any]:
        """解析pairtools stats输出文件|Parse pairtools stats output file

        Args:
            stats_file: 统计文件路径|Statistics file path

        Returns:
            Dict[str, Any]: 解析后的统计数据|Parsed statistics data
        """
        self.logger.info(f"解析统计文件|Parsing statistics file: {stats_file}")

        stats = {}
        summary = {}

        try:
            with open(stats_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    # 分割键值对|Split key-value pair
                    parts = line.split('\t')
                    if len(parts) != 2:
                        continue

                    key = parts[0]
                    value = parts[1]

                    # 尝试转换为数字|Try to convert to number
                    try:
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except ValueError:
                        pass  # 保持为字符串|Keep as string

                    # 提取关键统计指标|Extract key statistics
                    if key in ['total_unmapped', 'total_single_sided_mapped', 'total_mapped',
                               'total_dups', 'cis', 'trans']:
                        stats[key] = value

            # 计算百分比和比例|Calculate percentages and ratios
            total_pairs = sum([
                stats.get('total_unmapped', 0),
                stats.get('total_single_sided_mapped', 0),
                stats.get('total_mapped', 0)
            ])

            if total_pairs > 0:
                summary['total_pairs'] = total_pairs
                summary['total_unmapped'] = stats.get('total_unmapped', 0)
                summary['total_single_sided_mapped'] = stats.get('total_single_sided_mapped', 0)
                summary['total_mapped'] = stats.get('total_mapped', 0)
                summary['total_dups'] = stats.get('total_dups', 0)
                summary['cis'] = stats.get('cis', 0)
                summary['trans'] = stats.get('trans', 0)

                # 计算百分比|Calculate percentages
                summary['unmapped_rate'] = (summary['total_unmapped'] / total_pairs) * 100
                summary['single_sided_rate'] = (summary['total_single_sided_mapped'] / total_pairs) * 100
                summary['mapped_rate'] = (summary['total_mapped'] / total_pairs) * 100

                if summary['total_mapped'] > 0:
                    summary['dup_rate'] = (summary['total_dups'] / summary['total_mapped']) * 100
                else:
                    summary['dup_rate'] = 0.0

                # 计算cis/trans比例|Calculate cis/trans ratio
                if summary['trans'] > 0:
                    summary['cis_trans_ratio'] = summary['cis'] / summary['trans']
                else:
                    summary['cis_trans_ratio'] = float('inf') if summary['cis'] > 0 else 0.0

            self.stats_data = summary
            return summary

        except FileNotFoundError:
            self.logger.error(f"找不到统计文件|Statistics file not found: {stats_file}")
            raise
        except Exception as e:
            self.logger.error(f"解析统计文件时出错|Error parsing statistics file: {str(e)}")
            raise

    def assess_quality(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """评估数据质量|Assess data quality

        Args:
            stats: 统计数据|Statistics data

        Returns:
            Dict[str, Any]: 评估结果|Assessment results
        """
        self.logger.info("评估数据质量|Assessing data quality")

        thresholds = self.config.get_thresholds()
        results = {
            'stats': stats,
            'thresholds': thresholds,
            'assessment': {},
            'passed': True
        }

        # 评估各项指标|Assess each metric
        metrics = [
            ('unmapped_rate', stats.get('unmapped_rate', 0), thresholds['max_unmapped_rate'], 'le'),
            ('single_sided_rate', stats.get('single_sided_rate', 0), thresholds['max_single_sided_rate'], 'le'),
            ('mapped_rate', stats.get('mapped_rate', 0), thresholds['min_mapped_rate'], 'ge'),
            ('dup_rate', stats.get('dup_rate', 0), thresholds['max_dup_rate'], 'le'),
            ('cis_trans_ratio', stats.get('cis_trans_ratio', 0), thresholds['min_cis_trans_ratio'], 'ge')
        ]

        for metric_name, value, threshold, comparison in metrics:
            if comparison == 'le':
                passed = value <= threshold
                desc = f"{value:.2f}% (阈值|threshold <= {threshold}%)"
            else:  # ge
                passed = value >= threshold
                desc = f"{value:.2f}% (阈值|threshold >= {threshold}%)" if metric_name != 'cis_trans_ratio' else f"{value:.2f} (阈值|threshold >= {threshold})"

            results['assessment'][metric_name] = {
                'value': value,
                'threshold': threshold,
                'passed': passed,
                'description': desc
            }

            if not passed:
                results['passed'] = False

        return results

    def calculate(self) -> Dict[str, Any]:
        """执行完整的质量评估流程|Execute complete quality assessment pipeline

        Returns:
            Dict[str, Any]: 评估结果|Assessment results
        """
        # 获取输入文件（如果需要则转换BAM）|Get input file (convert BAM if needed)
        input_file = self._get_input_file()

        # 运行pairtools stats|Run pairtools stats
        stats_file = self.run_pairtools_stats(input_file)

        # 解析统计文件|Parse statistics file
        stats = self.parse_stats_file(stats_file)

        # 评估质量|Assess quality
        results = self.assess_quality(stats)

        return results

    def format_number(self, num: int) -> str:
        """格式化数字，添加千位分隔符|Format number with thousands separator

        Args:
            num: 要格式化的数字|Number to format

        Returns:
            str: 格式化后的字符串|Formatted string
        """
        return f"{num:,}"
