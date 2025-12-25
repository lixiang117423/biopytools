"""
SNP Index计算核心模块 | SNP Index Calculation Core Module
"""

import csv
import logging
from typing import List, Tuple, Optional
from .config import SNPIndexConfig
from .utils import (
    open_file, parse_vcf_line, extract_ad_values, calculate_snp_index,
    parse_quality_filters, get_sample_names_from_vcf, validate_vcf_file,
    format_number, ensure_dir_exists
)


class SNPIndexCalculator:
    """SNP Index计算器 | SNP Index Calculator"""

    def __init__(self, config: SNPIndexConfig):
        """
        初始化计算器 | Initialize calculator

        Args:
            config: 配置对象 | Configuration object
        """
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self.sample_names = []
        self.stats = {
            'total_variants': 0,
            'filtered_variants': 0,
            'processed_variants': 0,
            'low_depth_variants': 0,
            'low_quality_variants': 0
        }

    def validate_input(self) -> bool:
        """
        验证输入文件和参数 | Validate input files and parameters

        Returns:
            bool: 验证是否通过 | Whether validation passed
        """
        self.logger.info("验证输入文件和参数 | Validating input files and parameters")

        # 验证VCF文件 | Validate VCF file
        if not validate_vcf_file(self.config.input_vcf, min_samples=2):
            self.logger.error(f"VCF文件格式错误或样本数不足2个 | VCF file format error or insufficient samples: {self.config.input_vcf}")
            return False

        # 获取样本名称 | Get sample names
        self.sample_names = get_sample_names_from_vcf(self.config.input_vcf)
        if len(self.sample_names) < 2:
            self.logger.error(f"VCF文件样本数不足，需要至少2个样本，实际{len(self.sample_names)}个 | Insufficient samples in VCF file, need at least 2, found {len(self.sample_names)}")
            return False

        # 使用指定的样本名称或前两个样本 | Use specified sample names or first two samples
        if self.config.sample_names:
            # 验证指定的样本名是否存在 | Validate specified sample names exist
            for name in self.config.sample_names:
                if name not in self.sample_names:
                    self.logger.error(f"指定的样本名在VCF文件中不存在 | Specified sample name not found in VCF file: {name}")
                    return False
            self.sample_names = self.config.sample_names[:2]  # 只使用前两个 | Use only first two
        else:
            self.sample_names = self.sample_names[:2]  # 使用前两个样本 | Use first two samples

        self.logger.info(f"使用样本进行计算 | Using samples for calculation: {self.sample_names[0]} and {self.sample_names[1]}")

        # 创建输出目录 | Create output directory
        ensure_dir_exists(self.config.output_dir)

        return True

    def calculate(self) -> bool:
        """
        执行SNP index计算 | Perform SNP index calculation

        Returns:
            bool: 计算是否成功 | Whether calculation succeeded
        """
        self.logger.info("=" * 60)
        self.logger.info("开始SNP index计算 | Starting SNP index calculation")
        self.logger.info("=" * 60)
        self.logger.info(f"输入VCF文件 | Input VCF file: {self.config.input_vcf}")
        self.logger.info(f"输出文件 | Output file: {self.config.output_file}")
        self.logger.info(f"最小深度 | Minimum depth: {self.config.min_depth}")
        self.logger.info(f"最小质量 | Minimum quality: {self.config.min_quality}")

        # 验证输入 | Validate input
        if not self.validate_input():
            return False

        try:
            # 处理VCF文件 | Process VCF file
            success = self._process_vcf_file()
            if success:
                self._log_statistics()
                self.logger.info("SNP index计算完成 | SNP index calculation completed")
            else:
                self.logger.error("SNP index计算失败 | SNP index calculation failed")

            return success

        except Exception as e:
            self.logger.error(f"计算过程中发生错误 | Error occurred during calculation: {str(e)}", exc_info=True)
            return False

    def _process_vcf_file(self) -> bool:
        """
        处理VCF文件主逻辑 | Main logic for processing VCF file

        Returns:
            bool: 处理是否成功 | Whether processing succeeded
        """
        sample1_name = self.sample_names[0]
        sample2_name = self.sample_names[1]

        try:
            with open_file(self.config.input_vcf, 'r') as f, \
                    open(self.config.output_file, 'w', newline='', encoding='utf-8') as out_f:

                writer = csv.writer(out_f, delimiter='\t')

                # 写入表头 | Write header
                header = self._generate_header(sample1_name, sample2_name)
                writer.writerow(header)

                # 跳过注释行，找到数据开始 | Skip comment lines, find data start
                for line in f:
                    if line.startswith('#CHROM'):
                        break

                # 处理数据行 | Process data lines
                for line in f:
                    self.stats['total_variants'] += 1
                    self._process_vcf_line(line, writer, sample1_name, sample2_name)

                    # 进度报告 | Progress report
                    if self.stats['processed_variants'] > 0 and self.stats['processed_variants'] % 10000 == 0:
                        self.logger.info(f"已处理 | Processed: {self.stats['processed_variants']:,} 个变异位点")

            return True

        except Exception as e:
            self.logger.error(f"处理VCF文件时出错 | Error processing VCF file: {str(e)}")
            return False

    def _generate_header(self, sample1_name: str, sample2_name: str) -> List[str]:
        """
        生成输出文件表头 | Generate output file header

        Args:
            sample1_name: 样本1名称 | Sample 1 name
            sample2_name: 样本2名称 | Sample 2 name

        Returns:
            list: 表头列表 | Header list
        """
        return [
            'Chromosome', 'Position', 'Reference', 'Alternative',
            f'{sample1_name}_Ref_Depth', f'{sample1_name}_Alt_Depth', f'{sample1_name}_SNP_index',
            f'{sample2_name}_Ref_Depth', f'{sample2_name}_Alt_Depth', f'{sample2_name}_SNP_index',
            'Delta_SNP_index'
        ]

    def _process_vcf_line(self, line: str, writer: csv.writer, sample1_name: str, sample2_name: str) -> None:
        """
        处理单个VCF数据行 | Process single VCF data line

        Args:
            line: VCF行 | VCF line
            writer: CSV写入器 | CSV writer
            sample1_name: 样本1名称 | Sample 1 name
            sample2_name: 样本2名称 | Sample 2 name
        """
        try:
            # 解析VCF行 | Parse VCF line
            chrom, pos, ref, alt, samples = parse_vcf_line(line)

            # 获取样本索引 | Get sample indices
            sample1_idx = self.sample_names.index(sample1_name)
            sample2_idx = self.sample_names.index(sample2_name)

            if len(samples) <= max(sample1_idx, sample2_idx):
                return  # 样本数不足 | Insufficient samples

            # 提取质量信息 | Extract quality information
            qual, mq = parse_quality_filters(line)

            # 质量过滤 | Quality filtering
            if qual is not None and qual < self.config.min_quality:
                self.stats['low_quality_variants'] += 1
                return
            if mq is not None and mq < self.config.min_mapping_quality:
                self.stats['low_quality_variants'] += 1
                return

            # 提取AD值 | Extract AD values
            sample1_ref, sample1_alt = extract_ad_values(samples[sample1_idx])
            sample2_ref, sample2_alt = extract_ad_values(samples[sample2_idx])

            # 深度过滤 | Depth filtering
            sample1_total = sample1_ref + sample1_alt
            sample2_total = sample2_ref + sample2_alt

            if sample1_total < self.config.min_depth or sample2_total < self.config.min_depth:
                self.stats['low_depth_variants'] += 1
                return

            # 计算SNP index | Calculate SNP index
            snp_index1 = calculate_snp_index(sample1_ref, sample1_alt)
            snp_index2 = calculate_snp_index(sample2_ref, sample2_alt)
            delta_snp_index = snp_index1 - snp_index2

            # 写入结果 | Write results
            row = [
                chrom, pos, ref, alt,
                sample1_ref, sample1_alt, format_number(snp_index1),
                sample2_ref, sample2_alt, format_number(snp_index2),
                format_number(delta_snp_index)
            ]
            writer.writerow(row)

            self.stats['processed_variants'] += 1

        except Exception as e:
            self.logger.warning(f"处理VCF行时出错，跳过此行 | Error processing VCF line, skipping: {str(e)}")
            self.stats['filtered_variants'] += 1

    def _log_statistics(self) -> None:
        """记录统计信息 | Log statistics"""
        self.logger.info("=" * 60)
        self.logger.info("计算统计信息 | Calculation Statistics")
        self.logger.info("=" * 60)
        self.logger.info(f"总变异位点 | Total variants: {self.stats['total_variants']:,}")
        self.logger.info(f"已处理变异位点 | Processed variants: {self.stats['processed_variants']:,}")
        self.logger.info(f"低质量位点 | Low quality variants: {self.stats['low_quality_variants']:,}")
        self.logger.info(f"低深度位点 | Low depth variants: {self.stats['low_depth_variants']:,}")
        self.logger.info(f"过滤位点 | Filtered variants: {self.stats['filtered_variants']:,}")

        if self.stats['total_variants'] > 0:
            process_rate = (self.stats['processed_variants'] / self.stats['total_variants']) * 100
            self.logger.info(f"处理成功率 | Success rate: {process_rate:.1f}%")

    def get_statistics(self) -> dict:
        """
        获取统计信息 | Get statistics

        Returns:
            dict: 统计信息字典 | Statistics dictionary
        """
        return self.stats.copy()