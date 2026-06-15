"""
FASTQ文件统计计算器模块|FASTQ File Statistics Calculator Module
"""

import csv
import os
import subprocess
import tempfile
from typing import Dict, List, Optional
from .utils import format_number, format_size, get_file_size


class FastqStatsCalculator:
    """FASTQ文件统计计算器|FASTQ File Statistics Calculator"""

    def __init__(self, config, logger):
        """
        初始化计算器|Initialize calculator

        Args:
            config: FastqStatsConfig配置对象|FastqStatsConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def run_seqkit_stats(self, files: List[str]) -> Optional[Dict[str, Dict]]:
        """
        使用seqkit stats处理文件列表|Process files using seqkit stats

        Args:
            files: 文件路径列表|List of file paths

        Returns:
            统计信息字典|Statistics dictionary
        """
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as tmp:
            tmp_file = tmp.name

        try:
            # 构建seqkit命令|Build seqkit command
            cmd = ['seqkit', 'stats'] + files + ['-T', '-j', str(self.config.threads)]

            self.logger.info(f"运行seqkit命令|Running seqkit command with {len(files)} files")

            # 运行seqkit|Run seqkit
            with open(tmp_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                self.logger.error(f"seqkit运行失败|seqkit execution failed: {result.stderr}")
                return None

            # 读取结果|Read results
            stats_dict = {}
            with open(tmp_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    stats_dict[row['file']] = {
                        'format': row['format'],
                        'type': row['type'],
                        'num_seqs': int(row['num_seqs'].replace(',', '')),
                        'sum_len': int(row['sum_len'].replace(',', '')),
                        'min_len': int(row['min_len'].replace(',', '')),
                        'avg_len': float(row['avg_len'].replace(',', '')),
                        'max_len': int(row['max_len'].replace(',', ''))
                    }

            self.logger.info(f"seqkit成功处理{len(stats_dict)}个文件|seqkit successfully processed {len(stats_dict)} files")
            return stats_dict

        finally:
            # 清理临时文件|Clean up temp file
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

    def process_samples(self, samples: Dict[str, Dict[str, str]]) -> List[Dict]:
        """
        处理所有样品|Process all samples

        Args:
            samples: 样品字典|Sample dictionary {sample_name: {'R1': r1_path, 'R2': r2_path}}

        Returns:
            处理结果列表|List of processing results
        """
        # 收集所有需要处理的文件|Collect all files to process
        all_files = []
        file_to_sample = {}
        seen_files = set()  # 避免重复添加同一个文件|Avoid adding same file twice

        for sample, files_dict in samples.items():
            for read_type, filepath in files_dict.items():
                # 如果R1和R2是同一个文件，只处理一次
                # If R1 and R2 are the same file, only process once
                if filepath not in seen_files:
                    seen_files.add(filepath)
                    all_files.append(filepath)
                    file_to_sample[filepath] = (sample, read_type)
                else:
                    # 如果文件已经存在，检查是否需要更新read_type
                    # 优先使用R1而不是R2（单端数据可能被错误识别为双端）
                    # If file already exists, check if we need to update read_type
                    # Prefer R1 over R2 (single-end data might be misidentified as paired-end)
                    current_sample, current_read_type = file_to_sample[filepath]
                    if current_read_type == 'R2' and read_type == 'R1':
                        file_to_sample[filepath] = (sample, read_type)

        self.logger.info(f"正在处理{len(all_files)}个文件|Processing {len(all_files)} files")

        # 运行seqkit stats|Run seqkit stats
        stats_dict = self.run_seqkit_stats(all_files)

        if not stats_dict:
            self.logger.error("seqkit处理失败|seqkit processing failed")
            return []

        # 组织结果|Organize results
        results = {}
        for filepath, stats in stats_dict.items():
            sample, read_type = file_to_sample[filepath]

            if sample not in results:
                results[sample] = {
                    'sample_name': sample,
                    'R1_file': 'N/A',
                    'R2_file': 'N/A'
                }

            # 填充文件名|Fill file names
            results[sample][f'{read_type}_file'] = os.path.basename(filepath)

            # 填充统计信息|Fill statistics
            results[sample][f'{read_type}_reads'] = stats['num_seqs']
            results[sample][f'{read_type}_min_length'] = stats['min_len']
            results[sample][f'{read_type}_max_length'] = stats['max_len']
            results[sample][f'{read_type}_mean_length'] = round(stats['avg_len'], 2)
            results[sample][f'{read_type}_total_bases'] = stats['sum_len']

            # 获取文件大小|Get file size
            file_size, uncompressed_size = get_file_size(filepath)
            results[sample][f'{read_type}_file_size'] = file_size
            results[sample][f'{read_type}_file_size_formatted'] = format_size(file_size)
            if uncompressed_size is not None:
                results[sample][f'{read_type}_uncompressed_size'] = uncompressed_size
                results[sample][f'{read_type}_uncompressed_size_formatted'] = format_size(uncompressed_size)
            else:
                results[sample][f'{read_type}_uncompressed_size'] = 'N/A'
                results[sample][f'{read_type}_uncompressed_size_formatted'] = 'N/A'

        # 计算双末端总和|Calculate paired-end totals
        for sample, data in results.items():
            # 检查是否有真正的双末端数据（R1和R2是不同的文件）
            # Check if there's real paired-end data (R1 and R2 are different files)
            has_r1 = 'R1_reads' in data and data['R1_reads'] != 'N/A'
            has_r2 = 'R2_reads' in data and data['R2_reads'] != 'N/A'

            # 检查R1和R2文件是否不同
            # Check if R1 and R2 files are different
            r1_file = data.get('R1_file', '')
            r2_file = data.get('R2_file', '')
            is_different_pair = r1_file != r2_file

            if has_r1 and has_r2 and is_different_pair:
                # 真正的双末端数据|Real paired-end data
                data['total_reads'] = data['R1_reads'] + data['R2_reads']
                data['total_bases'] = data['R1_total_bases'] + data['R2_total_bases']
            elif has_r1:
                # 只有R1|Only R1
                data['total_reads'] = data['R1_reads']
                data['total_bases'] = data['R1_total_bases']
            elif has_r2:
                # 只有R2（单端数据被标记为R2）|Only R2 (single-end data marked as R2)
                data['total_reads'] = data['R2_reads']
                data['total_bases'] = data['R2_total_bases']
            else:
                data['total_reads'] = 0
                data['total_bases'] = 0

        return list(results.values())

    def write_csv(self, results: List[Dict], output_file: str):
        """
        写入CSV文件|Write to CSV file

        Args:
            results: 结果列表|List of results
            output_file: 输出文件路径|Output file path
        """
        if not results:
            self.logger.error("没有结果可写入|No results to write")
            return

        # 确定所有可能的字段|Determine all possible fields
        fieldnames = [
            'sample_name', 'R1_file', 'R2_file',
            'R1_reads', 'R1_min_length', 'R1_max_length',
            'R1_mean_length', 'R1_total_bases',
            'R1_file_size_formatted', 'R1_uncompressed_size_formatted'
        ]

        # 检查是否有R2数据|Check if R2 data exists
        if any('R2_reads' in r and r['R2_reads'] != 'N/A' for r in results):
            fieldnames.extend([
                'R2_reads', 'R2_min_length', 'R2_max_length',
                'R2_mean_length', 'R2_total_bases',
                'R2_file_size_formatted', 'R2_uncompressed_size_formatted',
                'total_reads', 'total_bases'
            ])

        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)

        self.logger.info(f"结果已保存到|Results saved to: {output_file}")

    def write_excel(self, results: List[Dict], output_file: str):
        """
        写入Excel文件|Write to Excel file

        Args:
            results: 结果列表|List of results
            output_file: 输出文件路径|Output file path
        """
        try:
            import pandas as pd
        except ImportError:
            self.logger.error("需要安装pandas和openpyxl库|pandas and openpyxl libraries are required")
            self.logger.error("请运行|Please run: pip install pandas openpyxl")
            raise

        df = pd.DataFrame(results)

        # 调整列顺序|Adjust column order
        cols = ['sample_name', 'R1_file', 'R2_file']
        r1_cols = [c for c in df.columns if c.startswith('R1_')]
        r2_cols = [c for c in df.columns if c.startswith('R2_')]
        total_cols = [c for c in df.columns if c.startswith('total_')]

        ordered_cols = cols + r1_cols + r2_cols + total_cols
        ordered_cols = [c for c in ordered_cols if c in df.columns]

        df = df[ordered_cols]

        df.to_excel(output_file, index=False, engine='openpyxl')
        self.logger.info(f"结果已保存到|Results saved to: {output_file}")

    def write_results(self, results: List[Dict]):
        """
        根据配置写入结果|Write results based on configuration

        Args:
            results: 结果列表|List of results
        """
        if self.config.output_format == 'excel':
            self.write_excel(results, self.config.output_file)
        else:
            self.write_csv(results, self.config.output_file)

    def print_summary(self, results: List[Dict]):
        """
        打印结果摘要|Print results summary

        Args:
            results: 结果列表|List of results
        """
        self.logger.info(f"\n统计摘要|Statistics Summary:")
        self.logger.info(f"{'样品名称|Sample Name':<20} {'R1 reads':<12} {'R2 reads':<12} {'总reads|Total':<12} {'R1大小|R1 Size':<12} {'R1解压|R1 Uncomp':<12}")
        self.logger.info(f"{'':<20} {'':<12} {'':<12} {'':<12} {'(压缩|compressed)':<12} {'(uncompressed)':<12}")

        for r in sorted(results, key=lambda x: x['sample_name']):
            sample_name = r['sample_name'][:20]
            r1_reads = format_number(r.get('R1_reads', 0))
            r2_reads = format_number(r.get('R2_reads', 0)) if r.get('R2_reads', 'N/A') != 'N/A' else 'N/A'
            total_reads = format_number(r.get('total_reads', 0))
            r1_size = r.get('R1_file_size_formatted', 'N/A')
            r1_uncomp = r.get('R1_uncompressed_size_formatted', 'N/A')

            self.logger.info(f"{sample_name:<20} {r1_reads:<12} {str(r2_reads):<12} {total_reads:<12} {r1_size:<12} {r1_uncomp:<12}")
