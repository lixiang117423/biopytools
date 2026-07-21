"""
HiTE 结果处理模块|HiTE Results Processing Module

提供HiTE和panHiTE的结果处理和汇总功能
Provides results processing and summary generation for HiTE and panHiTE
"""

import os
import glob
import json
from typing import Dict
from .utils import get_singularity_version


class HiteResultsProcessor:
    """
    HiTE 结果处理器|HiTE Results Processor

    扫描 01_hite/ 产物,生成 software_versions.yml 与 summary 到 00_pipeline_info/
    Scans 01_hite/ outputs, writes software_versions.yml + summary to 00_pipeline_info/
    """

    HITE_VERSION = "3.3.3"

    def __init__(self, config, logger):
        """
        初始化结果处理器|Initialize results processor

        Args:
            config: HiteConfig 实例|HiteConfig instance
            logger: 日志器|logger
        """
        self.config = config
        self.logger = logger
        self.output_files: Dict[str, str] = {}

    def process_results(self) -> Dict[str, str]:
        """
        扫描 HiTE 输出(在 hite_out_dir)|Scan HiTE outputs in hite_out_dir
        """
        self.logger.info("扫描HiTE输出文件|Scanning HiTE output files")

        expected_files = {
            'confident_te': 'confident_TE.cons.fa',
            'confident_ltr': 'confident_ltr_cut.fa.cons',
            'confident_tir': 'confident_tir*.fa',
            'confident_helitron': 'confident_helitron*.fa',
            'confident_non_ltr': 'confident_non_ltr*.fa',
            'all_te': 'all_TE.fa',
            'low_confident_te': 'low_confident_TE.cons.fa',
        }

        for name, pattern in expected_files.items():
            files = glob.glob(os.path.join(self.config.hite_out_dir, pattern))
            if files:
                self.output_files[name] = files[0]
                self.logger.info(f"  找到|Found {name}: {files[0]}")
            else:
                self.logger.warning(f"  未找到|Not found {name}: {pattern}")

        if self.config.annotate:
            for name, filename in {
                'hite_out': 'HiTE.out', 'hite_gff': 'HiTE.gff', 'hite_tbl': 'HiTE.tbl',
            }.items():
                fp = os.path.join(self.config.hite_out_dir, filename)
                if os.path.exists(fp):
                    self.output_files[name] = fp
                    self.logger.info(f"  找到|Found {name}: {fp}")
                else:
                    self.logger.warning(f"  未找到|Not found {name}: {filename}")

        if self.config.domain:
            domain_file = os.path.join(self.config.hite_out_dir, 'confident_TE.cons.fa.domain')
            if os.path.exists(domain_file):
                self.output_files['domain'] = domain_file
                self.logger.info(f"  找到|Found domain: {domain_file}")

        self.logger.info(
            f"共找到|Total found {len(self.output_files)} 个输出文件|output files"
        )
        return self.output_files

    def write_software_versions(self) -> None:
        """
        生成 software_versions.yml 到 00_pipeline_info/|Write software_versions.yml
        """
        import yaml
        singularity_ver = get_singularity_version(self.config.singularity_path)

        info = {
            'pipeline': {'name': 'biopytools hite', 'version': '1.0.0'},
            'tools': {
                'HiTE': {
                    'version': self.HITE_VERSION,
                    'path': self.config.sif_file,
                    'container': 'singularity',
                    'command': 'python /HiTE/main.py',
                },
                'singularity': {
                    'version': singularity_ver,
                    'path': self.config.singularity_path,
                },
            },
            'parameters': {
                'genome': self.config.genome,
                'threads': self.config.threads,
                'plant': self.config.plant,
                'annotate': self.config.annotate,
                'recover': self.config.recover,
                'domain': self.config.domain,
                'te_type': self.config.te_type,
                'chunk_size': self.config.chunk_size,
                'miu': self.config.miu,
                'min_te_len': self.config.min_te_len,
                'remove_nested': self.config.remove_nested,
            },
        }

        out_file = os.path.join(self.config.pipeline_info_dir, 'software_versions.yml')
        with open(out_file, 'w', encoding='utf-8') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"软件版本已保存|Software versions saved to: {out_file}")

    def generate_summary(self) -> Dict:
        """生成结果汇总 JSON 到 00_pipeline_info/|Generate summary JSON"""
        self.logger.info("生成结果汇总|Generating results summary")
        summary = {
            'analysis_type': 'HiTE single-genome analysis',
            'hite_out_dir': self.config.hite_out_dir,
            'output_files': self.output_files,
            'config': {
                'genome': self.config.genome,
                'threads': self.config.threads,
                'plant': self.config.plant,
                'annotate': self.config.annotate,
                'recover': self.config.recover,
                'domain': self.config.domain,
                'te_type': self.config.te_type,
            },
        }
        summary_file = os.path.join(self.config.pipeline_info_dir, 'hite_summary.json')
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        self.logger.info(f"结果汇总已保存|Summary saved to: {summary_file}")
        return summary

    def print_summary(self):
        """打印结果汇总|Print results summary"""
        self.logger.info("=" * 80)
        self.logger.info("HiTE分析结果汇总|HiTE Analysis Results Summary")
        self.logger.info("=" * 80)
        if self.output_files:
            self.logger.info(f"\n输出文件|Output files ({len(self.output_files)}):")
            for name, path in self.output_files.items():
                size = os.path.getsize(path) if os.path.exists(path) else 0
                self.logger.info(f"  {name:20s}: {path} ({size / 1048576:.2f} MB)")
        else:
            self.logger.warning("\n未找到输出文件|No output files found")
        self.logger.info("=" * 80)


# ===== PanHiteResultsProcessor 保留原样,勿动|Preserved as-is, do NOT touch =====


class PanHiteResultsProcessor:
    """
    panHiTE 结果处理器|panHiTE Results Processor

    处理panHiTE群体基因组分析的输出结果
    Process output results from panHiTE pan-genome analysis
    """

    def __init__(self, config, logger):
        """
        初始化结果处理器|Initialize results processor

        Args:
            config: PanHiteConfig配置对象|PanHiteConfig object
            logger: 日志器对象|Logger object
        """
        self.config = config
        self.logger = logger
        self.output_files: Dict[str, str] = {}

    def process_results(self) -> Dict[str, str]:
        """
        处理panHiTE输出结果|Process panHiTE output results

        扫描输出目录，查找所有预期的输出文件
        Scan output directory and find all expected output files

        Returns:
            Dict[str, str]: 输出文件字典|Dictionary of output files
        """
        self.logger.info("扫描panHiTE输出文件|Scanning panHiTE output files")

        # panHiTE特定的输出文件|panHiTE specific output files
        expected_files = {
            'pante_library': 'panTE_library.fa',
            'pante_annotation': 'panTE_annotation.gff',
        }

        # 查找输出文件|Find output files
        for name, filename in expected_files.items():
            filepath = os.path.join(self.config.output_dir, filename)
            if os.path.exists(filepath):
                self.output_files[name] = filepath
                self.logger.info(f"  找到|Found {name}: {filepath}")
            else:
                self.logger.warning(f"  未找到|Not found {name}: {filename}")

        # 查找可能的额外输出文件|Find possible additional output files
        additional_patterns = {
            'panTE_stats': '*_stats.txt',
            'panTE_cluster': '*_cluster.fa',
            'presence_absence': 'presence_absence*.txt'
        }

        for name, pattern in additional_patterns.items():
            files = glob.glob(os.path.join(self.config.output_dir, pattern))
            if files:
                self.output_files[name] = files[0]
                self.logger.info(f"  找到|Found {name}: {files[0]}")

        self.logger.info(
            f"共找到|Total found {len(self.output_files)} 个输出文件|output files"
        )

        return self.output_files

    def generate_summary(self) -> Dict:
        """
        生成结果汇总|Generate results summary

        创建JSON格式的结果汇总文件
        Create JSON format results summary file

        Returns:
            Dict: 结果汇总字典|Results summary dictionary
        """
        self.logger.info("生成结果汇总|Generating results summary")

        summary = {
            'analysis_type': 'panHiTE pan-genome analysis',
            'output_dir': self.config.output_dir,
            'output_files': self.output_files,
            'config': {
                'pan_genomes_dir': self.config.pan_genomes_dir,
                'genome_list': self.config.genome_list,
                'threads': self.config.threads,
                'skip_analyze': self.config.skip_analyze,
                'te_type': self.config.te_type
            }
        }

        # 添加可选配置|Add optional configuration
        if self.config.genes_dir:
            summary['config']['genes_dir'] = self.config.genes_dir
        if self.config.rna_dir:
            summary['config']['rna_dir'] = self.config.rna_dir

        # 保存为JSON文件|Save as JSON file
        summary_file = os.path.join(self.config.output_dir, 'panhite_summary.json')
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)

        self.logger.info(f"结果汇总已保存|Summary saved to: {summary_file}")

        return summary

    def print_summary(self):
        """
        打印结果汇总|Print results summary

        在控制台输出格式化的结果汇总
        Print formatted results summary to console
        """
        self.logger.info("="*80)
        self.logger.info("panHiTE分析结果汇总|panHiTE Analysis Results Summary")
        self.logger.info("="*80)

        if self.output_files:
            self.logger.info(f"\n输出文件|Output files ({len(self.output_files)}):")
            for name, path in self.output_files.items():
                file_size = os.path.getsize(path) if os.path.exists(path) else 0
                size_mb = file_size / (1024 * 1024)
                self.logger.info(f"  {name:20s}: {path} ({size_mb:.2f} MB)")
        else:
            self.logger.warning("\n未找到输出文件|No output files found")

        self.logger.info("="*80)
