"""
流程信息生成器模块|Pipeline Information Generator Module
"""

import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, Any


class PipelineInfoGenerator:
    """流程信息生成器|Pipeline Information Generator"""

    def __init__(self, logger, config):
        """
        初始化生成器|Initialize generator

        Args:
            logger: 日志器|Logger
            config: 配置对象|Configuration object
        """
        self.logger = logger
        self.config = config

    def get_software_versions(self) -> Dict[str, str]:
        """
        获取软件版本信息|Get software version information

        Returns:
            Dict[str, str]: 软件版本字典|Software version dictionary
        """
        versions = {}

        # 获取bowtie2版本|Get bowtie2 version
        try:
            result = subprocess.run(
                [self.config.bowtie2_path, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            bowtie2_version = result.stdout.strip().split('\n')[0]
            versions['bowtie2'] = bowtie2_version
        except Exception as e:
            self.logger.warning(f"无法获取bowtie2版本|Cannot get bowtie2 version: {e}")
            versions['bowtie2'] = "unknown"

        # 获取samtools版本|Get samtools version
        try:
            result = subprocess.run(
                [self.config.samtools_path, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            samtools_version = result.stdout.strip().split('\n')[0]
            versions['samtools'] = samtools_version
        except Exception as e:
            self.logger.warning(f"无法获取samtools版本|Cannot get samtools version: {e}")
            versions['samtools'] = "unknown"

        # 获取Python版本|Get Python version
        import sys
        versions['python'] = f"Python {sys.version.split()[0]}"

        return versions

    def get_pipeline_params(self) -> Dict[str, Any]:
        """
        获取流程参数信息|Get pipeline parameters information

        Returns:
            Dict[str, Any]: 流程参数字典|Pipeline parameters dictionary
        """
        params = {
            'analysis_name': 'Insertion Site Detection',
            'output_dir': self.config.output_dir,
            'genome': self.config.genome,
            'insert_sequence': self.config.insert_sequence,
            'threads': self.config.threads,
            'min_clip': self.config.min_clip,
            'min_support': self.config.min_support,
            'score_threshold': self.config.score_threshold,
            'skip_existing': self.config.skip_existing,
            'bowtie2_path': self.config.bowtie2_path,
            'samtools_path': self.config.samtools_path,
            'read1_suffix': self.config.read1_suffix,
            'read2_suffix': self.config.read2_suffix
        }
        return params

    def generate_pipeline_info(self, output_dir: str, sample_id: str = None):
        """
        生成流程信息文件|Generate pipeline information files

        Args:
            output_dir: 输出目录|Output directory
            sample_id: 样品ID（可选）|Sample ID (optional)
        """
        # 创建00_pipeline_info目录|Create 00_pipeline_info directory
        if sample_id:
            pipeline_info_dir = Path(output_dir) / sample_id / "00_pipeline_info"
        else:
            pipeline_info_dir = Path(output_dir) / "00_pipeline_info"

        pipeline_info_dir.mkdir(parents=True, exist_ok=True)

        # 获取软件版本|Get software versions
        software_versions = self.get_software_versions()

        # 获取流程参数|Get pipeline parameters
        pipeline_params = self.get_pipeline_params()

        # 添加时间戳|Add timestamp
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # 写入software_versions.yml|Write software_versions.yml
        versions_file = pipeline_info_dir / "software_versions.yml"
        with open(versions_file, 'w') as f:
            f.write(f"# Insert Detection Software Versions\n")
            f.write(f"# Generated at: {timestamp}\n\n")
            f.write("software_versions:\n")
            for software, version in software_versions.items():
                f.write(f"  {software}: \"{version}\"\n")

        # 写入pipeline_params.yml|Write pipeline_params.yml
        params_file = pipeline_info_dir / "pipeline_params.yml"
        with open(params_file, 'w') as f:
            f.write(f"# Insert Detection Pipeline Parameters\n")
            f.write(f"# Generated at: {timestamp}\n\n")
            f.write("pipeline_parameters:\n")
            for key, value in pipeline_params.items():
                if isinstance(value, str):
                    f.write(f"  {key}: \"{value}\"\n")
                else:
                    f.write(f"  {key}: {value}\n")

        self.logger.info(f"流程信息已生成|Pipeline info generated: {pipeline_info_dir}")
