"""
流程信息生成器模块|Pipeline Information Generator Module
"""

import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

import yaml


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
        self.start_time = None

    def get_software_versions(self) -> Dict[str, str]:
        """获取软件版本信息|Get software version information"""
        from .utils import build_conda_command

        versions = {}

        try:
            cmd = build_conda_command(self.config.orthofinder_path, ['-h'])
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            version_line = result.stdout.strip().split('\n')[0] if result.stdout.strip() else ""
            versions['orthofinder'] = version_line if "OrthoFinder" in version_line else "unknown"
        except Exception as e:
            self.logger.warning(f"无法获取orthofinder版本|Cannot get orthofinder version: {e}")
            versions['orthofinder'] = "unknown"

        import sys
        versions['python'] = f"Python {sys.version.split()[0]}"

        try:
            import pandas as pd
            versions['pandas'] = f"pandas {pd.__version__}"
        except ImportError:
            versions['pandas'] = "not installed"

        try:
            import numpy as np
            versions['numpy'] = f"numpy {np.__version__}"
        except ImportError:
            versions['numpy'] = "not installed"

        try:
            import matplotlib
            versions['matplotlib'] = f"matplotlib {matplotlib.__version__}"
        except ImportError:
            versions['matplotlib'] = "not installed"

        return versions

    def get_pipeline_params(self) -> Dict[str, Any]:
        """获取流程参数信息|Get pipeline parameters information"""
        params = {
            'analysis_name': 'OrthoFinder Pangenome Analysis',
            'project_name': self.config.project_name,
            'input_dir': self.config.input_dir,
            'output_dir': self.config.output_dir,
            'sequence_type': self.config.sequence_type,
            'threads': self.config.threads,
            'search_program': self.config.search_program,
            'mcl_inflation': self.config.mcl_inflation,
            'softcore_missing_threshold': self.config.softcore_missing_threshold,
            'dispensable_missing_threshold': self.config.dispensable_missing_threshold,
            'basic_analysis_only': self.config.basic_analysis_only,
            'generate_trees': self.config.generate_trees,
            'enable_single_copy_analysis': self.config.enable_single_copy_analysis,
            'extract_sequences': self.config.extract_sequences,
            'single_copy_output_format': self.config.single_copy_output_format,
            'enable_rarefaction': self.config.enable_rarefaction,
            'rarefaction_iterations': self.config.rarefaction_iterations,
            'generate_plots': self.config.generate_plots,
            'plot_format': self.config.plot_format,
            'orthofinder_path': self.config.orthofinder_path,
            'msa_program': self.config.msa_program,
            'tree_program': self.config.tree_program,
            'resume_from_existing': self.config.resume_from_existing,
            'skip_orthofinder': self.config.skip_orthofinder,
            'force_overwrite': self.config.force_overwrite
        }
        return params

    def generate_pipeline_info(self, output_dir: str):
        """生成流程信息文件|Generate pipeline information files"""
        self.start_time = datetime.now()

        pipeline_info_dir = Path(output_dir) / "00_pipeline_info"
        pipeline_info_dir.mkdir(parents=True, exist_ok=True)

        software_versions = self.get_software_versions()
        pipeline_params = self.get_pipeline_params()
        timestamp = self.start_time.strftime("%Y-%m-%d %H:%M:%S")

        versions_file = pipeline_info_dir / "software_versions.yml"
        versions_data = {
            'pipeline': {
                'name': 'biopytools orthofinder',
                'version': '1.0.0'
            },
            'tools': {
                name: {
                    'version': ver
                } for name, ver in software_versions.items()
            },
            'parameters': pipeline_params,
            'execution': {
                'start_time': timestamp,
                'end_time': '',
                'runtime_seconds': 0
            }
        }
        with open(versions_file, 'w', encoding='utf-8') as f:
            yaml.dump(versions_data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

        params_file = pipeline_info_dir / "pipeline_params.yml"
        params_data = {
            'generated_at': timestamp,
            'pipeline': 'OrthoFinder Pangenome Analysis v1.0',
            'parameters': pipeline_params
        }
        with open(params_file, 'w', encoding='utf-8') as f:
            yaml.dump(params_data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

        self.logger.info(f"流程信息已生成|Pipeline info generated: {pipeline_info_dir}")

    def finalize_pipeline_info(self, output_dir: str):
        """完成时更新流程信息|Update pipeline info on completion"""
        if not self.start_time:
            return

        end_time = datetime.now()
        runtime_seconds = int((end_time - self.start_time).total_seconds())

        pipeline_info_dir = Path(output_dir) / "00_pipeline_info"
        versions_file = pipeline_info_dir / "software_versions.yml"

        if versions_file.exists():
            try:
                with open(versions_file, 'r', encoding='utf-8') as f:
                    data = yaml.safe_load(f)

                data['execution'] = {
                    'start_time': self.start_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'runtime_seconds': runtime_seconds
                }

                with open(versions_file, 'w', encoding='utf-8') as f:
                    yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
            except Exception:
                pass
