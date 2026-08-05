"""PSVCP 主程序模块|PSVCP Main Module

PSVCPRunner 组合 config + logger + PangenomeConstructor;main() 为 argparse 入口。
|PSVCPRunner composes config + logger + PangenomeConstructor; main() is the argparse entry.
"""

import os
import sys

from .config import PSVCPConfig
from .pipeline import PangenomeConstructor
from .utils import PSVCPLogger, build_conda_command


class PSVCPRunner:
    """PSVCP 运行器|PSVCP Runner"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # log_file / verbose 不属于 config 字段,单独取出|separate non-config kwargs
        log_file = kwargs.pop('log_file', None)
        verbose = kwargs.pop('verbose', False)

        self.config = PSVCPConfig(**kwargs)
        self.config.validate()

        # 日志默认写 output_dir/psvcp.log|default log file
        os.makedirs(self.config.output_dir, exist_ok=True)
        log_file = log_file or os.path.join(self.config.output_dir, 'psvcp.log')
        self.logger_manager = PSVCPLogger(log_file=log_file, verbose=verbose)
        self.logger = self.logger_manager.get_logger()

        self.constructor = PangenomeConstructor(self.config, self.logger)
        self.logger.info("PSVCP 泛基因组构建运行器已初始化|PSVCP pangenome Runner initialized")
        self.logger.info(f"genome_dir   : {self.config.genome_dir}")
        self.logger.info(f"genome_list  : {self.config.genome_list}")
        self.logger.info(f"output_dir   : {self.config.output_dir}")
        self.logger.info(f"threads      : {self.config.threads}")

    def run(self) -> bool:
        """运行构建|Run construction"""
        try:
            success = self.constructor.run()
            if success:
                self._write_software_versions()
            return success
        except Exception as e:
            self.logger.error(f"运行失败|Run failed: {e}")
            return False

    def _write_software_versions(self):
        """写 software_versions.yml(§12.5)|write software_versions.yml"""
        try:
            import yaml
        except ImportError:
            self.logger.warning("未安装 pyyaml,跳过 software_versions.yml|pyyaml missing, skip")
            return

        info_dir = os.path.join(self.config.output_dir, '00_pipeline_info')
        os.makedirs(info_dir, exist_ok=True)
        tools = {}
        for name, path, ver_args in [
            ('nucmer', self.config.nucmer_path, ['--version']),
            ('assemblytics', self.config.assemblytics_path, ['--version']),
            ('samtools', self.config.samtools_path, ['--version']),
            ('bedtools', self.config.bedtools_path, ['--version']),
        ]:
            try:
                cmd = build_conda_command(path, ver_args)
                import subprocess
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                tools[name] = {'version': (r.stdout or r.stderr or '').strip() or 'unknown',
                               'path': path}
            except Exception:
                tools[name] = {'version': 'unknown', 'path': path}

        info = {
            'pipeline': {'name': 'biopytools psvcp', 'version': '1.0.0'},
            'tools': tools,
            'parameters': {
                'threads': self.config.threads,
                'assemblytics_unique_length': self.config.ASSEMBLYTICS_UNIQUE_LENGTH,
                'assemblytics_min_size': self.config.ASSEMBLYTICS_MIN_SIZE,
                'assemblytics_max_size': self.config.ASSEMBLYTICS_MAX_SIZE,
                'nucmer_maxgap': self.config.NUCMER_MAXGAP,
                'nucmer_mincluster': self.config.NUCMER_MINCLUSTER,
                'nucmer_diagdiff': self.config.NUCMER_DIAGDIFF,
                'min_insertion_size': self.config.MIN_INSERTION_SIZE,
                'genome_list': self.config.read_genome_list(),
            },
        }
        out_file = os.path.join(info_dir, 'software_versions.yml')
        with open(out_file, 'w', encoding='utf-8') as f:
            yaml.safe_dump(info, f, default_flow_style=False, sort_keys=False)
        self.logger.info(f"版本信息已保存|Version info saved: {out_file}")


def main():
    """命令行入口|CLI entry"""
    import argparse

    parser = argparse.ArgumentParser(
        description="PSVCP 线性泛基因组构建(MUMmer+Assemblytics,PAV 并入)|"
                    "PSVCP linear pangenome construction (MUMmer+Assemblytics, PAV incorporation)",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--genome-dir', required=True,
                        help='genome 目录(含 {name}.fa + {name}.gff/.gff3)|genome dir with {name}.fa + {name}.gff/.gff3')
    parser.add_argument('-l', '--genome-list', required=True,
                        help='genome_list 文本(行1=ref,其余=query)|genome_list (line1=ref, rest=queries)')
    parser.add_argument('-o', '--output-dir', default='~/psvcp_out',
                        help='输出目录|output directory (default: ~/psvcp_out)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|threads (default: 12)')
    parser.add_argument('--force', action='store_true',
                        help='忽略断点续传,强制重跑|ignore checkpoint, force rerun')
    parser.add_argument('--log-file', help='日志文件路径(默认 output_dir/psvcp.log)|log file path')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细输出|verbose output')

    args = parser.parse_args()

    try:
        runner = PSVCPRunner(
            genome_dir=args.genome_dir,
            genome_list=args.genome_list,
            output_dir=args.output_dir,
            threads=args.threads,
            force=args.force,
            log_file=args.log_file,
            verbose=args.verbose,
        )
        success = runner.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"程序执行失败|Program execution failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
