"""MGA共识基因组组装主模块|MGA consensus genome assembly main module"""

import argparse
import gzip
import subprocess
import sys
from datetime import datetime

from ..common.paths import expand_path
from .config import MGAConfig
from .utils import (
    MGALogger,
    build_mga_command,
    check_dependencies,
    write_software_versions,
)


class MGARunner:
    """MGA组装主类|MGA assembly main class"""

    def __init__(self, **kwargs):
        # 初始化配置(含validate)|Init config (with validate)
        self.config = MGAConfig(**kwargs)
        self.config.validate()
        # 初始化日志(写入99_logs)|Init logging (to 99_logs)
        self.logger_manager = MGALogger(self.config.logs_path)
        self.logger = self.logger_manager.get_logger()

    def _check_read_names(self):
        """扫描reads首行,read name含空格则WARNING(不自动改)|Warn on whitespace in read names"""
        reads = self.config.reads
        bad = False
        try:
            opener = gzip.open if reads.endswith(".gz") else open
            with opener(reads, "rt") as fh:
                checked = 0
                for line in fh:
                    if line.startswith("@") or line.startswith(">"):
                        if " " in line.rstrip("\n"):
                            bad = True
                        checked += 1
                        if checked >= 10:
                            break
        except Exception as e:
            self.logger.warning(f"无法扫描read name|Cannot scan read names: {e}")
            return
        if bad:
            self.logger.warning(
                "检测到read name含空格,MGA/LJA可能出错;建议先用trim_header去空格|"
                "Read names contain whitespace; MGA/LJA may fail. Trim with trim_header first."
            )

    def run(self) -> bool:
        """运行MGA组装(断点续传+dry-run)|Run MGA assembly (checkpoint + dry-run)"""
        cfg = self.config

        self.logger.info("=" * 60)
        self.logger.info("MGA共识基因组组装|MGA consensus genome assembly")
        self.logger.info("=" * 60)
        self.logger.info(f"输入reads|Input reads: {cfg.reads}")
        self.logger.info(f"输出目录|Output dir: {cfg.output_dir}")
        self.logger.info(f"线程数|Threads: {cfg.threads}")
        self.logger.info(f"conda环境|conda env: {cfg.conda_env}")

        # 断点续传:最终产物存在则整体跳过|Checkpoint: skip if final product exists
        if cfg.assembly_fasta.exists():
            self.logger.info(f"最终产物已存在,跳过|Final product exists, skip: {cfg.assembly_fasta}")
            return True

        # read name检查|read name check
        self._check_read_names()

        # 依赖检查(dry_run跳过,便于拿命令)|dependency check (skip in dry-run)
        versions = None
        if not cfg.dry_run:
            versions = check_dependencies(cfg, self.logger)
            if versions is None:
                self.logger.error("依赖检查失败|Dependency check failed")
                return False

        # 构建命令(显式conda run)|build command (explicit conda run)
        cmd = build_mga_command(cfg.mga_path, cfg.conda_env, cfg.reads,
                                cfg.output_dir, cfg.threads)
        self.logger.info(f"执行|Executing: MGA共识组装|MGA consensus assembly")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        if cfg.dry_run:
            self.logger.info("dry-run模式,不执行|dry-run mode, not executing")
            return True

        # 执行:stdout/stderr透传(实时显示,超算.out/.err捕获)|execute (passthrough)
        start_time = datetime.now()
        try:
            subprocess.run(cmd, shell=False, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MGA执行失败|MGA failed (exit {e.returncode})")
            write_software_versions(cfg, self.logger, start_time, datetime.now(), versions or {})
            return False
        end_time = datetime.now()

        # 版本记录|version info
        write_software_versions(cfg, self.logger, start_time, end_time, versions or {})

        # 结果检查|result check
        if cfg.assembly_fasta.exists():
            self.logger.info(f"组装完成|Assembly completed: {cfg.assembly_fasta}")
            self.logger.info("=" * 60)
            return True
        self.logger.error(f"未找到最终产物|Final product not found: {cfg.assembly_fasta}")
        return False


def main():
    """命令行入口|CLI entry"""
    parser = argparse.ArgumentParser(
        description="MGA共识基因组组装|MGA consensus genome assembly",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools mga -r reads.fastq.gz -o out_dir/",
    )
    parser.add_argument("-r", "--reads", required=True,
                        help="[FILE] HiFi reads(fasta/fastq,可gz)|HiFi reads (fasta/fastq, may be gz)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="[DIR] 输出目录|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=50,
                        help="线程数(默认50)|Threads (default 50)")
    parser.add_argument("--mga-path", default=None,
                        help="MGA二进制路径|MGA binary path")
    parser.add_argument("--conda-env", default="mga",
                        help="conda环境名(默认mga)|conda env name (default mga)")
    parser.add_argument("--dry-run", action="store_true",
                        help="只打印命令不执行|Print command without executing")
    args = parser.parse_args()

    kwargs = dict(
        reads=args.reads,
        output_dir=args.output_dir,
        threads=args.threads,
        conda_env=args.conda_env,
        dry_run=args.dry_run,
    )
    if args.mga_path:
        kwargs["mga_path"] = expand_path(args.mga_path)

    try:
        runner = MGARunner(**kwargs)
        success = runner.run()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        sys.exit(130)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"运行失败|Run failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
