"""trimal 主程序模块|trimal Main Module"""

import os
import sys
from typing import List

from .config import TrimalConfig
from .utils import TrimalLogger, run_command


class TrimalRunner:
    """trimal 运行器|trimal Runner"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        self.config = TrimalConfig(**kwargs)
        self.config.validate()

        self.logger_manager = TrimalLogger(
            log_file=self.config.log_file,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()
        self.logger.info("trimal 运行器已初始化|trimal Runner initialized")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"方法|Method: {self.config.method}")

    # ------------------------------------------------------------------ #
    # 命令构建(纯函数,便于测试)|Command building (pure, testable)
    # ------------------------------------------------------------------ #
    def _method_args(self) -> List[str]:
        """方法→trimal 标志|method→trimal flag"""
        m = self.config.method
        if m == 'automated1':
            return ['-automated1']
        if m == 'gappyout':
            return ['-gappyout']
        if m == 'strict':
            return ['-strict']
        if m == 'strictplus':
            return ['-strictplus']
        if m == 'gt':
            # 实测 v1.5: -gt <value>(非 -gapthreshold)|verified v1.5: -gt <value>
            return ['-gt', str(self.config.gt_threshold)]
        if m == 'cons':
            return ['-cons', str(self.config.cons_threshold)]
        return []

    def _format_arg(self) -> List[str]:
        """格式→trimal 标志|format→trimal flag"""
        f = self.config.output_format
        if f == 'keep':
            return []
        return [f'-{f}']

    def _output_ext(self) -> str:
        """输出扩展名|output extension"""
        if self.config.output_format == 'keep':
            ext = os.path.splitext(self.config.input_file)[1]
            return ext if ext else '.fasta'
        return f'.{self.config.output_format}'

    def _build_args(self, kind: str, out_file: str) -> List[str]:
        """
        构建 trimal 参数列表|Build trimal arg list

        kind: 'main' | 'complementary' | 'backtrans'
        """
        args = ['-in', self.config.input_file]

        # backtrans 需 CDS 文件,且其 -out 接收反向翻译的 NT 比对
        # |backtrans needs CDS file; its -out receives the backtranslated NT alignment
        if kind == 'backtrans' and self.config.backtrans_file:
            args += ['-backtrans', self.config.backtrans_file]

        args += self._method_args()

        # complementary 的 -out 接收互补比对|complementary -out receives the complementary alignment
        if kind == 'complementary':
            args += ['-complementary']

        args += self._format_arg()

        if self.config.keep_header:
            args += ['-keepheader']

        # colnumbering 走 stdout,可与主调用合并|colnumbering goes to stdout, combinable with main run
        if kind == 'main' and self.config.colnumbering:
            args += ['-colnumbering']

        args += ['-out', out_file]
        return args

    def _run_trimal(self, args: List[str], description: str = ""):
        """直接以绝对路径执行 trimal(独立二进制,无需 conda 包装)|Run trimal directly by abs path"""
        cmd = [self.config.trimal_path] + args
        return run_command(cmd, self.logger, description)

    # ------------------------------------------------------------------ #
    # 主流程|Main pipeline
    # ------------------------------------------------------------------ #
    def run_extraction(self) -> bool:
        """运行修剪流程|Run trimming pipeline"""
        try:
            # 创建输出目录|Create output dirs
            trimal_dir = os.path.join(self.config.output_dir, '01_trimal')
            info_dir = os.path.join(self.config.output_dir, '00_pipeline_info')
            log_dir = os.path.join(self.config.output_dir, '99_logs')
            for d in (trimal_dir, info_dir, log_dir):
                os.makedirs(d, exist_ok=True)

            self._setup_default_log(log_dir)

            ext = self._output_ext()
            sample = self.config.sample_name
            trimmed = os.path.join(trimal_dir, f'{sample}.trimmed{ext}')

            # 主修剪(断点续传)|Main trim (checkpoint resume)
            if os.path.exists(trimmed):
                self.logger.info(f"跳过已完成步骤|Skipping completed step: 主修剪|main trim")
            else:
                self.logger.info("开始主修剪|Starting main trimming")
                stdout, _ = self._run_trimal(
                    self._build_args('main', trimmed),
                    '主修剪|main trim'
                )
                self.logger.info(f"修剪比对已保存|Trimmed alignment saved: {trimmed}")

                # colnumbering 从主调用 stdout 写入|colnumbering written from main-run stdout
                if self.config.colnumbering:
                    cn_file = os.path.join(trimal_dir, f'{sample}.colnumbering.tsv')
                    with open(cn_file, 'w', encoding='utf-8') as f:
                        f.write(stdout)
                    self.logger.info(f"列号映射已保存|Column mapping saved: {cn_file}")

            # 互补比对(单独调用,-out 被其占用)|Complementary (separate run)
            if self.config.complementary:
                comp = os.path.join(trimal_dir, f'{sample}.complementary{ext}')
                self.logger.info("生成互补比对|Generating complementary alignment")
                self._run_trimal(
                    self._build_args('complementary', comp),
                    '互补比对|complementary'
                )
                self.logger.info(f"互补比对已保存|Complementary alignment saved: {comp}")

            # 反向翻译(单独调用,-out 被其占用)|Backtranslation (separate run)
            if self.config.backtrans_file:
                bt = os.path.join(trimal_dir, f'{sample}.backtrans{ext}')
                self.logger.info("生成反向翻译比对|Generating backtranslated alignment")
                self._run_trimal(
                    self._build_args('backtrans', bt),
                    '反向翻译|backtrans'
                )
                self.logger.info(f"反向翻译比对已保存|Backtranslated alignment saved: {bt}")

            # 版本信息|Software versions
            self._write_software_versions(info_dir)

            self.logger.info("trimal 流程完成|trimal pipeline completed")
            return True

        except Exception as e:
            self.logger.error(f"流程失败|Pipeline failed: {e}")
            return False

    def _setup_default_log(self, log_dir):
        """若无 log_file,默认写入 99_logs/trimal.log|Default log file if none specified"""
        if self.config.log_file is None:
            self.config.log_file = os.path.join(log_dir, 'trimal.log')
            import logging
            formatter = logging.Formatter(
                '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            fh = logging.FileHandler(self.config.log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

    def _get_trimal_version(self) -> str:
        """获取 trimal 版本|Get trimal version"""
        try:
            cmd = [self.config.trimal_path, '--version']
            out, err = run_command(cmd, self.logger, "trimal 版本|trimal version")
            return (out.strip() or err.strip()) or 'unknown'
        except Exception as e:
            self.logger.warning(f"无法获取 trimal 版本|Cannot get trimal version: {e}")
            return 'unknown'

    def _write_software_versions(self, info_dir):
        """写 software_versions.yml|Write software_versions.yml"""
        try:
            import yaml
        except ImportError:
            self.logger.warning("未安装 pyyaml,跳过 software_versions.yml|pyyaml missing, skip")
            return

        info = {
            'pipeline': {'name': 'biopytools trimal', 'version': '1.0.0'},
            'tools': {
                'trimal': {
                    'version': self._get_trimal_version(),
                    'path': self.config.trimal_path,
                    'command': 'trimal --version',
                }
            },
            'parameters': {
                'method': self.config.method,
                'gt_threshold': self.config.gt_threshold,
                'cons_threshold': self.config.cons_threshold,
                'output_format': self.config.output_format,
                'colnumbering': self.config.colnumbering,
                'complementary': self.config.complementary,
                'backtrans': bool(self.config.backtrans_file),
                'keep_header': self.config.keep_header,
                'sample_name': self.config.sample_name,
                'input_file': self.config.input_file,
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
        description="trimal 多序列比对自动修剪|MSA automated trimming with trimAl",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入比对文件|Input alignment file')
    parser.add_argument('-o', '--output-dir', default='./trimal_output',
                        help='输出目录|Output directory (default: ./trimal_output)')
    parser.add_argument('--method', default='automated1',
                        choices=TrimalConfig.VALID_METHODS,
                        help='修剪方法|Trimming method (default: automated1)')
    parser.add_argument('--gt-threshold', type=float, default=0.9,
                        help='gap 阈值[0,1],仅 method=gt|gap threshold [0,1], only method=gt (default: 0.9)')
    parser.add_argument('--cons-threshold', type=int, default=80,
                        help='保守度阈值[0,100],仅 method=cons|conservation [0,100], only method=cons (default: 80)')
    parser.add_argument('--format', default='keep',
                        choices=TrimalConfig.VALID_FORMATS,
                        help='输出格式|Output format (default: keep=沿用输入|input format)')
    parser.add_argument('--colnumbering', action='store_true',
                        help='输出新旧列号映射|Output old/new column mapping')
    parser.add_argument('--backtrans',
                        help='CDS 文件,AA→NT 反向翻译|CDS file for AA→NT backtranslation')
    parser.add_argument('--complementary', action='store_true',
                        help='输出被修剪列的互补比对|Output complementary alignment of trimmed columns')
    parser.add_argument('--keep-header', action='store_true',
                        help='保留完整 FASTA 头|Keep full FASTA headers')
    parser.add_argument('--sample-name',
                        help='输出文件名前缀|Output filename prefix (default: 输入 basename|input basename)')
    parser.add_argument('--log-file',
                        help='日志文件路径|Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出|Verbose output')

    args = parser.parse_args()

    try:
        runner = TrimalRunner(
            input_file=args.input,
            output_dir=args.output_dir,
            method=args.method,
            gt_threshold=args.gt_threshold,
            cons_threshold=args.cons_threshold,
            output_format=args.format,
            colnumbering=args.colnumbering,
            backtrans_file=args.backtrans,
            complementary=args.complementary,
            keep_header=args.keep_header,
            sample_name=args.sample_name,
            log_file=args.log_file,
            verbose=args.verbose
        )
        success = runner.run_extraction()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"程序执行失败|Program execution failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
