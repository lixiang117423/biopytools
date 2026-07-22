"""phylo_trim 主程序模块|phylo_trim Main Module

整合 mafft-fasttree + trimal:输入序列 → MAFFT → FastTree(trimal前) → trimal → FastTree(trimal后)
|Integrates mafft-fasttree + trimal: sequences → MAFFT → FastTree(before) → trimal → FastTree(after)
仅 import 复用,不修改 mafft_fasttree / trimal 两模块|Import-only reuse; does not modify either module
"""

import os
import sys
import subprocess
from pathlib import Path

from .config import PhyloTrimConfig
from .utils import PhyloTrimLogger
# 复用现有模块(只 import)|reuse existing modules (import only)
from ..mafft_fasttree import PhyloTreeBuilder, PhyloConfig
from ..mafft_fasttree.sequence_processor import SequenceProcessor
from ..mafft_fasttree.tree_builder import FastTreeBuilder
from ..mafft_fasttree.utils import CommandRunner
from ..trimal import TrimalRunner


class PhyloTrimRunner:
    """phylo_trim 运行器|phylo_trim Runner"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        self.config = PhyloTrimConfig(**kwargs)
        self.config.validate()

        self.logger_manager = PhyloTrimLogger(
            log_file=self.config.log_file, verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()
        self.logger.info("phylo_trim 运行器已初始化|phylo_trim Runner initialized")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"trimal|trimal: {'跳过|skip' if self.config.skip_trimal else self.config.trimal_method}")

    # ------------------------------------------------------------------ #
    # 辅助|Helpers
    # ------------------------------------------------------------------ #
    def _detect_seq_type(self) -> str:
        """检测序列类型(用户指定优先)|Detect seq type (user-specified takes priority)"""
        if self.config.seq_type:
            return self.config.seq_type
        proc = SequenceProcessor(self.config, self.logger)
        return proc.detect_sequence_type(self.config.input_file)

    def _trimmed_path(self, trimal_dir: str, mafft_fa: str) -> str:
        """推算 trimal 修剪输出路径(与 TrimalRunner 命名一致)|Compute trimmed output path"""
        if self.config.trimal_format == 'keep':
            ext = os.path.splitext(mafft_fa)[1] or '.fasta'
        else:
            ext = f'.{self.config.trimal_format}'
        return os.path.join(trimal_dir, '01_trimal', f'{self.config.sample_name}.trimmed{ext}')

    # ------------------------------------------------------------------ #
    # 子流程|Sub-flows
    # ------------------------------------------------------------------ #
    def _run_raw_flow(self, seq_type: str, before_dir: str):
        """前流程:mafft → fasttree(产出 trimal 前树)|raw flow: mafft → fasttree (before-tree)"""
        self.logger.info("运行前流程(mafft + fasttree)|Running raw flow (mafft + fasttree)")
        builder = PhyloTreeBuilder(
            input_file=self.config.input_file,
            output_dir=before_dir,
            seq_type=seq_type,
            threads=self.config.threads,
            mafft_params=self.config.mafft_params,
            fasttree_params=self.config.fasttree_params,
            mafft_path=self.config.mafft_path,
            fasttree_path=self.config.fasttree_path,
        )
        try:
            builder.run_pipeline()
        except SystemExit as e:
            # PhyloTreeBuilder.run_pipeline 失败时内部 sys.exit(1)|exits on failure
            if e.code:
                raise RuntimeError(f"前流程失败|Raw flow failed (exit={e.code})")

    def _run_trimal(self, mafft_fa: str, trimal_dir: str):
        """trimal 修剪|trimal trimming"""
        self.logger.info("运行 trimal 修剪|Running trimal trimming")
        runner = TrimalRunner(
            input_file=mafft_fa,
            output_dir=trimal_dir,
            method=self.config.trimal_method,
            gt_threshold=self.config.gt_threshold,
            cons_threshold=self.config.cons_threshold,
            output_format=self.config.trimal_format,
            sample_name=self.config.sample_name,
            trimal_path=self.config.trimal_path,
        )
        if not runner.run_extraction():
            raise RuntimeError("trimal 修剪失败|trimal failed")

    def _run_after_tree(self, trimmed: str, seq_type: str, after_dir: str, after_nwk: str):
        """后树:fasttree(trimal 后比对)|after-tree: fasttree on trimmed alignment"""
        self.logger.info("构建后树(fasttree on trimmed)|Building after-tree (fasttree on trimmed)")
        aft_config = PhyloConfig(
            input_file=trimmed,
            output_dir=after_dir,
            seq_type=seq_type,
            fasttree_path=self.config.fasttree_path,
            fasttree_params=self.config.fasttree_params,
        )
        cmd_runner = CommandRunner(self.logger, Path(after_dir))
        builder = FastTreeBuilder(aft_config, self.logger, cmd_runner)
        if not builder.build_tree(Path(trimmed), Path(after_nwk), seq_type):
            raise RuntimeError("后树构建失败|After-tree build failed")

    # ------------------------------------------------------------------ #
    # 主流程|Main pipeline
    # ------------------------------------------------------------------ #
    def run(self) -> bool:
        """运行整合流程|Run the integrated pipeline"""
        try:
            out = self.config.output_dir
            before_dir = os.path.join(out, '01_mafft_fasttree')
            trimal_dir = os.path.join(out, '02_trimal')
            after_dir = os.path.join(out, '03_fasttree_trimmed')
            info_dir = os.path.join(out, '00_pipeline_info')
            log_dir = os.path.join(out, '99_logs')
            for d in (before_dir, trimal_dir, after_dir, info_dir, log_dir):
                os.makedirs(d, exist_ok=True)

            self._setup_default_log(log_dir)

            base = self.config.sample_name
            mafft_fa = os.path.join(before_dir, f'{base}.mafft.fa')
            before_nwk = os.path.join(before_dir, f'{base}.fasttree.nwk')
            after_nwk = os.path.join(after_dir, f'{base}.fasttree.trimmed.nwk')

            seq_type = self._detect_seq_type()
            self.logger.info(f"序列类型|Sequence type: {seq_type}")

            # 前流程(断点续传)|raw flow (checkpoint)
            if os.path.exists(before_nwk):
                self.logger.info(f"跳过已完成步骤|Skipping completed step: 前流程|raw flow")
            else:
                self._run_raw_flow(seq_type, before_dir)

            # trimal + 后树|trimal + after-tree
            if not self.config.skip_trimal:
                trimmed = self._trimmed_path(trimal_dir, mafft_fa)
                if os.path.exists(trimmed):
                    self.logger.info(f"跳过已完成步骤|Skipping completed step: trimal")
                else:
                    self._run_trimal(mafft_fa, trimal_dir)

                if os.path.exists(after_nwk):
                    self.logger.info(f"跳过已完成步骤|Skipping completed step: 后树|after-tree")
                else:
                    self._run_after_tree(trimmed, seq_type, after_dir, after_nwk)

            self._write_software_versions(info_dir)
            self.logger.info("phylo_trim 流程完成|phylo_trim pipeline completed")
            return True

        except Exception as e:
            self.logger.error(f"流程失败|Pipeline failed: {e}")
            return False

    def _setup_default_log(self, log_dir):
        """若无 log_file,默认写 99_logs/phylo_trim.log|Default log file if none specified"""
        if self.config.log_file is None:
            import logging
            self.config.log_file = os.path.join(log_dir, 'phylo_trim.log')
            formatter = logging.Formatter(
                '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            fh = logging.FileHandler(self.config.log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

    def _get_tool_version(self, tool: str, args) -> str:
        """best-effort 获取版本|best-effort version capture

        取首行;若输出像帮助文本(FastTree -version 触发 usage)则返回 unknown
        |Take first line; return unknown if output looks like help text
        """
        try:
            r = subprocess.run([tool] + list(args), shell=False, capture_output=True,
                               text=True, timeout=15)
            out = (r.stdout or '').strip() or (r.stderr or '').strip()
            if not out:
                return 'unknown'
            if 'Unknown' in out or 'Common options' in out or len(out.splitlines()) > 5:
                return 'unknown'
            return out.splitlines()[0].strip() or 'unknown'
        except Exception:
            return 'unknown'

    def _write_software_versions(self, info_dir):
        """写 software_versions.yml|Write software_versions.yml"""
        try:
            import yaml
        except ImportError:
            self.logger.warning("未安装 pyyaml,跳过 software_versions.yml|pyyaml missing, skip")
            return

        info = {
            'pipeline': {'name': 'biopytools phylo-trim', 'version': '1.0.0'},
            'tools': {
                'mafft': {
                    'version': self._get_tool_version(self.config.mafft_path, ['--version']),
                    'path': self.config.mafft_path,
                },
                'fasttree': {
                    'version': self._get_tool_version(self.config.fasttree_path, ['-version']),
                    'path': self.config.fasttree_path,
                },
                'trimal': {
                    'version': self._get_tool_version(self.config.trimal_path, ['--version']),
                    'path': self.config.trimal_path,
                },
            },
            'parameters': {
                'skip_trimal': self.config.skip_trimal,
                'trimal_method': self.config.trimal_method,
                'gt_threshold': self.config.gt_threshold,
                'cons_threshold': self.config.cons_threshold,
                'trimal_format': self.config.trimal_format,
                'seq_type': seq_type_value(self.config),
                'threads': self.config.threads,
                'mafft_params': self.config.mafft_params,
                'fasttree_params': self.config.fasttree_params,
                'input_file': self.config.input_file,
            },
        }
        out_file = os.path.join(info_dir, 'software_versions.yml')
        with open(out_file, 'w', encoding='utf-8') as f:
            yaml.safe_dump(info, f, default_flow_style=False, sort_keys=False)
        self.logger.info(f"版本信息已保存|Version info saved: {out_file}")


def seq_type_value(config):
    """避免重复调用检测,仅记录用户指定或'auto'|record user value or 'auto'"""
    return config.seq_type or 'auto'


def main():
    """命令行入口|CLI entry"""
    import argparse

    parser = argparse.ArgumentParser(
        description="整合 mafft-fasttree + trimal,输出 trimal 前后两棵树"
                    "|Integrate mafft-fasttree + trimal, output before/after-trimal trees",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True,
                        help='输入原始序列 FASTA|Input raw sequence FASTA')
    parser.add_argument('-o', '--output-dir', default='./phylo_trim_output',
                        help='输出目录|Output directory (default: ./phylo_trim_output)')
    parser.add_argument('--skip-trimal', action='store_true',
                        help='跳过 trimal,只出 trimal 前树|Skip trimal, before-tree only')

    # trimal 参数|trimal params
    parser.add_argument('--trimal-method', default='automated1',
                        choices=PhyloTrimConfig.VALID_TRIMAL_METHODS,
                        help='trimal 方法|trimal method (default: automated1)')
    parser.add_argument('--gt-threshold', type=float, default=0.9,
                        help='trimal gap 阈值[0,1]|trimal gap threshold [0,1] (default: 0.9)')
    parser.add_argument('--cons-threshold', type=int, default=80,
                        help='trimal 保守度[0,100]|trimal conservation [0,100] (default: 80)')
    parser.add_argument('--trimal-format', default='keep',
                        choices=PhyloTrimConfig.VALID_TRIMAL_FORMATS,
                        help='trimal 输出格式|trimal output format (default: keep)')

    # mafft/fasttree 参数|mafft/fasttree params
    parser.add_argument('--seq-type', choices=['protein', 'nucleotide'],
                        help='序列类型(不指定自动检测)|Sequence type (auto-detect if not specified)')
    parser.add_argument('-t', '--threads', type=int, default=88,
                        help='MAFFT 线程数|MAFFT threads (default: 88)')
    parser.add_argument('--mafft-params', default='--auto',
                        help='MAFFT 额外参数|Additional MAFFT parameters (default: --auto)')
    parser.add_argument('--fasttree-params', default='',
                        help='FastTree 额外参数|Additional FastTree parameters')
    parser.add_argument('--mafft-path', default='mafft',
                        help='MAFFT 路径|MAFFT path (default: mafft)')
    parser.add_argument('--fasttree-path', default='fasttree',
                        help='FastTree 路径|FastTree path (default: fasttree)')

    parser.add_argument('--sample-name',
                        help='输出文件名前缀|Output filename prefix (default: 输入 basename|input basename)')
    parser.add_argument('--log-file',
                        help='日志文件路径|Log file path')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出|Verbose output')

    args = parser.parse_args()

    try:
        runner = PhyloTrimRunner(
            input_file=args.input,
            output_dir=args.output_dir,
            skip_trimal=args.skip_trimal,
            trimal_method=args.trimal_method,
            gt_threshold=args.gt_threshold,
            cons_threshold=args.cons_threshold,
            trimal_format=args.trimal_format,
            seq_type=args.seq_type,
            threads=args.threads,
            mafft_params=args.mafft_params,
            fasttree_params=args.fasttree_params,
            mafft_path=args.mafft_path,
            fasttree_path=args.fasttree_path,
            sample_name=args.sample_name,
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
