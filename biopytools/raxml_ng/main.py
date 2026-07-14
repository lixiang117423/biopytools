"""
RAxML-NG 系统发育树分析主模块|RAxML-NG Phylogenetic Tree Analysis Main Module
"""

import argparse
import sys
from datetime import datetime

from ..common.paths import expand_path
from .config import RaxmlNgConfig
from .utils import RaxmlNgLogger, CommandRunner, check_dependencies, build_conda_command


class RaxmlNgRunner:
    """RAxML-NG 系统发育树分析主类|Main RAxML-NG Phylogenetic Tree Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置(含 validate)|Init config (with validate)
        self.config = RaxmlNgConfig(**kwargs)
        self.config.validate()

        # 初始化日志(写入 99_logs)|Init logging (to 99_logs)
        self.logger_manager = RaxmlNgLogger(self.config.logs_path)
        self.logger = self.logger_manager.get_logger()

        # 命令执行器(工作目录=输出目录,RAxML-NG 把 {prefix}.raxml.* 写进 cwd)|Command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger) is not None

    def _build_command(self) -> list:
        """按 mode 构建 RAxML-NG 命令列表|Build RAxML-NG command list by mode"""
        cfg = self.config

        # 通用参数|Common args
        args = ['--prefix', cfg.prefix, '--threads', str(cfg.threads)]
        if cfg.seed is not None:
            args.extend(['--seed', str(cfg.seed)])
        if cfg.redo:
            args.append('--redo')

        if cfg.mode == 'all':
            args = ['--all'] + args
            args.extend(['--msa', cfg.input_file])
            if cfg.model:
                args.extend(['--model', cfg.model])
            args.extend(['--bs-trees', cfg.bs_trees])
            if cfg.bs_metric != 'fbp':
                args.extend(['--bs-metric', cfg.bs_metric])
            if cfg.outgroup:
                args.extend(['--outgroup', cfg.outgroup])
        elif cfg.mode == 'search':
            args = ['--search'] + args
            args.extend(['--msa', cfg.input_file])
            if cfg.model:
                args.extend(['--model', cfg.model])
            if cfg.tree:
                args.extend(['--tree', cfg.tree])
            if cfg.outgroup:
                args.extend(['--outgroup', cfg.outgroup])
        elif cfg.mode == 'support':
            args = ['--support'] + args
            args.extend(['--tree', cfg.tree])
            args.extend(['--bs-trees', cfg.bs_trees])
            if cfg.bs_metric != 'fbp':
                args.extend(['--bs-metric', cfg.bs_metric])

        # 用 build_conda_command 包装(静态二进制直接调用)|Wrap with build_conda_command
        return build_conda_command(cfg.raxml_ng_path, args)

    def _write_software_versions(self, start_time: datetime, end_time: datetime, version_str: str):
        """生成 00_pipeline_info/software_versions.yml|Generate software_versions.yml"""
        try:
            import yaml
        except ImportError:
            self.logger.warning("缺少 pyyaml,跳过 software_versions.yml|Missing pyyaml, skip software_versions.yml")
            return

        info = {
            'pipeline': {
                'name': 'biopytools raxml_ng',
                'version': '1.0.0',
            },
            'tools': {
                'raxml-ng': {
                    'version': version_str,
                    'path': self.config.raxml_ng_path,
                    'command': f"{self.config.raxml_ng_path} --version",
                },
            },
            'parameters': {
                'mode': self.config.mode,
                'model': self.config.model or 'auto',
                'prefix': self.config.prefix,
                'threads': self.config.threads,
                'bs_trees': self.config.bs_trees,
                'bs_metric': self.config.bs_metric,
                'tree': self.config.tree,
                'outgroup': self.config.outgroup,
                'seed': self.config.seed,
                'redo': self.config.redo,
            },
            'execution': {
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': int((end_time - start_time).total_seconds()),
            },
        }

        out_file = self.config.info_path / 'software_versions.yml'
        with open(out_file, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"版本信息已写入|Version info written: {out_file}")

    def run(self) -> bool:
        """运行分析流程|Run analysis pipeline"""
        cfg = self.config
        start_time = datetime.now()

        self.logger.info("=" * 60)
        self.logger.info("RAxML-NG 系统发育树分析流程启动|RAxML-NG Phylogenetic Analysis Pipeline Started")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input file: {cfg.input_file}")
        self.logger.info(f"输出目录|Output dir: {cfg.output_dir}")
        self.logger.info(f"前缀|Prefix: {cfg.prefix}")
        self.logger.info(f"模式|Mode: {cfg.mode}")
        self.logger.info(f"线程数|Threads: {cfg.threads}")
        if cfg.model:
            self.logger.info(f"模型|Model: {cfg.model}")
        else:
            self.logger.info("模型|Model: 自适应选择|auto (adaptive)")

        # 步骤1: 依赖检查(并取版本)|Step 1: dependency check (and version)
        self.logger.info("步骤1: 检查软件依赖|Step 1: Checking software dependencies")
        version_str = check_dependencies(self.config, self.logger)
        if version_str is None:
            self.logger.error("依赖检查失败|Dependency check failed")
            sys.exit(1)

        # 步骤2: 执行建树(RAxML-NG 原生断点续传:不传 --redo 即续传)|Step 2: run tree analysis
        self.logger.info(
            f"步骤2: 执行 RAxML-NG {cfg.mode} 分析|Step 2: Running RAxML-NG {cfg.mode} analysis"
        )
        cmd = self._build_command()
        if not self.cmd_runner.run(cmd, description=f"RAxML-NG {cfg.mode}"):
            self.logger.error(
                f"RAxML-NG 执行失败,详见|RAxML-NG failed, see: {cfg.prefix}.raxml.log"
            )
            self._write_software_versions(start_time, datetime.now(), version_str)
            sys.exit(1)

        # 步骤3: 版本信息|Step 3: version info
        end_time = datetime.now()
        self._write_software_versions(start_time, end_time, version_str)

        self.logger.info("=" * 60)
        self.logger.info("RAxML-NG 分析流程完成|RAxML-NG Analysis Pipeline Completed")
        self.logger.info("=" * 60)
        self.logger.info(f"结果保存在|Results saved in: {cfg.output_dir}")
        self.logger.info(f"最佳树|Best tree: {cfg.prefix}.raxml.bestTree")
        return True


def main():
    """命令行入口|CLI entry"""
    parser = argparse.ArgumentParser(
        description='RAxML-NG 系统发育树分析脚本|RAxML-NG Phylogenetic Tree Analysis Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i alignment.fasta -o tree_results -p my_tree
        """,
    )

    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input', required=True,
                          help='输入比对文件 (FASTA/PHYLIP)|Input alignment file (FASTA/PHYLIP)')
    required.add_argument('-o', '--output-dir', required=True,
                          help='输出目录|Output directory')

    core = parser.add_argument_group('核心参数|Core parameters')
    core.add_argument('-p', '--prefix', default=None,
                      help='输出文件前缀 (默认输入文件名)|Output prefix (default: input filename)')
    core.add_argument('--mode', choices=['all', 'search', 'support'], default='all',
                      help='分析模式 (默认: all)|Analysis mode (default: all)')
    core.add_argument('-m', '--model', default=None,
                      help='进化模型 (不指定则自适应)|Evolutionary model (auto if not specified)')
    core.add_argument('-t', '--threads', type=int, default=12,
                      help='线程数 (默认: 12)|Number of threads (default: 12)')

    boot = parser.add_argument_group('Bootstrap参数|Bootstrap parameters')
    boot.add_argument('-b', '--bs-trees', default='1000',
                      help='Bootstrap重复数或autoMRE{N} (默认: 1000)|Bootstrap replicates or autoMRE{N}')
    boot.add_argument('--bs-metric', default='fbp',
                      help='支持值类型: fbp/tbe/sh/ebg/rbs (默认: fbp)|Branch support metric')

    tree_opts = parser.add_argument_group('树选项|Tree options')
    tree_opts.add_argument('--tree', default=None,
                           help='起始树(search)/参考树(support)|Starting tree (search) / reference tree (support)')
    tree_opts.add_argument('--outgroup', default=None,
                           help='外群名称 (逗号分隔)|Outgroup taxon names (comma-separated)')

    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument('--seed', type=int, default=None, help='随机种子|Random seed')
    other.add_argument('--redo', action='store_true',
                       help='覆盖已有结果,忽略checkpoint|Overwrite results, ignore checkpoint')
    other.add_argument('--raxml-ng-path', default=None, help='RAxML-NG程序路径|RAxML-NG program path')

    args = parser.parse_args()

    kwargs = dict(
        input_file=args.input,
        output_dir=args.output_dir,
        prefix=args.prefix,
        mode=args.mode,
        model=args.model,
        threads=args.threads,
        bs_trees=args.bs_trees,
        bs_metric=args.bs_metric,
        tree=args.tree,
        outgroup=args.outgroup,
        seed=args.seed,
        redo=args.redo,
    )

    if args.raxml_ng_path:
        kwargs['raxml_ng_path'] = expand_path(args.raxml_ng_path)

    try:
        runner = RaxmlNgRunner(**kwargs)
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
