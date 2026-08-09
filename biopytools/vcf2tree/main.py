"""
VCF转树主程序模块|VCF to Tree Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import Vcf2TreeConfig
from .utils import Vcf2TreeLogger, CommandRunner, check_dependencies
from .vcf_converter import VcfToFastaConverter
from .fasttree_builder import FastTreeBuilder
from .iqtree_builder import IqtreeBuilder


class Vcf2TreeRunner:
    """VCF转树运行器|VCF to Tree Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = Vcf2TreeConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = Vcf2TreeLogger(self.config.log_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化各组件|Initialize components
        self.vcf_converter = VcfToFastaConverter(self.config, self.logger)

        if self.config.method == 'fasttree':
            self.tree_builder = FastTreeBuilder(
                self.config, self.logger, self.cmd_runner
            )
        else:
            self.tree_builder = IqtreeBuilder(
                self.config, self.logger, self.cmd_runner
            )

        self.logger.info(
            "VCF转树工具已初始化|VCF to Tree tool initialized"
        )
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"建树方法|Tree method: {self.config.method}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

    def _is_step_completed(self, output_file: Path) -> bool:
        """检查步骤是否已完成|Check if step is completed"""
        return output_file.exists() and output_file.stat().st_size > 0

    def run_pipeline(self):
        """运行完整的VCF→树流程|Run complete VCF→Tree pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info(
                "开始VCF转树流程|Starting VCF to Tree pipeline"
            )
            self.logger.info("=" * 60)

            # 检查依赖|Check dependencies
            check_dependencies(self.config, self.logger)

            # 步骤1: VCF→FASTA|Step 1: VCF→FASTA
            if self._is_step_completed(self.config.snps_fa):
                self.logger.info(
                    f"跳过已完成步骤|Skipping completed step: VCF→FASTA "
                    f"(断点续传|checkpoint resume)"
                )
            else:
                if not self.vcf_converter.convert():
                    raise RuntimeError(
                        "VCF→FASTA转换失败|VCF→FASTA conversion failed"
                    )

            # 步骤2: 建树|Step 2: Build tree
            if self._is_step_completed(self.config.tree_nwk):
                self.logger.info(
                    f"跳过已完成步骤|Skipping completed step: "
                    f"建树|Tree building (断点续传|checkpoint resume)"
                )
            else:
                if not self.tree_builder.build():
                    raise RuntimeError(
                        f"建树失败|Tree building failed ({self.config.method})"
                    )

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info(
                "VCF转树流程完成|VCF to Tree pipeline completed"
            )
            self.logger.info("=" * 60)
            self.logger.info("输出文件|Output files:")
            self.logger.info(
                f"  - SNP FASTA对齐矩阵|SNP FASTA alignment: {self.config.snps_fa}"
            )
            self.logger.info(
                f"  - 系统发育树|Phylogenetic tree: {self.config.tree_nwk}"
            )
            self.logger.info(
                f"  - 日志文件|Log file: {self.logger_manager.log_file}"
            )
            self.logger.info(
                f"输出目录|Output directory: {self.config.output_dir}"
            )

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF转系统发育树工具|VCF to Phylogenetic Tree Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i variants.vcf -o tree_output
  %(prog)s -i variants.vcf -o tree_output --method fasttree
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入VCF文件路径|Input VCF file path')

    # 建树方法|Tree method（默认IQ-TREE|default IQ-TREE）
    parser.add_argument('--method',
                       choices=['fasttree', 'iqtree'],
                       default='iqtree',
                       help='建树方法|Tree method: fasttree or iqtree (default: iqtree)')

    # 输出参数|Output arguments
    parser.add_argument('-o', '--output-dir',
                       default='./vcf2tree_output',
                       help='输出目录|Output directory')

    # 运行参数|Run parameters
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('-g', '--outgroup',
                       default='',
                       help='外群样本名称|Outgroup sample name')
    parser.add_argument('--min-samples-locus', type=int, default=4,
                       help='位点最少样本数|Minimum samples per locus')

    # FastTree参数|FastTree parameters
    fasttree_group = parser.add_argument_group('FastTree参数|FastTree parameters')
    fasttree_group.add_argument('--fasttree-path', default='fasttree',
                               help='FastTree软件路径|FastTree software path')
    fasttree_group.add_argument('--fasttree-params', default='',
                               help='FastTree额外参数|Additional FastTree parameters')

    # IQ-TREE参数|IQ-TREE parameters
    iqtree_group = parser.add_argument_group('IQ-TREE参数|IQ-TREE parameters')
    iqtree_group.add_argument('--iqtree-path', default=None,
                             help='IQ-TREE软件路径|IQ-TREE software path')
    iqtree_group.add_argument('--iqtree-bootstrap', type=int, default=1000,
                             help='IQ-TREE UFBoot重复次数|IQ-TREE UFBoot replicates')
    iqtree_group.add_argument('--iqtree-model', default=None,
                             help='IQ-TREE进化模型(默认ModelFinder自动)|IQ-TREE model (default: ModelFinder)')

    args = parser.parse_args()

    # 构建配置|Build config
    kwargs = dict(
        input_file=args.input,
        method=args.method,
        output_dir=args.output_dir,
        threads=args.threads,
        outgroup=args.outgroup,
        min_samples_locus=args.min_samples_locus,
        fasttree_path=args.fasttree_path,
        fasttree_params=args.fasttree_params,
        iqtree_bootstrap=args.iqtree_bootstrap,
    )

    if args.iqtree_path:
        from .config import _expand_path
        kwargs['iqtree_path'] = _expand_path(args.iqtree_path)
    if args.iqtree_model:
        kwargs['iqtree_model'] = args.iqtree_model

    runner = Vcf2TreeRunner(**kwargs)
    runner.run_pipeline()


if __name__ == "__main__":
    main()
