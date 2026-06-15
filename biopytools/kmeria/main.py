"""KMERIA主程序模块|KMERIA Main Module"""

import argparse
import sys
from pathlib import Path
from .config import (
    PipelineConfig,
    CountConfig,
    KctmConfig,
    FilterConfig,
    M2bConfig,
    AssoConfig
)
from .utils import KMERIALogger
from .processors import PipelineProcessor
from .processors.count import CountProcessor
from .processors.kctm import KctmProcessor
from .processors.filter import FilterProcessor
from .processors.m2b import M2bProcessor
from .processors.asso import AssoProcessor


class KMERIARunner:
    """KMERIA运行器|KMERIA Runner"""

    def __init__(self, command, **kwargs):
        self.command = command
        self.kwargs = kwargs

    def _filter_kwargs(self, config_class):
        """
        过滤kwargs，只保留config支持的参数|Filter kwargs to keep only config-supported parameters

        Args:
            config_class: 配置类|Config class

        Returns:
            过滤后的kwargs字典|Filtered kwargs dict
        """
        from dataclasses import fields

        # 参数名映射|Parameter name mapping
        param_mapping = {
            'samples': 'samples_file',
            'pheno_file': 'pheno_file',  # 保持原样|Keep as is
            'depth_file': 'depth_file',   # 保持原样|Keep as is
            # AssoConfig 映射|AssoConfig mapping
            'kin_method': 'kinship_method',
            'maf': 'minor_allele_freq',
            'miss': 'missing_threshold',
        }

        # 映射参数名|Map parameter names
        filtered_kwargs = {}
        for key, value in self.kwargs.items():
            # 处理特殊的布尔标志|Handle special boolean flags
            if key == 'no_normalize' and config_class == M2bConfig:
                filtered_kwargs['normalize'] = False
                continue
            elif key == 'quantile_norm' and config_class == M2bConfig:
                filtered_kwargs['quantile_norm'] = value
                continue
            elif key == 'use_gemma' and config_class == AssoConfig:
                # --use-gemma → use_bimbam_tools=False
                if value:
                    filtered_kwargs['use_bimbam_tools'] = False
                continue
            elif key == 'no_kinship' and config_class == AssoConfig:
                # --no-kinship → use_kinship=False
                if value:
                    filtered_kwargs['use_kinship'] = False
                continue

            new_key = param_mapping.get(key, key)
            filtered_kwargs[new_key] = value

        # 获取支持的参数|Get supported parameters
        supported_params = {f.name for f in fields(config_class)}

        # 过滤只保留支持的参数|Filter to keep only supported parameters
        final_kwargs = {k: v for k, v in filtered_kwargs.items() if k in supported_params}

        return final_kwargs

    def run(self):
        """运行KMERIA分析|Run KMERIA analysis"""
        # 初始化日志|Initialize logging
        log_file = self.kwargs.get('log_file')
        logger_manager = KMERIALogger(log_file)
        logger = logger_manager.get_logger()

        try:
            if self.command == 'pipeline':
                return self._run_pipeline(logger)
            elif self.command == 'count':
                return self._run_count(logger)
            elif self.command == 'kctm':
                return self._run_kctm(logger)
            elif self.command == 'filter':
                return self._run_filter(logger)
            elif self.command == 'm2b':
                return self._run_m2b(logger)
            elif self.command == 'asso':
                return self._run_asso(logger)
            else:
                logger.error(f"未知命令|Unknown command: {self.command}")
                return False

        except Exception as e:
            logger.error(f"执行出错|Execution error: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False

    def _run_pipeline(self, logger):
        """运行完整流程|Run complete pipeline"""
        from .utils import CommandRunner

        config = PipelineConfig(**self._filter_kwargs(PipelineConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = PipelineProcessor(config, logger, cmd_runner)

        return processor.run()

    def _run_count(self, logger):
        """运行k-mer计数|Run k-mer counting"""
        from .utils import CommandRunner

        config = CountConfig(**self._filter_kwargs(CountConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = CountProcessor(config, logger, cmd_runner)

        return processor.run()

    def _run_kctm(self, logger):
        """运行矩阵构建|Run matrix construction"""
        from .utils import CommandRunner

        config = KctmConfig(**self._filter_kwargs(KctmConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = KctmProcessor(config, logger, cmd_runner)

        return processor.run()

    def _run_filter(self, logger):
        """运行过滤|Run filtering"""
        from .utils import CommandRunner

        config = FilterConfig(**self._filter_kwargs(FilterConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = FilterProcessor(config, logger, cmd_runner)

        return processor.run()

    def _run_m2b(self, logger):
        """运行格式转换|Run format conversion"""
        from .utils import CommandRunner

        config = M2bConfig(**self._filter_kwargs(M2bConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = M2bProcessor(config, logger, cmd_runner)

        return processor.run()

    def _run_asso(self, logger):
        """运行关联分析|Run association analysis"""
        from .utils import CommandRunner

        config = AssoConfig(**self._filter_kwargs(AssoConfig))
        config.validate()

        cmd_runner = CommandRunner(logger, config.output_dir)
        processor = AssoProcessor(config, logger, cmd_runner)

        return processor.run()


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="K-mer GWAS分析工具|K-mer GWAS Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='子命令|Sub-command')

    # ========== pipeline 子命令 ==========
    pipeline_parser = subparsers.add_parser('pipeline', help='完整分析流程|Complete pipeline')

    pipeline_parser.add_argument('-i', '--fastq-dir', required=True,
                                help='FASTQ文件目录|FASTQ files directory')
    pipeline_parser.add_argument('--sample', '--samples', dest='samples', required=True,
                                help='样本列表文件|Sample list file')
    pipeline_parser.add_argument('-d', '--depth-file', required=True,
                                help='测序深度文件|Sequencing depth file')
    pipeline_parser.add_argument('-p', '--pheno-file', required=True,
                                help='表型文件|Phenotype file')
    pipeline_parser.add_argument('-o', '--output-dir', default='./kmeria_results',
                                help='输出目录|Output directory')

    pipeline_parser.add_argument('-k', '--kmer-size', type=int, default=31,
                                help='K-mer大小|K-mer size (default: 31)')
    pipeline_parser.add_argument('--max-abund', type=int, default=1000,
                                help='最大丰度|Maximum abundance (default: 1000)')
    pipeline_parser.add_argument('--missing-ratio', type=float, default=0.8,
                                help='缺失率（k-mer默认0.8，与GitHub一致）|Missing ratio (k-mer default 0.8, matches GitHub)')
    pipeline_parser.add_argument('--ploidy', type=int, default=4,
                                help='倍性（默认4，适配多倍体）|Ploidy (default 4, for polyploids)')

    pipeline_parser.add_argument('--step', choices=['count', 'kctm', 'filter', 'm2b', 'asso'],
                                help='从指定步骤开始|Start from specified step')
    pipeline_parser.add_argument('-t', '--threads', type=int, default=24,
                                help='线程数|Thread count (default: 24)')
    pipeline_parser.add_argument('--batch-size', type=int, default=4,
                                help='批处理大小|Batch size (default: 4)')

    pipeline_parser.add_argument('--pheno-col', type=int, default=1,
                                help='表型列|Phenotype column (default: 1)')
    pipeline_parser.add_argument('--kinship-file',
                                help='亲缘关系矩阵|Kinship matrix file')
    pipeline_parser.add_argument('--covar-file',
                                help='协变量文件|Covariate file')

    pipeline_parser.add_argument('--enable-qc', action='store_true', default=True,
                                help='启用质控|Enable QC (default: True)')

    pipeline_parser.add_argument('--genome-file',
                                help='参考基因组|Reference genome (for Post-GWAS analysis)')
    pipeline_parser.add_argument('--gff-file',
                                help='GFF注释文件|GFF annotation file (optional)')
    pipeline_parser.add_argument('--sample-ratio', type=float, default=0.1,
                                help='高p值位点抽样比例|Sampling ratio for high p-value loci (default: 0.1)')
    pipeline_parser.add_argument('--window-size', type=int, default=200000,
                                help='基因查找窗口大小|Gene search window size (default: 200000 = 200kb)')

    # 比对工具选择参数|Alignment tool selection parameters
    pipeline_parser.add_argument('--alignment-tool', type=str, default='bwa', choices=['bwa', 'blast'],
                                help='Post-GWAS比对工具选择 (默认: bwa)|Alignment tool for Post-GWAS (default: bwa)')
    pipeline_parser.add_argument('--bwa-k', type=int, default=9,
                                help='BWA mem -k 参数，最小种子长度 (默认: 9)|BWA mem -k parameter, minimum seed length (default: 9)')
    pipeline_parser.add_argument('--bwa-t-min-score', type=int, default=10, dest='bwa_T',
                                help='BWA mem -T 参数，最小输出分数 (默认: 10)|BWA mem -T parameter, minimum score to output (default: 10)')
    pipeline_parser.add_argument('--as-ratio', type=float, default=0.95,
                                help='BWA AS过滤阈值，保留AS >= 最高AS * ratio的所有比对 (默认: 0.95)|BWA AS filtering threshold, keep alignments with AS >= max_AS * ratio (default: 0.95)')

    pipeline_parser.add_argument('--log-file',
                                help='日志文件|Log file')
    pipeline_parser.add_argument('-f', '--force', action='store_true',
                                help='强制重新运行所有步骤|Force re-run all steps')

    # ========== count 子命令 ==========
    count_parser = subparsers.add_parser('count', help='k-mer计数|K-mer counting')

    count_parser.add_argument('-i', '--fastq-dir', required=True,
                             help='FASTQ文件目录|FASTQ files directory')
    count_parser.add_argument('--sample', '--samples', dest='samples', required=True,
                             help='样本列表文件|Sample list file')
    count_parser.add_argument('-o', '--output-dir', default='./01_kmer_counts',
                             help='输出目录|Output directory')

    count_parser.add_argument('-k', '--kmer-size', type=int, default=31,
                             help='K-mer大小|K-mer size (default: 31)')
    count_parser.add_argument('-t', '--threads', type=int, default=24,
                             help='线程数|Thread count (default: 24)')
    count_parser.add_argument('-b', '--batch-size', type=int, default=4,
                             help='批处理大小|Batch size (default: 4)')
    count_parser.add_argument('-C', '--count-separate-strands', action='store_true',
                             help='分别计数链|Count strands separately')
    count_parser.add_argument('-T', '--text-output', action='store_true',
                             help='文本输出|Text output')
    count_parser.add_argument('--log-file',
                             help='日志文件|Log file')
    count_parser.add_argument('-f', '--force', action='store_true',
                             help='强制重新运行|Force re-run')

    # ========== kctm 子命令 ==========
    kctm_parser = subparsers.add_parser('kctm', help='k-mer矩阵构建|K-mer matrix construction')

    kctm_parser.add_argument('-i', '--input-dir', required=True,
                            help='输入目录|Input directory')
    kctm_parser.add_argument('-o', '--output-dir', default='./02_kmer_matrices',
                            help='输出目录|Output directory')
    kctm_parser.add_argument('-t', '--threads', type=int, default=24,
                            help='线程数|Thread count (default: 24)')
    kctm_parser.add_argument('--log-file',
                            help='日志文件|Log file')

    # ========== filter 子命令 ==========
    filter_parser = subparsers.add_parser('filter', help='k-mer过滤|K-mer filtering')

    filter_parser.add_argument('-i', '--input-dir', required=True,
                              help='输入目录|Input directory')
    filter_parser.add_argument('-o', '--output-dir', default='./03_filtered_matrices',
                              help='输出目录|Output directory')
    filter_parser.add_argument('-d', '--depth-file', required=True,
                              help='测序深度文件|Sequencing depth file')
    filter_parser.add_argument('-c', '--max-abund', type=int, default=1000,
                              help='最大丰度|Maximum abundance (default: 1000)')
    filter_parser.add_argument('-s', '--missing-ratio', type=float, default=0.8,
                              help='缺失率（默认0.8，与GitHub一致）|Missing ratio (default 0.8, matches GitHub)')
    filter_parser.add_argument('-p', '--ploidy', type=int, default=4,
                              help='倍性（默认4）|Ploidy (default 4)')
    filter_parser.add_argument('-t', '--threads', type=int, default=24,
                              help='线程数|Thread count (default: 24)')
    filter_parser.add_argument('--log-file',
                              help='日志文件|Log file')

    # ========== m2b 子命令 ==========
    m2b_parser = subparsers.add_parser('m2b', help='转换为BIMBAM格式|Convert to BIMBAM format')

    m2b_parser.add_argument('--in', '-i', required=True, dest='input_dir',
                           help='输入目录|Input directory')
    m2b_parser.add_argument('--out', '-o', default='./04_bimbam', dest='output_dir',
                           help='输出目录|Output directory')
    m2b_parser.add_argument('-t', '--threads', type=int, default=24,
                           help='线程数|Thread count (default: 24)')
    m2b_parser.add_argument('--no-normalize', action='store_true',
                           help='不归一化|No normalization')
    m2b_parser.add_argument('--quantile-norm', action='store_true',
                           help='分位数归一化|Quantile normalization')
    m2b_parser.add_argument('--log-file',
                           help='日志文件|Log file')

    # ========== asso 子命令 ==========
    asso_parser = subparsers.add_parser('asso', help='k-mer关联分析|K-mer association analysis')

    asso_parser.add_argument('-i', '--input-dir', required=True,
                            help='输入目录|Input directory')
    asso_parser.add_argument('-p', '--pheno-file', required=True,
                            help='表型文件|Phenotype file')
    asso_parser.add_argument('-o', '--output-dir', default='./05_association',
                            help='输出目录|Output directory')
    asso_parser.add_argument('-n', '--pheno-col', type=int, default=1,
                            help='表型列|Phenotype column (default: 1)')
    asso_parser.add_argument('-c', '--covar-file',
                            help='协变量文件|Covariate file')
    asso_parser.add_argument('-k', '--kinship-file',
                            help='亲缘关系矩阵|Kinship matrix file')
    asso_parser.add_argument('-t', '--threads', type=int, default=64,
                            help='线程数|Thread count (default: 64)')

    # 工具与方法|Tool and method
    asso_parser.add_argument('--use-gemma', action='store_true',
                            help='使用gemma替代bimbamAsso|Use gemma instead of bimbamAsso')
    asso_parser.add_argument('-m', '--analysis-method',
                            help='分析方法（default/lm/lmm）|Analysis method')
    asso_parser.add_argument('--kin-method', type=int, default=3,
                            help='kinship计算方法（1=IBS均值, 2=IBS随机, 3=Balding-Nichols）|Kinship method (default: 3)')
    asso_parser.add_argument('--no-kinship', action='store_true',
                            help='不使用kinship矩阵|Do not use kinship matrix')
    asso_parser.add_argument('--disable-gls', action='store_true',
                            help='禁用GLS，使用OLS|Disable GLS, use OLS')
    asso_parser.add_argument('--write-eigen', action='store_true',
                            help='输出特征值/特征向量|Write eigenvalues/eigenvectors')

    # 精度|Precision
    asso_parser.add_argument('--kin-precision', type=int, default=10,
                            help='kinship精度|Kinship precision (default: 10)')
    asso_parser.add_argument('--out-precision', type=int, default=5,
                            help='输出精度|Output precision (default: 5)')

    # 质量控制|Quality control
    asso_parser.add_argument('--maf', type=float,
                            help='次等位基因频率过滤|Minor allele frequency filter')
    asso_parser.add_argument('--miss', type=float,
                            help='缺失阈值|Missing threshold')

    # 分块分析|Block analysis
    asso_parser.add_argument('--start-marker', type=int,
                            help='起始marker索引|Start marker index')
    asso_parser.add_argument('--end-marker', type=int,
                            help='结束marker索引|End marker index')

    # 其他选项|Other options
    asso_parser.add_argument('--generate-plots', action='store_true',
                            help='生成图表|Generate plots')
    asso_parser.add_argument('--compress', action='store_true',
                            help='压缩输出|Compress output')
    asso_parser.add_argument('--verbose', action='store_true',
                            help='详细输出|Verbose output')
    asso_parser.add_argument('--dry-run', action='store_true',
                            help='仅显示命令不执行|Show commands without executing')
    asso_parser.add_argument('--no-validate', action='store_true',
                            help='跳过输入验证|Skip input validation')
    asso_parser.add_argument('--no-check-deps', action='store_true',
                            help='跳过依赖检查|Skip dependency check')
    asso_parser.add_argument('--no-cleanup', action='store_true',
                            help='保留临时文件|Keep temporary files')

    asso_parser.add_argument('--log-file',
                            help='日志文件|Log file')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    if not args.command:
        print("请指定子命令|Please specify a sub-command")
        print("使用 -h|--help 查看帮助|Use -h|--help for help")
        sys.exit(1)

    # 转换args为字典|Convert args to dictionary
    kwargs = vars(args)
    command = kwargs.pop('command')  # 移除command避免重复传递|Remove command to avoid duplicate

    # 创建runner并运行|Create runner and run
    runner = KMERIARunner(command, **kwargs)

    success = runner.run()

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
