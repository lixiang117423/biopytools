"""
VG主程序模块|VG Main Module
"""

import sys
import argparse
from pathlib import Path

from .config import (
    VGConstructConfig,
    VGIndexConfig,
    VGGiraffeConfig,
    VGDeconstructConfig
)
from .utils import VGLogger, validate_vg_environment, run_vg_command


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='VG变异图分析工具|VG Variation Graph Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
子命令|Subcommands:
  construct    从VCF和参考基因组构建变异图|Build variation graph from VCF and reference
  index        创建图索引|Create graph indexes
  giraffe      快速序列比对|Fast read alignment
  deconstruct  从变异图生成VCF|Export VCF from variation graph

示例|Examples:
  # 构建变异图|Build variation graph
  biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg

  # 创建索引|Create indexes
  biopytools vg index -i graph.vg -o graph_prefix --xg --giraffe

  # 序列比对|Align reads
  biopytools vg giraffe -g graph_prefix -f reads.fq -o alignments.gam

  # 导出VCF|Export VCF
  biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf
        '''
    )

    parser.add_argument('--vg-env',
                       default='vg_v.1.7.0',
                       help='VG conda环境名称|VG conda environment name (default: vg_v.1.7.0)')
    parser.add_argument('--log-level',
                       choices=['DEBUG', 'INFO', 'WARN', 'ERROR'],
                       default='INFO',
                       help='日志级别|Log level (default: INFO)')

    subparsers = parser.add_subparsers(dest='command', help='子命令|Subcommand')

    # ========== construct命令|construct command ==========
    construct_parser = subparsers.add_parser(
        'construct',
        help='构建变异图|Build variation graph',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法|Basic usage
  biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg

  # 指定染色体区域|Specify chromosome region
  biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg -R chr1

  # 保存alt等位基因路径|Save alt allele paths
  biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg --alt-paths
        '''
    )
    construct_parser.add_argument('-r', '--reference', required=True,
                                  help='参考基因组FASTA文件|Reference FASTA file')
    construct_parser.add_argument('-v', '--vcf', required=True,
                                  help='VCF文件|VCF file')
    construct_parser.add_argument('-o', '--output', required=True,
                                  help='输出VG文件|Output VG file')
    construct_parser.add_argument('-R', '--region',
                                  help='指定染色体区域|Specify chromosome region')
    construct_parser.add_argument('-t', '--threads', type=int, default=12,
                                  help='线程数|Number of threads (default: 12)')
    construct_parser.add_argument('--alt-paths', action='store_true',
                                  help='保存alt等位基因路径|Save alt allele paths')
    construct_parser.add_argument('--no-progress', action='store_true',
                                  help='不显示进度|Do not show progress')

    # ========== index命令|index command ==========
    index_parser = subparsers.add_parser(
        'index',
        help='创建图索引|Create graph indexes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 创建XG索引|Create XG index
  biopytools vg index -i graph.vg -o graph --xg

  # 创建GIRAFFE索引（用于比对）|Create GIRAFFE indexes (for alignment)
  biopytools vg index -i graph.vg -o graph --giraffe

  # 创建所有索引|Create all indexes
  biopytools vg index -i graph.vg -o graph --xg --gcsa --gbwt --giraffe
        '''
    )
    index_parser.add_argument('-i', '--input', required=True,
                             help='输入图文件|Input graph file')
    index_parser.add_argument('-o', '--output', required=True,
                             help='输出前缀|Output prefix')
    index_parser.add_argument('--xg', action='store_true',
                             help='创建XG索引|Create XG index')
    index_parser.add_argument('--gcsa', action='store_true',
                             help='创建GCSA索引|Create GCSA index')
    index_parser.add_argument('--gbwt', action='store_true',
                             help='创建GBWT索引|Create GBWT index')
    index_parser.add_argument('--giraffe', action='store_true',
                             help='创建GIRAFFE索引|Create GIRAFFE indexes')
    index_parser.add_argument('-k', '--kmer-size', type=int, default=16,
                             help='GCSA k-mer大小|GCSA k-mer size (default: 16)')
    index_parser.add_argument('-t', '--threads', type=int, default=12,
                             help='线程数|Number of threads (default: 12)')
    index_parser.add_argument('--no-progress', action='store_true',
                             help='不显示进度|Do not show progress')

    # ========== giraffe命令|giraffe command ==========
    giraffe_parser = subparsers.add_parser(
        'giraffe',
        help='快速序列比对|Fast read alignment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 单端测序比对|Single-end alignment
  biopytools vg giraffe -g graph -f reads.fq -o alignments.gam

  # 双端测序比对|Paired-end alignment
  biopytools vg giraffe -g graph -f reads1.fq -f2 reads2.fq -o alignments.gam

  # 输出GAF格式|Output GAF format
  biopytools vg giraffe -g graph -f reads.fq -o alignments.gaf --format GAF

  # 指定片段长度|Specify fragment length
  biopytools vg giraffe -g graph -f reads.fq -o alignments.gam -l 500 -s 50
        '''
    )
    giraffe_parser.add_argument('-g', '--graph', required=True,
                               help='图文件前缀（索引）|Graph file prefix (indexed)')
    giraffe_parser.add_argument('-f', '--reads', required=True,
                               help='输入reads文件|Input reads file')
    giraffe_parser.add_argument('-o', '--output', required=True,
                               help='输出GAM文件|Output GAM file')
    giraffe_parser.add_argument('-f2', '--reads2',
                               help='第二个reads文件（双端测序）|Second reads file (paired-end)')
    giraffe_parser.add_argument('-t', '--threads', type=int, default=12,
                               help='线程数|Number of threads (default: 12)')
    giraffe_parser.add_argument('-l', '--fragment-length', type=int, default=0,
                               help='片段长度（0=自动检测）|Fragment length (0=auto)')
    giraffe_parser.add_argument('-s', '--fragment-std-dev', type=int, default=0,
                               help='片段长度标准差（0=自动检测）|Fragment length std dev (0=auto)')
    giraffe_parser.add_argument('--min-identity', type=float, default=0.0,
                               help='最小相似度|Min identity (default: 0.0)')
    giraffe_parser.add_argument('--format', choices=['GAM', 'GAF'], default='GAM',
                               help='输出格式|Output format (default: GAM)')
    giraffe_parser.add_argument('--progress', action='store_true',
                               help='显示进度|Show progress')

    # ========== deconstruct命令|deconstruct command ==========
    deconstruct_parser = subparsers.add_parser(
        'deconstruct',
        help='从变异图生成VCF|Export VCF from variation graph',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法|Basic usage
  biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf

  # 指定样本|Specify samples
  biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf -s sample1 -s sample2
        '''
    )
    deconstruct_parser.add_argument('-i', '--input', required=True,
                                   help='输入图文件|Input graph file')
    deconstruct_parser.add_argument('-o', '--output', required=True,
                                   help='输出VCF文件|Output VCF file')
    deconstruct_parser.add_argument('-r', '--reference-path', required=True,
                                   help='参考路径名称|Reference path name')
    deconstruct_parser.add_argument('-s', '--samples', action='append',
                                   help='样本列表（可多次使用）|Sample list (can be used multiple times)')
    deconstruct_parser.add_argument('-t', '--threads', type=int, default=12,
                                   help='线程数|Number of threads (default: 12)')
    deconstruct_parser.add_argument('--no-progress', action='store_true',
                                   help='不显示进度|Do not show progress')

    return parser.parse_args()


def cmd_construct(args, logger):
    """执行construct命令|Execute construct command"""
    from .construct import VGConstructRunner

    config = VGConstructConfig(
        vg_env=args.vg_env,
        log_level=args.log_level,
        reference=args.reference,
        vcf=args.vcf,
        output=args.output,
        region=args.region,
        threads=args.threads,
        alt_paths=args.alt_paths,
        progress=not args.no_progress
    )

    config.validate()
    runner = VGConstructRunner(config, logger)
    return runner.run()


def cmd_index(args, logger):
    """执行index命令|Execute index command"""
    from .index import VGIndexRunner

    config = VGIndexConfig(
        vg_env=args.vg_env,
        log_level=args.log_level,
        input_graph=args.input,
        output_prefix=args.output,
        xg=args.xg,
        gcsa=args.gcsa,
        gbwt=args.gbwt,
        giraffe=args.giraffe,
        kmer_size=args.kmer_size,
        threads=args.threads,
        progress=not args.no_progress
    )

    config.validate()
    runner = VGIndexRunner(config, logger)
    return runner.run()


def cmd_giraffe(args, logger):
    """执行giraffe命令|Execute giraffe command"""
    from .giraffe import VGGiraffeRunner

    config = VGGiraffeConfig(
        vg_env=args.vg_env,
        log_level=args.log_level,
        graph=args.graph,
        reads=args.reads,
        output=args.output,
        reads2=args.reads2,
        threads=args.threads,
        fragment_length=args.fragment_length,
        fragment_std_dev=args.fragment_std_dev,
        min_identity=args.min_identity,
        output_format=args.format,
        progress=args.progress
    )

    config.validate()
    runner = VGGiraffeRunner(config, logger)
    return runner.run()


def cmd_deconstruct(args, logger):
    """执行deconstruct命令|Execute deconstruct command"""
    from .deconstruct import VGDeconstructRunner

    config = VGDeconstructConfig(
        vg_env=args.vg_env,
        log_level=args.log_level,
        input_graph=args.input,
        output=args.output,
        reference_path=args.reference_path,
        samples=args.samples,
        threads=args.threads,
        progress=not args.no_progress
    )

    config.validate()
    runner = VGDeconstructRunner(config, logger)
    return runner.run()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 检查是否指定了子命令|Check if subcommand is specified
        if not args.command:
            print("错误|Error: 必须指定子命令|Must specify a subcommand")
            print("使用|Use: biopytools vg --help 查看帮助|to see help")
            sys.exit(1)

        # 设置日志|Setup logging
        log_file = Path(f"vg_{args.command}.log")
        logger_manager = VGLogger(log_file, args.log_level)
        logger = logger_manager.get_logger()

        # 打印头部信息|Print header information
        logger.info("")
        logger.info("=" * 60)
        logger.info(f"VG {args.command.upper()} 分析|VG {args.command.upper()} Analysis")
        logger.info("=" * 60)

        # 验证环境|Validate environment
        if not validate_vg_environment(args.vg_env, logger):
            logger.error("VG环境验证失败|VG environment validation failed")
            sys.exit(1)

        # 根据子命令执行相应功能|Execute corresponding function based on subcommand
        if args.command == 'construct':
            success = cmd_construct(args, logger)
        elif args.command == 'index':
            success = cmd_index(args, logger)
        elif args.command == 'giraffe':
            success = cmd_giraffe(args, logger)
        elif args.command == 'deconstruct':
            success = cmd_deconstruct(args, logger)
        else:
            logger.error(f"未知子命令|Unknown subcommand: {args.command}")
            sys.exit(1)

        # 打印结果|Print result
        if success:
            logger.info("")
            logger.info("=" * 60)
            logger.info("分析成功完成|Analysis completed successfully")
            logger.info("=" * 60)
            sys.exit(0)
        else:
            logger.error("")
            logger.error("=" * 60)
            logger.error("分析失败|Analysis failed")
            logger.error("=" * 60)
            sys.exit(1)

    except KeyboardInterrupt:
        sys.stderr.write("用户中断操作|User interrupted operation\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"程序错误|Program error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()
