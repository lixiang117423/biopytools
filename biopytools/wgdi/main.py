"""
WGDI主程序模块|WGDI Main Module
"""

import argparse
import sys
from .config import WGDIConfig, DotPlotConfig, CollinearityConfig, CalKsConfig
from .utils import WGDILogger, CommandRunner
from .processor import WGDIProcessor


def register_dotplot(subparsers):
    """注册dotplot子命令|Register dotplot subcommand"""
    parser = subparsers.add_parser(
        'dotplot',
        help='同源基因点图|Homologous gene dotplot'
    )

    # 必需文件|Required files
    parser.add_argument('-b', '--blast', required=True,
                       help='BLAST结果文件|BLAST result file')
    parser.add_argument('--gff1', required=True,
                       help='物种1 GFF文件|Species 1 GFF file')
    parser.add_argument('--gff2', required=True,
                       help='物种2 GFF文件|Species 2 GFF file')
    parser.add_argument('--lens1', required=True,
                       help='物种1 LENS文件|Species 1 LENS file')
    parser.add_argument('--lens2', required=True,
                       help='物种2 LENS文件|Species 2 LENS file')

    # 可选参数|Optional parameters
    parser.add_argument('--genome1-name', default='Genome1',
                       help='物种1名称|Species 1 name (default: Genome1)')
    parser.add_argument('--genome2-name', default='Genome2',
                       help='物种2名称|Species 2 name (default: Genome2)')
    parser.add_argument('--blast-reverse', action='store_true',
                       help='BLAST结果前两列互换|Swap first two columns of BLAST result')
    parser.add_argument('--multiple', type=int, default=1,
                       help='最佳同源基因数|Number of best homologous genes (default: 1)')
    parser.add_argument('--score', type=int, default=100,
                       help='BLAST分数阈值|BLAST score threshold (default: 100)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                       help='BLAST e-value阈值|BLAST e-value threshold (default: 1e-5)')
    parser.add_argument('--repeat-number', type=int, default=10,
                       help='重复基因最大数量|Max number of repetitive genes (default: 10)')
    parser.add_argument('--position', choices=['order', 'start', 'end'],
                       default='order',
                       help='基因位置类型|Gene position type (default: order)')
    parser.add_argument('--ancestor-left',
                       help='左侧物种祖先染色体区域|Ancestral chromosome region for left species')
    parser.add_argument('--ancestor-top',
                       help='顶部物种祖先染色体区域|Ancestral chromosome region for top species')
    parser.add_argument('--markersize', type=float, default=0.5,
                       help='点大小|Marker size (default: 0.5)')
    parser.add_argument('--figsize', default='10,10',
                       help='图像尺寸|Figure size (default: 10,10)')
    parser.add_argument('--savefig', default='dotplot.png',
                       help='输出图像文件|Output image file (default: dotplot.png)')

    # 通用参数|Common parameters
    parser.add_argument('-o', '--output-dir', default='./wgdi_output',
                       help='输出目录|Output directory (default: ./wgdi_output)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='线程数|Number of threads (default: 8)')

    parser.set_defaults(func='dotplot')


def register_collinearity(subparsers):
    """注册collinearity子命令|Register collinearity subcommand"""
    parser = subparsers.add_parser(
        'collinearity',
        help='共线性分析|Collinearity analysis'
    )

    # 必需文件|Required files
    parser.add_argument('-b', '--blast', required=True,
                       help='BLAST结果文件|BLAST result file')
    parser.add_argument('--gff1', required=True,
                       help='物种1 GFF文件|Species 1 GFF file')
    parser.add_argument('--gff2', required=True,
                       help='物种2 GFF文件|Species 2 GFF file')
    parser.add_argument('--lens1', required=True,
                       help='物种1 LENS文件|Species 1 LENS file')
    parser.add_argument('--lens2', required=True,
                       help='物种2 LENS文件|Species 2 LENS file')

    # 处理参数|Processing parameters
    parser.add_argument('--blast-reverse', action='store_true',
                       help='BLAST结果前两列互换|Swap first two columns of BLAST result')
    parser.add_argument('--comparison', choices=['genomes', 'chromosomes'],
                       default='genomes',
                       help='比较模式|Comparison mode (default: genomes)')
    parser.add_argument('--multiple', type=int, default=1,
                       help='最佳同源基因数|Number of best homologous genes (default: 1)')
    parser.add_argument('--process', type=int, default=8,
                       help='进程数|Number of processes (default: 8)')
    parser.add_argument('--score', type=int, default=100,
                       help='BLAST分数阈值|BLAST score threshold (default: 100)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                       help='BLAST e-value阈值|BLAST e-value threshold (default: 1e-5)')
    parser.add_argument('--grading', default='50,40,25',
                       help='评分参数(红,蓝,灰)|Scoring parameters for red,blue,gray (default: 50,40,25)')
    parser.add_argument('--mg', default='40,40',
                       help='最大gap值|Max gap values (default: 40,40)')
    parser.add_argument('--pvalue', type=float, default=1.0,
                       help='P值阈值|P-value threshold (default: 1.0)')
    parser.add_argument('--repeat-number', type=int, default=20,
                       help='重复基因最大数量|Max number of repetitive genes (default: 20)')
    parser.add_argument('--position', choices=['order', 'start', 'end'],
                       default='order',
                       help='基因位置类型|Gene position type (default: order)')

    # 输出参数|Output parameters
    parser.add_argument('--savefile', default='collinearity.txt',
                       help='输出文件|Output file (default: collinearity.txt)')

    # 通用参数|Common parameters
    parser.add_argument('-o', '--output-dir', default='./wgdi_output',
                       help='输出目录|Output directory (default: ./wgdi_output)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='线程数|Number of threads (default: 8)')

    parser.set_defaults(func='collinearity')


def register_calks(subparsers):
    """注册calks子命令|Register calks subcommand"""
    parser = subparsers.add_parser(
        'calks',
        help='计算Ka/Ks|Calculate Ka/Ks'
    )

    # 必需文件|Required files
    parser.add_argument('-c', '--collinearity', required=True,
                       help='共线性分析结果文件|Collinearity analysis result file')
    parser.add_argument('--fasta1', required=True,
                       help='物种1 FASTA文件|Species 1 FASTA file')
    parser.add_argument('--fasta2', required=True,
                       help='物种2 FASTA文件|Species 2 FASTA file')

    # 输出参数|Output parameters
    parser.add_argument('--savefile', default='ks_results.txt',
                       help='输出文件|Output file (default: ks_results.txt)')

    # 通用参数|Common parameters
    parser.add_argument('-o', '--output-dir', default='./wgdi_output',
                       help='输出目录|Output directory (default: ./wgdi_output)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                       help='线程数|Number of threads (default: 8)')

    parser.set_defaults(func='calks')


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='WGDI工具封装 - 比较基因组学分析|WGDI Wrapper - Comparative Genomics Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # DotPlot分析|DotPlot analysis
  %(prog)s dotplot -b blast.txt --gff1 species1.gff --gff2 species2.gff --lens1 species1.lens --lens2 species2.lens

  # Collinearity分析|Collinearity analysis
  %(prog)s collinearity -b blast.txt --gff1 species1.gff --gff2 species2.gff --lens1 species1.lens --lens2 species2.lens

  # CalKs分析|CalKs analysis
  %(prog)s calks -c collinearity.txt --fasta1 species1.fa --fasta2 species2.fa
        '''
    )

    # 添加全局参数|Add global arguments
    parser.add_argument('--wgdi-path',
                       help='WGDI软件路径|WGDI software path')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')

    # 添加子命令|Add subcommands
    subparsers = parser.add_subparsers(dest='command', help='子命令|Subcommand')
    register_dotplot(subparsers)
    register_collinearity(subparsers)
    register_calks(subparsers)

    # 解析参数|Parse arguments
    args = parser.parse_args()

    # 检查是否提供了子命令|Check if subcommand is provided
    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        # 创建通用配置|Create common configuration
        wgdi_path = getattr(args, 'wgdi_path', None)
        if wgdi_path:
            common_config = WGDIConfig(wgdi_path=wgdi_path)
        else:
            common_config = WGDIConfig()

        # 覆盖通用参数|Override common parameters
        if hasattr(args, 'threads'):
            common_config.threads = args.threads
        if hasattr(args, 'output_dir'):
            common_config.working_dir = args.output_dir
            common_config.threads = getattr(args, 'threads', 8)

        # 创建日志器|Create logger
        logger_manager = WGDILogger()
        logger = logger_manager.get_logger()

        # 创建命令执行器|Create command runner
        cmd_runner = CommandRunner(logger, common_config.working_dir)

        # 创建处理器|Create processor
        processor = WGDIProcessor(common_config, logger, cmd_runner)

        # 根据子命令执行|Execute based on subcommand
        if args.func == 'dotplot':
            config = DotPlotConfig(
                blast_file=args.blast,
                gff1_file=args.gff1,
                gff2_file=args.gff2,
                lens1_file=args.lens1,
                lens2_file=args.lens2,
                genome1_name=args.genome1_name,
                genome2_name=args.genome2_name,
                blast_reverse=args.blast_reverse,
                multiple=args.multiple,
                score=args.score,
                evalue=args.evalue,
                repeat_number=args.repeat_number,
                position=args.position,
                ancestor_left=args.ancestor_left,
                ancestor_top=args.ancestor_top,
                markersize=args.markersize,
                figsize=args.figsize,
                savefig=args.savefig,
                threads=args.threads,
                working_dir=args.output_dir
            )
            success = processor.run_dotplot(config)

        elif args.func == 'collinearity':
            config = CollinearityConfig(
                blast_file=args.blast,
                gff1_file=args.gff1,
                gff2_file=args.gff2,
                lens1_file=args.lens1,
                lens2_file=args.lens2,
                blast_reverse=args.blast_reverse,
                comparison=args.comparison,
                multiple=args.multiple,
                process=args.process,
                evalue=args.evalue,
                score=args.score,
                grading=args.grading,
                mg=args.mg,
                pvalue=args.pvalue,
                repeat_number=args.repeat_number,
                position=args.position,
                savefile=args.savefile,
                threads=args.threads,
                working_dir=args.output_dir
            )
            success = processor.run_collinearity(config)

        elif args.func == 'calks':
            config = CalKsConfig(
                collinearity_file=args.collinearity,
                fasta1_file=args.fasta1,
                fasta2_file=args.fasta2,
                savefile=args.savefile,
                threads=args.threads,
                working_dir=args.output_dir
            )
            success = processor.run_calks(config)

        # 根据结果退出|Exit based on result
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
