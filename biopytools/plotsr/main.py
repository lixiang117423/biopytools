"""
PlotSR主程序模块|PlotSR Main Module
"""

import argparse
import sys
import os
import glob
from pathlib import Path

from .config import PlotSRConfig
from .utils import PlotSRLogger, discover_genomes_in_folder
from .pipeline import PlotSRPipeline


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="多基因组共线性可视化工具|Multi-Genome Synteny Visualization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""示例|Examples: biopytools plotsr -i genome1.fa -i genome2.fa -o output/

Map文件格式|Map file format (支持指定名称和顺序|Support specifying names and order):
  name1\\tpath/to/genome1.fa
  name2\\tpath/to/genome2.fa

使用map文件|Using map file: biopytools plotsr -i genomes.map -o output/"""
    )

    parser.add_argument(
        '-i', '--input',
        action='append',
        dest='genomes',
        required=True,
        help='输入基因组FASTA文件（可多次使用）或包含基因组的文件夹|'
             'Input genome FASTA files (can be used multiple times) or folder containing genomes'
    )

    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='输出目录|Output directory'
    )

    parser.add_argument(
        '-n', '--names',
        help='基因组名称（逗号分隔）|Genome names (comma-separated, e.g., Col-0,Ler,Cvi)'
    )

    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads [default: 12]'
    )

    parser.add_argument(
        '--minimap2-preset',
        default='asm5',
        choices=['asm5', 'asm10', 'asm20'],
        help='minimap2预设参数|minimap2 preset [default: asm5]'
    )

    parser.add_argument(
        '-s', '--min-sr-size',
        type=int,
        default=10000,
        help='最小结构变异大小|Minimum structural variant size [default: 10000]'
    )

    parser.add_argument(
        '--output-format',
        default='pdf',
        choices=['pdf', 'png', 'svg'],
        help='输出格式|Output format [default: pdf]'
    )

    # PlotSR可视化参数|PlotSR visualization parameters
    parser.add_argument(
        '-f', '--font-size',
        type=int,
        default=6,
        help='字体大小|Font size [default: 6]'
    )

    parser.add_argument(
        '-d', '--dpi',
        type=int,
        default=300,
        help='图片DPI|Image DPI [default: 300]'
    )

    parser.add_argument(
        '--space-ratio',
        type=float,
        default=0.7,
        help='同源染色体间距(0.1-0.75)|Space for homologous chromosomes [default: 0.7]'
    )

    parser.add_argument(
        '-v', '--vertical',
        action='store_true',
        help='垂直排列染色体|Plot vertical chromosomes'
    )

    parser.add_argument(
        '--itx',
        action='store_true',
        help='染色体间交互模式|Inter-chromosomal plotting mode'
    )

    # 过滤参数|Filtering parameters
    parser.add_argument(
        '--nosyn',
        action='store_true',
        help='不绘制同源区域|Do not plot syntenic regions'
    )

    parser.add_argument(
        '--noinv',
        action='store_true',
        help='不绘制倒位|Do not plot inversions'
    )

    parser.add_argument(
        '--notr',
        action='store_true',
        help='不绘制易位|Do not plot translocations'
    )

    parser.add_argument(
        '--nodup',
        action='store_true',
        help='不绘制重复|Do not plot duplications'
    )

    # 染色体过滤参数|Chromosome filtering parameters
    parser.add_argument(
        '-c', '--chr',
        action='append',
        dest='chromosomes',
        help='指定要显示的染色体（可多次使用，支持数字如1或名称如Chr1）|'
             'Specify chromosomes to display (can be used multiple times, supports number like 1 or name like Chr1)'
    )

    # 流程控制参数|Pipeline control parameters
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        default=True,
        help='跳过已完成的步骤（默认启用）|Skip completed steps (default: enabled)'
    )

    parser.add_argument(
        '--force-run',
        action='store_false',
        dest='skip_existing',
        help='强制重新运行所有步骤|Force re-run all steps'
    )

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )

    return parser.parse_args()


def prepare_genomes(input_list: list) -> tuple:
    """
    处理输入参数，返回基因组文件列表和名称列表|Process input parameters, return genome files and names

    Args:
        input_list: 输入路径列表|List of input paths

    Returns:
        tuple: (基因组文件列表, 基因组名称列表)|(List of genome files, List of genome names)
    """
    genomes = []
    names = []

    for item in input_list:
        if os.path.isfile(item):
            # 检查是否是map文件|Check if it's a map file
            if item.endswith('.map') or item.endswith('.txt'):
                genome_list, name_list = parse_map_file(item)
                genomes.extend(genome_list)
                names.extend(name_list)
            else:
                # 直接是文件|Is a file
                genomes.append(item)
        elif os.path.isdir(item):
            # 是文件夹，自动发现|Is a folder, auto-discover
            discovered = discover_genomes_in_folder(item)
            genomes.extend(discovered)

    return genomes, names


def parse_map_file(map_file: str) -> tuple:
    """
    解析map文件|Parse map file

    文件格式|File format:
        genome_name1\tpath/to/genome1.fa
        genome_name2\tpath/to/genome2.fa

    Args:
        map_file: map文件路径|Map file path

    Returns:
        tuple: (基因组路径列表, 名称列表)|(List of genome paths, List of names)
    """
    genomes = []
    names = []

    with open(map_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                name = parts[0].strip()
                path = parts[1].strip()

                if os.path.exists(path):
                    names.append(name)
                    genomes.append(path)
                else:
                    print(f"警告|Warning: 文件不存在|File not found: {path}")

    if not genomes:
        raise ValueError(f"map文件中没有有效的基因组|No valid genomes in map file: {map_file}")

    return genomes, names


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 处理基因组输入，返回文件和名称|Process genome inputs, return files and names
        genomes, names_from_map = prepare_genomes(args.genomes)

        if len(genomes) < 2:
            print(f"错误|Error: 至少需要2个基因组|At least 2 genomes required, found: {len(genomes)}")
            sys.exit(1)

        # 处理基因组名称|Process genome names
        # 优先级：map文件 > 命令行参数 > 自动提取
        # Priority: map file > command line > auto-extract
        if names_from_map:
            names = names_from_map
        elif args.names:
            names = args.names.split(',')
        else:
            names = None

        # 处理染色体参数，支持逗号分隔|Process chromosome parameters, support comma-separated
        chromosomes = None
        if args.chromosomes:
            chromosomes = []
            for chr_item in args.chromosomes:
                # 支持逗号分隔|Support comma-separated
                if ',' in chr_item:
                    chromosomes.extend([c.strip() for c in chr_item.split(',')])
                else:
                    chromosomes.append(chr_item.strip())

        # 创建配置|Create configuration
        config = PlotSRConfig(
            genomes=genomes,
            output_dir=args.output_dir,
            names=names,
            threads=args.threads,
            minimap2_preset=args.minimap2_preset,
            min_sr_size=args.min_sr_size,
            output_format=args.output_format,
            font_size=args.font_size,
            dpi=args.dpi,
            space_ratio=args.space_ratio,
            vertical=args.vertical,
            itx=args.itx,
            nosyn=args.nosyn,
            noinv=args.noinv,
            notr=args.notr,
            nodup=args.nodup,
            chromosomes=chromosomes,
            skip_existing=args.skip_existing
        )

        config.validate()

        # 创建日志|Create logger
        log_file = os.path.join(config.output_dir, 'plotsr.log')
        logger_manager = PlotSRLogger(log_file=log_file)
        logger = logger_manager.get_logger()

        # 输出配置信息|Output configuration info
        logger.info(f"输入基因组|Input genomes: {len(config.genomes)}")
        for i, (genome, name) in enumerate(zip(config.genomes, config.names)):
            logger.info(f"  [{i+1}] {name}: {genome}")
        logger.info(f"输出目录|Output directory: {config.output_dir}")
        logger.info(f"线程数|Threads: {config.threads}")

        # 运行流程|Run pipeline
        pipeline = PlotSRPipeline(config, logger)
        success = pipeline.run()

        if success:
            logger.info("分析成功完成|Analysis completed successfully")
            sys.exit(0)
        else:
            logger.error("分析失败|Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
