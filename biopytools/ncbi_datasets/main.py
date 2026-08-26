"""NCBI datasets 下载主程序模块|NCBI datasets download main module"""

import argparse
import sys

from .config import INCLUDE_OPTIONS, NCBIDatasetsConfig
from .downloader import NCBIDatasetsDownloader
from .utils import ModuleLogger


def create_argument_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器|Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='NCBI taxon 基因组批量下载(datasets CLI)|NCBI taxon genome batch '
                    'download (datasets CLI)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -t 67593 -o ~/tmp/ncbi_download/
  %(prog)s -t 67593 --assembly-source refseq --assembly-level complete,chromosome \\
        --include-gff3 --include-protein -o ~/tmp/ncbi_download/
  %(prog)s -t 67593 --dry-run -o ~/tmp/ncbi_download/   # 只出清单不下载|manifest only
        '''
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '--taxon', '-t',
        required=True,
        help='NCBI taxon 编号|NCBI taxon ID'
    )

    # 输出设置|Output settings
    parser.add_argument(
        '--output-dir', '-o',
        default='./output',
        help='输出目录|Output directory'
    )

    # 筛选参数|Filter parameters
    parser.add_argument(
        '--assembly-source',
        choices=['refseq', 'genbank'],
        help='只下载 RefSeq 或 GenBank 来源|Only RefSeq or GenBank assemblies'
    )
    parser.add_argument(
        '--assembly-level',
        help='组装级别过滤(逗号分隔)|Assembly level filter (comma-separated): '
             'complete, chromosome, scaffold, contig'
    )
    parser.add_argument(
        '--reference',
        action='store_true',
        help='只下载参考基因组|Reference genomes only'
    )
    parser.add_argument(
        '--annotated',
        action='store_true',
        help='只下载有注释的基因组|Annotated genomes only'
    )

    # 下载内容|Download content
    parser.add_argument(
        '--include-gff3',
        action='store_true',
        help='额外下载 GFF3 基因注释|Also download GFF3 gene annotation'
    )
    parser.add_argument(
        '--include-protein',
        action='store_true',
        help='额外下载蛋白序列|Also download protein sequences'
    )
    parser.add_argument(
        '--include-cds',
        action='store_true',
        help='额外下载 CDS 序列|Also download CDS sequences'
    )
    parser.add_argument(
        '--include-seq-report',
        action='store_true',
        help='额外下载 seq-report 汇总|Also download seq-report summary'
    )

    # 运行模式|Run mode
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='只查询 assembly 清单,不下载|Query manifest only, no download'
    )
    parser.add_argument(
        '--no-organize',
        action='store_true',
        help='下载后不整理到 02_organized|Do not organize into 02_organized after download'
    )

    # 高级选项|Advanced options
    parser.add_argument(
        '--datasets-path',
        help='datasets 工具路径(默认走环境变量/配置/~bin)|datasets tool path '
             '(default: env var / config / ~/bin)'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='日志级别|Log level'
    )

    return parser


def main():
    """主函数|Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()

    try:
        config = NCBIDatasetsConfig(
            taxon=args.taxon,
            output_dir=args.output_dir,
            assembly_source=args.assembly_source,
            assembly_level=args.assembly_level,
            reference=args.reference,
            annotated=args.annotated,
            include_gff3=args.include_gff3,
            include_protein=args.include_protein,
            include_cds=args.include_cds,
            include_seq_report=args.include_seq_report,
            dry_run=args.dry_run,
            organize=not args.no_organize,
            datasets_path=args.datasets_path,
            log_level=args.log_level,
        )
        config.validate()

        # 初始化日志|Initialize logging
        logger_manager = ModuleLogger(str(config.log_file), config.log_level)
        logger = logger_manager.get_logger()

        # 运行下载流程|Run download pipeline
        downloader = NCBIDatasetsDownloader(config, logger)
        success = downloader.run()

        sys.exit(0 if success else 1)
    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
