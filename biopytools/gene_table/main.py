"""gene-table 命令行入口|gene-table command-line entry"""

import argparse
import sys

from ..common.paths import expand_path
from .config import GeneTableConfig
from .extractor import GeneTableExtractor
from .utils import GeneTableLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse CLI args"""
    parser = argparse.ArgumentParser(
        prog='gene-table',
        description='基因信息+序列合并表|Gene info + sequence merged table (gene DNA + CDS + Protein)')
    parser.add_argument('-g', '--genome', required=True,
                        help='基因组 FASTA|Genome FASTA')
    parser.add_argument('-f', '--gff', required=True,
                        help='GFF3 注释(支持 .gz)|GFF3 annotation (gz supported)')
    parser.add_argument('-o', '--output', required=True,
                        help='输出表路径(或目录)|Output table path (or directory)')
    parser.add_argument('--prefix', default=None,
                        help='输出文件前缀 + Sample 列(默认取 GFF 文件名)|Output prefix + Sample column')
    parser.add_argument('--longest-only', action='store_true',
                        help='每基因仅保留最长转录本(默认全部)|Keep only longest transcript per gene')
    parser.add_argument('--transcript-types', nargs='+', default=['mRNA', 'transcript'],
                        help='视为转录本的 feature 类型|Feature types treated as transcripts')
    parser.add_argument('--gene-type', default='gene', help='基因 feature 类型|Gene feature type')
    parser.add_argument('--min-length', type=int, default=0,
                        help='基因 DNA 最小长度过滤(0=不过滤)|Min gene-DNA length filter')
    parser.add_argument('--gffread', default=None,
                        help='gffread 路径(默认自动检测)|gffread path (auto-detected by default)')
    parser.add_argument('--log-file', default=None, help='日志文件|Log file')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='日志级别|Log level')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细日志|Verbose')
    return parser.parse_args(argv)


def main():
    """主函数|Main entry"""
    args = parse_arguments()
    try:
        config = GeneTableConfig(
            genome_file=args.genome, gff_file=args.gff, output=args.output,
            prefix=args.prefix, longest_only=args.longest_only,
            transcript_types=tuple(args.transcript_types), gene_type=args.gene_type,
            min_length=args.min_length, log_file=args.log_file, log_level=args.log_level,
            verbose=args.verbose)
        if args.gffread:
            # 用户显式指定优先,且需要展开 ~|explicit user override, expand ~
            config.gffread_path = expand_path(args.gffread)
        config.validate()

        logger = GeneTableLogger(config.log_file, config.log_level, config.verbose).get_logger()
        logger.info(f"输出文件|Outputs: TSV={config.tsv_path} gene.fa={config.gene_fa} "
                    f"cds.fa={config.cds_fa} pep.fa={config.pep_fa}")
        GeneTableExtractor(config, logger).run()
        sys.exit(0)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
