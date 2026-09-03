"""全长转录本分析命令行入口|Full-length RNA analysis CLI entry"""

import argparse
import sys

from .config import RnaIsoConfig
from .pipeline import RnaIsoPipeline
from .utils import RnaIsoLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='全长转录本分析(IsoSeq3+IsoQuant)|Full-length transcript analysis (IsoSeq3+IsoQuant)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='示例|Examples: biopytools rna-iso -i sample.fq.gz --data-type ont -g genome.fa -o output/'
    )

    parser.add_argument('-i', '--reads', nargs='+', required=True,
                        help='输入文件(可多个):subreads.bam/ccs.bam/fasta/fastq(±gz)|Input files (multiple allowed)')
    parser.add_argument('--data-type', choices=['pacbio', 'ont'], default=None,
                        help='reads文件(fasta/fastq)时必填;BAM输入自动嗅探|Required for fasta/fastq inputs; auto-sniffed for BAM')
    parser.add_argument('-g', '--reference', default=None,
                        help='参考基因组FASTA(isoquant引擎必填)|Reference genome FASTA (required for isoquant engine)')
    parser.add_argument('--genedb', default=None,
                        help='参考注释GTF/GFF(可选,提高精确度)|Reference annotation GTF/GFF (optional, improves precision)')
    parser.add_argument('--engine', choices=['isoquant', 'isoseq3', 'both'], default='isoquant',
                        help='转录本重建引擎(默认isoquant)|Transcript engine (default: isoquant)')
    parser.add_argument('--primers', default=None,
                        help='引物fasta(默认内置Clontech SMARTer)|Primer fasta (default: built-in Clontech SMARTer)')
    parser.add_argument('--min-passes', type=int, default=1,
                        help='ccs最小pass数(默认1,Iso-Seq官方推荐)|ccs min passes (default: 1)')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数(默认12)|Number of threads (default: 12)')
    parser.add_argument('-p', '--prefix', default='rna_sample',
                        help='样本前缀(默认rna_sample)|Sample prefix')
    parser.add_argument('-o', '--output-dir', default='./rna_iso_output',
                        help='输出目录(默认./rna_iso_output)|Output directory')

    return parser.parse_args(argv)


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = RnaIsoConfig(**vars(args))
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error:\n{e}", file=sys.stderr)
        sys.exit(1)

    logger_manager = RnaIsoLogger(config.log_dir)
    logger = logger_manager.get_logger()

    pipeline = RnaIsoPipeline(config)
    success = pipeline.run(logger)

    if not success:
        logger.error("流程失败|Pipeline failed")
        sys.exit(1)


if __name__ == '__main__':
    main()
