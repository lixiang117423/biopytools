"""
EviAnn基因组注释主程序模块|EviAnn Genome Annotation Main Module
"""

import argparse
import os
import sys

from .calculator import EviAnnRunner
from .config import EviAnnConfig
from .utils import ModuleLogger


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='EviAnn基因组注释工具(证据驱动的真核基因组注释)|'
                    'EviAnn Genome Annotation Tool (evidence-based)')

    # 必需参数|Required arguments
    parser.add_argument('-g', '--genome',
                        required=True,
                        help='基因组FASTA文件|Genome FASTA file')
    parser.add_argument('-o', '--output-dir',
                        required=True,
                        help='输出目录|Output directory')

    # 数据输入(三选一或-e)|Data inputs (one of three, or -e)
    parser.add_argument('--rnaseq-data',
                        help='转录组数据文件或目录(逗号分隔多个),'
                             '自动识别二代/三代|RNA-seq file(s) or dir(s), '
                             'comma-separated; auto-classified')
    parser.add_argument('--sample-sheet',
                        help='样本清单TSV(sample_id/r1/r2/long_reads,'
                             '逗号分隔多文件)|Sample sheet TSV')
    parser.add_argument('-r', '--rnaseq',
                        help='EviAnn原生-r描述文件(透传)|'
                             'EviAnn native -r description file')
    parser.add_argument('-e', '--transcripts',
                        help='近缘物种转录本FASTA|Transcripts FASTA from '
                             'related species')
    parser.add_argument('-p', '--proteins',
                        help='近缘物种蛋白质FASTA|Proteins FASTA from '
                             'related species')

    # 可选参数|Optional arguments
    parser.add_argument('-s', '--uniprot',
                        help='UniProt-SwissProt FASTA|UniProt-SwissProt FASTA')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('-m', '--max-intron',
                        type=int,
                        help='最大内含子长度|Maximum intron length '
                             '(default: auto)')
    parser.add_argument('-d', '--ploidy',
                        type=int,
                        default=2,
                        help='基因组倍性|Genome ploidy (default: 2)')
    parser.add_argument('-c', '--cds-gff',
                        help='含现有CDS的GFF|GFF with existing CDS')
    parser.add_argument('--lncrna-tpm',
                        type=float,
                        default=1.0,
                        help='lncRNA最小TPM|Minimum TPM for lncRNA '
                             '(default: 1.0)')
    parser.add_argument('--min-prot',
                        type=int,
                        help='无同源证据时ab initio ORF最小蛋白长度(aa)|'
                             'Min protein length for ab initio ORF '
                             '(default: 75)')

    # 布尔选项|Boolean options
    parser.add_argument('--partial',
                        action='store_true',
                        help='包含部分CDS|Include partial CDS')
    parser.add_argument('-f', '--functional',
                        action='store_true',
                        help='执行功能注释|Perform functional annotation')
    parser.add_argument('--mito-contigs',
                        help='线粒体contig列表文件|File with mitochondrial '
                             'contigs')
    parser.add_argument('--extra-gff',
                        help='额外GFF特征|Extra features from external GFF')
    parser.add_argument('--debug',
                        action='store_true',
                        help='保留中间文件|Keep intermediate files')
    parser.add_argument('--verbose',
                        action='store_true',
                        help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = EviAnnConfig(**vars(args))
        config.validate()

        # 日志目录先建好(FileHandler需要)|Create logs dir first (FileHandler)
        log_dir = os.path.join(config.output_dir, '99_logs')
        os.makedirs(log_dir, exist_ok=True)
        logger_manager = ModuleLogger(
            log_file=os.path.join(log_dir, 'eviann.log'))
        logger = logger_manager.logger

        logger.info("=" * 80)
        logger.info("EviAnn基因组注释|EviAnn Genome Annotation")
        logger.info(f"  基因组|Genome: {config.genome}")
        logger.info(f"  输出目录|Output dir: {config.output_dir}")
        logger.info("=" * 80)

        runner = EviAnnRunner(config, logger)
        success = runner.run()
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
