"""
QIIME2主程序模块|QIIME2 Main Program Module
"""

import argparse
import os
import sys


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='QIIME2微生物组多样性分析|QIIME2 Microbiome Diversity Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i raw_reads/ -o qiime2_output
  %(prog)s -i raw_reads/ --amplicon its --skip-cutadapt
  %(prog)s -i raw_reads/ --method otu --trunc-len-f 220 --trunc-len-r 200
  %(prog)s -i raw_reads/ --classifier ~/db/silva_classifier.qza

截断长度建议|Truncation length guide (V3-V4):
  2x250 bp: --trunc-len-f 220 --trunc-len-r 200
  2x300 bp: --trunc-len-f 230 --trunc-len-r 210
  (默认0=不截断,依据demux质量图调整|default 0=no truncation, adjust by demux QC)
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input-dir',
                       required=True,
                       help='双端FASTQ输入目录|Input directory of paired-end FASTQ')

    # 输出|Output
    parser.add_argument('-o', '--output-dir',
                       default='./qiime2_output',
                       help='输出目录|Output directory (default: ./qiime2_output)')

    # 分析类型|Analysis type
    parser.add_argument('--amplicon',
                       choices=['16s', 'its'],
                       default='16s',
                       help='扩增子类型|Amplicon type (default: 16s)')
    parser.add_argument('--method',
                       choices=['asv', 'otu'],
                       default='asv',
                       help='聚类方法|Method: ASV(DADA2) or OTU(vsearch) (default: asv)')

    # 引物与截断|Primers and truncation
    parser.add_argument('--fwd-primer',
                       default='CCTACGGGNGGCWGCAG',
                       help='正向引物序列(IUPAC)|Forward primer (default: 341F)')
    parser.add_argument('--rev-primer',
                       default='GACTACHVGGGTATCTAATCC',
                       help='反向引物序列(IUPAC)|Reverse primer (default: 806R)')
    parser.add_argument('--trunc-len-f', type=int, default=0,
                       help='R1截断长度(0=不截断)|R1 truncation length (0=none)')
    parser.add_argument('--trunc-len-r', type=int, default=0,
                       help='R2截断长度(0=不截断)|R2 truncation length (0=none)')
    parser.add_argument('--trim-left-f', type=int, default=0,
                       help='R1左侧裁剪|R1 trim left')
    parser.add_argument('--trim-left-r', type=int, default=0,
                       help='R2左侧裁剪|R2 trim left')

    # 抽样深度|Sampling depth
    parser.add_argument('--sampling-depth', type=int, default=0,
                       help='抽平深度(0=自动取第10百分位)|Rarefaction depth (0=auto)')

    # 聚类/分类参数|Clustering/classification params
    parser.add_argument('--perc-identity', type=float, default=0.97,
                       help='OTU聚类相似度|OTU identity (default: 0.97)')
    parser.add_argument('--confidence', type=float, default=0.7,
                       help='classify-sklearn置信度|classification confidence (default: 0.7)')
    parser.add_argument('--min-length', type=int, default=50,
                       help='extract-reads最小长度|extract-reads min length (default: 50)')
    parser.add_argument('--max-length', type=int, default=0,
                       help='extract-reads最大长度(0=不限)|extract-reads max length (0=none)')

    # 运行控制|Run control
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads (default: 12)')
    parser.add_argument('--validate-level', choices=['min', 'max'], default='min',
                       help='tools import校验级别|import validate level (default: min)')

    # 分类器配置|Classifier configuration
    parser.add_argument('--classifier',
                       help='预训练分类器(.qza),省略则自动训练|Pre-trained classifier (.qza), auto-train if omitted')

    # 路径配置|Path configuration
    parser.add_argument('--database-dir',
                       help='原始参考库目录(SILVA/UNITE)|Raw reference DB directory (default: ~/database/qiime2)')
    parser.add_argument('--qiime-path',
                       help='qiime可执行文件路径|qiime executable path (default: ~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime)')
    parser.add_argument('--classifier-cache-dir',
                       help='分类器缓存目录|Classifier cache directory (default: <db>/classifier_cache)')

    # 样品命名|Sample naming
    parser.add_argument('--r1-suffix', default='_1.clean.fq.gz',
                       help='R1文件后缀|R1 file suffix (default: _1.clean.fq.gz)')
    parser.add_argument('--r2-suffix', default='_2.clean.fq.gz',
                       help='R2文件后缀|R2 file suffix (default: _2.clean.fq.gz)')

    # 跳过控制|Skip control
    parser.add_argument('--skip-cutadapt', action='store_true',
                       help='跳过引物切除(数据已去引物)|Skip primer trimming')
    parser.add_argument('--skip-phylogeny', action='store_true',
                       help='跳过系统发育建树(ITS自动跳过)|Skip phylogeny (auto for ITS)')

    # 元数据|Metadata
    parser.add_argument('--metadata-file',
                       help='样品元数据TSV(可选)|Sample metadata TSV (optional)')

    # 杂项|Misc
    parser.add_argument('-f', '--force', action='store_true',
                       help='覆盖已有输出|Overwrite existing outputs')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        from .config import Qiime2Config
        from .pipeline import Qiime2Pipeline
        from .classifier import ClassifierTrainer
        from .utils import Qiime2Logger, CommandRunner

        # 构建配置|Build configuration
        config_kwargs = {
            'input_dir': args.input_dir,
            'output_dir': args.output_dir,
            'amplicon': args.amplicon,
            'method': args.method,
            'fwd_primer': args.fwd_primer,
            'rev_primer': args.rev_primer,
            'trunc_len_f': args.trunc_len_f,
            'trunc_len_r': args.trunc_len_r,
            'trim_left_f': args.trim_left_f,
            'trim_left_r': args.trim_left_r,
            'sampling_depth': args.sampling_depth,
            'perc_identity': args.perc_identity,
            'confidence': args.confidence,
            'min_length': args.min_length,
            'max_length': args.max_length,
            'threads': args.threads,
            'validate_level': args.validate_level,
            'r1_suffix': args.r1_suffix,
            'r2_suffix': args.r2_suffix,
            'skip_cutadapt': args.skip_cutadapt,
            'skip_phylogeny': args.skip_phylogeny,
            'force': args.force,
            'verbose': args.verbose,
        }

        # 路径配置(仅当用户指定时覆盖默认值)|Path config (override defaults only when specified)
        if args.classifier:
            config_kwargs['classifier'] = args.classifier
        if args.database_dir:
            config_kwargs['database_dir'] = args.database_dir
        if args.qiime_path:
            config_kwargs['qiime_bin'] = args.qiime_path
        if args.classifier_cache_dir:
            config_kwargs['classifier_cache_dir'] = args.classifier_cache_dir
        if args.metadata_file:
            config_kwargs['metadata_file'] = args.metadata_file

        config = Qiime2Config(**config_kwargs)
        config.validate()

        # ITS使用默认16S引物且未跳过cutadapt时提示|Warn if ITS uses default 16S primers without skip
        if config.amplicon == 'its' and not config.skip_cutadapt:
            print("提示|Note: ITS建议使用ITS专用引物或加--skip-cutadapt|"
                  "ITS: use ITS primers or --skip-cutadapt", file=sys.stderr)

        # 创建日志|Create logger
        log_file = os.path.join(config.log_dir, 'qiime2_pipeline.log')
        logger_manager = Qiime2Logger(log_file)
        logger = logger_manager.get_logger()

        # 获取或训练分类器(物种注释前提)|Get or train classifier (prerequisite for taxonomy)
        cmd_runner = CommandRunner(logger, config.output_dir)
        trainer = ClassifierTrainer(config, logger, cmd_runner)
        classifier_qza, classifier_source = trainer.get_or_train()

        # 运行流程|Run pipeline
        pipeline = Qiime2Pipeline(config, logger, classifier_qza, classifier_source)
        success = pipeline.run_pipeline()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"运行错误|Runtime error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"未预期的错误|Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
