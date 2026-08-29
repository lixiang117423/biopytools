"""pathorepeat命令行入口|pathorepeat CLI Entry

示例|Examples: biopytools pathorepeat -i genome.fa -o out_dir/
"""

import argparse
import os
import sys

from .config import TESORTER_DBS, PathorepeatConfig
from .pipeline import PathorepeatPipeline
from .utils import PathorepeatLogger, collect_genomes


def parse_arguments(argv=None):
    """解析命令行参数|Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        prog='biopytools pathorepeat',
        description='病原菌重复序列注释:RepeatModeler2从头建模+RepeatMasker软屏蔽'
                    '+TEsorter分类精修,EDTA的非植物替代|Pathogen repeat annotation: '
                    'RepeatModeler2 de novo + RepeatMasker soft mask + TEsorter, '
                    'a non-EDTA alternative',
        epilog='示例|Examples: biopytools pathorepeat -i genome.fa -o out_dir/')
    parser.add_argument('-i', '--input', required=True,
                        help='基因组FASTA或文件夹(批量)|Genome FASTA or directory '
                             '(batch)')
    parser.add_argument('-o', '--output-dir', default='./pathorepeat_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Thread count')
    parser.add_argument('--masking-mode', default='xsmall',
                        choices=['xsmall', 'soft', 'hard', 'x'],
                        help='屏蔽模式(xsmall=小写软屏蔽,保留序列信息,病原菌默认)'
                             '|Masking mode (xsmall=lowercase soft mask, default '
                             'for pathogens)')
    parser.add_argument('--ltr-struct', dest='ltr_struct', action='store_true',
                        default=True, help=argparse.SUPPRESS)
    parser.add_argument('--no-ltr-struct', dest='ltr_struct', action='store_false',
                        help='关闭RepeatModeler -LTRStruct(更快但LTR建库变差)'
                             '|Disable -LTRStruct (faster but worse LTR library)')
    parser.add_argument('--tesorter-db', default='rexdb', choices=TESORTER_DBS,
                        help='TEsorter数据库(REXdb植物/动物lineage为主,卵菌/原生生物'
                             '可试gydb)|TEsorter db (REXdb is plant/metazoa-heavy; '
                             'try gydb for oomycetes/protists)')
    parser.add_argument('--db-hmm', default=None,
                        help='自定义TEsorter HMM文件(优先于--tesorter-db)|Custom '
                             'TEsorter HMM file (overrides --tesorter-db)')
    parser.add_argument('--famdb-dir', default=None,
                        help='Dfam famdb数据目录(含famdb.py与*.h5;设置后注入FAMDB_DIR,'
                             '启用RM2自带分类;不设则分类失败时自动降级)|Dfam famdb dir '
                             '(famdb.py + *.h5; injected as FAMDB_DIR to enable RM2 '
                             'classification; auto-degrades if unset)')
    parser.add_argument('--effector-bed', default=None,
                        help='effector候选区BED(仅单文件模式)|Effector regions BED '
                             '(single-sample mode only)')
    parser.add_argument('--effector-gff', default=None,
                        help='effector候选区GFF3(仅单文件模式)|Effector regions '
                             'GFF3 (single-sample mode only)')
    parser.add_argument('--genome-name', default=None,
                        help='输出文件前缀(仅单文件模式,默认basename剥后缀)|Output '
                             'prefix (single-sample only)')
    parser.add_argument('--skip-completed', dest='skip_completed',
                        action='store_true', default=True, help=argparse.SUPPRESS)
    parser.add_argument('--no-skip-completed', dest='skip_completed',
                        action='store_false',
                        help='全部重跑(忽略已完成步骤)|Rerun all (ignore done steps)')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='日志级别|Log level')
    return parser.parse_args(argv)


def main():
    """主函数|Main function"""
    args = parse_arguments()
    try:
        config = PathorepeatConfig(
            input=args.input, output_dir=args.output_dir,
            threads=args.threads, masking_mode=args.masking_mode,
            ltr_struct=args.ltr_struct, tesorter_db=args.tesorter_db,
            db_hmm=args.db_hmm, famdb_dir=args.famdb_dir,
            effector_bed=args.effector_bed,
            effector_gff=args.effector_gff, genome_name=args.genome_name,
            skip_completed=args.skip_completed, log_level=args.log_level)
        config.validate()

        genomes = collect_genomes(config.input)
        # 单文件 --genome-name 前缀覆盖由 config.sample_name() 生效
        # |Single-sample --genome-name takes effect via config.sample_name()

        logs_dir = os.path.join(config.output_dir, '99_logs')
        logger = PathorepeatLogger(logs_dir, config.log_level).get_logger()
        logger.info(f"pathorepeat模块启动|pathorepeat started "
                    f"(samples={len(genomes)}, batch={config.is_batch}, "
                    f"threads={config.threads}, masking={config.masking_mode})")

        pipeline = PathorepeatPipeline(config, logger)
        ok = pipeline.run(genomes)
        sys.exit(0 if ok else 1)
    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
