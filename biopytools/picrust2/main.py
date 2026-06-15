"""
PICRUSt2主程序模块|PICRUSt2 Main Module
"""

import argparse
import os
import sys


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='PICRUSt2微生物群落功能丰度预测|PICRUSt2 Microbial Community Functional Abundance Prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -s study_seqs.fna -i seqabun.biom -o picrust2_out
  %(prog)s -s study_seqs.fna -i seqabun.tsv -o picrust2_out -t 24 --stratified
  %(prog)s -s study_seqs.fna -i seqabun.biom -o picrust2_out --in_traits EC,KO,GO
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-s', '--study-fasta',
                       required=True,
                       help='代表序列FASTA文件|FASTA of unaligned study sequences')
    parser.add_argument('-i', '--input',
                       required=True,
                       help='特征表文件(自动识别BIOM/TSV/Excel/Mothur)|Input table of sequence abundances (auto-detect BIOM/TSV/Excel/Mothur)')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir',
                       default='./picrust2_output',
                       help='输出目录|Output directory (default: ./picrust2_output)')

    # 流程参数|Pipeline parameters
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=12,
                       help='线程数|Number of threads (default: 12)')
    parser.add_argument('--max-nsti',
                       type=float,
                       default=2.0,
                       help='最大NSTI阈值|Maximum NSTI value (default: 2.0)')
    parser.add_argument('--stratified',
                       action='store_true',
                       help='生成分层输出表|Generate stratified output tables')
    parser.add_argument('--in-traits',
                       default='EC,KO',
                       help='功能数据库(EC,KO,GO,PFAM,BIGG,CAZY)|Gene families to predict (default: EC,KO)')
    parser.add_argument('--placement-tool',
                       choices=['epa-ng', 'sepp'],
                       default='epa-ng',
                       help='序列放置工具|Placement tool (default: epa-ng)')
    parser.add_argument('--hsp-method',
                       choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'],
                       default='mp',
                       help='隐状态预测方法|HSP method (default: mp)')
    parser.add_argument('--edge-exponent',
                       type=float,
                       default=0.5,
                       help='HSP edge exponent (default: 0.5)')

    # 流程选择|Pipeline selection
    parser.add_argument('--pipeline',
                       choices=['auto', 'split', 'single'],
                       default='auto',
                       help='流程类型: auto自动检测, split双域, single单参考|Pipeline type (default: auto)')

    # 过滤参数|Filter parameters
    parser.add_argument('--min-align',
                       type=float,
                       default=0.8,
                       help='最小比对比例|Minimum alignment proportion (default: 0.8)')
    parser.add_argument('--min-reads',
                       type=int,
                       default=1,
                       help='每ASV最小reads数|Minimum reads per ASV (default: 1)')
    parser.add_argument('--min-samples',
                       type=int,
                       default=1,
                       help='每ASV最小样本数|Minimum samples per ASV (default: 1)')

    # 通路控制|Pathway control
    parser.add_argument('--no-pathways',
                       action='store_true',
                       help='跳过通路推断|Skip pathway inference')
    parser.add_argument('--coverage',
                       action='store_true',
                       help='计算通路覆盖度|Calculate pathway coverages')
    parser.add_argument('--skip-minpath',
                       action='store_true',
                       help='跳过MinPath|Skip MinPath')
    parser.add_argument('--no-gap-fill',
                       action='store_true',
                       help='跳过gap filling|Skip gap filling')
    parser.add_argument('--per-sequence-contrib',
                       action='store_true',
                       help='逐序列运行MinPath|Run MinPath per sequence')

    # 其他选项|Other options
    parser.add_argument('--skip-norm',
                       action='store_true',
                       help='跳过marker gene拷贝数归一化|Skip normalization by marker gene copies')
    parser.add_argument('--remove-intermediate',
                       action='store_true',
                       help='移除中间文件|Remove intermediate files')
    parser.add_argument('--verbose',
                       action='store_true',
                       help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        from .config import Picrust2Config
        from .pipeline import Picrust2Pipeline
        from .utils import Picrust2Logger

        log_file = os.path.join(args.output_dir, "99_logs", "picrust2_pipeline.log")
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        logger_manager = Picrust2Logger(log_file)
        logger = logger_manager.get_logger()

        config = Picrust2Config(
            study_fasta=args.study_fasta,
            input_table=args.input,
            output_dir=args.output_dir,
            threads=args.threads,
            max_nsti=args.max_nsti,
            stratified=args.stratified,
            in_traits=args.in_traits,
            placement_tool=args.placement_tool,
            hsp_method=args.hsp_method,
            edge_exponent=args.edge_exponent,
            min_align=args.min_align,
            min_reads=args.min_reads,
            min_samples=args.min_samples,
            pipeline=args.pipeline,
            skip_pathways=args.no_pathways,
            coverage=args.coverage,
            skip_norm=args.skip_norm,
            remove_intermediate=args.remove_intermediate,
            verbose=args.verbose,
            per_sequence_contrib=args.per_sequence_contrib,
            skip_minpath=args.skip_minpath,
            no_gap_fill=args.no_gap_fill,
        )

        config.validate()

        pipeline = Picrust2Pipeline(config, logger)
        pipeline.run_pipeline()

        sys.exit(0)

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
