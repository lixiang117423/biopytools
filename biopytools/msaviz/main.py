"""
主程序模块|Main Module
"""

import argparse
import os
import sys
from pathlib import Path
import tempfile

from .msaviz import MsaViz
from .logger import MsaVizLogger
from .aligner import MAFFTAligner


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=' MSA可视化工具|MSA Visualization Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples: %(prog)s -i sequences.fa -o output.png
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--infile',
                       required=True,
                       help=' 输入序列文件或MSA文件|Input sequences or MSA file')
    parser.add_argument('-o', '--outfile',
                       required=True,
                       help=' 输出可视化文件|Output visualization file (*.png|*.jpg|*.svg|*.pdf)')

    # 比对参数|Alignment parameters
    parser.add_argument('--skip-align',
                       action='store_true',
                       help=' 跳过MAFFT比对（输入已是比对结果）|Skip MAFFT alignment (input is already aligned)')
    parser.add_argument('--mafft-path',
                       default='mafft',
                       help=' MAFFT可执行文件路径|MAFFT executable path (default: mafft)')
    parser.add_argument('--mafft-params',
                       default='--auto',
                       help=' MAFFT参数|MAFFT parameters (default: --auto)')
    parser.add_argument('--threads',
                       type=int,
                       default=4,
                       help=' MAFFT线程数|MAFFT threads (default: 4)')
    parser.add_argument('--keep-alignment',
                       action='store_true',
                       default=True,
                       help=' 保留比对结果文件|Keep alignment result file (default: True)')
    parser.add_argument('--no-keep-alignment',
                       action='store_false',
                       dest='keep_alignment',
                       help=' 不保留比对结果文件|Do not keep alignment result file')

    # 格式参数|Format parameters
    parser.add_argument('--format',
                       default='fasta',
                       help=' MSA文件格式|MSA file format (default: fasta)')

    # 颜色方案|Color scheme
    parser.add_argument('--color-scheme',
                       default='Zappo',
                       help=' 颜色方案|Color scheme (default: Zappo)')

    # 区域参数|Region parameters
    parser.add_argument('--start',
                       type=int,
                       default=1,
                       help=' 起始位置(1-based)|Start position (1-based, default: 1)')
    parser.add_argument('--end',
                       type=int,
                       help=' 结束位置(1-based)|End position (1-based, default: MSA length)')

    # 显示参数|Display parameters
    parser.add_argument('--wrap-length',
                       type=int,
                       default=60,
                       help=' 换行长度|Wrap length (default: 60)')
    parser.add_argument('--wrap-space-size',
                       type=float,
                       default=3.0,
                       help=' 换行间距|Wrap space size (default: 3.0)')
    parser.add_argument('--label-type',
                       choices=['id', 'description'],
                       default='id',
                       help=' 标签类型|Label type (default: id)')
    parser.add_argument('--show-label',
                       action='store_true',
                       default=True,
                       help=' 显示序列标签|Show sequence labels (default: True)')
    parser.add_argument('--no-show-label',
                       action='store_false',
                       dest='show_label',
                       help=' 不显示序列标签|Do not show sequence labels')
    parser.add_argument('--show-grid',
                       action='store_true',
                       help=' 显示网格|Show grid')
    parser.add_argument('--show-count',
                       action='store_true',
                       help=' 显示字符统计|Show sequence character count')
    parser.add_argument('--show-consensus',
                       action='store_true',
                       help=' 显示consensus序列|Show consensus sequence')
    parser.add_argument('--consensus-color',
                       default='#1f77b4',
                       help=' Consensus颜色|Consensus color (default: #1f77b4)')
    parser.add_argument('--consensus-size',
                       type=float,
                       default=2.0,
                       help=' Consensus大小|Consensus size (default: 2.0)')

    # 排序|Sorting
    parser.add_argument('--sort',
                       action='store_true',
                       help=' 按NJ树排序|Sort by NJ tree')

    # 输出参数|Output parameters
    parser.add_argument('--dpi',
                       type=int,
                       default=300,
                       help=' 图像DPI|Figure DPI (default: 300)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 检查输入文件是否存在|Check if input file exists
        infile = Path(args.infile)
        if not infile.exists():
            print(f"错误|Error: 输入文件不存在|Input file does not exist: {args.infile}")
            sys.exit(1)

        # 创建输出目录|Create output directory
        outfile = Path(args.outfile)
        outfile.parent.mkdir(parents=True, exist_ok=True)

        # 初始化日志|Initialize logging
        logger_manager = MsaVizLogger()
        logger = logger_manager.get_logger()

        logger.info("=" * 80)
        logger.info(" MSA可视化|MSA Visualization")
        logger.info("=" * 80)

        # 确定输入文件（可能是原始序列或比对结果）|Determine input file
        input_file = infile

        # 默认行为：运行比对（除非用户明确跳过）|Default: run alignment (unless user explicitly skips)
        if not args.skip_align:
            logger.info(" 步骤1/2: MAFFT多序列比对|Step 1/2: MAFFT Multiple Sequence Alignment")

            # 初始化MAFFT比对器|Initialize MAFFT aligner
            aligner = MAFFTAligner(
                mafft_path=args.mafft_path,
                threads=args.threads,
                logger=logger
            )

            # 检查MAFFT是否可用|Check if MAFFT is available
            if not aligner.check_mafft():
                logger.error(f" MAFFT未找到或不可用|MAFFT not found or not available: {args.mafft_path}")
                logger.error(f" 请安装MAFFT或使用--mafft-path指定路径|Please install MAFFT or specify path with --mafft-path")
                logger.error(f" 或者使用--skip-align选项如果输入已是比对结果|Or use --skip-align if input is already aligned")
                sys.exit(1)

            # 检测序列类型|Detect sequence type
            seq_type = aligner.detect_sequence_type(input_file)
            if seq_type:
                logger.info(f" 检测到序列类型|Detected sequence type: {seq_type}")

            # 确定比对输出文件|Determine alignment output file
            if args.keep_alignment:
                # 保留比对文件到输出目录|Keep alignment file in output directory
                alignment_file = outfile.parent / f"{outfile.stem}.aligned.fa"
            else:
                # 使用临时文件|Use temp file
                alignment_file = Path(tempfile.gettempdir()) / f"{infile.stem}.aligned.fa"

            # 运行MAFFT比对|Run MAFFT alignment
            if not aligner.run_alignment(
                input_file,
                alignment_file,
                mafft_params=args.mafft_params
            ):
                logger.error(f" 比对失败|Alignment failed")
                sys.exit(1)

            # 使用比对结果进行可视化|Use alignment result for visualization
            input_file = alignment_file
            if args.keep_alignment:
                logger.info(f" 比对结果已保存|Alignment saved: {alignment_file}")
            else:
                logger.info(f" 比对完成|Alignment completed")

            logger.info(" 步骤2/2: 生成可视化|Step 2/2: Generating visualization")
        else:
            logger.info(f" 跳过比对，直接可视化|Skip alignment, visualize directly")
            logger.info(f" 输入文件|Input file: {args.infile}")

        logger.info(f" 输出文件|Output file: {args.outfile}")

        # 创建MsaViz对象|Create MsaViz object
        mv = MsaViz(
            msa=str(input_file),
            format=args.format,
            color_scheme=args.color_scheme,
            start=args.start,
            end=args.end,
            wrap_length=args.wrap_length,
            wrap_space_size=args.wrap_space_size,
            show_label=args.show_label,
            label_type=args.label_type,
            show_grid=args.show_grid,
            show_count=args.show_count,
            show_consensus=args.show_consensus,
            consensus_color=args.consensus_color,
            consensus_size=args.consensus_size,
            sort=args.sort,
        )

        # 保存图像|Save figure
        logger.info(" 生成可视化图像|Generating visualization...")
        mv.savefig(args.outfile, dpi=args.dpi)

        logger.info(f" 完成|Completed: {args.outfile}")
        logger.info("=" * 80)
        sys.exit(0)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
