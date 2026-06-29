"""a-liner pipeline 命令行入口|aliner pipeline CLI entry"""

import argparse
import sys

from .config import AlinerConfig
from .pipeline import AlinerPipeline


def parse_args(argv=None):
    """解析命令行参数|Parse command line arguments"""
    p = argparse.ArgumentParser(
        description="a-liner共线性可视化pipeline（FASTA->minimap2->a-liner）|"
                    "a-liner synteny pipeline (FASTA->minimap2->a-liner)")
    p.add_argument('--ref', required=True, help='参考基因组FASTA|Reference genome FASTA')
    p.add_argument('--query', required=True, help='查询基因组FASTA|Query genome FASTA')
    p.add_argument('--ref-seqs', required=True, type=lambda s: s.split(','),
                   help='ref侧序列（逗号分隔，如 chrZ,chrW 或 chrZ:1-30000000）|'
                        'ref-side seqs (comma-separated)')
    p.add_argument('--query-seqs', required=True, type=lambda s: s.split(','),
                   help='query侧序列（逗号分隔，与ref等长）|query-side seqs (comma-separated, equal length to ref)')
    p.add_argument('-o', '--output-dir', default='./aliner_output', help='输出目录|Output directory')
    p.add_argument('--out-prefix', default='synteny', help='输出文件前缀|Output prefix')
    p.add_argument('--preset', default='asm5', choices=['asm5', 'asm10', 'asm20'],
                   help='minimap2预设（近缘asm5/远缘asm10,asm20）|minimap2 preset')
    p.add_argument('--min-identity', type=int, default=70, help='a-liner identity阈值(%%)|identity threshold')
    p.add_argument('--min-alignment-len', type=int, default=1000, help='最小比对长度|min alignment length')
    p.add_argument('--colormap', type=int, default=5, choices=[0, 1, 2, 3, 4, 5], help='配色(0-5)|colormap')
    p.add_argument('--figure-size', type=float, nargs=2, default=[6, 0], metavar=('W', 'H'),
                   help='图尺寸[宽 高]，高0自适应|figure size [w h], h=0 auto')
    p.add_argument('--threads', '-t', type=int, default=12, help='线程数|Threads')
    p.add_argument('--extra-args', default='', help='透传给a-liner的额外参数|pass-through args to a-liner')
    return p.parse_args(argv)


def main():
    """主函数|Main function"""
    args = parse_args()
    try:
        config = AlinerConfig(
            ref_fasta=args.ref, query_fasta=args.query,
            ref_seqs=args.ref_seqs, query_seqs=args.query_seqs,
            output_dir=args.output_dir, out_prefix=args.out_prefix,
            preset=args.preset, min_identity=args.min_identity,
            min_alignment_len=args.min_alignment_len, colormap=args.colormap,
            figure_size=tuple(args.figure_size), threads=args.threads,
            extra_args=args.extra_args,
        )
        AlinerPipeline(config).run()
        sys.exit(0)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
