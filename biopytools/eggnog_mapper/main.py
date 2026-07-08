"""
eggnog-mapper 功能注释主程序|eggnog-mapper functional annotation main entry.
"""

import argparse
import sys

from .config import EggnogMapperConfig
from .utils import EggnogMapperLogger, EggnogMapperRunner


def parse_arguments():
    """解析命令行参数|Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="eggnog-mapper功能注释工具|eggnog-mapper functional annotation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools eggnog-mapper -i proteins.faa -o out/",
    )

    parser.add_argument("-i", "--input", required=True,
                        help="输入FASTA(蛋白/CDS/基因组)|Input FASTA")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="输出目录|Output directory")
    parser.add_argument("--itype", default="proteins",
                        choices=["proteins", "CDS", "genome", "metagenome"],
                        help="输入类型|Input type (default: proteins)")
    parser.add_argument("--translate", action="store_true",
                        help="CDS翻译为蛋白(itype=CDS/genome/metagenome)|Translate CDS")
    parser.add_argument("-m", "--mode", default="mmseqs",
                        choices=["mmseqs", "diamond", "hmmer", "no_search", "cache"],
                        help="搜索模式|Search mode (default: mmseqs)")
    parser.add_argument("--cpu", type=int, default=12,
                        help="线程数|Threads (default: 12)")
    parser.add_argument("--sensmode", default="sensitive",
                        help="灵敏度|Sensitivity (default: sensitive)")
    parser.add_argument("--seed-ortholog-evalue", type=float, default=0.001,
                        help="seed ortholog E值|evalue (default: 0.001)")
    parser.add_argument("--data-dir", default=None,
                        help="DB目录|DB directory (default: ~/database/eggnog)")
    parser.add_argument("--prefix", default=None,
                        help="输出前缀|Output prefix (default: input stem)")
    parser.add_argument("--emapper-path", default=None,
                        help="emapper.py路径|emapper.py path override")
    parser.add_argument("--resume", action="store_true", help="续传|Resume")
    parser.add_argument("--override", action="store_true",
                        help="覆盖|Override existing output")
    parser.add_argument("--no-format", action="store_true",
                        help="跳过重排版,只留原生产物|Skip reformat")

    return parser.parse_args()


def main():
    """主函数|Main entry."""
    args = parse_arguments()
    try:
        kwargs = dict(
            input_file=args.input,
            output_dir=args.output_dir,
            itype=args.itype,
            translate=args.translate,
            mode=args.mode,
            cpu=args.cpu,
            sensmode=args.sensmode,
            seed_ortholog_evalue=args.seed_ortholog_evalue,
            prefix=args.prefix,
            resume=args.resume,
            override=args.override,
            no_format=args.no_format,
        )
        if args.data_dir:
            kwargs["data_dir"] = args.data_dir
        if args.emapper_path:
            kwargs["emapper_path"] = args.emapper_path

        config = EggnogMapperConfig(**kwargs)
        config.validate()

        logger = EggnogMapperLogger(config.output_path).get_logger()
        ok = EggnogMapperRunner(config, logger).run()
        sys.exit(0 if ok else 1)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
