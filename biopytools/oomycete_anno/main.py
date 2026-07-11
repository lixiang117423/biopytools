"""疫霉菌基因组注释命令行入口|Oomycete annotation CLI entry"""

import argparse
import sys

from .config import OomyceteAnnoConfig
from .pipeline import OomyceteAnnoRunner
from .utils import OomyceteAnnoLogger


def parse_arguments():
    """解析命令行参数|Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="疫霉菌基因组注释(T2T Augustus流程)|Oomycete genome annotation (T2T Augustus pipeline)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools oomycete-anno -g genome.fa -s phytophthora --rnaseq-dirs rna1/ rna2/ -o out/ -t 24",
    )

    # 必需|required
    parser.add_argument("-g", "--genome", required=True, help="基因组 FASTA|Genome FASTA")
    parser.add_argument("-s", "--species", required=True,
                        help="Augustus 物种名(简单字母, 如 phytophthora)|Augustus species name (simple, e.g. phytophthora)")

    # 可选证据|optional evidence
    parser.add_argument("--rnaseq-dirs", nargs="+", default=None,
                        help="二代 RNA-seq 目录(可多个)|Short RNA-seq dir(s)")
    parser.add_argument("--prot-seq", default=None,
                        help="同源蛋白文件(Phase2)|Homologous proteins (P2)")
    parser.add_argument("--isoseq", default=None,
                        help="三代转录本文件(Phase2)|Long-read transcripts (P2)")

    # 文件模式/链特异性|file patterns / strandness
    parser.add_argument("--read1-pattern", default="_1.clean.fq.gz",
                        help="R1 文件后缀模式|R1 suffix pattern")
    parser.add_argument("--read2-pattern", default="_2.clean.fq.gz",
                        help="R2 文件后缀模式|R2 suffix pattern")
    parser.add_argument("--rna-strandness", default="",
                        help="链特异性: ''(非链特异性,默认) / FR / RF|Strandness: ''(unstranded,default)/FR/RF")

    # 输出/流程|output & pipeline
    parser.add_argument("-o", "--output-dir", default="./oomycete_anno_output",
                        help="输出目录|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=12, help="线程数|Threads")
    parser.add_argument("--no-soft-masking", dest="soft_masking", action="store_false",
                        help="禁用软屏蔽(改用硬屏蔽)|Disable soft masking (use hard masking)")

    # GeneMark 特殊覆盖(非 conda 路径)|GeneMark special overrides (non-conda)
    parser.add_argument("--gmes-petap-path", default=None,
                        help="GeneMark gmes_petap.pl 路径(默认 ~/software/GeneMark/...)|GeneMark gmes_petap.pl path")
    parser.add_argument("--genemark-perl-env", default=None,
                        help="GeneMark perl 提供环境(默认 braker_v.3.0.8)|GeneMark perl provider env")

    # 步骤跳过|skip flags
    parser.add_argument("--skip-repeat", action="store_true", help="跳过重复屏蔽|Skip repeat masking")
    parser.add_argument("--skip-rna", action="store_true", help="跳过 RNA-seq 比对|Skip RNA-seq alignment")
    parser.add_argument("--skip-iso", action="store_true", help="跳过三代转录本(Phase2)|Skip long-read (P2)")
    parser.add_argument("--skip-protein", action="store_true", help="跳过蛋白证据(Phase2)|Skip protein hints (P2)")
    parser.add_argument("--skip-ltr", action="store_true", help="跳过 LTR 注解(Phase2)|Skip LTR annotation (P2)")

    return parser.parse_args()


def main():
    """主函数|Main function."""
    args = parse_arguments()

    try:
        # 构建 config, 仅在用户显式覆盖时传入特殊路径|build config
        kwargs = dict(
            genome=args.genome,
            species=args.species,
            rnaseq_dirs=args.rnaseq_dirs,
            prot_seq=args.prot_seq,
            isoseq=args.isoseq,
            read1_pattern=args.read1_pattern,
            read2_pattern=args.read2_pattern,
            rna_strandness=args.rna_strandness,
            output_dir=args.output_dir,
            threads=args.threads,
            soft_masking=args.soft_masking,
            skip_repeat=args.skip_repeat,
            skip_rna=args.skip_rna,
            skip_iso=args.skip_iso,
            skip_protein=args.skip_protein,
            skip_ltr=args.skip_ltr,
        )
        if args.gmes_petap_path:
            kwargs["gmes_petap_path"] = args.gmes_petap_path
        if args.genemark_perl_env:
            kwargs["genemark_perl_env"] = args.genemark_perl_env

        config = OomyceteAnnoConfig(**kwargs)
        config.validate()

        logger = OomyceteAnnoLogger(config.log_dir).get_logger()
        ok = OomyceteAnnoRunner(config, logger).run()
        sys.exit(0 if ok else 1)

    except ValueError as e:
        print(f"参数错误|Config error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
