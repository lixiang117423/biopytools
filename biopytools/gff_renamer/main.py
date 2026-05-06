"""
GFF重命名主程序模块|GFF Renamer Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import GFFRenamerConfig
from .utils import GFFRenamerLogger
from .renamer import GFFRenamer


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="GFF文件ID规范化工具|GFF File ID Standardization Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""示例|Examples: python -m biopytools.gff_renamer -i input.gff -o output.gff -p CDRT -s Ov"""
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='输入GFF文件|Input GFF file')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='输出GFF文件|Output GFF file')

    parser.add_argument('-p', '--prefix',
                        required=True,
                        help='ID前缀|ID prefix (e.g., CDRT, AGIS)')

    parser.add_argument('-s', '--species',
                        required=True,
                        help='物种缩写|Species abbreviation (e.g., Ov, Os)')

    parser.add_argument('-t', '--threads',
                        type=int,
                        default=12,
                        help='并行线程数|Number of parallel threads')

    parser.add_argument('--output-mrna-mapping',
                        action='store_true',
                        help='输出mRNA映射文件|Output mRNA mapping file')

    parser.add_argument('--mrna-mapping-file',
                        help='mRNA映射文件路径（可选，默认为输出文件名_mrna_mapping.tsv）|mRNA mapping file path (optional, defaults to output_filename_mrna_mapping.tsv)')

    parser.add_argument('--chr-mapping',
                        help='染色体映射文件路径|Chromosome mapping file path')

    parser.add_argument('--naming-format',
                        choices=['standard', 'simple', 'compact'],
                        default='standard',
                        help='命名格式|Naming format (standard/simple/compact, default: standard)')

    parser.add_argument('--include-utr',
                        action='store_true',
                        help='包含UTR特征|Include UTR features (five_prime_UTR, three_prime_UTR)')

    parser.add_argument('--skip-clean',
                        action='store_true',
                        help='跳过AGAT清洗步骤|Skip AGAT GFF cleaning step')

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 1.0.0')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = GFFRenamerConfig(
            input_file=args.input,
            output_file=args.output,
            prefix=args.prefix,
            species=args.species,
            threads=args.threads,
            output_mrna_mapping=args.output_mrna_mapping,
            mrna_mapping_file=args.mrna_mapping_file,
            chr_mapping_file=args.chr_mapping,
            naming_format=args.naming_format,
            include_utr=args.include_utr,
            skip_gff_clean=args.skip_clean
        )
        config.validate()

        # 创建日志|Create logger
        logger_manager = GFFRenamerLogger()
        logger = logger_manager.get_logger()

        # 运行重命名|Run renaming
        renamer = GFFRenamer(config, logger)
        success = renamer.rename_gff()

        if success:
            logger.info("重命名成功完成|Renaming completed successfully")
            sys.exit(0)
        else:
            logger.error("重命名失败|Renaming failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
