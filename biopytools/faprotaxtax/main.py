"""
FAPROTAX主程序模块|FAPROTAX Main Program Module
"""

import argparse
import os
import sys


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='FAPROTAX微生物群落功能注释|FAPROTAX Microbial Community Functional Annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i otu_table.biom -o faprotaxtax_output
  %(prog)s -i otu_table.tsv -o faprotaxtax_output --collapse-by-metadata taxonomy
  %(prog)s -i otu_table.biom -o faprotaxtax_output --group-leftovers-as other -n columns_before_collapsing
        """
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input-table',
                       required=True,
                       help='输入OTU/ASV表（BIOM或TSV格式）|Input OTU/ASV table (BIOM or TSV format)')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir',
                       default='./faprotaxtax_output',
                       help='输出目录|Output directory (default: ./faprotaxtax_output)')

    # 工具路径|Tool paths
    parser.add_argument('-g', '--groups-file',
                       help='FAPROTAX功能组数据库文件路径|Path to FAPROTAX functional groups database file')

    parser.add_argument('--collapse-table-path',
                       help='collapse_table.py脚本路径|Path to collapse_table.py script')

    # 核心参数|Core parameters
    parser.add_argument('--collapse-by-metadata',
                       help='用于功能注释的BIOM元数据字段名（如: taxonomy）|'
                            'BIOM metadata field for functional annotation (e.g., taxonomy)')

    parser.add_argument('--group-leftovers-as',
                       help='未匹配到功能组的记录归为此组名|'
                            'Group name for records not matching any functional group')

    parser.add_argument('-n', '--normalize',
                       choices=[
                           'none', 'columns_before_collapsing', 'rows_before_collapsing',
                           'columns_after_collapsing', 'rows_after_collapsing',
                           'columns_before_collapsing_excluding_unassigned',
                           'rows_before_collapsing_excluding_unassigned'
                       ],
                       default='none',
                       help='标准化方式|Normalization method (default: none)')

    parser.add_argument('--average',
                       choices=['none', 'across_records', 'across_group_members',
                                'across_used_group_members', 'maximum', 'minimum',
                                'minimum_across_records'],
                       default='none',
                       help='组内聚合方式|Aggregation method within groups (default: none)')

    parser.add_argument('--row-names-are-in-column',
                       help='包含行名的列名或索引（仅经典TSV表格）|'
                            'Column name or index containing row names (for TSV tables only)')

    parser.add_argument('--output-format',
                       choices=['auto', 'BIOM', 'classical'],
                       default='auto',
                       help='输出格式|Output format (default: auto)')

    # 控制参数|Control parameters
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=1,
                       help='线程数|Number of threads (default: 1)')

    parser.add_argument('-f', '--force',
                       action='store_true',
                       help='覆盖已存在的输出文件|Overwrite existing output files')

    parser.add_argument('-v', '--verbose',
                       action='store_true',
                       help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        from .config import FaprotaxtaxConfig
        from .pipeline import FaprotaxtaxPipeline
        from .utils import FaprotaxtaxLogger

        # 构建配置|Build configuration
        config_kwargs = {
            'input_table': args.input_table,
            'output_dir': args.output_dir,
            'collapse_by_metadata': args.collapse_by_metadata,
            'group_leftovers_as': args.group_leftovers_as,
            'normalize': args.normalize,
            'average': args.average,
            'row_names_are_in_column': args.row_names_are_in_column,
            'output_format': args.output_format,
            'threads': args.threads,
            'force': args.force,
            'verbose': args.verbose,
        }

        # 工具路径（仅当用户指定时覆盖默认值）|Tool paths (override defaults only when specified)
        if args.groups_file:
            config_kwargs['groups_file'] = args.groups_file
        if args.collapse_table_path:
            config_kwargs['collapse_table_path'] = args.collapse_table_path

        config = FaprotaxtaxConfig(**config_kwargs)
        config.validate()

        # 创建日志|Create logger
        log_file = os.path.join(config.log_dir, "faprotaxtax_pipeline.log")
        logger_manager = FaprotaxtaxLogger(log_file)
        logger = logger_manager.get_logger()

        # 运行流程|Run pipeline
        pipeline = FaprotaxtaxPipeline(config, logger)
        success = pipeline.run_pipeline()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"未预期的错误|Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
