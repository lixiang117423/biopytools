"""
K-mer工具主程序模块|K-mer Tools Main Module
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from .config import (
    BuildConfig, QueryConfig, SplitFastaConfig,
    GenFofConfig, ImportDBConfig, ExtractConfig
)
from .utils import (
    KmerToolsLogger, KmtricksChecker, BgzipChecker, RocksDBChecker,
    generate_fof_file, split_fasta_file, format_number,
    ensure_pipeline_dirs, generate_software_versions_yml
)
from .extract import KmerExtractor
from .build import KmerBuildPipeline
from .query import KmerQueryPipeline
from .rocksdb.importer import RocksDBImporter


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='K-mer工具集 - 建库、查询和分析|K-mer Tools - Build, Query and Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
子命令|Subcommands:
  build       构建k-mer库|Build k-mer database
  compare     比较两个kmer矩阵|Compare two kmer matrix files
  count       K-mer丰度分析|K-mer abundance analysis
  extract     从FASTA提取k-mer|Extract k-mers from FASTA
  gen-fof     生成FOF文件|Generate FOF file
  import-db   导入RocksDB|Import to RocksDB
  intersect   提取目标kmer丰度|Extract target kmer abundance
  kmer2vcf    Kmer丰度转VCF|Kmer abundance to VCF converter
  query       查询k-mer库|Query k-mer database
  split-fasta 分割FASTA文件|Split FASTA file

示例|Examples:
  # 构建k-mer库|Build k-mer database
  %(prog)s build -i ./fastq_dir -o ./kmer_db

  # 比较两个kmer矩阵|Compare two kmer matrices
  %(prog)s compare -f1 matrix1.txt -f2 matrix2.txt -o comparison

  # K-mer丰度分析|K-mer abundance analysis
  %(prog)s count -i ./fastq_dir -p "*.fq.gz" -k kmer_lib.fa -o ./count_results

  # 从FASTA提取k-mer|Extract k-mers from FASTA
  %(prog)s extract -i input.fa -o output_dir

  # 生成FOF文件|Generate FOF file
  %(prog)s gen-fof --dir ./fastq_dir -o fof.txt

  # 导入RocksDB|Import to RocksDB
  %(prog)s import-db -i matrix.txt.gz -o rocksdb_dir

  # 提取目标kmer丰度|Extract target kmer abundance
  %(prog)s intersect -m kmer_matrix.txt -k target_kmers.fa -o results.txt

  # Kmer丰度转VCF|Kmer abundance to VCF
  %(prog)s kmer2vcf -i kmer_matrix.txt -o output.vcf

  # 查询k-mer库|Query k-mer database
  %(prog)s query -d ./kmer_db/rocksdb -q query.fa -o ./results

  # 分割FASTA文件|Split FASTA file
  %(prog)s split-fasta input.fa output_dir
        '''
    )

    subparsers = parser.add_subparsers(dest='command', help='子命令|Subcommand')

    # ========== build命令|build command ==========
    build_parser = subparsers.add_parser('build', help='构建k-mer库|Build k-mer database')
    build_parser.add_argument('-i', '--input', required=True,
                             help='输入目录(包含FASTQ/FASTA文件)|Input directory (containing FASTQ/FASTA files)')
    build_parser.add_argument('-o', '--output', required=True,
                             help='输出目录|Output directory')

    # 模式选择|Mode selection
    build_parser.add_argument('--use-kmindex', action='store_true',
                             help='使用kmindex模式 (默认: 使用kmtricks)|Use kmindex mode (default: use kmtricks)')

    # k-mer参数|K-mer parameters
    build_parser.add_argument('-k', '--kmer-size', type=int, default=51,
                             help='k-mer大小 (默认: 51)|K-mer size (default: 51)')
    build_parser.add_argument('--hard-min', type=int, default=2,
                             help='最小丰度 (默认: 2)|Minimum abundance (default: 2)')
    build_parser.add_argument('--recurrence-min', type=int, default=1,
                             help='最小重现次数 (默认: 1)|Minimum recurrence (default: 1)')

    # kmtricks参数|kmtricks parameters
    build_parser.add_argument('--mode', default='kmer:pa:bin',
                             help='kmtricks模式 (默认: kmer:pa:bin)|kmtricks mode (default: kmer:pa:bin)')
    build_parser.add_argument('--minimizer-size', type=int, default=10,
                             help='minimizer大小 (默认: 10)|Minimizer size (default: 10)')
    build_parser.add_argument('--tmp-dir', default='',
                             help='kmtricks临时目录 (默认: 输出目录/tmp)|kmtricks temporary directory (default: output_dir/tmp)')
    build_parser.add_argument('--nb-partitions', type=int, default=0,
                             help='分区数 (默认: 0=自动计算, -1=kmtricks默认)|Partitions (default: 0=auto, -1=kmtricks default)')
    build_parser.add_argument('--fof-file', default='',
                             help='预存的FOF文件路径|Pre-existing FOF file path')
    build_parser.add_argument('--header-file', default='',
                             help='预存的header文件路径|Pre-existing header file path')

    # kmindex参数|kmindex parameters
    build_parser.add_argument('--index-name', default='kmer_index',
                             help='kmindex索引名称 (默认: kmer_index)|kmindex index name (default: kmer_index)')
    build_parser.add_argument('--bloom-size', type=int, default=1000000000000,
                             help='kmindex布隆过滤器大小 (默认: 1000000000000)|kmindex bloom filter size (default: 1000000000000)')

    # 通用参数|Common parameters
    build_parser.add_argument('-t', '--threads', type=int, default=64,
                             help='线程数 (默认: 64)|Thread count (default: 64)')
    build_parser.add_argument('--kmtricks-path', default=None,
                             help='kmtricks路径 (默认按 KMTRICKS_PATH环境变量>配置文件>内置默认 解析)|kmtricks path (resolved via KMTRICKS_PATH env>config>built-in)')
    build_parser.add_argument('--kmindex-path', default=None,
                             help='kmindex路径 (默认按 KMINDEX_PATH环境变量>配置文件>内置默认 解析)|kmindex path (resolved via KMINDEX_PATH env>config>built-in)')
    build_parser.add_argument('--bgzip-path', default=None,
                             help='bgzip路径 (默认按 BGZIP_PATH环境变量>配置文件>内置默认 解析)|bgzip path (resolved via BGZIP_PATH env>config>built-in)')
    build_parser.add_argument('--fof-suffix-1', default='_1.clean.fq.gz',
                             help='R1文件后缀 (默认: _1.clean.fq.gz)|R1 file suffix (default: _1.clean.fq.gz)')
    build_parser.add_argument('--fof-suffix-2', default='_2.clean.fq.gz',
                             help='R2文件后缀 (默认: _2.clean.fq.gz)|R2 file suffix (default: _2.clean.fq.gz)')

    # ========== query命令|query command ==========
    query_parser = subparsers.add_parser('query', help='查询k-mer库|Query k-mer database')

    # 模式选择|Mode selection
    query_parser.add_argument('--use-kmindex', action='store_true',
                             help='使用kmindex模式 (默认: 使用kmtricks)|Use kmindex mode (default: use kmtricks)')

    # 通用参数|Common parameters
    query_parser.add_argument('-d', '--database',
                             help='RocksDB数据库目录 (kmtricks模式必需，默认模式)|RocksDB database directory (required for kmtricks mode, default)')
    query_parser.add_argument('-i', '--index',
                             help='kmindex索引目录 (kmindex模式必需，需加--use-kmindex)|kmindex index directory (required for kmindex mode, requires --use-kmindex)')
    query_parser.add_argument('-q', '--query', required=True,
                             help='查询FASTA文件|Query FASTA file')
    query_parser.add_argument('-o', '--output', required=True,
                             help='输出目录|Output directory')
    query_parser.add_argument('-k', '--kmer-size', type=int, default=51,
                             help='k-mer大小，必须与建库时一致 (默认: 51)|K-mer size, must match build (default: 51)')

    # kmtricks参数|kmtricks parameters
    query_parser.add_argument('--header-db-key', default='kmer_header',
                             help='数据库中的header key (默认: kmer_header)|Header key in database (default: kmer_header)')
    query_parser.add_argument('--bed',
                             help='BED文件路径（用于生成位置丰度文件）|BED file path (for generating position-abundance file)')

    # kmindex参数|kmindex parameters
    query_parser.add_argument('--index-name', default='kmer_index',
                             help='kmindex索引名称 (默认: kmer_index)|kmindex index name (default: kmer_index)')
    query_parser.add_argument('--zvalue', type=int, default=0,
                             help='findere算法z值 (默认: 0)|findere z-value (default: 0)')
    query_parser.add_argument('--threshold', type=float, default=0.0,
                             help='共享k-mer阈值 (默认: 0.0)|Shared k-mer threshold (default: 0.0)')
    query_parser.add_argument('--format', default='matrix', choices=['json', 'matrix'],
                             help='输出格式 (默认: matrix)|Output format (default: matrix)')

    # 通用参数|Common parameters
    query_parser.add_argument('-t', '--threads', type=int, default=64,
                             help='线程数 (默认: 64)|Thread count (default: 64)')
    query_parser.add_argument('--kmindex-path', default=None,
                             help='kmindex路径 (默认按 KMINDEX_PATH环境变量>配置文件>内置默认 解析)|kmindex path (resolved via KMINDEX_PATH env>config>built-in)')

    # ========== split-fasta命令|split-fasta command ==========
    split_parser = subparsers.add_parser('split-fasta', help='分割FASTA文件|Split FASTA file')
    split_parser.add_argument('input_fasta', help='输入FASTA文件|Input FASTA file')
    split_parser.add_argument('output_dir', help='输出目录|Output directory')

    # ========== gen-fof命令|gen-fof command ==========
    fof_parser = subparsers.add_parser('gen-fof', help='生成FOF文件|Generate FOF file')
    fof_parser.add_argument('--dir', required=True,
                           help='输入目录|Input directory')
    fof_parser.add_argument('-o', '--output', required=True,
                           help='输出FOF文件|Output FOF file')
    fof_parser.add_argument('--suffix-1', default='_1.clean.fq.gz',
                           help='R1文件后缀 (默认: _1.clean.fq.gz)|R1 file suffix (default: _1.clean.fq.gz)')
    fof_parser.add_argument('--suffix-2', default='_2.clean.fq.gz',
                           help='R2文件后缀 (默认: _2.clean.fq.gz)|R2 file suffix (default: _2.clean.fq.gz)')

    # ========== import-db命令|import-db command ==========
    import_parser = subparsers.add_parser('import-db', help='导入RocksDB|Import to RocksDB')
    import_parser.add_argument('-i', '--input', required=True,
                             help='输入矩阵文件(可以是gzip压缩)|Input matrix file (can be gzipped)')
    import_parser.add_argument('-o', '--output', required=True,
                             help='输出RocksDB目录|Output RocksDB directory')
    import_parser.add_argument('--input-delimiter', default='\t',
                             help='输入分隔符 (默认: tab)|Input delimiter (default: tab)')
    import_parser.add_argument('--batch-size', type=int, default=20000,
                             help='批量写入大小 (默认: 20000)|Batch write size (default: 20000)')
    import_parser.add_argument('--bloom-bits', type=int, default=15,
                             help='Bloom filter位数 (默认: 15)|Bloom filter bits per key (default: 15)')
    import_parser.add_argument('--force-overwrite', action='store_true',
                             help='强制覆盖已存在的数据库|Force overwrite existing database')
    import_parser.add_argument('--header-file',
                             help='Header文件路径|Header file path')
    import_parser.add_argument('--header-db-key', default='kmer_header',
                             help='数据库中的header key (默认: kmer_header)|Header key in database (default: kmer_header)')

    # ========== extract命令|extract command ==========
    extract_parser = subparsers.add_parser('extract', help='从FASTA提取k-mer|Extract k-mers from FASTA')
    extract_parser.add_argument('-i', '--input', required=True,
                               help='输入FASTA文件|Input FASTA file')
    extract_parser.add_argument('-o', '--output', required=True,
                               help='输出目录|Output directory')
    extract_parser.add_argument('-k', '--kmer-size', type=int, default=51,
                               help='k-mer大小 (默认: 51)|K-mer size (default: 51)')
    extract_parser.add_argument('--method', default='unikmer', choices=['unikmer', 'pyfastx'],
                               help='提取方法 (默认: unikmer)|Extraction method (default: unikmer)')
    extract_parser.add_argument('--unikmer-path', default=None,
                               help='unikmer路径 (默认按 UNIKMER_PATH环境变量>配置文件>内置默认 解析)|unikmer path (resolved via UNIKMER_PATH env>config>built-in)')
    extract_parser.add_argument('--kmer-output',
                               help='kmer FASTA文件 (默认: output_dir/basename_kmer_k.fa)|Kmer FASTA file')
    extract_parser.add_argument('--kmer-pos-output',
                               help='kmer位置文件 (默认: output_dir/basename_kmer_k_pos.txt)|Kmer position file')
    extract_parser.add_argument('--no-bed', action='store_true',
                               help='不输出BED文件|Do not output BED file')

    # 通用参数|Common parameters
    extract_parser.add_argument('-t', '--threads', type=int, default=64,
                               help='线程数 (默认: 64)|Thread count (default: 64)')

    # ========== count命令|count command ==========
    count_parser = subparsers.add_parser('count', help='K-mer丰度分析|K-mer abundance analysis')

    # 必需参数|Required parameters
    count_parser.add_argument('-i', '--input', required=True,
                             help='输入文件目录|Input files directory')
    count_parser.add_argument('-p', '--pattern', required=True,
                             help='文件模式，支持FASTQ和FASTA格式|File pattern (FASTQ/FASTA)')
    count_parser.add_argument('-k', '--kmer-lib', required=True,
                             help='K-mer库文件(FASTA格式)|K-mer library file (FASTA format)')
    count_parser.add_argument('-o', '--output', required=True,
                             help='输出目录|Output directory')

    # 可选参数|Optional parameters
    count_parser.add_argument('-b', '--bed-file',
                             help='BED文件路径|BED file path')
    count_parser.add_argument('-m', '--kmer-size', type=int, default=51,
                             help='K-mer长度 (默认: 51)|K-mer size (default: 51)')
    count_parser.add_argument('-s', '--hash-size', default='1000M',
                             help='哈希表大小 (默认: 1000M)|Hash table size (default: 1000M)')
    count_parser.add_argument('-t', '--threads', type=int, default=8,
                             help='线程数 (默认: 8)|Thread count (default: 8)')
    count_parser.add_argument('-w', '--window-size', type=int, default=500000,
                             help='滑动窗口大小bp (默认: 500000)|Window size in bp (default: 500000)')
    count_parser.add_argument('--step-size', type=int,
                             help='滑动窗口步长bp (默认: window-size/5)|Step size in bp (default: window-size/5)')
    count_parser.add_argument('-C', '--canonical', action='store_true',
                             help='统计正向和反向互补链|Count both strands')
    count_parser.add_argument('--keep-temp', action='store_true',
                             help='保留临时文件|Keep temporary files')
    count_parser.add_argument('--keep-binary', action='store_true',
                             help='保留0/1存在缺失矩阵|Keep 0/1 matrix')
    count_parser.add_argument('--jellyfish-path', default=None,
                             help='Jellyfish路径 (默认按 JELLYFISH_PATH环境变量>配置文件>内置默认 解析)|Jellyfish path (resolved via JELLYFISH_PATH env>config>built-in)')
    count_parser.add_argument('-v', '--verbose', action='store_true',
                             help='详细输出|Verbose output')

    # kmer2vcf子命令|kmer2vcf subcommand
    kmer2vcf_parser = subparsers.add_parser('kmer2vcf', help='Kmer丰度转VCF|Kmer abundance to VCF converter')

    # 必需参数|Required parameters
    kmer2vcf_parser.add_argument('-i', '--input-matrix', required=True,
                           help=' 输入kmer丰度矩阵文件|Input kmer abundance matrix file (TSV format)')
    kmer2vcf_parser.add_argument('-o', '--output-vcf', required=True,
                           help=' 输出VCF文件路径|Output VCF file path (.vcf or .vcf.gz)')

    # 模式选择参数|Mode selection arguments
    mode_group = kmer2vcf_parser.add_mutually_exclusive_group()
    mode_group.add_argument('--fast-mode', action='store_true', default=True,
                           help=' 快速模式（单次遍历，默认）|Fast mode (single pass, default)')
    mode_group.add_argument('--standard-mode', action='store_true',
                           help=' 标准模式（3遍处理+排序）|Standard mode (3-pass processing + sorting)')

    # 快速模式参数|Fast mode parameters
    kmer2vcf_parser.add_argument('--chr-length', type=int, default=100000000,
                           help=' 每条染色体长度（快速模式，默认100M）|Chromosome length for fast mode (default: 100M)')
    kmer2vcf_parser.add_argument('--chr-number', type=int, default=0,
                           help=' 染色体数量（如果设置则优先使用）|Number of chromosomes (if set, takes priority over --chr-length)')
    kmer2vcf_parser.add_argument('--min-freq', type=int, default=0,
                           help=' 最小出现频次过滤（快速模式）|Minimum frequency filter for fast mode (default: 0=no filter)')
    kmer2vcf_parser.add_argument('--kmer-length', type=int, default=51,
                           help=' Kmer长度，用于VCF INFO字段|Kmer length for VCF INFO field (default: 51)')
    kmer2vcf_parser.add_argument('--no-header', action='store_true',
                           help=' 输入文件没有header行（第一行就是样本名）|Input file has no header line (first line is sample names)')

    # 标准模式参数|Standard mode parameters
    kmer2vcf_parser.add_argument('-m', '--min-agg-count', type=int, default=3,
                           help=' 最小聚合频次阈值（标准模式）|Minimum aggregated count threshold for standard mode (default: 3)')

    # 通用参数|Common parameters
    kmer2vcf_parser.add_argument('-t', '--threads', type=int, default=12,
                           help=' 线程数（标准模式）|Number of threads for standard mode (default: 12)')
    kmer2vcf_parser.add_argument('-T', '--temp-dir', default=None,
                           help=' 临时文件目录（标准模式）|Temporary directory for standard mode (default: ./temp)')

    # intersect子命令|intersect subcommand
    intersect_parser = subparsers.add_parser('intersect', help='提取目标kmer的丰度信息|Extract abundance for target kmers')

    # 必需参数|Required parameters
    intersect_parser.add_argument('-m', '--kmer-matrix', required=True,
                           help='kmer矩阵文件|Kmer matrix file (TSV format)')
    intersect_parser.add_argument('-k', '--kmer-fasta', required=True,
                           help='目标kmer的fasta文件|Target kmer fasta file')
    intersect_parser.add_argument('-o', '--output-file', required=True,
                           help='输出文件路径|Output file path')

    # 可选参数|Optional parameters
    intersect_parser.add_argument('-t', '--threads', type=int, default=12,
                           help='线程数 (默认: 12)|Number of threads (default: 12)')
    intersect_parser.add_argument('-r', '--use-reverse-complement', action='store_true', default=True,
                           help='使用反向互补查询 (默认: 启用)|Use reverse complement query (default: enabled)')
    intersect_parser.add_argument('--no-reverse-complement', dest='use_reverse_complement',
                           action='store_false',
                           help='不使用反向互补查询|Do not use reverse complement query')
    intersect_parser.add_argument('--keep-not-found', action='store_true', default=True,
                           help='保留未找到的kmer (默认: 保留)|Keep kmers not found (default: enabled)')
    intersect_parser.add_argument('--no-keep-not-found', dest='keep_not_found',
                           action='store_false',
                           help='不保留未找到的kmer|Do not keep kmers not found')
    intersect_parser.add_argument('-f', '--output-format', choices=['tsv', 'csv'], default='tsv',
                           help='输出格式 (默认: tsv)|Output format (default: tsv)')
    intersect_parser.add_argument('-w', '--window-size', type=int, default=100000,
                           help='窗口大小：每N个kmer统计一次 (默认: 100000)|Window size for statistics per N kmers (default: 100000)')

    # compare子命令|compare subcommand
    compare_parser = subparsers.add_parser('compare', help='比较两个kmer矩阵文件|Compare two kmer matrix files')

    # 必需参数|Required parameters
    compare_parser.add_argument('-f1', '--file1', required=True,
                           help='第一个kmer矩阵文件|First kmer matrix file')
    compare_parser.add_argument('-f2', '--file2', required=True,
                           help='第二个kmer矩阵文件|Second kmer matrix file')
    compare_parser.add_argument('-o', '--output-prefix', required=True,
                           help='输出文件前缀|Output file prefix')

    # 可选参数|Optional parameters
    compare_parser.add_argument('-w', '--window-size', type=int, default=100000,
                           help='窗口大小（行数） (默认: 100000)|Window size in lines (default: 100000)')

    return parser.parse_args()


def cmd_build(args):
    """执行build命令|Execute build command"""
    # 创建配置|Create configuration
    config = BuildConfig(
        input_dir=args.input,
        output_dir=args.output,
        use_kmtricks=not args.use_kmindex,
        kmer_size=args.kmer_size,
        hard_min=args.hard_min,
        recurrence_min=args.recurrence_min,
        mode=args.mode,
        minimizer_size=args.minimizer_size,
        index_name=args.index_name,
        bloom_size=args.bloom_size,
        threads=args.threads,
        kmtricks_path=args.kmtricks_path,
        kmindex_path=args.kmindex_path,
        bgzip_path=args.bgzip_path,
        fof_suffix_1=args.fof_suffix_1,
        fof_suffix_2=args.fof_suffix_2,
        tmp_dir=args.tmp_dir,
        nb_partitions=args.nb_partitions,
        fof_file=args.fof_file,
        header_file=args.header_file,
        run_dir=str(Path(args.output) / "kmtricks_run")
    )

    config.validate()

    # 记录开始时间|Record start time
    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    logs_dir = ensure_pipeline_dirs(config.output_path)
    logger_manager = KmerToolsLogger(logs_dir / "build.log")
    logger = logger_manager.get_logger()

    # 检查依赖|Check dependencies
    if config.use_kmtricks:
        checker = KmtricksChecker(logger, config.kmtricks_path)
        if not checker.check_kmtricks():
            sys.exit(1)

        checker = BgzipChecker(logger, config.bgzip_path)
        if not checker.check_bgzip():
            sys.exit(1)
    else:
        from .utils import KmindexChecker
        checker = KmindexChecker(logger, config.kmindex_path)
        if not checker.check_kmindex():
            sys.exit(1)

    # 运行建库流程|Run build pipeline
    pipeline = KmerBuildPipeline(config, logger)
    results = pipeline.run()

    if results['success']:
        # 记录软件版本与参数(§12.5)|Record software versions and params
        generate_software_versions_yml(
            config.output_path, "kmertools build",
            tools={'kmtricks': config.kmtricks_path, 'kmindex': config.kmindex_path, 'bgzip': config.bgzip_path},
            params={
                'use_kmtricks': config.use_kmtricks, 'kmer_size': config.kmer_size,
                'hard_min': config.hard_min, 'recurrence_min': config.recurrence_min,
                'mode': config.mode, 'minimizer_size': config.minimizer_size,
                'threads': config.threads, 'input_dir': str(config.input_path),
                'output_dir': str(config.output_path),
            },
            start_time=start_time,
        )
        logger.info("建库成功|Build completed successfully")
        sys.exit(0)
    else:
        logger.error(f"建库失败|Build failed: {results.get('error', 'Unknown error')}")
        sys.exit(1)


def cmd_query(args):
    """执行query命令|Execute query command"""
    # 创建配置|Create configuration
    config = QueryConfig(
        rocksdb_dir=args.database if args.database else "",
        query_fasta=args.query,
        output_dir=args.output,
        index_dir=args.index if args.index else "",
        use_kmtricks=not args.use_kmindex,
        kmer_size=args.kmer_size,
        header_db_key=args.header_db_key,
        bed_file=args.bed if args.bed else "",
        index_name=args.index_name,
        zvalue=args.zvalue,
        threshold=args.threshold,
        output_format=args.format,
        threads=args.threads,
        kmindex_path=args.kmindex_path
    )

    config.validate()

    # 记录开始时间|Record start time
    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    logs_dir = ensure_pipeline_dirs(config.output_path)
    logger_manager = KmerToolsLogger(logs_dir / "query.log")
    logger = logger_manager.get_logger()

    # 检查依赖|Check dependencies
    if config.use_kmtricks:
        checker = RocksDBChecker(logger)
        if not checker.check_rocksdb():
            sys.exit(1)
    else:
        from .utils import KmindexChecker
        checker = KmindexChecker(logger, config.kmindex_path)
        if not checker.check_kmindex():
            sys.exit(1)

    # 运行查询流程|Run query pipeline
    pipeline = KmerQueryPipeline(config, logger)
    results = pipeline.run()

    if results['success']:
        generate_software_versions_yml(
            config.output_path, "kmertools query",
            tools={'kmindex': config.kmindex_path},
            params={
                'use_kmtricks': config.use_kmtricks, 'kmer_size': config.kmer_size,
                'threshold': config.threshold, 'zvalue': config.zvalue,
                'output_format': config.output_format, 'threads': config.threads,
                'rocksdb_dir': str(getattr(config, 'rocksdb_path', '')),
                'query_fasta': str(getattr(config, 'query_fasta_path', '')),
            },
            start_time=start_time,
        )
        logger.info("查询成功|Query completed successfully")
        sys.exit(0)
    else:
        logger.error(f"查询失败|Query failed: {results.get('error', 'Unknown error')}")
        sys.exit(1)


def cmd_split_fasta(args):
    """执行split-fasta命令|Execute split-fasta command"""
    # 创建配置|Create configuration
    config = SplitFastaConfig(
        input_fasta=args.input_fasta,
        output_dir=args.output_dir
    )

    config.validate()

    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    logs_dir = ensure_pipeline_dirs(config.output_path)
    logger_manager = KmerToolsLogger(logs_dir / "split-fasta.log")
    logger = logger_manager.get_logger()

    # 执行分割|Execute split
    count = split_fasta_file(config.input_fasta_path, config.output_path, logger)

    if count > 0:
        generate_software_versions_yml(
            config.output_path, "kmertools split-fasta",
            tools={},
            params={'input_fasta': str(config.input_fasta_path), 'split_count': count},
            start_time=start_time,
        )
        logger.info(f"成功分割|Successfully split {count} 个序列|sequences")
        sys.exit(0)
    else:
        logger.error("分割失败|Split failed")
        sys.exit(1)


def cmd_gen_fof(args):
    """执行gen-fof命令|Execute gen-fof command"""
    # 创建配置|Create configuration
    config = GenFofConfig(
        input_dir=args.dir,
        output_file=args.output,
        suffix_1=args.suffix_1,
        suffix_2=args.suffix_2
    )

    config.validate()

    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    output_base = Path(config.output_file).parent
    logs_dir = ensure_pipeline_dirs(output_base)
    logger_manager = KmerToolsLogger(logs_dir / "gen-fof.log")
    logger = logger_manager.get_logger()

    # 生成FOF文件|Generate FOF file
    success = generate_fof_file(
        config.input_path,
        Path(config.output_file),
        config.suffix_1,
        config.suffix_2,
        logger
    )

    if success:
        generate_software_versions_yml(
            output_base, "kmertools gen-fof",
            tools={},
            params={'input_dir': str(config.input_path), 'output_file': str(config.output_file),
                    'suffix_1': config.suffix_1, 'suffix_2': config.suffix_2},
            start_time=start_time,
        )
        logger.info("FOF文件生成成功|FOF file generated successfully")
        sys.exit(0)
    else:
        logger.error("FOF文件生成失败|FOF file generation failed")
        sys.exit(1)


def cmd_import_db(args):
    """执行import-db命令|Execute import-db command"""
    # 创建配置|Create configuration
    config = ImportDBConfig(
        input_matrix=args.input,
        output_db=args.output,
        input_delimiter=args.input_delimiter,
        batch_size=args.batch_size,
        bloom_bits=args.bloom_bits,
        force_overwrite=args.force_overwrite,
        header_file=args.header_file,
        header_db_key=args.header_db_key
    )

    config.validate()

    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    output_base = Path(config.output_db).parent
    logs_dir = ensure_pipeline_dirs(output_base)
    logger_manager = KmerToolsLogger(logs_dir / "import-db.log")
    logger = logger_manager.get_logger()

    # 检查依赖|Check dependencies
    checker = RocksDBChecker(logger)
    if not checker.check_rocksdb():
        sys.exit(1)

    # 创建导入器|Create importer
    importer = RocksDBImporter(
        matrix_file=config.input_matrix,
        db_path=config.output_db,
        input_delimiter=config.input_delimiter,
        batch_size=config.batch_size,
        bloom_bits=config.bloom_bits,
        force_overwrite=config.force_overwrite,
        header_file=config.header_file,
        header_db_key=config.header_db_key,
        logger=logger
    )

    # 执行导入|Execute import
    success = importer.import_to_rocksdb()

    if success:
        generate_software_versions_yml(
            output_base, "kmertools import-db",
            tools={},
            params={'input_matrix': str(config.input_matrix_path),
                    'output_db': str(config.output_db_path),
                    'input_delimiter': repr(config.input_delimiter),
                    'batch_size': config.batch_size, 'bloom_bits': config.bloom_bits,
                    'header_db_key': config.header_db_key},
            start_time=start_time,
        )
        logger.info("RocksDB导入成功|RocksDB import completed successfully")
        sys.exit(0)
    else:
        logger.error("RocksDB导入失败|RocksDB import failed")
        sys.exit(1)


def cmd_extract(args):
    """执行extract命令|Execute extract command"""
    import shutil

    # 查找unikmer的完整路径|Find full path of unikmer
    unikmer_path = args.unikmer_path
    if args.method == "unikmer":
        full_path = shutil.which(args.unikmer_path)
        if full_path:
            unikmer_path = full_path

    # 创建配置|Create configuration
    config = ExtractConfig(
        fasta_file=args.input,
        output_dir=args.output,
        kmer_size=args.kmer_size,
        extract_method=args.method,
        unikmer_path=unikmer_path,
        kmer_output=args.kmer_output if args.kmer_output else "",
        kmer_pos_output=args.kmer_pos_output if args.kmer_pos_output else "",
        extract_output_bed=not args.no_bed,
        threads=args.threads
    )

    config.validate()

    start_time = datetime.now()

    # 初始化日志(§12: 日志集中到 99_logs/)|Init logging (§12: logs centralized to 99_logs/)
    logs_dir = ensure_pipeline_dirs(config.output_path)
    logger_manager = KmerToolsLogger(logs_dir / "extract.log")
    logger = logger_manager.get_logger()

    # 检查依赖（如果使用unikmer）|Check dependencies (if using unikmer)
    if config.extract_method == "unikmer":
        try:
            # unikmer使用'version'子命令而不是'--version'标志|unikmer uses 'version' subcommand
            result = subprocess.run([config.unikmer_path, "version"],
                                  capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"未找到unikmer或执行失败|unikmer not found or failed: {config.unikmer_path}")
                sys.exit(1)
        except FileNotFoundError:
            logger.error(f"未找到unikmer|unikmer not found: {config.unikmer_path}")
            sys.exit(1)
        except Exception as e:
            logger.error(f"检查unikmer时出错|Error checking unikmer: {e}")
            sys.exit(1)

    # 创建提取器|Create extractor
    extractor = KmerExtractor(config, logger)

    # 执行提取|Execute extraction
    success = extractor.extract_kmers_from_fasta()

    if success:
        generate_software_versions_yml(
            config.output_path, "kmertools extract",
            tools={'unikmer': config.unikmer_path} if config.extract_method == 'unikmer' else {},
            params={'fasta_file': str(config.fasta_path), 'kmer_size': config.kmer_size,
                    'extract_method': config.extract_method, 'threads': config.threads},
            start_time=start_time,
        )
        logger.info("K-mer提取成功|K-mer extraction completed successfully")
        sys.exit(0)
    else:
        logger.error("K-mer提取失败|K-mer extraction failed")
        sys.exit(1)


def cmd_count(args):
    """执行count命令|Execute count command"""
    try:
        from .kmer_count.config import KmerCountConfig
        from .kmer_count.main import KmerCountAnalyzer
    except ImportError as e:
        print(f"导入kmer_count模块失败|Failed to import kmer_count module: {e}")
        sys.exit(1)

    # 创建配置对象（模仿argparse的Namespace）|Create config object (mimicking argparse Namespace)
    class ArgsNamespace:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

    args_namespace = ArgsNamespace(
        input=args.input,
        pattern=args.pattern,
        kmer_lib=args.kmer_lib,
        output=args.output,
        bed_file=args.bed_file if args.bed_file else "",
        kmer_size=args.kmer_size,
        hash_size=args.hash_size,
        threads=args.threads,
        window_size=args.window_size,
        step_size=args.step_size,
        canonical=args.canonical,
        keep_temp=args.keep_temp,
        keep_binary=args.keep_binary,
        jellyfish_path=args.jellyfish_path,
        verbose=args.verbose
    )

    try:
        # 创建配置|Create configuration
        config = KmerCountConfig.from_args(args_namespace)
        config.validate()

        # 创建分析器并运行|Create analyzer and run
        analyzer = KmerCountAnalyzer(config)
        analyzer.run_analysis()

        # 记录软件版本与参数(§12.5)|Record software versions and params
        generate_software_versions_yml(
            config.output_dir, "kmertools count",
            tools={'jellyfish': config.jellyfish_path},
            params={'kmer_size': config.kmer_size, 'hash_size': config.hash_size,
                    'threads': config.threads, 'canonical': config.canonical,
                    'window_size': config.window_size, 'step_size': config.step_size,
                    'pattern': config.pattern, 'input_dir': str(config.input_dir)},
        )

    except KeyboardInterrupt:
        print("分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def cmd_vcf(args):
    """执行vcf命令|Execute vcf command"""
    try:
        from .kmer2vcf.config import Kmer2VcfConfig
        from .kmer2vcf.main import Kmer2VcfAnalyzer
    except ImportError as e:
        print(f"导入kmer2vcf模块失败|Failed to import kmer2vcf module: {e}")
        sys.exit(1)

    try:
        # 创建配置|Create configuration
        config = Kmer2VcfConfig(
            input_matrix=args.input_matrix,
            output_vcf=args.output_vcf,
            fast_mode=args.fast_mode if hasattr(args, 'fast_mode') else True,
            standard_mode=args.standard_mode if hasattr(args, 'standard_mode') else False,
            chr_length=args.chr_length if hasattr(args, 'chr_length') else 100000000,
            chr_number=args.chr_number if hasattr(args, 'chr_number') else 0,
            min_freq=args.min_freq if hasattr(args, 'min_freq') else 0,
            kmer_length=args.kmer_length if hasattr(args, 'kmer_length') else 51,
            no_header=args.no_header if hasattr(args, 'no_header') else False,
            min_agg_count=args.min_agg_count,
            threads=args.threads,
            temp_dir=args.temp_dir
        )
        config.validate()

        # 创建分析器并运行|Create analyzer and run
        analyzer = Kmer2VcfAnalyzer(
            input_matrix=args.input_matrix,
            output_vcf=args.output_vcf,
            fast_mode=args.fast_mode if hasattr(args, 'fast_mode') else True,
            standard_mode=args.standard_mode if hasattr(args, 'standard_mode') else False,
            chr_length=args.chr_length if hasattr(args, 'chr_length') else 100000000,
            chr_number=args.chr_number if hasattr(args, 'chr_number') else 0,
            min_freq=args.min_freq if hasattr(args, 'min_freq') else 0,
            kmer_length=args.kmer_length if hasattr(args, 'kmer_length') else 51,
            no_header=args.no_header if hasattr(args, 'no_header') else False,
            min_agg_count=args.min_agg_count,
            threads=args.threads,
            temp_dir=args.temp_dir
        )
        analyzer.run_analysis()

    except KeyboardInterrupt:
        print("分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def cmd_intersect(args):
    """执行intersect命令|Execute intersect command"""
    try:
        from .intersect.main import main as intersect_main
    except ImportError as e:
        print(f"导入intersect模块失败|Failed to import intersect module: {e}")
        sys.exit(1)

    # 构建参数列表|Build argument list
    sys.argv = [
        'intersect',
        '-m', args.kmer_matrix,
        '-k', args.kmer_fasta,
        '-o', args.output_file,
        '-t', str(args.threads),
    ]

    # 添加可选参数|Add optional arguments
    if not args.use_reverse_complement:
        sys.argv.append('--no-reverse-complement')
    if not args.keep_not_found:
        sys.argv.append('--no-keep-not-found')
    if args.output_format != 'tsv':
        sys.argv.extend(['-f', args.output_format])
    if args.window_size != 100000:
        sys.argv.extend(['-w', str(args.window_size)])

    try:
        intersect_main()
    except SystemExit as e:
        return e.code
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def cmd_compare(args):
    """执行compare命令|Execute compare command"""
    try:
        from .compare.main import main as compare_main
    except ImportError as e:
        print(f"导入compare模块失败|Failed to import compare module: {e}")
        sys.exit(1)

    # 构建参数列表|Build argument list
    sys.argv = [
        'compare.py',
        '-f1', args.file1,
        '-f2', args.file2,
        '-o', args.output_prefix,
    ]

    # 添加可选参数|Add optional arguments
    if args.window_size != 100000:
        sys.argv.extend(['-w', str(args.window_size)])

    try:
        compare_main()
    except SystemExit as e:
        return e.code
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    """主函数|Main function"""
    args = parse_arguments()

    if not args.command:
        print("请指定子命令|Please specify a subcommand")
        print("使用 -h/--help 查看帮助|Use -h/--help for help")
        sys.exit(1)

    try:
        if args.command == 'build':
            cmd_build(args)
        elif args.command == 'query':
            cmd_query(args)
        elif args.command == 'extract':
            cmd_extract(args)
        elif args.command == 'count':
            cmd_count(args)
        elif args.command == 'split-fasta':
            cmd_split_fasta(args)
        elif args.command == 'gen-fof':
            cmd_gen_fof(args)
        elif args.command == 'import-db':
            cmd_import_db(args)
        elif args.command == 'kmer2vcf':
            cmd_vcf(args)
        elif args.command == 'intersect':
            cmd_intersect(args)
        elif args.command == 'compare':
            cmd_compare(args)
        else:
            print(f"未知命令|Unknown command: {args.command}")
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"程序执行出错|Program execution error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
