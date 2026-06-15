"""
KMC主程序模块|KMC Main Module
"""

import argparse
import sys
import json
from pathlib import Path

from .config import KMCConfig
from .utils import KMCLogger, format_number
from .kmer_counter import KMCCounter
from .matrix_builder import KMCMatrixBuilder
from .kmer_query import KMCQuery
from .matrix_exporter import KMCMatrixExporter


class KMCManager:
    """KMC管理器|KMC Manager"""

    def __init__(self, **kwargs):
        """初始化KMC管理器|Initialize KMC manager"""
        # 初始化配置|Initialize configuration
        self.config = KMCConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = KMCLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化各模块|Initialize modules
        self.counter = KMCCounter(self.config, self.logger_manager)
        self.matrix_builder = KMCMatrixBuilder(self.config, self.logger_manager)
        self.query = KMCQuery(self.config, self.logger_manager)
        self.exporter = KMCMatrixExporter(self.config, self.logger_manager)

    def run_count(self):
        """运行k-mer统计|Run k-mer counting"""
        self.logger.info("开始k-mer统计模式|Starting k-mer counting mode")

        results = self.counter.run()

        # 输出统计结果|Output statistics
        success_count = sum(1 for v in results.values() if v)
        self.logger.info(
            f"统计完成|Counting completed: "
            f"{success_count}/{len(results)} 个样本成功|samples successful"
        )

        # 保存metadata文件（包含正确的kmer_size）|Save metadata file (with correct kmer_size)
        if all(results.values()):
            self._save_count_metadata()

        return all(results.values())

    def _save_count_metadata(self):
        """保存count步骤的metadata|Save metadata from count step"""
        import json
        from pathlib import Path

        metadata_file = self.config.output_path / 'kmer_metadata.json'
        metadata = {
            'n_kmers': 0,  # 未知，将在matrix步骤填充|Unknown, will fill in matrix step
            'kmer_file': None,  # 未知，将在matrix步骤填充|Unknown, will fill in matrix step
            'kmer_size': self.config.kmer_size,  # 保存k-mer大小|Save k-mer size
            'index_type': None,
            'index_path': None
        }

        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        self.logger.info(f"Metadata文件已保存|Metadata file saved: {metadata_file}")
        self.logger.info(f"k-mer大小|K-mer size: {metadata['kmer_size']}")

    def run_matrix(self):
        """运行丰度矩阵构建|Run abundance matrix building"""
        self.logger.info("开始丰度矩阵构建模式|Starting abundance matrix building mode")

        # 如果指定了input_dir，使用input_dir/kmc_databases作为数据库目录
        # If input_dir specified, use input_dir/kmc_databases as database directory
        if self.config.input_dir:
            from pathlib import Path
            custom_db_path = Path(self.config.input_dir) / 'kmc_databases'
            self.logger.info(f"使用自定义输入目录|Using custom input directory: {custom_db_path}")
            # 临时修改kmc_db_path|Temporarily modify kmc_db_path
            original_db_path = self.matrix_builder.config.kmc_db_path
            self.matrix_builder.config.kmc_db_path = custom_db_path

        # 不传递 sample_names，让它自动扫描|Don't pass sample_names, let it auto-scan
        success = self.matrix_builder.build_matrix()

        # 如果修改了kmc_db_path，恢复原值|If modified kmc_db_path, restore original value
        if self.config.input_dir:
            self.matrix_builder.config.kmc_db_path = original_db_path

        if success:
            # 输出矩阵统计信息|Output matrix statistics
            stats = self.matrix_builder.get_matrix_statistics()

            self.logger.info("矩阵统计|Matrix statistics:")
            self.logger.info(f"  k-mer数量|K-mer count: {format_number(stats.get('n_kmers', 0))}")
            self.logger.info(f"  样本数量|Sample count: {stats.get('n_samples', 0)}")
            self.logger.info(f"  k-mer长度|K-mer size: {stats.get('kmer_size', 0)}")
            self.logger.info(f"  存储类型|Storage type: {stats.get('storage_type', 'unknown')}")

            if 'sparsity_percent' in stats:
                self.logger.info(f"  稀疏度|Sparsity: {stats['sparsity_percent']}")

            # 自动导出TSV文件（默认开启，除非指定 --no-export）|Auto export TSV files (default on, unless --no-export specified)
            no_export = getattr(self, 'no_export', False)
            if not no_export:
                self.logger.info("自动导出TSV文件|Auto exporting TSV files")
                abundance_file = str(self.config.output_path / 'abundance_matrix.tsv')
                presence_file = str(self.config.output_path / 'presence_matrix.tsv')
                self._auto_export_tsv(abundance_file, presence_file)

        return success

    def _auto_export_tsv(self, abundance_file: str, presence_file: str):
        """自动导出为TSV格式（简单版，用于kmc matrix）|Auto export to TSV format (simple, for kmc matrix)

        Args:
            abundance_file: 丰度文件路径|Abundance file path
            presence_file: 存在文件路径|Presence file path (0/1)
        """
        from .matrix_exporter import KMCMatrixExporter
        import h5py
        from pathlib import Path

        exporter = KMCMatrixExporter(self.config, self.logger_manager)

        # 导出丰度矩阵|Export abundance matrix
        self.logger.info(f"导出丰度文件|Exporting abundance file: {abundance_file}")
        exporter.export_to_tsv(
            output_file=abundance_file,
            format='full',
            min_abundance=1  # 不过滤，导出所有|No filtering, export all
        )

        # 导出存在/缺失矩阵 (0/1)|Export presence/absence matrix (0/1)
        self.logger.info(f"导出存在文件|Exporting presence file: {presence_file}")
        self._export_presence_matrix(presence_file)

        self.logger.info("TSV文件导出完成|TSV files exported")

    def _export_presence_matrix(self, output_file: str):
        """导出存在/缺失矩阵 (0/1)|Export presence/absence matrix (0/1)

        Args:
            output_file: 输出文件路径|Output file path
        """
        import h5py
        from pathlib import Path

        matrix_file = self.config.output_path / 'abundance_matrix.h5'
        kmer_dict_file = self.config.output_path / 'kmer_dictionary.h5'

        if not matrix_file.exists():
            self.logger.error(f"矩阵文件不存在|Matrix file not exists: {matrix_file}")
            return

        # 创建输出目录|Create output directory
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(matrix_file, 'r') as f:
            # 获取样本名称|Get sample names
            sample_names = [s.decode('utf-8') if isinstance(s, bytes) else s
                           for s in f['sample_names'][:]]

            # 获取存储类型|Get storage type
            storage_type = f.attrs['storage_type']
            if isinstance(storage_type, bytes):
                storage_type = storage_type.decode('utf-8')

            # 获取 k-mer 列表|Get k-mer list
            with h5py.File(kmer_dict_file, 'r') as kf:
                kmer_list = [k.decode('utf-8') if isinstance(k, bytes) else k
                             for k in kf['kmer'][:]]

            n_kmers = len(kmer_list)
            n_samples = len(sample_names)

            self.logger.info(f"导出 {n_kmers} k-mers x {n_samples} samples 到存在矩阵|"
                           f"Exporting {n_kmers} k-mers x {n_samples} samples to presence matrix")

            # 写入TSV|Write to TSV
            with open(output_file, 'w') as out:
                # 写入表头|Write header
                out.write('kmer\t' + '\t'.join(sample_names) + '\n')

                if storage_type == 'sparse':
                    # 稀疏矩阵：构建存在矩阵|Sparse matrix: build presence matrix
                    # 初始化全0矩阵|Initialize all-zero matrix
                    presence = [[0] * n_samples for _ in range(n_kmers)]

                    # 填充存在的k-mer|Fill existing k-mers
                    kmer_ids = f['kmer_id'][:]
                    sample_ids = f['sample_id'][:]

                    for kid, sid in zip(kmer_ids, sample_ids):
                        if kid < n_kmers and sid < n_samples:
                            presence[kid][sid] = 1

                    # 写入数据|Write data
                    for i, kmer in enumerate(kmer_list):
                        row = [kmer] + [str(presence[i][j]) for j in range(n_samples)]
                        out.write('\t'.join(row) + '\n')

                else:
                    # 密集矩阵|Dense matrix
                    abundance = f['abundance'][:]

                    for i, kmer in enumerate(kmer_list):
                        row = [kmer] + ['1' if abundance[i, j] > 0 else '0'
                                        for j in range(n_samples)]
                        out.write('\t'.join(row) + '\n')

        self.logger.info(f"存在文件已导出|Presence file exported: {output_file}")

    def run_query(self, query_fasta: str, output_file: str = None):
        """运行k-mer查询|Run k-mer query

        Args:
            query_fasta: FASTA文件路径|FASTA file path
            output_file: 输出文件路径（可选）|Output file path (optional)
        """
        self.logger.info("开始k-mer查询模式|Starting k-mer query mode")

        # FASTA批量查询（优化版，支持百万级）|FASTA batch query (optimized, supports millions)
        self.logger.info(f"从FASTA文件批量查询|Batch query from FASTA: {query_fasta}")

        # 使用优化的批量查询方法|Use optimized batch query method
        results = self.query.query_kmers_from_fasta(query_fasta)

        # 输出结果|Output results
        if output_file:
            self._save_batch_query_results(results, output_file)
        else:
            self._print_batch_query_results(results)

        return True

    def _save_batch_query_results(self, results: dict, output_file: str):
        """保存批量查询结果到文件|Save batch query results to file"""
        self.logger.info(f"保存结果到文件|Saving results to file: {output_file}")

        with open(output_file, 'w') as f:
            # 获取样本名称（排除 'sequence' 键）|Get sample names (exclude 'sequence' key)
            if results:
                first_result = next(iter(results.values()))
                sample_names = [k for k in first_result.keys() if k != 'sequence']
            else:
                sample_names = []

            # 写入表头：kmer_name + sequence + sample_names|Write header
            f.write("kmer_name\tsequence\t" + "\t".join(sample_names) + "\n")

            # 写入每个k-mer的结果|Write results for each k-mer
            for kmer_name, result in results.items():
                sequence = result.get('sequence', 'N/A')
                row = [kmer_name, sequence] + [str(result.get(s, 0)) for s in sample_names]
                f.write("\t".join(row) + "\n")

        self.logger.info(f"批量查询结果已保存|Batch query results saved: {len(results)} k-mers")

    def _print_batch_query_results(self, results: dict):
        """打印批量查询结果|Print batch query results"""
        self.logger.info(f"批量查询结果|Batch query results ({len(results)} k-mers):")

        for kmer_name, result in results.items():
            sequence = result.get('sequence', 'N/A')

            # 提取丰度数据（排除 'sequence' 键）|Extract abundance data (exclude 'sequence' key)
            abundances = {k: v for k, v in result.items() if k != 'sequence'}

            if abundances:
                sorted_abundances = sorted(abundances.items(), key=lambda x: x[1], reverse=True)
                self.logger.info(f"{kmer_name} ({sequence}):")
                for sample, abundance in sorted_abundances:
                    self.logger.info(f"  {sample}: {abundance}")
            else:
                self.logger.warning(f"{kmer_name} ({sequence}): 未找到|Not found")

    def run_add(self):
        """运行添加新样本模式|Run add new samples mode"""
        self.logger.info("开始添加新样本模式（增量式）|Starting add new samples mode (incremental)")

        # 从已有数据库读取k-mer大小|Read k-mer size from existing database
        import json
        metadata_file = self.config.output_path / 'kmer_metadata.json'

        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
                original_kmer_size = self.config.kmer_size
                self.config.kmer_size = metadata.get('kmer_size', self.config.kmer_size)
                if original_kmer_size != self.config.kmer_size:
                    self.logger.info(f"从metadata读取k-mer大小|Read k-mer size from metadata: "
                                   f"{original_kmer_size} → {self.config.kmer_size}")
                else:
                    self.logger.info(f"使用k-mer大小|Using k-mer size: {self.config.kmer_size}")
        else:
            self.logger.warning(f"元数据文件不存在，使用配置的k-mer大小|Metadata file not exists, using configured k-mer size: {self.config.kmer_size}")

        # 处理目录输入|Handle directory input
        if self.config.input_dir:
            self.logger.info(f"检测到目录输入模式|Detected directory input mode: {self.config.input_dir}")
            self.logger.info(f"数据模式|Data mode: {'单末端|Single-end' if self.config.single_end else '双末端|Paired-end'}")
            self.logger.info(f"文件模式|File pattern: *{self.config.read1_suffix}" +
                            (f", *{self.config.read2_suffix}" if not self.config.single_end else ""))

            # 使用KMCSampleFinder查找样本|Use KMCSampleFinder to find samples
            from .sample_finder import KMCSampleFinder
            sample_finder = KMCSampleFinder(self.config, self.logger)

            samples = sample_finder.find_samples()

            if not samples:
                self.logger.error("未找到有效的样本文件|No valid sample files found")
                return False

            # 验证样本|Validate samples
            if not sample_finder.validate_samples(samples):
                return False

            # 提取样本名和文件路径|Extract sample names and file paths
            new_samples = []
            new_sample_names = []

            for sample_name, file_paths in samples:
                # KMC可以接受多个输入文件，用空格分隔|KMC accepts multiple input files, space-separated
                input_str = ' '.join(str(f) for f in file_paths)
                new_samples.append(input_str)
                new_sample_names.append(sample_name)

            self.logger.info(f"找到 {len(new_samples)} 个新样本|Found {len(new_samples)} new samples: {new_sample_names}")

        else:
            # 使用命令行指定的文件列表|Use file list from command line
            new_samples = self.config.input_files if self.config.input_files else []
            new_sample_names = self.config.sample_names if self.config.sample_names else []

            if new_samples:
                self.logger.info(f"使用文件列表模式|Using file list mode: {len(new_samples)} 个样本|samples")
            else:
                self.logger.info("未指定新样本，仅更新矩阵|No new samples specified, only updating matrix")

        # 使用增量式添加（优化版）|Use incremental add (optimized)
        success = self.matrix_builder.add_samples_incremental(
            new_samples,
            new_sample_names
        )

        if success:
            # 输出更新后的统计信息|Output updated statistics
            stats = self.matrix_builder.get_matrix_statistics()

            self.logger.info("更新后的矩阵统计|Updated matrix statistics:")
            self.logger.info(f"  k-mer数量|K-mer count: {format_number(stats.get('n_kmers', 0))}")
            self.logger.info(f"  样本数量|Sample count: {stats.get('n_samples', 0)}")
            self.logger.info(f"  k-mer长度|K-mer size: {stats.get('kmer_size', 0)}")
            self.logger.info(f"  存储类型|Storage type: {stats.get('storage_type', 'unknown')}")

            if 'sparsity_percent' in stats:
                self.logger.info(f"  稀疏度|Sparsity: {stats['sparsity_percent']}")

            # 自动导出TSV文件（默认开启，除非指定 --no-export）|Auto export TSV files (default on, unless --no-export specified)
            no_export = getattr(self, 'no_export', False)
            if not no_export:
                self.logger.info("自动导出TSV文件|Auto exporting TSV files")
                abundance_file = str(self.config.output_path / 'abundance_matrix.tsv')
                presence_file = str(self.config.output_path / 'presence_matrix.tsv')
                self._auto_export_tsv(abundance_file, presence_file)

        return success

    def run_export(self):
        """运行矩阵导出模式|Run matrix export mode"""
        self.logger.info("开始矩阵导出模式|Starting matrix export mode")

        # 获取导出参数|Get export parameters
        output_file = getattr(self, 'export_output_file', 'abundance_matrix.tsv')
        export_format = getattr(self, 'export_format', 'sparse')
        min_abundance = getattr(self, 'export_min_abundance', 1)

        self.logger.info(f"输出文件|Output file: {output_file}")
        self.logger.info(f"导出格式|Export format: {export_format}")
        self.logger.info(f"最小丰度|Min abundance: {min_abundance}")

        # 如果指定了input_dir，使用它作为矩阵所在目录
        # If input_dir specified, use it as matrix directory
        if self.config.input_dir:
            from pathlib import Path
            custom_matrix_path = Path(self.config.input_dir)
            self.logger.info(f"使用自定义矩阵目录|Using custom matrix directory: {custom_matrix_path}")
            # 临时修改output_path|Temporarily modify output_path
            original_output_path = self.exporter.config.output_path
            self.exporter.config.output_path = custom_matrix_path
            self.exporter.matrix_file = custom_matrix_path / 'abundance_matrix.h5'
            self.exporter.kmer_dict_file = custom_matrix_path / 'kmer_dictionary.h5'

        success = self.exporter.export_to_tsv(output_file, export_format, min_abundance)

        # 如果修改了output_path，恢复原值|If modified output_path, restore original value
        if self.config.input_dir:
            self.exporter.config.output_path = original_output_path
            self.exporter.matrix_file = original_output_path / 'abundance_matrix.h5'
            self.exporter.kmer_dict_file = original_output_path / 'kmer_dictionary.h5'

        return success

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            if self.config.mode == 'count':
                success = self.run_count()
            elif self.config.mode == 'matrix':
                success = self.run_matrix()
            elif self.config.mode == 'query':
                query_fasta = getattr(self, 'query_fasta', None)
                output_file = getattr(self, 'output_file', None)

                if not query_fasta:
                    self.logger.error("查询模式需要指定FASTA文件|Query mode requires FASTA file")
                    return False

                success = self.run_query(query_fasta, output_file)
            elif self.config.mode == 'add':
                success = self.run_add()
            elif self.config.mode == 'export':
                success = self.run_export()
            else:
                self.logger.error(f"无效的操作模式|Invalid operation mode: {self.config.mode}")
                return False

            if not success:
                sys.exit(1)

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='KMC k-mer分析工具(模块化版本)|KMC K-mer Analysis Tool (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 操作模式|Operation mode
    parser.add_argument('-m', '--mode',
                       choices=['count', 'matrix', 'query', 'add', 'export'],
                       default='count',
                       help='操作模式|Operation mode:\n'
                            'count: 统计k-mer|count k-mers\n'
                            'matrix: 构建丰度矩阵|build abundance matrix\n'
                            'query: 查询k-mer|query k-mer\n'
                            'add: 添加新样本|add new samples\n'
                            'export: 导出矩阵为TSV|export matrix to TSV')

    # 输入参数|Input parameters
    parser.add_argument('-d', '--input-dir',
                       help='输入目录(自动识别双末端测序)|Input directory (auto-detect paired-end)')
    parser.add_argument('-i', '--input',
                       nargs='+',
                       help='输入文件列表(FASTQ/FASTA)|Input files list (FASTQ/FASTA)')
    parser.add_argument('-n', '--sample-names',
                       nargs='+',
                       help='样本名称(默认使用文件名)|Sample names (default: use filename)')

    # 双末端测序参数|Paired-end parameters
    parser.add_argument('--read1-suffix',
                       default='_1.clean.fq.gz',
                       help='Read1文件后缀|Read1 file suffix')
    parser.add_argument('--read2-suffix',
                       default='_2.clean.fq.gz',
                       help='Read2文件后缀|Read2 file suffix')
    parser.add_argument('--single-end',
                       action='store_true',
                       help='单末端测序模式|Single-end sequencing mode')

    # 核心参数|Core parameters
    parser.add_argument('-k', '--kmer-size',
                       type=int,
                       default=21,
                       help='k-mer大小|k-mer size')
    parser.add_argument('--min-count',
                       type=int,
                       default=2,
                       help='最小计数阈值|Minimum count threshold')
    parser.add_argument('--max-count',
                       type=int,
                       help='最大计数阈值|Maximum count threshold')

    # 路径参数|Path parameters
    parser.add_argument('--kmc-path',
                       default='~/miniforge3/envs/kmc_v.3.2.4/bin',
                       help='KMC软件路径|KMC software path')
    parser.add_argument('-o', '--output-dir',
                       default='./kmc_output',
                       help='输出目录|Output directory')
    parser.add_argument('--tmp-dir',
                       default='./kmc_tmp',
                       help='临时文件目录|Temporary directory')

    # 处理参数|Processing parameters
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=12,
                       help='线程数|Number of threads')
    parser.add_argument('--memory-limit',
                       help='内存限制(如: 12G)|Memory limit (e.g., 12G)')
    parser.add_argument('--max-memory',
                       type=int,
                       default=500,
                       help='最大内存使用量(GB)|Maximum memory usage (GB)')

    # 矩阵参数|Matrix parameters
    parser.add_argument('--matrix-format',
                       choices=['hdf5', 'tsv', 'sqlite'],
                       default='hdf5',
                       help='矩阵存储格式|Matrix storage format')
    parser.add_argument('--dense-matrix',
                       action='store_true',
                       help='使用密集矩阵(默认稀疏)|Use dense matrix (default: sparse)')
    parser.add_argument('--no-export',
                       action='store_true',
                       help='不自动导出TSV文件(默认自动导出丰度和存在文件)|Do not auto export TSV files (default: auto export abundance and presence files)')
    parser.add_argument('--no-keep-dump',
                       action='store_true',
                       help='不保留dump文件(默认保留到dump_files目录)|Do not keep dump files (default: keep in dump_files directory)')

    # k-mer索引参数|K-mer index parameters
    parser.add_argument('--index-mode',
                       choices=['auto', 'memory', 'db'],
                       default='auto',
                       help='k-mer索引模式|K-mer index mode:\\n'
                            'auto: 自动选择(默认)|auto: Auto select based on file size (default)\\n'
                            'memory: 使用内存索引(快但占用大内存)|memory: Use in-memory index (fast but high memory)\\n'
                            'db: 使用数据库索引(省内存但稍慢)|db: Use database index (low memory but slightly slower)')
    parser.add_argument('--index-threshold',
                       type=float,
                       default=1.0,
                       help='自动选择索引的文件大小阈值(GB)|File size threshold for auto index selection (GB)')

    # 查询参数|Query parameters
    parser.add_argument('-f', '--input-fasta',
                       required=False,
                       help='查询的FASTA文件(批量查询)|FASTA file for batch query')
    parser.add_argument('--output-file',
                       help='查询结果输出文件/导出TSV文件|Query result output file / Export TSV file')

    # 导出参数|Export parameters
    parser.add_argument('--format',
                       choices=['full', 'sparse'],
                       default='sparse',
                       help='导出格式|Export format: full (完整矩阵|full matrix) or sparse (稀疏格式|sparse format)')
    parser.add_argument('--min-abundance',
                       type=int,
                       default=1,
                       help='最小丰度阈值|Minimum abundance threshold')

    args = parser.parse_args()

    # 验证输入参数|Validate input parameters
    if args.mode == 'count' and not args.input_dir and not args.input:
        parser.error("count模式需要指定-d/--input-dir或-i/--input|Count mode requires -d/--input-dir or -i/--input")

    # 验证查询参数|Validate query parameters
    if args.mode == 'query' and not args.input_fasta:
        parser.error("查询模式需要指定-i/--input-fasta|Query mode requires -i/--input-fasta")

    # 创建管理器并运行|Create manager and run
    manager = KMCManager(
        mode=args.mode,
        input_dir=args.input_dir,
        input_files=args.input,
        sample_names=args.sample_names,
        read1_suffix=args.read1_suffix,
        read2_suffix=args.read2_suffix,
        single_end=args.single_end,
        kmer_size=args.kmer_size,
        min_count=args.min_count,
        max_count=args.max_count,
        kmc_path=args.kmc_path,
        output_dir=args.output_dir,
        tmp_dir=args.tmp_dir,
        threads=args.threads,
        memory_limit=args.memory_limit,
        max_memory=args.max_memory,
        matrix_format=args.matrix_format,
        sparse_storage=not args.dense_matrix,
        index_mode=args.index_mode,
        index_threshold_gb=args.index_threshold,
        keep_dump=not args.no_keep_dump  # 默认True，--no-keep-dump时为False
    )

    # 设置查询k-mer|Set query k-mer
    if args.mode == 'query':
        if args.input_fasta:
            manager.query_fasta = args.input_fasta

        if args.output_file:
            manager.output_file = args.output_file

    # 设置导出参数|Set export parameters
    if args.mode == 'export':
        manager.export_output_file = args.output_file
        manager.export_format = args.format
        manager.export_min_abundance = args.min_abundance

    # 设置矩阵自动导出TSV参数|Set matrix auto export TSV parameter
    if args.mode in ['matrix', 'add'] and args.no_export:
        manager.no_export = args.no_export

    manager.run_analysis()


if __name__ == "__main__":
    main()
