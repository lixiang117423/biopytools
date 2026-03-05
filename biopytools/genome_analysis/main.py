"""
基因组分析主模块|Genome Analysis Main Module
"""

import os
import sys
import time
import subprocess
import argparse
import tempfile
import yaml
from pathlib import Path
from datetime import datetime
from .config import GenomeAnalysisConfig
from .utils import GenomeAnalysisLogger, GenomeScopeRunner, SmudgeplotRunner, SampleFinder


def find_fastq_files(input_path: str, logger) -> list:
    """
    查找FASTQ文件|Find FASTQ files

    Args:
        input_path: 输入路径(文件或目录)|Input path (file or directory)
        logger: 日志对象|Logger object

    Returns:
        FASTQ文件列表|List of FASTQ files
    """
    # 检查输入是文件还是目录
    if os.path.isfile(input_path):
        logger.info(f"输入为单个FASTQ文件: {input_path}|Input is a single FASTQ file: {input_path}")
        return [input_path]

    # 如果是目录，递归搜索所有FASTQ文件
    logger.info(f"正在目录中搜索FASTQ文件: {input_path}|Searching FASTQ files in directory: {input_path}")

    fastq_files = list(Path(input_path).rglob('*.fastq')) + \
                  list(Path(input_path).rglob('*.fq')) + \
                  list(Path(input_path).rglob('*.fastq.gz')) + \
                  list(Path(input_path).rglob('*.fq.gz'))

    if not fastq_files:
        logger.error(f"未找到任何FASTQ文件|No FASTQ files found: {input_path}")
        return []

    fastq_files = [str(f) for f in fastq_files]
    logger.info(f"找到{len(fastq_files)}个FASTQ文件|Found {len(fastq_files)} FASTQ files")

    for f in fastq_files[:5]:  # 只显示前5个|Only show first 5
        logger.info(f"  - {f}")
    if len(fastq_files) > 5:
        logger.info(f"  ... 还有{len(fastq_files) - 5}个文件|... and {len(fastq_files) - 5} more files")

    return fastq_files


def setup_temp_dir(output_dir: str, logger) -> str:
    """
    设置临时目录|Setup temporary directory

    Args:
        output_dir: 输出目录|Output directory
        logger: 日志对象|Logger object

    Returns:
        临时目录路径|Temporary directory path
    """
    # 在输出目录下创建临时文件夹
    temp_dir = os.path.join(output_dir, 'tmp')
    os.makedirs(temp_dir, exist_ok=True)

    # 设置环境变量|Set environment variables
    os.environ['TMPDIR'] = temp_dir
    os.environ['TMP'] = temp_dir
    os.environ['TEMP'] = temp_dir

    logger.info(f"设置临时目录: {temp_dir}|Set temporary directory: {temp_dir}")

    return temp_dir


def check_dependencies(logger, run_smudgeplot=False) -> bool:
    """
    检查依赖环境|Check dependencies

    Args:
        logger: 日志对象|Logger object
        run_smudgeplot: 是否运行Smudgeplot|Whether to run Smudgeplot

    Returns:
        是否所有依赖都可用|Whether all dependencies are available
    """
    logger.info("=" * 60)
    logger.info("检查依赖环境|Checking dependencies")
    logger.info("=" * 60)

    dependencies = [
        ('jellyfish', 'jellyfish'),
    ]

    # 如果需要运行Smudgeplot，添加相关依赖检查|Add Smudgeplot related dependencies if needed
    if run_smudgeplot:
        dependencies.append(('FastK', 'FastK'))
        dependencies.append(('smudgeplot', 'smudgeplot'))

    all_ok = True
    for cmd, name in dependencies:
        try:
            # 检查命令|Check command
            subprocess.run(['which', cmd], capture_output=True, check=True)
            logger.info(f"[OK] {name} 已找到|{name} found")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"[ERROR] {name} 未找到|{name} not found")
            all_ok = False

    return all_ok


def generate_software_versions_yml(output_dir: str, config: GenomeAnalysisConfig, logger):
    """
    生成software_versions.yml文件|Generate software_versions.yml file

    Args:
        output_dir: 输出目录|Output directory
        config: 配置对象|Configuration object
        logger: 日志对象|Logger object
    """
    logger.info("生成软件版本信息文件|Generating software versions file")

    # 定义工具列表|Define tools list
    tools = {}
    tool_commands = {
        'jellyfish': ['jellyfish', '--version'],
        'python3': ['python3', '--version'],
    }

    # 如果运行Smudgeplot，添加相关工具|Add Smudgeplot related tools if running
    if not config.skip_smudgeplot:
        tool_commands.update({
            'FastK': ['FastK', '-v'],  # FastK -v 不输出版本，会返回空
            'smudgeplot': ['smudgeplot', '--version'],
        })

    # 获取工具版本|Get tool versions
    for tool_name, cmd in tool_commands.items():
        try:
            # 同时捕获 stdout 和 stderr|Capture both stdout and stderr
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=10
            )

            # 优先从 stdout 获取，其次 stderr|Prefer stdout, fallback to stderr
            output = result.stdout.strip() if result.stdout.strip() else result.stderr.strip()
            version = output.split('\n')[0] if output else 'unknown'

            # 尝试获取路径|Try to get path
            path_result = subprocess.run(
                ['which', tool_name],
                capture_output=True,
                text=True,
                timeout=5
            )
            tool_path = path_result.stdout.strip() if path_result.returncode == 0 else 'not found'

            tools[tool_name] = {
                'version': version,
                'path': tool_path
            }

            # 记录调试信息|Log debug info
            if version == 'unknown':
                logger.debug(f"无法获取 {tool_name} 版本: stdout='{result.stdout.strip()}', stderr='{result.stderr.strip()}'")

        except subprocess.TimeoutExpired:
            tools[tool_name] = {
                'version': 'timeout',
                'path': 'not found'
            }
            logger.warning(f"获取 {tool_name} 版本超时")
        except FileNotFoundError:
            tools[tool_name] = {
                'version': 'not installed',
                'path': 'not found'
            }
            logger.warning(f"{tool_name} 未安装")
        except Exception as e:
            tools[tool_name] = {
                'version': 'error',
                'path': 'not found'
            }
            logger.warning(f"获取 {tool_name} 版本时出错: {str(e)}")

    # 构建版本信息|Build version info
    version_info = {
        'pipeline': {
            'name': 'biopytools genomescope',
            'version': '1.0.0',
            'command': ' '.join(sys.argv)
        },
        'tools': tools,
        'parameters': {
            'kmer_size': config.kmer_size,
            'read_length': config.read_length,
            'threads': config.threads,
            'hash_size': config.hash_size,
            'max_kmer_cov': config.max_kmer_cov,
            'skip_smudgeplot': config.skip_smudgeplot
        },
        'execution': {
            'timestamp': datetime.now().isoformat(),
            'hostname': os.uname().nodename
        }
    }

    # 创建输出目录|Create output directory
    # 使用00前缀表示流程元数据目录（与99_logs对应）|Use 00 prefix for pipeline metadata (mirrors 99_logs)
    pipeline_info_dir = os.path.join(output_dir, '00_pipeline_info')
    os.makedirs(pipeline_info_dir, exist_ok=True)

    # 写入YAML文件|Write YAML file
    output_file = os.path.join(pipeline_info_dir, 'software_versions.yml')
    with open(output_file, 'w') as f:
        yaml.dump(version_info, f, default_flow_style=False, sort_keys=False)

    logger.info(f"软件版本信息已保存|Software versions saved to: {output_file}")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='基因组分析工具|Genome Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i fastq_dir -o output_dir -l 150
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='输入FASTQ文件或目录|Input FASTQ file or directory')
    required.add_argument('-o', '--output-dir', required=True,
                         help='输出目录|Output directory')
    required.add_argument('-l', '--read-length', required=True, type=int,
                         help='测序读长|Read length')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|optional arguments')
    optional.add_argument('-k', '--kmer-size', type=int, default=21,
                         help='K-mer大小|K-mer size')
    optional.add_argument('-t', '--threads', type=int, default=64,
                         help='线程数|Number of threads')
    optional.add_argument('-s', '--hash-size', default='10G',
                         help='Jellyfish哈希表大小|Jellyfish hash size')
    optional.add_argument('-c', '--max-kmer-cov', type=int, default=1000,
                         help='最大k-mer覆盖度|Max k-mer coverage')
    optional.add_argument('--skip-smudgeplot', action='store_true',
                         help='跳过Smudgeplot倍性分析|Skip Smudgeplot ploidy analysis')
    optional.add_argument('--fastk-table',
                         default='',
                         help='FastK表文件路径|FastK table file path')
    optional.add_argument('--fastk-memory', default='16G',
                         help='FastK内存大小|FastK memory size')
    optional.add_argument('--read1-suffix', default='*_1.clean.fq.gz',
                         help='Read1文件后缀模式|Read1 file suffix pattern')

    args = parser.parse_args()

    # 创建输出目录|Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # 设置顶层控制台日志（样本循环内会重新设置）|Setup top-level console logging (will reset in sample loop)
    # INFO→stdout→.out, WARNING→stderr→.err
    import logging
    logger = logging.getLogger("GenomeAnalysis")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # stdout handler - INFO级别|stdout handler - INFO level
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    # stderr handler - WARNING及以上|stderr handler - WARNING and above
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Genome Analysis Pipeline")
        logger.info("=" * 60)
        logger.info(f"输入路径: {args.input}")
        logger.info(f"输出目录: {args.output_dir}")
        logger.info(f"读长: {args.read_length}")
        logger.info(f"K-mer大小: {args.kmer_size}")
        logger.info(f"线程数: {args.threads}")

        # 检查依赖|Check dependencies
        if not check_dependencies(logger, not args.skip_smudgeplot):
            logger.error("依赖检查失败，请安装所需软件")
            sys.exit(1)

        # 初始化配置|Initialize config
        config = GenomeAnalysisConfig(
            input_dir=args.input,
            output_dir=args.output_dir,
            read_length=args.read_length,
            kmer_size=args.kmer_size,
            threads=args.threads,
            hash_size=args.hash_size,
            max_kmer_cov=args.max_kmer_cov,
            read1_suffix=args.read1_suffix,
            skip_smudgeplot=args.skip_smudgeplot
        )

        config.validate()

        # 使用SampleFinder查找样品|Use SampleFinder to find samples
        sample_finder = SampleFinder(logger, config.read1_suffix)
        samples = sample_finder.find_samples(args.input)

        if not samples:
            logger.error("未找到任何样品")
            sys.exit(1)

        # 处理每个样品|Process each sample
        successful_count = 0
        failed_count = 0

        for sample_name, fastq_files in samples:
            # 为每个样品创建子目录|Create subdirectory for each sample
            sample_output_dir = os.path.join(args.output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            # 为每个样本设置独立的日志系统|Setup independent logging for each sample
            logger_manager = GenomeAnalysisLogger(sample_output_dir)
            logger = logger_manager.get_logger()

            # 生成样本级的软件版本信息文件|Generate sample-level software versions file
            generate_software_versions_yml(sample_output_dir, config, logger)

            logger.info("")
            logger.info("=" * 60)
            logger.info(f"处理样品|Processing sample: {sample_name}")
            logger.info("=" * 60)

            # 设置临时目录到样品目录下|Setup temp directory under sample directory
            setup_temp_dir(sample_output_dir, logger)

            # 使用相对路径，避免包含父目录的中文字符导致R无法处理|Use relative path to avoid Chinese characters in parent directory
            # 不使用 os.path.abspath()，因为它会生成包含中文的绝对路径|Don't use os.path.abspath() as it generates absolute paths with Chinese
            output_dir_abs = sample_output_dir

            # 创建各步骤的输出目录（添加数字前缀表示流程顺序）|Create output directories (with number prefix)
            jellyfish_dir = os.path.join(output_dir_abs, "01_jellyfish")
            os.makedirs(jellyfish_dir, exist_ok=True)
            # 文件命名：样本名.工具名|File naming: sample.tool
            output_prefix = os.path.join(jellyfish_dir, f"{sample_name}.jellyfish")

            # FastK目录：无编号前缀，因为：
            # 1) 仅在运行Smudgeplot时需要（条件性步骤）
            # 2) 支持通过--fastk-table使用预生成表
            # 3) 属于Smudgeplot预处理，非独立分析步骤
            # FastK directory: No number prefix because:
            # 1) Only needed when running Smudgeplot (conditional step)
            # 2) Supports pre-generated tables via --fastk-table
            # 3) Preprocessing for Smudgeplot, not an independent analysis step
            fastk_dir = os.path.join(output_dir_abs, "fastk")
            # FastK目录会在需要时创建|FastK directory will be created when needed

            # 运行GenomeScope流程|Run GenomeScope pipeline
            gs_runner = GenomeScopeRunner(logger)

            # 步骤1: Jellyfish count
            if not gs_runner.run_jellyfish_count(
                [str(f) for f in fastq_files], output_prefix, config.kmer_size,
                config.hash_size, config.threads
            ):
                logger.error(f"样品 {sample_name} Jellyfish count失败")
                failed_count += 1
                continue

            # 步骤2: Jellyfish histo
            jf_file = f"{output_prefix}.jf"
            histo_file = f"{output_prefix}.histo"
            if not gs_runner.run_jellyfish_histo(jf_file, histo_file, config.threads):
                logger.error(f"样品 {sample_name} Jellyfish histo失败")
                failed_count += 1
                continue

            # 步骤3: GenomeScope
            # histo_file 已在上面定义|histo_file is defined above
            gs_output_dir = os.path.join(output_dir_abs, "02_genomescope")
            os.makedirs(gs_output_dir, exist_ok=True)

            kcov = gs_runner.run_genomescope(
                histo_file, config.kmer_size, config.read_length,
                gs_output_dir, config.max_kmer_cov
            )

            if kcov is None:
                logger.warning(f"样品 {sample_name} 未能获取kcov值，将尝试使用默认值继续")
                kcov = 50.0  # 默认值

            config.kcov = kcov

            # 步骤4-5: Smudgeplot (默认运行，可使用--skip-smudgeplot跳过)
            sm_runner = SmudgeplotRunner(logger)

            if not args.skip_smudgeplot:
                logger.info("")
                logger.info("=" * 60)
                logger.info(f"样品 {sample_name} 运行Smudgeplot倍性分析")
                logger.info("=" * 60)

                # 确定FastK表路径|Determine FastK table path
                if args.fastk_table:
                    fastk_table = args.fastk_table
                    logger.info(f"使用已存在的FastK表: {fastk_table}")
                else:
                    # 运行FastK生成FastK表
                    os.makedirs(fastk_dir, exist_ok=True)
                    fastk_table = os.path.join(fastk_dir, "fastk_table")
                    logger.info(f"将生成FastK表: {fastk_table}")

                    if not sm_runner.run_fastk(
                        [str(f) for f in fastq_files], fastk_table, config.kmer_size,
                        config.threads, args.fastk_memory, sample_name
                    ):
                        logger.error(f"样品 {sample_name} FastK步骤失败")
                        failed_count += 1
                        continue

                # 运行Smudgeplot
                smudgeplot_output = os.path.join(output_dir_abs, "03_smudgeplot")
                if not sm_runner.run_smudgeplot(
                    fastk_table, smudgeplot_output, sample_name, kcov,
                    config.kmer_size, config.threads
                ):
                    logger.warning(f"样品 {sample_name} Smudgeplot步骤失败（可能是纯合二倍体，属于正常现象）")
                    # Smudgeplot失败不影响整体成功状态
                    # Smudgeplot failure doesn't affect overall success

                else:
                    logger.info(f"样品 {sample_name} Smudgeplot输出: {smudgeplot_output}")
            else:
                logger.info("")
                logger.info("=" * 60)
                logger.info(f"样品 {sample_name} 跳过Smudgeplot倍性分析")
                logger.info("=" * 60)

            logger.info(f"样品 {sample_name} 处理完成")
            successful_count += 1

        # 完成统计|Completion statistics
        elapsed_time = time.time() - start_time

        logger.info("")
        logger.info("=" * 60)
        logger.info("Pipeline Summary")
        logger.info("=" * 60)
        logger.info(f"总运行时间: {elapsed_time:.2f} seconds")
        logger.info(f"总样品数: {len(samples)}")
        logger.info(f"成功处理: {successful_count}")
        logger.info(f"失败样品: {failed_count}")
        if len(samples) > 0:
            logger.info(f"成功率: {(successful_count/len(samples))*100:.1f}%")
        logger.info("Pipeline completed successfully")

    except KeyboardInterrupt:
        logger.warning("用户中断操作")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
        sys.exit(1)


class GenomeAnalysis:
    """基因组分析类|Genome Analysis Class"""

    def __init__(self, **kwargs):
        """
        初始化分析器|Initialize analyzer

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        self.config = GenomeAnalysisConfig(**kwargs)
        self.config.validate()

        # 创建输出目录|Create output directory
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 初始化日志|Initialize logging
        self.logger_manager = GenomeAnalysisLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 生成软件版本信息文件|Generate software versions file
        generate_software_versions_yml(self.config.output_dir, self.config, self.logger)

        # 初始化工具类|Initialize utility classes
        self.gs_runner = GenomeScopeRunner(self.logger)
        self.sm_runner = SmudgeplotRunner(self.logger)

    def run(self):
        """
        运行分析流程|Run analysis pipeline

        Returns:
            是否成功|Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("Genome Analysis Pipeline")
            self.logger.info("=" * 60)
            self.logger.info(f"配置信息:\n{self.config}")

            # 检查依赖|Check dependencies
            if not check_dependencies(self.logger):
                self.logger.error("依赖检查失败")
                return False

            # 查找FASTQ文件|Find FASTQ files
            fastq_files = find_fastq_files(self.config.input_dir, self.logger)
            if not fastq_files:
                self.logger.error("未找到FASTQ文件")
                return False

            # 获取绝对路径|Get absolute paths
            output_dir_abs = os.path.abspath(self.config.output_dir)

            # 从第一个FASTQ文件中提取样本名|Extract sample name from first FASTQ file
            import os
            first_fastq = fastq_files[0]
            sample_name = Path(first_fastq).stem
            # 移除常见的后缀|Remove common suffixes
            for suffix in ['_R1', '_R2', '_1', '_2', '.fastq', '.fq']:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)]
                    break

            # 创建各步骤的输出目录（添加数字前缀表示流程顺序）|Create output directories (with number prefix)
            jellyfish_dir = os.path.join(output_dir_abs, "01_jellyfish")
            os.makedirs(jellyfish_dir, exist_ok=True)
            # 文件命名：样本名.工具名|File naming: sample.tool
            output_prefix = os.path.join(jellyfish_dir, f"{sample_name}.jellyfish")

            # 运行GenomeScope流程|Run GenomeScope pipeline
            # 步骤1: Jellyfish count
            if not self.gs_runner.run_jellyfish_count(
                fastq_files, output_prefix, self.config.kmer_size,
                self.config.hash_size, self.config.threads
            ):
                return False

            # 步骤2: Jellyfish histo
            jf_file = f"{output_prefix}.jf"
            histo_file = f"{output_prefix}.histo"
            if not self.gs_runner.run_jellyfish_histo(jf_file, histo_file, self.config.threads):
                return False

            # 步骤3: GenomeScope
            # histo_file 已在上面定义|histo_file is defined above
            gs_output_dir = os.path.join(output_dir_abs, "02_genomescope")
            os.makedirs(gs_output_dir, exist_ok=True)

            kcov = self.gs_runner.run_genomescope(
                histo_file, self.config.kmer_size, self.config.read_length,
                gs_output_dir, self.config.max_kmer_cov
            )

            if kcov is None:
                self.logger.warning("未能获取kcov值")
                return False

            self.config.kcov = kcov

            # 完成统计|Completion statistics
            elapsed_time = time.time() - start_time

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间: {elapsed_time:.2f} seconds")
            self.logger.info(f"K-mer coverage: {kcov}")
            self.logger.info(f"GenomeScope输出: {gs_output_dir}")
            self.logger.info("Pipeline completed successfully")

            return True

        except Exception as e:
            self.logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
            return False


if __name__ == "__main__":
    main()
