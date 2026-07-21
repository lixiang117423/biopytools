"""
基因组分析主模块|Genome Analysis Main Module
"""

import os
import sys
import time
from typing import Optional
import subprocess
import argparse
import tempfile
import yaml
import glob
import shutil
from pathlib import Path
from datetime import datetime
try:
    from openpyxl import Workbook, Font
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

from .config import GenomeAnalysisConfig
from .utils import GenomeAnalysisLogger, GenomeScopeRunner, SmudgeplotRunner, SampleFinder, get_conda_env, build_conda_command


def _parse_genomescope_model(model_file: str) -> Optional[float]:
    """
    从GenomeScope model.txt中提取杂合度点估计值(r参数)|Extract heterozygosity point estimate (r param) from model.txt

    Args:
        model_file: model.txt 文件路径

    Returns:
        杂合度值(小数)或None|Heterozygosity value (decimal) or None
    """
    if not os.path.exists(model_file):
        return None
    try:
        import re
        with open(model_file, 'r') as f:
            for line in f:
                # 匹配 r 或 r1 参数行: "r1      1.742e-01  ..."
                if re.match(r'^r\d?\s', line.strip()):
                    match = re.search(r'^r\d?\s+([\d.eE+-]+)', line.strip())
                    if match:
                        return float(match.group(1))
        return None
    except Exception:
        return None


def _parse_genomescope_summary(summary_file: str) -> dict:
    """
    从GenomeScope summary.txt中解析关键指标|Parse key metrics from GenomeScope summary.txt

    Args:
        summary_file: summary.txt 文件路径

    Returns:
        包含各指标的字典或None|Dict of metrics or None
    """
    if not os.path.exists(summary_file):
        return None

    metrics = {}
    try:
        import re
        with open(summary_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('GenomeScope') or line.startswith('input') or \
                   line.startswith('output') or line.startswith('p =') or line.startswith('k =') or \
                   line.startswith('max_kmercov') or line.startswith('property'):
                    continue

                # 找到最后两个以数字开头的token作为min和max值
                # Find the last two tokens starting with a digit as min and max values
                parts = line.split()
                value_indices = []
                for i, t in enumerate(parts):
                    if re.match(r'^[\d]+', t):
                        value_indices.append(i)

                if len(value_indices) < 2:
                    continue

                min_idx = value_indices[-2]
                max_idx = value_indices[-1]
                prop_name = ' '.join(parts[:min_idx])

                # 去除数值中的逗号和单位|Remove commas and units from values
                min_val = re.sub(r'^([\d,.]+).*$', r'\1', parts[min_idx])
                max_val = re.sub(r'^([\d,.]+).*$', r'\1', parts[max_idx])

                key = prop_name.lower().replace(' ', '_').replace('(', '').replace(')', '')
                metrics[f'{key}_min'] = min_val
                metrics[f'{key}_max'] = max_val
        return metrics if metrics else None
    except Exception:
        return None


def generate_summary_table(samples, output_dir: str, sample_kcov_dict: dict, logger):
    """
    生成多样本汇总表|Generate multi-sample summary table

    Args:
        samples: 样本列表 [(sample_name, fastq_files), ...]
        output_dir: 输出目录
        sample_kcov_dict: {sample_name: kcov} 字典
        logger: 日志对象
    """
    logger.info("")
    logger.info("=" * 60)
    logger.info("生成汇总表|Generating summary table")
    logger.info("=" * 60)

    all_dir = os.path.join(output_dir, "all")
    os.makedirs(all_dir, exist_ok=True)
    summary_tsv = os.path.join(all_dir, "summary.tsv")

    header = [
        "Sample",
        "Genome_Haploid_Length(bp)",
        "Genome_Repeat_Length(bp)",
        "Genome_Unique_Length(bp)",
        "Heterozygous",
        "Read_Error_Rate",
        "Model_Fit(%)",
        "Kmer_Coverage",
    ]

    rows = []
    parsed_count = 0
    for sample_tuple in samples:
        sample_name = sample_tuple[0]
        summary_file = os.path.join(output_dir, sample_name, "02_genomescope", "summary.txt")
        metrics = _parse_genomescope_summary(summary_file)

        if metrics:
            parsed_count += 1
            kcov = sample_kcov_dict.get(sample_name, '')

            # 定义范围格式化:必须在 het_str 之前定义,否则 het_point 为 None 时
            # 会引用尚未定义的 fmt_range -> UnboundLocalError(原 bug)
            # Define range formatter BEFORE het_str: otherwise het_point=None references
            # an undefined fmt_range -> UnboundLocalError (the original bug)
            def fmt_range(key_base):
                mn = metrics.get(f'{key_base}_min', '')
                mx = metrics.get(f'{key_base}_max', '')
                if mn == mx:
                    return mn
                return f"{mn}-{mx}" if mn and mx else ''

            # 从model.txt提取杂合度点估计值(r参数)|Extract heterozygosity point estimate (r param) from model.txt
            model_file = os.path.join(output_dir, sample_name, "02_genomescope", "model.txt")
            het_point = _parse_genomescope_model(model_file)
            het_str = f"{het_point * 100:.4f}%" if het_point is not None else fmt_range('heterozygous_ab')

            row = [
                sample_name,
                fmt_range('genome_haploid_length'),
                fmt_range('genome_repeat_length'),
                fmt_range('genome_unique_length'),
                het_str,
                fmt_range('read_error_rate'),
                fmt_range('model_fit'),
                f"{kcov:.1f}" if kcov else '',
            ]
        else:
            row = [sample_name] + [''] * (len(header) - 1)
            logger.warning(f"  未找到summary.txt或解析失败|summary.txt not found or parse failed: {sample_name}")
        rows.append(row)

    with open(summary_tsv, 'w', encoding='utf-8') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

    logger.info(f"汇总表已保存|Summary table saved to: {summary_tsv}")

    # 同时输出xlsx|Also output xlsx
    if HAS_OPENPYXL:
        xlsx_file = summary_tsv.replace('.tsv', '.xlsx')
        wb = Workbook()
        ws = wb.active
        ws.title = "Summary"
        for cell in ws[1]:
            # openpyxl 新版 .copy() 已弃用,改用 Font 相加|newer openpyxl deprecates .copy(), use Font addition
            cell.font = cell.font + Font(bold=True)
        ws.append(header)
        for row in rows:
            ws.append(row)
        # 自动调整列宽|Auto-adjust column widths
        for col in ws.columns:
            max_len = max(len(str(c.value or '')) for c in col)
            ws.column_dimensions[col[0].column_letter].width = max_len + 2
        wb.save(xlsx_file)
        logger.info(f"汇总表已保存|Summary table saved to: {xlsx_file}")
    else:
        logger.warning("openpyxl未安装，跳过xlsx输出|openpyxl not installed, skipping xlsx output")

    logger.info(f"成功解析: {parsed_count}/{len(samples)} 个样本|Successfully parsed: {parsed_count}/{len(samples)} samples")


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
        ('genomescope2', 'GenomeScope 2.0'),
    ]

    # 如果需要运行Smudgeplot，添加相关依赖检查|Add Smudgeplot related dependencies if needed
    if run_smudgeplot:
        dependencies.append(('FastK', 'FastK'))
        dependencies.append(('smudgeplot', 'smudgeplot'))

    all_ok = True
    for cmd, name in dependencies:
        # 统一检测：conda环境优先，PATH兜底，与运行时调用方式(conda run/直接调用)一致
        # Unified detection: conda env first, PATH fallback, consistent with runtime
        env = get_conda_env(cmd)
        if env:
            logger.info(f"[OK] {name} 已找到 (conda环境: {env})|{name} found (conda env: {env})")
        elif shutil.which(cmd):
            logger.info(f"[OK] {name} 已找到 (PATH)|{name} found (in PATH)")
        else:
            logger.error(f"[ERROR] {name} 未找到|{name} not found")
            logger.error(f"请安装|Please install: conda create -n <env> -c bioconda {cmd}")
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
        'R': ['R', '--version'],
        'Rscript': ['Rscript', '--version'],
    }

    # 如果运行Smudgeplot，添加相关工具|Add Smudgeplot related tools if running
    if not config.skip_smudgeplot:
        tool_commands.update({
            'FastK': ['FastK', '-v'],  # FastK -v 不输出版本，会返回空
            'smudgeplot': ['smudgeplot', '--version'],
        })

    # 获取工具版本(经build_conda_command,保证仅在conda环境里的工具也能取到版本)|
    # Get versions via build_conda_command so conda-only tools are also detected
    for tool_name, cmd in tool_commands.items():
        try:
            wrapped_cmd = build_conda_command(cmd[0], cmd[1:])
            # 同时捕获 stdout 和 stderr|Capture both stdout and stderr
            result = subprocess.run(
                wrapped_cmd,
                capture_output=True,
                text=True,
                timeout=10
            )

            # 优先从 stdout 获取，其次 stderr|Prefer stdout, fallback to stderr
            output = result.stdout.strip() if result.stdout.strip() else result.stderr.strip()
            version = output.split('\n')[0] if output else 'unknown'

            # 路径:优先conda env,其次PATH|Path: prefer conda env, then PATH
            env = get_conda_env(cmd[0])
            if env:
                tool_path = f"conda:{env}"
            else:
                tool_path = shutil.which(cmd[0]) or 'not found'

            tools[tool_name] = {
                'version': version,
                'path': tool_path
            }

            # 记录调试信息|Log debug info
            if version == 'unknown':
                logger.debug(f"无法获取 {tool_name} 版本|Failed to get {tool_name} version: stdout='{result.stdout.strip()}', stderr='{result.stderr.strip()}'")

        except subprocess.TimeoutExpired:
            tools[tool_name] = {
                'version': 'timeout',
                'path': 'not found'
            }
            logger.warning(f"获取 {tool_name} 版本超时|Timeout while getting {tool_name} version")
        except FileNotFoundError:
            tools[tool_name] = {
                'version': 'not installed',
                'path': 'not found'
            }
            logger.warning(f"{tool_name} 未安装|{tool_name} not installed")
        except Exception as e:
            tools[tool_name] = {
                'version': 'error',
                'path': 'not found'
            }
            logger.warning(f"获取 {tool_name} 版本时出错|Error getting {tool_name} version: {str(e)}")

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
    optional.add_argument('--ploidy', type=int, default=2,
                         help='基因组倍性 1-6 (默认: 2，由Smudgeplot自动推断)|Genome ploidy level 1-6 (default: 2, auto-inferred by Smudgeplot)')
    optional.add_argument('--genomescope-env', default='genomescope_v.2.0.1',
                         help='GenomeScope conda环境名称 (默认: genomescope_v.2.0.1)|GenomeScope conda env name (default: genomescope_v.2.0.1)')
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
    # INFO→stdout→.out, WARNING→stderr→.err (§2.3 超算日志分离|§2.3 job scheduler log separation)
    import logging
    logger = logging.getLogger("GenomeAnalysis")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False  # 避免传播到root logger导致重复|Avoid propagation to root logger

    # 规范格式：点号分隔毫秒(§2.1)|Standard format: dot-separated ms (§2.1)
    formatter = logging.Formatter(
        '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

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
        logger.info(f"输入路径: {args.input}|Input path: {args.input}")
        logger.info(f"输出目录: {args.output_dir}|Output directory: {args.output_dir}")
        logger.info(f"读长: {args.read_length}|Read length: {args.read_length}")
        logger.info(f"K-mer大小: {args.kmer_size}|K-mer size: {args.kmer_size}")
        logger.info(f"线程数: {args.threads}|Threads: {args.threads}")

        # 检查依赖|Check dependencies
        if not check_dependencies(logger, not args.skip_smudgeplot):
            logger.error("依赖检查失败，请安装所需软件|Dependency check failed, please install required software")
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
            skip_smudgeplot=args.skip_smudgeplot,
            ploidy=args.ploidy,
            genomescope_env=args.genomescope_env
        )

        config.validate()

        # 使用SampleFinder查找样品|Use SampleFinder to find samples
        sample_finder = SampleFinder(logger, config.read1_suffix)
        samples = sample_finder.find_samples(args.input)

        if not samples:
            logger.error("未找到任何样品|No samples found")
            sys.exit(1)

        # 处理每个样品|Process each sample
        successful_count = 0
        failed_count = 0
        sample_kcov_dict = {}  # 收集每个样本的kcov|Collect kcov for each sample

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
            gs_runner = GenomeScopeRunner(logger, config)

            # 步骤1: Jellyfish count
            if not gs_runner.run_jellyfish_count(
                [str(f) for f in fastq_files], output_prefix, config.kmer_size,
                config.hash_size, config.threads
            ):
                logger.error(f"样品 {sample_name} Jellyfish count失败|Sample {sample_name} Jellyfish count failed")
                failed_count += 1
                continue

            # 步骤2: Jellyfish histo
            jf_file = f"{output_prefix}.jf"
            histo_file = f"{output_prefix}.histo"
            if not gs_runner.run_jellyfish_histo(jf_file, output_prefix, config.threads):
                logger.error(f"样品 {sample_name} Jellyfish histo失败|Sample {sample_name} Jellyfish histo failed")
                failed_count += 1
                continue

            # 步骤3: GenomeScope
            # histo_file 已在上面定义|histo_file is defined above
            gs_output_dir = os.path.join(output_dir_abs, "02_genomescope")
            os.makedirs(gs_output_dir, exist_ok=True)

            kcov = gs_runner.run_genomescope(
                histo_file, config.kmer_size, config.read_length,
                gs_output_dir, config.max_kmer_cov, config.ploidy
            )

            if kcov is None:
                logger.warning(f"样品 {sample_name} 未能获取kcov值，将尝试使用默认值继续|Sample {sample_name} failed to get kcov value, will try to continue with default value")
                kcov = 50.0  # 默认值

            config.kcov = kcov
            sample_kcov_dict[sample_name] = kcov

            # 步骤4-5: Smudgeplot (默认运行，可使用--skip-smudgeplot跳过)
            sm_runner = SmudgeplotRunner(logger)

            if not config.skip_smudgeplot:
                logger.info("")
                logger.info("=" * 60)
                logger.info(f"样品 {sample_name} 运行Smudgeplot倍性分析|Sample {sample_name} running Smudgeplot ploidy analysis")
                logger.info("=" * 60)

                # 检查Smudgeplot最终输出是否已存在|Check if Smudgeplot final output already exists
                smudgeplot_output = os.path.join(output_dir_abs, "03_smudgeplot")
                existing_plots = glob.glob(os.path.join(smudgeplot_output, f"{sample_name}*.pdf"))

                if existing_plots:
                    logger.info(f"Smudgeplot输出已存在，跳过Smudgeplot流程|Smudgeplot output exists, skipping Smudgeplot pipeline")
                    logger.info(f"  输出文件|Output files: {', '.join([os.path.basename(f) for f in existing_plots])}")
                    logger.info("如需重新运行，请删除现有输出文件|To re-run, delete existing output files")
                else:
                    # Smudgeplot输出不存在，运行完整流程|Smudgeplot output doesn't exist, run full pipeline

                    # 确定FastK表路径|Determine FastK table path
                    if args.fastk_table:
                        fastk_table = args.fastk_table
                        logger.info(f"使用已存在的FastK表: {fastk_table}|Using existing FastK table: {fastk_table}")
                    else:
                        # 运行FastK生成FastK表
                        os.makedirs(fastk_dir, exist_ok=True)
                        fastk_table = os.path.join(fastk_dir, "fastk_table")
                        logger.info(f"将生成FastK表: {fastk_table}|Will generate FastK table: {fastk_table}")

                        if not sm_runner.run_fastk(
                            [str(f) for f in fastq_files], fastk_table, config.kmer_size,
                            config.threads, args.fastk_memory, sample_name
                        ):
                            logger.error(f"样品 {sample_name} FastK步骤失败|Sample {sample_name} FastK step failed")
                            failed_count += 1
                            continue

                    # 运行Smudgeplot
                    if not sm_runner.run_smudgeplot(
                        fastk_table, smudgeplot_output, sample_name, kcov,
                        config.kmer_size, config.threads
                    ):
                        logger.warning(f"样品 {sample_name} Smudgeplot步骤失败（可能是纯合二倍体，属于正常现象）|Sample {sample_name} Smudgeplot step failed (may be homozygous diploid, which is normal)")
                        # Smudgeplot失败不影响整体成功状态
                        # Smudgeplot failure doesn't affect overall success
                    else:
                        logger.info(f"样品 {sample_name} Smudgeplot输出: {smudgeplot_output}|Sample {sample_name} Smudgeplot output: {smudgeplot_output}")
            else:
                logger.info("")
                logger.info("=" * 60)
                logger.info(f"样品 {sample_name} 跳过Smudgeplot倍性分析|Sample {sample_name} skipping Smudgeplot ploidy analysis")
                logger.info("=" * 60)

            logger.info(f"样品 {sample_name} 处理完成|Sample {sample_name} processing completed")
            successful_count += 1

        # 汇总结果文件|Collect result files
        logger.info("")
        logger.info("=" * 60)
        logger.info("汇总结果文件|Collecting result files")
        logger.info("=" * 60)

        # 使用总输出目录而不是样品目录|Use total output dir instead of sample dir
        all_dir = os.path.join(args.output_dir, "all")
        os.makedirs(all_dir, exist_ok=True)
        logger.info(f"创建汇总目录|Created summary directory: {all_dir}")

        collected_count = 0
        for sample_tuple in samples:
            # samples是tuple列表: (sample_name, fastq_files)|samples is list of tuples: (sample_name, fastq_files)
            sample_name = sample_tuple[0]
            # 拷贝GenomeScope的linear_plot.png|Copy GenomeScope linear_plot.png
            gs_linear_plot = os.path.join(args.output_dir, sample_name, "02_genomescope", "linear_plot.png")
            if os.path.exists(gs_linear_plot):
                target_file = os.path.join(all_dir, f"{sample_name}_genomescope.png")
                shutil.copy2(gs_linear_plot, target_file)
                logger.info(f"  拷贝|Copied: {sample_name}/linear_plot.png -> all/{sample_name}_genomescope.png")
                collected_count += 1
            else:
                logger.warning(f"  文件不存在|File not found: {sample_name}/02_genomescope/linear_plot.png")

            # 拷贝Smudgeplot的smudgeplot.png|Copy Smudgeplot smudgeplot.png
            sm_smudgeplot = os.path.join(args.output_dir, sample_name, "03_smudgeplot", f"{sample_name}_smudgeplot.png")
            if os.path.exists(sm_smudgeplot):
                target_file = os.path.join(all_dir, f"{sample_name}_smudgeplot.png")
                shutil.copy2(sm_smudgeplot, target_file)
                logger.info(f"  拷贝|Copied: {sample_name}/{sample_name}_smudgeplot.png -> all/{sample_name}_smudgeplot.png")
                collected_count += 1
            else:
                logger.warning(f"  文件不存在|File not found: {sample_name}/03_smudgeplot/{sample_name}_smudgeplot.png")

        logger.info(f"汇总完成|Collection completed: {collected_count} 个文件|files")

        # 生成汇总表|Generate summary table
        generate_summary_table(samples, args.output_dir, sample_kcov_dict, logger)

        # 完成统计|Completion statistics
        elapsed_time = time.time() - start_time

        logger.info("")
        logger.info("=" * 60)
        logger.info("Pipeline Summary")
        logger.info("=" * 60)
        logger.info(f"总运行时间: {elapsed_time:.2f} seconds|Total runtime: {elapsed_time:.2f} seconds")
        logger.info(f"总样品数: {len(samples)}|Total samples: {len(samples)}")
        logger.info(f"成功处理: {successful_count}|Successfully processed: {successful_count}")
        logger.info(f"失败样品: {failed_count}|Failed samples: {failed_count}")
        if len(samples) > 0:
            logger.info(f"成功率: {(successful_count/len(samples))*100:.1f}%|Success rate: {(successful_count/len(samples))*100:.1f}%")
        logger.info("Pipeline completed successfully")

    except KeyboardInterrupt:
        logger.warning("用户中断操作|User interrupted operation")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Pipeline失败|Pipeline failed: {str(e)}", exc_info=True)
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
        self.gs_runner = GenomeScopeRunner(self.logger, self.config)
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
            self.logger.info(f"配置信息|Configuration info:\n{self.config}")

            # 检查依赖|Check dependencies
            if not check_dependencies(self.logger):
                self.logger.error("依赖检查失败|Dependency check failed")
                return False

            # 查找FASTQ文件|Find FASTQ files
            fastq_files = find_fastq_files(self.config.input_dir, self.logger)
            if not fastq_files:
                self.logger.error("未找到FASTQ文件|No FASTQ files found")
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
            if not self.gs_runner.run_jellyfish_histo(jf_file, output_prefix, self.config.threads):
                return False

            # 步骤3: GenomeScope
            # histo_file 已在上面定义|histo_file is defined above
            gs_output_dir = os.path.join(output_dir_abs, "02_genomescope")
            os.makedirs(gs_output_dir, exist_ok=True)

            kcov = self.gs_runner.run_genomescope(
                histo_file, self.config.kmer_size, self.config.read_length,
                gs_output_dir, self.config.max_kmer_cov, self.config.ploidy
            )

            if kcov is None:
                self.logger.warning("未能获取kcov值|Failed to get kcov value")
                return False

            self.config.kcov = kcov

            # 完成统计|Completion statistics
            elapsed_time = time.time() - start_time

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间: {elapsed_time:.2f} seconds|Total runtime: {elapsed_time:.2f} seconds")
            self.logger.info(f"K-mer coverage: {kcov}")
            self.logger.info(f"GenomeScope输出: {gs_output_dir}|GenomeScope output: {gs_output_dir}")
            self.logger.info("Pipeline completed successfully")

            return True

        except Exception as e:
            self.logger.error(f"Pipeline失败|Pipeline failed: {str(e)}", exc_info=True)
            return False


if __name__ == "__main__":
    main()
