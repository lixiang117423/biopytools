"""泛基因组Block构建工具入口|Pan-Genome Block Construction Tool Entry Point"""

import argparse
import re
import sys
import os
from pathlib import Path
from datetime import datetime

FASTA_EXTENSIONS = {'.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'}


def _natural_sort_key(s: str):
    """自然排序key：字母排序，数字按数值排序|Natural sort key: alphabetical with numeric values"""
    parts = re.split(r'(\d+)', s)
    return [(not p.isdigit(), p.lower()) if not p.isdigit() else (False, int(p)) for p in parts]


def _extract_genome_name(filename: str) -> str:
    """从FASTA文件名提取基因组名称|Extract genome name from FASTA filename"""
    name = filename
    for ext in sorted(FASTA_EXTENSIONS, key=len, reverse=True):
        if name.lower().endswith(ext):
            name = name[:-len(ext)]
            break
    return name


def prepare_genome_list(input_dir: str, output_file: str):
    """扫描目录中的FASTA文件，生成genome_list.txt|Scan FASTA files in directory and generate genome_list.txt"""
    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        print(f"错误|Error: 目录不存在|Directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    fasta_files = []
    for fname in os.listdir(input_dir):
        lower = fname.lower()
        if any(lower.endswith(ext) for ext in FASTA_EXTENSIONS):
            full_path = os.path.join(input_dir, fname)
            if os.path.isfile(full_path):
                fasta_files.append(fname)

    if not fasta_files:
        print(f"错误|Error: 目录中未找到FASTA文件|No FASTA files found in: {input_dir}", file=sys.stderr)
        sys.exit(1)

    fasta_files.sort(key=_natural_sort_key)

    output_file = os.path.abspath(output_file)
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)

    with open(output_file, 'w') as f:
        for fname in fasta_files:
            name = _extract_genome_name(fname)
            full_path = os.path.join(input_dir, fname)
            f.write(f"{name}\t{full_path}\n")

    print(f"基因组列表已生成|Genome list generated: {output_file}")
    print(f"共|Total: {len(fasta_files)} 个基因组|genomes")


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="泛基因组Block构建工具|Pan-Genome Block Construction Tool"
    )

    parser.add_argument('-i', '--genome-list',
                        help='基因组列表文件|Genome list file (name<TAB>path)')
    parser.add_argument('-o', '--output-dir', default='./pan_blocks_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads (default: 12)')
    parser.add_argument('--parallel-alignments', type=int, default=4,
                        help='并行比对数|Number of parallel alignments (default: 4)')
    parser.add_argument('--min-alignment-length', type=int, default=10000,
                        help='最小比对长度|Minimum alignment length for delta-filter (default: 10000)')
    parser.add_argument('--genome-order',
                        help='基因组优先级顺序文件|Genome priority order file (one name per line)')
    parser.add_argument('--chromosome',
                        help='指定染色体|Specific chromosome to process')
    parser.add_argument('--step', choices=['align', 'build', 'plot'],
                        help='执行特定步骤|Run specific step')
    parser.add_argument('--nucmer', help='nucmer可执行文件路径|nucmer executable path')
    parser.add_argument('--delta-filter', help='delta-filter可执行文件路径|delta-filter executable path')
    parser.add_argument('--show-coords', help='show-coords可执行文件路径|show-coords executable path')
    parser.add_argument('--bedtools', help='bedtools可执行文件路径|bedtools executable path')
    parser.add_argument('--minimap2', help='minimap2可执行文件路径|minimap2 executable path')
    parser.add_argument('--plot-format', choices=['svg', 'png'], default='svg',
                        help='绘图格式|Plot format (default: svg)')
    parser.add_argument('--plot-width', type=int, default=20,
                        help='绘图宽度|Plot width (default: 20)')
    parser.add_argument('--plot-height', type=int, default=10,
                        help='绘图高度|Plot height (default: 10)')
    parser.add_argument('--prepare', metavar='DIR',
                        help='从目录自动生成genome_list.txt|Auto-generate genome_list.txt from FASTA directory')
    parser.add_argument('--prepare-output', default='./genome_list.txt',
                        help='prepare输出文件路径|Output file path for prepare (default: ./genome_list.txt)')

    return parser.parse_args()


def generate_software_versions_yml(output_dir: str, config,
                                   start_time: datetime):
    """生成software_versions.yml|Generate software_versions.yml"""
    import subprocess
    from pathlib import Path

    end_time = datetime.now()
    runtime = int((end_time - start_time).total_seconds())

    tools = {}
    for name, path in [('nucmer', config.nucmer_path), ('delta-filter', config.delta_filter_path),
                       ('show-coords', config.show_coords_path), ('bedtools', config.bedtools_path),
                       ('minimap2', config.minimap2_path)]:
        try:
            result = subprocess.run([path, '--version'], capture_output=True, text=True, timeout=10)
            version = result.stdout.strip().split('\n')[0] if result.returncode == 0 else 'unknown'
            tools[name] = {'version': version, 'path': path}
        except Exception:
            tools[name] = {'version': 'unknown', 'path': path}

    info = {
        'pipeline': {'name': 'biopytools pan_blocks', 'version': '1.0.0'},
        'tools': tools,
        'parameters': {
            'threads': config.threads,
            'parallel_alignments': config.parallel_alignments,
            'min_alignment_length': config.min_alignment_length,
            'genome_count': len(config.genomes),
            'genome_order': config.genome_order_list,
        },
        'execution': {
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'runtime_seconds': runtime,
        }
    }

    output_file = Path(output_dir) / '00_pipeline_info' / 'software_versions.yml'
    output_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        import yaml
        with open(output_file, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
    except ImportError:
        import json
        with open(output_file.with_suffix('.json'), 'w') as f:
            json.dump(info, f, indent=2)


def main():
    """主函数|Main function"""
    args = parse_arguments()

    if args.prepare:
        prepare_genome_list(args.prepare, args.prepare_output)
        sys.exit(0)

    if not args.genome_list:
        print("错误|Error: 缺少必需参数 -i/--genome-list|Missing required argument -i/--genome-list", file=sys.stderr)
        sys.exit(1)

    start_time = datetime.now()

    try:
        from .config import PanBlocksConfig
        from .utils import PanBlocksLogger, check_dependencies

        config_kwargs = {
            'genome_list': args.genome_list,
            'output_dir': args.output_dir,
            'threads': args.threads,
            'parallel_alignments': args.parallel_alignments,
            'min_alignment_length': args.min_alignment_length,
            'plot_format': args.plot_format,
            'plot_width': args.plot_width,
            'plot_height': args.plot_height,
        }

        if args.genome_order:
            config_kwargs['genome_order_file'] = args.genome_order
        if args.chromosome:
            config_kwargs['chromosome'] = args.chromosome
        if args.nucmer:
            config_kwargs['nucmer_path'] = args.nucmer
        if args.delta_filter:
            config_kwargs['delta_filter_path'] = args.delta_filter
        if args.show_coords:
            config_kwargs['show_coords_path'] = args.show_coords
        if args.bedtools:
            config_kwargs['bedtools_path'] = args.bedtools
        if args.minimap2:
            config_kwargs['minimap2_path'] = args.minimap2

        config = PanBlocksConfig(**config_kwargs)
        config.validate()

        logger_manager = PanBlocksLogger(config.output_dir)
        logger = logger_manager.get_logger()

        logger.info("泛基因组Block构建工具启动|Pan-Blocks Construction Tool started")
        logger.info(f"基因组数量|Genome count: {len(config.genomes)}")
        logger.info(f"基因组顺序|Genome order: {', '.join(config.genome_order_list)}")

        check_dependencies(config, logger)

        # Step 1: 两两比对
        if args.step in [None, 'align']:
            from .aligner import PairwiseAligner
            aligner = PairwiseAligner(config, logger)
            if not aligner.run_all_alignments():
                logger.error("两两比对失败|Pairwise alignment failed")
                sys.exit(1)

        # Step 2: 构建Pan-Blocks
        if args.step in [None, 'build']:
            from .builder import PanBlockBuilder
            builder = PanBlockBuilder(config, logger)
            if not builder.build_all_chromosomes():
                logger.error("Pan-Blocks构建失败|Pan-blocks construction failed")
                sys.exit(1)

        # Step 3: 可视化
        if args.step in [None, 'plot']:
            from .plotter import PanBlocksPlotter
            plotter = PanBlocksPlotter(config, logger)
            if not plotter.plot_all_chromosomes():
                logger.error("可视化失败|Plotting failed")
                sys.exit(1)

        generate_software_versions_yml(config.output_dir, config, start_time)

        logger.info(f"泛基因组Block构建完成|Pan-Blocks construction completed")
        logger.info(f"输出目录|Output directory: {config.output_dir}")
        sys.exit(0)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
