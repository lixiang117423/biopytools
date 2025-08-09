"""
命令行接口模块 | Command Line Interface Module 💻⌨️
"""

import argparse
import sys
import os
import json
from pathlib import Path
from typing import List, Optional

from .config import KmerConfig, AssignmentStrategy
from .analyzer import KmerAnalyzer

def create_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器 🛠️"""
    
    parser = argparse.ArgumentParser(
        prog='run_kmer_analysis',
        description="""
Universal K-mer Analysis Toolkit | 通用K-mer分析工具包 🧬🔬

Supports flexible k-mer analysis for FASTA/FASTQ files with automatic role detection
支持FASTA/FASTQ文件的灵活k-mer分析，具备自动角色检测功能

Examples | 使用示例:
  # Auto-assignment mode | 自动分配模式 ✨
  run_kmer_analysis --input /data/genomes/*.fa /data/samples/*.fq.gz \\
                 --output results/ --kmer-size 51
  
  # Explicit assignment | 明确指定角色 👉
  run_kmer_analysis --kmer-sources /ref/*.fa \\
                 --query-targets /samples/*.fastq.gz \\
                 --output results/
  
  # Interactive mode | 交互模式 💬
  run_kmer_analysis --input /data/mixed/ --interactive
  
  # Configuration file | 配置文件模式 📄
  run_kmer_analysis --config analysis.yaml
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # === 输入输出参数 | Input/Output Parameters ===
    io_group = parser.add_argument_group('Input/Output | 输入输出 📂')
    
    io_group.add_argument(
        '--input', '-i',
        nargs='+',
        help='Input file paths (files, directories, wildcards) | 输入文件路径（文件、目录、通配符） 📁'
    )
    
    io_group.add_argument(
        '--kmer-sources', '--ks',
        nargs='+',
        help='K-mer source files (explicit mode) | K-mer库来源文件（明确模式） 📚'
    )
    
    io_group.add_argument(
        '--query-targets', '--qt',
        nargs='+',
        help='Query target files (explicit mode) | 查询目标文件（明确模式） 🎯'
    )
    
    io_group.add_argument(
        '--output', '-o',
        default='kmer_analysis_results',
        help='Output directory (default: kmer_analysis_results) | 输出目录（默认：kmer_analysis_results） 📂'
    )
    
    io_group.add_argument(
        '--config', '-c',
        help='Configuration file (YAML format) | 配置文件（YAML格式） 📄'
    )

    io_group.add_argument(
    '--fastq-pattern', '--fp',
    help='FASTQ file naming pattern, e.g., "*_1.clean.fq.gz" | FASTQ文件命名模式，如"*_1.clean.fq.gz"'
    )

    io_group.add_argument(
        '--fasta-pattern', '--fap', 
        help='FASTA file naming pattern, e.g., "*.fa" | FASTA文件命名模式，如"*.fa"'
    )
        
    # === K-mer参数 | K-mer Parameters ===
    kmer_group = parser.add_argument_group('K-mer Parameters | K-mer参数 🧬')
    
    kmer_group.add_argument(
        '--kmer-size', '-k',
        type=int,
        default=51,
        help='K-mer size (default: 51) | K-mer长度（默认：51） 📏'
    )
    
    kmer_group.add_argument(
        '--min-count', '--ci',
        type=int,
        default=1,
        help='Minimum k-mer count (default: 1) | 最小k-mer计数（默认：1） 📉'
    )
    
    kmer_group.add_argument(
        '--max-count', '--cx',
        type=int,
        default=int(1e9),
        help='Maximum k-mer count (default: 1e9) | 最大k-mer计数（默认：1e9） 📈'
    )
    
    kmer_group.add_argument(
        '--no-canonical', '--nc',
        action='store_true',
        help='Disable canonical k-mer form | 禁用k-mer标准形式 🔄'
    )
    
    # === 系统资源参数 | System Resource Parameters ===
    sys_group = parser.add_argument_group('System Resources | 系统资源 🖥️')
    
    sys_group.add_argument(
        '--threads', '-t',
        type=int,
        default=88,
        help='Number of threads (default: 88) | 线程数（默认：88） 🧵'
    )
    
    sys_group.add_argument(
        '--memory', '-m',
        type=int,
        default=880,
        help='Memory limit in GB (default: 880) | 内存限制（GB，默认：880） 🧠'
    )
    
    sys_group.add_argument(
        '--temp-dir', '--tmp',
        default='/tmp/kmer_universal',
        help='Temporary directory (default: /tmp/kmer_universal) | 临时目录（默认：/tmp/kmer_universal） 🗑️'
    )
    
    # === 分析参数 | Analysis Parameters ===
    analysis_group = parser.add_argument_group('Analysis Parameters | 分析参数 🔬')
    
    analysis_group.add_argument(
        '--assignment-strategy', '--as',
        choices=['explicit', 'size_based', 'type_based', 'intelligent', 'interactive'],
        default='intelligent',
        help='Role assignment strategy (default: intelligent) | 角色分配策略（默认：intelligent） 🤔'
    )
    
    analysis_group.add_argument(
        '--window-sizes', '--ws',
        nargs='+',
        type=int,
        default=[500000],
        help='Sliding window sizes (default: 500000) | 滑窗大小（默认：500000） 🪟'
    )
    
    analysis_group.add_argument(
        '--size-threshold', '--st',
        type=float,
        default=1.0,
        help='File size threshold in GB for role assignment (default: 1.0) | 角色分配的文件大小阈值（GB，默认：1.0） ⚖️'
    )
    
    analysis_group.add_argument(
        '--no-positions', '--np',
        action='store_true',
        help='Disable position tracking for FASTA files | 禁用FASTA文件的位置追踪 📍'
    )
    
    # === 输出格式参数 | Output Format Parameters ===
    output_group = parser.add_argument_group('Output Formats | 输出格式 📝')
    
    output_group.add_argument(
        '--output-formats', '--of',
        nargs='+',
        choices=['fasta', 'csv', 'txt'],
        default=['fasta', 'csv', 'txt'],
        help='Output formats (default: fasta csv txt) | 输出格式（默认：fasta csv txt） 📄'
    )
    
    output_group.add_argument(
        '--no-library', '--nl',
        action='store_true',
        help='Skip k-mer library output | 跳过k-mer库输出 📚'
    )
    
    output_group.add_argument(
        '--no-sliding-window', '--nsw',
        action='store_true',
        help='Skip sliding window analysis | 跳过滑窗分析 🪟'
    )
    
    # === 高级选项 | Advanced Options ===
    advanced_group = parser.add_argument_group('Advanced Options | 高级选项 🚀')
    
    advanced_group.add_argument(
        '--chunk-size', '--cs',
        type=float,
        default=2.0,
        help='File chunk size in GB (default: 2.0) | 文件分片大小（GB，默认：2.0） 🍰'
    )
    
    advanced_group.add_argument(
        '--cache-size', '--cache',
        type=float,
        default=50.0,
        help='Cache size in GB (default: 50.0) | 缓存大小（GB，默认：50.0） 💨'
    )
    
    advanced_group.add_argument(
        '--signature-length', '--sl',
        type=int,
        choices=range(5, 12),
        default=9,
        help='KMC signature length (default: 9) | KMC签名长度（默认：9） ✍️'
    )
    
    advanced_group.add_argument(
        '--ram-only', '--ro',
        action='store_true',
        help='Use RAM-only mode for KMC | KMC使用仅RAM模式 ⚡️'
    )
    
    advanced_group.add_argument(
        '--no-strict-memory', '--nsm',
        action='store_true',
        help='Disable strict memory mode | 禁用严格内存模式 🔓'
    )
    
    advanced_group.add_argument(
        '--enable-compression', '--ec',
        action='store_true',
        default=True,
        help='Enable compression (default: True) | 启用压缩（默认：True） 📦'
    )
    
    advanced_group.add_argument(
        '--keep-intermediate', '--ki',
        action='store_true',
        help='Keep intermediate files | 保留中间文件 🖇️'
    )
    
    advanced_group.add_argument(
        '--resume', '-r',
        action='store_true',
        help='Resume from previous run | 从之前的运行恢复 ⏯️'
    )
    
    # === 运行模式 | Run Modes ===
    mode_group = parser.add_argument_group('Run Modes | 运行模式 🏃')
    
    mode_group.add_argument(
        '--interactive', '--int',
        action='store_true',
        help='Interactive mode for file role assignment | 交互模式进行文件角色分配 💬'
    )
    
    mode_group.add_argument(
        '--dry-run', '--dr',
        action='store_true',
        help='Dry run mode (show what would be done) | 干运行模式（显示将要执行的操作） 💨'
    )
    
    mode_group.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output | 详细输出 🗣️'
    )
    
    mode_group.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Quiet mode (minimal output) | 静默模式（最少输出） 🤫'
    )
    
    # === 其他选项 | Other Options ===
    parser.add_argument(
        '--version',
        action='version',
        version='K-mer Universal Analysis Toolkit v1.0.0 ℹ️'
    )
    
    parser.add_argument(
        '--save-config', '--sc',
        help='Save current configuration to file | 保存当前配置到文件 💾'
    )
    
    parser.add_argument(
        '--export-results', '--er',
        help='Export results in JSON format | 以JSON格式导出结果 📤'
    )
    
    return parser

def validate_args(args) -> None:
    """验证命令行参数 ✅"""
    
    # 检查输入参数
    if not args.config:
        if not (args.input or (args.kmer_sources and args.query_targets)):
            raise ValueError(
                "❌ Must provide either --input for auto mode or "
                "--kmer-sources and --query-targets for explicit mode"
            )
        
        if args.input and (args.kmer_sources or args.query_targets):
            raise ValueError(
                "❌ Cannot use --input with --kmer-sources/--query-targets simultaneously"
            )
    
    # 检查输出目录
    output_dir = Path(args.output)
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"❌ Output path exists but is not a directory: {args.output}")
    
    # 检查系统资源
    if args.threads < 1:
        raise ValueError("❌ Number of threads must be >= 1")
    
    if args.memory < 1:
        raise ValueError("❌ Memory limit must be >= 1 GB")
    
    # 检查k-mer参数
    if not (1 <= args.kmer_size <= 256):
        raise ValueError("❌ K-mer size must be between 1 and 256")
    
    if args.min_count < 0:
        raise ValueError("❌ Minimum count must be >= 0")
    
    if args.max_count <= args.min_count:
        raise ValueError("❌ Maximum count must be > minimum count")
    
    # 检查窗口大小
    if any(ws <= 0 for ws in args.window_sizes):
        raise ValueError("❌ Window sizes must be positive")

def create_config_from_args(args) -> KmerConfig:
    """从命令行参数创建配置对象 ⚙️"""
    
    config_dict = {
        'kmer_size': args.kmer_size,
        'threads': args.threads,
        'memory_gb': args.memory,
        'min_count': args.min_count,
        'max_count': args.max_count,
        'canonical_form': not args.no_canonical,
        'signature_length': args.signature_length,
        'ram_only_mode': args.ram_only,
        'strict_memory': not args.no_strict_memory,
        'output_dir': args.output,
        'temp_dir': args.temp_dir,
        'fastq_pattern': args.fastq_pattern,
        'fasta_pattern': args.fasta_pattern,
        'assignment_strategy': AssignmentStrategy(args.assignment_strategy),
        'size_threshold_gb': args.size_threshold,
        'window_sizes': [] if args.no_sliding_window else args.window_sizes,
        'include_positions': not args.no_positions,
        'output_formats': [] if args.no_library else args.output_formats,
        'chunk_size_gb': args.chunk_size,
        'enable_compression': args.enable_compression,
        'cache_size_gb': args.cache_size,
        'keep_intermediate': args.keep_intermediate,
        'verbose': args.verbose and not args.quiet,
        'dry_run': args.dry_run,
        'resume': args.resume
    }
    
    # 设置输入路径
    if args.input:
        config_dict['input_paths'] = args.input
    if args.kmer_sources:
        config_dict['kmer_source_paths'] = args.kmer_sources
    if args.query_targets:
        config_dict['query_target_paths'] = args.query_targets
    
    return KmerConfig(**config_dict)

def main():
    """主函数 🚀"""
    
    try:
        # 解析命令行参数
        parser = create_parser()
        args = parser.parse_args()
        
        # 从配置文件加载（如果提供）
        if args.config:
            if not os.path.exists(args.config):
                print(f"❌ Error: Configuration file not found: {args.config}", file=sys.stderr)
                sys.exit(1)
            
            config = KmerConfig.from_yaml(args.config)
            # 命令行参数覆盖配置文件
            if args.output != 'kmer_analysis_results':  # 默认值检查
                config.output_dir = args.output
            if args.verbose:
                config.verbose = True
            if args.quiet:
                config.verbose = False
        else:
            # 验证参数
            validate_args(args)
            # 从命令行参数创建配置
            config = create_config_from_args(args)
        
        # 保存配置（如果需要）
        if args.save_config:
            config.to_yaml(args.save_config)
            print(f"💾 Configuration saved to: {args.save_config}")
            if args.dry_run:
                return
        
        # 干运行模式
        if args.dry_run:
            print("💨 === Dry Run Mode ===")
            print(f"K-mer size: {config.kmer_size}")
            print(f"Threads: {config.threads}")
            print(f"Memory: {config.memory_gb}GB")
            print(f"Output directory: {config.output_dir}")
            print(f"Assignment strategy: {config.assignment_strategy.value}")
            if config.input_paths:
                print(f"Input paths: {config.input_paths}")
            if config.kmer_source_paths:
                print(f"K-mer sources: {config.kmer_source_paths}")
            if config.query_target_paths:
                print(f"Query targets: {config.query_target_paths}")
            print("Would execute analysis with above configuration")
            return
        
        # 创建分析器
        analyzer = KmerAnalyzer(config)
        
        # 执行分析
        if config.input_paths:
            # 自动模式
            if args.interactive:
                config.assignment_strategy = AssignmentStrategy.INTERACTIVE
            
            results = analyzer.auto_analyze(
                input_paths=config.input_paths,
                output_dir=config.output_dir,
                kmer_size=config.kmer_size,
                assignment_strategy=config.assignment_strategy
            )
        else:
            # 明确模式
            results = analyzer.explicit_analyze(
                kmer_sources=config.kmer_source_paths,
                query_targets=config.query_target_paths,
                output_dir=config.output_dir,
                kmer_size=config.kmer_size
            )
        
        # 显示结果
        if not args.quiet:
            print("\n🎉 === Analysis Completed Successfully ===")
            print(f"📚 K-mer library size: {results['kmer_library_size']:,}")
            print(f"🎯 Target samples: {results['target_samples']}")
            print(f"⏱️ Runtime: {results['runtime_seconds']:.1f} seconds")
            print(f"📂 Output directory: {config.output_dir}")
            
            if results.get('output_files'):
                print("\n📄 Output files:")
                for file_type, file_path in results['output_files'].items():
                    print(f"  {file_type}: {file_path}")
        
        # 导出结果（如果需要）
        if args.export_results:
            # 移除不可序列化的对象
            export_results = {k: v for k, v in results.items() if k != 'config'}
            export_results['config_summary'] = {
                'kmer_size': config.kmer_size,
                'threads': config.threads,
                'memory_gb': config.memory_gb,
                'assignment_strategy': config.assignment_strategy.value
            }
            
            with open(args.export_results, 'w') as f:
                json.dump(export_results, f, indent=2)
            print(f"📤 Results exported to: {args.export_results}")
        
    except KeyboardInterrupt:
        print("\n🛑 Analysis interrupted by user", file=sys.stderr)
        sys.exit(1)
    
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()