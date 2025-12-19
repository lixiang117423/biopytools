"""
TASSEL GWAS分析主程序 | TASSEL GWAS Analysis Main Program
"""

import argparse
import sys
from pathlib import Path

# 确保可以导入模块 | Ensure module can be imported
sys.path.insert(0, str(Path(__file__).parent.parent))

from .core import TASSELGWASAnalyzer
from .batch_processor import TASSELBatchProcessor
from .utils import TASSELLogger


def parse_arguments():
    """解析命令行参数 | Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="TASSEL GWAS分析工具 🌾 | TASSEL GWAS Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:

  # 基本GWAS分析 - 自动识别并处理所有表型
  python -m biopytools.tassel_gwas -i input.vcf.gz -p traits.txt -o results

  # 指定模型和参数
  python -m biopytools.tassel_gwas -i input.vcf.gz -p traits.txt -o results \\
      --model MLM --memory 200g --maf 0.05 --miss 0.1

  # 使用Q矩阵和Kinship矩阵
  python -m biopytools.tassel_gwas -i input.vcf.gz -p traits.txt -o results \\
      --q-matrix Q.txt --kinship kinship.txt

  # 并行处理多个表型
  python -m biopytools.tassel_gwas -i input.vcf.gz -p traits.txt -o results \\
      --parallel --workers 8
        """
    )

    # 必需参数 | Required parameters
    parser.add_argument('-i', '--vcf', required=True,
                        help='VCF基因型文件路径 | VCF genotype file path')
    parser.add_argument('-p', '--pheno', required=True,
                        help='表型文件路径 | Phenotype file path (可包含多个表型)')
    parser.add_argument('-o', '--output', required=True,
                        help='输出目录路径 | Output directory path')

    # 分析参数 | Analysis parameters
    parser.add_argument('--model', choices=['GLM', 'MLM', 'BOTH'], default='MLM',
                        help='GWAS模型选择 | GWAS model selection (default: MLM)')
    parser.add_argument('--memory', default='100g',
                        help='Java最大内存 | Java maximum memory (default: 100g)')
    parser.add_argument('--threads', type=int, default=4,
                        help='并行线程数 | Number of parallel threads (default: 4)')
    parser.add_argument('--maf', type=float,
                        help='最小等位基因频率过滤 | Minimum allele frequency filter')
    parser.add_argument('--miss', type=float,
                        help='最大缺失率过滤 | Maximum missing rate filter')
    parser.add_argument('--pca-components', type=int, default=5,
                        help='PCA主成分数量（用作MLM协变量）| Number of PCA components for MLM covariates (default: 5)')

    # 矩阵文件 | Matrix files
    parser.add_argument('--q-matrix',
                        help='群体结构Q矩阵文件 | Population structure Q matrix file')
    parser.add_argument('--kinship',
                        help='亲缘关系K矩阵文件 | Kinship K matrix file')

    # 其他选项 | Other options
    parser.add_argument('--tassel-path',
                        help='TASSEL安装路径 | TASSEL installation path')
    parser.add_argument('--skip-sort', action='store_true',
                        help='跳过VCF排序 | Skip VCF sorting')
    parser.add_argument('--keep-temp', action='store_true',
                        help='保留临时文件 | Keep temporary files')

    # 并行处理选项 | Parallel processing options
    parser.add_argument('--parallel', action='store_true',
                        help='并行处理多个表型 | Parallel process multiple traits')
    parser.add_argument('--workers', type=int, default=4,
                        help='并行工作线程数 | Number of parallel workers (default: 4)')

    # 日志选项 | Logging options
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARN', 'ERROR'],
                        default='INFO', help='日志级别 | Log level (default: INFO)')

    return parser.parse_args()


def main():
    """主函数 | Main function"""
    try:
        # 解析参数 | Parse arguments
        args = parse_arguments()

        # 验证输入文件 | Validate input files
        vcf_file = Path(args.vcf).expanduser()
        pheno_file = Path(args.pheno).expanduser()
        output_dir = Path(args.output).expanduser()

        if not vcf_file.exists():
            print(f"❌ VCF文件不存在 | VCF file does not exist: {vcf_file}")
            sys.exit(1)

        if not pheno_file.exists():
            print(f"❌ 表型文件不存在 | Phenotype file does not exist: {pheno_file}")
            sys.exit(1)

        # 创建输出目录 | Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # 自动识别表型数量并处理 | Automatically detect and process all traits
        print(f"🌾 开始GWAS分析 | Starting GWAS analysis")
        print(f"📊 VCF文件 | VCF file: {vcf_file}")
        print(f"📋 表型文件 | Phenotype file: {pheno_file}")
        print(f"📁 输出目录 | Output directory: {output_dir}")

        # 设置日志 | Setup logging
        logger_manager = TASSELLogger(output_dir, "gwas.log")
        logger = logger_manager.get_logger()

        # 配置批量处理器 | Configure batch processor
        global_config = {
            'model': args.model,
            'memory_max': args.memory,
            'threads': args.threads,
            'maf_filter': args.maf,
            'miss_filter': args.miss,
            'skip_sort': args.skip_sort,
            'pca_components': args.pca_components,  # 添加PCA数量参数
            'q_matrix': args.q_matrix,
            'kinship': args.kinship,
            'keep_temp': args.keep_temp,
            'log_level': args.log_level
        }

        if args.tassel_path:
            global_config['tassel_path'] = Path(args.tassel_path).expanduser()

        # 创建批量处理器 | Create batch processor
        processor = TASSELBatchProcessor(global_config, logger)

        # 运行GWAS分析 | Run GWAS analysis
        success = processor.process_batch(
            pheno_file, vcf_file, output_dir, global_config,
            parallel=args.parallel, max_workers=args.workers
        )

        # 输出摘要 | Print summary
        summary = processor.get_summary()
        print(f"\n📊 分析完成 | Analysis completed")
        print(f"   总表型数 | Total traits: {summary['total']}")
        print(f"   成功处理 | Successful: {summary['successful']}")
        print(f"   处理失败 | Failed: {summary['failed']}")
        print(f"   成功率 | Success rate: {summary['success_rate']:.1f}%")

        if summary['failed_traits']:
            print(f"\n❌ 失败的表型 | Failed traits:")
            for trait in summary['failed_traits']:
                print(f"   • {trait['trait_name']}: {trait['error']}")

        sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        print("\n🛑 用户中断操作 | User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 程序错误 | Program error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()