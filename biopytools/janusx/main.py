"""
JanusX主程序模块|JanusX Main Module

命令行接口入口|Command-line interface entry point
"""

import argparse
import sys
from typing import List

from .config import JanusXGWASConfig, JanusXGSConfig, JanusXPCAConfig, JanusXPostGWASConfig
from .gwas_runner import JanusXGWASRunner
from .gs_runner import JanusXGSRunner
from .pca_runner import JanusXPCARunner
from .postgwas_runner import JanusXPostGWASRunner


def parse_gwas_arguments(subparser):
    """解析GWAS参数|Parse GWAS arguments"""
    parser = subparser.add_parser(
        'gwas',
        help='全基因组关联分析|Genome-Wide Association Study (GWAS)'
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--genotype', required=True,
                       help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
    parser.add_argument('-p', '--pheno', required=True,
                       help='表型文件|Phenotype file')

    # 可选参数|Optional parameters
    parser.add_argument('-t', '--type', default='vcf', choices=['vcf', 'bfile'],
                       help='基因型类型|Genotype type (default: vcf)')
    parser.add_argument('-m', '--models', nargs='+', default=['lmm'],
                       choices=['lm', 'lmm', 'fastlmm', 'farmcpu'],
                       help='GWAS模型|GWAS models (default: lmm)')
    parser.add_argument('--maf', type=float, default=0.05,
                       help='最小等位基因频率阈值|Minor allele frequency threshold (default: 0.05)')
    parser.add_argument('--geno', type=float, default=0.0,
                       help='缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)')
    parser.add_argument('-k', '--grm', default='1',
                       help='亲缘关系矩阵方法|GRM method (default: 1)')
    parser.add_argument('-q', '--qcov', type=int, default=0,
                       help='PCA数量或Q矩阵文件路径|Number of PCs or path to Q matrix file (default: 0, no Q matrix)')
    parser.add_argument('-c', '--cov',
                       help='协变量文件|Covariate file')
    parser.add_argument('-n', '--ncol', nargs='+', type=int,
                       help='表型列索引(零基)|Phenotype column indices (zero-based)')
    parser.add_argument('--plot', action='store_true',
                       help='生成图表|Generate plots')
    parser.add_argument('--chunksize', type=int, default=100000,
                       help='SNP分块大小|SNP chunk size (default: 100000)')
    parser.add_argument('--mmap-limit', type=int,
                       help='内存映射限制|Memory map limit')
    parser.add_argument('-th', '--threads', type=int, default=12,
                       help='线程数|Number of threads (default: 12)')
    parser.add_argument('-o', '--output-dir', default='./janusx_gwas_output',
                       help='输出目录|Output directory (default: ./janusx_gwas_output)')
    parser.add_argument('--prefix',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--janusx-path',
                       help='JanusX可执行文件路径|JanusX executable path')


def parse_gs_arguments(subparser):
    """解析GS参数|Parse GS arguments"""
    parser = subparser.add_parser(
        'gs',
        help='基因组选择分析|Genomic Selection (GS)'
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--genotype', required=True,
                       help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
    parser.add_argument('-p', '--pheno', required=True,
                       help='表型文件|Phenotype file')

    # 可选参数|Optional parameters
    parser.add_argument('-t', '--type', default='vcf', choices=['vcf', 'bfile'],
                       help='基因型类型|Genotype type (default: vcf)')
    parser.add_argument('-m', '--models', nargs='+', default=['GBLUP'],
                       choices=['GBLUP', 'rrBLUP', 'BayesA', 'BayesB', 'BayesCpi'],
                       help='GS模型|GS models (default: GBLUP)')
    parser.add_argument('--maf', type=float, default=0.05,
                       help='最小等位基因频率阈值|Minor allele frequency threshold (default: 0.05)')
    parser.add_argument('--geno', type=float, default=0.0,
                       help='缺失率阈值(0=不过滤|Missing rate threshold, 0=no filtering)')
    parser.add_argument('--pcd', action='store_true',
                       help='启用PCA降维|Enable PCA-based dimensionality reduction')
    parser.add_argument('-n', '--ncol', type=int,
                       help='表型列索引(零基)|Phenotype column index (zero-based)')
    parser.add_argument('--cv', type=int,
                       help='交叉验证折数|K-fold cross-validation')
    parser.add_argument('--plot', action='store_true',
                       help='生成图表|Generate plots')
    parser.add_argument('-o', '--output-dir', default='./janusx_gs_output',
                       help='输出目录|Output directory (default: ./janusx_gs_output)')
    parser.add_argument('--prefix',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--janusx-path',
                       help='JanusX可执行文件路径|JanusX executable path')


def parse_pca_arguments(subparser):
    """解析PCA参数|Parse PCA arguments"""
    parser = subparser.add_parser(
        'pca',
        help='主成分分析|Principal Component Analysis (PCA)'
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--genotype',
                       help='基因型文件(VCF或PLINK前缀)|Genotype file (VCF or PLINK prefix)')
    parser.add_argument('-t', '--type', default='vcf', choices=['vcf', 'bfile'],
                       help='基因型类型|Genotype type (default: vcf)')
    parser.add_argument('--grm',
                       help='预计算的GRM前缀|Precomputed GRM prefix')
    parser.add_argument('--pcfile',
                       help='预计算的PCA文件|Precomputed PCA file')

    # 可选参数|Optional parameters
    parser.add_argument('-d', '--dim', type=int, default=3,
                       help='输出主成分数量|Number of PCs to output (default: 3)')
    parser.add_argument('--plot', action='store_true',
                       help='生成2D散点图|Generate 2D scatter plot')
    parser.add_argument('--plot3d', action='store_true',
                       help='生成3D旋转GIF|Generate 3D rotating GIF')
    parser.add_argument('-g', '--group',
                       help='分组文件路径|Group file path')
    parser.add_argument('--color', type=int, default=1,
                       help='调色板索引(0-6)|Color palette index (default: 1)')
    parser.add_argument('-o', '--output-dir', default='./janusx_pca_output',
                       help='输出目录|Output directory (default: ./janusx_pca_output)')
    parser.add_argument('--prefix',
                       help='输出文件前缀|Output file prefix')
    parser.add_argument('--janusx-path',
                       help='JanusX可执行文件路径|JanusX executable path')


def parse_postgwas_arguments(subparser):
    """解析PostGWAS参数|Parse PostGWAS arguments"""
    parser = subparser.add_parser(
        'postgwas',
        help='GWAS结果可视化和注释|GWAS result visualization and annotation'
    )

    # 必需参数|Required parameters
    parser.add_argument('-f', '--files', nargs='+', required=True,
                       help='GWAS结果文件列表|List of GWAS result files')

    # 可选参数|Optional parameters
    parser.add_argument('--chr-col', default='#CHROM',
                       help='染色体列名|Chromosome column name (default: #CHROM)')
    parser.add_argument('--pos-col', default='POS',
                       help='位置列名|Position column name (default: POS)')
    parser.add_argument('--pvalue-col', default='p',
                       help='P值列名|P-value column name (default: p)')
    parser.add_argument('--threshold', type=float,
                       help='显著性阈值|Significance threshold')
    parser.add_argument('--noplot', action='store_true',
                       help='禁用绘图|Disable plotting')
    parser.add_argument('--color', type=int, default=0,
                       help='颜色方案索引(0-6, -1为auto)|Color scheme index (default: 0)')
    parser.add_argument('--hl', '--highlight',
                       help='高亮区域BED文件|Highlight regions BED file')
    parser.add_argument('--format', default='png', choices=['pdf', 'png', 'svg', 'tif'],
                       help='输出格式|Output format (default: png)')
    parser.add_argument('-a', '--anno',
                       help='注释文件路径|Annotation file path')
    parser.add_argument('--ab', '--anno-broaden', type=int,
                       help='注释窗口|Annotation window (kb)')
    parser.add_argument('--desc-item', default='description',
                       help='GFF描述键|GFF description key (default: description)')
    parser.add_argument('-o', '--output-dir', default='./janusx_postgwas_output',
                       help='输出目录|Output directory (default: ./janusx_postgwas_output)')
    parser.add_argument('--prefix', default='JanusX',
                       help='输出文件前缀|Output file prefix (default: JanusX)')
    parser.add_argument('-th', '--threads', type=int, default=12,
                       help='线程数|Number of threads (default: 12)')
    parser.add_argument('--janusx-path',
                       help='JanusX可执行文件路径|JanusX executable path')


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='JanusX GWAS和基因组选择分析工具|JanusX GWAS and Genomic Selection Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # GWAS分析|GWAS analysis
  janusx gwas -i data.vcf.gz -p pheno.txt -m lmm -o output

  # GS分析|GS analysis
  janusx gs -i data.vcf.gz -p pheno.txt -m GBLUP rrBLUP -o output

  # PCA分析|PCA analysis
  janusx pca -i data.vcf.gz -d 5 --plot --plot3d -o output

  # PostGWAS可视化|PostGWAS visualization
  janusx postgwas -f result.lmm.tsv --threshold 1e-6 -o output
        '''
    )

    # 子命令|Subcommands
    subparsers = parser.add_subparsers(dest='mode', help='分析模式|Analysis mode')

    # 添加各模块参数|Add module parameters
    parse_gwas_arguments(subparsers)
    parse_gs_arguments(subparsers)
    parse_pca_arguments(subparsers)
    parse_postgwas_arguments(subparsers)

    args = parser.parse_args()

    # 检查模式|Check mode
    if not args.mode:
        parser.print_help()
        sys.exit(1)

    try:
        if args.mode == 'gwas':
            # 创建GWAS配置|Create GWAS config
            config = JanusXGWASConfig(
                genotype=args.genotype,
                pheno=args.pheno,
                genotype_type=args.type,
                models=args.models,
                maf=args.maf,
                geno=args.geno,
                grm=args.grm,
                qcov=args.qcov,
                cov=args.cov,
                ncol=args.ncol,
                plot=args.plot,
                chunksize=args.chunksize,
                mmap_limit=args.mmap_limit,
                threads=args.threads,
                output_dir=args.output_dir,
                prefix=args.prefix,
                janusx_path=args.janusx_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = JanusXGWASRunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

        elif args.mode == 'gs':
            # 创建GS配置|Create GS config
            config = JanusXGSConfig(
                genotype=args.genotype,
                pheno=args.pheno,
                genotype_type=args.type,
                models=args.models,
                maf=args.maf,
                geno=args.geno,
                pcd=args.pcd,
                ncol=args.ncol,
                cv=args.cv,
                plot=args.plot,
                output_dir=args.output_dir,
                prefix=args.prefix,
                janusx_path=args.janusx_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = JanusXGSRunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

        elif args.mode == 'pca':
            # 创建PCA配置|Create PCA config
            config = JanusXPCAConfig(
                genotype=args.genotype,
                genotype_type=args.type,
                grm=args.grm,
                pcfile=args.pcfile,
                dim=args.dim,
                plot=args.plot,
                plot3d=args.plot3d,
                group=args.group,
                color=args.color,
                output_dir=args.output_dir,
                prefix=args.prefix,
                janusx_path=args.janusx_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = JanusXPCARunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

        elif args.mode == 'postgwas':
            # 创建PostGWAS配置|Create PostGWAS config
            config = JanusXPostGWASConfig(
                files=args.files,
                chr_col=args.chr_col,
                pos_col=args.pos_col,
                pvalue_col=args.pvalue_col,
                threshold=args.threshold,
                noplot=args.noplot,
                color=args.color,
                highlight=getattr(args, 'hl', None),  # argparse使用hl作为属性名
                format=args.format,
                anno=args.anno,
                anno_broaden=getattr(args, 'ab', None),  # argparse使用ab作为属性名
                desc_item=args.desc_item,
                output_dir=args.output_dir,
                prefix=args.prefix,
                threads=args.threads,
                janusx_path=args.janusx_path
            )

            # 验证配置|Validate config
            config.validate()

            # 创建运行器并运行|Create runner and run
            runner = JanusXPostGWASRunner(config)
            success = runner.run()

            sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
