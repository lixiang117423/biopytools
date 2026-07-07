"""
连锁不平衡热图分析命令|LD Heatmap Analysis Command
"""

import click
import sys
import os


def _lazy_import_ldblockshow_main():
    """延迟加载ldblockshow主函数|Lazy load ldblockshow main function"""
    try:
        from ...ldblockshow.main import main as ldblockshow_main
        return ldblockshow_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(ctx, param, value):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and value and not os.path.exists(value):
        raise click.BadParameter(f"文件不存在|File does not exist: {value}")
    return value


@click.command(
    short_help='连锁不平衡热图分析工具|LD heatmap analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
# 必需参数|Required parameters
@click.option('-i', '--vcf-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value),
              help='VCF变异文件路径|VCF variant file path')
@click.option('-o', '--output-dir',
              required=True,
              help='输出目录(自动创建)；每 region 产物落在 目录/<label>.* '
              '|Output directory (auto-created); per-region outputs land in dir/<label>.*')
@click.option('-r', '--region',
              help='分析区域，格式chr:start-end|Analysis region, format chr:start-end')
@click.option('-b', '--bed',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='基因组BED文件(每行 chrom start end [name])，等价多个 -r 批量出图'
              '|Genomic BED (cols: chrom start end [name]), equivalent to multiple -r')
# 输入文件|Input files
@click.option('--in-genotype',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='SNP Genotype格式文件路径|SNP Genotype format file path')
@click.option('--in-plink',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='Plink文件前缀(bed+bim+fam或ped+map)|Plink file prefix')
# LD度量|LD statistic
@click.option('--sele-var',
              type=click.Choice(['1', '2', '3', '4']),
              default='1',
              show_default=True,
              help='选择LD度量统计量|Select LD statistic (1=D\', 2=R², 3/4=Both)')
# 过滤参数|Filter parameters
@click.option('--maf', type=float, default=0.05, show_default=True,
              help='最小次要等位基因频率|Minimum minor allele frequency')
@click.option('--miss', type=float, default=0.25, show_default=True,
              help='最大缺失率|Maximum missing ratio')
@click.option('--hwe', type=float, default=0.0, show_default=True,
              help='Hardy-Weinberg平衡P值阈值|Hardy-Weinberg equilibrium P-value threshold')
@click.option('--het', type=float, default=1.0, show_default=True,
              help='最大杂合率|Maximum heterozygosity ratio')
@click.option('--enable-oth-var', is_flag=True,
              help='允许indel/SV/CNV变异|Allow bi-indel bi-sv bi-cnv variants')
# Block检测|Block detection
@click.option('--block-type',
              type=click.Choice(['1', '2', '3', '4', '5']),
              default='1',
              show_default=True,
              help='Block检测方法 (1=Gabriel, 2=SolidSpine, 3=BlockCut, 4=FixBlock, 5=NoBlock)')
@click.option('--block-cut', default="0.85:0.90", show_default=True,
              help='BlockType3的cutoff|Cutoff for BlockType3')
@click.option('--fix-block',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='固定block文件路径|Fixed block file path')
# 可视化参数|Visualization parameters
@click.option('--in-gwas',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='GWAS文件，3列: chr pos pvalue，chr名须与VCF一致|GWAS file, 3 cols: chr pos pvalue')
@click.option('--in-gff',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='GFF3注释文件路径|GFF3 annotation file path')
@click.option('--mer-min-snp-num', type=int, default=50, show_default=True,
              help='合并网格的最小SNP数|Minimum SNP number to merge grids')
# 输出格式|Output format
@click.option('--no-out-png', is_flag=True,
              help='不输出PNG格式图像|Do not output PNG format image')
@click.option('--out-pdf', is_flag=True,
              help='输出PDF格式图像|Output PDF format image')
# 其他参数|Other parameters
@click.option('--sub-pop',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='亚群样本文件路径|Subgroup sample file path')
@click.option('--tag-snp-cut', type=float, default=0.80, show_default=True,
              help='TagSNP的LD cutoff|LD cutoff for TagSNP')
# 绘图参数|Drawing parameters (ShowLDSVG)
@click.option('--cutline', type=float, default=5.0, show_default=True,
              help='GWAS P值显著性阈值(-log10)|GWAS P-value significance cutoff (-log10)')
@click.option('--point-size', type=int,
              help='GWAS散点大小|GWAS point size')
@click.option('--top-site',
              help='指定GWAS峰值位点(chr:pos)|Specify GWAS peak site (chr:pos)')
@click.option('--no-log-p', is_flag=True,
              help='不对P值取-log10|Do not -log10 transform P-value')
@click.option('--no-gene-name', is_flag=True,
              help='不显示基因名|Do not show gene names')
@click.option('--show-num', is_flag=True,
              help='在热图中显示R²/D\'值|Show R²/D\' values in heatmap')
@click.option('--spe-snp-name',
              callback=lambda ctx, param, value: _validate_file_exists(ctx, param, value) if value else None,
              help='特殊SNP名称文件(chr site Name)|Special SNP name file')
@click.option('--show-gwas-spe-snp', is_flag=True,
              help='在GWAS图中显示特殊SNP名称|Show special SNP names in GWAS plot')
@click.option('--resize-h', type=int,
              help='图像高度，宽度按比例自动调整|Image height, width auto-adjusted')
@click.option('--no-show-ldist', type=int,
              help='超过此距离的SNP对不显示LD|NoShow pairwise LD over this distance')
def ldblockshow(vcf_file, output_dir, region, bed, in_genotype, in_plink,
                sele_var, maf, miss, hwe, het, enable_oth_var,
                block_type, block_cut, fix_block,
                in_gwas, in_gff, mer_min_snp_num,
                no_out_png, out_pdf,
                sub_pop, tag_snp_cut,
                cutline, point_size, top_site, no_log_p, no_gene_name,
                show_num, spe_snp_name, show_gwas_spe_snp, resize_h, no_show_ldist):
    """
    连锁不平衡热图分析工具|LD Heatmap Analysis Tool

    基于VCF文件生成连锁不平衡(LD)热图，支持GWAS和基因注释可视化|Generate linkage disequilibrium (LD) heatmap from VCF files, support GWAS and gene annotation visualization

    示例|Example: biopytools ldblockshow -i variants.vcf.gz -o ld_result -r chr1:1000000-2000000
    """

    ldblockshow_main = _lazy_import_ldblockshow_main()

    args = ['ldblockshow.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf_file])
    args.extend(['-o', output_dir])
    # region 与 bed 二选一(由 main 校验)|region XOR bed (validated by main)
    if region:
        args.extend(['-r', region])
    if bed:
        args.extend(['-b', bed])

    # 输入文件|Input files
    if in_genotype:
        args.extend(['--in-genotype', in_genotype])
    if in_plink:
        args.extend(['--in-plink', in_plink])

    # LD度量选择|LD statistic selection
    if sele_var != '1':
        args.extend(['--sele-var', sele_var])

    # 过滤参数|Filter parameters
    if maf != 0.05:
        args.extend(['--maf', str(maf)])
    if miss != 0.25:
        args.extend(['--miss', str(miss)])
    if hwe != 0.0:
        args.extend(['--hwe', str(hwe)])
    if het != 1.0:
        args.extend(['--het', str(het)])
    if enable_oth_var:
        args.append('--enable-oth-var')

    # Block检测参数|Block detection parameters
    if block_type != '1':
        args.extend(['--block-type', block_type])
    if block_type == '3' and block_cut != "0.85:0.90":
        args.extend(['--block-cut', block_cut])
    if block_type == '4' and fix_block:
        args.extend(['--fix-block', fix_block])

    # 可视化参数|Visualization parameters
    if in_gwas:
        args.extend(['--in-gwas', in_gwas])
    if in_gff:
        args.extend(['--in-gff', in_gff])
    if mer_min_snp_num != 50:
        args.extend(['--mer-min-snp-num', str(mer_min_snp_num)])

    # 输出格式|Output format
    if no_out_png:
        args.append('--no-out-png')
    if out_pdf:
        args.append('--out-pdf')

    # 其他参数|Other parameters
    if sub_pop:
        args.extend(['--sub-pop', sub_pop])
    if tag_snp_cut != 0.80:
        args.extend(['--tag-snp-cut', str(tag_snp_cut)])

    # 绘图参数|Drawing parameters
    if cutline != 5.0:
        args.extend(['--cutline', str(cutline)])
    if point_size:
        args.extend(['--point-size', str(point_size)])
    if top_site:
        args.extend(['--top-site', top_site])
    if no_log_p:
        args.append('--no-log-p')
    if no_gene_name:
        args.append('--no-gene-name')
    if show_num:
        args.append('--show-num')
    if spe_snp_name:
        args.extend(['--spe-snp-name', spe_snp_name])
    if show_gwas_spe_snp:
        args.append('--show-gwas-spe-snp')
    if resize_h:
        args.extend(['--resize-h', str(resize_h)])
    if no_show_ldist:
        args.extend(['--no-show-ldist', str(no_show_ldist)])

    original_argv = sys.argv
    sys.argv = args

    try:
        ldblockshow_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
