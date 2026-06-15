"""
PopLDdecay CLI命令|PopLDdecay CLI Command
"""

import click


@click.command()
@click.option('-i', '--input', required=True, type=click.Path(exists=True), help='输入VCF或Genotype文件|Input VCF or Genotype file')
@click.option('-o', '--output', required=True, type=click.Path(), help='输出文件前缀|Output file prefix')
@click.option('-t', '--type', type=click.Choice(['vcf', 'genotype']), default='vcf', help='输入文件类型|Input file type (vcf/genotype) [default: vcf]')
@click.option('-d', '--max-dist', default=300, type=int, help='最大距离(kb)|Max distance in kb [default: 300]')
@click.option('-m', '--min-maf', default=0.005, type=float, help='最小等位基因频率|Min minor allele frequency [default: 0.005]')
@click.option('--max-het', default=0.88, type=float, help='最大杂合率|Max heterozygous rate [default: 0.88]')
@click.option('--max-miss', default=0.25, type=float, help='最大缺失率|Max missing rate [default: 0.25]')
@click.option('-s', '--subpop', type=click.Path(exists=True), help='子群体样本列表文件|Subpopulation sample list file')
@click.option('--out-type', default=1, type=int, help='输出类型(1:R^2, 2:R^2&D\', 3:Pairwise LD)|Output type [default: 1]')
@click.option('--bin1', default=10, type=int, help='短距离bin大小|Bin size for short distance [default: 10]')
@click.option('--bin2', default=100, type=int, help='长距离bin大小|Bin size for long distance [default: 100]')
@click.option('--break-point', default=100, type=int, help='短/长距离分界点|Break point [default: 100]')
@click.option('--max-x', type=int, help='最大X坐标(kb)|Max X coordinate in kb')
@click.option('--measure', type=click.Choice(['r2', 'D', 'both']), default='r2', help='LD度量方法|LD measure method [default: r2]')
@click.option('--method', type=click.Choice(['MeanBin', 'HW', 'MedianBin', 'PercentileBin']), default='MeanBin', help='绘图方法|Plotting method [default: MeanBin]')
@click.option('--percentile', default=0.5, type=float, help='百分位数|Percentile for PercentileBin [default: 0.5]')
@click.option('--no-plot', is_flag=True, help='不绘制图像|Do not plot figure')
@click.option('--no-recommend-threshold', is_flag=True, help='不推荐LD阈值|Do not recommend LD threshold')
def poplddecay(input, output, type, max_dist, min_maf, max_het, max_miss, subpop, out_type, bin1, bin2, break_point, max_x, measure, method, percentile, no_plot, no_recommend_threshold):
    """
    连锁不平衡衰减分析工具|Linkage Disequilibrium Decay Analysis Tool

    示例|Example: biopytools poplddecay -i variants.vcf -o output
    """
    # 延迟加载|Lazy loading
    from biopytools.poplddecay.main import main as poplddecay_main

    # 创建args对象|Create args object
    class Args:
        pass

    args = Args()
    args.input = input
    args.output = output
    args.type = type
    args.max_dist = max_dist
    args.min_maf = min_maf
    args.max_het = max_het
    args.max_miss = max_miss
    args.subpop = subpop
    args.out_type = out_type
    args.bin1 = bin1
    args.bin2 = bin2
    args.break_point = break_point
    args.max_x = max_x
    args.measure = measure
    args.method = method
    args.percentile = percentile
    args.no_plot = no_plot
    args.recommend_threshold = not no_recommend_threshold  # 默认开启

    return poplddecay_main(args)
