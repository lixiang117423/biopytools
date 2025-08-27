"""
群体遗传分析命令 | Population Genetics Analysis Command
"""

import click
import sys
from ...popgen_toolkit.main import main as popgen_main


@click.command(short_help='群体遗传学多样性分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径 | Input VCF file path')
@click.option('--output', '-o',
              default='./popgen_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./popgen_output)')
@click.option('--groups', '-g',
              type=click.Path(exists=True),
              help='分组信息文件 | Group information file')
# @click.option('--skip-qc',
#               is_flag=True,
#               help='跳过质量控制过滤，直接使用输入VCF文件 | Skip quality control filtering, use input VCF file directly')
@click.option('--skip-qc',
              is_flag=True,
              default=True,  # 改为True，与config.py一致
              help='跳过质量控制过滤 | Skip quality control filtering (default)')
# @click.option('--all',
#               is_flag=True,
#               default=False,  # 改为False
#               help='计算所有参数 | Calculate all parameters')
@click.option('--fst',
              is_flag=True,
              help='计算Fst | Calculate Fst')
@click.option('--pi',
              is_flag=True,
              help='计算π | Calculate π')
@click.option('--theta-w',
              is_flag=True,
              help='计算θw | Calculate θw')
@click.option('--tajima-d',
              is_flag=True,
              help='计算Tajima\'s D | Calculate Tajima\'s D')
@click.option('--ibd',
              is_flag=True,
              help='计算IBD | Calculate IBD')
@click.option('--ld',
              is_flag=True,
              help='计算LD | Calculate LD')
# @click.option('--ne',
#               is_flag=True,
#               help='计算有效群体大小 | Calculate effective population size')
@click.option('--windows', '-w',
              multiple=True,
              type=int,
              default=[10000, 100000, 500000],
              help='滑动窗口大小 (bp) | Sliding window sizes (bp) (default: 10000,100000,500000)')
@click.option('--overlap',
              default=0.9,
              type=float,
              help='窗口重叠率 | Window overlap rate (default: 0.9)')
@click.option('--maf', '-m',
              default=0.01,
              type=float,
              help='MAF阈值 | MAF threshold (default: 0.01)')
@click.option('--missing', '-M',
              default=0.1,
              type=float,
              help='缺失率阈值 | Missing rate threshold (default: 0.1)')
@click.option('--hwe', '-H',
              default=1e-6,
              type=float,
              help='HWE p值阈值 | HWE p-value threshold (default: 1e-6)')
@click.option('--min-dp',
              default=10,
              type=int,
              help='最小测序深度 | Minimum depth (default: 10)')
@click.option('--max-dp',
              default=100,
              type=int,
              help='最大测序深度 | Maximum depth (default: 100)')
@click.option('--format', '-f',
              default='txt',
              type=click.Choice(['txt', 'csv', 'tsv', 'json']),
              help='输出文件格式 | Output file format (default: txt)')
@click.option('--threads', '-t',
              default=4,
              type=int,
              help='线程数 | Number of threads (default: 4)')
@click.option('--vcftools-path',
              default='vcftools',
              type=str,
              help='VCFtools路径 | VCFtools path (default: vcftools)')
@click.option('--plink-path',
              default='plink',
              type=str,
              help='PLINK路径 | PLINK path (default: plink)')
@click.option('--bcftools-path',
              default='bcftools',
              type=str,
              help='BCFtools路径 | BCFtools path (default: bcftools)')
@click.option('--smcpp-path',
              default='smc++',
              type=str,
              help='SMC++路径 | SMC++ path (default: smc++)')
def popgen(vcf, output, groups, all, fst, pi, theta_w, tajima_d, ibd, ld, ne,
           windows, overlap, maf, missing, hwe, min_dp, max_dp, format,
           threads, vcftools_path, plink_path, bcftools_path, smcpp_path, skip_qc):
    """
    群体遗传分析工具 (模块化版本)
    
    基于VCF文件进行全面的群体遗传学分析，计算多种遗传多样性参数，
    包括核苷酸多样性、群体分化、连锁不平衡、有效群体大小等指标。
    
    功能特点 | Features:
    - 多种遗传多样性参数计算
    - 滑动窗口分析支持
    - 群体分组比较分析
    - 质量控制和数据过滤
    - 多种输出格式支持
    - 并行计算加速
    
    输出文件 | Output Files:
    - diversity_*.txt: 各种多样性参数结果
    - fst_results.txt: 群体分化指数结果
    - ibd_results.txt: 同源性分析结果
    - ld_results.txt: 连锁不平衡分析结果
    - ne_results.txt: 有效群体大小结果
    - summary_report.txt: 综合分析报告
    
    示例 | Examples:
    
    \b
    # 基本分析（计算所有参数）
    biopytools popgen -v variants.vcf.gz -o popgen_results
    
    \b
    # 带分组信息的分析
    biopytools popgen -v data.vcf.gz -o results -g groups.txt --fst --pi --ld
    
    \b
    # 自定义参数的全面分析
    biopytools popgen -v variants.vcf -o results -g samples_groups.txt \\
        --all -t 8 --format csv
    
    \b
    # 严格质控的多样性分析
    biopytools popgen -v large_dataset.vcf.gz -o strict_analysis \\
        --maf 0.05 --missing 0.05 --hwe 1e-8 \\
        --min-dp 20 --max-dp 200 -t 32
    
    \b
    # 自定义滑动窗口分析
    biopytools popgen -v genome.vcf.gz -o window_analysis \\
        -w 50000 -w 200000 -w 1000000 --overlap 0.8 \\
        --pi --theta-w --tajima-d
    
    \b
    # 指定外部工具路径
    biopytools popgen -v data.vcf.gz -o custom_tools \\
        --vcftools-path /opt/vcftools/bin/vcftools \\
        --plink-path /usr/local/bin/plink \\
        --smcpp-path /home/user/smc++/smc++
    
    分组文件格式 | Group File Format:
    sample1    group1
    sample2    group1
    sample3    group2
    sample4    group2
    
    参数说明 | Parameters Description:
    
    遗传多样性参数:
    - π (pi): 核苷酸多样性，衡量序列变异程度
    - θw (theta-w): Watterson's theta，基于分离位点数的多样性估计
    - Tajima's D: 检测选择压力和群体历史的统计量
    - Fst: 群体分化指数，衡量群体间遗传差异
    - IBD: 同源性分析，个体间亲缘关系
    - LD: 连锁不平衡，基因座间的关联程度
    - Ne: 有效群体大小，影响遗传漂变的个体数
    
    质控参数:
    - MAF: 最小等位基因频率，过滤稀有变异
    - Missing: 缺失率阈值，移除数据质量差的位点
    - HWE: Hardy-Weinberg平衡检验，检测基因型频率偏离
    - Depth: 测序深度范围，确保数据可靠性
    
    滑动窗口分析:
    - 支持多种窗口大小同时分析
    - 可调节窗口重叠率
    - 适用于基因组尺度的变异模式分析
    
    分析流程 | Analysis Pipeline:
    1. VCF文件质量控制和预处理
    2. 样本分组信息加载（如果提供）
    3. 多样性参数计算（π, θw, Tajima's D）
    4. 群体分化分析（Fst计算）
    5. 同源性分析（IBD计算）
    6. 连锁不平衡分析（LD计算）
    7. 有效群体大小估算（SMC++）
    8. 结果整合和报告生成
    
    外部依赖工具 | External Dependencies:
    - VCFtools: VCF文件处理和统计计算
    - PLINK: 群体遗传学分析
    - BCFtools: VCF文件操作
    - SMC++: 有效群体大小估算
    
    性能建议 | Performance Tips:
    - 大数据集建议增加线程数和内存
    - 根据研究目标选择合适的窗口大小
    - 严格的质控参数有助于提高结果可靠性
    - 分组分析需要平衡的样本设计
    """
    
    # 构建参数列表传递给原始main函数
    args = ['popgen.py']
    
    # 必需参数
    args.extend(['-v', vcf])
    
    # 可选参数（只在非默认值时添加）
    if output != './popgen_output':
        args.extend(['-o', output])
    
    if groups:
        args.extend(['-g', groups])
    
    # 分析选择参数
    # if not all:
    #     if fst:
    #         args.append('--fst')
    #     if pi:
    #         args.append('--pi')
    #     if theta_w:
    #         args.append('--theta-w')
    #     if tajima_d:
    #         args.append('--tajima-d')
    #     if ibd:
    #         args.append('--ibd')
    #     if ld:
    #         args.append('--ld')
    #     if ne:
    #         args.append('--ne')
    # elif all and not (fst or pi or theta_w or tajima_d or ibd or ld or ne):
    #     args.append('--all')
    # 分析选择参数
    specific_params = [fst, pi, theta_w, tajima_d, ibd, ld, ne]
    has_specific_params = any(specific_params)

    if all and not has_specific_params:
        # 只指定了--all，没有指定具体参数
        args.append('--all')
    elif has_specific_params:
        # 指定了具体的参数，忽略--all标志
        if fst:
            args.append('--fst')
        if pi:
            args.append('--pi')
        if theta_w:
            args.append('--theta-w')
        if tajima_d:
            args.append('--tajima-d')
        if ibd:
            args.append('--ibd')
        if ld:
            args.append('--ld')
        if ne:
            args.append('--ne')
    else:
        # 什么都没指定，默认计算所有
         args.append('--all')

    # 质控选项
    if skip_qc:
        args.append('--skip-qc')
    
    # if skip_qc:  # 用户明确指定要跳过
    #     args.append('--skip-qc')
    # elif not skip_qc:  # 用户没有指定或明确指定不跳过
    #     # 什么都不做，让main.py使用其默认值
    #     pass
    
    # 滑动窗口参数
    if windows and set(windows) != {10000, 100000, 500000}:
        for window in windows:
            args.extend(['-w', str(window)])
    
    if overlap != 0.9:
        args.extend(['--overlap', str(overlap)])
    
    # 质控参数
    if maf != 0.01:
        args.extend(['-m', str(maf)])
    
    if missing != 0.1:
        args.extend(['-M', str(missing)])
    
    if hwe != 1e-6:
        args.extend(['-H', str(hwe)])
    
    if min_dp != 10:
        args.extend(['--min-dp', str(min_dp)])
    
    if max_dp != 100:
        args.extend(['--max-dp', str(max_dp)])
    
    # 输出格式
    if format != 'txt':
        args.extend(['-f', format])
    
    # 计算资源
    if threads != 4:
        args.extend(['-t', str(threads)])
    
    # 工具路径
    if vcftools_path != 'vcftools':
        args.extend(['--vcftools-path', vcftools_path])
    
    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])
    
    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])
    
    if smcpp_path != 'smc++':
        args.extend(['--smcpp-path', smcpp_path])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        popgen_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"分析失败 | Analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv