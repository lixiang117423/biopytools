"""TreeMix CLI组定义|TreeMix CLI Group Definition

biopytools treemix {prepare,run,all}
"""

import sys
import click


@click.group(
    short_help='TreeMix群体历史与基因流分析|TreeMix Population History & Gene Flow',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
    invoke_without_command=True,
)
def treemix():
    """TreeMix群体历史与基因流分析 - 输入准备/运行/完整流程|TreeMix Population History & Gene Flow

    \b
    示例|Examples:
        biopytools treemix all -i input.vcf.gz -o output --cluster cluster.txt
        biopytools treemix prepare -i input.vcf.gz -o output
        biopytools treemix run -i input.frq.gz -o output --m-max 5 --root outgroup
    """
    pass


@treemix.command(short_help='输入准备|Input Preparation')
@click.option('-i', '--input', 'vcf_file', required=True,
              help='[FILE] VCF输入文件|VCF input file')
@click.option('-o', '--output', 'output_dir', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('--cluster', default=None,
              help='[FILE] 样本分组文件 (sample_id population); 不提供则自动推断|Cluster file')
@click.option('--pop-delimiter', default='_',
              help='[STR] 自动推断群体时的分隔符 (default: _)')
@click.option('--ld-window', type=int, default=50,
              help='[INT] LD窗口SNP数 (default: 50)')
@click.option('--ld-step', type=int, default=10,
              help='[INT] LD步长 (default: 10)')
@click.option('--ld-r2', type=float, default=0.2,
              help='[FLOAT] LD r2阈值 (default: 0.2)')
@click.option('--treemix-path', default=None,
              help='[FILE] treemix路径|treemix binary path')
@click.option('--plink-path', default=None,
              help='[FILE] plink路径|plink binary path')
@click.option('--bcftools-path', default=None,
              help='[FILE] bcftools路径|bcftools binary path')
def prepare(vcf_file, output_dir, **kwargs):
    """输入准备 - VCF转TreeMix格式|Prepare TreeMix input from VCF

    \b
    流程: 按分组文件筛选样品 → LD过滤 → plink频率计算 → TreeMix格式转换

    示例|Examples: biopytools treemix prepare -i input.vcf.gz -o output --cluster cluster.txt
    """
    from .config import TreemixConfig
    from .runner import TreemixRunner

    config = TreemixConfig(
        vcf_file=vcf_file,
        output_dir=output_dir,
        cluster_file=kwargs.get('cluster') or '',
        pop_delimiter=kwargs['pop_delimiter'],
        ld_window=kwargs['ld_window'],
        ld_step=kwargs['ld_step'],
        ld_r2=kwargs['ld_r2'],
        treemix_path=kwargs['treemix_path'] if kwargs['treemix_path'] else TreemixConfig().treemix_path,
        plink_path=kwargs['plink_path'] if kwargs['plink_path'] else TreemixConfig().plink_path,
        bcftools_path=kwargs['bcftools_path'] if kwargs['bcftools_path'] else TreemixConfig().bcftools_path,
    )
    runner = TreemixRunner(config)
    success = runner.run_prepare()
    sys.exit(0 if success else 1)


@treemix.command(short_help='运行TreeMix|Run TreeMix')
@click.option('-i', '--input', 'treemix_input', required=True,
              help='[FILE] TreeMix输入文件 (.frq.gz)|TreeMix input file')
@click.option('-o', '--output', 'output_dir', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('--m-max', type=int, default=10,
              help='[INT] 测试m=0..m_max (default: 10)')
@click.option('-m', 'm_value', type=int, default=None,
              help='[INT] 指定m值 (单独运行或绘图)|Specific m value')
@click.option('--root', default=None,
              help='[STR] 外群群体名|Outgroup population name')
@click.option('--bootstrap', type=int, default=1000,
              help='[INT] bootstrap重复次数 (default: 1000)')
@click.option('--replicates', type=int, default=10,
              help='[INT] 每个m值的重复次数 (default: 10)')
@click.option('--noss', is_flag=True, default=False,
              help='关闭样本量校正|Turn off sample size correction')
@click.option('--global', 'global_opt', is_flag=True, default=False,
              help='全局拓扑重排 (默认开启, 加此选项强制开启)|Global tree rearrangements')
@click.option('--no-global', 'no_global_opt', is_flag=True, default=False,
              help='关闭全局拓扑重排|Disable global tree rearrangements')
@click.option('--se', is_flag=True, default=False,
              help='计算迁移权重标准误|Calculate standard errors')
@click.option('-k', 'k_value', type=int, default=500,
              help='[INT] SNP block大小 (default: 500)')
@click.option('-t', '--threads', type=int, default=12,
              help='[INT] 并行线程数 (default: 12)')
@click.option('--seed', type=int, default=None,
              help='[INT] 随机种子|Random seed')
@click.option('--treemix-path', default=None,
              help='[FILE] treemix路径|treemix binary path')
@click.option('--r-path', default=None,
              help='[FILE] R脚本路径|R binary path')
@click.option('--plotting-funcs', default=None,
              help='[FILE] plotting_funcs.R路径|plotting_funcs.R path')
def run(treemix_input, output_dir, **kwargs):
    """运行TreeMix - scan + OptM + 绘图|Run TreeMix: scan + OptM + plot

    \b
    完整流程:
      1. Scan: m=0..m_max, 每个m运行N次重复 (含bootstrap)
      2. OptM: 用R包OptM确定最优m值
      3. Plot: 绘制最优m的树图和残差图

    示例|Examples: biopytools treemix run -i input.frq.gz -o output --root outgroup
    """
    from .config import TreemixConfig
    from .runner import TreemixRunner

    global_opt = kwargs['global_opt'] or (not kwargs['no_global_opt'])

    config = TreemixConfig(
        treemix_input=treemix_input,
        output_dir=output_dir,
        m_max=kwargs['m_max'],
        m=kwargs['m_value'],
        root=kwargs.get('root') or '',
        bootstrap=kwargs['bootstrap'],
        replicates=kwargs['replicates'],
        noss=kwargs['noss'],
        global_opt=global_opt,
        se=kwargs['se'],
        k=kwargs['k_value'],
        threads=kwargs['threads'],
        seed=kwargs['seed'],
        treemix_path=kwargs['treemix_path'] if kwargs['treemix_path'] else TreemixConfig().treemix_path,
        r_path=kwargs['r_path'] if kwargs['r_path'] else TreemixConfig().r_path,
        plotting_funcs_r=kwargs['plotting_funcs'] or '',
    )
    runner = TreemixRunner(config)
    success = runner.run_treemix()
    sys.exit(0 if success else 1)


@treemix.command(short_help='完整流程|Full Pipeline')
@click.option('-i', '--input', 'vcf_file', required=True,
              help='[FILE] VCF输入文件|VCF input file')
@click.option('-o', '--output', 'output_dir', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('--cluster', default=None,
              help='[FILE] 样本分组文件 (sample_id population); 不提供则自动推断|Cluster file')
@click.option('--pop-delimiter', default='_',
              help='[STR] 自动推断群体时的分隔符 (default: _)')
@click.option('--m-max', type=int, default=10,
              help='[INT] 测试m=0..m_max (default: 10)')
@click.option('-m', 'm_value', type=int, default=None,
              help='[INT] 指定m值 (单独运行或绘图)|Specific m value')
@click.option('--root', default=None,
              help='[STR] 外群群体名|Outgroup population name')
@click.option('--bootstrap', type=int, default=1000,
              help='[INT] bootstrap重复次数 (default: 1000)')
@click.option('--replicates', type=int, default=10,
              help='[INT] 每个m值的重复次数 (default: 10)')
@click.option('--noss', is_flag=True, default=False,
              help='关闭样本量校正|Turn off sample size correction')
@click.option('--global', 'global_opt', is_flag=True, default=False,
              help='全局拓扑重排 (默认开启, 加此选项强制开启)|Global tree rearrangements')
@click.option('--no-global', 'no_global_opt', is_flag=True, default=False,
              help='关闭全局拓扑重排|Disable global tree rearrangements')
@click.option('--se', is_flag=True, default=False,
              help='计算迁移权重标准误|Calculate standard errors')
@click.option('-k', 'k_value', type=int, default=500,
              help='[INT] SNP block大小 (default: 500)')
@click.option('-t', '--threads', type=int, default=12,
              help='[INT] 并行线程数 (default: 12)')
@click.option('--seed', type=int, default=None,
              help='[INT] 随机种子|Random seed')
@click.option('--treemix-path', default=None,
              help='[FILE] treemix路径|treemix binary path')
@click.option('--plink-path', default=None,
              help='[FILE] plink路径|plink binary path')
@click.option('--r-path', default=None,
              help='[FILE] R脚本路径|R binary path')
@click.option('--plotting-funcs', default=None,
              help='[FILE] plotting_funcs.R路径|plotting_funcs.R path')
@click.option('--bcftools-path', default=None,
              help='[FILE] bcftools路径|bcftools binary path')
def all_cmd(vcf_file, output_dir, **kwargs):
    """完整流程 - 输入准备 + scan + OptM + 绘图|Full pipeline: prepare + scan + OptM + plot

    \b
    完整流程:
      1. Prepare: 样品筛选 → LD过滤 → plink频率计算 → TreeMix格式转换
      2. Scan: m=0..m_max, 每个m运行N次重复 (含bootstrap)
      3. OptM: 用R包OptM确定最优m值
      4. Plot: 绘制m=0和最优m的树图和残差图

    示例|Examples: biopytools treemix all -i input.vcf.gz -o output --cluster cluster.txt
    """
    from .config import TreemixConfig
    from .runner import TreemixRunner

    global_opt = kwargs['global_opt'] or (not kwargs['no_global_opt'])

    config = TreemixConfig(
        vcf_file=vcf_file,
        output_dir=output_dir,
        cluster_file=kwargs.get('cluster') or '',
        pop_delimiter=kwargs['pop_delimiter'],
        m_max=kwargs['m_max'],
        m=kwargs['m_value'],
        root=kwargs.get('root') or '',
        bootstrap=kwargs['bootstrap'],
        replicates=kwargs['replicates'],
        noss=kwargs['noss'],
        global_opt=global_opt,
        se=kwargs['se'],
        k=kwargs['k_value'],
        threads=kwargs['threads'],
        seed=kwargs['seed'],
        treemix_path=kwargs['treemix_path'] if kwargs['treemix_path'] else TreemixConfig().treemix_path,
        plink_path=kwargs['plink_path'] if kwargs['plink_path'] else TreemixConfig().plink_path,
        r_path=kwargs['r_path'] if kwargs['r_path'] else TreemixConfig().r_path,
        plotting_funcs_r=kwargs['plotting_funcs'] or '',
        bcftools_path=kwargs['bcftools_path'] if kwargs['bcftools_path'] else TreemixConfig().bcftools_path,
    )
    runner = TreemixRunner(config)
    success = runner.run_all()
    sys.exit(0 if success else 1)
