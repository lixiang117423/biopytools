"""泛基因组变异分析命令|Pan-genome Variant Analysis Command"""

import sys
import click


@click.command(
    short_help='泛基因组变异分析|Pan-genome Variant Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
)
@click.option('-i', '--input', 'input_file', required=True,
              help='[FILE] GBZ图文件或VCF文件|GBZ graph or VCF file')
@click.option('-o', '--output', 'output_dir', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('-P', '--ref-path', default=None,
              help='[STR] 参考路径前缀，GBZ输入时必需 (如T2T)|Reference path prefix (required for GBZ)')
@click.option('-t', '--threads', type=int, default=12,
              help='[INT] 线程数 (default: 12)')
@click.option('--ref-size', type=float, default=0.0,
              help='[FLOAT] 参考基因组大小Mb (0=自动推断, default: 0)')
@click.option('--permutations', type=int, default=100,
              help='[INT] 增长曲线随机置换次数 (default: 100)')
@click.option('--vg-env', default=None,
              help='[STR] vg conda环境名 (default: vg_v.1.7.0)')
@click.option('--r-path', default=None,
              help='[FILE] Rscript路径|Rscript binary path')
def panvar(input_file, output_dir, **kwargs):
    """泛基因组变异分析 - vg deconstruct + 变异统计 + 增长曲线+gamma|Pan-genome Variant Analysis

    \b
    自动识别输入类型: .gbz执行deconstruct+统计, .vcf直接统计
    Auto-detect input: .gbz runs deconstruct+summary, .vcf runs summary only

    示例|Examples: biopytools panvar -i cactus_output.gbz -P T2T -o output/
    """
    from ..panvar.config import PanvarConfig
    from ..panvar.runner import PanvarRunner

    config = PanvarConfig(
        input_file=input_file,
        output_dir=output_dir,
        ref_path=kwargs['ref_path'] or '',
        threads=kwargs['threads'],
        ref_size_mb=kwargs['ref_size'],
        n_permutations=kwargs['permutations'],
        vg_env=kwargs['vg_env'] if kwargs['vg_env'] else PanvarConfig().vg_env,
        r_path=kwargs['r_path'] if kwargs['r_path'] else PanvarConfig().r_path,
    )
    runner = PanvarRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)


__all__ = ['panvar']
