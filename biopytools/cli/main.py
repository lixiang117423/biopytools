#!/usr/bin/env python3
"""
BioPyTools 统一CLI入口点 | BioPyTools Unified CLI Entry Point
"""

import click
from .._version import __version__

def import_available_commands():
    """动态导入可用的命令模块"""
    available_commands = {}
    
    command_modules = [
        # (模块文件名, 命令名, 描述文本)
        ('admixture', 'admixture', '🧬 ADMIXTURE群体结构分析'),
        ('annovar', 'annovar', '📝 ANNOVAR变异注释'),
        ('blast', 'blast', '🧬 BLAST序列比对分析'),
        ('coverage', 'coverage', '📊 BAM覆盖度分析'),
        ('ena_downloader', 'ena-downloader', '📥 ENA数据下载工具'),
        ('fastp', 'fastp', '🧹 FASTQ数据质量控制'),
        ('genomesyn', 'genomesyn', '🗺️ 基因组共线性分析'),
        ('geneinfo', 'geneinfo', '📄 从GFF文件提取基因信息'),
        ('gtx', 'gtx', '🔬 运行GTX WGS流程'),
        ('hifiasm', 'hifiasm', '🧩 运行hifiasm基因组组装'),
        ('kaks', 'kaks', '🧮 Ka/Ks计算'),
        ('kmer_count', 'kmer-count', '🔢 K-mer丰度矩阵计算'),
        ('kmer_query', 'kmer-query', '✂️ K-mer提取'),
        ('longestmrna', 'longest-mrna', '📜 提取最长转录本'),
        ('minimap2', 'minimap2', '🔗 Minimap2比对与区域提取'),
        ('plinkgwas', 'plink-gwas', '📈 PLINK GWAS分析'),
        ('popgen', 'popgen', '🌍 群体遗传学多样性分析'),
        ('rnaseq', 'rnaseq', 'RNA-seq表达定量流程'),
        ('split_fasta_id', 'split-fasta-id', '🔪 分割FASTA文件ID'),
        ('vcf_filter', 'vcf-filter', '🩸 VCF文件筛选'),
        ('vcf_genotype', 'vcf-genotype', '🔬 VCF基因型提取'),
        ('vcf_pca', 'vcf-pca', '📊 VCF主成分分析 (PCA)'),
        ('vcf_nj_tree', 'vcf-nj-tree', '🌳 VCF构建NJ进化树'),
        ('vcf_sample_hete', 'vcf-sample-hete', '📈 VCF样本基因型统计'),
        ('vcf_sequence', 'vcf-sequence', '🧬 从基因组和VCF提取序列'),
        ('bismark', 'bismark', '🧬 全基因组甲基化'),
        ('transcriptome_prediction', 'mrna-prediction', '🧬 基于转录组的基因预测'),
    ]
    
    for module_name, command_name, description in command_modules:
        try:
            module = __import__(f'biopytools.cli.commands.{module_name}', 
                              fromlist=[module_name])
            command_func = getattr(module, module_name)
            available_commands[command_name] = (command_func, description)
        except (ImportError, AttributeError) as e:
            # --- 关键修改在这里 ---
            # 不再沉默！打印出错误警告，帮助我们调试。
            click.secho(f" [!] 警告: 无法加载命令 '{command_name}'. 错误: {e}", fg='yellow', err=True)
            pass
    
    return available_commands


# --- 后续代码保持不变 ---

@click.group(invoke_without_command=True)
@click.version_option(version=__version__, prog_name='biopytools')
@click.pass_context
def cli(ctx):
    """
    BioPyTools - 生物信息学分析工具包
    ... (文档字符串保持不变) ...
    """
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        # 您的自定义帮助显示逻辑保持不变
        available_commands = import_available_commands()
        if available_commands:
            click.echo(f"\n当前可用的命令 | Currently available commands:")
            for cmd_name in sorted(available_commands.keys()):
                _, description = available_commands[cmd_name]
                click.echo(f"  {cmd_name:<20} {description}")
        # ... (您的“即将推出”部分也保持不变) ...

def register_commands():
    """注册所有可用的子命令"""
    available_commands = import_available_commands()
    for command_name, (command_func, description) in available_commands.items():
        # 这里我们恢复使用 docstring 的方式，因为那是您最初正确的方式
        # 如果您想强制使用列表里的描述，可以取消下面一行的注释
        # command_func.short_help = description 
        cli.add_command(command_func, name=command_name)

def main():
    """主入口函数"""
    register_commands()
    cli()

if __name__ == '__main__':
    main()