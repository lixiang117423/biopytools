# #!/usr/bin/env python3
# """
# BioPyTools 统一CLI入口点 | BioPyTools Unified CLI Entry Point
# """

# import click
# from .._version import __version__

# def import_available_commands():
#     """动态导入可用的命令模块"""
#     available_commands = {}
    
#     command_modules = [
#         # (模块文件名, 命令名, 描述文本)
#         ('admixture', 'admixture', '🧬 ADMIXTURE群体结构分析'),
#         ('annovar', 'annovar', '📝 ANNOVAR变异注释'),
#         ('blast', 'blast', '🧬 BLAST序列比对分析'),
#         ('coverage', 'coverage', '📊 BAM覆盖度分析'),
#         ('ena_downloader', 'ena-downloader', '📥 ENA数据下载工具'),
#         ('fastp', 'fastp', '🧹 FASTQ数据质量控制'),
#         ('genomesyn', 'genomesyn', '🗺️ 基因组共线性分析'),
#         ('geneinfo', 'geneinfo', '📄 从GFF文件提取基因信息'),
#         ('gtx', 'gtx', '🔬 运行GTX WGS流程'),
#         ('hifiasm', 'hifiasm', '🧩 运行hifiasm基因组组装'),
#         ('kaks', 'kaks', '🧮 Ka/Ks计算'),
#         ('kmer_count', 'kmer-count', '🔢 K-mer丰度矩阵计算'),
#         ('kmer_query', 'kmer-extractor', '✂️ K-mer提取'),
#         ('longestmrna', 'longest-mrna', '📜 提取最长转录本'),
#         ('minimap2', 'minimap2', '🔗 Minimap2比对与区域提取'),
#         ('plinkgwas', 'plink-gwas', '📈 PLINK GWAS分析'),
#         ('popgen', 'popgen', '🌍 群体遗传学多样性分析'),
#         ('rnaseq', 'rnaseq', 'RNA-seq表达定量流程'),
#         ('split_fasta_id', 'split-fasta-id', '🔪 分割FASTA文件ID'),
#         ('vcf_filter', 'vcf-filter', '🩸 VCF文件筛选'),
#         ('vcf_genotype', 'vcf-genotype', '🔬 VCF基因型提取'),
#         ('vcf_pca', 'vcf-pca', '📊 VCF主成分分析 (PCA)'),
#         ('vcf_nj_tree', 'vcf-nj-tree', '🌳 VCF构建NJ进化树'),
#         ('vcf_sample_hete', 'vcf-sample-hete', '📈 VCF样本基因型统计'),
#         ('vcf_sequence', 'vcf-sequence', '🧬 从基因组和VCF提取序列'),
#         ('bismark', 'bismark', '🧬 全基因组甲基化'),
#         ('transcriptome_prediction', 'mrna-prediction', '🧬 基于转录组的基因预测'),
#         ('parabricks','parabricks','机遇GPU的全基因组流程'),
#         ('raxml', 'raxml', 'RAxML系统发育树'),
#         ('vcf2phylip','vcf2phylip','vcf转phylip格式'),
#         ('repeat_analyzer','repeat-analyzer','重复序列分析模块'),
#         ('edta','edta','EDTA模块'),
#         ('genomethreader','genome-threader','genomethreader预测基因结构')
#     ]
    
#     for module_name, command_name, description in command_modules:
#         try:
#             module = __import__(f'biopytools.cli.commands.{module_name}', 
#                               fromlist=[module_name])
#             command_func = getattr(module, module_name)
#             available_commands[command_name] = (command_func, description)
#         except (ImportError, AttributeError) as e:
#             # --- 关键修改在这里 ---
#             # 不再沉默！打印出错误警告，帮助我们调试。
#             click.secho(f" [!] 警告: 无法加载命令 '{command_name}'. 错误: {e}", fg='yellow', err=True)
#             pass
    
#     return available_commands


# # --- 后续代码保持不变 ---

# # @click.group(invoke_without_command=True)
# # @click.version_option(version=__version__, prog_name='biopytools')
# # @click.pass_context
# # def cli(ctx):
# #     """
# #     BioPyTools - 生物信息学分析工具包
# #     ... (文档字符串保持不变) ...
# #     """
# #     if ctx.invoked_subcommand is None:
# #         click.echo(ctx.get_help())
# #         # 您的自定义帮助显示逻辑保持不变
# #         available_commands = import_available_commands()
# #         if available_commands:
# #             click.echo(f"\n当前可用的命令 | Currently available commands:")
# #             for cmd_name in sorted(available_commands.keys()):
# #                 _, description = available_commands[cmd_name]
# #                 click.echo(f"  {cmd_name:<20} {description}")
# #         # ... (您的“即将推出”部分也保持不变) ...

# @click.group(
#     # 使用 context_settings 统一配置
#     context_settings=dict(
#         # 显式定义触发帮助文档的选项，确保 -h 和 --help 都可用
#         help_option_names=['-h', '--help'],
#         # 设置帮助文档的最大内容宽度，防止不必要的换行
#         max_content_width=120  
#     ),
#     invoke_without_command=True
# )
# # @click.version_option(version=__version__, prog_name='biopytools')
# @click.version_option(__version__, '-v', '--version', prog_name='biopytools', message='%(prog)s, version %(version)s')
# @click.pass_context
# def cli(ctx):
#     """
#     BioPyTools - 生物信息学分析工具包

    
#     要查看特定命令的帮助，请运行：biopytools <命令> -h/--help, 如biopytools fastp -h
#     """
#     # 这段代码只在用户直接运行 `biopytools` (不带任何子命令或-h) 时执行
#     if ctx.invoked_subcommand is None:
#         # 显示click生成的标准帮助信息
#         click.echo(ctx.get_help())
        
#         # 附加你自定义的命令列表，使其更美观
#         available_commands = import_available_commands()
#         if available_commands:
#             click.echo(f"\n当前可用的命令 | Currently available commands:")
#             # 计算最长命令名的长度以对齐
#             max_len = max(len(name) for name in available_commands.keys()) if available_commands else 20
            
#             for cmd_name in sorted(available_commands.keys()):
#                 _, description = available_commands[cmd_name]
#                 # 使用 f-string 进行格式化对齐
#                 click.echo(f"  {cmd_name:<{max_len + 2}} {description}")
        
#         # click.echo("\n即将推出 | Coming soon:")
#         # click.echo("  vcf-anno             - VCF变异注释")
#         # click.echo("  gsea                 - 基因富集分析")
#         # click.echo("\n项目地址 | Project URL: https://github.com/your-repo")

# def register_commands():
#     """注册所有可用的子命令"""
#     available_commands = import_available_commands()
#     for command_name, (command_func, description) in available_commands.items():
#         # 这里我们恢复使用 docstring 的方式，因为那是您最初正确的方式
#         # 如果您想强制使用列表里的描述，可以取消下面一行的注释
#         # command_func.short_help = description 
#         cli.add_command(command_func, name=command_name)

# def main():
#     """主入口函数"""
#     register_commands()
#     cli()

# if __name__ == '__main__':
#     main()

#!/usr/bin/env python3
"""
BioPyTools 统一CLI入口点 | BioPyTools Unified CLI Entry Point
"""

import click
from .._version import __version__

# 硬编码所有命令信息，用于快速显示帮助
COMMAND_REGISTRY = [
    # (模块文件名, 命令名, 描述文本)
    ('admixture', 'admixture', '🧬 ADMIXTURE群体结构分析'),
    ('annovar', 'annovar', '📝 ANNOVAR变异注释'),
    ('blast', 'blast', '🧬 BLAST序列比对分析'),
    ('coverage', 'coverage', '📊 BAM覆盖度分析'),
    ('ena_downloader', 'ena-downloader', '📥 ENA数据下载工具'),
    ('fastp', 'fastp', '🧹 FASTQ数据质量控制'),
    ('genomesyn', 'genomesyn', '🗺️  基因组共线性分析'),
    ('geneinfo', 'geneinfo', '📄 从GFF文件提取基因信息'),
    ('gtx', 'gtx', '🔬 运行GTX WGS流程'),
    ('hifiasm', 'hifiasm', '🧩 运行hifiasm基因组组装'),
    ('kaks', 'kaks', '🧮 Ka/Ks计算'),
    ('kmer_count', 'kmer-count', '🔢 K-mer丰度矩阵计算'),
    ('kmer_query', 'kmer-extractor', '✂️  K-mer提取'),
    ('longestmrna', 'longest-mrna', '📜 提取最长转录本'),
    ('minimap2', 'minimap2', '🔗 Minimap2比对与区域提取'),
    ('plinkgwas', 'plink-gwas', '📈 PLINK GWAS分析'),
    ('popgen', 'popgen', '🌍 群体遗传学多样性分析'),
    ('rnaseq', 'rnaseq', '🧬 RNA-seq表达定量流程'),
    ('split_fasta_id', 'split-fasta-id', '🔪 分割FASTA文件ID'),
    ('vcf_filter', 'vcf-filter', '🩸 VCF文件筛选'),
    ('vcf_genotype', 'vcf-genotype', '🔬 VCF基因型提取'),
    ('vcf_pca', 'vcf-pca', '📊 VCF主成分分析 (PCA)'),
    ('vcf_nj_tree', 'vcf-nj-tree', '🌳 VCF构建NJ进化树'),
    ('vcf_sample_hete', 'vcf-sample-hete', '📈 VCF样本基因型统计'),
    ('vcf_sequence', 'vcf-sequence', '🧬 从基因组和VCF提取序列'),
    ('bismark', 'bismark', '🧬 全基因组甲基化'),
    ('transcriptome_prediction', 'mrna-prediction', '🧬 基于转录组的基因预测'),
    ('parabricks', 'parabricks', '🧬 基于GPU的全基因组流程'),
    ('raxml', 'raxml', '🌳 RAxML系统发育树'),
    ('vcf2phylip', 'vcf2phylip', '🔄 vcf转phylip格式'),
    ('repeat_analyzer', 'repeat-analyzer', '🔄 重复序列分析模块'),
    ('edta', 'edta', '🧬 EDTA重复元件注释'),
    ('genomethreader', 'genome-threader', '🔬 GenomeThreader预测基因结构'),
    ('orthofinder', 'orthofinder', '🧬 OrthoFinder泛基因组分析工具包'),
    ('genomeasm', 'genomeasm', '🧬 三代基因组组装流程'),
    ('gffconverter', 'renamegff', '✂️  GFF文件整理工具'),
    ('indelpav', 'indelpav', '🧬 INDEL PAV分析工具'),
    ('busco', 'busco', '🧬 BUSCO质量评估分析工具'),
    ('genebank2fasta', 'genebank2fasta','🧬 GenBank序列提取工具'),
    ('parse_seq', 'parse-seq','🧬 核酸或蛋白序列提取工具'),
    ('parse_gene_dna', 'parse-gene-dna','🧬 基因DNA序列提取工具'),
    ('bwa', 'bwa','🧬 全基因组比对工具'),
    ('mafft_fasttree', 'mafft-fasttree','🌳 系统发育树构建工具'),
    ('bwa_gatk', 'bwa-gatk','🧬 全基因组比对和编译检测工具'),
    ('iqtree', 'iqtree','🌲 IQ-TREE系统发育树分析工具'),
    ('msa', 'msa','🧬 多序列比对分析工具'),
    ('sra2fastq', 'sra2fastq','🧬 SRA转FASTQ转换工具')
]

# 将硬编码信息转换为字典，方便查询
COMMAND_INFO = {cmd_name: description for _, cmd_name, description in COMMAND_REGISTRY}

class LazyGroup(click.Group):
    """懒加载组类 - 只在需要时才导入命令模块"""
    
    def get_command(self, ctx, cmd_name):
        """获取命令时才导入对应模块"""
        # 查找对应的模块名
        module_name = None
        for mod_name, command_name, _ in COMMAND_REGISTRY:
            if command_name == cmd_name:
                module_name = mod_name
                break
        
        if module_name is None:
            return None
            
        try:
            module = __import__(f'biopytools.cli.commands.{module_name}', 
                              fromlist=[module_name])
            command_func = getattr(module, module_name)
            return command_func
        except (ImportError, AttributeError) as e:
            click.secho(f" [!] 错误: 无法加载命令 '{cmd_name}'. 错误: {e}", fg='red', err=True)
            return None
    
    def list_commands(self, ctx):
        """返回所有可用命令列表"""
        return [cmd_name for _, cmd_name, _ in COMMAND_REGISTRY]
    
    def format_commands(self, ctx, formatter):
        """自定义命令列表格式化，使用硬编码的emoji描述"""
        commands = []
        for _, cmd_name, description in sorted(COMMAND_REGISTRY, key=lambda x: x[1]):
            commands.append((cmd_name, description))
        
        if commands:
            with formatter.section("Commands"):
                formatter.write_dl(commands)

@click.group(
    cls=LazyGroup,
    context_settings=dict(
        help_option_names=['-h', '--help'],
        max_content_width=120  
    ),
    invoke_without_command=True
)
@click.version_option(__version__, '-v', '--version', prog_name='biopytools', message='%(prog)s, version %(version)s')
@click.pass_context
def cli(ctx):
    """
    BioPyTools - 生物信息学分析工具包

    
    要查看特定命令的帮助，请运行：biopytools <命令> -h/--help, 如biopytools fastp -h
    """
    if ctx.invoked_subcommand is None:
        # 显示click生成的标准帮助信息（包含我们自定义的Commands部分）
        click.echo(ctx.get_help())

def main():
    """主入口函数"""
    cli()

if __name__ == '__main__':
    main()