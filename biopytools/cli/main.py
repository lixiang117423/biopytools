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
    ('bam_stats', 'bam-stats', '📊 BAM文件批量统计分析'),
    ('bam_cov', 'bam-cov', '📊 BAM覆盖度统计'),
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
    ('sra2fastq', 'sra2fastq','🧬 SRA转FASTQ转换工具'),
    ('filter_snp_indel', 'filter-snp-indel','🧬 SNP和INDEL过滤工具'),
    ('gatk_joint', 'gatk-joint','🧬 GATK Joint Genotyping工具'),
    ('filter_annovar', 'filter-annovar','🧬 基因区域变异提取工具'),
    ('metagraph_kmer', 'metagraph-kmer','🧬 K-mer库构建与查询分析工具'),
    ('fastq2vcf_parabricks', 'fastq2vcf-parabricks','🧬 Fastq到VCF (Parabricks) 全流程分析'),
    ('fastq2vcf_gtx', 'fastq2vcf-gtx','🧬 Fastq到VCF (GTX) 全流程分析'),
    ('hifi_hic', 'hifi-hic',"🧬 使用HiFi和Hi-C数据进行基因组组装"),
    ('allhic', 'allhic',"'🧬 使用ALLHiC进行染色体挂载'"),
    ('get_link_from_CNCB', 'get-link-from-CNCB','📥 从CNCB批量获取测序数据下载链接'),
    ('gwas_lambda', 'gwas-lambda', '📊 GWAS Lambda GC计算工具'),
    ('mcyc', 'mcyc', '🧬 甲烷循环基因丰度分析工具'),
    ('haphic', 'haphic', '🧬 HapHiC基因组scaffolding工具'),
    ('subseq', 'subseq', '🧬 序列子集提取工具'),
    ('tassel_gwas', 'tassel-gwas', '🌾 TASSEL GWAS分析工具'),
    ('snp_index', 'snp-index', '🧬 SNP index计算和分析工具'),
    ('gtx_joint', 'gtx-joint', '🧬 GTX Joint Calling命令生成工具'),
    ('vcf_renamer', 'vcf-renamer', '🏷️ VCF样品名称重命名工具'),
    ('egapx_batch', 'egapx-batch', '🧬 EGAPx批量运行配置生成工具'),
    ('genome_analysis', 'genome-analysis', '🧬 GenomeScope2基因组评估工具'),
    ('rename_chromosomes', 'rename-chromosomes', '🧬 染色体重命名工具'),
    ('dsuite', 'dsuite', '🧬 Dsuite D统计量分析工具'),
    ('ngenomesyn', 'ngenomesyn', '🧬 NGenomeSyn可视化工具'),
    ('assembly2agp', 'assembly2agp', '🧬 Assembly文件转AGP格式工具'),
    ('gemma_gwas', 'gemma-gwas', '🧬 GEMMA GWAS批量分析工具'),
    ('vcf_merger', 'vcf-merger', '🧬 VCF按染色体合并工具'),
    ('hicanu', 'hicanu', '🧬 HiCanu基因组组装工具'),
    ('find_telomere', 'find-telomere', '🧬 端粒识别分析工具'),
    ('get_plastome', 'get-plastome', '🧬 叶绿体基因组组装工具'),
    ('vcf_sampler', 'vcf-sampler', '🎲 VCF文件SNP抽样工具')
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