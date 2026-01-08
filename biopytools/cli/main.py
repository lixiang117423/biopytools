#!/usr/bin/env python3
"""
BioPyTools 统一CLI入口点|BioPyTools Unified CLI Entry Point
"""

import click
from .._version import __version__

# 硬编码所有命令信息，用于快速显示帮助
COMMAND_REGISTRY = [
    # (模块文件名, 命令名, 描述文本)
    ('admixture', 'admixture', 'ADMIXTURE群体结构分析|ADMIXTURE Population Structure Analysis'),
    ('annovar', 'annovar', 'ANNOVAR变异注释|ANNOVAR Variant Annotation'),
    ('allhic', 'allhic', '使用ALLHiC进行染色体挂载|Use ALLHiC for chromosome scaffolding'),
    ('assembly2agp', 'assembly2agp', 'Assembly文件转AGP格式工具|Assembly to AGP format converter'),
    ('assembly_stats', 'assembly-stats', '基因组装配统计|Genome Assembly Statistics'),
    ('assembly_qv', 'assembly-qv', '装配质量QV值计算|Assembly Quality QV Calculation'),
    ('bam_stats', 'bam-stats', 'BAM文件批量统计分析|BAM File Batch Statistics Analysis'),
    ('bam_cov', 'bam-cov', 'BAM覆盖度统计|BAM Coverage Statistics'),
    ('bam2fastq', 'bam2fastq', 'BAM to FASTQ批量转换工具|BAM to FASTQ batch conversion tool'),
    ('bismark', 'bismark', '全基因组甲基化分析|Whole genome methylation analysis'),
    ('blast', 'blast', 'BLAST序列比对分析|BLAST Sequence Alignment Analysis'),
    ('bwa', 'bwa', '全基因组比对工具|Whole genome alignment tool'),
    ('bwa_gatk', 'bwa-gatk', '全基因组比对和变异检测|Whole genome alignment and variant detection'),
    ('busco', 'busco', 'BUSCO质量评估分析工具|BUSCO quality assessment tool'),
    ('dsuite', 'dsuite', 'Dsuite D统计量分析工具|Dsuite D-statistics analysis tool'),
    ('edta', 'edta', 'EDTA重复元件注释|EDTA repeat element annotation'),
    ('egapx_batch', 'egapx-batch', 'EGAPx批量运行配置生成工具|EGAPx batch run configuration generator'),
    ('ena_downloader', 'ena-downloader', 'ENA数据下载工具|ENA data download tool'),
    ('fastp', 'fastp', 'FASTQ数据质量控制|FASTQ data quality control'),
    ('fastq2vcf_gtx', 'fastq2vcf-gtx', 'Fastq到VCF (GTX) 全流程分析|Fastq to VCF (GTX) pipeline'),
    ('fastq2vcf_parabricks', 'fastq2vcf-parabricks', 'Fastq到VCF (Parabricks) 全流程分析|Fastq to VCF (Parabricks) pipeline'),
    # ('filter_annovar', 'filter-annovar', '基因区域变异提取工具|Gene region variant extraction tool'),
    ('filter_snp_indel', 'filter-snp-indel', 'SNP和INDEL过滤工具|SNP and INDEL filtering tool'),
    ('find_telomere', 'find-telomere', '端粒识别分析工具|Telomere identification analysis tool'),
    ('gatk_joint', 'gatk-joint', 'GATK Joint Genotyping工具|GATK Joint Genotyping tool'),
    ('gemma_gwas', 'gemma-gwas', 'GEMMA GWAS批量分析工具|GEMMA GWAS batch analysis tool'),
    ('genebank2fasta', 'genebank2fasta', 'GenBank序列提取工具|GenBank sequence extraction tool'),
    ('geneinfo', 'geneinfo', '从GFF文件提取基因信息|Extract gene information from GFF file'),
    ('gff_renamer', 'gff-renamer', 'GFF文件ID规范化工具|GFF file ID standardization tool'),
    ('genome_analysis', 'genome-analysis', 'GenomeScope2基因组评估工具|GenomeScope2 genome evaluation tool'),
    # ('genome_threader', 'genome-threader', 'GenomeThreader预测基因结构|GenomeThreader gene structure prediction'),
    ('genomeasm', 'genomeasm', '三代基因组组装流程|Third generation genome assembly pipeline'),
    ('genomesyn', 'genomesyn', '基因组共线性分析|Genome collinearity analysis'),
    ('get_link_from_CNCB', 'get-link-from-CNCB', '从CNCB批量获取测序数据下载链接|Batch download links from CNCB'),
    ('get_plastome', 'get-plastome', '叶绿体基因组组装工具|Chloroplast genome assembly tool'),
    ('gffconverter', 'renamegff', 'GFF文件整理工具|GFF file organization tool'),
    ('gtx', 'gtx', '运行GTX WGS流程|Run GTX WGS pipeline'),
    ('gtx_joint', 'gtx-joint', 'GTX Joint Calling命令生成工具|GTX Joint Calling command generator'),
    ('gwas_lambda', 'gwas-lambda', 'GWAS Lambda GC计算工具|GWAS Lambda GC calculation tool'),
    ('haphic', 'haphic', 'HapHiC基因组scaffolding工具|HapHiC genome scaffolding tool'),
    ('hicanu', 'hicanu', 'HiCanu基因组组装工具|HiCanu genome assembly tool'),
    ('hifi_hic', 'hifi-hic', '使用HiFi和Hi-C数据进行基因组组装|Genome assembly using HiFi and Hi-C data'),
    ('hifiasm', 'hifiasm', '运行hifiasm基因组组装|Run hifiasm genome assembly'),
    ('indelpav', 'indelpav', 'INDEL PAV分析工具|INDEL PAV analysis tool'),
    ('iqtree', 'iqtree', 'IQ-TREE系统发育树分析工具|IQ-TREE phylogenetic tree analysis tool'),
    ('kaks', 'kaks', 'Ka/Ks计算|Ka/Ks calculation'),
    ('kmer_count', 'kmer-count', 'K-mer丰度矩阵计算|K-mer abundance matrix calculation'),
    ('kmertools', 'kmertools', 'K-mer工具集 - 建库、查询和分析|K-mer Tools - Build, Query and Analysis'),
    ('kmer_query', 'kmer-extractor', 'K-mer提取|K-mer extraction'),
    ('lai', 'lai', 'LAI组装质量指数计算工具|LAI Assembly Index calculator'),
    ('longestmrna', 'longest-mrna', '提取最长转录本|Extract longest transcript'),
    ('longrnaseq', 'longrnaseq', '三代转录组比对工具|Long RNA-seq alignment tool'),
    ('mafft_fasttree', 'mafft-fasttree', '系统发育树构建工具|Phylogenetic tree construction tool'),
    ('mcyc', 'mcyc', '甲烷循环基因丰度分析工具|Methane cycle gene abundance analysis tool'),
    # ('metagraph_kmer', 'metagraph-kmer', 'K-mer库构建与查询分析工具|K-mer library construction and query tool'),
    ('minimap2', 'minimap2', 'Minimap2比对与区域提取|Minimap2 alignment and region extraction'),
    ('msa', 'msa', '多序列比对分析工具|Multiple sequence alignment analysis tool'),
    ('ngenomesyn', 'ngenomesyn', 'NGenomeSyn可视化工具|NGenomeSyn visualization tool'),
    ('orthofinder', 'orthofinder', 'OrthoFinder泛基因组分析工具包|OrthoFinder pan-genome analysis toolkit'),
    ('parabricks', 'parabricks', '基于GPU的全基因组流程|GPU-based whole genome pipeline'),
    ('parse_gene_dna', 'parse-gene-dna', '基因DNA序列提取工具|Gene DNA sequence extraction tool'),
    ('parse_seq', 'parse-seq', '核酸或蛋白序列提取工具|Nucleotide or protein sequence extraction tool'),
    ('plinkgwas', 'plink-gwas', 'PLINK GWAS分析|PLINK GWAS analysis'),
    # ('popgen', 'popgen', '群体遗传学多样性分析|Population genetics diversity analysis'),
    ('raxml', 'raxml', 'RAxML系统发育树|RAxML phylogenetic tree'),
    ('rename_chromosomes', 'rename-chromosomes', '染色体重命名工具|Chromosome renaming tool'),
    ('repeat_analyzer', 'repeat-analyzer', '重复序列分析模块|Repeat sequence analysis module'),
    ('rnaseq', 'rnaseq', 'RNA-seq表达定量流程|RNA-seq expression quantification pipeline'),
    ('snp_index', 'snp-index', 'SNP index计算和分析工具|SNP index calculation and analysis tool'),
    ('sra2fastq', 'sra2fastq', 'SRA转FASTQ转换工具|SRA to FASTQ conversion tool'),
    ('split_fasta_id', 'split-fasta-id', '分割FASTA文件ID|Split FASTA file ID'),
    ('subseq', 'subseq', '序列子集提取工具|Sequence subset extraction tool'),
    ('tassel_gwas', 'tassel-gwas', 'TASSEL GWAS分析工具|TASSEL GWAS analysis tool'),
    ('vcf2phylip', 'vcf2phylip', 'VCF转phylip格式|VCF to phylip format conversion'),
    ('vcf_filter', 'vcf-filter', 'VCF文件筛选|VCF file filtering'),
    ('vcf_genotype', 'vcf-genotype', 'VCF基因型提取|VCF genotype extraction'),
    ('vcf_merger', 'vcf-merger', 'VCF按染色体合并工具|VCF chromosome merge tool'),
    ('vcf_nj_tree', 'vcf-nj-tree', 'VCF构建NJ进化树|VCF NJ phylogenetic tree construction'),
    ('vcf_pca', 'vcf-pca', 'VCF主成分分析 (PCA)|VCF Principal Component Analysis (PCA)'),
    ('vcf_renamer', 'vcf-renamer', 'VCF样品名称重命名工具|VCF sample name renaming tool'),
    ('vcf_sampler', 'vcf-sampler', 'VCF文件SNP抽样工具|VCF file SNP sampling tool'),
    ('vcf_sample_hete', 'vcf-sample-hete', 'VCF样本基因型统计|VCF sample genotype statistics'),
    ('vcf_sequence', 'vcf-sequence', '从基因组和VCF提取序列|Extract sequences from genome and VCF')
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
        """自定义命令列表格式化，使用分隔符连接命令和描述"""
        # 计算最长命令名长度，用于对齐
        max_cmd_len = max(len(cmd_name) for _, cmd_name, _ in COMMAND_REGISTRY)

        with formatter.section("Commands"):
            for _, cmd_name, description in sorted(COMMAND_REGISTRY, key=lambda x: x[1]):
                # 创建带分隔符的格式：command -------- description
                separator = '-' * (max_cmd_len - len(cmd_name) + 4)
                line = f"  {cmd_name} {separator} {description}"
                formatter.write_text(line)

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