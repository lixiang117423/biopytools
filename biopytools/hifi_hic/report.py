"""
组装结果报告模块 |Assembly Results Report Module
"""

import os
from pathlib import Path
from datetime import datetime

try:
    from utils import get_fasta_stats, format_time
except ImportError:
    from biopytools.genome_assembler.utils import get_fasta_stats, format_time

class ReportGenerator:
    """报告生成器|Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, fasta_results: dict, start_time: datetime, end_time: datetime):
        """生成组装摘要报告|Generate assembly summary report"""
        report_file = os.path.join(self.config.stat_dir, f"{self.config.prefix}_assembly_statistics.txt")
        
        duration = int((end_time - start_time).total_seconds())
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write(" 基因组组装统计报告|Genome Assembly Statistics Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f" 样本名称|Sample Name: {self.config.prefix}\n")
            f.write(f" 报告日期|Report Date: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f" 总耗时|Total Duration: {format_time(duration)}\n\n")
            
            f.write(f"  组装参数|Assembly Parameters:\n")
            f.write(f"  - 基因组大小|Genome Size: {self.config.genome_size}\n")
            f.write(f"  - 倍性|Ploidy: {self.config.n_hap}\n")
            f.write(f"  - 线程数|Threads: {self.config.threads}\n\n")
            
            f.write(f" 输入文件|Input Files:\n")
            f.write(f"  - HiFi数据|HiFi Data: {os.path.basename(self.config.hifi_data)}\n")
            if self.config.has_hic:
                f.write(f"  - Hi-C R1|Hi-C R1: {os.path.basename(self.config.hic_r1)}\n")
                f.write(f"  - Hi-C R2|Hi-C R2: {os.path.basename(self.config.hic_r2)}\n")
            else:
                f.write(f"  - Hi-C数据|Hi-C Data: 未使用|Not used\n")
            f.write("\n")
            
            f.write("=" * 50 + "\n")
            f.write("组装结果|Assembly Results\n")
            f.write("=" * 50 + "\n\n")
            
            for suffix in ["hap1", "hap2", "primary", "alternate"]:
                fasta_name = f"{self.config.prefix}.{suffix}.fa"
                if fasta_name in fasta_results:
                    stats = fasta_results[fasta_name]
                    if stats:
                        f.write(f" {fasta_name}\n")
                        f.write("-" * 40 + "\n")
                        
                        total_gb = stats['total_len'] / 1_000_000_000
                        f.write(f"序列数量|Number of Sequences: {stats['num_seqs']:,}\n")
                        f.write(f"总长度|Total Length: {stats['total_len']:,} bp ({total_gb:.2f} Gb)\n")
                        
                        # 获取完整统计信息|Get full statistics
                        fasta_path = os.path.join(self.config.fasta_dir, fasta_name)
                        full_stats = get_fasta_stats(fasta_path)
                        
                        if 'max_len' in full_stats:
                            f.write(f"最长序列|Longest Sequence: {full_stats['max_len']:,} bp\n")
                        
                        if 'n50' in full_stats:
                            f.write(f"N50: {full_stats['n50']:,} bp\n")
                        
                        f.write("\n")
            
            f.write("=" * 50 + "\n")
            f.write("文件位置|File Locations\n")
            f.write("=" * 50 + "\n")
            f.write(f"原始输出|Raw Output: {self.config.raw_dir}\n")
            f.write(f"FASTA文件|FASTA Files: {self.config.fasta_dir}\n")
            f.write(f"日志文件|Log Files: {self.config.log_dir}\n")
            f.write(f"统计文件|Statistics: {self.config.stat_dir}\n")
            f.write("=" * 80 + "\n")
        
        self.logger.info(f" 统计报告已生成|Statistics report generated: {report_file}")
    
    def print_file_tree(self):
        """打印输出文件结构|Print output file structure"""
        self.logger.info("输出文件结构|Output File Structure:")

        if self.config.has_ngs:
            # 有NGS数据的目录结构|Directory structure with NGS data
            if self.config.has_hic:
                gfa_suffix = "hic"
                tree = f"""
 {self.config.work_dir}/
   │
   ├──  01.raw_output/     (初次组装原始输出|Initial assembly raw output)
   │   ├──  {self.config.prefix}.hic.hap1.p_ctg.gfa
   │   ├──  {self.config.prefix}.hic.hap2.p_ctg.gfa
   │   ├──  {self.config.prefix}.hic.p_ctg.gfa
   │   └──  其他中间文件|Other intermediate files
   │
   ├──  02.fasta/          (初次组装FASTA|Initial assembly FASTA)
   │   ├──  {self.config.prefix}.hap1.fa
   │   ├──  {self.config.prefix}.hap2.fa
   │   ├──  {self.config.prefix}.primary.fa
   │   ├──  {self.config.prefix}.alternate.fa
   │   ├──  {self.config.prefix}.p_ctg.contig_reads.tsv (contig-reads映射|mapping)
   │   ├──  {self.config.prefix}.hap1.p_ctg.contig_reads.tsv
   │   └──  {self.config.prefix}.hap2.p_ctg.contig_reads.tsv
   │
   ├──  03.ngs_polish/     (NGS筛选和重新组装|NGS filtering & reassembly)
   │   ├──  01.bwa_alignment/              (BWA比对结果|BWA alignment)
   │   │   └──  bam/{self.config.prefix}.bam
   │   ├──  02.coverage_filter/            (覆盖度过滤结果|Coverage filter)
   │   │   └──  {self.config.prefix}_high_quality.list
   │   ├──  03.filtered_reads/             (筛选的reads|Filtered reads)
   │   │   └──  {self.config.prefix}_high_quality_reads.fq.gz
   │   ├──  04.reassembly/                 (重新组装结果|Reassembly)
   │   │   ├──  01.raw_output/
   │   │   └──  02.fasta/
   │   └──  {self.config.prefix}.polished.fa (最终polished基因组|Final polished genome)
   │
   ├──  04.statistics/     (统计信息|Statistics)
   │   └──  {self.config.prefix}_assembly_statistics.txt
   │
   └──  05.logs/           (日志文件|Log files)
       └──  assembly.log
            """
            else:
                tree = f"""
 {self.config.work_dir}/
   │
   ├──  01.raw_output/     (初次组装原始输出|Initial assembly raw output)
   │   ├──  {self.config.prefix}.hap1.p_ctg.gfa
   │   ├──  {self.config.prefix}.hap2.p_ctg.gfa
   │   ├──  {self.config.prefix}.p_ctg.gfa
   │   └──  其他中间文件|Other intermediate files
   │
   ├──  02.fasta/          (初次组装FASTA|Initial assembly FASTA)
   │   ├──  {self.config.prefix}.hap1.fa
   │   ├──  {self.config.prefix}.hap2.fa
   │   ├──  {self.config.prefix}.primary.fa
   │   ├──  {self.config.prefix}.alternate.fa
   │   ├──  {self.config.prefix}.p_ctg.contig_reads.tsv (contig-reads映射|mapping)
   │   ├──  {self.config.prefix}.hap1.p_ctg.contig_reads.tsv
   │   └──  {self.config.prefix}.hap2.p_ctg.contig_reads.tsv
   │
   ├──  03.ngs_polish/     (NGS筛选和重新组装|NGS filtering & reassembly)
   │   ├──  01.bwa_alignment/              (BWA比对结果|BWA alignment)
   │   │   └──  bam/{self.config.prefix}.bam
   │   ├──  02.coverage_filter/            (覆盖度过滤结果|Coverage filter)
   │   │   └──  {self.config.prefix}_high_quality.list
   │   ├──  03.filtered_reads/             (筛选的reads|Filtered reads)
   │   │   └──  {self.config.prefix}_high_quality_reads.fq.gz
   │   ├──  04.reassembly/                 (重新组装结果|Reassembly)
   │   │   ├──  01.raw_output/
   │   │   └──  02.fasta/
   │   └──  {self.config.prefix}.polished.fa (最终polished基因组|Final polished genome)
   │
   ├──  04.statistics/     (统计信息|Statistics)
   │   └──  {self.config.prefix}_assembly_statistics.txt
   │
   └──  05.logs/           (日志文件|Log files)
       └──  assembly.log
            """
        else:
            # 无NGS数据的目录结构|Directory structure without NGS data
            if self.config.has_hic:
                tree = f"""
 {self.config.work_dir}/
   │
   ├──  01.raw_output/     (原始输出文件|Raw output files)
   │   ├──  {self.config.prefix}.hic.hap1.p_ctg.gfa
   │   ├──  {self.config.prefix}.hic.hap2.p_ctg.gfa
   │   ├──  {self.config.prefix}.hic.p_ctg.gfa
   │   ├──  {self.config.prefix}.hic.a_ctg.gfa (如果存在|if present)
   │   └──  其他中间文件|Other intermediate files
   │
   ├──  02.fasta/          (FASTA格式|FASTA format)
   │   ├──  {self.config.prefix}.hap1.fa
   │   ├──  {self.config.prefix}.hap2.fa
   │   ├──  {self.config.prefix}.primary.fa
   │   ├──  {self.config.prefix}.alternate.fa
   │   ├──  {self.config.prefix}.p_ctg.contig_reads.tsv (contig-reads映射|contig-reads mapping)
   │   ├──  {self.config.prefix}.hap1.p_ctg.contig_reads.tsv
   │   └──  {self.config.prefix}.hap2.p_ctg.contig_reads.tsv
   │
   ├──  03.statistics/     (统计信息|Statistics)
   │   └──  {self.config.prefix}_assembly_statistics.txt
   │
   └──  04.logs/           (日志文件|Log files)
       └──  assembly.log
            """
            else:
                tree = f"""
 {self.config.work_dir}/
   │
   ├──  01.raw_output/     (原始输出文件|Raw output files)
   │   ├──  {self.config.prefix}.hap1.p_ctg.gfa
   │   ├──  {self.config.prefix}.hap2.p_ctg.gfa
   │   ├──  {self.config.prefix}.p_ctg.gfa
   │   ├──  {self.config.prefix}.a_ctg.gfa (如果存在|if present)
   │   └──  其他中间文件|Other intermediate files
   │
   ├──  02.fasta/          (FASTA格式|FASTA format)
   │   ├──  {self.config.prefix}.hap1.fa
   │   ├──  {self.config.prefix}.hap2.fa
   │   ├──  {self.config.prefix}.primary.fa
   │   ├──  {self.config.prefix}.alternate.fa
   │   ├──  {self.config.prefix}.p_ctg.contig_reads.tsv (contig-reads映射|contig-reads mapping)
   │   ├──  {self.config.prefix}.hap1.p_ctg.contig_reads.tsv
   │   └──  {self.config.prefix}.hap2.p_ctg.contig_reads.tsv
   │
   ├──  03.statistics/     (统计信息|Statistics)
   │   └──  {self.config.prefix}_assembly_statistics.txt
   │
   └──  04.logs/           (日志文件|Log files)
       └──  assembly.log
            """

        self.logger.info(tree)

        self.logger.info(" 下一步建议|Next Steps:")
        if self.config.has_ngs:
            self.logger.info("   1. 查看polished基因组|1. View polished genome")
            self.logger.info(f"      {self.config.ngs_polish_dir}/{self.config.prefix}.polished.fa")
        self.logger.info("   2. 查看统计报告|2. View statistics report")
        self.logger.info("   3. 检查组装质量|3. Check assembly quality (QUAST/BUSCO)")
        self.logger.info("   4. 可视化组装|4. Visualize assembly (Bandage)")
        self.logger.info("   5. 查看contig-reads映射|5. View contig-reads mapping")
