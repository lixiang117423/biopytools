"""
组装结果报告模块 |Assembly Results Report Module
"""

import os
from pathlib import Path
from datetime import datetime

from .utils import get_fasta_stats, format_time

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
            
            # 根据n_hap动态构建后缀列表|Build suffix list dynamically based on n_hap
            suffixes = [f"hap{i}" for i in range(1, self.config.n_hap + 1)] + ["primary", "alternate"]

            for suffix in suffixes:
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
    
    def _build_file_tree(self) -> str:
        """
        根据配置动态构建文件树字符串|Build file tree string dynamically based on config

        Returns:
            str: 文件树字符串|File tree string
        """
        prefix = self.config.prefix
        hic = ".hic" if self.config.has_hic else ""
        n_hap = self.config.n_hap

        # 构建GFA文件列表|Build GFA file list
        gfa_lines = []
        for i in range(1, n_hap + 1):
            gfa_lines.append(f"│   ├──  {prefix}{hic}.hap{i}.p_ctg.gfa")
        gfa_lines.append(f"│   ├──  {prefix}{hic}.p_ctg.gfa")
        if not self.config.has_ngs:
            gfa_lines.append(f"│   ├──  {prefix}{hic}.a_ctg.gfa (如果存在|if present)")
        gfa_lines.append("│   └──  其他中间文件|Other intermediate files")
        gfa_block = "\n".join(gfa_lines)

        # 构建FASTA文件列表|Build FASTA file list
        fasta_lines = []
        for i in range(1, n_hap + 1):
            fasta_lines.append(f"│   ├──  {prefix}.hap{i}.fa")
        fasta_lines.append(f"│   ├──  {prefix}.primary.fa")
        fasta_lines.append(f"│   ├──  {prefix}.alternate.fa")
        fasta_lines.append(f"│   ├──  {prefix}.p_ctg.contig_reads.tsv (contig-reads映射|mapping)")
        for i in range(1, n_hap + 1):
            fasta_lines.append(f"│   ├──  {prefix}.hap{i}.p_ctg.contig_reads.tsv")
        # 移除最后一个├──改为└──|Replace last ├── with └──
        fasta_lines[-1] = fasta_lines[-1].replace("├──", "└──")
        fasta_block = "\n".join(fasta_lines)

        if self.config.has_ngs:
            ngs_block = f"""
   ├──  03.ngs_polish/     (NGS筛选和重新组装|NGS filtering & reassembly)
   │   ├──  01.bwa_alignment/              (BWA比对结果|BWA alignment)
   │   │   └──  bam/{prefix}.bam
   │   ├──  02.coverage_filter/            (覆盖度过滤结果|Coverage filter)
   │   │   └──  {prefix}_high_quality.list
   │   ├──  03.filtered_reads/             (筛选的reads|Filtered reads)
   │   │   └──  {prefix}_high_quality_reads.fq.gz
   │   ├──  04.reassembly/                 (重新组装结果|Reassembly)
   │   │   ├──  01.raw_output/
   │   │   └──  02.fasta/
   │   └──  {prefix}.polished.fa (最终polished基因组|Final polished genome)
   │
   ├──  04.statistics/     (统计信息|Statistics)
   │   └──  {prefix}_assembly_statistics.txt
   │
   └──  05.logs/           (日志文件|Log files)
       └──  assembly.log"""
        else:
            ngs_block = f"""
   ├──  03.statistics/     (统计信息|Statistics)
   │   └──  {prefix}_assembly_statistics.txt
   │
   └──  04.logs/           (日志文件|Log files)
       └──  assembly.log"""

        tree = f"""
 {self.config.work_dir}/
   │
   ├──  01.raw_output/     (初次组装原始输出|Initial assembly raw output)
{gfa_block}
   │
   ├──  02.fasta/          (初次组装FASTA|Initial assembly FASTA)
{fasta_block}
   │
{ngs_block}
            """
        return tree

    def print_file_tree(self):
        """打印输出文件结构|Print output file structure"""
        self.logger.info("输出文件结构|Output File Structure:")

        tree = self._build_file_tree()

        self.logger.info(tree)

        self.logger.info(" 下一步建议|Next Steps:")
        if self.config.has_ngs:
            self.logger.info("   1. 查看polished基因组|1. View polished genome")
            self.logger.info(f"      {self.config.ngs_polish_dir}/{self.config.prefix}.polished.fa")
        self.logger.info("   2. 查看统计报告|2. View statistics report")
        self.logger.info("   3. 检查组装质量|3. Check assembly quality (QUAST/BUSCO)")
        self.logger.info("   4. 可视化组装|4. Visualize assembly (Bandage)")
        self.logger.info("   5. 查看contig-reads映射|5. View contig-reads mapping")
