"""
📊 统计信息生成模块 | Statistics Generation Module
"""

import subprocess
from pathlib import Path
from typing import List, Dict

class StatisticsGenerator:
    """统计信息生成器 | Statistics Generator"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def generate_statistics(self, samples: List[dict], bam_files: List[Path], 
                           filtered_vcfs: Dict[str, Path]):
        """生成统计信息 | Generate statistics"""
        self.logger.info("=" * 80)
        self.logger.info("📊 生成统计信息 | Generating statistics")
        self.logger.info("=" * 80)
        
        stats_file = self.config.stats_dir / "pipeline_statistics.txt"
        
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("🧬 BWA-GATK变异检测流程统计报告\n")
            f.write("BWA-GATK Variant Calling Pipeline Statistics Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 1. 样本信息 | Sample information
            f.write("📋 样本信息 | Sample Information\n")
            f.write("-" * 80 + "\n")
            f.write(f"总样本数 | Total samples: {len(samples)}\n")
            f.write(f"测序类型 | Sequencing type: {samples[0]['seq_type']}\n")
            f.write(f"倍性 | Ploidy: {self.config.ploidy}\n")
            if self.config.intervals:
                f.write(f"目标区间 | Target intervals: {Path(self.config.intervals).name}\n")
            f.write("\n")
            
            # 2. 比对统计 | Alignment statistics
            f.write("🎯 比对统计 | Alignment Statistics\n")
            f.write("-" * 80 + "\n")
            
            for sample, bam_file in zip(samples, bam_files):
                stats = self._get_alignment_stats(bam_file)
                f.write(f"\n样本 | Sample: {sample['name']}\n")
                f.write(f"  总reads数 | Total reads: {stats.get('total_reads', 'N/A')}\n")
                f.write(f"  比对reads数 | Mapped reads: {stats.get('mapped_reads', 'N/A')}\n")
                f.write(f"  比对率 | Mapping rate: {stats.get('mapping_rate', 'N/A')}\n")
                f.write(f"  正确配对率 | Properly paired: {stats.get('properly_paired', 'N/A')}\n")
                f.write(f"  重复率 | Duplicate rate: {stats.get('duplicate_rate', 'N/A')}\n")
            
            f.write("\n")
            
            # 3. 变异统计 | Variant statistics
            f.write("🔬 变异统计 | Variant Statistics\n")
            f.write("-" * 80 + "\n")
            
            # 原始变异数 | Raw variants
            raw_vcf = self.config.vcf_dir / "raw_variants.vcf.gz"
            if raw_vcf.exists():
                raw_stats = self._count_variants(raw_vcf)
                f.write(f"\n原始变异 | Raw variants:\n")
                f.write(f"  总变异数 | Total variants: {raw_stats['total']}\n")
                f.write(f"  SNP数量 | SNPs: {raw_stats.get('SNP', 0)}\n")
                f.write(f"  InDel数量 | InDels: {raw_stats.get('INDEL', 0)}\n")
            
            # 过滤后变异数 | Filtered variants
            f.write(f"\n硬过滤后变异 | Hard filtered variants:\n")
            if filtered_vcfs['snp_hard'].exists():
                snp_stats = self._count_variants(filtered_vcfs['snp_hard'])
                f.write(f"  SNP数量 | SNPs: {snp_stats['total']}\n")
            
            if filtered_vcfs['indel_hard'].exists():
                indel_stats = self._count_variants(filtered_vcfs['indel_hard'])
                f.write(f"  InDel数量 | InDels: {indel_stats['total']}\n")
            
            # 过滤统计 | Filtering statistics
            raw_snp_vcf = self.config.vcf_dir / "raw_snp.vcf.gz"
            raw_indel_vcf = self.config.vcf_dir / "raw_indel.vcf.gz"
            
            if raw_snp_vcf.exists() and filtered_vcfs['snp_hard'].exists():
                raw_snp_count = self._count_variants(raw_snp_vcf)['total']
                filtered_snp_count = self._count_variants(filtered_vcfs['snp_hard'])['total']
                snp_filtered = raw_snp_count - filtered_snp_count
                snp_filter_rate = (snp_filtered / raw_snp_count * 100) if raw_snp_count > 0 else 0
                
                f.write(f"\n过滤统计 | Filtering statistics:\n")
                f.write(f"  SNP过滤数量 | SNPs filtered: {snp_filtered}\n")
                f.write(f"  SNP过滤率 | SNP filter rate: {snp_filter_rate:.2f}%\n")
            
            if raw_indel_vcf.exists() and filtered_vcfs['indel_hard'].exists():
                raw_indel_count = self._count_variants(raw_indel_vcf)['total']
                filtered_indel_count = self._count_variants(filtered_vcfs['indel_hard'])['total']
                indel_filtered = raw_indel_count - filtered_indel_count
                indel_filter_rate = (indel_filtered / raw_indel_count * 100) if raw_indel_count > 0 else 0
                
                f.write(f"  InDel过滤数量 | InDels filtered: {indel_filtered}\n")
                f.write(f"  InDel过滤率 | InDel filter rate: {indel_filter_rate:.2f}%\n")
            
            f.write("\n")
            
            # 4. 输出文件 | Output files
            f.write("📁 输出文件 | Output Files\n")
            f.write("-" * 80 + "\n")
            f.write(f"\nBAM文件目录 | BAM directory: {self.config.bam_dir}\n")
            f.write(f"GVCF文件目录 | GVCF directory: {self.config.gvcf_dir}\n")
            f.write(f"VCF文件目录 | VCF directory: {self.config.vcf_dir}\n")
            
            f.write(f"\n最终VCF文件 | Final VCF files:\n")
            f.write(f"  硬过滤SNP | Hard filtered SNPs: {filtered_vcfs['snp_hard'].name}\n")
            f.write(f"  硬过滤InDel | Hard filtered InDels: {filtered_vcfs['indel_hard'].name}\n")
            f.write(f"  软过滤SNP | Soft filtered SNPs: {filtered_vcfs['snp_soft'].name}\n")
            f.write(f"  软过滤InDel | Soft filtered InDels: {filtered_vcfs['indel_soft'].name}\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        self.logger.info(f"✅ 统计报告已生成 | Statistics report generated: {stats_file}")
    
    def _get_alignment_stats(self, bam_file: Path) -> dict:
        """获取比对统计 | Get alignment statistics"""
        try:
            # samtools flagstat
            result = subprocess.run(
                [self.config.samtools_path, 'flagstat', str(bam_file)],
                capture_output=True, text=True, check=True
            )
            
            stats = {}
            lines = result.stdout.strip().split('\n')
            
            for line in lines:
                if 'in total' in line:
                    stats['total_reads'] = line.split()[0]
                elif 'mapped (' in line and 'primary' not in line:
                    parts = line.split()
                    stats['mapped_reads'] = parts[0]
                    stats['mapping_rate'] = parts[4].strip('()')
                elif 'properly paired' in line:
                    stats['properly_paired'] = line.split()[4].strip('()')
            
            # 获取重复率 | Get duplicate rate
            metrics_file = self.config.stats_dir / f"{bam_file.stem.replace('.dedup', '')}.dedup_metrics.txt"
            if metrics_file.exists():
                with open(metrics_file) as f:
                    for line in f:
                        if line.startswith('LIBRARY'):
                            next_line = next(f)
                            parts = next_line.strip().split('\t')
                            if len(parts) > 8:
                                stats['duplicate_rate'] = f"{float(parts[8])*100:.2f}%"
                            break
            
            return stats
        except Exception as e:
            self.logger.warning(f"⚠️  无法获取比对统计 | Cannot get alignment stats: {e}")
            return {}
    
    def _count_variants(self, vcf_file: Path) -> dict:
        """统计变异数量 | Count variants"""
        try:
            # 使用bcftools或grep统计 | Use bcftools or grep to count
            result = subprocess.run(
                f"zgrep -v '^#' {vcf_file} | wc -l",
                shell=True, capture_output=True, text=True, check=True
            )
            total = int(result.stdout.strip())
            
            return {'total': total}
        except Exception as e:
            self.logger.warning(f"⚠️  无法统计变异 | Cannot count variants: {e}")
            return {'total': 0}
