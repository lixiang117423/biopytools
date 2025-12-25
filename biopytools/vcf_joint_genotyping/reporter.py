"""
VCF联合分型报告生成模块 | VCF Joint Genotyping Report Generation Module
"""

import subprocess
from pathlib import Path
from typing import Dict, Any
from datetime import datetime

class VCFJointGenotypingReporter:
    """VCF联合分型报告生成器 | VCF Joint Genotyping Reporter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_paths = config.get_output_paths()
    
    def generate_report(self, stats: Dict[str, Any]) -> bool:
        """生成处理报告 | Generate processing report"""
        self.logger.info("生成处理报告 | Generating processing report")
        
        try:
            with open(self.output_paths['report_file'], 'w', encoding='utf-8') as f:
                self._write_header(f)
                self._write_basic_info(f)
                self._write_filtering_parameters(f)
                self._write_output_files(f)
                self._write_statistics(f, stats)
                self._write_chromosome_distribution(f)
            
            self.logger.info(f"处理报告已生成 | Processing report generated: {self.output_paths['report_file']}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成处理报告失败 | Failed to generate processing report: {e}")
            return False
    
    def _write_header(self, f):
        """写入报告头部 | Write report header"""
        f.write("VCF联合分型和过滤处理报告\n")
        f.write("VCF Joint Genotyping and Filtering Report\n")
        f.write("=" * 60 + "\n\n")
        
        # 添加处理时间 | Add processing time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"处理时间 | Processing time: {current_time}\n")
        f.write(f"报告生成时间 | Report generation time: {current_time}\n\n")
    
    def _write_basic_info(self, f):
        """写入基本信息 | Write basic information"""
        f.write("基本信息 | Basic Information\n")
        f.write("-" * 30 + "\n")
        f.write(f"输入目录 | Input directory: {self.config.vcf_input_dir}\n")
        f.write(f"输出目录 | Output directory: {self.config.output_dir}\n")
        f.write(f"参考基因组 | Reference genome: {self.config.reference_genome}\n")
        f.write(f"输入VCF文件数 | Input VCF files: {len(list(Path(self.config.vcf_input_dir).glob('*.vcf.gz')))}\n")
        f.write("\n")
    
    def _write_filtering_parameters(self, f):
        """写入过滤参数 | Write filtering parameters"""
        f.write("过滤参数 | Filtering Parameters\n")
        f.write("-" * 30 + "\n")
        
        f.write("SNP过滤参数 | SNP filtering parameters:\n")
        f.write(f"  - MAF >= {self.config.snp_maf}\n")
        f.write(f"  - 最大缺失率 | Max missing rate: {self.config.snp_max_missing}\n")
        f.write(f"  - Hardy-Weinberg平衡 | Hardy-Weinberg equilibrium: p >= {self.config.snp_hwe_pvalue}\n")
        f.write(f"  - 平均深度 | Mean depth: {self.config.snp_min_mean_dp}-{self.config.snp_max_mean_dp}\n")
        
        f.write("\nINDEL过滤参数 | INDEL filtering parameters:\n")
        f.write(f"  - MAF >= {self.config.indel_maf}\n")
        f.write(f"  - 最大缺失率 | Max missing rate: {self.config.indel_max_missing}\n")
        f.write(f"  - Hardy-Weinberg平衡 | Hardy-Weinberg equilibrium: p >= {self.config.indel_hwe_pvalue}\n")
        f.write(f"  - 平均深度 | Mean depth: {self.config.indel_min_mean_dp}-{self.config.indel_max_mean_dp}\n")
        f.write("\n")
    
    def _write_output_files(self, f):
        """写入输出文件信息 | Write output files information"""
        f.write("输出文件 | Output Files\n")
        f.write("-" * 30 + "\n")
        
        output_files = [
            ("样本映射文件 | Sample mapping file", self.output_paths['sample_map_file']),
            ("合并VCF | Merged VCF", self.output_paths['merged_vcf']),
            ("SNP VCF", self.output_paths['snp_vcf']),
            ("INDEL VCF", self.output_paths['indel_vcf']),
            ("最终过滤SNP VCF | Final filtered SNP VCF", self.output_paths['final_snp_vcf']),
            ("最终过滤INDEL VCF | Final filtered INDEL VCF", self.output_paths['final_indel_vcf']),
            ("SNP索引文件 | SNP index file", Path(f"{self.output_paths['final_snp_vcf']}.tbi")),
            ("INDEL索引文件 | INDEL index file", Path(f"{self.output_paths['final_indel_vcf']}.tbi"))
        ]
        
        for description, file_path in output_files:
            if file_path.exists():
                file_size = self._get_file_size(file_path)
                f.write(f"  - {description}: {file_path.name} ({file_size})\n")
            else:
                f.write(f"  - {description}: {file_path.name} (不存在 | Not exists)\n")
        
        f.write("\n")
    
    def _write_statistics(self, f, stats: Dict[str, Any]):
        """写入统计信息 | Write statistics"""
        f.write("最终统计 | Final Statistics\n")
        f.write("-" * 30 + "\n")
        
        if 'snp_count' in stats:
            f.write(f"过滤后SNP总数 | Filtered SNP count: {stats['snp_count']}\n")
            if 'snp_filter_rate' in stats:
                f.write(f"SNP过滤保留率 | SNP filter retention rate: {stats['snp_filter_rate']:.2f}%\n")
        
        if 'indel_count' in stats:
            f.write(f"过滤后INDEL总数 | Filtered INDEL count: {stats['indel_count']}\n")
            if 'indel_filter_rate' in stats:
                f.write(f"INDEL过滤保留率 | INDEL filter retention rate: {stats['indel_filter_rate']:.2f}%\n")
        
        if 'total_variants' in stats:
            f.write(f"过滤后变异总数 | Total filtered variants: {stats['total_variants']}\n")
        
        f.write("\n")
    
    def _write_chromosome_distribution(self, f):
        """写入染色体分布信息 | Write chromosome distribution"""
        f.write("染色体分布 | Chromosome Distribution\n")
        f.write("-" * 30 + "\n")
        
        # SNP染色体分布 | SNP chromosome distribution
        if self.output_paths['final_snp_vcf'].exists():
            f.write("SNP各染色体分布 | SNP chromosome distribution:\n")
            chr_dist = self._get_chromosome_distribution(self.output_paths['final_snp_vcf'])
            for chr_name, count in chr_dist[:10]:  # 只显示前10个
                f.write(f"  {chr_name}: {count} SNPs\n")
            if len(chr_dist) > 10:
                f.write(f"  ... (还有 | and {len(chr_dist) - 10} 个染色体 | more chromosomes)\n")
            f.write("\n")
        
        # INDEL染色体分布 | INDEL chromosome distribution
        if self.output_paths['final_indel_vcf'].exists():
            f.write("INDEL各染色体分布 | INDEL chromosome distribution:\n")
            chr_dist = self._get_chromosome_distribution(self.output_paths['final_indel_vcf'])
            for chr_name, count in chr_dist[:10]:  # 只显示前10个
                f.write(f"  {chr_name}: {count} INDELs\n")
            if len(chr_dist) > 10:
                f.write(f"  ... (还有 | and {len(chr_dist) - 10} 个染色体 | more chromosomes)\n")
    
    def _get_file_size(self, file_path: Path) -> str:
        """获取文件大小 | Get file size"""
        if file_path.exists():
            try:
                result = subprocess.run(f"du -h {file_path}", shell=True, capture_output=True, text=True)
                return result.stdout.split()[0]
            except:
                return "未知 | Unknown"
        return "不存在 | Not exists"
    
    def _get_chromosome_distribution(self, vcf_file: Path):
        """获取染色体分布 | Get chromosome distribution"""
        try:
            cmd = f"zcat {vcf_file} | grep -v '^#' | cut -f1 | sort | uniq -c | sort -nr"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            chr_dist = []
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        count = int(parts[0])
                        chr_name = parts[1]
                        chr_dist.append((chr_name, count))
            
            return chr_dist
        except:
            return []
