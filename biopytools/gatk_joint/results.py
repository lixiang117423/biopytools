"""
结果汇总模块 | Results Summary Module
"""

import subprocess
from pathlib import Path

class ResultsSummary:
    """结果汇总生成器 | Results Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self):
        """生成结果汇总 | Generate summary report"""
        self.logger.info("📊 生成结果汇总报告 | Generating summary report")
        
        report_file = self.config.output_path / "analysis_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("GATK Joint Genotyping 分析汇总报告\n")
            f.write("GATK Joint Genotyping Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 输入信息 | Input information
            f.write("📥 输入信息 | Input Information:\n")
            f.write(f"  - 输入目录 | Input directory: {self.config.input_dir}\n")
            f.write(f"  - 参考基因组 | Reference genome: {self.config.reference}\n")
            f.write(f"  - 文件类型 | File type: {self.config.file_type.upper()}\n")
            if self.config.intervals:
                f.write(f"  - 分析区间 | Analysis intervals: {self.config.intervals}\n")
            f.write(f"\n")
            
            # 计算参数 | Computing parameters
            f.write("⚙️ 计算参数 | Computing Parameters:\n")
            f.write(f"  - 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  - 内存设置 | Memory: {self.config.memory}\n")
            f.write(f"\n")
            
            # 过滤参数 | Filtering parameters
            f.write("🔬 过滤参数 | Filtering Parameters:\n")
            f.write(f"  SNP过滤 | SNP filters:\n")
            f.write(f"    QD < {self.config.snp_qd}, FS > {self.config.snp_fs}, MQ < {self.config.snp_mq}\n")
            f.write(f"    MQRankSum < {self.config.snp_mqrs}, ReadPosRankSum < {self.config.snp_rprs}\n")
            f.write(f"    SOR > {self.config.snp_sor}\n")
            f.write(f"  INDEL过滤 | INDEL filters:\n")
            f.write(f"    QD < {self.config.indel_qd}, FS > {self.config.indel_fs}\n")
            f.write(f"    ReadPosRankSum < {self.config.indel_rprs}, SOR > {self.config.indel_sor}\n")
            f.write(f"\n")
            
            # 输出文件 | Output files
            f.write("📤 输出文件 | Output Files:\n")
            f.write(f"  - 原始VCF | Raw VCF: {self.config.raw_vcf}\n")
            f.write(f"  - SNP过滤VCF | SNP filtered VCF: {self.config.snp_filtered_vcf}\n")
            f.write(f"  - INDEL过滤VCF | INDEL filtered VCF: {self.config.indel_filtered_vcf}\n")
            f.write(f"  - 合并过滤VCF | Merged filtered VCF: {self.config.merged_vcf}\n")
            f.write(f"\n")
            
            # 统计信息 | Statistics
            f.write("📈 变异统计 | Variant Statistics:\n")
            self._add_variant_stats(f)
            
            f.write(f"\n输出目录 | Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"✅ 汇总报告已生成 | Summary report generated: {report_file}")
    
    def _add_variant_stats(self, file_handle):
        """添加变异统计信息 | Add variant statistics"""
        try:
            # 统计原始VCF | Count raw VCF
            if hasattr(self.config, 'raw_vcf'):
                raw_count = self._count_variants(self.config.raw_vcf)
                file_handle.write(f"  - 原始变异数 | Raw variants: {raw_count}\n")
            
            # 统计SNP | Count SNPs
            if hasattr(self.config, 'snp_filtered_vcf'):
                snp_count = self._count_variants(self.config.snp_filtered_vcf, pass_only=True)
                file_handle.write(f"  - 过滤后SNP数 | Filtered SNPs (PASS): {snp_count}\n")
            
            # 统计INDEL | Count INDELs
            if hasattr(self.config, 'indel_filtered_vcf'):
                indel_count = self._count_variants(self.config.indel_filtered_vcf, pass_only=True)
                file_handle.write(f"  - 过滤后INDEL数 | Filtered INDELs (PASS): {indel_count}\n")
            
            # 统计合并结果 | Count merged results
            if hasattr(self.config, 'merged_vcf'):
                merged_count = self._count_variants(self.config.merged_vcf, pass_only=True)
                file_handle.write(f"  - 合并变异数 | Merged variants (PASS): {merged_count}\n")
                
        except Exception as e:
            file_handle.write(f"  统计信息获取失败 | Failed to get statistics: {e}\n")
    
    def _count_variants(self, vcf_file, pass_only=False):
        """统计VCF中的变异数 | Count variants in VCF"""
        try:
            if pass_only:
                cmd = f"{self.config.bcftools_path} view -f PASS {vcf_file} | grep -v '^#' | wc -l"
            else:
                cmd = f"zgrep -v '^#' {vcf_file} | wc -l"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return result.stdout.strip()
        except:
            return "N/A"
