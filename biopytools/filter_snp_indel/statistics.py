"""
VCF统计模块 | VCF Statistics Module
"""

import subprocess
from pathlib import Path

class VCFStatistics:
    """VCF文件统计器 | VCF File Statistics"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    # def count_variants(self, vcf_file: Path) -> int:
    #     """统计变异数量 | Count variants"""
    #     try:
    #         cmd = f"{self.config.bcftools_path} view -H {vcf_file} | wc -l"
    #         result = subprocess.run(
    #             cmd,
    #             shell=True,
    #             capture_output=True,
    #             text=True,
    #             check=True
    #         )
    #         count = int(result.stdout.strip())
    #         return count
    #     except Exception as e:
    #         self.logger.error(f"❌ 统计失败 | Count failed: {e}")
    #         return 0

    def count_variants(self, vcf_file_name: Path) -> int:
        """统计变异数量 | Count variants"""
        try:
            # 关键修正：构建文件的完整路径
            # The key fix: construct the full path to the file
            full_vcf_path = self.config.output_path / vcf_file_name

            # 在执行命令前，先检查文件是否存在，增加健壮性
            if not full_vcf_path.exists():
                self.logger.warning(f"⚠️ 统计时文件未找到 | File not found during statistics: {full_vcf_path}")
                return 0

            # 使用文件的完整路径执行命令
            cmd = f"{self.config.bcftools_path} view -H {full_vcf_path} | wc -l"
            
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            count = int(result.stdout.strip())
            return count
        except Exception as e:
            self.logger.error(f"❌ 统计失败 | Count failed for file {vcf_file_name}: {e}")
            return 0
    
    def generate_statistics_report(self, separator, filter_obj) -> bool:
        """生成统计报告 | Generate statistics report"""
        self.logger.info("=" * 60)
        self.logger.info("📊 步骤4: 生成统计报告 | Step 4: Generating statistics report")
        self.logger.info("=" * 60)
        
        # 统计各个文件的变异数量 | Count variants in each file
        stats = {
            "原始SNP | Raw SNPs": self.count_variants(separator.raw_snp_file),
            "过滤后SNP | Filtered SNPs": self.count_variants(filter_obj.filtered_snp_file),
            "原始INDEL | Raw INDELs": self.count_variants(separator.raw_indel_file),
            "过滤后INDEL | Filtered INDELs": self.count_variants(filter_obj.filtered_indel_file),
            "合并后变异 | Merged variants": self.count_variants(filter_obj.merged_file)
        }
        
        # 计算过滤率 | Calculate filtering rates
        snp_rate = 0.0
        indel_rate = 0.0
        
        if stats["原始SNP | Raw SNPs"] > 0:
            snp_rate = (stats["过滤后SNP | Filtered SNPs"] / stats["原始SNP | Raw SNPs"]) * 100
        
        if stats["原始INDEL | Raw INDELs"] > 0:
            indel_rate = (stats["过滤后INDEL | Filtered INDELs"] / stats["原始INDEL | Raw INDELs"]) * 100
        
        # 生成报告文件 | Generate report file
        report_file = self.config.output_path / "filtering_statistics.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("🧬 VCF SNP/INDEL 过滤统计报告 | VCF SNP/INDEL Filtering Statistics Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("📁 输入文件 | Input file:\n")
            f.write(f"   {self.config.vcf_file}\n\n")
            
            f.write("📂 输出目录 | Output directory:\n")
            f.write(f"   {self.config.output_dir}\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("📊 变异统计 | Variant Statistics\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"🔹 原始SNP数量 | Raw SNPs:              {stats['原始SNP | Raw SNPs']:>12,}\n")
            f.write(f"🔸 过滤后SNP数量 | Filtered SNPs:       {stats['过滤后SNP | Filtered SNPs']:>12,}\n")
            f.write(f"📈 SNP保留率 | SNP retention rate:    {snp_rate:>11.2f}%\n\n")
            
            f.write(f"🔹 原始INDEL数量 | Raw INDELs:          {stats['原始INDEL | Raw INDELs']:>12,}\n")
            f.write(f"🔸 过滤后INDEL数量 | Filtered INDELs:   {stats['过滤后INDEL | Filtered INDELs']:>12,}\n")
            f.write(f"📈 INDEL保留率 | INDEL retention rate: {indel_rate:>11.2f}%\n\n")
            
            f.write(f"🔗 合并后变异总数 | Total merged:       {stats['合并后变异 | Merged variants']:>12,}\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("⚙️ 过滤参数 | Filtering Parameters\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("SNP过滤参数 | SNP Filtering Parameters:\n")
            f.write(f"  • 最小质量值 | Min QUAL:              {self.config.snp_qual}\n")
            f.write(f"  • 最小深度 | Min DP:                  {self.config.snp_dp}\n")
            f.write(f"  • 最小比对质量 | Min MQ:              {self.config.snp_mq}\n")
            f.write(f"  • 质量/深度比 | Min QD:               {self.config.snp_qd}\n")
            f.write(f"  • 最大FisherStrand | Max FS:          {self.config.snp_fs}\n")
            f.write(f"  • 最大SOR | Max SOR:                  {self.config.snp_sor}\n")
            f.write(f"  • 最小MQRankSum:                      {self.config.snp_mqrs}\n")
            f.write(f"  • 最小ReadPosRankSum:                 {self.config.snp_rprs}\n\n")
            f.write(f"  • 最小MAF | Min MAF:                  {self.config.snp_maf}\n\n")  # [新增] 报告中显示MAF
            
            f.write("INDEL过滤参数 | INDEL Filtering Parameters:\n")
            f.write(f"  • 最小质量值 | Min QUAL:              {self.config.indel_qual}\n")
            f.write(f"  • 最小深度 | Min DP:                  {self.config.indel_dp}\n")
            f.write(f"  • 最小比对质量 | Min MQ:              {self.config.indel_mq}\n")
            f.write(f"  • 质量/深度比 | Min QD:               {self.config.indel_qd}\n")
            f.write(f"  • 最大FisherStrand | Max FS:          {self.config.indel_fs}\n")
            f.write(f"  • 最大SOR | Max SOR:                  {self.config.indel_sor}\n")
            f.write(f"  • 最小ReadPosRankSum:                 {self.config.indel_rprs}\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("📦 输出文件列表 | Output Files\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"1. 原始SNP文件 | Raw SNP file:\n")
            f.write(f"   {separator.raw_snp_file.name}\n\n")
            
            f.write(f"2. 过滤后SNP文件 | Filtered SNP file:\n")
            f.write(f"   {filter_obj.filtered_snp_file.name}\n\n")
            
            f.write(f"3. 原始INDEL文件 | Raw INDEL file:\n")
            f.write(f"   {separator.raw_indel_file.name}\n\n")
            
            f.write(f"4. 过滤后INDEL文件 | Filtered INDEL file:\n")
            f.write(f"   {filter_obj.filtered_indel_file.name}\n\n")
            
            f.write(f"5. 合并后文件 | Merged file:\n")
            f.write(f"   {filter_obj.merged_file.name}\n\n")
            
            f.write("=" * 80 + "\n")
        
        # 打印到日志 | Print to log
        self.logger.info("=" * 60)
        self.logger.info("📊 过滤统计 | Filtering Statistics")
        self.logger.info("=" * 60)
        
        for key, value in stats.items():
            self.logger.info(f"{key}: {value:,}")
        
        self.logger.info(f"📈 SNP保留率 | SNP retention rate: {snp_rate:.2f}%")
        self.logger.info(f"📈 INDEL保留率 | INDEL retention rate: {indel_rate:.2f}%")
        
        self.logger.info(f"📄 统计报告已保存 | Statistics report saved: {report_file}")
        
        return True
