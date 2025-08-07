"""
PacBio HiFi结构变异报告生成模块 | PacBio HiFi Structural Variant Report Generation Module
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Dict
from .utils import SVCommandRunner, get_vcf_stats


class SVReporter:
    """结构变异报告生成器 | Structural Variant Reporter"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def generate_report(self) -> bool:
        """生成统计报告 | Generate statistical report"""
        self.runner.logger.info("步骤8: 生成统计报告 | Step 8: Generating statistical report...")
        
        report_file = self.config.results_dir / f"{self.config.sample_name}_SV_summary_report.txt"
        
        # 收集统计信息 | Collect statistics
        stats = self._collect_statistics()
        
        # 生成报告内容 | Generate report content
        report_content = self._generate_report_content(stats)
        
        # 写入报告文件 | Write report file
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.runner.logger.info("统计报告完成 | Statistical report completed.")
        return True
    
    def _collect_statistics(self) -> Dict:
        """收集统计信息 | Collect statistics"""
        stats = {}
        
        # 各工具检测结果 | Results from each tool
        callers = ["pbsv", "sniffles2", "cutesv"]
        for caller in callers:
            vcf_file = self.config.results_dir / f"{self.config.sample_name}_{caller}_filtered.vcf.gz"
            if vcf_file.exists():
                stats[f"{caller}_count"] = get_vcf_stats(vcf_file, self.runner.logger)['total_sv']
            else:
                stats[f"{caller}_count"] = 0
        
        # 合并后结果 | Merged results
        merged_vcf = self.config.results_dir / f"{self.config.sample_name}_merged_sv.vcf"
        if merged_vcf.exists():
            stats['merged'] = get_vcf_stats(merged_vcf, self.runner.logger)
        else:
            stats['merged'] = {key: 0 for key in ['total_sv', 'del', 'ins', 'dup', 'inv', 'bnd', 'tra']}
        
        # 大片段SV统计 | Large SV statistics
        sv_types = ["large_deletions", "large_insertions", "duplications", 
                   "inversions", "translocations", "very_large_sv"]
        for sv_type in sv_types:
            vcf_file = self.config.results_dir / f"{self.config.sample_name}_{sv_type}.vcf"
            if vcf_file.exists():
                stats[sv_type] = get_vcf_stats(vcf_file, self.runner.logger)['total_sv']
            else:
                stats[sv_type] = 0
        
        return stats
    
    def _generate_report_content(self, stats: Dict) -> str:
        """生成报告内容 | Generate report content"""
        return f"""========================================
PacBio HiFi 结构变异检测报告
========================================
样本名称: {self.config.sample_name}
分析日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
参考基因组: {self.config.ref_genome}

========================================
分析参数:
========================================
最小SV长度: {self.config.min_sv_length} bp
最小支持读数: {self.config.min_support}
质量阈值: {self.config.quality_threshold}
大片段SV阈值: {self.config.large_sv_threshold} bp
超大片段SV阈值: {self.config.very_large_sv_threshold} bp

========================================
各工具检测结果统计:
========================================

pbsv检测的SV数量: {stats.get('pbsv_count', 0)}
Sniffles2检测的SV数量: {stats.get('sniffles2_count', 0)}
cuteSV检测的SV数量: {stats.get('cutesv_count', 0)}

========================================
合并后SV类型统计:
========================================

总SV数量: {stats['merged']['total_sv']}
删除 (DEL): {stats['merged']['del']}
插入 (INS): {stats['merged']['ins']}
重复序列 (DUP): {stats['merged']['dup']}
倒位 (INV): {stats['merged']['inv']}
断点 (BND): {stats['merged']['bnd']}
易位 (TRA): {stats['merged']['tra']}

========================================
大片段SV统计:
========================================

大片段删除 (>{self.config.large_sv_threshold}bp): {stats.get('large_deletions', 0)}
大片段插入 (>{self.config.large_sv_threshold}bp): {stats.get('large_insertions', 0)}
重复序列: {stats.get('duplications', 0)}
倒位: {stats.get('inversions', 0)}
易位: {stats.get('translocations', 0)}
超大片段SV (>{self.config.very_large_sv_threshold}bp): {stats.get('very_large_sv', 0)}

========================================
文件输出位置:
========================================

主要结果文件:
- 合并SV结果: {self.config.results_dir}/{self.config.sample_name}_merged_sv.vcf
- 大片段删除: {self.config.results_dir}/{self.config.sample_name}_large_deletions.vcf
- 大片段插入: {self.config.results_dir}/{self.config.sample_name}_large_insertions.vcf
- 重复序列: {self.config.results_dir}/{self.config.sample_name}_duplications.vcf
- 倒位: {self.config.results_dir}/{self.config.sample_name}_inversions.vcf
- 易位: {self.config.results_dir}/{self.config.sample_name}_translocations.vcf
- 超大片段SV: {self.config.results_dir}/{self.config.sample_name}_very_large_sv.vcf

日志文件: {self.config.logs_dir}/
临时文件: {self.config.temp_dir}/
"""
    
    def generate_visualization_data(self) -> bool:
        """生成可视化数据 | Generate visualization data"""
        self.runner.logger.info("步骤9: 准备可视化数据 | Step 9: Preparing visualization data...")
        
        merged_vcf = self.config.results_dir / f"{self.config.sample_name}_merged_sv.vcf"
        if not merged_vcf.exists():
            self.runner.logger.warning("合并VCF文件不存在，跳过可视化数据生成 | Merged VCF not found, skipping visualization data generation")
            return True
        
        # SV长度分布 | SV length distribution
        lengths_file = self.config.results_dir / f"{self.config.sample_name}_sv_lengths.txt"
        self.runner.run(
            f"{self.config.bcftools_path} query -f '%SVLEN\\n' {merged_vcf} | "
            f"grep -v '\\.' | sed 's/-//' > {lengths_file}",
            "生成SV长度分布数据"
        )
        
        # SV类型分布 | SV type distribution
        types_file = self.config.results_dir / f"{self.config.sample_name}_sv_types.txt"
        self.runner.run(
            f"{self.config.bcftools_path} query -f '%SVTYPE\\n' {merged_vcf} | "
            f"sort | uniq -c > {types_file}",
            "生成SV类型分布数据"
        )
        
        # 染色体SV分布 | Chromosome SV distribution
        chrom_file = self.config.results_dir / f"{self.config.sample_name}_sv_by_chrom.txt"
        self.runner.run(
            f"{self.config.bcftools_path} query -f '%CHROM\\n' {merged_vcf} | "
            f"sort | uniq -c > {chrom_file}",
            "生成染色体SV分布数据"
        )
        
        self.runner.logger.info("可视化数据准备完成 | Visualization data preparation completed.")
        return True
