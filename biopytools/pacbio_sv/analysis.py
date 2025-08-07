"""
PacBio HiFi结构变异分析模块 | PacBio HiFi Structural Variant Analysis Module
"""

import os
from pathlib import Path
from typing import Dict, List
from .utils import SVCommandRunner, get_vcf_stats


class SVFilter:
    """结构变异过滤器 | Structural Variant Filter"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def filter_vcf_files(self) -> bool:
        """过滤VCF文件 | Filter VCF files"""
        self.runner.logger.info("步骤5: 过滤和质控结果 | Step 5: Filtering and QC results...")
        
        callers = ["pbsv", "sniffles2", "cutesv"]
        filter_cmd = f"QUAL>={self.config.quality_threshold} && abs(INFO/SVLEN)>={self.config.min_sv_length}"
        
        for caller in callers:
            input_vcf = self.config.results_dir / f"{self.config.sample_name}_{caller}.vcf"
            output_vcf = self.config.results_dir / f"{self.config.sample_name}_{caller}_filtered.vcf"
            
            if not input_vcf.exists():
                self.runner.logger.warning(f"{caller}结果文件不存在 | {caller} result file not found: {input_vcf}")
                continue
            
            self.runner.logger.info(f"过滤{caller}结果 | Filtering {caller} results...")
            if not self.runner.run(
                f"{self.config.bcftools_path} filter -i '{filter_cmd}' {input_vcf} > {output_vcf}",
                f"过滤{caller}结果"
            ):
                return False
            
            # 压缩和索引 | Compress and index
            if not self.runner.run(f"{self.config.bgzip_path} -f {output_vcf}", f"压缩{caller}结果"):
                return False
            if not self.runner.run(f"{self.config.tabix_path} -f -p vcf {output_vcf}.gz", f"索引{caller}结果"):
                return False
        
        self.runner.logger.info("过滤完成 | Filtering completed.")
        return True


class SVMerger:
    """结构变异合并器 | Structural Variant Merger"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def merge_results(self) -> bool:
        """合并多个工具的结果 | Merge results from multiple tools"""
        self.runner.logger.info("步骤6: 合并多个工具的结果 | Step 6: Merging results from multiple tools...")
        
        # 创建文件列表 | Create file list
        file_list = self.config.temp_dir / "sv_files.list"
        
        vcf_files = []
        callers = ["pbsv", "sniffles2", "cutesv"]
        for caller in callers:
            vcf_file = self.config.results_dir / f"{self.config.sample_name}_{caller}_filtered.vcf.gz"
            if vcf_file.exists():
                vcf_files.append(str(vcf_file))
        
        if len(vcf_files) < 2:
            self.runner.logger.warning("可用的VCF文件少于2个，跳过合并步骤 | Less than 2 VCF files available, skipping merge step")
            return True
        
        with open(file_list, 'w') as f:
            f.write('\n'.join(vcf_files))
        
        # 使用SURVIVOR合并 | Merge using SURVIVOR
        merged_vcf = self.config.results_dir / f"{self.config.sample_name}_merged_sv.vcf"
        if not self.runner.run(
            f"{self.config.survivor_path} merge {file_list} {self.config.survivor_distance} "
            f"{self.config.survivor_min_callers} 1 1 0 {self.config.min_sv_length} {merged_vcf}",
            "SURVIVOR合并"
        ):
            return False
        
        self.runner.logger.info("结果合并完成 | Results merging completed.")
        return True


class LargeSVAnalyzer:
    """大片段结构变异分析器 | Large Structural Variant Analyzer"""
    
    def __init__(self, config, runner: SVCommandRunner):
        self.config = config
        self.runner = runner
    
    def analyze_large_svs(self) -> bool:
        """分析大片段结构变异 | Analyze large structural variants"""
        self.runner.logger.info("步骤7: 分析大片段结构变异 | Step 7: Analyzing large structural variants...")
        
        merged_vcf = self.config.results_dir / f"{self.config.sample_name}_merged_sv.vcf"
        if not merged_vcf.exists():
            self.runner.logger.warning("合并VCF文件不存在，跳过大片段SV分析 | Merged VCF not found, skipping large SV analysis")
            return True
        
        # 定义SV类型和过滤条件 | Define SV types and filter conditions
        sv_analyses = {
            "large_deletions": f'INFO/SVTYPE=="DEL" && abs(INFO/SVLEN)>{self.config.large_sv_threshold}',
            "large_insertions": f'INFO/SVTYPE=="INS" && abs(INFO/SVLEN)>{self.config.large_sv_threshold}',
            "duplications": 'INFO/SVTYPE=="DUP"',
            "inversions": 'INFO/SVTYPE=="INV"',
            "translocations": 'INFO/SVTYPE=="BND" || INFO/SVTYPE=="TRA"',
            "very_large_sv": f'abs(INFO/SVLEN)>{self.config.very_large_sv_threshold}'
        }
        
        for sv_type, filter_condition in sv_analyses.items():
            output_vcf = self.config.results_dir / f"{self.config.sample_name}_{sv_type}.vcf"
            self.runner.logger.info(f"提取{sv_type} | Extracting {sv_type}...")
            
            if not self.runner.run(
                f"{self.config.bcftools_path} view -i '{filter_condition}' {merged_vcf} > {output_vcf}",
                f"提取{sv_type}"
            ):
                return False
        
        self.runner.logger.info("大片段SV分析完成 | Large SV analysis completed.")
        return True
