"""
重复序列统计分析模块 |Duplicate Statistics Analysis Module
"""

from pathlib import Path
from typing import Dict, Any
from .utils import CommandRunner

class DuplicateStatsAnalyzer:
    """重复序列统计分析器 |Duplicate Statistics Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def get_duplicate_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取重复序列统计|Get duplicate statistics"""
        self.logger.info(f" 分析重复序列|Analyzing duplicates: {bam_file}")
        
        # 检查是否已标记重复|Check if duplicates are already marked
        marked_duplicates = self._check_marked_duplicates(bam_file)
        
        if marked_duplicates:
            return self._get_flagged_duplicate_stats(bam_file)
        else:
            return self._estimate_duplicate_rate(bam_file)
    
    def _check_marked_duplicates(self, bam_file: str) -> bool:
        """检查BAM文件是否已标记重复 |Check if BAM file has marked duplicates"""
        cmd = f"{self.config.samtools_path} view -f 1024 {bam_file}|head -1"
        success, output = self.cmd_runner.run(cmd, "Check marked duplicates", use_threads=True)
        
        return success and output.strip() != ""
    
    def _get_flagged_duplicate_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取已标记重复的统计 |Get statistics for flagged duplicates"""
        # 总读段数（支持多线程）| Total reads (with multi-threading)
        cmd_total = f"{self.config.samtools_path} view -c {bam_file}"
        success_total, total_output = self.cmd_runner.run(cmd_total, "Count total reads", use_threads=True)
        
        # 重复读段数（支持多线程）| Duplicate reads (with multi-threading)
        cmd_dup = f"{self.config.samtools_path} view -c -f 1024 {bam_file}"
        success_dup, dup_output = self.cmd_runner.run(cmd_dup, "Count duplicate reads", use_threads=True)
        
        if not (success_total and success_dup):
            return {}
        
        try:
            total_reads = int(total_output.strip())
            duplicate_reads = int(dup_output.strip())
            duplicate_rate = (duplicate_reads / total_reads * 100) if total_reads > 0 else 0
            
            return {
                'duplicates_marked': True,
                'total_reads': total_reads,
                'duplicate_reads': duplicate_reads,
                'duplicate_rate': duplicate_rate,
                'unique_reads': total_reads - duplicate_reads
            }
        except ValueError:
            return {}
    
    def _estimate_duplicate_rate(self, bam_file: str) -> Dict[str, Any]:
        """估算重复率|Estimate duplicate rate"""
        self.logger.warning("BAM文件未标记重复，进行简单估算|BAM file has no marked duplicates, performing simple estimation")
        
        # 通过坐标统计可能的重复|Count potential duplicates by coordinates
        cmd = f"{self.config.samtools_path} view -F 4 {bam_file}|head -100000|awk '{{print $3,$4,$6}}'|sort|uniq -c|awk '$1>1{{sum+=$1}} END{{print sum+0}}'"
        success, output = self.cmd_runner.run(cmd, "Estimate duplicates")
        
        cmd_total = f"{self.config.samtools_path} view -c -F 4 {bam_file}|head -100000"
        success_total, total_output = self.cmd_runner.run(cmd_total, "Count mapped reads for estimation")
        
        if not (success and success_total):
            return {'duplicates_marked': False, 'estimation_available': False}
        
        try:
            estimated_duplicates = int(output.strip()) if output.strip() else 0
            total_sampled = min(100000, int(total_output.strip()))
            estimated_rate = (estimated_duplicates / total_sampled * 100) if total_sampled > 0 else 0
            
            return {
                'duplicates_marked': False,
                'estimation_available': True,
                'estimated_duplicate_rate': estimated_rate,
                'sampled_reads': total_sampled,
                'note': 'This is an rough estimation based on coordinate analysis'
            }
        except ValueError:
            return {'duplicates_marked': False, 'estimation_available': False}
