"""
转换结果报告模块 |Conversion Results Report Module
"""

from pathlib import Path
from datetime import datetime

class ReportGenerator:
    """报告生成器|Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, results: dict, start_time: datetime, end_time: datetime):
        """生成转换摘要报告|Generate conversion summary report"""
        report_file = self.config.output_path / "conversion_summary.txt"
        
        total = results['total']
        success_count = len(results['success'])
        failed_count = len(results['failed'])
        duration = (end_time - start_time).total_seconds()
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("SRA转FASTQ转换摘要报告|SRA to FASTQ Conversion Summary Report\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"开始时间|Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"结束时间|End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总耗时|Total Duration: {duration:.1f} 秒|seconds\n\n")
            
            f.write(f"输入路径|Input Path: {self.config.input_path}\n")
            f.write(f"输出目录|Output Directory: {self.config.output_dir}\n\n")

            f.write(f"转换配置|Conversion Settings:\n")
            tool_name = "parallel-fastq-dump " if self.config.use_parallel else "fastq-dump"
            f.write(f"  - 转换工具|Tool: {tool_name}\n")
            f.write(f"  - 线程数|Threads: {self.config.threads}\n")
            f.write(f"  - 压缩输出|Compress: {'是|Yes' if self.config.compress else '否|No'}\n")
            f.write(f"  - 拆分双端|Split Files: {'是|Yes' if self.config.split_files else '否|No'}\n")
            
            if self.config.tmpdir:
                f.write(f"  - 临时目录|Temp Dir: {self.config.tmpdir}\n")
            if self.config.min_read_len > 0:
                f.write(f"  - 最小读长|Min Read Length: {self.config.min_read_len}\n")
            
            f.write(f"转换统计|Conversion Statistics:\n")
            f.write(f"  - 总文件数|Total Files: {total}\n")
            f.write(f"  - 成功|Success: {success_count}\n")
            f.write(f"  - 失败|Failed: {failed_count}\n")
            f.write(f"  - 成功率|Success Rate: {success_count/total*100:.1f}%\n")
            
            if success_count > 0:
                avg_time = duration / success_count
                f.write(f"  - 平均速度|Avg Speed: {avg_time:.1f} 秒/文件|sec/file\n")
            
            f.write(f"\n")
            
            if results['success']:
                f.write(f"成功转换的文件 ({len(results['success'])})|Successfully Converted Files:\n")
                for sra in results['success']:
                    f.write(f"  - {Path(sra).name}\n")
                f.write("\n")

            if results['failed']:
                f.write(f"转换失败的文件 ({len(results['failed'])})|Failed Conversions:\n")
                for sra in results['failed']:
                    f.write(f"  - {Path(sra).name}\n")
                f.write("\n")

            f.write(f"详细日志|Detailed Log: {self.config.output_path / 'sra_conversion.log'}\n")
            f.write("=" * 80 + "\n")

        self.logger.info(f"转换摘要已生成|Summary report generated: {report_file}")
