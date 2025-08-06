"""
📋 分析报告生成模块 | Analysis Report Generation Module
"""

import os
from datetime import datetime
from pathlib import Path

class AnalysisReporter:
    """📝 分析报告生成器 | Analysis Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_report(self):
        """📄 生成分析报告 | Generate analysis report"""
        self.logger.info("📝 生成分析报告 | Generating analysis report")
        
        # 读取样本数量 | Read sample count
        sample_count = 0
        if self.config.sample_list_file.exists():
            with open(self.config.sample_list_file, 'r') as f:
                sample_count = len([line for line in f if line.strip()])
        
        # 生成报告内容 | Generate report content
        report_content = self._create_report_content(sample_count)
        
        # 保存报告 | Save report
        report_file = self.config.output_path / "analysis_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        self.logger.info(f"✅ 分析报告已生成: {report_file} | Analysis report generated: {report_file}")
        
        # 显示主要结果文件 | Show main result files
        self._show_result_files()
    
    def _create_report_content(self, sample_count):
        """📄 创建报告内容 | Create report content"""
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        report = f"""🧬 K-mer PAV分析报告 | K-mer PAV Analysis Report
==========================================

📅 分析时间 | Analysis Time: {current_time}
🧬 基因组文件 | Genome File: {self.config.genome_file}
📂 FASTQ文件目录 | FASTQ Directory: {self.config.fastq_dir}
📁 输出目录 | Output Directory: {self.config.output_dir}

⚙️ 分析参数 | Analysis Parameters:
------------------------
🧮 K-mer长度 | K-mer Size: {self.config.kmer_size}bp
🔧 使用线程数 | Threads: {self.config.threads}
✅ Canonical K-mer: {self.config.canonical}
📊 排序输出 | Sort Output: {self.config.sort_output}
🎯 FASTQ文件模式 | FASTQ Pattern: {self.config.fastq_pattern}
📏 窗口大小 | Window Size: {self.config.window_size:,} bp
{'👣 步长 | Step Size: ' + str(self.config.step_size) + ' bp (重叠窗口 | Overlapping)' if self.config.overlapping else '📦 非重叠窗口 | Non-overlapping windows'}

🧪 样本信息 | Sample Information:
-------------------------
📊 总样本数 | Total Samples: {sample_count}

📁 输出文件 | Output Files:
--------------------
1. genome_{self.config.kmer_size}mers.unik - 🧬 基因组k-mer二进制文件 | Genome k-mer binary file
2. genome_{self.config.kmer_size}mers_positions.bed - 📍 k-mer位置信息 | K-mer position information
3. fastq_kmers/ - 📂 各样本k-mer文件目录 | Sample k-mer files directory
4. presence_results/ - 🎯 存在性检查结果目录 | Presence check results directory
5. kmer{self.config.kmer_size}_presence_matrix.csv - 📊 最终存在性矩阵 | Final presence matrix
6. kmer{self.config.kmer_size}_sample_stats.csv - 📈 样本统计信息 | Sample statistics
7. window_analysis/ - 🪟 窗口分析结果目录 | Window analysis results directory
   - kmer{self.config.kmer_size}_window_analysis_complete.csv - 📊 完整窗口分析结果 | Complete window analysis results
   - kmer{self.config.kmer_size}_window_presence_ratios.bed - 🛏️ BED格式窗口存在比例 | BED format window presence ratios
   - kmer{self.config.kmer_size}_window_sample_stats.csv - 📈 窗口样本统计 | Window sample statistics
8. logs/ - 📝 所有运行日志 | All running logs

✅ 分析状态 | Analysis Status: 完成 | Completed

📖 使用说明 | Usage Instructions:
------------------------
1. 📊 存在性矩阵文件可用于后续的群体遗传学分析
   The presence matrix file can be used for subsequent population genetics analysis
   
2. 📈 样本统计文件包含每个样本的k-mer检出率信息
   Sample statistics file contains k-mer detection rate information for each sample
   
3. 📍 位置信息文件提供每个k-mer在基因组中的坐标
   Position information file provides genomic coordinates for each k-mer

4. 🪟 窗口分析结果:
   Window analysis results:
   - 🛏️ BED格式文件包含每个窗口中每个样本的k-mer存在比例
     BED format file contains k-mer presence ratios for each sample in each window
   - 📊 可用于基因组区域水平的变异分析
     Can be used for genomic region-level variation analysis
   - 🔗 适合与其他基因组注释信息结合分析
     Suitable for combined analysis with other genomic annotation information

🆘 技术支持 | Technical Support:
---------------------------
如有问题，请检查logs目录中的详细日志文件
For issues, please check detailed log files in the logs directory
"""
        return report
    
    def _show_result_files(self):
        """📁 显示结果文件 | Show result files"""
        self.logger.info("📋 主要结果文件 | Main result files:")
        
        result_files = [
            f"kmer{self.config.kmer_size}_presence_matrix.csv",
            f"kmer{self.config.kmer_size}_sample_stats.csv",
            "analysis_report.txt"
        ]
        
        # 添加窗口分析结果文件 | Add window analysis result files
        window_files = [
            f"window_analysis/kmer{self.config.kmer_size}_window_analysis_complete.csv",
            f"window_analysis/kmer{self.config.kmer_size}_window_presence_ratios.bed",
            f"window_analysis/kmer{self.config.kmer_size}_window_sample_stats.csv"
        ]
        
        for filename in result_files:
            filepath = self.config.output_path / filename
            if filepath.exists():
                size = filepath.stat().st_size
                self.logger.info(f"  ✅ {filename} ({size:,} bytes)")
            else:
                self.logger.warning(f"  ❌ {filename} (文件不存在 | File not found)")
        
        # 检查窗口分析结果文件 | Check window analysis result files
        self.logger.info("🪟 窗口分析结果文件 | Window analysis result files:")
        for filename in window_files:
            filepath = self.config.output_path / filename
            if filepath.exists():
                size = filepath.stat().st_size
                self.logger.info(f"  ✅ {filename} ({size:,} bytes)")
            else:
                self.logger.warning(f"  ❌ {filename} (文件不存在 | File not found)")
