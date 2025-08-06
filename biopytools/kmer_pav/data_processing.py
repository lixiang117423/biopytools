"""
🧮 K-mer数据处理模块 | K-mer Data Processing Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, SampleDiscovery

class GenomeKmerExtractor:
    """🧬 基因组K-mer提取器 | Genome K-mer Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def extract_genome_kmers(self):
        """🔍 从基因组提取k-mer | Extract k-mers from genome"""
        self.logger.info(f"🧬 从基因组提取{self.config.kmer_size}bp的k-mer | Extracting {self.config.kmer_size}bp k-mers from genome")
        
        # 构建命令 | Build command
        cmd_parts = [
            self.config.unikmer_path,
            "count",
            f"-k {self.config.kmer_size}",
            f"-j {self.config.threads}"
        ]
        
        if self.config.canonical:
            cmd_parts.append("-K --canonical")
        
        if self.config.sort_output:
            cmd_parts.append("-s --sort")
        
        output_prefix = self.config.output_path / f"genome_{self.config.kmer_size}mers"
        cmd_parts.extend([
            f"-o {output_prefix}",
            f'"{self.config.genome_file}"'
        ])
        
        cmd = " ".join(cmd_parts)
        log_file = self.config.logs_dir / "genome_kmer_extraction.log"
        
        success = self.cmd_runner.run(cmd, 
                                    f"🧬 基因组k-mer提取 | Genome k-mer extraction", 
                                    str(log_file))
        
        if success and self.config.genome_kmers_file.exists():
            self.logger.info("✅ 基因组k-mer提取完成 | Genome k-mer extraction completed")
            return True
        else:
            self.logger.error("❌ 基因组k-mer提取失败 | Genome k-mer extraction failed")
            return False
    
    def locate_kmers(self):
        """📍 获取k-mer在基因组中的位置信息 | Get k-mer positions in genome"""
        self.logger.info("📍 获取k-mer在基因组中的位置信息 | Getting k-mer positions in genome")
        
        cmd = (
            f"{self.config.unikmer_path} locate "
            f"-g \"{self.config.genome_file}\" "
            f"-o \"{self.config.genome_positions_file}\" "
            f"\"{self.config.genome_kmers_file}\""
        )
        
        log_file = self.config.logs_dir / "kmer_location.log"
        success = self.cmd_runner.run(cmd, 
                                    "📍 k-mer位置定位 | K-mer position location", 
                                    str(log_file))
        
        if success and self.config.genome_positions_file.exists():
            self.logger.info("✅ k-mer位置信息获取完成 | K-mer position information completed")
            return True
        else:
            self.logger.error("❌ k-mer位置定位失败 | K-mer position location failed")
            return False

class SampleKmerExtractor:
    """🧪 样本K-mer提取器 | Sample K-mer Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.sample_discovery = SampleDiscovery(config.fastq_dir, config.fastq_pattern, logger)
    
    def discover_and_process_samples(self):
        """🔍 发现并处理样本 | Discover and process samples"""
        # 发现样本 | Discover samples
        samples = self.sample_discovery.discover_samples()
        
        # 保存样本列表 | Save sample list
        with open(self.config.sample_list_file, 'w') as f:
            for sample in samples:
                f.write(f"{sample}\n")
        
        # 处理每个样本 | Process each sample
        failed_samples = []
        for i, sample in enumerate(samples, 1):
            self.logger.info(f"🧪 处理样本 {i}/{len(samples)}: {sample} | Processing sample {i}/{len(samples)}: {sample}")
            
            try:
                if self.extract_sample_kmers(sample):
                    self.logger.info(f"✅ 样本 {sample} k-mer提取完成 | Sample {sample} k-mer extraction completed")
                else:
                    failed_samples.append(sample)
                    self.logger.error(f"❌ 样本 {sample} k-mer提取失败 | Sample {sample} k-mer extraction failed")
            except Exception as e:
                failed_samples.append(sample)
                self.logger.error(f"❌ 样本 {sample} 处理出错: {e} | Sample {sample} processing error: {e}")
        
        if failed_samples:
            self.logger.warning(f"⚠️ 以下样本处理失败 | The following samples failed: {', '.join(failed_samples)}")
        
        return samples, failed_samples
    
    def extract_sample_kmers(self, sample_name: str):
        """🧮 从单个样本提取k-mer | Extract k-mers from single sample"""
        try:
            # 获取样本文件路径 | Get sample file paths
            r1_path, r2_path = self.sample_discovery.get_sample_files(sample_name)
            
            # 构建命令 | Build command
            cmd_parts = [
                self.config.unikmer_path,
                "count",
                f"-k {self.config.kmer_size}",
                f"-j {self.config.threads}"
            ]
            
            if self.config.canonical:
                cmd_parts.append("-K --canonical")
            
            if self.config.sort_output:
                cmd_parts.append("-s --sort")
            
            output_prefix = self.config.fastq_kmers_dir / f"{sample_name}_{self.config.kmer_size}mers"
            cmd_parts.extend([
                f"-o {output_prefix}",
                f'"{r1_path}" "{r2_path}"'
            ])
            
            cmd = " ".join(cmd_parts)
            log_file = self.config.logs_dir / f"{sample_name}_kmer_extraction.log"
            
            return self.cmd_runner.run(cmd, 
                                     f"🧪 样本{sample_name} k-mer提取 | Sample {sample_name} k-mer extraction", 
                                     str(log_file))
        
        except FileNotFoundError as e:
            self.logger.error(f"❌ 样本文件未找到: {e} | Sample files not found: {e}")
            return False

class PresenceAnalyzer:
    """🔍 存在性分析器 | Presence Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def analyze_presence(self, samples: list):
        """🎯 分析基因组k-mer在样本中的存在情况 | Analyze presence of genome k-mers in samples"""
        self.logger.info("🔍 检查基因组k-mer在样本中的存在情况 | Checking presence of genome k-mers in samples")
        
        failed_samples = []
        for i, sample in enumerate(samples, 1):
            self.logger.info(f"🧪 检查样本 {i}/{len(samples)}: {sample} | Checking sample {i}/{len(samples)}: {sample}")
            
            if self.check_sample_presence(sample):
                self.logger.info(f"✅ 样本 {sample} 存在性检查完成 | Sample {sample} presence check completed")
            else:
                failed_samples.append(sample)
                self.logger.error(f"❌ 样本 {sample} 存在性检查失败 | Sample {sample} presence check failed")
        
        if failed_samples:
            self.logger.warning(f"⚠️ 以下样本存在性检查失败 | The following samples failed presence check: {', '.join(failed_samples)}")
        
        return failed_samples
    
    def check_sample_presence(self, sample_name: str):
        """🔎 检查单个样本的存在性 | Check presence for single sample"""
        sample_kmer_file = self.config.fastq_kmers_dir / f"{sample_name}_{self.config.kmer_size}mers.unik"
        
        if not sample_kmer_file.exists():
            self.logger.error(f"❌ 样本k-mer文件不存在: {sample_kmer_file} | Sample k-mer file does not exist: {sample_kmer_file}")
            # 创建空的存在性文件 | Create empty presence file
            empty_file = self.config.presence_results_dir / f"{sample_name}_present.txt"
            empty_file.touch()
            return False
        
        # 计算交集 | Calculate intersection
        output_prefix = self.config.presence_results_dir / f"{sample_name}_present"
        cmd = (
            f"{self.config.unikmer_path} inter "
            f"\"{self.config.genome_kmers_file}\" "
            f"\"{sample_kmer_file}\" "
            f"-o \"{output_prefix}\""
        )
        
        log_file = self.config.logs_dir / f"{sample_name}_intersection.log"
        success = self.cmd_runner.run(cmd, 
                                    f"🔄 样本{sample_name}交集计算 | Sample {sample_name} intersection calculation", 
                                    str(log_file))
        
        # 转换为文本格式 | Convert to text format
        if success:
            unik_file = f"{output_prefix}.unik"
            txt_file = self.config.presence_results_dir / f"{sample_name}_present.txt"
            
            if os.path.exists(unik_file):
                view_cmd = f"{self.config.unikmer_path} view \"{unik_file}\" > \"{txt_file}\""
                return self.cmd_runner.run(view_cmd, f"🔄 转换样本{sample_name}结果为文本 | Convert sample {sample_name} result to text")
            else:
                # 创建空文件表示没有交集 | Create empty file indicating no intersection
                txt_file.touch()
                return True
        
        return False
