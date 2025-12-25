"""
🚀 K-mer提取主程序模块 | K-mer Extraction Main Module
"""

import argparse
import sys
import os
import glob
from .config import KmerConfig
from .utils import KmerLogger, CommandRunner, check_dependencies
from .file_handler import FileTypeDetector, FastqFileMatcher, FastqFileMerger
from .extractor_core import UnikmrExtractor, JellyfishExtractor
from .output_formatter import KmerOutputFormatter

class KmerExtractor:
    """🧬 K-mer提取主类 | Main K-mer Extractor Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = KmerConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = KmerLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.file_detector = FileTypeDetector(self.logger)
        self.fastq_matcher = FastqFileMatcher(self.logger)
        self.fastq_merger = FastqFileMerger(self.logger, self.config.output_path)
        self.output_formatter = KmerOutputFormatter(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """🔍 检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def detect_and_process_files(self):
        """🔎 检测和处理输入文件 | Detect and process input files"""
        self.logger.info(f"🔎 检测输入文件类型 | Detecting input file types")
        
        # 处理输入路径（可能是文件或目录）| Process input paths
        all_input_files = []
        for input_path in self.config.input_files:
            if os.path.isdir(input_path):
                self.logger.info(f"📂 检测到目录输入 | Detected directory input: {input_path}")
                
                # 根据文件类型搜索相应的文件
                if self.config.file_type == 'fastq' or not self.config.file_type:
                    patterns = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
                elif self.config.file_type == 'fasta':
                    patterns = ['*.fasta', '*.fa', '*.fas', '*.fasta.gz', '*.fa.gz', '*.fas.gz']
                else:
                    patterns = ['*']
                
                dir_files = []
                for pattern in patterns:
                    dir_files.extend(glob.glob(os.path.join(input_path, pattern)))
                
                self.logger.info(f"📄 在目录中找到 {len(dir_files)} 个文件 | Found {len(dir_files)} files in directory")
                all_input_files.extend(dir_files)
            else:
                all_input_files.append(input_path)
        
        # 更新配置中的文件列表
        self.config.input_files = all_input_files
        self.logger.info(f"📋 总共处理 {len(all_input_files)} 个输入文件 | Processing {len(all_input_files)} input files total")
        
        # 自动检测文件类型（如果未指定）| Auto-detect file type
        if not self.config.file_type and all_input_files:
            first_file = all_input_files[0]
            detected_type = self.file_detector.detect_file_type(first_file)
            self.config.file_type = detected_type
            self.logger.info(f"🎯 自动检测文件类型 | Auto-detected file type: {detected_type}")
        
        if self.config.file_type == 'fastq':
            return self._process_fastq_files()
        elif self.config.file_type == 'fasta':
            return self._process_fasta_files()
        else:
            raise ValueError(f"❌ 不支持的文件类型 | Unsupported file type: {self.config.file_type}")
    
    def _process_fastq_files(self):
        """📄 处理FASTQ文件 | Process FASTQ files"""
        self.logger.info(f"📄 处理FASTQ文件 | Processing FASTQ files")
        
        if self.config.fastq_pattern:
            # 使用模式匹配 | Use pattern matching
            self.logger.info(f"🔍 使用模式匹配 | Using pattern matching: {self.config.fastq_pattern}")
            fastq_samples = self.fastq_matcher.match_paired_files(self.config.input_files, self.config.fastq_pattern)
        else:
            # 单文件处理 | Single file processing
            self.logger.info(f"📄 单文件模式处理 | Single file mode processing")
            fastq_samples = {}
            for i, file_path in enumerate(self.config.input_files):
                sample_name = f"sample_{i+1}"
                fastq_samples[sample_name] = (file_path, None)
                self.logger.info(f"📄 单文件样品 | Single file sample: {sample_name} -> {os.path.basename(file_path)}")
        
        if not fastq_samples:
            raise ValueError(f"❌ 未找到任何FASTQ样品文件 | No FASTQ sample files found. 请检查文件路径和模式 | Please check file paths and pattern.")
        
        return 'fastq', fastq_samples, None
    
    def _process_fasta_files(self):
        """🧬 处理FASTA文件 | Process FASTA files"""
        self.logger.info(f"🧬 处理FASTA文件 | Processing FASTA files")
        return 'fasta', None, self.config.input_files
    
    def run_extraction(self):
        """🚀 运行完整的K-mer提取流程 | Run complete K-mer extraction pipeline"""
        try:
            self.logger.info(f"🚀 开始K-mer提取流程 | Starting K-mer extraction pipeline")
            self.logger.info(f"{'=' * 60}")
            
            # 1. 检查依赖 | Check dependencies
            self.logger.info(f"🔍 第1步：检查依赖软件 | Step 1: Checking dependencies")
            self.check_dependencies()
            
            # 2. 检测文件类型 | Detect file types
            self.logger.info(f"🔎 第2步：检测文件类型 | Step 2: Detecting file types")
            file_type, fastq_samples, fasta_files = self.detect_and_process_files()
            
            # 3. 提取k-mer | Extract k-mers
            self.logger.info(f"🧬 第3步：提取k-mer | Step 3: Extracting k-mers")
            
            if file_type == 'fastq':
                # 使用Jellyfish处理FASTQ文件 | Use Jellyfish for FASTQ files
                self.logger.info(f"🐟 使用Jellyfish处理FASTQ文件 | Using Jellyfish for FASTQ files")
                
                # 3.1 合并FASTQ文件 | Merge FASTQ files
                merged_r1, merged_r2 = self.fastq_merger.merge_fastq_files(fastq_samples)
                
                # 3.2 使用Jellyfish提取k-mer | Extract k-mers using Jellyfish
                jellyfish_extractor = JellyfishExtractor(self.config, self.logger, self.cmd_runner)
                binary_file = jellyfish_extractor.extract_kmers_from_fastq(merged_r1, merged_r2)
                
                # 3.3 获取统计信息 | Get statistics
                self.logger.info(f"📊 第4步：获取统计信息 | Step 4: Getting statistics")
                stats = jellyfish_extractor.get_kmer_statistics(binary_file)
                
                # 3.4 转换为FASTA格式 | Convert to FASTA format
                self.logger.info(f"📤 第5步：格式化输出 | Step 5: Formatting output")
                output_fasta = self.config.output_path / f"{self.config.base_name}.fasta"
                
                success = self.output_formatter.convert_jellyfish_to_fasta(binary_file, str(output_fasta))
                
                if not success:
                    raise RuntimeError("❌ FASTA文件生成失败 | FASTA file generation failed")
                
                # 清理临时合并文件 | Clean up temporary merged files
                if os.path.exists(merged_r1):
                    os.remove(merged_r1)
                    self.logger.info(f"🗑️ 清理临时文件 | Cleaned up temporary file: {os.path.basename(merged_r1)}")
                if merged_r2 and os.path.exists(merged_r2):
                    os.remove(merged_r2)
                    self.logger.info(f"🗑️ 清理临时文件 | Cleaned up temporary file: {os.path.basename(merged_r2)}")
                
                # 处理二进制文件 | Handle binary file
                if self.config.keep_binary:
                    self.logger.info(f"🗃️ 保留二进制文件 | Keeping binary file: {os.path.basename(binary_file)}")
                else:
                    if os.path.exists(binary_file):
                        os.remove(binary_file)
                        self.logger.info(f"🗑️ 删除二进制文件 | Removed binary file: {os.path.basename(binary_file)}")
                
            else:
                # 使用Unikmer处理FASTA文件 | Use Unikmer for FASTA files
                self.logger.info(f"🧬 使用Unikmer处理FASTA文件 | Using Unikmer for FASTA files")
                
                unikmer_extractor = UnikmrExtractor(self.config, self.logger, self.cmd_runner)
                unik_file = unikmer_extractor.extract_kmers_from_fasta(fasta_files)
                
                # 4. 获取统计信息 | Get statistics
                self.logger.info(f"📊 第4步：获取统计信息 | Step 4: Getting statistics")
                stats = unikmer_extractor.get_kmer_statistics(unik_file)
                
                # 5. 格式化输出 | Format output
                self.logger.info(f"📤 第5步：格式化输出 | Step 5: Formatting output")
                output_fasta = self.config.output_path / f"{self.config.base_name}.fasta"
                
                # 生成FASTA文件 | Generate FASTA file
                success = self.output_formatter.convert_unik_to_fasta(
                    unik_file, str(output_fasta), 
                    is_fastq_input=False,
                    fastq_samples=None,
                    fasta_files=fasta_files
                )
                
                if not success:
                    raise RuntimeError("❌ FASTA文件生成失败 | FASTA file generation failed")
                
                # 生成BED文件（仅用于FASTA输入且用户要求时）| Generate BED file
                if self.config.output_bed:
                    self.logger.info(f"🎯 用户要求生成BED文件 | User requested BED file generation")
                    output_bed = self.config.output_path / f"{self.config.base_name}.bed"
                    bed_success = self.output_formatter.generate_bed_file(unik_file, str(output_bed), fasta_files)
                    if not bed_success:
                        self.logger.error(f"❌ BED文件生成失败 | BED file generation failed")
                else:
                    self.logger.info(f"ℹ️ 跳过BED文件生成 (使用--output-bed启用) | Skipping BED file generation (use --output-bed to enable)")
            
            # 6. 完成 | Complete
            self.logger.info(f"✅ K-mer提取完成 | K-mer extraction completed")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"📊 提取统计 | Extraction Statistics:")
            self.logger.info(f"   🧬 K-mer长度 | K-mer length: {self.config.kmer_length}")
            self.logger.info(f"   🔧 使用工具 | Tool used: {'Jellyfish' if file_type == 'fastq' else 'Unikmer'}")
            
            if file_type == 'fastq':
                self.logger.info(f"   📈 不同K-mer数量 | Distinct k-mers: {stats.get('total_kmers', 'Unknown')}")
                self.logger.info(f"   🎯 唯一K-mer数量 | Unique k-mers: {stats.get('unique_kmers', 'Unknown')}")
            else:
                self.logger.info(f"   📈 K-mer数量 | Number of k-mers: {stats.get('total_kmers', 'Unknown')}")
            
            self.logger.info(f"   📁 输出文件 | Output files:")
            self.logger.info(f"     📄 FASTA: {output_fasta}")
            
            if file_type == 'fasta' and self.config.output_bed:
                self.logger.info(f"     📋 BED: {self.config.output_path / f'{self.config.base_name}.bed'}")
            elif file_type == 'fasta':
                self.logger.info(f"     ℹ️ BED文件未生成 (使用--output-bed启用) | BED file not generated (use --output-bed to enable)")
            
            if file_type == 'fastq' and self.config.keep_binary:
                self.logger.info(f"     🗃️ 二进制文件: {self.config.output_path / f'{self.config.base_name}.jf'}")
            
            self.logger.info(f"📂 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ 提取流程在执行过程中意外终止 | Extraction pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """🎯 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 K-mer提取工具 (FASTA使用unikmer，FASTQ使用jellyfish) | K-mer Extraction Tool (unikmer for FASTA, jellyfish for FASTQ)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🌟 使用示例 | Examples:
  %(prog)s -i data.fasta -o results
  %(prog)s -i sample1.fastq sample2.fastq -o results -k 31 -t 16
  %(prog)s -i *.fastq.gz -o results --fastq-pattern "*_1.clean.fq.gz" -k 21
  %(prog)s -i data.fasta -o results --output-bed -k 25
  %(prog)s --input-files file1.fa file2.fa -o results --threads 8 --memory 100 --output-bed
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input-files', nargs='+', required=True,
                       help='🎯 输入文件路径 (FASTA/FASTQ，支持压缩格式) | Input file paths (FASTA/FASTQ, supports compressed formats)')
    
    # 输出参数 | Output arguments
    parser.add_argument('-o', '--output-dir', default='./kmer_output',
                       help='📂 输出目录 | Output directory')
    
    # K-mer参数 | K-mer parameters
    parser.add_argument('-k', '--kmer-length', type=int, default=51,
                       help='🧬 K-mer长度 (1-64) | K-mer length (1-64)')
    
    # 性能参数 | Performance parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🚀 线程数 | Number of threads')
    parser.add_argument('-m', '--memory', type=int, default=880,
                       help='💾 内存限制(GB) | Memory limit (GB)')
    
    # 文件类型和模式 | File type and pattern
    parser.add_argument('--file-type', choices=['fasta', 'fastq'],
                       help='📁 文件类型 (自动检测如果未指定) | File type (auto-detect if not specified)')
    parser.add_argument('--fastq-pattern',
                       help='🔗 FASTQ文件匹配模式 (例如: "*_1.clean.fq.gz") | FASTQ file matching pattern (e.g., "*_1.clean.fq.gz")')
    
    # 处理选项 | Processing options
    parser.add_argument('--no-canonical', action='store_true',
                       help='🔄 不使用canonical k-mer | Do not use canonical k-mers')
    parser.add_argument('--no-compress', action='store_true',
                       help='🗜️ 不压缩输出文件 | Do not compress output files')
    parser.add_argument('--output-bed', action='store_true',
                       help='📋 输出BED格式文件 (仅适用于FASTA输入) | Output BED format file (only for FASTA input)')
    parser.add_argument('--no-keep-binary', action='store_true',
                       help='🗑️ 不保留二进制文件 (默认保留) | Do not keep binary files (keep by default)')
    
    # 工具路径 | Tool paths
    parser.add_argument('--unikmer-path', default='unikmer',
                       help='⚙️ Unikmer软件路径 | Unikmer software path')
    parser.add_argument('--jellyfish-path', default='jellyfish',
                       help='🐟 Jellyfish软件路径 | Jellyfish software path')
    
    # Jellyfish特有参数 | Jellyfish-specific parameters
    parser.add_argument('--jellyfish-hash-size', default='10000M',
                       help='🗂️ Jellyfish哈希表大小 | Jellyfish hash table size')
    
    args = parser.parse_args()
    
    # 创建提取器并运行 | Create extractor and run
    extractor = KmerExtractor(
        input_files=args.input_files,
        output_dir=args.output_dir,
        kmer_length=args.kmer_length,
        threads=args.threads,
        memory_gb=args.memory,
        file_type=args.file_type,
        fastq_pattern=args.fastq_pattern,
        canonical=not args.no_canonical,
        compress_output=not args.no_compress,
        output_bed=args.output_bed,
        keep_binary=not args.no_keep_binary,
        unikmer_path=args.unikmer_path,
        jellyfish_path=args.jellyfish_path,
        jellyfish_hash_size=args.jellyfish_hash_size
    )
    
    extractor.run_extraction()

if __name__ == "__main__":
    main()
