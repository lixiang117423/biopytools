"""
🧬 FASTP质控批处理命令 | FASTP Quality Control Batch Command
"""

import click
import sys
# 假设原始 main.py 文件在 ...fastp.main 路径下
# from ...fastp.main import main as fastp_main

# --- 为了让此代码块能独立运行，我们在此处包含原始 main 函数 ---
# --- 在您的项目中，您应该使用上面的 import 语句，并删除这部分 ---
# START: In-place original main function for demonstration
def get_original_main_for_demo():
    """
    这是一个辅助函数，用于将原始代码包含进来，以便此示例可以独立运行。
    在实际项目中，您应该删除这个函数，并使用 `from ...fastp.main import main as fastp_main`。
    """
    import argparse
    # from .config import FastpConfig
    # from .utils import FastpLogger, CommandRunner
    # from .data_processing import SampleFinder
    # from .processing import FastpCore
    # from .results import SummaryGenerator
    # NOTE: 下面的类是原始代码的简化模拟，以使 main 函数能够执行
    class MockConfig:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
            self.report_path = f"{kwargs['output_dir']}/reports"
        def validate(self): pass
    class MockLogger:
        def get_logger(self): return self
        def info(self, msg): print(f"INFO: {msg}")
        def error(self, msg): print(f"ERROR: {msg}", file=sys.stderr)
    class MockRunner:
        def __init__(self, logger): pass
    class MockFinder:
        def __init__(self, config, logger): pass
        def find_sample_pairs(self): print("INFO: 模拟查找样本..."); return [('sample1', 's1_1.fq.gz', 's1_2.fq.gz')]
        def validate_sample_pairs(self, pairs): return True
    class MockCore:
        def __init__(self, config, logger, runner): pass
        def validate_fastp(self): print("INFO: 模拟检查 fastp..."); return True
        def create_output_directories(self): print("INFO: 模拟创建输出目录...")
        def process_sample(self, *args): print(f"INFO: 模拟处理样本: {args[0]}..."); return True
    class MockSummary:
        def __init__(self, config, logger): pass
        def generate_summary_report(self, *args): print("INFO: 模拟生成总结报告...")

    class FastpProcessor:
        def __init__(self, **kwargs):
            self.config = MockConfig(**kwargs)
            self.config.validate()
            self.logger_manager = MockLogger()
            self.logger = self.logger_manager.get_logger()
            self.cmd_runner = MockRunner(self.logger)
            self.sample_finder = MockFinder(self.config, self.logger)
            self.fastp_core = MockCore(self.config, self.logger, self.cmd_runner)
            self.summary_generator = MockSummary(self.config, self.logger)
        def run_batch_processing(self):
            self.logger.info("="*60)
            self.logger.info("开始FASTQ数据质控批处理 | Starting FASTQ data quality control batch processing")
            self.logger.info("="*60)
            if not self.fastp_core.validate_fastp(): sys.exit(1)
            self.fastp_core.create_output_directories()
            sample_pairs = self.sample_finder.find_sample_pairs()
            if not self.sample_finder.validate_sample_pairs(sample_pairs): sys.exit(1)
            self.logger.info(f"找到 {len(sample_pairs)} 个有效样本配对...")
            successful_count = sum(1 for sp in sample_pairs if self.fastp_core.process_sample(*sp))
            failed_count = len(sample_pairs) - successful_count
            self.summary_generator.generate_summary_report(successful_count, failed_count, len(sample_pairs))
            self.logger.info("="*60)
            self.logger.info("FASTQ质控批处理完成！")
            self.logger.info(f"成功处理: {successful_count}, 失败: {failed_count}")
            self.logger.info("="*60)

    def main():
        parser = argparse.ArgumentParser(description="FASTQ Data Quality Control Batch Processing Script")
        parser.add_argument("-i", "--input-dir", required=True)
        parser.add_argument("-o", "--output-dir", required=True)
        parser.add_argument("--fastp-path", default="fastp")
        parser.add_argument("-t", "--threads", type=int, default=12)
        parser.add_argument("-q", "--quality-threshold", type=int, default=30)
        parser.add_argument("-l", "--min-length", type=int, default=50)
        parser.add_argument("-u", "--unqualified-percent", type=int, default=40)
        parser.add_argument("-n", "--n-base-limit", type=int, default=10)
        parser.add_argument("--read1-suffix", default="_1.fq.gz")
        parser.add_argument("--read2-suffix", default="_2.fq.gz")
        args = parser.parse_args()
        processor = FastpProcessor(
            input_dir=args.input_dir, output_dir=args.output_dir, fastp_path=args.fastp_path,
            threads=args.threads, quality_threshold=args.quality_threshold, min_length=args.min_length,
            unqualified_percent=args.unqualified_percent, n_base_limit=args.n_base_limit,
            read1_suffix=args.read1_suffix, read2_suffix=args.read2_suffix
        )
        processor.run_batch_processing()
    return main
fastp_main = get_original_main_for_demo()
# END: In-place original main function

@click.command(short_help = "使用FASTP进行FASTQ数据质控的批处理工具")
@click.option('--input-dir', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='📥 输入原始FASTQ数据目录 | Input raw FASTQ data directory')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='📤 输出清洁FASTQ数据目录 | Output clean FASTQ data directory')
@click.option('--fastp-path',
              default='fastp',
              show_default=True,
              help='🔧 fastp可执行文件路径 | fastp executable path')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='⚡ 线程数 | Number of threads')
@click.option('--quality-threshold', '-q',
              type=int,
              default=30,
              show_default=True,
              help='✅ 质量阈值 | Quality threshold')
@click.option('--min-length', '-l',
              type=int,
              default=50,
              show_default=True,
              help='📏 最小长度 | Minimum length')
@click.option('--unqualified-percent', '-u',
              type=int,
              default=40,
              show_default=True,
              help='🚮 不合格碱基百分比阈值 | Unqualified base percentage threshold')
@click.option('--n-base-limit', '-n',
              type=int,
              default=10,
              show_default=True,
              help='🚫 N碱基数量限制 | N base count limit')
@click.option('--read1-suffix',
              default='_1.fq.gz',
              show_default=True,
              help='🧬 Read1文件后缀 | Read1 file suffix')
@click.option('--read2-suffix',
              default='_2.fq.gz',
              show_default=True,
              help='🧬 Read2文件后缀 | Read2 file suffix')
def fastp(input_dir, output_dir, fastp_path, threads, quality_threshold,
          min_length, unqualified_percent, n_base_limit, read1_suffix,
          read2_suffix):
    """
    使用FASTP进行FASTQ数据质控的批处理工具.

    自动查找输入目录中的配对FASTQ文件，并对每个样本执行FASTP质控。
    
    示例 | Examples:
    
    \b
    # 🎯 基本用法 (指定输入和输出目录)
    biopytools fastp -i ./raw_data -o ./clean_data
    
    \b
    # 🔧 自定义线程数和质控标准
    biopytools fastp -i ./raw_data -o ./clean_data -t 24 -q 25 -l 75
    
    \b
    # 📄 指定不同的文件后缀
    biopytools fastp -i ./raw_data -o ./clean_data \\
        --read1-suffix _R1.fastq.gz \\
        --read2-suffix _R2.fastq.gz
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'fastp']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    
    # 可选参数 (仅在值不为默认值时添加，保持命令行简洁) ⚙️ | Optional parameters (add only when non-default)
    if fastp_path != 'fastp':
        args.extend(['--fastp-path', fastp_path])
    
    if threads != 12:
        args.extend(['-t', str(threads)])
        
    if quality_threshold != 30:
        args.extend(['-q', str(quality_threshold)])
        
    if min_length != 50:
        args.extend(['-l', str(min_length)])
        
    if unqualified_percent != 40:
        args.extend(['-u', str(unqualified_percent)])
        
    if n_base_limit != 10:
        args.extend(['-n', str(n_base_limit)])
        
    if read1_suffix != '_1.fq.gz':
        args.extend(['--read1-suffix', read1_suffix])
        
    if read2_suffix != '_2.fq.gz':
        args.extend(['--read2-suffix', read2_suffix])
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        fastp_main()
    except SystemExit as e:
        # 处理程序正常退出 (如 argparse 的 --help) ✅ | Handle normal program exit
        if e.code != 0:
            click.secho(f"❌ 脚本执行被终止，退出码: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        click.secho(f"❌ 发生未知错误 | An unexpected error occurred: {e}", fg='red', err=True)
        sys.exit(1)
    finally:
        # 无论如何都要恢复原始的 sys.argv | Restore original sys.argv regardless of outcome
        sys.argv = original_argv

# 如果直接运行此文件用于测试 | If running this file directly for testing
if __name__ == '__main__':
    # 模拟命令行调用，例如: python your_script_name.py -i . -o ./out
    fastp()