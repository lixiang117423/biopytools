"""
SRA转换主程序模块|SRA Conversion Main Module
"""

import argparse
import sys
from datetime import datetime
from .config import ConvertConfig
from .utils import ConvertLogger, CommandRunner, check_dependencies
from .processor import SRAProcessor
from .report import ReportGenerator

class SRAConverter:
    """SRA转换主类|Main SRA Converter Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ConvertConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = ConvertLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化处理器|Initialize processors
        self.processor = SRAProcessor(self.config, self.logger, self.cmd_runner)
        self.report_generator = ReportGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_conversion(self):
        """运行完整的转换流程|Run complete conversion pipeline"""
        start_time = datetime.now()
        
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始SRA转FASTQ转换流程|Starting SRA to FASTQ Conversion Pipeline")
            self.logger.info("=" * 80)
            
            # 1. 检查依赖|Check dependencies
            self.check_dependencies()
            
            # 2. 显示配置信息|Display configuration
            self.logger.info(f"转换配置|Conversion Configuration:")
            self.logger.info(f"  输入路径|Input: {self.config.input_path}")
            self.logger.info(f"  输出目录|Output: {self.config.output_dir}")
            self.logger.info(f"  文件数量|Files: {len(self.config.input_files)}")
            
            tool_name = "parallel-fastq-dump (推荐)" if self.config.use_parallel else "fastq-dump (备选)"
            self.logger.info(f"  转换工具|Tool: {tool_name}")
            self.logger.info(f"  线程数|Threads: {self.config.threads}")
            self.logger.info(f"  压缩输出|Compress: {self.config.compress}")
            
            if self.config.tmpdir:
                self.logger.info(f"  临时目录|Temp Dir: {self.config.tmpdir}")
            
            # 3. 执行转换|Execute conversion
            self.logger.info(f"开始批量转换|Starting batch conversion...")
            results = self.processor.convert_all_files()
            
            end_time = datetime.now()
            
            # 4. 生成报告|Generate report
            self.logger.info(f"生成转换报告|Generating conversion report...")
            self.report_generator.generate_summary(results, start_time, end_time)
            
            # 5. 输出结果|Output results
            duration = (end_time - start_time).total_seconds()
            
            self.logger.info(f"=" * 80)
            self.logger.info(f"转换流程完成|Conversion Pipeline Completed")
            self.logger.info(f"=" * 80)
            self.logger.info(f"总耗时|Duration: {duration:.1f} 秒|seconds")
            self.logger.info(f"成功: {len(results['success'])}|Success: {len(results['success'])}")
            self.logger.info(f"失败: {len(results['failed'])}|Failed: {len(results['failed'])}")
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"转换流程在执行过程中意外终止|Conversion pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='SRA转FASTQ高速转换工具(parallel-fastq-dump版本)|SRA to FASTQ High-Speed Conversion Tool (parallel-fastq-dump version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i sra_dir/ -o fastq_output
        """
    )
    
    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True, 
                       help='输入SRA文件或文件夹路径|Input SRA file or folder path')
    
    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output', default='./fastq_output', 
                       help='输出目录|Output directory')
    
    # 转换参数|Conversion parameters
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    
    parser.add_argument('-c', '--compress', action='store_true', default=True,
                       help='压缩输出为.gz格式|Compress output to .gz format')
    parser.add_argument('--no-compress', dest='compress', action='store_false',
                       help='不压缩输出|Do not compress output')
    
    parser.add_argument('--split', action='store_true', default=True,
                       help='拆分双端测序文件|Split paired-end reads')
    parser.add_argument('--no-split', dest='split_files', action='store_false',
                       help='不拆分文件|Do not split files')
    
    # 高级参数|Advanced parameters
    parser.add_argument('--tmpdir', type=str, default=None,
                       help='临时目录 (可以设置为高速存储以加速)|Temporary directory (can be set to fast storage for acceleration)')
    
    # 过滤参数|Filtering parameters
    parser.add_argument('--min-len', type=int, default=0, 
                       help='最小读长过滤|Minimum read length filter')
    
    parser.add_argument('--clip', action='store_true',
                       help='剪切adapters|Clip adapters')
    
    args = parser.parse_args()
    
    # 创建转换器并运行|Create converter and run
    converter = SRAConverter(
        input_path=args.input,
        output_dir=args.output,
        compress=args.compress,
        threads=args.threads,
        split_files=args.split if hasattr(args, 'split') else args.split_files,
        skip_technical=True,
        clip=args.clip,
        min_read_len=args.min_len,
        tmpdir=args.tmpdir
    )
    
    converter.run_conversion()

if __name__ == "__main__":
    main()
