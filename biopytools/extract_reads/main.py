"""
基于contig-reads对应关系提取fastq reads主程序|Extract fastq reads by contig-reads mapping Main Module
"""

import argparse
import sys
from .config import ExtractReadsConfig
from .utils import ExtractReadsLogger
from .extractor import FastqReadsExtractor


class ReadsExtractor:
    """Reads提取器主类|Reads Extractor Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ExtractReadsConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = None
        self.logger_manager = ExtractReadsLogger(log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化提取器|Initialize extractor
        self.extractor = FastqReadsExtractor(self.config, self.logger)

    def run(self):
        """运行提取|Run extraction"""
        self.logger.info("开始Reads提取流程|Starting reads extraction pipeline")

        try:
            success = self.extractor.extract()

            if success:
                self.logger.info("Reads提取成功|Reads extraction completed successfully")
                return 0
            else:
                self.logger.error("Reads提取失败|Reads extraction failed")
                return 1

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return 1


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='基于contig-reads对应关系从fastq提取指定reads|Extract specified reads from fastq by contig-reads mapping',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-m', '--mapping',
                       required=True,
                       help='contig-reads对应关系文件(TSV格式)|contig-reads mapping file (TSV format)')
    parser.add_argument('-i', '--input',
                       required=True,
                       help='输入FASTQ文件(支持gzip压缩)|Input FASTQ file (gzip supported)')
    parser.add_argument('-o', '--output',
                       required=True,
                       help='输出文件|Output file')

    # 可选参数|Optional arguments
    parser.add_argument('--no-compress',
                       action='store_true',
                       help='不压缩输出文件|Do not compress output files')

    args = parser.parse_args()

    # 创建提取器并运行|Create extractor and run
    extractor = ReadsExtractor(
        mapping_file=args.mapping,
        fastq_file=args.input,
        output_file=args.output,
        compress_output=not args.no_compress
    )

    sys.exit(extractor.run())


if __name__ == "__main__":
    main()
