"""
ENA下载工具主程序模块|ENA Downloader Main Module
"""

import argparse
import sys
from .config import DownloadConfig
from .utils import DownloadLogger, CommandRunner, check_dependencies
from .metadata import MetadataDownloader
from .downloader import FastqDownloader
from .results import ResultsSummary

class ENADownloader:
    """ENA下载工具主类|Main ENA Downloader Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = DownloadConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = DownloadLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化各个处理器|Initialize processors
        self.metadata_downloader = MetadataDownloader(self.config, self.logger)
        self.fastq_downloader = FastqDownloader(self.config, self.logger, self.cmd_runner)
        self.results_summary = ResultsSummary(self.config, self.logger)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def download_metadata(self):
        """仅下载元数据|Download metadata only"""
        self.logger.info("开始元数据下载流程|Starting metadata download pipeline")

        try:
            # 下载元数据|Download metadata
            metadata_file = self.metadata_downloader.download_metadata()

            if metadata_file:
                self.logger.info("元数据下载完成|Metadata download completed")

                # 生成汇总报告|Generate summary report
                self.results_summary.generate_summary(metadata_file)
                return metadata_file
            else:
                self.logger.error("元数据下载失败|Metadata download failed")
                return None

        except Exception as e:
            self.logger.error(f"元数据下载过程中发生错误|Error during metadata download: {str(e)}")
            return None

    def run_full_pipeline(self):
        """运行完整的下载流程|Run full download pipeline"""
        self.logger.info("开始完整下载流程|Starting full download pipeline")

        try:
            # 检查依赖|Check dependencies (仅当需要下载FASTQ时)
            if self.config.method in ["save", "run"]:
                if not self.check_dependencies():
                    self.logger.error("依赖检查失败，无法继续|Dependency check failed, cannot continue")
                    return False

            # 下载元数据|Download metadata
            self.logger.info("步骤1/3: 下载元数据|Step 1/3: Download metadata")
            metadata_file = self.metadata_downloader.download_metadata()

            if not metadata_file:
                self.logger.error("元数据下载失败，流程终止|Metadata download failed, pipeline terminated")
                return False

            # 提取下载链接|Extract download links
            self.logger.info("步骤2/3: 提取下载链接|Step 2/3: Extract download links")
            download_links = self.metadata_downloader.get_download_links(metadata_file)

            if not download_links:
                self.logger.warning("未找到下载链接|No download links found")
                # 仍然生成汇总报告|Still generate summary report
                self.results_summary.generate_summary(metadata_file, download_links_count=0)
                return True

            # 处理下载|Process downloads
            self.logger.info("步骤3/3: 处理下载|Step 3/3: Process downloads")
            script_file = None

            if self.config.method == "save":
                # 生成下载脚本|Generate download script
                script_file = self.fastq_downloader.generate_download_script(download_links)
                if script_file:
                    self.logger.info(f"下载脚本已生成|Download script generated: {script_file}")
                    print(f'\n\033[32m下载脚本已生成|Download script generated: {script_file.name}\033[0m')
                    print(f'请运行以下命令开始下载|Please run the following command to start download:')
                    print(f'bash {script_file.name}')
                else:
                    self.logger.error("下载脚本生成失败|Download script generation failed")

            elif self.config.method == "run":
                # 直接执行下载|Execute downloads directly
                success = self.fastq_downloader.execute_downloads(download_links)
                if success:
                    self.logger.info("直接下载完成|Direct download completed")
                else:
                    self.logger.error("直接下载失败|Direct download failed")

            # 生成汇总报告|Generate summary report
            self.results_summary.generate_summary(
                metadata_file,
                script_file,
                len(download_links)
            )

            self.logger.info("完整下载流程结束|Full download pipeline completed")
            return True

        except Exception as e:
            self.logger.error(f"下载流程中发生错误|Error in download pipeline: {str(e)}")
            return False

def create_argument_parser():
    """创建命令行参数解析器|Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='ENA数据下载工具|ENA Data Download Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -a PRJNA661210
        '''
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '--accession', '-a',
        required=True,
        help='ENA项目编号|ENA accession number'
    )

    # 输出设置|Output settings
    parser.add_argument(
        '--output-dir', '-o',
        help='输出目录|Output directory'
    )

    parser.add_argument(
        '--create-dir', '-d',
        action='store_true',
        help='创建专门的输出目录|Create dedicated output directory'
    )

    parser.add_argument(
        '--metadata-format', '-f',
        choices=['tsv', 'csv', 'xlsx'],
        default='tsv',
        help='元数据文件格式|Metadata file format'
    )

    # 下载协议设置|Download protocol settings
    parser.add_argument(
        '--protocol', '-p',
        choices=['ftp', 'aspera'],
        default='ftp',
        help='下载协议类型|Download protocol type'
    )

    parser.add_argument(
        '--aspera-key', '-k',
        help='Aspera私钥路径|Path to aspera private key'
    )

    parser.add_argument(
        '--method', '-m',
        choices=['save', 'run'],
        default='save',
        help='执行模式|Execution mode'
    )

    # 特殊模式|Special modes
    parser.add_argument(
        '--metadata-only', '-M',
        action='store_true',
        help='仅下载元数据，不处理FASTQ文件|Only download metadata, do not process FASTQ files'
    )

    # 高级选项|Advanced options
    parser.add_argument(
        '--fields', '-F',
        nargs='+',
        help='自定义元数据字段|Custom metadata fields'
    )

    parser.add_argument(
        '--max-retries', '-r',
        type=int,
        default=3,
        help='API请求最大重试次数|Maximum API request retries'
    )

    return parser

def main():
    """主函数|Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()

    try:
        # 准备配置参数|Prepare configuration parameters
        config_kwargs = {
            'accession': args.accession,
            'output_dir': args.output_dir,
            'create_dir': args.create_dir,
            'protocol': args.protocol,
            'method': args.method,
            'aspera_key': args.aspera_key,
            'fields': args.fields,
            'metadata_format': args.metadata_format,
            'max_retries': args.max_retries
        }

        # 创建下载器|Create downloader
        downloader = ENADownloader(**config_kwargs)

        # 选择执行模式|Choose execution mode
        if args.metadata_only:
            # 仅下载元数据模式|Metadata-only mode
            result = downloader.download_metadata()
            if result:
                print(f'\n\033[32m元数据下载完成|Metadata download completed\033[0m')
                print(f'元数据文件|Metadata file: {result}')
            else:
                print('\n\033[31m元数据下载失败|Metadata download failed\033[0m')
                sys.exit(1)
        else:
            # 完整流程模式|Full pipeline mode
            success = downloader.run_full_pipeline()
            if success:
                print(f'\n\033[32m流程执行完成|Pipeline execution completed\033[0m')
            else:
                print('\n\033[31m流程执行失败|Pipeline execution failed\033[0m')
                sys.exit(1)

    except KeyboardInterrupt:
        print('\n\033[33m用户中断操作|User interrupted operation\033[0m')
        sys.exit(1)
    except Exception as e:
        print(f'\n\033[31m发生错误|Error occurred: {str(e)}\033[0m')
        sys.exit(1)

if __name__ == '__main__':
    main()
