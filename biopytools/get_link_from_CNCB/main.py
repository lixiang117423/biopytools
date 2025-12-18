"""
CNCB链接提取器主程序模块 | CNCB Link Extractor Main Module
"""

import os
import sys
import time
import logging
from typing import Dict, List, Tuple, Optional
from pathlib import Path

from .config import CNCBConfig
from .utils import (
    CNCBLogger, InputFileParser, FTPConnectionManager,
    PathCache, FileDownloader
)
from .ftp_searcher import CNCBFTPSearcher


class CNCLinkExtractor:
    """CNCB链接提取器主类 | Main CNCB Link Extractor Class"""

    def __init__(self, **kwargs):
        """初始化链接提取器 | Initialize link extractor"""
        self.config = CNCBConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = CNCBLogger(
            log_file=self.config.log_file,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化组件 | Initialize components
        self.path_cache = PathCache()
        self.ftp_manager = FTPConnectionManager(
            host=self.config.ftp_host,
            timeout=self.config.ftp_timeout,
            retry_attempts=self.config.retry_attempts
        )

        self.logger.info("🚀 CNCB链接提取器已初始化 | CNCB Link Extractor initialized")
        self.logger.info(f"📁 输入文件 | Input file: {self.config.input_file}")
        self.logger.info(f"💾 输出文件 | Output files: {self.config.get_output_info()}")

    def run_extraction(self) -> bool:
        """运行链接提取流程 | Run link extraction pipeline"""
        try:
            start_time = time.time()

            # 读取和验证输入文件 | Read and validate input file
            projects = self._read_input_file()
            if not projects:
                return False

            # 连接FTP服务器 | Connect to FTP server
            if not self._connect_ftp():
                return False

            # 提取链接 | Extract links
            success_urls, failed_ids = self._extract_links(projects)

            # 保存结果 | Save results
            self._save_results(success_urls, failed_ids, projects)

            # 生成报告 | Generate report
            self._generate_summary_report(projects, success_urls, failed_ids)

            # 输出统计信息 | Output statistics
            self._log_statistics(projects, success_urls, failed_ids, start_time)

            return True

        except Exception as e:
            self.logger.error(f"💥 链接提取过程中发生错误 | Error during link extraction: {e}")
            return False

        finally:
            self._cleanup()

    def _read_input_file(self) -> Optional[Dict[str, List[str]]]:
        """读取输入文件 | Read input file"""
        self.logger.info("📖 正在读取输入文件 | Reading input file...")

        # 验证文件格式 | Validate file format
        is_valid, message = InputFileParser.validate_input_file(self.config.input_file)
        if not is_valid:
            self.logger.error(f"❌ {message}")
            return None

        self.logger.info(f"✅ {message}")

        # 解析文件内容 | Parse file content
        projects = InputFileParser.read_and_group_by_project(self.config.input_file)
        if not projects:
            self.logger.error("❌ 无法解析输入文件 | Failed to parse input file")
            return None

        return projects

    def _connect_ftp(self) -> bool:
        """连接FTP服务器 | Connect to FTP server"""
        self.logger.info("🌐 正在连接到CNCB FTP服务器 | Connecting to CNCB FTP server...")

        if not self.ftp_manager.connect():
            self.logger.error("❌ FTP连接失败 | FTP connection failed")
            return False

        # 测试访问权限 | Test access permissions
        ftp_client = self.ftp_manager.get_ftp()
        searcher = CNCBFTPSearcher(
            ftp_client=ftp_client,
            base_dirs=self.config.base_dirs,
            path_cache=self.path_cache
        )

        if not searcher.test_ftp_access():
            self.logger.error("❌ FTP访问权限不足 | Insufficient FTP access permissions")
            return False

        # 显示服务器信息 | Show server information
        server_info = searcher.get_server_info()
        if server_info.get("base_dirs_accessible"):
            self.logger.info(f"📂 可访问的基础目录 | Accessible base directories: {server_info['base_dirs_accessible']}")

        return True

    def _extract_links(self, projects: Dict[str, List[str]]) -> Tuple[List[str], List[Tuple[str, str]]]:
        """提取链接 | Extract links"""
        self.logger.info("🔍 开始提取下载链接 | Starting link extraction...")

        all_found_urls = []
        all_failed_ids = []
        archive_map = self.config.get_archive_map()
        filename_templates = self.config.get_filename_templates()

        ftp_client = self.ftp_manager.get_ftp()
        searcher = CNCBFTPSearcher(
            ftp_client=ftp_client,
            base_dirs=self.config.base_dirs,
            path_cache=self.path_cache
        )

        total_projects = len(projects)
        project_count = 0

        for project_id, run_id_list in projects.items():
            project_count += 1
            self.logger.info("")
            self.logger.info(f"📁 [{project_count}/{total_projects}] 处理项目 | Processing project: {project_id} ({len(run_id_list)} 个ID)")

            for run_id in run_id_list:
                self.logger.info(f"  🔄 处理Run ID | Processing Run ID: {run_id}")

                # 确定数据库类型 | Determine archive type
                prefix = run_id[:3]
                archive = archive_map.get(prefix)
                if not archive:
                    self.logger.error(f"    ❌ 未知前缀 | Unknown prefix: {prefix}")
                    all_failed_ids.append((project_id, run_id))
                    continue

                # 查找基础路径 | Find base path
                base_path = searcher.find_base_path_for_id(archive, run_id)
                if not base_path:
                    self.logger.warning(f"    ⚠️ 未找到路径 | Path not found for: {run_id}")
                    all_failed_ids.append((project_id, run_id))
                    continue

                # 查找文件 | Find files
                found_urls = searcher.find_files_for_run(
                    run_id=run_id,
                    base_path=base_path,
                    filename_templates=filename_templates
                )

                if found_urls:
                    self.logger.info(f"    ✅ 找到 {len(found_urls)} 个文件 | Found {len(found_urls)} files")
                    all_found_urls.extend(found_urls)
                else:
                    self.logger.warning(f"    ⚠️ 未找到文件 | No files found for: {run_id}")
                    all_failed_ids.append((project_id, run_id))

        return all_found_urls, all_failed_ids

    def _save_results(self, success_urls: List[str], failed_ids: List[Tuple[str, str]], projects: Dict[str, List[str]]):
        """保存结果 | Save results"""
        self.logger.info("💾 正在保存结果 | Saving results...")

        # 保存成功的链接 | Save successful links
        try:
            with open(self.config.output_file, 'w', encoding='utf-8') as f:
                for url in sorted(success_urls):
                    f.write(url + '\n')
            self.logger.info(f"✅ 链接已保存 | Links saved to: {self.config.output_file}")
        except Exception as e:
            self.logger.error(f"❌ 保存链接失败 | Failed to save links: {e}")

        # 保存失败的ID | Save failed IDs
        if failed_ids:
            try:
                with open(self.config.failed_file, 'w', encoding='utf-8') as f:
                    for project, run_id in sorted(failed_ids):
                        f.write(f"{project}\t{run_id}\n")
                self.logger.info(f"📝 失败记录已保存 | Failed records saved to: {self.config.failed_file}")
            except Exception as e:
                self.logger.error(f"❌ 保存失败记录失败 | Failed to save failed records: {e}")

        # 生成下载脚本 | Generate download script
        if self.config.generate_download_script and success_urls:
            script_dir = os.path.dirname(os.path.abspath(self.config.output_file))
            script_path = os.path.join(script_dir, self.config.download_script)

            if FileDownloader.generate_download_script(
                urls=success_urls,
                script_file=script_path,
                make_executable=self.config.script_executable
            ):
                self.logger.info(f"📜 下载脚本已生成 | Download script generated: {script_path}")
                self.logger.info("💡 使用方法 | Usage:")
                self.logger.info(f"   bash {self.config.download_script}")
                self.logger.info("   # 或者后台运行 | Or run in background:")
                self.logger.info(f"   nohup bash {self.config.download_script} &")

    def _generate_summary_report(self, projects: Dict[str, List[str]], success_urls: List[str], failed_ids: List[Tuple[str, str]]):
        """生成总结报告 | Generate summary report"""
        output_dir = os.path.dirname(os.path.abspath(self.config.output_file))
        report_file = FileDownloader.generate_summary_report(
            projects=projects,
            successful_urls=success_urls,
            failed_ids=failed_ids,
            output_dir=output_dir
        )

        if report_file and os.path.exists(report_file):
            self.logger.info(f"📊 总结报告已生成 | Summary report generated: {report_file}")

    def _log_statistics(self, projects: Dict[str, List[str]], success_urls: List[str], failed_ids: List[Tuple[str, str]], start_time: float):
        """记录统计信息 | Log statistics"""
        total_ids = sum(len(ids) for ids in projects.values())
        successful_ids = total_ids - len(failed_ids)
        success_rate = (successful_ids / total_ids * 100) if total_ids > 0 else 0

        # 缓存统计 | Cache statistics
        cache_stats = self.path_cache.get_stats()

        # 时间统计 | Time statistics
        elapsed_time = time.time() - start_time

        self.logger.info("")
        self.logger.info("📊 处理统计 | Processing Statistics:")
        self.logger.info(f"  📁 总项目数 | Total Projects: {len(projects)}")
        self.logger.info(f"  🧬 总Run ID数 | Total Run IDs: {total_ids}")
        self.logger.info(f"  ✅ 成功ID数 | Successful IDs: {successful_ids}")
        self.logger.info(f"  ❌ 失败ID数 | Failed IDs: {len(failed_ids)}")
        self.logger.info(f"  📈 成功率 | Success Rate: {success_rate:.2f}%")
        self.logger.info(f"  📄 总链接数 | Total URLs: {len(success_urls)}")
        self.logger.info(f"  ⏱️  总耗时 | Total Time: {elapsed_time:.2f} seconds")
        self.logger.info(f"  🎯 缓存命中率 | Cache Hit Rate: {cache_stats['hit_rate_percent']}%")

        # 输出文件信息 | Output file information
        self.logger.info("")
        self.logger.info("📁 输出文件 | Output Files:")
        self.logger.info(f"  📄 链接文件 | Links File: {self.config.output_file}")
        if failed_ids:
            self.logger.info(f"  ❌ 失败记录 | Failed Records: {self.config.failed_file}")
        if self.config.generate_download_script:
            self.logger.info(f"  📜 下载脚本 | Download Script: {self.config.download_script}")

    def _cleanup(self):
        """清理资源 | Cleanup resources"""
        if self.ftp_manager:
            self.ftp_manager.disconnect()

        self.logger.info("🧹 资源清理完成 | Resource cleanup completed")
        self.logger.info("🎉 CNCB链接提取完成 | CNCB link extraction completed!")


def main():
    """主函数 | Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="从CNCB数据库批量获取测序数据下载链接",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 必需参数 | Required parameters
    parser.add_argument("input_file",
                       help="输入文件路径 (ProjectID和RunID，Tab分隔) | Input file path (ProjectID and Run ID, tab separated)")

    # 输出配置 | Output configuration
    parser.add_argument("-o", "--output",
                       help="输出文件路径 | Output file path (default: [input]_links.txt)")
    parser.add_argument("-f", "--failed",
                       help="失败记录文件路径 | Failed records file path (default: [input]_failed.txt)")
    parser.add_argument("--download-script",
                       help="下载脚本文件名 | Download script filename (default: download.sh)")

    # FTP配置 | FTP configuration
    parser.add_argument("--ftp-host",
                       default="download2.cncb.ac.cn",
                       help="FTP服务器地址 | FTP server host (default: download2.cncb.ac.cn)")
    parser.add_argument("--ftp-timeout", type=int, default=60,
                       help="FTP连接超时时间 | FTP connection timeout in seconds (default: 60)")

    # 性能配置 | Performance configuration
    parser.add_argument("--retry-attempts", type=int, default=3,
                       help="FTP连接重试次数 | FTP connection retry attempts (default: 3)")
    parser.add_argument("--threads", type=int, default=4,
                       help="并发线程数 | Concurrent threads (default: 4)")

    # 日志配置 | Logging configuration
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="详细输出模式 | Verbose output mode")
    parser.add_argument("--log-file",
                       help="日志文件路径 | Log file path")

    # 输出选项 | Output options
    parser.add_argument("--no-download-script", action="store_true",
                       help="不生成下载脚本 | Don't generate download script")
    parser.add_argument("--no-executable", action="store_true",
                       help="不设置脚本执行权限 | Don't make script executable")

    args = parser.parse_args()

    try:
        # 创建提取器 | Create extractor
        extractor = CNCLinkExtractor(
            input_file=args.input_file,
            output_file=args.output,
            failed_file=args.failed,
            download_script=args.download_script,
            ftp_host=args.ftp_host,
            ftp_timeout=args.ftp_timeout,
            retry_attempts=args.retry_attempts,
            max_threads=args.threads,
            verbose=args.verbose,
            log_file=args.log_file,
            generate_download_script=not args.no_download_script,
            script_executable=not args.no_executable
        )

        # 运行提取 | Run extraction
        success = extractor.run_extraction()

        # 退出 | Exit
        sys.exit(0 if success else 1)

    except Exception as e:
        logging.error(f"💥 程序执行失败 | Program execution failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()