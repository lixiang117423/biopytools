"""
CNCB链接提取器主程序模块|CNCB Link Extractor Main Module
"""

import os
import re
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
from .ena_searcher import ENALinkSearcher
from .gsa_searcher import GSAHTTPSearcher


class CNCLinkExtractor:
    """CNCB链接提取器主类|Main CNCB Link Extractor Class"""

    def __init__(self, **kwargs):
        """初始化链接提取器|Initialize link extractor"""
        self.config = CNCBConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = CNCBLogger(
            log_file=self.config.log_file,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化组件|Initialize components
        self.path_cache = PathCache()
        self.ftp_manager = FTPConnectionManager(
            host=self.config.ftp_host,
            timeout=self.config.ftp_timeout,
            retry_attempts=self.config.retry_attempts,
            logger=self.logger
        )
        # GSA HTTP搜索器(CRR前缀,不依赖FTP)|GSA HTTP searcher (CRR prefix, FTP-independent)
        self.gsa_searcher = GSAHTTPSearcher(
            logger=self.logger,
            retry_attempts=self.config.retry_attempts
        )

        self.logger.info("CNCB链接提取器已初始化|CNCB Link Extractor initialized")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出文件|Output files: {self.config.get_output_info()}")

        # ENA回退命中计数|ENA fallback hit counter
        self.ena_found_count = 0

    def run_extraction(self) -> bool:
        """运行链接提取流程|Run link extraction pipeline"""
        try:
            start_time = time.time()

            # 读取和验证输入文件|Read and validate input file
            projects = self._read_input_file()
            if not projects:
                return False

            # 连接FTP服务器|Connect to FTP server
            # 连接失败不退出:INSDC ID交给ENA回退,CRR(GSA)不依赖FTP
            # |Failure does not abort: INSDC IDs fall back to ENA; CRR (GSA) is FTP-independent
            ftp_connected = self._connect_ftp()
            if not ftp_connected:
                if self.config.ena_fallback:
                    self.logger.warning(
                        "FTP连接失败,所有INSDC ID将交给ENA回退查询|FTP connection failed, "
                        "all INSDC IDs will be attempted via ENA fallback"
                    )
                else:
                    self.logger.warning(
                        "FTP连接失败且ENA回退已关闭,所有INSDC ID将记入失败列表|FTP connection failed "
                        "and ENA fallback disabled, all INSDC IDs will be recorded as failed"
                    )

            # 提取链接|Extract links
            success_urls, failed_ids = self._extract_links(projects, ftp_connected)

            # CNCB未找到的ID回退ENA查询|Fall back to ENA for IDs missed by CNCB
            if self.config.ena_fallback and failed_ids:
                success_urls, failed_ids = self._fallback_ena(success_urls, failed_ids)

            # 保存结果|Save results
            self._save_results(success_urls, failed_ids, projects)

            # 生成报告|Generate report
            self._generate_summary_report(projects, success_urls, failed_ids)

            # 输出统计信息|Output statistics
            self._log_statistics(projects, success_urls, failed_ids, start_time)

            return True

        except Exception as e:
            self.logger.error(f"链接提取过程中发生错误|Error during link extraction: {e}")
            return False

        finally:
            self._cleanup()

    def _read_input_file(self) -> Optional[Dict[str, List[str]]]:
        """读取输入文件|Read input file"""
        self.logger.info("正在读取输入文件|Reading input file...")

        # 验证文件格式|Validate file format
        is_valid, message = InputFileParser.validate_input_file(self.config.input_file)
        if not is_valid:
            self.logger.error(f"{message}")
            return None

        self.logger.info(f"{message}")

        # 解析文件内容|Parse file content
        projects = InputFileParser.read_and_group_by_project(self.config.input_file, logger=self.logger)
        if not projects:
            self.logger.error("无法解析输入文件|Failed to parse input file")
            return None

        return projects

    def _connect_ftp(self) -> bool:
        """连接FTP服务器|Connect to FTP server"""
        self.logger.info("正在连接到CNCB FTP服务器|Connecting to CNCB FTP server...")

        if not self.ftp_manager.connect():
            self.logger.error("FTP连接失败|FTP connection failed")
            return False

        # 测试访问权限|Test access permissions
        ftp_client = self.ftp_manager.get_ftp()
        searcher = CNCBFTPSearcher(
            ftp_client=ftp_client,
            base_dirs=self.config.base_dirs,
            path_cache=self.path_cache,
            logger=self.logger
        )

        if not searcher.test_ftp_access():
            self.logger.error("FTP访问权限不足|Insufficient FTP access permissions")
            return False

        # 显示服务器信息|Show server information
        server_info = searcher.get_server_info()
        if server_info.get("base_dirs_accessible"):
            self.logger.info(f"可访问的基础目录|Accessible base directories: {server_info['base_dirs_accessible']}")

        return True

    def _extract_links(self, projects: Dict[str, List[str]],
                       ftp_connected: bool) -> Tuple[List[str], List[Tuple[str, str]]]:
        """
        提取链接|Extract links

        CRR(GSA原生)走HTTPS通道,不依赖FTP;SRR/ERR/DRR走FTP INSDC镜像
        |CRR (GSA-native) uses the HTTPS channel, FTP-independent; SRR/ERR/DRR use the FTP INSDC mirror
        """
        self.logger.info("开始提取下载链接|Starting link extraction...")

        all_found_urls = []
        all_failed_ids = []
        archive_map = self.config.get_archive_map()
        filename_templates = self.config.get_filename_templates()

        searcher = None
        if ftp_connected:
            ftp_client = self.ftp_manager.get_ftp()
            searcher = CNCBFTPSearcher(
                ftp_client=ftp_client,
                base_dirs=self.config.base_dirs,
                path_cache=self.path_cache,
                logger=self.logger
            )

        total_projects = len(projects)
        project_count = 0

        for project_id, run_id_list in projects.items():
            project_count += 1
            self.logger.info("")
            self.logger.info(f"[{project_count}/{total_projects}] 处理项目|Processing project: {project_id} ({len(run_id_list)} 个IDs|IDs)")

            for run_id in run_id_list:
                self.logger.info(f"处理Run ID|Processing Run ID: {run_id}")

                # 确定数据库类型|Determine archive type
                prefix = run_id[:3]
                archive = archive_map.get(prefix)

                # GSA原生Run ID走HTTPS通道|GSA-native run IDs use the HTTPS channel
                if archive == "GSA":
                    self._handle_gsa_run(project_id, run_id, all_found_urls, all_failed_ids)
                    continue

                # FTP未连接时INSDC ID直接记失败,留待ENA回退|Without FTP, INSDC IDs fail now for ENA fallback
                if not ftp_connected:
                    all_failed_ids.append((project_id, run_id))
                    continue

                if not archive:
                    self.logger.error(f"未知前缀|Unknown prefix: {prefix}")
                    all_failed_ids.append((project_id, run_id))
                    continue

                # 查找基础路径|Find base path
                base_path = searcher.find_base_path_for_id(archive, run_id)
                if not base_path:
                    self.logger.warning(f"未找到路径|Path not found for: {run_id}")
                    all_failed_ids.append((project_id, run_id))
                    continue

                # 查找文件|Find files
                found_urls = searcher.find_files_for_run(
                    run_id=run_id,
                    base_path=base_path,
                    filename_templates=filename_templates
                )

                if found_urls:
                    self.logger.info(f"找到 {len(found_urls)} 个文件|Found {len(found_urls)} files")
                    all_found_urls.extend(found_urls)
                else:
                    self.logger.warning(f"未找到文件|No files found for: {run_id}")
                    all_failed_ids.append((project_id, run_id))

        return all_found_urls, all_failed_ids

    def _handle_gsa_run(self, project_id: str, run_id: str,
                        found_urls: List[str], failed_ids: List[Tuple[str, str]]):
        """
        处理GSA原生Run ID|Handle a GSA-native run ID

        项目列必须是CRA编号;PRJCA等BioProject编号无法自动映射到GSA目录
        |The project column must be the CRA accession; BioProject IDs like PRJCA
        cannot be auto-mapped to GSA directories
        """
        if not re.match(r'^CRA\d+$', project_id):
            self.logger.warning(
                f"CRR是GSA原生Run ID,项目列请用CRA编号(如CRA010060),PRJCA无法自动映射"
                f"|CRR is a GSA-native run ID; use its CRA accession in the project column, "
                f"PRJCA cannot be auto-mapped: {project_id} -> {run_id}"
            )
            failed_ids.append((project_id, run_id))
            return

        found = self.gsa_searcher.search_files(project_id, run_id)
        if found:
            found_urls.extend(found)
        else:
            failed_ids.append((project_id, run_id))

    def _fallback_ena(self, success_urls: List[str],
                      failed_ids: List[Tuple[str, str]]) -> Tuple[List[str], List[Tuple[str, str]]]:
        """
        对CNCB未找到的ID批量回退ENA查询|Batch-query ENA for IDs missed by CNCB

        Args:
            success_urls: 已有的成功链接|Existing successful URLs
            failed_ids: CNCB未找到的(项目ID, Run ID)列表|CNCB-missed (project, run) pairs

        Returns:
            (合并后的链接列表, 仍失败的ID列表)|(merged URLs, still-failed IDs)
        """
        # 只查INSDC前缀(SRR/ERR/DRR);CRR是GSA原生ID,ENA没有,跳过以免误导
        # |Query only INSDC prefixes; CRR is GSA-native and absent from ENA, skip to avoid noise
        ena_candidates = [(p, r) for p, r in failed_ids if r[:3] in ("SRR", "ERR", "DRR")]
        non_ena_failed = [(p, r) for p, r in failed_ids if r[:3] not in ("SRR", "ERR", "DRR")]

        self.logger.info(f"开始ENA回退查询|Starting ENA fallback query: {len(ena_candidates)} 个IDs|IDs")
        if not ena_candidates:
            return success_urls, failed_ids

        searcher = ENALinkSearcher(
            logger=self.logger,
            retry_attempts=self.config.retry_attempts
        )
        run_ids = [run_id for _, run_id in ena_candidates]
        links_by_run = searcher.search_links(run_ids)

        still_failed = list(non_ena_failed)
        for project_id, run_id in ena_candidates:
            urls = links_by_run.get(run_id)
            if urls:
                self.logger.info(f"ENA回退命中|ENA fallback hit: {run_id} -> {len(urls)} 个文件|files")
                success_urls.extend(urls)
                self.ena_found_count += 1
            else:
                still_failed.append((project_id, run_id))

        self.logger.info(
            f"ENA回退完成|ENA fallback completed: 命中|hit {self.ena_found_count}, "
            f"未命中|missed {len(still_failed)}"
        )
        return success_urls, still_failed

    def _save_results(self, success_urls: List[str], failed_ids: List[Tuple[str, str]], projects: Dict[str, List[str]]):
        """保存结果|Save results"""
        self.logger.info("正在保存结果|Saving results...")

        # 保存成功的链接|Save successful links
        try:
            with open(self.config.output_file, 'w', encoding='utf-8') as f:
                for url in sorted(success_urls):
                    f.write(url + '\n')
            self.logger.info(f"链接已保存|Links saved to: {self.config.output_file}")
        except Exception as e:
            self.logger.error(f"保存链接失败|Failed to save links: {e}")

        # 保存失败的ID|Save failed IDs
        # 无条件写入(空文件也写):重跑全部成功时会清掉上次运行残留的陈旧失败记录
        # |Always write (even empty): a fully successful rerun clears stale failed records
        try:
            with open(self.config.failed_file, 'w', encoding='utf-8') as f:
                for project, run_id in sorted(failed_ids):
                    f.write(f"{project}\t{run_id}\n")
            self.logger.info(
                f"失败记录已保存|Failed records saved to: {self.config.failed_file} "
                f"({len(failed_ids)} 条|records)"
            )
        except Exception as e:
            self.logger.error(f"保存失败记录失败|Failed to save failed records: {e}")

        # 生成下载脚本|Generate download script
        if self.config.generate_download_script and success_urls:
            script_dir = os.path.dirname(os.path.abspath(self.config.output_file))
            script_path = os.path.join(script_dir, self.config.download_script)

            if FileDownloader.generate_download_script(
                urls=success_urls,
                script_file=script_path,
                make_executable=self.config.script_executable,
                logger=self.logger
            ):
                self.logger.info(f"下载脚本已生成|Download script generated: {script_path}")
                self.logger.info("使用方法|Usage:")
                self.logger.info(f"   bash {self.config.download_script}")
                self.logger.info("   # 或者后台运行|Or run in background:")
                self.logger.info(f"   nohup bash {self.config.download_script} &")

    def _generate_summary_report(self, projects: Dict[str, List[str]], success_urls: List[str], failed_ids: List[Tuple[str, str]]):
        """生成总结报告|Generate summary report"""
        output_dir = os.path.dirname(os.path.abspath(self.config.output_file))
        report_file = FileDownloader.generate_summary_report(
            projects=projects,
            successful_urls=success_urls,
            failed_ids=failed_ids,
            output_dir=output_dir,
            logger=self.logger
        )

        if report_file and os.path.exists(report_file):
            self.logger.info(f" 总结报告已生成|Summary report generated: {report_file}")

    def _log_statistics(self, projects: Dict[str, List[str]], success_urls: List[str], failed_ids: List[Tuple[str, str]], start_time: float):
        """记录统计信息|Log statistics"""
        total_ids = sum(len(ids) for ids in projects.values())
        successful_ids = total_ids - len(failed_ids)
        success_rate = (successful_ids / total_ids * 100) if total_ids > 0 else 0

        # 缓存统计|Cache statistics
        cache_stats = self.path_cache.get_stats()

        # 时间统计|Time statistics
        elapsed_time = time.time() - start_time

        self.logger.info("")
        self.logger.info("处理统计|Processing Statistics:")
        self.logger.info(f"   总项目数|Total Projects: {len(projects)}")
        self.logger.info(f"   总Run ID数|Total Run IDs: {total_ids}")
        self.logger.info(f"   成功ID数|Successful IDs: {successful_ids}")
        self.logger.info(f"   失败ID数|Failed IDs: {len(failed_ids)}")
        self.logger.info(f"   成功率|Success Rate: {success_rate:.2f}%")
        self.logger.info(f"   总链接数|Total URLs: {len(success_urls)}")
        if self.config.ena_fallback and self.ena_found_count:
            self.logger.info(f"   ENA回退命中数|ENA Fallback Hits: {self.ena_found_count}")
        self.logger.info(f"   总耗时|Total Time: {elapsed_time:.2f} seconds")
        self.logger.info(f"   缓存命中率|Cache Hit Rate: {cache_stats['hit_rate_percent']}%")

        # 输出文件信息|Output file information
        self.logger.info("")
        self.logger.info("输出文件|Output Files:")
        self.logger.info(f"   链接文件|Links File: {self.config.output_file}")
        if failed_ids:
            self.logger.info(f"   失败记录|Failed Records: {self.config.failed_file}")
        if self.config.generate_download_script:
            self.logger.info(f"   下载脚本|Download Script: {self.config.download_script}")

    def _cleanup(self):
        """清理资源|Cleanup resources"""
        if self.ftp_manager:
            self.ftp_manager.disconnect()

        self.logger.info("资源清理完成|Resource cleanup completed")
        self.logger.info("CNCB链接提取完成|CNCB link extraction completed!")


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="从CNCB数据库批量获取测序数据下载链接",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument("input_file",
                       help="输入文件路径 (ProjectID和RunID，Tab分隔;GSA项目(CRR前缀)项目列用CRA编号)|"
                            "Input file path (ProjectID and Run ID, tab separated; "
                            "use the CRA accession in the project column for GSA projects with CRR runs)")

    # 输出配置|Output configuration
    parser.add_argument("-o", "--output",
                       help="输出文件路径|Output file path (default: [input]_links.txt)")
    parser.add_argument("-f", "--failed",
                       help="失败记录文件路径|Failed records file path (default: [input]_failed.txt)")
    parser.add_argument("--download-script",
                       help="下载脚本文件名|Download script filename (default: download.sh)")

    # FTP配置|FTP configuration
    parser.add_argument("--ftp-host",
                       default="download2.cncb.ac.cn",
                       help="FTP服务器地址|FTP server host (default: download2.cncb.ac.cn)")
    parser.add_argument("--ftp-timeout", type=int, default=60,
                       help="FTP连接超时时间|FTP connection timeout in seconds (default: 60)")

    # 性能配置|Performance configuration
    parser.add_argument("--retry-attempts", type=int, default=3,
                       help="FTP连接重试次数|FTP connection retry attempts (default: 3)")

    # 日志配置|Logging configuration
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="详细输出模式|Verbose output mode")
    parser.add_argument("--log-file",
                       help="日志文件路径|Log file path")

    # 输出选项|Output options
    parser.add_argument("--no-download-script", action="store_true",
                       help="不生成下载脚本|Don't generate download script")
    parser.add_argument("--no-executable", action="store_true",
                       help="不设置脚本执行权限|Don't make script executable")

    # ENA回退选项|ENA fallback options
    parser.add_argument("--no-ena-fallback", action="store_true",
                       help="关闭ENA回退查询(纯CNCB模式)|Disable ENA fallback query (pure CNCB mode)")

    args = parser.parse_args()

    try:
        # 创建提取器|Create extractor
        extractor = CNCLinkExtractor(
            input_file=args.input_file,
            output_file=args.output,
            failed_file=args.failed,
            download_script=args.download_script,
            ftp_host=args.ftp_host,
            ftp_timeout=args.ftp_timeout,
            retry_attempts=args.retry_attempts,
            verbose=args.verbose,
            log_file=args.log_file,
            generate_download_script=not args.no_download_script,
            script_executable=not args.no_executable,
            ena_fallback=not args.no_ena_fallback
        )

        # 运行提取|Run extraction
        success = extractor.run_extraction()

        # 退出|Exit
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"程序执行失败|Program execution failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()