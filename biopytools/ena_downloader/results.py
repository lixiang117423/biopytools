"""
 ENA下载结果汇总模块|ENA Download Results Summary Module
"""

from pathlib import Path
from datetime import datetime

def format_number(num: int) -> str:
    """ 大数字格式化: 超过1百万用M单位保留2位小数|Format large numbers: >= 1 million as M with 2 decimals"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)

class ResultsSummary:
    """ 结果汇总生成器|Results Summary Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary(self, accession: str, metadata_file: Path, script_file: Path = None, download_links_count: int = 0):
        """ 生成单个编号的汇总报告|Generate summary report for one accession"""
        summary_file = self.config.output_path / f'{accession}.download_summary.txt'

        try:
            with summary_file.open('w', encoding='utf-8') as f:
                f.write("="*60 + "\n")
                f.write(f" ENA数据下载汇总报告|ENA Data Download Summary Report: {accession}\n")
                f.write("="*60 + "\n\n")

                #  基本信息|Basic information
                f.write(" 基本信息|Basic Information:\n")
                f.write(f"  -  编号|Accession: {accession}\n")
                f.write(f"  -  下载协议|Protocol: {self.config.protocol}\n")
                f.write(f"  -  执行方式|Method: {self.config.method}\n")
                f.write(f"  -  输出目录|Output Directory: {self.config.output_dir}\n")
                f.write(f"  -  生成时间|Generated Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("\n")

                #  元数据信息|Metadata information
                f.write(" 元数据信息|Metadata Information:\n")
                if metadata_file and metadata_file.exists():
                    f.write(f"  -  元数据文件|Metadata File: {metadata_file.name}\n")
                    f.write(f"  -  文件大小|File Size: {format_number(metadata_file.stat().st_size)} bytes\n")
                    f.write(f"  -  格式|Format: {self.config.metadata_format.upper()}\n")
                else:
                    f.write("  -  状态|Status: 元数据下载失败|Metadata download failed\n")
                f.write("\n")

                #  下载信息|Download information
                f.write(" 下载信息|Download Information:\n")
                f.write(f"  -  发现的FASTQ文件数量|FASTQ Files Found: {format_number(download_links_count)}\n")

                if script_file and script_file.exists():
                    f.write(f"  -  下载脚本|Download Script: {script_file.name}\n")
                    f.write(f"  -  脚本大小|Script Size: {format_number(script_file.stat().st_size)} bytes\n")

                    if self.config.method == "save":
                        f.write("\n")
                        f.write(" 下一步操作|Next Steps:\n")
                        f.write(f"   执行以下命令开始下载|Run the following command to start download:\n")
                        f.write(f"  bash {script_file.name}\n")
                else:
                    if self.config.method == "save":
                        f.write("  -  状态|Status: 下载脚本生成失败|Download script generation failed\n")
                    else:
                        f.write("  -  状态|Status: 直接下载已执行|Direct download executed\n")

                f.write("\n")

                #  文件列表(仅该编号相关文件)|File list (only files related to this accession)
                f.write(" 输出文件|Output Files:\n")
                for file_path in sorted(self.config.output_path.iterdir()):
                    if file_path.is_file() and (accession in file_path.name or file_path.name.endswith('.log')):
                        f.write(f"  -  {file_path.name}\n")

                f.write("\n")
                f.write("="*60 + "\n")

            self.logger.info(f" 汇总报告已生成|Summary report generated: {summary_file}")

        except Exception as e:
            self.logger.error(f" 生成汇总报告失败|Failed to generate summary report: {str(e)}")

    def generate_batch_overview(self, accessions):
        """ 批量场景总览: 每个编号的产物清单与状态|Batch overview: artifact list and status per accession"""
        overview_file = self.config.output_path / 'batch_overview.txt'
        if len(accessions) <= 1:
            # 单编号不生成总览, 保持输出精简|No overview for single accession to keep output lean
            return

        try:
            with overview_file.open('w', encoding='utf-8') as f:
                f.write("="*60 + "\n")
                f.write(" ENA批量下载总览|ENA Batch Download Overview\n")
                f.write("="*60 + "\n\n")
                f.write(f" 编号总数|Total accessions: {len(accessions)}\n")
                f.write(f" 生成时间|Generated Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                ok = 0
                for accession in accessions:
                    meta = self.config.output_path / f"{accession}.meta.{self.config.metadata_format}"
                    status = "成功|ok"
                    if not _has_data_rows_file(meta):
                        status = "无元数据|no metadata"
                    else:
                        ok += 1
                    f.write(f" [{status}] {accession}\n")

                f.write(f"\n 成功|Successful: {ok}/{len(accessions)}\n")
                f.write("="*60 + "\n")

            self.logger.info(f" 批量总览已生成|Batch overview generated: {overview_file}")

        except Exception as e:
            self.logger.error(f" 生成批量总览失败|Failed to generate batch overview: {str(e)}")


def _has_data_rows_file(path: Path) -> bool:
    """ 判断TSV元数据文件是否含数据行|Check whether a TSV metadata file has data rows"""
    if not path.exists():
        return False
    try:
        lines = [line for line in path.read_text(encoding='utf-8').strip().split('\n') if line.strip()]
    except IOError:
        return False
    return len(lines) > 1
