"""
主程序模块 | Main Program Module
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

# 注意: 在模块化模式下，提取脚本会添加以下导入语句:
# from .config import Config
# from .utils import Logger, PathValidator
# from .file_handler import FileHandler
# from .ftp_client import FTPClient
# from .link_finder import LinkFinder

class LinkFinderApp:
    """链接查找应用 | Link Finder Application"""
    
    def __init__(self, input_file: str, output_file: str = None, failed_file: str = None):
        # 初始化配置
        self.config = Config()
        self.config.validate()
        
        # 初始化日志
        self.logger = Logger("cncb_ftp_finder")
        
        # 验证输入文件
        PathValidator.validate_input_file(input_file)
        
        # 设置默认输出文件
        if not output_file:
            base_name = Path(input_file).stem
            output_file = f"{base_name}_links.txt"
        
        if not failed_file:
            base_name = Path(input_file).stem
            failed_file = f"{base_name}_failed.txt"
        
        # 确保输出目录存在
        self.output_file = PathValidator.ensure_output_dir(output_file)
        self.failed_file = PathValidator.ensure_output_dir(failed_file)
        
        # 初始化组件
        self.file_handler = FileHandler(self.logger)
        self.link_finder = LinkFinder(self.config, self.logger)
        
        # 存储输入文件路径
        self.input_file = input_file
    
    def run_search(self) -> bool:
        """运行搜索流程 | Run search process"""
        try:
            self.logger.info("🚀 --- 开始查找下载链接 (模块化版本) ---")
            
            # 读取项目数据
            projects = self.file_handler.read_and_group_by_project(self.input_file)
            if not projects:
                return False
            
            total_ids = sum(len(ids) for ids in projects.values())
            self.logger.info(f"📄 从 {self.input_file} 中读取到 {len(projects)} 个项目，共计 {total_ids} 个ID")
            
            # 连接FTP并搜索
            with FTPClient(self.config.ftp_host, self.config.ftp_timeout, self.logger) as ftp_client:
                if not ftp_client.connected:
                    self.logger.error("FTP连接失败，无法继续")
                    return False
                
                self.link_finder.ftp_client = ftp_client
                found_urls, failed_ids = self.link_finder.search_project_ids(ftp_client, projects)
            
            # 保存结果
            success = self.file_handler.write_links(str(self.output_file), found_urls)
            success &= self.file_handler.write_failed_ids(str(self.failed_file), failed_ids)
            
            # 输出统计信息
            self._print_summary(total_ids, len(found_urls), len(failed_ids))
            
            return success
            
        except Exception as e:
            self.logger.error(f"💥 发生严重错误，程序中断: {e}")
            return False
    
    def _print_summary(self, total_ids: int, found_count: int, failed_count: int):
        """打印统计摘要 | Print summary statistics"""
        self.logger.info("\n🎉 --- 查找完成 ---")
        self.logger.info(f"📊 总计处理: {total_ids} 个ID")
        self.logger.info(f"✔️ 成功找到文件的ID: {found_count} 个")
        self.logger.info(f"✖️ 失败 (未找到): {failed_count} 个ID，已保存至 {self.failed_file}")
        self.logger.info(f"💾 所有找到的下载链接已保存至 {self.output_file}")
        self.logger.info("\n👉 下一步：请使用aria2c进行高速下载。")

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🔍 从CNCB的FTP服务器上为给定的Project/Run ID列表查找下载链接 (模块化版本)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  %(prog)s input_file.txt
  
  # 指定输出文件 | Specify output files  
  %(prog)s input_file.txt -o output_links.txt -f failed_ids.txt

输入文件格式 | Input File Format:
  ProjectID\\tRunID
  P1\\tSRR123456
  P1\\tSRR123457
  P2\\tERR987654

输出文件 | Output Files:
  - 成功链接: [输出文件名]_links.txt
  - 失败ID: [输出文件名]_failed.txt

注意 | Note:
  这是完整的模块化版本，支持提取为独立模块文件。
        """
    )
    
    # 必需参数
    parser.add_argument("input_file", help="输入文件路径 | Input file path")
    
    # 可选参数
    parser.add_argument("-o", "--output", dest="output_file", 
                       help="保存成功找到的下载链接的输出文件 | Output file for successful download links")
    parser.add_argument("-f", "--failed", dest="failed_file",
                       help="保存未能找到的ID的输出文件 | Output file for failed IDs")
    
    args = parser.parse_args()
    
    # 创建并运行应用
    app = LinkFinderApp(
        input_file=args.input_file,
        output_file=args.output_file,
        failed_file=args.failed_file
    )
    
    success = app.run_search()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
