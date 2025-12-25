"""
CNCB FTP链接查找工具包 | CNCB FTP Link Finder Toolkit
功能: 从CNCB FTP服务器查找生物数据下载链接 | 
Features: Find biological data download links from CNCB FTP server
作者 | Author: Claude
版本 | Version: v2.0 - 模块化版本 | Modular version
日期 | Date: 2025-11-04

使用示例 | Usage Examples:
    from biopytools.find_cncb_link import LinkFinder, Config
    
    # 创建链接查找器 | Create link finder
    finder = LinkFinder(
        input_file="project_run_ids.txt",
        output_file="links.txt",
        failed_file="failed.txt"
    )
    
    # 执行查找 | Execute search
    finder.run_search()
"""

__version__ = "2.0.0"
__author__ = "Claude"

# 注意: 这个文件在单体模式下不会执行导入
# 在模块化模式下，提取脚本会添加正确的导入语句
# from .main import LinkFinder
# from .config import Config

__all__ = ['LinkFinder', 'Config']
