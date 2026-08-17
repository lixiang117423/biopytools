"""
从CNCB获取测序数据下载链接模块|CNCB Sequencing Data Download Links Module

基于FTP协议从CNCB数据库批量获取测序数据下载链接的工具模块。
支持SRA、ERA、DRA数据库的Run ID查找和链接生成。
CNCB未找到的Run ID自动回退到ENA Portal API查询下载链接。

FTP-based tool module for batch downloading sequencing data links from CNCB database.
Supports Run ID lookup and link generation for SRA, ERA, DRA databases.
Run IDs missed by CNCB automatically fall back to the ENA Portal API.
"""

__version__ = "1.1.0"
__author__ = "biopytools development team"

from .main import CNCLinkExtractor
from .config import CNCBConfig
from .ena_searcher import ENALinkSearcher

__all__ = ['CNCLinkExtractor', 'CNCBConfig', 'ENALinkSearcher']
