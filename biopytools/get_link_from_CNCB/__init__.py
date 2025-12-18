"""
从CNCB获取测序数据下载链接模块 | CNCB Sequencing Data Download Links Module

基于FTP协议从CNCB数据库批量获取测序数据下载链接的工具模块。
支持SRA、ERA、DRA数据库的Run ID查找和链接生成。

FTP-based tool module for batch downloading sequencing data links from CNCB database.
Supports Run ID lookup and link generation for SRA, ERA, DRA databases.
"""

__version__ = "1.0.0"
__author__ = "biopytools development team"

from .main import CNCLinkExtractor
from .config import CNCBConfig

__all__ = ['CNCLinkExtractor', 'CNCBConfig']