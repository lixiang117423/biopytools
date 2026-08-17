"""
从CNCB获取测序数据下载链接模块|CNCB Sequencing Data Download Links Module

基于FTP协议从CNCB数据库批量获取测序数据下载链接的工具模块。
支持SRA、ERA、DRA数据库的Run ID查找和链接生成。
CNCB未找到的Run ID自动回退到ENA Portal API查询下载链接。
GSA原生Run ID(CRR前缀)通过HTTPS通道查询(项目列需为CRA编号)。

FTP-based tool module for batch downloading sequencing data links from CNCB database.
Supports Run ID lookup and link generation for SRA, ERA, DRA databases.
Run IDs missed by CNCB automatically fall back to the ENA Portal API.
GSA-native run IDs (CRR prefix) are resolved via HTTPS (project column must be the CRA accession).
"""

__version__ = "1.2.1"
__author__ = "biopytools development team"

from .main import CNCLinkExtractor
from .config import CNCBConfig
from .ena_searcher import ENALinkSearcher
from .gsa_searcher import GSAHTTPSearcher

__all__ = ['CNCLinkExtractor', 'CNCBConfig', 'ENALinkSearcher', 'GSAHTTPSearcher']
