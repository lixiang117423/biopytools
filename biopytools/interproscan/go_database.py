"""
GO术语数据库模块|GO Term Database Module
加载和管理GO术语信息（ID、名称、ontology分类）|Load and manage GO term information (ID, name, ontology)
优先使用内置数据库，可选加载外部JSON文件|Use built-in database by default, optionally load external JSON
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional
from threading import Lock

# 导入内置GO数据库|Import built-in GO database
from .go_data import GO_DATABASE, GO_STATS


class GODatabase:
    """GO术语数据库|GO Term Database"""

    _instance = None
    _lock = Lock()

    def __new__(cls, *args, **kwargs):
        """单例模式|Singleton pattern"""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, go_json_path: Optional[str] = None, logger: Optional[logging.Logger] = None,
                 use_builtin: bool = True):
        """
        初始化GO数据库|Initialize GO database

        Args:
            go_json_path: 外部GO basic JSON文件路径（可选）|External GO basic JSON file path (optional)
            logger: 日志器|Logger
            use_builtin: 是否使用内置数据库（默认True）|Whether to use built-in database (default True)
        """
        # 如果已经初始化过，跳过|Skip if already initialized
        if hasattr(self, '_initialized') and self._initialized:
            return

        self.logger = logger or logging.getLogger(__name__)
        self.go_map: Dict[str, Dict[str, str]] = {}
        self._source = "unknown"
        self._initialized = False

        # 优先级：外部JSON > 内置数据库|Priority: External JSON > Built-in database
        if go_json_path and Path(go_json_path).exists():
            self._load_external_json(go_json_path)
        elif use_builtin:
            self._load_builtin_database()
        else:
            self.logger.warning("未指定GO数据库来源|No GO database source specified")

    def _load_builtin_database(self):
        """加载内置GO数据库|Load built-in GO database"""
        self.logger.info("使用内置GO数据库|Using built-in GO database")

        # 转换内置数据格式|Convert built-in data format
        for go_id, (name, ontology) in GO_DATABASE.items():
            # 映射ontology缩写为全称|Map ontology abbreviations to full names
            ontology_map = {
                'BP': 'Biological Process',
                'MF': 'Molecular Function',
                'CC': 'Cellular Component',
                'Other': 'Other'
            }
            ontology_full = ontology_map.get(ontology, ontology)

            self.go_map[go_id] = {
                'name': name,
                'ontology': ontology_full
            }

        self._source = "built-in"
        self._initialized = True

        stats = self.get_stats()
        self.logger.info(f"内置GO数据库加载完成|Built-in GO database loaded: {stats['total']} GO terms "
                        f"(BP: {stats['biological_process']}, MF: {stats['molecular_function']}, "
                        f"CC: {stats['cellular_component']})")

    def _load_external_json(self, go_json_path: str):
        """
        加载外部GO数据库|Load GO database from external JSON file

        Args:
            go_json_path: GO basic JSON文件路径|GO basic JSON file path
        """
        self.logger.info(f"加载外部GO数据库|Loading external GO database: {go_json_path}")

        try:
            with open(go_json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            nodes = data.get('graphs', [{}])[0].get('nodes', [])

            ontology_map = {
                'biological_process': 'Biological Process',
                'molecular_function': 'Molecular Function',
                'cellular_component': 'Cellular Component'
            }

            for node in nodes:
                # 提取GO ID（从URI中提取）|Extract GO ID from URI
                go_uri = node.get('id', '')
                if not go_uri.startswith('http://purl.obolibrary.org/obo/GO_'):
                    continue

                go_id = 'GO:' + go_uri.split('GO_')[-1]

                # 提取名称|Extract name
                name = node.get('lbl', '')

                # 提取ontology（namespace）|Extract ontology (namespace)
                ontology = 'Unknown'
                if 'meta' in node and 'basicPropertyValues' in node['meta']:
                    for prop in node['meta']['basicPropertyValues']:
                        if 'hasOBONamespace' in prop.get('pred', ''):
                            ontology = prop.get('val', 'Unknown')
                            break

                ontology_full = ontology_map.get(ontology, ontology)

                self.go_map[go_id] = {
                    'name': name,
                    'ontology': ontology_full
                }

            self._source = "external"
            self._initialized = True

            stats = self.get_stats()
            self.logger.info(f"外部GO数据库加载完成|External GO database loaded: {stats['total']} GO terms "
                            f"(BP: {stats['biological_process']}, MF: {stats['molecular_function']}, "
                            f"CC: {stats['cellular_component']})")

        except Exception as e:
            self.logger.error(f"加载外部GO数据库失败，尝试使用内置数据库|Failed to load external GO database, trying built-in: {str(e)}")
            self._load_builtin_database()

    def get_go_info(self, go_id: str) -> Dict[str, str]:
        """
        获取GO术语信息|Get GO term information

        Args:
            go_id: GO ID (e.g., "GO:0003700")

        Returns:
            包含name和ontology的字典|Dictionary containing name and ontology
        """
        # 清理GO ID格式|Clean GO ID format
        go_id_clean = go_id.strip()
        if not go_id_clean.startswith('GO:'):
            # 如果是从TSV中解析的格式如 "GO:0003700(PANTHER)"
            go_id_clean = go_id_clean.split('(')[0].split('|')[0].strip()

        return self.go_map.get(go_id_clean, {'name': '', 'ontology': ''})

    def get_go_name(self, go_id: str) -> str:
        """获取GO术语名称|Get GO term name"""
        return self.get_go_info(go_id).get('name', '')

    def get_go_ontology(self, go_id: str) -> str:
        """获取GO术语ontology|Get GO term ontology"""
        return self.get_go_info(go_id).get('ontology', '')

    def is_loaded(self) -> bool:
        """检查GO数据库是否已加载|Check if GO database is loaded"""
        return self._initialized and len(self.go_map) > 0

    def get_source(self) -> str:
        """获取数据库来源|Get database source (built-in or external)"""
        return self._source

    def get_stats(self) -> Dict[str, int]:
        """获取GO数据库统计信息|Get GO database statistics"""
        stats = {
            'total': len(self.go_map),
            'biological_process': 0,
            'molecular_function': 0,
            'cellular_component': 0,
            'other': 0
        }

        for go_info in self.go_map.values():
            ontology = go_info.get('ontology', 'Unknown')
            if 'Biological Process' in ontology:
                stats['biological_process'] += 1
            elif 'Molecular Function' in ontology:
                stats['molecular_function'] += 1
            elif 'Cellular Component' in ontology:
                stats['cellular_component'] += 1
            else:
                stats['other'] += 1

        return stats
