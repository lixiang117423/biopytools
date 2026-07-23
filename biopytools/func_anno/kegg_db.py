"""
KEGG pathway 静态映射查表|KEGG pathway static map lookup.

加载内置 kegg_data.py(由 scripts/build_kegg_pathway_map.py 生成, 含 name + category_A + category_B),
提供 ko_id -> (name, category_A, category_B) 查询. 纯离线, 计算节点无网可用.
|Load built-in kegg_data.py, provide ko_id -> (name, category_A, category_B) lookup. Offline.

category_A(7 大类, 用于过滤非物种通路)|category_A (7 top-level, for filtering):
    Metabolism / Genetic Information Processing / Environmental Information Processing /
    Cellular Processes / Organismal Systems / Human Diseases / Drug Development
"""

import logging
from pathlib import Path
from threading import Lock
from typing import Dict, Optional

from .kegg_data import KEGG_PATHWAY


class KEGGDatabase:
    """KEGG pathway 查表(单例)|KEGG pathway lookup (singleton)."""

    _instance = None
    _lock = Lock()

    def __new__(cls, *args, **kwargs):
        """单例模式|Singleton pattern."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, kegg_map_file: Optional[str] = None,
                 logger: Optional[logging.Logger] = None):
        """
        初始化|Initialize.

        Args:
            kegg_map_file: 可选外部映射 TSV(ko_id\\tname[\\tcat_a[\\tcat_b]]), 覆盖/补全内置.
            logger: 日志器|Logger.
        """
        if hasattr(self, "_initialized") and self._initialized:
            if kegg_map_file:
                self._load_external(kegg_map_file)
            return

        self.logger = logger or logging.getLogger(__name__)
        self.pathway_map: Dict[str, Dict[str, str]] = {}

        for kid, val in KEGG_PATHWAY.items():
            name, cat_a, cat_b = self._normalize(val)
            self.pathway_map[kid] = {"name": name, "category_a": cat_a, "category_b": cat_b}

        total = len(self.pathway_map)
        with_cat = sum(1 for v in self.pathway_map.values() if v["category_a"])
        self.logger.info(
            f"内置 KEGG 映射加载|Loaded KEGG map: {total} 条|entries "
            f"(category 命中|with category: {with_cat})"
        )

        if kegg_map_file:
            self._load_external(kegg_map_file)

        self._initialized = True

    @staticmethod
    def _normalize(val) -> tuple:
        """归一化 kegg_data 的值 → (name, category_a, category_b)|Normalize."""
        if isinstance(val, (tuple, list)):
            name = val[0] if val else ""
            cat_a = val[1] if len(val) > 1 else ""
            cat_b = val[2] if len(val) > 2 else ""
            return str(name), str(cat_a), str(cat_b)
        # 旧格式: 纯字符串 name|legacy bare string
        return str(val), "", ""

    def _load_external(self, path: str):
        """加载外部 TSV 覆盖/补全|Load external TSV. Format: ko_id\\tname[\\tcat_a[\\tcat_b]]."""
        p = Path(path)
        if not p.exists():
            self.logger.warning(f"外部 KEGG 映射不存在|External KEGG map not found: {path}")
            return
        n = 0
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                kid = parts[0].strip()
                name = parts[1].strip()
                cat_a = parts[2].strip() if len(parts) > 2 else ""
                cat_b = parts[3].strip() if len(parts) > 3 else ""
                existing = self.pathway_map.get(kid, {"name": "", "category_a": "", "category_b": ""})
                if name:
                    existing["name"] = name
                if cat_a:
                    existing["category_a"] = cat_a
                if cat_b:
                    existing["category_b"] = cat_b
                self.pathway_map[kid] = existing
                n += 1
        self.logger.info(f"外部 KEGG 映射合并|Merged external KEGG map: {n} 条|entries ({path})")

    def get_pathway_info(self, ko_id: str) -> Dict[str, str]:
        """
        查 pathway 信息|Look up pathway info.

        Returns:
            {"name":..., "category_a":..., "category_b":...}; 未命中返回空串.
        """
        kid = (ko_id or "").strip()
        if not kid or kid == "-":
            return {"name": "", "category_a": "", "category_b": ""}
        if not kid.startswith("ko") and kid[:5].isdigit():
            kid = "ko" + kid
        return self.pathway_map.get(kid, {"name": "", "category_a": "", "category_b": ""})

    def get_name(self, ko_id: str) -> str:
        """获取 pathway 描述|Get pathway description."""
        return self.get_pathway_info(ko_id).get("name", "")

    def get_category_a(self, ko_id: str) -> str:
        """获取 A 级分类(7 大类)|Get A-level category."""
        return self.get_pathway_info(ko_id).get("category_a", "")

    def get_category_b(self, ko_id: str) -> str:
        """获取 B 级分类(子类)|Get B-level category."""
        return self.get_pathway_info(ko_id).get("category_b", "")

    def get_category(self, ko_id: str) -> str:
        """兼容旧接口: 返回 A 级分类|Compat: return A-level category."""
        return self.get_category_a(ko_id)
