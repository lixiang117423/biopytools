"""确保biopytools可导入|Ensure biopytools importable from repo root"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def pytest_configure(config):
    """注册e2e marker|Register e2e marker(登录节点默认 -m "not e2e" 跳过)"""
    config.addinivalue_line(
        'markers', 'e2e: 计算节点e2e|compute-node end-to-end')
