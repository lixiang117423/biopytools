"""确保biopytools可导入|Ensure biopytools importable from repo root"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
