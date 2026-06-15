"""微观共线性可视化配置|Microsynteny Visualization Configuration"""

from dataclasses import dataclass
from typing import Optional

from ..config import JcviBaseConfig


@dataclass
class MicroConfig(JcviBaseConfig):
    """微观共线性可视化配置|Microsynteny Visualization Configuration

    继承JcviBaseConfig获得mcscan管道字段, 新增区间和基因过滤字段
    Inherits JcviBaseConfig for mcscan pipeline fields, adds region and gene filtering
    """

    region_a: str = ""
    region_b: str = ""
    genes_a: Optional[str] = None       # 基因列表文件路径|gene list file path
    genes_b: Optional[str] = None       # 基因列表文件路径|gene list file path
    glyph_style: str = "arrow"           # 基因 glyph 样式: arrow, box等|gene glyph style
    shadestyle: str = ""                # 阴影样式: line, gradient等|shade style
    liftover: bool = False              # 是否启用liftover推断|enable liftover inference
    extend_blocks: int = 5
    iter_count: int = 1

    def validate(self):
        super().validate()
        errors = []
        if not self.pairs or not any(
            isinstance(self.pairs, list) and self.pairs
        ):
            errors.append(
                "需要 --pairs A,B 参数|Parameter --pairs A,B is required"
            )
        if not self.region_a:
            errors.append(
                "需要 --region-a chr:start-end 参数|Parameter --region-a is required"
            )
        if not self.region_b:
            errors.append(
                "需要 --region-b chr:start-end 参数|Parameter --region-b is required"
            )
        for label, region in [("region_a", self.region_a), ("region_b", self.region_b)]:
            if region and ":" not in region:
                errors.append(
                    f"{label} 格式错误, 应为 chr:start-end|"
                    f"{label} format invalid, expected chr:start-end: {region}"
                )
        if errors:
            raise ValueError("\n".join(errors))
        return True
