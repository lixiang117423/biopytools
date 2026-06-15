"""宏观共线性可视化配置|Macrosynteny Visualization Configuration"""

from dataclasses import dataclass

from ..config import JcviBaseConfig


@dataclass
class MacroConfig(JcviBaseConfig):
    """宏观共线性可视化配置|Macrosynteny Visualization Configuration

    继承JcviBaseConfig获得mcscan管道字段, 新增karyotype可视化字段
    Inherits JcviBaseConfig for mcscan pipeline fields, adds karyotype visualization
    """

    species: str = ""

    # 筛选参数|Filter parameters
    minspan: int = 30             # screen最小跨度|minimal span for screen
    min_chr_genes: int = 20       # 染色体最小基因数|min genes per chromosome

    # 可视化参数|Visualization parameters
    figsize: str = ""             # 画布大小, 如"14x10"|figure size, e.g. "14x10"
    shadestyle: str = "line"      # 阴影样式: line, gradient等|shade style
    chrstyle: str = "rect"        # 染色体样式: rect, roundrect等|chromosome style

    # 断点续传|Checkpoint resume
    replot: bool = False          # 仅重新绘图|re-plot only

    @property
    def species_list(self) -> list:
        """解析物种列表为有序列表|Parse species string to ordered list"""
        return [s.strip() for s in self.species.split(',') if s.strip()]

    @property
    def is_multi_species(self) -> bool:
        """是否为多物种模式(>=3)|Whether multi-species mode (>=3)"""
        return len(self.species_list) >= 3

    def validate(self):
        super().validate()
        errors = []
        if not self.species or len(self.species_list) < 2:
            errors.append(
                f"至少需要2个物种, 当前: {self.species}|"
                f"At least 2 species required, got: {self.species}"
            )
        if errors:
            raise ValueError("\n".join(errors))
        return True
