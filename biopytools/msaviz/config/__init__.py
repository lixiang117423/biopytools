from __future__ import annotations

import csv
from enum import IntEnum, auto
from pathlib import Path

###########################################################
# Color Schemes Config
###########################################################


def get_color_schemes() -> dict[str, dict[str, str]]:
    """Get color schemes

    Returns
    -------
    name2color_scheme : dict[str, dict[str, str]]
        Color schemes dict
    """
    COLOR_SCHEMES_FILE = Path(__file__).parent / "color_schemes.tsv"
    name2color_scheme = {}
    with open(COLOR_SCHEMES_FILE) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        letters = header[1:]
        for row in reader:
            name, colors = row[0], row[1:]
            color_scheme = {}
            for letter, color in zip(letters, colors):
                color_scheme[letter] = color
            name2color_scheme[name] = color_scheme
    return name2color_scheme


COLOR_SCHEMES = get_color_schemes()

###########################################################
# Plot Config
###########################################################


class AxesType(IntEnum):
    """Plot axes type enum|绘图坐标轴类型枚举"""

    MSA = auto()
    CONSENSUS = auto()
    SPACE = auto()
    WRAP_SPACE = auto()
