#!/usr/bin/env python3
"""
测试DeepBSA绘图功能|Test DeepBSA Plotting
"""

from pathlib import Path
from plot_generator import plot_deepbsa_results

# 测试数据文件路径（使用你已有的数据）
data_file = Path("/share/org/YZWL/yzwl_lixg/project/18.小花糖芥/17.BSA重新比对到68-1基因组/05.deepbsa/merged_results/plot_data_for_R.csv")

# 输出文件路径
output_file = Path("/tmp/deepbsa_test_plot.png")

if data_file.exists():
    print(f"读取数据|Reading data: {data_file}")

    # 调试信息|Debug info
    from plot_generator import DeepBSAPlotter
    print(f"COLOR_POINT值|COLOR_POINT value: {DeepBSAPlotter.COLOR_POINT}")

    # 调用绘图函数|Call plotting function
    success = plot_deepbsa_results(
        data_file=data_file,
        output_file=output_file,
        width=12,
        height=8,
        dpi=300,
        color_threshold="#1F77B4FF",  # 蓝色|Blue
        color_smooth="#FF7F0EFF",      # 橙色|Orange
        color_point="grey80",          # 灰色|Grey
        logger=None
    )

    if success:
        print(f"✓ 绘图成功|Plotting successful!")
        print(f"  输出文件|Output file: {output_file}")
        print(f"  文件大小|File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    else:
        print("✗ 绘图失败|Plotting failed")
else:
    print(f"✗ 数据文件不存在|Data file not found: {data_file}")
    print("提示|Hint: 请先运行deepbsa生成绘图数据|Please run deepbsa first to generate plot data")
