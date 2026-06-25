#!/usr/bin/env python3
"""
测试DeepBSA绘图功能|Test DeepBSA Plotting

示例|Examples: python test_plot.py -i plot_data_for_R.csv -o /tmp/deepbsa_test_plot.png
"""

import argparse
from pathlib import Path
from plot_generator import plot_deepbsa_results


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(description="测试DeepBSA绘图功能|Test DeepBSA plotting")
    parser.add_argument('-i', '--input', required=True,
                        help='绘图数据CSV文件|Plot data CSV file')
    parser.add_argument('-o', '--output', default='/tmp/deepbsa_test_plot.png',
                        help='输出图片路径|Output image path')
    args = parser.parse_args()

    # 数据文件与输出文件|Data file and output file
    data_file = Path(args.input)
    output_file = Path(args.output)

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
            print(f"绘图成功|Plotting successful!")
            print(f"输出文件|Output file: {output_file}")
            print(f"文件大小|File size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
        else:
            print("绘图失败|Plotting failed")
    else:
        print(f"数据文件不存在|Data file not found: {data_file}")
        print("提示|Hint: 请先运行deepbsa生成绘图数据|Please run deepbsa first to generate plot data")


if __name__ == '__main__':
    main()
