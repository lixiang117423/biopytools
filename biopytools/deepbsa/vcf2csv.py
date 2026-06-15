"""
VCF转CSV模块|VCF to CSV Module

从VCF文件中提取AD（Allele Depth）字段，转换为DeepBSA可用的CSV格式。
Extract AD (Allele Depth) field from VCF file and convert to DeepBSA-compatible CSV format.

CSV输出格式（无header）|CSV output format (no header):
    CHROM, POS, REF, ALT, AD_ref_s1, AD_alt_s1, AD_ref_s2, AD_alt_s2, ...
"""

import os
import sys
import logging
from pathlib import Path

import pandas as pd


class Record:
    """单行VCF记录|Single VCF record"""

    def __init__(self, line: str):
        info = line.split("\t")
        self.CHROM = info[0]
        self.POS = info[1]
        self.ID = info[2]
        self.REF = info[3]
        self.ALT = info[4]
        self.QUAL = info[5]
        self.FILTER = info[6]
        self.INFO = info[7]
        self.FORMAT = info[8].split(":")
        self.sample_num = len(info) - 9
        self.GT = []
        for i in range(self.sample_num):
            GT_value = info[8 + i + 1].split(":")
            GT_dict = {}
            for g in range(len(GT_value)):
                GT_dict[self.FORMAT[g]] = GT_value[g]
            self.GT.append(GT_dict)


class VCF:
    """VCF文件迭代器|VCF file iterator"""

    def __init__(self, vcf_path: str):
        self.header = []
        self.reader = open(vcf_path, 'r')
        self.line = self.reader.readline().strip()
        while self.line.startswith('#'):
            self.header.append(self.line)
            self.line = self.reader.readline().strip()
        if self.line:
            self.record = Record(self.line)
        else:
            self.record = None

    def __iter__(self):
        return self

    def __next__(self):
        self.line = self.reader.readline().strip()
        if self.line:
            self.record = Record(self.line)
            return self.record
        else:
            self.reader.close()
            raise StopIteration()

    def close(self):
        self.reader.close()


def vcf2csv(input_file: Path, output_file: Path, logger: logging.Logger) -> Path:
    """将VCF文件转换为DeepBSA格式的CSV|Convert VCF file to DeepBSA-format CSV

    Args:
        input_file: 输入VCF文件路径|Input VCF file path
        output_file: 输出CSV文件路径|Output CSV file path
        logger: 日志器|Logger

    Returns:
        Path: 输出CSV文件路径|Output CSV file path
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_csv = output_file

    logger.info(f"输入文件|Input file: {input_file}")
    logger.info(f"输出文件|Output file: {output_csv}")

    vcf = VCF(str(input_file))
    data = []
    error_count = 0
    total = 0

    for r in vcf:
        total += 1
        try:
            row = [
                str(r.CHROM),
                int(r.POS),
                r.REF,
                r.ALT,
            ]
            for sample_gt in r.GT:
                ad_values = sample_gt["AD"].split(",")
                if len(ad_values) >= 2:
                    row.append(int(ad_values[0]))
                    row.append(int(ad_values[1]))
                else:
                    raise ValueError("AD field incomplete")
            data.append(row)
        except Exception:
            error_count += 1

    vcf.close()

    df = pd.DataFrame(data)
    df.to_csv(output_csv, header=False, index=False, sep=',')

    sample_count = (df.shape[1] - 4) // 2
    logger.info(f"总行数|Total VCF records: {total}")
    logger.info(f"有效行数|Valid records: {len(df)}")
    logger.info(f"跳过行数|Skipped records: {error_count}")
    logger.info(f"样本数|Sample count: {sample_count}")
    logger.info(f"CSV形状|CSV shape: {df.shape}")
    logger.info(f"CSV保存到|CSV saved to: {output_csv}")

    return output_csv


def main():
    """主入口函数|Main entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='VCF转CSV - 为DeepBSA准备输入数据|VCF to CSV - Prepare input data for DeepBSA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i input.vcf -o ./
  %(prog)s -i input.vcf -o ./deepbsa_input/
"""
    )

    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-file',
        required=True,
        help='输入VCF文件路径|Input VCF file path'
    )
    required.add_argument(
        '-o', '--output-file',
        required=True,
        help='输出CSV文件路径|Output CSV file path'
    )

    args = parser.parse_args()

    input_path = Path(args.input_file).expanduser().resolve()
    output_path = Path(args.output_file).expanduser().resolve()

    if not input_path.exists():
        print(f"错误|Error: 输入文件不存在|Input file not found: {input_path}")
        sys.exit(1)

    # 配置简单日志|Setup simple logger
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger('deepbsa.vcf2csv')

    try:
        csv_path = vcf2csv(input_path, output_path, logger)
        print(f"\n转换完成|Conversion completed: {csv_path}")
    except Exception as e:
        logger.error(f"转换失败|Conversion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
