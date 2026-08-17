"""
BRAKER4Phyto疫霉菌基因组注释主程序模块|BRAKER4Phyto Genome Annotation Main Module

薄包装复用 annorefine 的端到端流程(BRAKER 注释 + 同源查漏补缺 → 整合 GFF3),
唯一差异: 重复序列屏蔽默认关闭(疫霉效应子多位于 repeat 区)
|Thin wrapper reusing annorefine's end-to-end flow (BRAKER + homology
gap-filling → integrated GFF3); only difference: repeat masking off by default
"""

from ..annorefine.main import (
    parse_arguments as _annorefine_parse_arguments,
    run_end_to_end,
)


def parse_arguments():
    """解析命令行参数(疫霉默认: 不屏蔽重复)|Parse command line arguments with phyto defaults"""
    return _annorefine_parse_arguments(skip_repeat_default=True,
                                       prog_name="braker4phyto")


def main():
    """主入口|Main entry: 与 annorefine 相同,仅默认跳过 repeat 屏蔽"""
    args = parse_arguments()
    run_end_to_end(args, prog_name="braker4phyto")


if __name__ == '__main__':
    main()
