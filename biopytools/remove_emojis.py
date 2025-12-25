#!/usr/bin/env python3
"""
移除代码和文档中的emoji字符
根据开发规范要求，移除所有用户不友好的emoji符号
"""

import os
import re
import sys
from pathlib import Path
from typing import List, Tuple

# 定义需要移除的emoji字符列表
EMOJI_PATTERN = re.compile(
    '['
    '\U0001F300-\U0001F5FF'  # Symbols & pictographs
    '\U0001F600-\U0001F64F'  # Emoticons
    '\U0001F680-\U0001F6FF'  # Transport & map symbols
    '\U0001F700-\U0001F77F'  # Alchemical symbols
    '\U0001F780-\U0001F7FF'  # Geometric shapes
    '\U0001F800-\U0001F8FF'  # Supplemental arrows-C
    '\U0001F900-\U0001F9FF'  # Supplemental symbols and pictographs
    '\U0001FA00-\U0001FA6F'  # Chess symbols
    '\U0001FA70-\U0001FAFF'  # Symbols and pictographs Extended-A
    '\U00002702-\U000027B0'  # Dingbats
    '\U000024C2-\U0001F251'  # Enclosed characters
    ']+'
)

# 常见emoji及其替换文本的映射（可选）
EMOJI_REPLACEMENTS = {
    '🧬': '',  # DNA相关
    '📊': '',  # 数据/图表
    '📄': '',  # 文件
    '🔢': '',  # 数字
    '📂': '',  # 目录
    '📝': '',  # 文档
    '🔄': '',  # 循环/处理
    '🎯': '',  # 目标
    '🔧': '',  # 工具
    '💾': '',  # 存储
    '🚀': '',  # 快速
    '⏭️': '',  # 跳过
    '📈': '',  # 图表
    '🥤': '',  # Juicebox
    '🛠️': '',  # 工具
    '🔗': '',  # 链接
    '🔍': '',  # 搜索
    '🧪': '',  # 测试
    '💥': '',  # 错误
    '⚠️': '',  # 警告
    '✅': '',  # 成功
    '❌': '',  # 失败
    '🗄️': '',  # 数据库
    '⚡': '',  # 快速
    '📖': '',  # 文档
    '✨': '',  # 特性
    '🔬': '',  # 分析
    '🛡️': '',  # 保护
}


def remove_emojis_from_file(file_path: Path, backup: bool = True) -> Tuple[bool, int]:
    """
    从文件中移除emoji

    Args:
        file_path: 文件路径
        backup: 是否创建备份

    Returns:
        (是否修改, 移除的emoji数量)
    """
    try:
        # 读取原始内容
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # 备份原始文件
        if backup:
            backup_path = file_path.with_suffix(f'{file_path.suffix}.emoji_backup')
            with open(backup_path, 'w', encoding='utf-8') as f:
                f.write(content)

        # 统计原始emoji数量
        original_emojis = EMOJI_PATTERN.findall(content)
        emoji_count = len(original_emojis)

        if emoji_count == 0:
            return False, 0

        # 移除emoji - 使用更精确的替换
        cleaned_content = content
        for emoji, replacement in EMOJI_REPLACEMENTS.items():
            cleaned_content = cleaned_content.replace(emoji, replacement)

        # 额外清理：移除emoji后可能留下的多余空格
        cleaned_content = re.sub(r'\s+', ' ', cleaned_content)
        cleaned_content = re.sub(r'\n\s*\n', '\n\n', cleaned_content)

        # 写入清理后的内容
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(cleaned_content)

        return True, emoji_count

    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {e}", file=sys.stderr)
        return False, 0


def find_emoji_files(directory: Path, extensions: List[str]) -> List[Path]:
    """
    查找包含emoji的文件

    Args:
        directory: 搜索目录
        extensions: 文件扩展名列表

    Returns:
        包含emoji的文件路径列表
    """
    emoji_files = []

    for ext in extensions:
        for file_path in directory.rglob(f'*{ext}'):
            if file_path.is_file():
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read(1024)  # 只读取前1024字符进行快速检测
                        if EMOJI_PATTERN.search(content):
                            emoji_files.append(file_path)
                except Exception:
                    # 忽略无法读取的文件
                    continue

    return emoji_files


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(
        description='移除代码和文档中的emoji字符',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  %(prog)s --cli-only          # 只处理CLI包装器文件
  %(prog)s --docs-only         # 只处理文档文件
  %(prog)s --source-only       # 只处理源码文件（除CLI外的所有Python文件）
  %(prog)s --all               # 处理所有文件（CLI、文档、源码）
  %(prog)s --dry-run           # 预览模式，不实际修改文件
  %(prog)s --no-backup         # 不创建备份文件
        '''
    )

    parser.add_argument('--cli-only', action='store_true',
                       help='只处理CLI包装器文件')
    parser.add_argument('--docs-only', action='store_true',
                       help='只处理文档文件')
    parser.add_argument('--source-only', action='store_true',
                       help='只处理源码文件（除CLI外的所有Python文件）')
    parser.add_argument('--all', action='store_true',
                       help='处理所有文件（CLI、文档、源码）')
    parser.add_argument('--dry-run', action='store_true',
                       help='预览模式，只显示将要修改的文件')
    parser.add_argument('--no-backup', action='store_true',
                       help='不创建备份文件')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='详细输出')

    args = parser.parse_args()

    # 设置基础路径
    base_dir = Path(__file__).parent
    cli_dir = base_dir / 'cli' / 'commands'
    docs_dir = base_dir.parent / 'docs'
    source_dir = base_dir

    files_to_process = []

    # 收集要处理的文件
    if args.cli_only:
        # 只处理CLI文件
        if cli_dir.exists():
            cli_files = find_emoji_files(cli_dir, ['.py'])
            files_to_process.extend(cli_files)
            if args.verbose:
                print(f"找到 {len(cli_files)} 个包含emoji的CLI文件")
    elif args.docs_only:
        # 只处理文档文件
        if docs_dir.exists():
            docs_files = find_emoji_files(docs_dir, ['.md'])
            files_to_process.extend(docs_files)
            if args.verbose:
                print(f"找到 {len(docs_files)} 个包含emoji的文档文件")
    elif args.source_only:
        # 只处理源码文件（排除CLI目录）
        if source_dir.exists():
            # 递归查找所有Python文件，但排除CLI目录
            for py_file in source_dir.rglob('*.py'):
                if py_file.is_file() and cli_dir not in py_file.parents:
                    try:
                        with open(py_file, 'r', encoding='utf-8') as f:
                            content = f.read(1024)
                            if EMOJI_PATTERN.search(content):
                                files_to_process.append(py_file)
                    except Exception:
                        continue
            if args.verbose:
                print(f"找到 {len([f for f in files_to_process if f.suffix == '.py'])} 个包含emoji的源码文件")
    elif args.all:
        # 处理所有文件
        # CLI文件
        if cli_dir.exists():
            cli_files = find_emoji_files(cli_dir, ['.py'])
            files_to_process.extend(cli_files)
            if args.verbose:
                print(f"找到 {len(cli_files)} 个包含emoji的CLI文件")

        # 文档文件
        if docs_dir.exists():
            docs_files = find_emoji_files(docs_dir, ['.md'])
            files_to_process.extend(docs_files)
            if args.verbose:
                print(f"找到 {len(docs_files)} 个包含emoji的文档文件")

        # 其他源码文件（排除CLI）
        for py_file in source_dir.rglob('*.py'):
            if py_file.is_file() and cli_dir not in py_file.parents:
                try:
                    with open(py_file, 'r', encoding='utf-8') as f:
                        content = f.read(1024)
                        if EMOJI_PATTERN.search(content):
                            files_to_process.append(py_file)
                except Exception:
                    continue
        if args.verbose:
            source_count = len([f for f in files_to_process if f.suffix == '.py']) - len([f for f in files_to_process if cli_dir in f.parents])
            print(f"找到 {source_count} 个包含emoji的其他源码文件")
    else:
        # 默认行为：处理CLI和文档
        if cli_dir.exists():
            cli_files = find_emoji_files(cli_dir, ['.py'])
            files_to_process.extend(cli_files)
            if args.verbose:
                print(f"找到 {len(cli_files)} 个包含emoji的CLI文件")

        if docs_dir.exists():
            docs_files = find_emoji_files(docs_dir, ['.md'])
            files_to_process.extend(docs_files)
            if args.verbose:
                print(f"找到 {len(docs_files)} 个包含emoji的文档文件")

    if not files_to_process:
        print("未找到包含emoji的文件")
        return 0

    print(f"\n找到 {len(files_to_process)} 个包含emoji的文件:")
    for file_path in files_to_process:
        relative_path = file_path.relative_to(base_dir.parent)
        print(f"  - {relative_path}")

    if args.dry_run:
        print(f"\n预览模式：将处理以上 {len(files_to_process)} 个文件")
        return 0

    # 确认处理
    response = input(f"\n确认处理这 {len(files_to_process)} 个文件吗? (y/N): ")
    if response.lower() not in ['y', 'yes']:
        print("操作已取消")
        return 0

    # 处理文件
    total_modified = 0
    total_emojis_removed = 0

    for file_path in files_to_process:
        relative_path = file_path.relative_to(base_dir.parent)
        modified, emoji_count = remove_emojis_from_file(file_path, backup=not args.no_backup)

        if modified:
            total_modified += 1
            total_emojis_removed += emoji_count
            print(f"✓ {relative_path} (移除了 {emoji_count} 个emoji)")
        else:
            if args.verbose:
                print(f"- {relative_path} (无需修改)")

    # 输出总结
    print(f"\n处理完成:")
    print(f"  修改文件数: {total_modified}")
    print(f"  移除emoji数: {total_emojis_removed}")

    if not args.no_backup:
        print(f"  备份文件: 已创建 .emoji_backup 文件")

    return 0


if __name__ == '__main__':
    sys.exit(main())