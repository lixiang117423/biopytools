"""
FOF文件生成器模块|FOF File Generator Module
"""

import os
from pathlib import Path
from typing import List, Tuple

from .utils import extract_sample_name, format_number


class FofGenerator:
    """FOF文件生成器|FOF File Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def _collect_files(self) -> List[Path]:
        """
        收集需要处理的文件列表|Collect files to process

        Returns:
            文件路径列表（已排序）|Sorted list of file paths
        """
        input_path = self.config.input_path_obj

        if self.config.is_single_file:
            self.logger.info(f"输入为单个文件|Input is a single file: {input_path.name}")
            return [input_path]

        # 目录模式|Directory mode
        scan_mode = "递归|recursive" if self.config.recursive else "非递归|non-recursive"
        self.logger.info(f"扫描目录|Scanning directory: {input_path} ({scan_mode})")

        if self.config.recursive:
            pattern = '**/*' if not self.config.suffixes else None
            if pattern:
                files = list(input_path.glob(pattern))
            else:
                files = []
                for suffix in self.config.suffixes:
                    files.extend(input_path.glob(f'**/*{suffix}'))
                # 去重并保持顺序|Deduplicate and maintain order
                seen = set()
                files = [f for f in files if not (f in seen or seen.add(f))]
        else:
            if self.config.suffixes:
                files = []
                for suffix in self.config.suffixes:
                    files.extend(input_path.glob(f'*{suffix}'))
                seen = set()
                files = [f for f in files if not (f in seen or seen.add(f))]
            else:
                files = list(input_path.glob('*'))

        # 过滤掉目录，只保留文件|Filter out directories, keep only files
        files = [f for f in files if f.is_file() or (f.is_symlink() and f.resolve().is_file())]
        # 排除输出文件自身，避免自引用|Exclude output file itself to avoid self-reference
        output_resolved = self.config.output_file_obj.resolve()
        files = [f for f in files if f.resolve() != output_resolved]
        # 按文件名排序|Sort by filename
        files.sort(key=lambda f: f.name)

        self.logger.info(f"找到文件数|Files found: {format_number(len(files))}")
        return files

    def _generate_entries(self, files: List[Path]) -> List[Tuple[str, str]]:
        """
        生成FOF条目列表|Generate FOF entry list

        Args:
            files: 文件路径列表|List of file paths

        Returns:
            (样品名, 绝对路径) 元组列表|List of (sample_name, absolute_path) tuples
        """
        entries = []
        for filepath in files:
            sample_name = extract_sample_name(filepath, self.config.suffixes)
            abs_path = str(filepath.resolve())
            entries.append((sample_name, abs_path))
            self.logger.debug(f"样品|Sample: {sample_name} -> {abs_path}")
        return entries

    def _write_fof(self, entries: List[Tuple[str, str]]) -> int:
        """
        写入FOF文件|Write FOF file

        Args:
            entries: (样品名, 绝对路径) 元组列表|List of (sample_name, absolute_path) tuples

        Returns:
            写入行数|Number of lines written
        """
        output_file = self.config.output_file_obj
        self.logger.info(f"写入FOF文件|Writing FOF file: {output_file}")

        with open(output_file, 'w', encoding='utf-8') as f:
            for sample_name, abs_path in entries:
                f.write(f"{sample_name}\t{abs_path}\n")

        line_count = len(entries)
        self.logger.info(f"写入行数|Lines written: {format_number(line_count)}")
        return line_count

    def generate(self) -> bool:
        """
        执行FOF生成|Execute FOF generation

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 50)
        self.logger.info("FOF文件生成|FOF File Generation")
        self.logger.info("=" * 50)

        # 收集文件|Collect files
        files = self._collect_files()

        if not files:
            self.logger.warning("未找到任何文件|No files found")
            # 仍然创建空的FOF文件|Still create an empty FOF file
            self._write_fof([])
            return True

        # 生成条目|Generate entries
        entries = self._generate_entries(files)

        # 写入FOF文件|Write FOF file
        self._write_fof(entries)

        # 检查重复样品名|Check for duplicate sample names
        sample_names = [e[0] for e in entries]
        duplicates = [n for n in sample_names if sample_names.count(n) > 1]
        if duplicates:
            unique_dups = list(set(duplicates))
            self.logger.warning(
                f"检测到重复样品名|Duplicate sample names detected: {', '.join(unique_dups)}"
            )

        self.logger.info("FOF文件生成完成|FOF file generation completed")
        return True
