"""
BAM to FASTQ转换核心逻辑模块|BAM to FASTQ Conversion Core Logic
"""

import os
import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Tuple, List

from .utils import build_conda_command

# 用于从输出文件路径提取前缀|Pattern to strip output extension for bam2fastq prefix
_STRIP_OUTPUT_EXT = re.compile(r'\.(fq|fastq)(\.gz)?$', re.IGNORECASE)


class BAMConverter:
    """BAM文件转换器|BAM File Converter"""

    def __init__(self, config, logger):
        """
        初始化转换器|Initialize converter

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def _resolve_symlink(self, file_path: Path) -> Path:
        """
        正确解析软链接路径|Correctly resolve symlink path

        Args:
            file_path: 文件路径（可能是软链接）|File path (may be symlink)

        Returns:
            Path: 解析后的真实文件路径|Resolved real file path
        """
        if not file_path.is_symlink():
            return file_path

        try:
            # 读取软链接目标|Read symlink target
            link_target = os.readlink(file_path)
            self.logger.debug(f"检测到软链接|Detected symlink: {file_path} -> {link_target}")

            # 如果是绝对路径，直接使用|If absolute path, use directly
            if os.path.isabs(link_target):
                resolved = Path(link_target)
                if resolved.exists():
                    return resolved
                else:
                    self.logger.error(f"软链接指向的文件不存在|Symlink target does not exist: {resolved}")
                    return file_path

            # 使用realpath命令解析相对路径|Use realpath command to resolve relative path
            symlink_dir = file_path.parent

            # 尝试在不同的目录层级执行realpath|Try executing realpath from different directory levels
            test_dirs = [
                symlink_dir,                                    # 软链接所在目录|Symlink directory
                symlink_dir.parent,                             # 父目录|Parent directory (hifi/)
                symlink_dir.parent.parent,                      # 父父目录|Grandparent directory (raw/)
                symlink_dir.parent.parent.parent,               # 更上级|Higher level (01.data/)
                symlink_dir.parent.parent.parent.parent,        # 项目根目录|Project root
            ]

            for i, test_dir in enumerate(test_dirs):
                try:
                    result = subprocess.run(
                        ['realpath', link_target],
                        capture_output=True,
                        text=True,
                        cwd=test_dir,
                        check=False  # 不抛出异常|Don't raise exception
                    )

                    if result.returncode == 0:
                        resolved_path = Path(result.stdout.strip())
                        if resolved_path.exists():
                            self.logger.debug(f"使用realpath解析成功(层级{i})|Resolved using realpath (level {i}) from {test_dir}")
                            return resolved_path
                except Exception:
                    continue

            self.logger.error(f"realpath命令在所有目录层级都失败|realpath failed at all directory levels")

            # 如果realpath失败，尝试手动解析|If realpath fails, try manual resolution
            # 方法1: 基于软链接所在目录解析|Method 1: Resolve based on symlink's directory
            resolved_method1 = (symlink_dir / link_target).resolve()
            if resolved_method1.exists():
                self.logger.debug(f"解析成功(方法1)|Resolved (method 1): {resolved_method1}")
                return resolved_method1

            # 方法2: 基于软链接目录的父目录解析|Method 2: Resolve based on symlink's parent directory
            resolved_method2 = (symlink_dir.parent / link_target).resolve()
            if resolved_method2.exists():
                self.logger.debug(f"解析成功(方法2)|Resolved (method 2): {resolved_method2}")
                return resolved_method2

            # 方法3: 基于软链接目录的父父目录解析|Method 3: Resolve based on symlink's grandparent directory
            resolved_method3 = (symlink_dir.parent.parent / link_target).resolve()
            if resolved_method3.exists():
                self.logger.debug(f"解析成功(方法3)|Resolved (method 3): {resolved_method3}")
                return resolved_method3

            # 所有方法都失败|All methods failed
            self.logger.error(f"软链接解析失败|Symlink resolution failed: {file_path}")
            return file_path

        except Exception as e:
            self.logger.warning(f"无法解析软链接|Failed to resolve symlink: {file_path}, 错误|error: {e}")
            return file_path

    def convert_single_bam(self, bam_file: Path) -> Tuple[bool, str]:
        """
        转换单个BAM文件为FASTQ格式|Convert single BAM file to FASTQ format

        Args:
            bam_file: BAM文件路径|BAM file path

        Returns:
            Tuple[bool, str]: (是否成功|Success, 文件名|Filename)
        """
        base_name = bam_file.stem

        if self.config.output_is_file:
            # 文件输出模式: -o指定了输出文件路径|File output mode: -o specifies output file path
            # 从目标文件名推导bam2fastq的prefix|Derive bam2fastq prefix from target filename
            target_file = self.config.output_path
            parent_dir = self.config.output_dir_path
            prefix_stem = _STRIP_OUTPUT_EXT.sub('', self.config.output_filename)
            output_prefix = parent_dir / prefix_stem
        else:
            # 目录输出模式: 保持原有行为|Directory output mode: keep original behavior
            output_prefix = self.config.output_path / base_name
            target_file = None

        self.logger.info(f"正在处理|Processing: {bam_file.name}")

        # 获取真实文件路径，正确处理软链接|Get real file path, properly handle symlinks
        bam_file_to_use = self._resolve_symlink(bam_file)

        # 检查并生成PBI索引|Check and generate PBI index
        pbi_file = Path(str(bam_file_to_use) + '.pbi')
        if not pbi_file.exists():
            self.logger.info(f"生成PBI索引|Generating PBI index: {bam_file_to_use.name}")
            pbindex_cmd = build_conda_command('pbindex', [str(bam_file_to_use)])
            self.logger.info(f"命令|Command: {' '.join(pbindex_cmd)}")
            try:
                subprocess.run(
                    pbindex_cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                self.logger.error(f"PBI索引生成失败|PBI index generation failed for {bam_file.name}: {e.stderr}")
                return False, bam_file.name

        try:
            # 使用bam2fastq命令转换|Use bam2fastq to convert
            cmd_args = [
                '-o', str(output_prefix),
                '-j', str(self.config.threads),
                '-c', '6',  # 压缩级别|Compression level
                str(bam_file_to_use)
            ]
            cmd = build_conda_command(self.config.bam2fastq_path, cmd_args)
            self.logger.info(f"命令|Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # 文件输出模式：将bam2fastq生成的文件重命名为用户指定的文件名
            if self.config.output_is_file and target_file is not None:
                self._rename_output(prefix_stem, parent_dir, target_file)

            self.logger.info(f"完成|Completed: {bam_file.name}")
            return True, bam_file.name

        except subprocess.CalledProcessError as e:
            self.logger.error(f"转换失败|Conversion failed for {bam_file.name}: {e.stderr}")
            return False, bam_file.name

        except Exception as e:
            self.logger.error(f"处理出错|Error processing {bam_file.name}: {str(e)}")
            return False, bam_file.name

    def _rename_output(self, prefix_stem: str, parent_dir: Path, target_file: Path):
        """
        将bam2fastq生成的文件重命名为用户指定文件名|Rename bam2fastq output to user-specified filename

        bam2fastq产生的文件命名不稳定，可能为 {prefix}.fastq.gz 或 {prefix}.fq.gz
        此函数找到生成的文件并重命名|bam2fastq output naming is inconsistent,
        this function finds and renames the generated file.

        Args:
            prefix_stem: bam2fastq的prefix(stem)|bam2fastq prefix (stem only)
            parent_dir: 输出目录|Output directory
            target_file: 用户指定的目标文件路径|User-specified target file path
        """
        # 查找以prefix_stem开头的输出文件|Find output files starting with prefix_stem
        generated_files = sorted(parent_dir.glob(f"{prefix_stem}*"))
        if not generated_files:
            self.logger.warning(
                f"未找到bam2fastq输出文件|No bam2fastq output files found with prefix: {prefix_stem}"
            )
            return

        if len(generated_files) == 1:
            src = generated_files[0]
            # 如果目标文件已存在，先删除|If target exists, remove first
            if target_file.exists():
                target_file.unlink()
            src.rename(target_file)
            self.logger.info(f"输出文件|Output file: {target_file}")
        else:
            # 多个输出文件(多个read group)，合并后重命名|Multiple output files (multiple read groups), concatenate
            self.logger.info(
                f"合并{len(generated_files)}个输出文件|Merging {len(generated_files)} output files"
            )
            if target_file.exists():
                target_file.unlink()
            with open(target_file, 'wb') as dst:
                for src in generated_files:
                    with open(src, 'rb') as f:
                        dst.write(f.read())
                    src.unlink()  # 删除临时文件|Remove temp file
            self.logger.info(f"合并输出|Merged output: {target_file}")

    def convert_all_bams(self) -> dict:
        """
        转换所有BAM文件|Convert all BAM files

        Returns:
            dict: 转换结果统计|Conversion result statistics
        """
        # 根据输入类型获取BAM文件列表|Get BAM file list based on input type
        if self.config.is_single_file:
            # 单个文件|Single file
            bam_files = [self.config.input_path]
        else:
            # 目录：查找所有BAM文件|Directory: find all BAM files
            bam_files = list(self.config.input_path.glob("*.bam"))

            # 目录输入+文件输出模式冲突|Directory input + file output mode conflict
            if self.config.output_is_file and len(bam_files) > 1:
                self.logger.warning(
                    f"检测到{len(bam_files)}个BAM文件但指定了单文件输出|"
                    f"Found {len(bam_files)} BAM files but single output file specified: "
                    f"{self.config.output_filename}"
                )
                self.logger.warning(
                    "回退到目录输出模式|Falling back to directory output mode"
                )
                self.config.output_is_file = False
                self.config.output_dir_path = self.config.output_path
                self.config.output_dir_path.mkdir(parents=True, exist_ok=True)

        if not bam_files:
            self.logger.error(
                f"未找到BAM文件|No BAM files found in {self.config.input_dir}"
            )
            return {
                'total': 0,
                'success': 0,
                'failed': 0,
                'success_rate': 0.0
            }

        total_count = len(bam_files)
        self.logger.info(f"找到|Found {total_count} 个BAM文件|BAM files")
        self.logger.info(f"每个文件使用|Each file uses {self.config.threads} 个线程|threads")
        self.logger.info(f"并行处理|Parallel processing {self.config.jobs} 个文件|files")

        # 并行处理BAM文件|Process BAM files in parallel
        success_count = 0
        failed_count = 0
        failed_files = []

        with ThreadPoolExecutor(max_workers=self.config.jobs) as executor:
            futures = {
                executor.submit(self.convert_single_bam, bam_file): bam_file
                for bam_file in bam_files
            }

            for future in as_completed(futures):
                success, filename = future.result()
                if success:
                    success_count += 1
                else:
                    failed_count += 1
                    failed_files.append(filename)

        # 计算成功率|Calculate success rate
        success_rate = (success_count / total_count * 100) if total_count > 0 else 0.0

        # 记录转换结果|Log conversion results
        self.logger.info("=" * 50)
        self.logger.info(f"转换完成|Conversion completed!")
        self.logger.info(f"总数|Total: {total_count}")
        self.logger.info(f"成功|Success: {success_count}")
        self.logger.info(f"失败|Failed: {failed_count}")
        self.logger.info(f"成功率|Success rate: {success_rate:.2f}%")
        self.logger.info("=" * 50)

        if failed_files:
            self.logger.warning(f"失败的文件|Failed files: {', '.join(failed_files)}")

        return {
            'total': total_count,
            'success': success_count,
            'failed': failed_count,
            'success_rate': success_rate,
            'failed_files': failed_files
        }
