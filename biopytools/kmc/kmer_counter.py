"""
KMC k-mer统计模块|KMC K-mer Counting Module
"""

import os
import shutil
from pathlib import Path
from typing import List, Optional

from .config import KMCConfig
from .utils import KMCLogger, CommandRunner
from .sample_finder import KMCSampleFinder


class KMCCounter:
    """KMC k-mer统计器|KMC K-mer Counter"""

    def __init__(self, config: KMCConfig, logger: Optional[KMCLogger] = None):
        """初始化KMC统计器|Initialize KMC counter"""
        self.config = config
        self.config.validate()

        if logger is None:
            self.logger_manager = KMCLogger(self.config.output_path)
            self.logger = self.logger_manager.get_logger()
        else:
            self.logger_manager = logger
            self.logger = logger.get_logger()

        # CommandRunner的工作目录设置为kmc_databases，这样KMC生成的文件会直接在正确的位置|Set CommandRunner working dir to kmc_databases so KMC outputs are in correct place
        self.cmd_runner = CommandRunner(self.logger, self.config.kmc_db_path)
        self.sample_finder = KMCSampleFinder(self.config, self.logger)

    def count_sample(self, input_file: str, sample_name: str, tmp_dir: Optional[str] = None) -> bool:
        """统计单个样本的k-mer|Count k-mers for a single sample"""
        self.logger.info(f"开始统计样本|Starting to count sample: {sample_name}")

        if tmp_dir is None:
            tmp_dir = str(self.config.tmp_path / sample_name)

        # 创建临时目录|Create temporary directory
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)

        # KMC输出文件名（只是文件名，不是路径）|KMC output file name (filename only, not path)
        kmc_output_name = sample_name

        # 最终输出目录|Final output directory
        output_dir = self.config.kmc_db_path
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建输入文件列表（KMC需要用@file.list格式传递多个文件）|Create input file list
        input_list_file = Path(tmp_dir) / "input_files.txt"
        input_files = input_file.split()  # 输入可能是多个文件用空格分隔

        # 将输入文件转换为绝对路径|Convert input files to absolute paths
        with open(input_list_file, 'w') as f:
            for file_path in input_files:
                abs_path = str(Path(file_path).resolve())
                f.write(f"{abs_path}\n")

        self.logger.debug(f"创建输入文件列表|Created input file list: {input_list_file}")

        # 构建KMC命令|Build KMC command
        # KMC用法: kmc [options] <@input_list> <output_name> <working_directory>
        # 注意：使用绝对路径以避免工作目录问题|Note: Use absolute paths to avoid working directory issues
        cmd_parts = [
            self.config.get_kmc_bin(),
            f'-k{self.config.kmer_size}',
            f'-ci{self.config.min_count}',
            f'-t{self.config.threads}'
        ]

        # 添加最大计数参数|Add max count parameter
        if self.config.max_count is not None:
            cmd_parts.append(f'-cs{self.config.max_count}')

        # 添加内存限制|Add memory limit
        if self.config.memory_limit:
            cmd_parts.append(f'-m{self.config.memory_limit}')

        # 输入文件列表（使用绝对路径和@前缀）|Input file list (use absolute path and @ prefix)
        abs_input_list = str(input_list_file.resolve())
        cmd_parts.append(f"@{abs_input_list}")

        # 输出文件名（只是文件名）|Output file name (filename only)
        cmd_parts.append(kmc_output_name)

        # 工作目录（使用绝对路径）|Working directory (use absolute path)
        abs_tmp_dir = str(Path(tmp_dir).resolve())
        cmd_parts.append(abs_tmp_dir)

        cmd = ' '.join(cmd_parts)

        # 执行命令|Execute command
        success = self.cmd_runner.run(cmd, f"统计样本|Counting sample {sample_name}")

        if success:
            # 移动生成的文件到目标目录|Move generated files to target directory
            self._move_kmc_output(tmp_dir, output_dir, sample_name)
            self.logger.info(f"样本统计完成|Sample counting completed: {sample_name}")
            # 清理临时文件|Clean temporary files
            self._cleanup_tmp(tmp_dir)
        else:
            self.logger.error(f"样本统计失败|Sample counting failed: {sample_name}")

        return success

    def count_samples(self) -> dict:
        """统计多个样本的k-mer|Count k-mers for multiple samples"""
        # 使用sample_finder查找样本|Use sample_finder to find samples
        samples = self.sample_finder.find_samples()

        if not samples:
            self.logger.error("未找到样本文件|No sample files found")
            return {}

        # 验证样本|Validate samples
        if not self.sample_finder.validate_samples(samples):
            return {}

        self.logger.info(f"开始统计 {len(samples)} 个样本|Starting to count {len(samples)} samples")

        results = {}

        for i, (sample_name, file_paths) in enumerate(samples, 1):
            self.logger.info(f"处理进度|Progress: {i}/{len(samples)}: {sample_name}")

            # KMC支持多个输入文件|KMC supports multiple input files
            input_files_str = ' '.join(str(f) for f in file_paths)

            success = self.count_sample(input_files_str, sample_name)
            results[sample_name] = success

        # 统计结果|Statistics
        success_count = sum(1 for v in results.values() if v)
        self.logger.info(
            f"统计完成|Counting completed: "
            f"{success_count}/{len(results)} 个样本成功|samples successful"
        )

        return results

    def run(self) -> dict:
        """运行k-mer统计|Run k-mer counting"""
        return self.count_samples()

    def union_databases(self, db_names: List[str], output_name: str = 'global_kmers') -> bool:
        """合并多个KMC数据库|Union multiple KMC databases

        使用pipeline方式逐个合并：db1 ∪ db2 ∪ db3 ... ∪ dbN
        Uses pipeline approach to union databases one by one: db1 ∪ db2 ∪ db3 ... ∪ dbN

        注意：CommandRunner的工作目录是kmc_db_path，所以使用相对路径（只有文件名）
        Note: CommandRunner working dir is kmc_db_path, so use relative paths (filenames only)
        """
        self.logger.info(f"开始合并 {len(db_names)} 个数据库|Starting to union {len(db_names)} databases")

        if len(db_names) < 2:
            self.logger.error("至少需要2个数据库才能合并|At least 2 databases required for union")
            return False

        # 使用文件名（相对路径），因为CommandRunner的working_dir已经是kmc_db_path
        # Use filenames (relative paths) because CommandRunner's working_dir is already kmc_db_path
        db_names_only = db_names

        # 线程参数|Thread parameter
        threads_param = f"-t{self.config.threads}"

        # 使用pipeline方式逐个合并|Use pipeline approach to union one by one
        # 初始：使用前两个数据库|Initial: union first two databases
        current_union = f".temp_union_0"

        cmd = f"{self.config.get_kmc_tools_bin()} {threads_param} simple {db_names_only[0]} {db_names_only[1]} union {current_union}"
        if not self.cmd_runner.run(cmd, f"合并数据库 1-2|Union databases 1-2"):
            return False

        # 逐个合并剩余数据库|Union remaining databases one by one
        for i in range(2, len(db_names_only)):
            next_union = f".temp_union_{i}"
            cmd = f"{self.config.get_kmc_tools_bin()} {threads_param} simple {current_union} {db_names_only[i]} union {next_union}"

            if not self.cmd_runner.run(cmd, f"合并数据库 {i+1}|Union database {i+1}"):
                # 清理临时文件|Clean temporary files
                self._cleanup_temp_unions(i)
                return False

            # 删除旧的临时文件|Remove old temp file
            try:
                pre_file = str(self.config.kmc_db_path / f"{current_union}.kmc_pre")
                suf_file = str(self.config.kmc_db_path / f"{current_union}.kmc_suf")
                if os.path.exists(pre_file):
                    os.remove(pre_file)
                if os.path.exists(suf_file):
                    os.remove(suf_file)
            except Exception as e:
                self.logger.warning(f"删除临时文件失败|Failed to remove temp file: {e}")

            current_union = next_union

        # 最终输出|Final output
        final_output = str(self.config.output_path / output_name)

        # 移动最终结果到目标位置|Move final result to target location
        try:
            current_pre = str(self.config.kmc_db_path / f"{current_union}.kmc_pre")
            current_suf = str(self.config.kmc_db_path / f"{current_union}.kmc_suf")
            final_pre = final_output + ".kmc_pre"
            final_suf = final_output + ".kmc_suf"

            if os.path.exists(current_pre):
                shutil.move(current_pre, final_pre)
            if os.path.exists(current_suf):
                shutil.move(current_suf, final_suf)
        except Exception as e:
            self.logger.error(f"移动最终结果失败|Failed to move final result: {e}")
            return False

        self.logger.info(f"数据库合并完成|Database union completed: {output_name}")
        return True

    def union_two_databases(self, db1: str, db2: str, output_name: str) -> bool:
        """增量式合并两个数据库|Incremental union of two databases

        将 db2 合并到 db1 中，生成新的全局数据库
        Merge db2 into db1, generating new global database

        Args:
            db1: 已有的全局数据库（如 global_kmers）|Existing global database
            db2: 新样本数据库|New sample database
            output_name: 输出数据库名称|Output database name

        Returns:
            是否成功|Success or not
        """
        self.logger.info(f"增量式合并数据库|Incremental union: {db1} + {db2} → {output_name}")

        # 线程参数|Thread parameter
        threads_param = f"-t{self.config.threads}"

        # 处理路径：global_kmers在output_path，样本数据库在kmc_db_path
        # Handle paths: global_kmers in output_path, sample databases in kmc_db_path
        if db1 == 'global_kmers':
            # global_kmers在父目录（output_path），使用相对路径|global_kmers in parent directory (output_path), use relative path
            db1_path = '../global_kmers'
        elif db1.startswith('global_kmers_temp_'):
            # 临时数据库也在父目录|Temp database also in parent directory
            db1_path = f'../{db1}'
        else:
            # 样本数据库在当前目录（kmc_db_path）|Sample database in current directory (kmc_db_path)
            db1_path = db1

        # db2是样本名，在当前目录|db2 is sample name, in current directory
        db2_path = db2

        # 输出到当前目录|Output to current directory
        output_path = output_name

        # 执行合并命令|Execute union command
        cmd = f"{self.config.get_kmc_tools_bin()} {threads_param} simple {db1_path} {db2_path} union {output_path}"

        if not self.cmd_runner.run(cmd, f"增量合并数据库|Incremental union: {db1} + {db2}"):
            return False

        # 移动最终结果到目标位置|Move final result to target location
        final_output = str(self.config.output_path / output_name)

        try:
            current_pre = str(self.config.kmc_db_path / f"{output_name}.kmc_pre")
            current_suf = str(self.config.kmc_db_path / f"{output_name}.kmc_suf")
            final_pre = f"{final_output}.kmc_pre"
            final_suf = f"{final_output}.kmc_suf"

            if os.path.exists(current_pre):
                shutil.move(current_pre, final_pre)
            if os.path.exists(current_suf):
                shutil.move(current_suf, final_suf)
        except Exception as e:
            self.logger.error(f"移动最终结果失败|Failed to move final result: {e}")
            return False

        self.logger.info(f"增量合并完成|Incremental union completed: {output_name}")
        return True

    def _cleanup_temp_unions(self, last_index: int):
        """清理临时合并文件|Clean temporary union files"""
        try:
            for i in range(last_index):
                temp_path = str(self.config.kmc_db_path / f".temp_union_{i}")
                if os.path.exists(temp_path + ".kmc_pre"):
                    os.remove(temp_path + ".kmc_pre")
                if os.path.exists(temp_path + ".kmc_suf"):
                    os.remove(temp_path + ".kmc_suf")
        except Exception as e:
            self.logger.warning(f"清理临时文件失败|Failed to clean temp files: {e}")

    def _move_kmc_output(self, tmp_dir: str, output_dir: Path, sample_name: str):
        """验证KMC输出文件|Verify KMC output files

        注意：KMC会将输出文件生成在当前工作目录中，而不是在working_directory中。
        由于CommandRunner设置了cwd=output_dir，文件会直接生成在output_dir中。
        因此这里只需要验证文件是否存在，不需要移动。

        Note: KMC generates output files in current working directory, not in working_directory.
        Since CommandRunner sets cwd=output_dir, files are generated directly in output_dir.
        So this method only verifies files exist, no moving needed.
        """
        # 验证文件是否已经生成在输出目录中|Verify files are already in output directory
        kmc_pre = output_dir / f"{sample_name}.kmc_pre"
        kmc_suf = output_dir / f"{sample_name}.kmc_suf"

        if kmc_pre.exists() and kmc_suf.exists():
            self.logger.debug(f"KMC输出文件已生成|KMC output files generated: {kmc_pre}, {kmc_suf}")
        else:
            self.logger.error(f"KMC输出文件不完整|KMC output files incomplete:")
            if not kmc_pre.exists():
                self.logger.error(f"  缺少|Missing: {kmc_pre}")
            if not kmc_suf.exists():
                self.logger.error(f"  缺少|Missing: {kmc_suf}")
            raise FileNotFoundError(f"KMC output files not found in {output_dir}")

    def _cleanup_tmp(self, tmp_dir: str):
        """清理临时文件|Clean temporary files"""
        try:
            import shutil
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
                self.logger.debug(f"已清理临时文件|Cleaned temporary files: {tmp_dir}")
        except Exception as e:
            self.logger.warning(f"清理临时文件失败|Failed to clean temporary files: {e}")

    def get_database_info(self, sample_name: str) -> dict:
        """获取KMC数据库信息|Get KMC database information"""
        db_path = self.config.get_sample_db_path(sample_name)

        info = {
            'sample_name': sample_name,
            'database_exists': False,
            'kmer_count': 0,
            'total_kmers': 0
        }

        # 检查数据库文件是否存在|Check if database files exist
        pre_file = f"{db_path}.kmc_pre"
        suf_file = f"{db_path}.kmc_suf"

        if os.path.exists(pre_file) and os.path.exists(suf_file):
            info['database_exists'] = True

            # 使用kmc_tools获取统计信息|Use kmc_tools to get statistics
            # 这里需要调用KMC API或解析输出|Need to call KMC API or parse output
            # 暂时标记为需要实现|Mark as to be implemented
            pass

        return info
