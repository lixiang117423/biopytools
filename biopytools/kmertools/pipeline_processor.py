"""KmerPAV流程处理器|KmerPAV Pipeline Processor

处理kmtricks流程：pipeline -> aggregate -> import to rocksdb
Process kmtricks pipeline: pipeline -> aggregate -> import to rocksdb
"""

import os
import re
import gzip
import shutil
from pathlib import Path
from typing import Optional, Tuple, Dict
from . import rocksdb_importer
from .utils import build_conda_command


class PipelineProcessor:
    """流程处理器|Pipeline Processor

    负责执行kmtricks pipeline、aggregate和rocksdb导入
    Responsible for executing kmtricks pipeline, aggregate and rocksdb import
    """

    def __init__(self, config, logger, cmd_runner):
        """初始化流程处理器|Initialize pipeline processor

        Args:
            config: KmerPAVConfig配置对象|KmerPAVConfig object
            logger: 日志器|Logger instance
            cmd_runner: 命令执行器|Command runner instance
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    @staticmethod
    def has_invalid_chars(sample_name: str) -> bool:
        """检测样品名是否包含kmtricks不支持的特殊字符|Check if sample name contains invalid characters

        kmtricks FOF格式只支持 [A-Za-z0-9_-] 这些字符
        kmtricks FOF format only supports [A-Za-z0-9_-]

        Args:
            sample_name: 样品名|Sample name

        Returns:
            bool: 如果包含特殊字符返回True|Return True if contains invalid characters
        """
        # 检查是否包含非字母数字下划线连字符的字符
        # Check for characters other than alphanumeric, underscore, hyphen
        return bool(re.search(r'[^A-Za-z0-9_-]', sample_name))

    @staticmethod
    def sanitize_sample_name(sample_name: str) -> str:
        """清理样品名，将特殊字符替换为下划线|Sanitize sample name by replacing special chars with underscore

        Args:
            sample_name: 原始样品名|Original sample name

        Returns:
            str: 清理后的样品名|Sanitized sample name
        """
        # 先替换一些常见的特殊字符为更易读的形式
        # First replace some common special characters with more readable forms
        char_map = {
            '①': '1', '②': '2', '③': '3', '④': '4', '⑤': '5',
            '⑥': '6', '⑦': '7', '⑧': '8', '⑨': '9', '⑩': '10',
            '⑪': '11', '⑫': '12', '⑬': '13', '⑭': '14', '⑮': '15',
            '⑯': '16', '⑰': '17', '⑱': '18', '⑲': '19', '⑳': '20',
            '⓵': '1', '⓶': '2', '⓷': '3', '⓸': '4', '⓹': '5',
            '⓺': '6', '⓻': '7', '⓼': '8', '⓽': '9', '⓾': '10',
        }
        for old_char, new_char in char_map.items():
            sample_name = sample_name.replace(old_char, new_char)

        # 将其他不属于 [A-Za-z0-9_-] 的字符替换为下划线
        # Replace other characters not in [A-Za-z0-9_-] with underscore
        sample_name = re.sub(r'[^A-Za-z0-9_-]', '_', sample_name)

        return sample_name

    def generate_fof_file(self) -> Optional[str]:
        """生成FOF文件|Generate File-of-Files file

        支持单末端和双末端测序数据
        Support both single-end and paired-end sequencing data
        自动处理包含特殊字符的样品名|Automatically handle sample names with special characters

        Returns:
            str: FOF文件路径|FOF file path, 失败返回None|None if failed
        """
        self.logger.info("生成FOF文件|Generating FOF file")

        input_dir = Path(self.config.input_dir)
        fof_path = self.config.output_path / "input.fof"
        cleaned_dir = self.config.output_path / "cleaned_data"
        mapping_file = self.config.output_path / "sample_mapping.txt"

        if not input_dir.exists():
            self.logger.error(f"输入目录不存在|Input directory not found: {input_dir}")
            return None

        # 首先尝试查找双末端数据（paired-end）|Try to find paired-end data first
        r1_files = sorted(input_dir.glob("*_1.clean.fq.gz"))
        if not r1_files:
            r1_files = sorted(input_dir.glob("*_1.filter.fq.gz"))
        if not r1_files:
            r1_files = sorted(input_dir.glob("*_1.fq.gz"))
        if not r1_files:
            r1_files = sorted(input_dir.glob("*_1.fastq.gz"))

        # 如果找到双末端数据|If paired-end data found
        if r1_files:
            self.logger.info(f"检测到双末端数据|Detected paired-end data: {len(r1_files)} samples")

            # 扫描并分类样品：正常 vs 异常|Scan and classify samples: normal vs invalid
            invalid_samples = []
            sample_info = []  # 存储 (原始样品名, r1路径, r2路径, 清理后样品名)

            for r1_file in r1_files:
                # 确定对应的R2文件|Determine corresponding R2 file
                r2_file = r1_file.name.replace("_1.", "_2.")
                r2_path = r1_file.parent / r2_file

                if not r2_path.exists():
                    self.logger.warning(f"未找到对应的R2文件|R2 file not found: {r2_path}")
                    continue

                # 样本名（去除后缀）|Sample name (remove suffix)
                sample_name = r1_file.name.replace("_1.clean.fq.gz", "")
                sample_name = sample_name.replace("_1.filter.fq.gz", "")
                sample_name = sample_name.replace("_1.fq.gz", "")
                sample_name = sample_name.replace("_1.fastq.gz", "")

                sanitized_name = self.sanitize_sample_name(sample_name)

                if self.has_invalid_chars(sample_name):
                    invalid_samples.append({
                        'original': sample_name,
                        'sanitized': sanitized_name,
                        'r1': r1_file,
                        'r2': r2_path
                    })
                    self.logger.warning(f"检测到异常样品名|Detected invalid sample name: {sample_name} -> {sanitized_name}")

                sample_info.append((sample_name, r1_file, r2_path, sanitized_name))

            # 如果有异常样品，创建 cleaned_data 目录并拷贝文件|If invalid samples exist, create cleaned_data and copy files
            mapping_lines = []
            if invalid_samples:
                self.logger.info(f"发现 {len(invalid_samples)} 个异常样品，进行处理中|Found {len(invalid_samples)} invalid samples, processing...")
                cleaned_dir.mkdir(exist_ok=True)

                for sample in invalid_samples:
                    orig_name = sample['original']
                    sane_name = sample['sanitized']
                    r1_src = sample['r1']
                    r2_src = sample['r2']

                    # 拷贝并重命名 R1 文件|Copy and rename R1 file
                    r1_name_new = r1_src.name.replace(orig_name, sane_name)
                    r1_dst = cleaned_dir / r1_name_new
                    shutil.copy2(r1_src, r1_dst)
                    self.logger.debug(f"拷贝文件|Copying file: {r1_src} -> {r1_dst}")

                    # 拷贝并重命名 R2 文件|Copy and rename R2 file
                    r2_name_new = r2_src.name.replace(orig_name, sane_name)
                    r2_dst = cleaned_dir / r2_name_new
                    shutil.copy2(r2_src, r2_dst)
                    self.logger.debug(f"拷贝文件|Copying file: {r2_src} -> {r2_dst}")

                    # 记录映射信息|Record mapping info
                    r1_rel_orig = os.path.relpath(r1_src, fof_path.parent)
                    r2_rel_orig = os.path.relpath(r2_src, fof_path.parent)
                    r1_rel_new = f"cleaned_data/{r1_name_new}"
                    r2_rel_new = f"cleaned_data/{r2_name_new}"
                    mapping_lines.append(f"{sane_name}\t{orig_name}\tR1_original: {r1_rel_orig}\tR1_cleaned: {r1_rel_new}\n")
                    mapping_lines.append(f"{sane_name}\t{orig_name}\tR2_original: {r2_rel_orig}\tR2_cleaned: {r2_rel_new}\n")

                # 生成映射文件|Generate mapping file
                with open(mapping_file, 'w') as mf:
                    mf.write("# 样品名映射文件|Sample name mapping file\n")
                    mf.write("# 格式|Format: 清理后名称\t原始名称\tR1原始路径\tR1清理后路径\n")
                    mf.write("# Format|Format: Sanitized_Name\tOriginal_Name\tR1_Original_Path\tR1_Cleaned_Path\n")
                    mf.write("".join(mapping_lines))

                self.logger.info(f"样品映射文件已生成|Sample mapping file generated: {mapping_file}")

            # 生成 FOF 文件|Generate FOF file
            with open(fof_path, 'w') as f:
                for sample_name, r1_file, r2_path, sanitized_name in sample_info:
                    if self.has_invalid_chars(sample_name):
                        # 使用清理后的文件路径|Use cleaned file path
                        r1_name_new = r1_file.name.replace(sample_name, sanitized_name)
                        r2_name_new = r2_path.name.replace(sample_name, sanitized_name)
                        r1_rel = f"cleaned_data/{r1_name_new}"
                        r2_rel = f"cleaned_data/{r2_name_new}"
                        f.write(f"{sanitized_name}: {r1_rel} ; {r2_rel}\n")
                    else:
                        # 使用相对路径（相对于当前工作目录）|Use relative path (relative to current working directory)
                        # kmtricks 不支持中文路径，所以使用相对路径避免中文目录名|kmtricks doesn't support Chinese paths, use relative path
                        r1_rel = os.path.relpath(r1_file, os.getcwd())
                        r2_rel = os.path.relpath(r2_path, os.getcwd())
                        f.write(f"{sample_name}: {r1_rel} ; {r2_rel}\n")

            self.logger.info(f"FOF文件已生成（双末端）|FOF file generated (paired-end): {fof_path}")
            # 保存样本数到配置|Save sample count to config
            self.config._num_samples = len(r1_files)
            return str(fof_path)

        # 如果没有找到双末端数据，尝试查找单末端数据（single-end）|If no paired-end data, try single-end
        self.logger.info("未检测到双末端数据，尝试查找单末端数据|No paired-end data detected, trying single-end data")

        # 查找单末端fastq文件|Find single-end fastq files
        # 排除已经是双末端的文件（包含_1或_2）|Exclude paired-end files (containing _1 or _2)
        all_gz_files = sorted(input_dir.glob("*.fq.gz"))
        single_end_files = [f for f in all_gz_files if "_1." not in f.name and "_2." not in f.name]

        if not single_end_files:
            all_gz_files = sorted(input_dir.glob("*.fastq.gz"))
            single_end_files = [f for f in all_gz_files if "_1." not in f.name and "_2." not in f.name]

        if not single_end_files:
            self.logger.error("未找到FASTQ文件（既未找到双末端也未找到单末端数据）|No FASTQ files found (neither paired-end nor single-end)")
            return None

        self.logger.info(f"检测到单末端数据|Detected single-end data: {len(single_end_files)} samples")

        # 扫描并分类样品|Scan and classify samples
        invalid_samples = []
        sample_info = []

        for fastq_file in single_end_files:
            # 样本名（去除后缀）|Sample name (remove suffix)
            sample_name = fastq_file.name.replace(".fq.gz", "").replace(".fastq.gz", "")
            sanitized_name = self.sanitize_sample_name(sample_name)

            if self.has_invalid_chars(sample_name):
                invalid_samples.append({
                    'original': sample_name,
                    'sanitized': sanitized_name,
                    'file': fastq_file
                })
                self.logger.warning(f"检测到异常样品名|Detected invalid sample name: {sample_name} -> {sanitized_name}")

            sample_info.append((sample_name, fastq_file, sanitized_name))

        # 处理异常样品（单末端）|Process invalid samples (single-end)
        mapping_lines = []
        if invalid_samples:
            self.logger.info(f"发现 {len(invalid_samples)} 个异常样品，进行处理中|Found {len(invalid_samples)} invalid samples, processing...")
            cleaned_dir.mkdir(exist_ok=True)

            for sample in invalid_samples:
                orig_name = sample['original']
                sane_name = sample['sanitized']
                file_src = sample['file']

                # 拷贝并重命名文件|Copy and rename file
                file_name_new = file_src.name.replace(orig_name, sane_name)
                file_dst = cleaned_dir / file_name_new
                shutil.copy2(file_src, file_dst)
                self.logger.debug(f"拷贝文件|Copying file: {file_src} -> {file_dst}")

                # 记录映射信息|Record mapping info
                file_rel_orig = os.path.relpath(file_src, fof_path.parent)
                file_rel_new = f"cleaned_data/{file_name_new}"
                mapping_lines.append(f"{sane_name}\t{orig_name}\tOriginal: {file_rel_orig}\tCleaned: {file_rel_new}\n")

            # 生成映射文件|Generate mapping file
            with open(mapping_file, 'w') as mf:
                mf.write("# 样品名映射文件|Sample name mapping file\n")
                mf.write("# 格式|Format: 清理后名称\t原始名称\t原始路径\t清理后路径\n")
                mf.write("# Format|Format: Sanitized_Name\tOriginal_Name\tOriginal_Path\tCleaned_Path\n")
                mf.write("".join(mapping_lines))

            self.logger.info(f"样品映射文件已生成|Sample mapping file generated: {mapping_file}")

        # 生成 FOF 文件（单末端）|Generate FOF file (single-end)
        with open(fof_path, 'w') as f:
            for sample_name, fastq_file, sanitized_name in sample_info:
                if self.has_invalid_chars(sample_name):
                    # 使用清理后的文件路径|Use cleaned file path
                    file_name_new = fastq_file.name.replace(sample_name, sanitized_name)
                    file_rel = f"cleaned_data/{file_name_new}"
                    f.write(f"{sanitized_name}: {file_rel}\n")
                else:
                    # 使用相对路径（相对于当前工作目录）|Use relative path (relative to current working directory)
                    # kmtricks 不支持中文路径，所以使用相对路径避免中文目录名|kmtricks doesn't support Chinese paths, use relative path
                    fastq_rel = os.path.relpath(fastq_file, os.getcwd())
                    f.write(f"{sample_name}: {fastq_rel}\n")

        self.logger.info(f"FOF文件已生成（单末端）|FOF file generated (single-end): {fof_path}")
        # 保存样本数到配置|Save sample count to config
        self.config._num_samples = len(single_end_files)
        return str(fof_path)

    def run_kmtricks_pipeline(self, fof_file: str) -> bool:
        """运行kmtricks pipeline|Run kmtricks pipeline

        Args:
            fof_file: FOF文件路径|FOF file path

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("运行kmtricks pipeline|Running kmtricks pipeline")

        # 自动计算分区数（如果需要）|Auto calculate partitions if needed
        if self.config.nb_partitions == 0 and self.config._num_samples > 0:
            optimal_partitions = self.config.calculate_nb_partitions(self.config._num_samples)
            if optimal_partitions > 0:
                self.logger.info(f"自动计算分区数|Auto-calculated partitions: {optimal_partitions}")
                self.logger.info(f"样本数|Samples: {self.config._num_samples}, 线程数|Threads: {self.config.threads}")
                # 临时覆盖配置值，用于本次运行
                self.config.nb_partitions = optimal_partitions

        # 构建命令|Build command
        run_dir = Path(self.config.run_dir)

        # 删除已存在的运行目录|Remove existing run directory
        if run_dir.exists():
            self.logger.info(f"删除已存在的运行目录|Removing existing run directory: {run_dir}")
            shutil.rmtree(run_dir)

        # 构建参数列表|Build argument list
        args = [
            "pipeline",
            "--file", fof_file,
            "--run-dir", str(run_dir),
            "--kmer-size", str(self.config.kmer_size),
            "--hard-min", str(self.config.hard_min),
            "--recurrence-min", str(self.config.recurrence_min),
            "--mode", "kmer:pa:bin",  # presence/absence binary
            "--cpr",  # 启用压缩|Enable compression
            "-t", str(self.config.threads)
        ]

        # 如果指定了分区数，添加参数|Add nb-partitions if specified
        if self.config.nb_partitions > 0:
            args.extend(["--nb-partitions", str(self.config.nb_partitions)])
            self.logger.info(f"使用分区数|Using partitions: {self.config.nb_partitions}")
        elif self.config.nb_partitions == -1:
            self.logger.info("使用kmtricks默认分区数|Using kmtricks default partitions")

        # 使用conda命令包装器|Use conda command wrapper
        cmd = build_conda_command(self.config.kmtricks_path, args)

        # 执行命令|Execute command
        success = self.cmd_runner.run_command(
            cmd,
            description="执行kmtricks pipeline|Executing kmtricks pipeline"
        )

        if success:
            self.logger.info("kmtricks pipeline完成|kmtricks pipeline completed")
        else:
            self.logger.error("kmtricks pipeline失败|kmtricks pipeline failed")

        return success

    def run_kmtricks_aggregate(self) -> Optional[str]:
        """运行kmtricks aggregate|Run kmtricks aggregate

        Returns:
            str: 聚合文件路径|Aggregated file path, 失败返回None|None if failed
        """
        self.logger.info("运行kmtricks aggregate|Running kmtricks aggregate")

        run_dir = Path(self.config.run_dir)
        output_file = self.config.output_path / "kmer_matrix.txt"

        # 构建参数列表|Build argument list
        args = [
            "aggregate",
            "--run-dir", str(run_dir),
            "--pa-matrix", "kmer",
            "--format", "text",
            "--cpr-in",
            "-t", str(self.config.threads),
            "--output", str(output_file)
        ]

        # 使用conda命令包装器|Use conda command wrapper
        cmd = build_conda_command(self.config.kmtricks_path, args)

        success = self.cmd_runner.run_command(
            cmd,
            description="执行kmtricks aggregate|Executing kmtricks aggregate"
        )

        if success:
            self.logger.info(f"kmtricks aggregate完成|kmtricks aggregate completed: {output_file}")
            return str(output_file)
        else:
            self.logger.error("kmtricks aggregate失败|kmtricks aggregate failed")
            return None

    def compress_matrix_file(self, matrix_file: str) -> str:
        """压缩矩阵文件|Compress matrix file

        Args:
            matrix_file: 矩阵文件路径|Matrix file path

        Returns:
            str: 压缩文件路径|Compressed file path
        """
        self.logger.info("压缩矩阵文件|Compressing matrix file")

        gz_file = f"{matrix_file}.gz"

        # 如果已存在压缩文件，先删除|Remove existing compressed file if it exists
        if os.path.exists(gz_file):
            self.logger.info(f"删除已存在的压缩文件|Removing existing compressed file: {gz_file}")
            os.remove(gz_file)

        # 使用conda命令包装器构建bgzip命令|Use conda command wrapper to build bgzip command
        args = ["-k", matrix_file, "-@", str(self.config.threads)]
        cmd = build_conda_command("bgzip", args)

        # 检查bgzip是否可用|Check if bgzip is available
        if "conda" not in cmd and not self.cmd_runner.check_command_exists("bgzip"):
            self.logger.warning("bgzip未找到，使用gzip压缩|bgzip not found, using gzip")
            cmd = ["gzip", "-c", matrix_file]
            with open(gz_file, 'wb') as f_out:
                subprocess.run(cmd, stdout=f_out, check=True)
        else:
            self.cmd_runner.run_command(cmd, description="压缩矩阵文件|Compressing matrix file")

        self.logger.info(f"矩阵文件已压缩|Matrix file compressed: {gz_file}")
        return gz_file

    def extract_header_from_fof(self, fof_file: str) -> Optional[str]:
        """从FOF文件提取样本名|Extract sample names from FOF file

        Args:
            fof_file: FOF文件路径|FOF file path

        Returns:
            str: header文件路径|Header file path, 失败返回None|None if failed
        """
        self.logger.info("提取样本名|Extracting sample names")

        header_file = self.config.output_path / "header.txt"

        try:
            with open(fof_file, 'r') as f_in, open(header_file, 'w') as f_out:
                for line in f_in:
                    # 格式: sample_name: r1_file ; r2_file
                    # Format: sample_name: r1_file ; r2_file
                    sample_name = line.split(':')[0].strip()
                    f_out.write(f"{sample_name}\n")

            self.logger.info(f"样本名已提取|Sample names extracted: {header_file}")
            return str(header_file)

        except Exception as e:
            self.logger.error(f"提取样本名失败|Failed to extract sample names: {e}")
            return None

    def import_to_rocksdb(self, matrix_gz: str, header_file: str) -> bool:
        """导入矩阵到RocksDB|Import matrix to RocksDB

        Args:
            matrix_gz: 压缩的矩阵文件|Compressed matrix file
            header_file: 样本名文件|Sample name file

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("导入矩阵到RocksDB|Importing matrix to RocksDB")

        # 检查是否已有rocksdb python库|Check if rocksdb python library exists
        try:
            import rocksdb
        except ImportError:
            self.logger.error("rocksdb Python库未安装|rocksdb Python library not installed")
            self.logger.error("请安装: pip install python-rocksdb|Please install: pip install python-rocksdb")
            return False

        rocksdb_path = str(self.config.rocksdb_path)

        # 调用rocksdb_importer模块进行导入|Call rocksdb_importer module to import
        try:
            self.logger.info(f"开始导入kmer矩阵到RocksDB|Starting kmer matrix import to RocksDB: {rocksdb_path}")

            # 调用导入函数|Call import function
            success = rocksdb_importer.import_gz_to_rocksdb(
                gz_file_path=str(matrix_gz),
                db_path=rocksdb_path,
                input_delimiter=' ',  # kmer_matrix使用空格分隔|kmer_matrix uses space delimiter
                value_join_char='',  # 直接连接|Direct concatenation
                batch_size=20000,
                bloom_bits=15,
                force_overwrite=True,  # 覆盖已存在的数据库|Overwrite existing database
                header_file_path=str(header_file),
                header_db_key=self.config.header_db_key,
                header_storage_delimiter='\t'
            )

            if success:
                self.logger.info(f"RocksDB导入成功|RocksDB import successful: {rocksdb_path}")
                return True
            else:
                self.logger.error("RocksDB导入失败|RocksDB import failed")
                return False

        except Exception as e:
            self.logger.error(f"RocksDB导入异常|RocksDB import error: {e}")
            return False

    def generate_matrix_with_header(self, matrix_file: str, header_file: str) -> Optional[str]:
        """生成带表头的矩阵文件|Generate matrix file with header

        Args:
            matrix_file: 原始矩阵文件路径|Original matrix file path
            header_file: 样本名文件路径|Sample name file path

        Returns:
            str: 带表头的矩阵文件路径|Matrix file path with header
        """
        self.logger.info("生成带表头的矩阵文件|Generating matrix file with header")

        output_file = self.config.output_path / "kmer_matrix_with_header.txt"

        try:
            # 读取header并合并成一行|Read header and merge into one line
            with open(header_file, 'r') as f:
                sample_names = [line.strip() for line in f if line.strip()]

            # 第一列是KMER，后面是样本名（用空格分隔）|First column is KMER, followed by sample names (space-separated)
            header_line = "KMER " + " ".join(sample_names) + "\n"

            # 写入输出文件|Write to output file
            with open(output_file, 'w') as f_out:
                f_out.write(header_line)

                # 追加原始矩阵内容|Append original matrix content
                with open(matrix_file, 'r') as f_in:
                    # 直接写入原始内容（空格分隔）|Write original content as-is (space-delimited)
                    for line in f_in:
                        f_out.write(line)

            self.logger.info(f"带表头的矩阵文件已生成|Matrix file with header generated: {output_file}")
            return str(output_file)

        except Exception as e:
            self.logger.error(f"生成带表头的矩阵文件失败|Failed to generate matrix file with header: {e}")
            return None

    def run_full_pipeline(self) -> bool:
        """运行完整流程|Run full pipeline

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("开始KmerPAV流程|Starting KmerPAV pipeline")

        # 1. 生成FOF文件|Generate FOF file
        if self.config.fof_file:
            fof_file = self.config.fof_file
        else:
            fof_file = self.generate_fof_file()
            if not fof_file:
                return False

        # 2. 运行kmtricks pipeline|Run kmtricks pipeline
        if not self.run_kmtricks_pipeline(fof_file):
            return False

        # 3. 运行kmtricks aggregate|Run kmtricks aggregate
        matrix_file = self.run_kmtricks_aggregate()
        if not matrix_file:
            return False

        # 4. 压缩矩阵文件|Compress matrix file
        matrix_gz = self.compress_matrix_file(matrix_file)

        # 5. 提取样本名|Extract sample names
        if self.config.header_file:
            header_file = self.config.header_file
        else:
            header_file = self.extract_header_from_fof(fof_file)
            if not header_file:
                return False

        # 6. 导入到RocksDB|Import to RocksDB
        if not self.import_to_rocksdb(matrix_gz, header_file):
            return False

        # 7. 生成带表头的矩阵文件|Generate matrix file with header
        self.generate_matrix_with_header(matrix_file, header_file)

        self.logger.info("KmerPAV流程完成|KmerPAV pipeline completed")
        return True
