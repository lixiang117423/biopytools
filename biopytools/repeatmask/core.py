"""
重复序列屏蔽核心逻辑模块|Repeat Masking Core Logic Module
"""

import os
import shutil
from pathlib import Path
from typing import Optional


class RepeatModelerRunner:
    """RepeatModeler运行器|RepeatModeler Runner"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.db_name = None
        self.db_prefix = None  # 数据库文件前缀（完整路径）
        self.output_lib = None
        # 记录运行开始时的临时目录，用于清理
        self.temp_dirs_before = set()

    def _check_database_exists(self) -> bool:
        """检查数据库是否已存在|Check if database already exists"""
        # BuildDatabase 生成的文件列表
        db_extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.nog', '.njs', '.nni']

        for ext in db_extensions:
            if not (self.db_prefix.with_suffix(ext)).exists():
                return False
        return True

    def _cleanup_temp_dirs(self):
        """清理RepeatModeler生成的临时目录|Cleanup RepeatModeler temporary directories"""
        import glob

        # 查找当前目录下所有 RM_* 开头的目录
        rm_dirs = []
        for pattern in ['RM_*', 'RM_*.']:
            rm_dirs.extend(glob.glob(pattern))

        if not rm_dirs:
            return

        # 清理临时目录
        cleaned_count = 0
        for rm_dir in rm_dirs:
            try:
                if Path(rm_dir).is_dir():
                    shutil.rmtree(rm_dir)
                    cleaned_count += 1
                    self.logger.debug(f"清理临时目录|Cleaned temp directory: {rm_dir}")
            except Exception as e:
                self.logger.warning(f"清理临时目录失败|Failed to clean temp directory {rm_dir}: {e}")

        if cleaned_count > 0:
            self.logger.info(f"清理RepeatModeler临时目录|Cleaned {cleaned_count} RepeatModeler temporary directories")

    def build_database(self) -> bool:
        """构建RepeatModeler数据库|Build RepeatModeler database"""
        # 使用完整路径，确保数据库文件输出到output_path目录
        self.db_name = str(self.config.output_path / f"{self.config.base_name}_db")
        self.db_prefix = self.config.output_path / f"{self.config.base_name}_db"

        # 检查数据库是否已存在（断点续传）
        if self._check_database_exists():
            self.logger.info(f"数据库已存在，跳过构建|Database already exists, skipping: {self.db_prefix}")
            return True

        cmd = f"{self.config.builddatabase_path} -name {self.db_name} {self.config.genome}"

        log_file = self.config.output_path / "builddatabase.log"
        success = self.cmd_runner.run_with_log(
            cmd,
            log_file,
            f"构建RepeatModeler数据库|Building RepeatModeler database: {self.db_name}",
            timeout=None
        )

        if success:
            self.logger.info(f"RepeatModeler数据库创建成功|RepeatModeler database created: {self.db_prefix}")

        return success

    def run_repeatmodeler(self) -> bool:
        """运行RepeatModeler|Run RepeatModeler"""
        if not self.db_name:
            self.logger.error("数据库未创建，无法运行RepeatModeler|Database not created, cannot run RepeatModeler")
            return False

        # 库文件在数据库同目录下生成
        expected_lib = Path(str(self.db_prefix) + "-families.fa")

        # 检查库文件是否已存在（断点续传）
        if expected_lib.exists():
            self.logger.info(f"库文件已存在，跳过RepeatModeler|Library file already exists, skipping RepeatModeler: {expected_lib}")
            self.output_lib = expected_lib
            # 清理之前运行可能残留的临时目录|Cleanup potential residual temp directories from previous runs
            self._cleanup_temp_dirs()
            return True

        # 构建命令|Build command
        cmd = f"{self.config.repeatmodeler_path} -database {self.db_name} -threads {self.config.threads}"

        if self.config.use_ltr:
            cmd += " -LTRStruct"

        if self.config.modeler_quick:
            cmd += " -quick"

        log_file = self.config.output_path / "repeatmodeler.log"
        success = self.cmd_runner.run_with_log(
            cmd,
            log_file,
            "运行RepeatModeler识别重复序列|Running RepeatModeler to identify repeats",
            timeout=None
        )

        # 检查是否生成了库文件|Check if library file was generated
        # 库文件在数据库同目录下生成
        if expected_lib.exists():
            self.output_lib = expected_lib
            self.logger.info(f"RepeatModeler库文件生成|RepeatModeler library file generated: {expected_lib}")
            # 清理临时目录|Cleanup temporary directories
            self._cleanup_temp_dirs()
            return True
        else:
            self.logger.warning(f"未找到预期的库文件|Expected library file not found: {expected_lib}")
            return False

    def get_library_path(self) -> Optional[Path]:
        """获取生成的库文件路径|Get generated library file path"""
        return self.output_lib


class RepeatMaskerRunner:
    """RepeatMasker运行器|RepeatMasker Runner"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_files = {}

    def run_repeatmasker(self, library_path: Optional[Path] = None) -> bool:
        """
        运行RepeatMasker|Run RepeatMasker

        Args:
            library_path: 自定义重复序列库路径|Custom repeat library path
        """
        # RepeatMasker输出目录|RepeatMasker output directory
        output_dir = self.config.output_path / "repeatmasker_output"
        output_dir.mkdir(exist_ok=True)

        # 复制基因组文件到输出目录|Copy genome file to output directory
        genome_copy = output_dir / Path(self.config.genome).name
        if not genome_copy.exists():
            shutil.copy2(self.config.genome, genome_copy)

        # 检查是否已完成（断点续传）
        genome_name = Path(self.config.genome).name

        # 检查 masked 文件（支持新旧格式）|Check masked file (support both old and new formats)
        if genome_name.endswith('.fa'):
            base_name = genome_name[:-3]
            masked_file = (output_dir / f"{base_name}.masked.fa") if (output_dir / f"{base_name}.masked.fa").exists() else (output_dir / f"{genome_name}.masked")
        else:
            masked_file = output_dir / f"{genome_name}.masked"

        out_file = output_dir / f"{genome_name}.out"

        expected_outputs = {'masked': masked_file, 'out': out_file}

        # 检查关键输出文件是否都存在
        all_exist = all(f.exists() for f in expected_outputs.values())
        if all_exist:
            self.logger.info(f"RepeatMasker输出文件已存在，跳过屏蔽|RepeatMasker output files already exist, skipping masking")
            self._collect_output_files(output_dir)
            return True

        # 构建命令|Build command
        cmd_parts = [
            f"cd {output_dir} &&",
            f"{self.config.repeatmasker_path}"
        ]

        # 添加屏蔽模式参数|Add masking mode parameter
        if self.config.masking_mode == 'soft':
            cmd_parts.append("-xsmall")  # 小写|lowercase
        elif self.config.masking_mode == 'x':
            cmd_parts.append("-x")  # X字母|letter X
        # hard模式使用默认的N屏蔽|hard mode uses default N masking

        # 添加物种或库参数|Add species or library parameter
        if self.config.use_dfam and self.config.species:
            # 使用Dfam/Repbase数据库|Use Dfam/Repbase database
            cmd_parts.append(f"-species {self.config.species}")

        if library_path and library_path.exists():
            # 使用自定义库|Use custom library
            # 转换为绝对路径，避免cd后相对路径失效
            library_abs_path = library_path.resolve()
            cmd_parts.append(f"-lib {library_abs_path}")

        # 添加其他参数|Add other parameters
        cmd_parts.extend([
            f"-pa {self.config.threads}",
            "-gff",  # 生成GFF格式输出|Generate GFF format output
            "-dir", ".",  # 输出到当前目录（已cd到output_dir）
            genome_copy.name
        ])

        cmd = " ".join(cmd_parts)

        log_file = output_dir / "repeatmasker_run.log"
        success = self.cmd_runner.run_with_log(
            cmd,
            log_file,
            "运行RepeatMasker屏蔽重复序列|Running RepeatMasker to mask repeats",
            timeout=None
        )

        if success:
            # 重命名输出文件为更标准的格式|Rename output files to standard format
            self._rename_output_files(output_dir)
            # 收集输出文件|Collect output files
            self._collect_output_files(output_dir)

        return success

    def _rename_output_files(self, output_dir: Path):
        """复制输出文件为标准格式（保留原文件）|Copy output files to standard format (keep original)"""
        genome_name = Path(self.config.genome).name

        # 复制 .fa.masked 为 .masked.fa（保留原始文件）|Copy .fa.masked to .masked.fa (keep original)
        masked_file_old = output_dir / f"{genome_name}.masked"
        if genome_name.endswith('.fa'):
            # 去掉 .fa 后缀，添加 .masked.fa|Remove .fa suffix, add .masked.fa
            base_name = genome_name[:-3]
            masked_file_new = output_dir / f"{base_name}.masked.fa"

            if masked_file_old.exists() and not masked_file_new.exists():
                shutil.copy2(str(masked_file_old), str(masked_file_new))
                self.logger.info(f"创建屏蔽文件副本|Created masked file copy: {masked_file_old} -> {masked_file_new}")
            elif masked_file_new.exists():
                # 副本已存在|Copy already exists
                pass

    def _collect_output_files(self, output_dir: Path):
        """收集输出文件|Collect output files"""
        genome_name = Path(self.config.genome).name

        # 构建可能的输出文件路径|Build possible output file paths
        possible_outputs = {
            'out': output_dir / f"{genome_name}.out",
            'gff': output_dir / f"{genome_name}.out.gff",
            'tbl': output_dir / f"{genome_name}.tbl"
        }

        # 处理 masked 文件的特殊命名规则|Handle special naming for masked file
        if genome_name.endswith('.fa'):
            base_name = genome_name[:-3]
            # 优先使用新格式 .masked.fa，兼容旧格式 .fa.masked|Prefer new format .masked.fa, fallback to old .fa.masked
            masked_new = output_dir / f"{base_name}.masked.fa"
            masked_old = output_dir / f"{genome_name}.masked"
            possible_outputs['masked'] = masked_new if masked_new.exists() else masked_old
        else:
            possible_outputs['masked'] = output_dir / f"{genome_name}.masked"

        for file_type, file_path in possible_outputs.items():
            if file_path.exists():
                self.output_files[file_type] = file_path
                self.logger.info(f"RepeatMasker输出文件生成|RepeatMasker output file generated: {file_type} -> {file_path}")

    def get_output_files(self) -> dict:
        """获取输出文件字典|Get output files dictionary"""
        return self.output_files

    def get_masked_genome(self) -> Optional[Path]:
        """获取屏蔽基因组文件|Get masked genome file"""
        return self.output_files.get('masked')

    def get_annotation_file(self) -> Optional[Path]:
        """获取注释文件|Get annotation file"""
        return self.output_files.get('out')

    def get_gff_file(self) -> Optional[Path]:
        """获取GFF文件|Get GFF file"""
        return self.output_files.get('gff')


class RepeatMaskPipeline:
    """重复序列屏蔽流程主类|Repeat Masking Pipeline Main Class"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 初始化各个运行器|Initialize runners
        self.modeler_runner = RepeatModelerRunner(config, logger, cmd_runner)
        self.masker_runner = RepeatMaskerRunner(config, logger, cmd_runner)

        # 存储结果|Store results
        self.results = {
            'genome_stats': {},
            'modeler_library': None,
            'repeatmasker_files': {}
        }

    def build_denovo_library(self) -> Optional[Path]:
        """构建从头重复库|Build de novo repeat library"""
        if self.config.skip_modeler:
            self.logger.info("跳过RepeatModeler步骤|Skipping RepeatModeler step")
            return None

        self.logger.info("开始RepeatModeler流程|Starting RepeatModeler pipeline")

        # 构建数据库|Build database
        if not self.modeler_runner.build_database():
            self.logger.error("RepeatModeler数据库构建失败|RepeatModeler database construction failed")
            return None

        # 运行RepeatModeler|Run RepeatModeler
        if not self.modeler_runner.run_repeatmodeler():
            self.logger.error("RepeatModeler运行失败|RepeatModeler execution failed")
            return None

        library_path = self.modeler_runner.get_library_path()
        if library_path:
            self.results['modeler_library'] = library_path
            self.logger.info(f"从头重复库构建完成|De novo repeat library constructed: {library_path}")

        return library_path

    def run_masking(self, library_path: Optional[Path] = None) -> bool:
        """
        运行重复序列屏蔽|Run repeat masking

        Args:
            library_path: 重复序列库路径|Repeat library path
        """
        self.logger.info("开始RepeatMasker屏蔽流程|Starting RepeatMasker masking pipeline")

        if not self.masker_runner.run_repeatmasker(library_path):
            self.logger.error("RepeatMasker运行失败|RepeatMasker execution failed")
            return False

        self.results['repeatmasker_files'] = self.masker_runner.get_output_files()

        # 输出结果摘要|Output results summary
        self._print_summary()

        return True

    def _print_summary(self):
        """打印结果摘要|Print results summary"""
        self.logger.info("=" * 60)
        self.logger.info("重复序列屏蔽完成|Repeat masking completed")
        self.logger.info("=" * 60)

        if self.results['modeler_library']:
            self.logger.info(f"RepeatModeler库|RepeatModeler library: {self.results['modeler_library']}")

        masked_genome = self.masker_runner.get_masked_genome()
        if masked_genome:
            self.logger.info(f"屏蔽基因组|Masked genome: {masked_genome}")

        annotation = self.masker_runner.get_annotation_file()
        if annotation:
            self.logger.info(f"注释文件|Annotation file: {annotation}")

        gff_file = self.masker_runner.get_gff_file()
        if gff_file:
            self.logger.info(f"GFF文件|GFF file: {gff_file}")

        self.logger.info("=" * 60)

    def get_results(self) -> dict:
        """获取结果|Get results"""
        return self.results
