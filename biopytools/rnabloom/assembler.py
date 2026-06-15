"""
RNA-Bloom转录组组装器|RNA-Bloom Transcriptome Assembler
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, List
from .utils import CommandRunner, DependencyChecker


class TranscriptomeAssembler:
    """转录组从头组装器|De Novo Transcriptome Assembler"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """初始化组装器|Initialize assembler

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 检查依赖|Check dependencies
        self.dep_checker = DependencyChecker(logger)
        self._check_dependencies()

        # 确定RNA-Bloom JAR文件路径|Determine RNA-Bloom JAR file path
        self.rnabloom_jar = self._locate_rnabloom_jar()

    def _check_dependencies(self):
        """检查依赖软件|Check dependency software"""
        self.logger.info("检查依赖软件|Checking dependency software")

        # Java必需|Java is required
        if not self.dep_checker.check_java():
            raise RuntimeError("Java未安装，请安装Java 11或17|Java is not installed. Please install Java 11 or 17")

        # minimap2必需|minimap2 is required
        if not self.dep_checker.check_minimap2():
            raise RuntimeError(
                "minimap2未安装，请安装: conda install -c bioconda minimap2|"
                "minimap2 is not installed. Please install: conda install -c bioconda minimap2"
            )

        # ntCard可选|ntCard is optional
        self.dep_checker.check_ntcard()

    def _locate_rnabloom_jar(self) -> str:
        """定位RNA-Bloom JAR文件|Locate RNA-Bloom JAR file

        Returns:
            str: JAR文件路径或rnabloom命令|JAR file path or rnabloom command
        """
        rnabloom_path = self.config.rnabloom_path

        # 如果是默认值"rnabloom"，使用命令|If default "rnabloom", use command
        if rnabloom_path == "rnabloom":
            self.logger.info("使用RNA-Bloom命令|Using RNA-Bloom command")
            return "rnabloom"

        # 如果直接指向JAR文件|If directly pointing to JAR file
        if rnabloom_path.endswith('.jar'):
            if os.path.exists(rnabloom_path):
                self.logger.info(f"使用RNA-Bloom JAR文件|Using RNA-Bloom JAR file: {rnabloom_path}")
                return rnabloom_path
            else:
                raise FileNotFoundError(f"RNA-Bloom JAR文件不存在|JAR file not found: {rnabloom_path}")

        # 如果是目录，查找JAR文件|If directory, find JAR file
        if os.path.isdir(rnabloom_path):
            jar_path = os.path.join(rnabloom_path, "RNA-Bloom.jar")
            if os.path.exists(jar_path):
                self.logger.info(f"找到RNA-Bloom JAR文件|Found RNA-Bloom JAR file: {jar_path}")
                return jar_path
            else:
                raise FileNotFoundError(
                    f"在目录中未找到RNA-Bloom JAR文件|"
                    f"RNA-Bloom JAR not found in directory: {rnabloom_path}"
                )

        # 如果是可执行文件|If executable file
        if os.path.isfile(rnabloom_path):
            self.logger.info(f"使用RNA-Bloom可执行文件|Using RNA-Bloom executable: {rnabloom_path}")
            return rnabloom_path

        raise RuntimeError(f"无法确定RNA-Bloom路径|Cannot determine RNA-Bloom path: {rnabloom_path}")

    def build_command_args(self) -> List[str]:
        """构建RNA-Bloom命令参数|Build RNA-Bloom command arguments

        Returns:
            List[str]: 参数列表|Argument list
        """
        args = []

        # 输入文件参数|Input file parameters
        if self.config.cell_list:
            # 单细胞混合模式|Single-cell pooled mode
            args.extend(["-pool", self.config.cell_list])
        else:
            # 常规模式|Normal mode
            if self.config.long_reads:
                args.extend(["-long", self.config.long_reads])

                # PacBio数据|PacBio data
                if self.config.is_pacbio:
                    args.append("-lrpb")

            # 短reads|Short reads
            if self.config.left_reads and self.config.right_reads:
                args.extend(["-left", self.config.left_reads])
                args.extend(["-right", self.config.right_reads])

            if self.config.single_end_forward:
                args.extend(["-sef", self.config.single_end_forward])

            if self.config.single_end_reverse:
                args.extend(["-ser", self.config.single_end_reverse])

            # 链特异性|Strand-specific
            if self.config.stranded:
                args.append("-stranded")

            # 反向互补|Reverse-complement
            if self.config.revcomp_left:
                args.append("-revcomp-left")

            if self.config.revcomp_right:
                args.append("-revcomp-right")

        # 参考转录本|Reference transcripts
        if self.config.reference_transcripts:
            args.extend(["-ref", self.config.reference_transcripts])

        # Bloom filter配置|Bloom filter configuration
        if self.config.memory_gb:
            args.extend(["-mem", str(self.config.memory_gb)])

        if self.config.false_positive_rate is not None:
            args.extend(["-fpr", str(self.config.false_positive_rate)])

        if self.config.num_kmers is not None:
            args.extend(["-nk", str(self.config.num_kmers)])

        # 其他参数|Other parameters
        args.extend(["-t", str(self.config.threads)])
        args.extend(["-outdir", self.config.output_dir])

        # 停止阶段|Stop at stage
        if self.config.stage is not None:
            args.extend(["-stage", str(self.config.stage)])

        # 输出选项|Output options
        if self.config.write_uracil:
            args.append("-uracil")

        return args

    def run_assembly(self) -> bool:
        """运行转录组组装|Run transcriptome assembly

        Returns:
            bool: 成功返回True，失败返回False|True if successful, False otherwise
        """
        self.logger.info("=" * 60)
        self.logger.info("开始RNA-Bloom转录组组装|Starting RNA-Bloom transcriptome assembly")
        self.logger.info("=" * 60)
        self.logger.info(f"组装模式|Assembly mode: {self.config.get_assembly_mode()}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        # 构建命令参数|Build command arguments
        args = self.build_command_args()

        # 执行RNA-Bloom|Execute RNA-Bloom
        if self.rnabloom_jar == "rnabloom":
            # 使用conda安装的rnabloom命令|Use conda-installed rnabloom command
            cmd = "rnabloom " + " ".join(args)
            success = self.cmd_runner.run(cmd, "转录组组装|Transcriptome assembly")
        else:
            # 使用JAR文件|Use JAR file
            success = self.cmd_runner.run_java_jar(
                self.rnabloom_jar,
                args,
                "转录组组装|Transcriptome assembly"
            )

        if success:
            self.logger.info("=" * 60)
            self.logger.info("转录组组装完成|Transcriptome assembly completed")
            self.logger.info("=" * 60)
            self._post_assembly_summary()
        else:
            self.logger.error("=" * 60)
            self.logger.error("转录组组装失败|Transcriptome assembly failed")
            self.logger.error("=" * 60)

        return success

    def _post_assembly_summary(self):
        """组装后摘要|Post-assembly summary"""
        self.logger.info("组装结果摘要|Assembly result summary:")

        # 检查输出文件|Check output files
        output_files = [
            ("rnabloom.transcripts.fa", "主要转录本|Main transcripts"),
            ("rnabloom.transcripts.short.fa", "短转录本|Short transcripts"),
            ("rnabloom.transcripts.nr.fa", "去冗余转录本|Non-redundant transcripts")
        ]

        for file_name, description in output_files:
            file_path = os.path.join(self.config.output_dir, file_name)
            if os.path.exists(file_path):
                # 统计序列数量|Count sequences
                try:
                    result = subprocess.run(
                        ["grep", "-c", "^>", file_path],
                        capture_output=True,
                        text=True,
                        check=False
                    )
                    count = result.stdout.strip() if result.returncode == 0 else "未知|unknown"
                    self.logger.info(f"  {description}|{description}: {count} 条|sequences")
                except Exception:
                    self.logger.info(f"  {description}|{description}: 已生成|generated")
            else:
                self.logger.warning(f"  {description}不存在|{description} not found: {file_path}")

    def get_output_files(self) -> dict:
        """获取输出文件列表|Get list of output files

        Returns:
            dict: 输出文件路径字典|Dictionary of output file paths
        """
        output_files = {
            "transcripts": os.path.join(self.config.output_dir, "rnabloom.transcripts.fa"),
            "transcripts_short": os.path.join(self.config.output_dir, "rnabloom.transcripts.short.fa"),
            "transcripts_nr": os.path.join(self.config.output_dir, "rnabloom.transcripts.nr.fa"),
        }

        # 只返回存在的文件|Only return existing files
        return {k: v for k, v in output_files.items() if os.path.exists(v)}

    def validate_output(self) -> bool:
        """验证输出结果|Validate output results

        Returns:
            bool: 验证通过返回True|True if validation passed
        """
        self.logger.info("验证输出结果|Validating output results")

        output_files = self.get_output_files()

        if not output_files:
            self.logger.error("未找到任何输出文件|No output files found")
            return False

        # 检查主要转录本文件|Check main transcripts file
        if "transcripts" not in output_files:
            self.logger.error("主要转录本文件不存在|Main transcripts file not found")
            return False

        # 统计转录本数量|Count transcripts
        try:
            result = subprocess.run(
                ["grep", "-c", "^>", output_files["transcripts"]],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                count = int(result.stdout.strip())
                self.logger.info(f"组装了{count}条转录本|Assembled {count} transcripts")

                if count == 0:
                    self.logger.warning("警告：未组装出任何转录本|Warning: No transcripts assembled")
                    return False
        except Exception as e:
            self.logger.warning(f"无法统计转录本数量|Cannot count transcripts: {e}")

        self.logger.info("输出验证通过|Output validation passed")
        return True
