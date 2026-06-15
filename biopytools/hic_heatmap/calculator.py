"""Hi-C热图分析流程|Hi-C heatmap analysis pipeline (Juicer + PlotHiC)"""

import os
import subprocess
from pathlib import Path
from typing import List
import shutil

from .utils import build_conda_command


class HiCPipeline:
    """Hi-C热图分析流程类|Hi-C heatmap analysis pipeline class (Juicer + PlotHiC)"""

    def __init__(self, config, logger):
        """初始化流程|Initialize pipeline

        Args:
            config: HiCConfig配置对象|HiCConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 获取运行模式|Get run mode
        self.mode = self.config.get_mode()

        # 输出文件路径|Output file paths
        # 使用genome_id作为前缀|Use genome_id as prefix
        if self.config.genome_id:
            self.base_name = self.config.genome_id
        else:
            self.base_name = Path(self.config.genome).stem

        # 根据模式设置.hic文件路径|Set .hic file path based on mode
        if self.mode == 'hic':
            self.hic_file = Path(self.config.hic_file)
        else:
            self.hic_file = self.config.juicer_aligned_dir / f"inter_30.hic"

        self.heatmap_file = self.config.plot_dir / f"{self.base_name}_hic_heatmap.{self.config.output_format}"

    def run(self) -> bool:
        """运行完整流程|Run complete pipeline

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        mode_name = {'fastq': 'FASTQ', 'bam': 'BAM', 'hic': 'HiC'}[self.mode]
        self.logger.info(f"开始Hi-C热图分析流程 (Juicer + PlotHiC) - 模式|Mode: {mode_name}|Starting Hi-C heatmap analysis pipeline (Juicer + PlotHiC) - Mode: {mode_name}")
        self.logger.info("=" * 60)

        # 根据模式执行不同流程|Execute different workflows based on mode
        if self.mode == 'fastq':
            # 模式1：从FASTQ开始的完整流程|Mode 1: Complete workflow from FASTQ
            if not self._prepare_juicer_inputs():
                return False
            if not self._run_juicer():
                return False

        elif self.mode == 'bam':
            # 模式2：从BAM开始，跳过比对|Mode 2: Start from BAM, skip alignment
            self.logger.info("模式|Mode: BAM模式（跳过比对）|BAM mode (skip alignment)")
            if not self._prepare_juicer_inputs():
                return False
            if not self._run_juicer_from_bam():
                return False

        elif self.mode == 'hic':
            # 模式3：已有.hic文件，直接绘图|Mode 3: Existing .hic file, plot directly
            self.logger.info("模式|Mode: HiC模式（直接绘图）|HiC mode (direct plotting)")
            self.logger.info(f"使用现有Hi-C文件|Using existing Hi-C file: {self.hic_file}")

        # 所有模式都需要绘图|All modes need plotting
        if not self._plot_heatmap():
            return False

        self.logger.info("=" * 60)
        self.logger.info("Hi-C热图分析完成|Hi-C heatmap analysis completed")
        self.logger.info(f"结果保存在|Results saved to: {self.config.output_dir}")
        self.logger.info("=" * 60)

        return True

    def _prepare_juicer_inputs(self) -> bool:
        """准备Juicer输入文件|Prepare Juicer input files

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤1: 准备Juicer输入文件|Step 1: Preparing Juicer input files")

        # 创建references和restriction_sites目录|Create references and restriction_sites directories
        refs_dir = self.config.output_path / "references"
        sites_dir = self.config.output_path / "restriction_sites"
        refs_dir.mkdir(exist_ok=True)
        sites_dir.mkdir(exist_ok=True)

        # 1. 生成染色体大小文件|Generate chromosome sizes file
        chrom_sizes_file = refs_dir / f"{self.config.genome_id}.chrom.sizes"
        if not chrom_sizes_file.exists() or self.config.force_juicer:
            self.logger.info(f"生成染色体大小文件|Generating chromosome sizes: {chrom_sizes_file}")
            if not self._generate_chrom_sizes(chrom_sizes_file):
                return False
        else:
            self.logger.info(f"染色体大小文件已存在|Chromosome sizes file exists: {chrom_sizes_file}")

        self.config.chrom_sizes_file = chrom_sizes_file

        # 2. 生成酶切位点文件|Generate restriction sites file
        site_file = sites_dir / f"{self.config.genome_id}_{self.config.restriction_enzyme}.txt"
        if not site_file.exists() or self.config.force_juicer:
            self.logger.info(f"生成酶切位点文件|Generating restriction sites: {site_file}")
            if not self._generate_restriction_sites(site_file):
                return False
        else:
            self.logger.info(f"酶切位点文件已存在|Restriction sites file exists: {site_file}")

        self.config.site_file = site_file

        # 3. 检查BWA索引|Check BWA index
        if not os.path.exists(f"{self.config.genome}.bwt"):
            self.logger.info("构建BWA索引|Building BWA index")
            # 使用conda run wrapper|Use conda run wrapper
            cmd = build_conda_command(
                self.config.bwa_path,
                ["index", self.config.genome]
            )
            self._log_command(cmd, "BWA index")
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                self.logger.info("BWA索引已建立|BWA index built")
            except subprocess.CalledProcessError as e:
                self.logger.error(f"BWA索引失败|Failed to build BWA index: {e}")
                return False
        else:
            self.logger.info("BWA索引已存在|BWA index exists")

        return True

    def _generate_chrom_sizes(self, output_file: Path) -> bool:
        """生成染色体大小文件|Generate chromosome sizes file

        Args:
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            # 使用samtools faidx生成索引|Generate index using samtools faidx
            # 使用conda run wrapper|Use conda run wrapper
            cmd = build_conda_command(
                self.config.samtools_path,
                ["faidx", self.config.genome]
            )
            self._log_command(cmd, "samtools faidx")
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # 从.fai文件生成chrom.sizes|Generate chrom.sizes from .fai file
            fai_file = Path(f"{self.config.genome}.fai")
            with open(fai_file, 'r') as f_in, open(output_file, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        f_out.write(f"{fields[0]}\t{fields[1]}\n")

            self.logger.info(f"染色体大小文件已生成|Chromosome sizes file created: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成染色体大小文件失败|Failed to generate chromosome sizes: {e}")
            return False

    def _generate_restriction_sites(self, output_file: Path) -> bool:
        """生成酶切位点文件|Generate restriction sites file

        Args:
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 是否成功|Whether successful
        """
        # 检查是否已存在|Check if already exists
        if output_file.exists() and not self.config.force_juicer:
            self.logger.info(f"酶切位点文件已存在|Restriction sites file exists: {output_file}")
            self.config.site_file = output_file
            return True

        try:
            # 使用Juicer的generate_site_positions.py脚本
            script_path = self.config.juicer_dir / "misc" / "generate_site_positions.py"

            # 检查脚本是否存在|Check if script exists
            if not script_path.exists():
                self.logger.warning(f"Juicer脚本不存在|Juicer script not found: {script_path}")
                self.logger.warning("将尝试使用简化的酶切位点生成方式|Will try simplified restriction site generation")
                return self._generate_restriction_sites_simple(output_file)

            # 检查Python2是否可用|Check if Python2 is available
            python2_cmd = self._find_python2()
            if not python2_cmd:
                self.logger.warning("未找到Python2|Python2 not found")
                self.logger.warning("将尝试使用Python3运行脚本|Will try to run script with Python3")
                python2_cmd = "python3"

            # 临时文件|Temporary file
            temp_file = self.config.output_path / f"{self.config.genome_id}_{self.config.restriction_enzyme}.txt"

            cmd = [
                python2_cmd,
                str(script_path),
                self.config.restriction_enzyme,
                self.config.genome_id,
                self.config.genome
            ]

            self._log_command(cmd, "generate_site_positions.py")
            self.logger.info("生成酶切位点中|Generating restriction sites...")

            result = subprocess.run(
                cmd,
                cwd=self.config.output_path,  # 在输出目录运行|Run in output directory
                capture_output=True,
                text=True,
                timeout=600  # 10分钟超时|10 minutes timeout
            )

            # 如果Python2失败，尝试Python3
            if result.returncode != 0 and python2_cmd == "python2":
                self.logger.warning("Python2运行失败，尝试Python3|Python2 failed, trying Python3")
                cmd[0] = "python3"
                result = subprocess.run(
                    cmd,
                    cwd=self.config.output_path,
                    capture_output=True,
                    text=True,
                    timeout=600
                )

            if result.returncode != 0:
                self.logger.warning(f"Juicer脚本失败|Juicer script failed: {result.stderr}")
                self.logger.info("尝试使用简化方式生成酶切位点|Trying simplified restriction site generation")
                return self._generate_restriction_sites_simple(output_file)

            # 排序并移动到最终位置|Sort and move to final location
            self.logger.info("排序酶切位点文件|Sorting restriction sites file...")
            cmd_sort = ["sort", "-k1,1V", "-k2,2n", str(temp_file)]
            with open(output_file, 'w') as f_out:
                subprocess.run(cmd_sort, stdout=f_out, check=True)

            # 删除临时文件|Remove temporary file
            temp_file.unlink()

            self.logger.info(f"酶切位点文件已生成|Restriction sites file created: {output_file}")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error("生成酶切位点超时|Timeout generating restriction sites")
            return False
        except Exception as e:
            self.logger.error(f"生成酶切位点失败|Failed to generate restriction sites: {e}")
            return False

    def _generate_restriction_sites_simple(self, output_file: Path) -> bool:
        """使用简化方式生成酶切位点文件|Generate restriction sites file using simplified method

        Args:
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq

            self.logger.info("使用简化方式生成酶切位点|Using simplified method to generate restriction sites")

            # 简单的酶切位点识别|Simple restriction site recognition
            enzyme_sites = {
                'MboI': 'GATC',
                'HindIII': 'AAGCTT',
                'NcoI': 'CCATGG',
                'EcoRI': 'GAATTC',
                'BamHI': 'GGATCC'
            }

            if self.config.restriction_enzyme not in enzyme_sites:
                self.logger.error(f"不支持的酶|Unsupported enzyme: {self.config.restriction_enzyme}")
                self.logger.error(f"支持的酶|Supported enzymes: {', '.join(enzyme_sites.keys())}")
                return False

            site_sequence = enzyme_sites[self.config.restriction_enzyme]

            self.logger.info(f"扫描基因组中的酶切位点|Scanning genome for {self.config.restriction_enzyme} sites ({site_sequence})...")

            with open(output_file, 'w') as f_out:
                for chrom in SeqIO.parse(self.config.genome, "fasta"):
                    chrom_name = chrom.id
                    seq = str(chrom.seq)

                    # 查找所有酶切位点|Find all restriction sites
                    site_positions = []
                    pos = 0
                    while True:
                        pos = seq.find(site_sequence, pos)
                        if pos == -1:
                            break
                        # 1-based position for Juicer
                        site_positions.append(pos + 1)
                        pos += 1

                    # 写入位点|Write sites
                    for site_pos in site_positions:
                        f_out.write(f"{chrom_name}\t{site_pos}\t{site_pos}\t{self.config.restriction_enzyme}\n")

                    self.logger.debug(f"{chrom_name}: 找到{len(site_positions)}个位点|Found {len(site_positions)} sites")

            self.logger.info(f"酶切位点文件已生成|Restriction sites file created: {output_file}")
            return True

        except ImportError:
            self.logger.error("需要Biopython包|Biopython package required")
            self.logger.error("请安装: pip install biopython|Please install: pip install biopython")
            return False
        except Exception as e:
            self.logger.error(f"简化方式生成失败|Simplified generation failed: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _find_python2(self) -> str:
        """查找Python2可执行文件|Find Python2 executable

        Returns:
            str: Python2路径或空字符串|Python2 path or empty string
        """
        # 常见Python2路径|Common Python2 paths
        python2_candidates = [
            "python2",
            "python2.7",
            "~/miniforge3/envs/juicer_v.1.6/bin/python2",
        ]

        for python_cmd in python2_candidates:
            try:
                result = subprocess.run(
                    [python_cmd, "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                # 检查是否是Python 2.x|Check if it's Python 2.x
                if result.returncode == 0 and ("2.7" in result.stdout or "2." in result.stdout):
                    return python_cmd
            except:
                continue

        return ""

    def _find_plothic(self) -> str:
        """查找plothic可执行文件|Find plothic executable

        Returns:
            str: plothic路径或空字符串|plothic path or empty string
        """
        # 首先检查配置中指定的plothic路径|First check configured plothic path
        if self.config.plothic_path and os.path.exists(self.config.plothic_path):
            try:
                result = subprocess.run(
                    [self.config.plothic_path, "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    return str(self.config.plothic_path)
            except:
                pass

        # 常见plothic路径|Common plothic paths
        plothic_candidates = [
            "plothic",
        ]

        for plothic_cmd in plothic_candidates:
            try:
                result = subprocess.run(
                    [plothic_cmd, "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    return plothic_cmd
            except:
                continue

        return ""

    def _run_juicer(self) -> bool:
        """运行Juicer流程|Run Juicer pipeline

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤2: 运行Juicer|Step 2: Running Juicer")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.hic_file.exists():
            self.logger.info(f"跳过|Skipping: Hi-C文件已存在|Hi-C file exists: {self.hic_file}")
            return True

        # 清理旧的输出|Clean old outputs
        # Juicer要求aligned目录不存在才能运行|Juicer requires aligned directory to not exist
        self.logger.info("清理旧的Juicer输出|Cleaning old Juicer outputs...")
        for dir_name in ["aligned", "splits", "duplicates"]:
            dir_path = self.config.output_path / dir_name
            if dir_path.exists():
                shutil.rmtree(dir_path)
                self.logger.info(f"删除目录|Deleted: {dir_path}")

        # 准备FASTQ文件|Prepare FASTQ files
        # Juicer期望FASTQ文件在output/fastq/目录中，文件名格式为*_R*.fastq*|Juicer expects FASTQ files in output/fastq/ with naming *_R*.fastq*
        fastq_dir = self.config.output_path / "fastq"
        fastq_dir.mkdir(exist_ok=True)

        # 检查FASTQ文件|Check FASTQ files
        if not self.config.fastq_r1 or not self.config.fastq_r2:
            self.logger.error("FASTQ模式需要fastq_r1和fastq_r2|FASTQ mode requires fastq_r1 and fastq_r2")
            return False

        # 复制或链接FASTQ文件到fastq目录，并重命名为符合Juicer要求的格式|Copy or link FASTQ files to fastq directory with Juicer-compatible naming
        # Juicer要求文件名包含_R1和_R2|Juicer requires _R1 and _R2 in filenames
        import os
        import re

        # 提取文件名的基本部分（去掉_R1, _R2, _1, _2等后缀）|Extract base name (remove _R1, _R2, _1, _2 suffixes)
        def get_fastq_base(filename):
            """提取FASTQ文件的基本名称|Extract FASTQ base name"""
            name = Path(filename).stem  # 去掉.gz等扩展名|Remove .gz extension
            # 去掉常见的read标识符|Remove common read identifiers
            name = re.sub(r'[_\.][Rr]?[12]$', '', name)  # _R1, _R2, _1, _2, .r1, .r2
            name = re.sub(r'[_\.]read[12]$', '', name, flags=re.IGNORECASE)  # _read1, _read2
            return name

        # 获取共同的基本名称|Get common base name
        base_name = get_fastq_base(self.config.fastq_r1)

        # 构建符合Juicer要求的文件名|Build Juicer-compatible filenames
        # 格式: basename_R1.fastq.gz | Format: basename_R1.fastq.gz
        r1_target_name = f"{base_name}_R1.fastq.gz"
        r2_target_name = f"{base_name}_R2.fastq.gz"

        r1_target = fastq_dir / r1_target_name
        r2_target = fastq_dir / r2_target_name

        # 复制FASTQ文件并重命名|Copy and rename FASTQ files
        # 直接复制文件而不是创建软链接，避免symlink路径解析问题
        # Copy files directly instead of creating symlinks to avoid path resolution issues
        if not r1_target.exists() or not os.path.samefile(self.config.fastq_r1, r1_target):
            if r1_target.exists():
                r1_target.unlink()
            shutil.copy2(self.config.fastq_r1, r1_target)
            self.logger.info(f"复制|Copied: {self.config.fastq_r1} -> {r1_target}")

        if not r2_target.exists() or not os.path.samefile(self.config.fastq_r2, r2_target):
            if r2_target.exists():
                r2_target.unlink()
            shutil.copy2(self.config.fastq_r2, r2_target)
            self.logger.info(f"复制|Copied: {self.config.fastq_r2} -> {r2_target}")

        # 构建Juicer命令|Build Juicer command
        # 将所有路径转换为绝对路径，因为会切换到output目录运行|Convert all paths to absolute since we'll cd to output directory
        genome_abs = os.path.abspath(self.config.genome)
        chrom_sizes_abs = os.path.abspath(self.config.chrom_sizes_file)
        site_abs = os.path.abspath(self.config.site_file)
        output_abs = os.path.abspath(self.config.output_path)
        juicer_dir_abs = os.path.abspath(self.config.juicer_dir)

        cmd = [
            str(self.config.juicer_sh),
            "-g", self.config.genome_id,
            "-d", output_abs,
            "-s", self.config.restriction_enzyme,
            "-p", chrom_sizes_abs,
            "-y", site_abs,
            "-z", genome_abs,
            "-D", juicer_dir_abs,
            "-t", str(self.config.threads)
        ]

        self._log_command(cmd, "Juicer")
        self.logger.info("运行Juicer流程中|Running Juicer pipeline...")
        self.logger.info("这可能需要较长时间，请耐心等待|This may take a while, please be patient...")

        try:
            # 设置环境变量|Set environment variables
            env = os.environ.copy()
            env["PATH"] = f"{os.path.dirname(self.config.java_path)}:{env['PATH']}"

            # 切换到输出目录运行Juicer|Change to output directory to run Juicer
            # 这样Juicer创建的symlinks路径才能正确解析
            # This ensures Juicer-created symlinks resolve correctly
            original_dir = os.getcwd()
            os.chdir(str(self.config.output_path))

            try:
                result = subprocess.run(
                    cmd,
                    env=env,
                    check=True,
                    timeout=72000,  # 2小时超时|2 hours timeout
                    text=True
                )
            finally:
                # 恢复原始工作目录|Restore original working directory
                os.chdir(original_dir)

            self.logger.info(f"Hi-C文件已生成|Hi-C file created: {self.hic_file}")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error("Juicer运行超时|Juicer timeout")
            return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Juicer运行失败|Juicer failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"Juicer异常|Juicer error: {e}")
            return False

    def _run_juicer_from_bam(self) -> bool:
        """从BAM文件运行Juicer后处理|Run Juicer post-processing from BAM file

        注意：Juicer官方脚本不支持直接从BAM开始
        此方法将BAM复制到aligned目录并使用Juicer的CPP模块生成.hic
        Note: Juicer official script doesn't support starting from BAM directly
        This method copies BAM to aligned directory and uses Juicer's CPP module to generate .hic

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤2: 从BAM运行Juicer后处理|Step 2: Running Juicer post-processing from BAM")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.hic_file.exists():
            self.logger.info(f"跳过|Skipping: Hi-C文件已存在|Hi-C file exists: {self.hic_file}")
            return True

        self.logger.warning("注意：Juicer不原生支持从BAM开始|Warning: Juicer doesn't natively support starting from BAM")
        self.logger.info("推荐方案1：使用.hic模式（如果已有.hic文件）|Recommended option 1: Use .hic mode (if .hic file exists)")
        self.logger.info("推荐方案2：使用HiC-Pro从BAM生成.hic|Recommended option 2: Use HiC-Pro to generate .hic from BAM")
        self.logger.info("当前方案：尝试从BAM提取paired-end reads并生成Juicer格式|Current approach: Try to extract paired-end reads from BAM and generate Juicer format")

        # TODO: 实现BAM到Juicer格式的转换
        # 这需要：
        # 1. 从BAM提取paired-end reads
        # 2. 转换为Juicer的merged_nodup格式
        # 3. 使用Juicer的CPP模块生成.hic

        self.logger.error("BAM模式暂未实现|BAM mode not yet implemented")
        self.logger.error("请使用以下替代方案|Please use the following alternatives:")
        self.logger.error("  1. 使用-h/--hic-file参数直接提供.hic文件|  1. Use -h/--hic-file to provide .hic file directly")
        self.logger.error("  2. 使用FASTQ模式运行完整Juicer流程|  2. Use FASTQ mode to run complete Juicer pipeline")
        self.logger.error("  3. 使用其他工具（如HiC-Pro）从BAM生成.hic|  3. Use other tools (e.g., HiC-Pro) to generate .hic from BAM")

        return False

    def _plot_heatmap(self) -> bool:
        """绘制Hi-C热图|Plot Hi-C heatmap using plothic command

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤3: 绘制全基因组热图|Step 3: Plotting whole genome heatmap")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.heatmap_file.exists() and self.heatmap_file.stat().st_size > 0:
            self.logger.info(f"跳过|Skipping: 热图已存在|Heatmap exists: {self.heatmap_file}")
            return True

        # 方法1：尝试使用plothic（需要3D-DNA处理过的.hic文件）|Method 1: Try plothic (requires 3D-DNA processed .hic)
        try:
            # 需要先生成染色体信息文件|Generate chromosome info file first
            chr_info_file = self.config.plot_dir / "chr_info.txt"
            self._generate_chr_info(chr_info_file)

            # 查找plothic可执行文件|Find plothic executable
            plothic_cmd = self._find_plothic()
            if not plothic_cmd:
                self.logger.error("未找到plothic命令|plothic command not found")
                return False

            # 构建plothic命令|Build plothic command
            # 使用conda run wrapper|Use conda run wrapper
            cmd = build_conda_command(
                plothic_cmd,
                [
                    "-hic", str(self.hic_file),
                    "-chr", str(chr_info_file),
                    "-o", str(self.config.plot_dir),
                    "-r", str(self.config.resolution),
                    "-cmap", self.config.color_map,
                    "-format", self.config.output_format,
                    "-dpi", str(self.config.dpi),
                    "-g", self.config.genome_id
                ]
            )

            self._log_command(cmd, "plothic")
            self.logger.info("绘制热图中|Plotting heatmap...")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800  # 30分钟超时|30 minutes timeout
            )

            if result.returncode != 0:
                self.logger.warning(f"plothic失败|plothic failed，尝试备用方法|trying fallback method")
                self.logger.debug(f"plothic stderr: {result.stderr}")
                self.logger.debug(f"plothic stdout: {result.stdout}")
                return self._plot_heatmap_python()

            # plothic输出的文件名格式可能是：genomeID_resolution.{format}
            # 查找生成的文件|Find generated file
            import glob
            pattern = str(self.config.plot_dir / f"{self.config.genome_id}_*.{self.config.output_format}")
            generated_files = glob.glob(pattern)

            if not generated_files:
                # 也可能只是heatmap.{format}
                pattern2 = str(self.config.plot_dir / f"*.{self.config.output_format}")
                generated_files = glob.glob(pattern2)

            if generated_files:
                # 重命名到目标文件|Rename to target file
                generated_file = Path(generated_files[0])
                if generated_file != self.heatmap_file:
                    import shutil
                    shutil.move(str(generated_file), str(self.heatmap_file))
                    self.logger.info(f"重命名|Renamed: {generated_file} -> {self.heatmap_file}")
                self.logger.info(f"热图已生成|Heatmap created: {self.heatmap_file}")
                return True
            else:
                # plothic成功但没有生成文件，触发备用方案
                self.logger.warning("plothic成功但未生成输出文件|plothic succeeded but no output file，尝试备用方法|trying fallback method")
                return self._plot_heatmap_python()

        except subprocess.TimeoutExpired:
            self.logger.warning("plothic绘制超时|plothic timeout，尝试备用方法|trying fallback method")
            return self._plot_heatmap_python()
        except Exception as e:
            self.logger.warning(f"plothic绘制失败|plothic failed: {e}，尝试备用方法|trying fallback method")
            import traceback
            traceback.print_exc()
            return self._plot_heatmap_python()

    def _plot_heatmap_python(self):
        """使用Python原生绘图（备用方案）|Plot using Python native (fallback method)

        使用hicstraw + matplotlib直接读取.hic文件并绘制热图
        Use hicstraw + matplotlib to read .hic file and plot heatmap directly

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            import hicstraw
            import matplotlib
            matplotlib.use('Agg')  # 非交互式后端|Non-interactive backend
            import matplotlib.pyplot as plt
            import numpy as np

            self.logger.info("使用Python原生绘图|Using Python native plotting")

            # 读取.hic文件|Read .hic file
            hic_obj = hicstraw.HiCFile(str(self.hic_file))

            # 获取可用分辨率|Get available resolutions
            resolutions = hic_obj.getResolutions()
            self.logger.info(f"可用分辨率|Available resolutions: {resolutions}")

            # 选择最接近的分辨率|Select closest resolution
            target_res = self.config.resolution
            if target_res not in resolutions:
                # 找最接近的|Find closest
                closest = min(resolutions, key=lambda x: abs(x - target_res))
                self.logger.info(f"使用最接近的分辨率|Using closest resolution: {closest}")
                target_res = closest

            # 获取染色体信息|Get chromosome info
            chromosomes = [chrom for chrom in hic_obj.getChromosomes()
                         if chrom.name not in ['All', 'all', 'assembly']]

            if not chromosomes:
                self.logger.error("未找到有效染色体|No valid chromosomes found")
                return False

            # 只绘制前几个染色体（避免内存问题）|Only plot first few chromosomes
            max_chroms = 8
            chromosomes = chromosomes[:max_chroms]

            self.logger.info(f"绘制染色体|Plotting chromosomes: {[c.name for c in chromosomes]}")

            # 计算累积位置|Calculate cumulative positions
            chr_starts = {}
            cumulative = 0
            for chrom in chromosomes:
                chr_starts[chrom.name] = cumulative
                cumulative += chrom.length

            total_length = cumulative

            # 创建矩阵|Create matrix
            matrix_size = total_length // target_res + 1
            contact_matrix = np.zeros((matrix_size, matrix_size))

            # 逐对染色体提取数据|Extract data for each chromosome pair
            for i, chrom1 in enumerate(chromosomes):
                for j, chrom2 in enumerate(chromosomes):
                    if i > j:  # 只处理上三角|Only process upper triangle
                        continue

                    try:
                        matrix_obj = hic_obj.getMatrixZoomData(
                            chrom1.name, chrom2.name,
                            "observed", "NONE", "BP", target_res
                        )

                        # 获取数据|Get data
                        records = matrix_obj.getRecordsAsMatrix(
                            0, chrom1.length,
                            0, chrom2.length
                        )

                        # 处理NaN值：将NaN替换为0|Handle NaN values: replace NaN with 0
                        records = np.nan_to_num(records, nan=0.0, posinf=0.0, neginf=0.0)

                        # 计算在全局矩阵中的位置|Calculate position in global matrix
                        start_row = chr_starts[chrom1.name] // target_res
                        start_col = chr_starts[chrom2.name] // target_res
                        end_row = start_row + records.shape[0]
                        end_col = start_col + records.shape[1]

                        # 确保不越界|Ensure no out of bounds
                        end_row = min(end_row, matrix_size)
                        end_col = min(end_col, matrix_size)

                        # 填充矩阵|Fill matrix
                        contact_matrix[start_row:end_row, start_col:end_col] = records[:end_row-start_row, :end_col-start_col]

                        # 对称填充|Symmetric fill
                        if i != j:
                            contact_matrix[start_col:end_col, start_row:end_row] = records.T[:end_col-start_col, :end_row-start_row]

                    except Exception as e:
                        self.logger.warning(f"无法提取{chrom1.name}-{chrom2.name}|Failed to extract {chrom1.name}-{chrom2.name}: {e}")
                        continue

            # 绘制热图|Plot heatmap
            fig, ax = plt.subplots(figsize=(12, 12))

            # 对数变换|Log transform
            # 添加一个小值避免log(0)|Add small value to avoid log(0)
            contact_matrix = np.log10(contact_matrix + 1)

            # 处理对数变换后的NaN/inf（虽然上面已经处理过，但再保险一次）
            # Handle NaN/inf after log transform (already handled above, but double-check)
            contact_matrix = np.nan_to_num(contact_matrix, nan=0.0, posinf=0.0, neginf=0.0)

            # 使用vmin和vmax设置合理的颜色范围|Use vmin and vmax to set reasonable color range
            # 排除对角线附近的值来计算颜色范围|Exclude values near diagonal to calculate color range
            # 这样可以更好地展示染色体间的相互作用|This better shows inter-chromosomal interactions
            mask = np.triu_indices_from(contact_matrix, k=10)  # 排除对角线附近10个bin|Exclude 10 bins near diagonal
            valid_values = contact_matrix[mask]
            valid_values = valid_values[valid_values > 0]  # 只考虑非零值|Only consider non-zero values

            if len(valid_values) > 0:
                vmax = np.percentile(valid_values, 99)  # 使用99百分位数作为最大值|Use 99th percentile as max
                vmin = np.percentile(valid_values, 1)   # 使用1百分位数作为最小值|Use 1st percentile as min
            else:
                vmax = np.max(contact_matrix)
                vmin = np.min(contact_matrix[contact_matrix > 0]) if np.any(contact_matrix > 0) else 0

            im = ax.imshow(contact_matrix, cmap=self.config.color_map, aspect='equal', vmin=vmin, vmax=vmax)

            # 添加染色体边界线|Add chromosome boundary lines
            cumulative = 0
            for chrom in chromosomes[:-1]:
                cumulative += chrom.length
                pos = cumulative // target_res
                ax.axhline(pos, color='white', linewidth=0.5)
                ax.axvline(pos, color='white', linewidth=0.5)

            # 添加颜色条|Add colorbar
            plt.colorbar(im, ax=ax, label='log10(contacts + 1)')

            # 设置标题|Set title
            ax.set_title(f"{self.config.genome_id} Hi-C Heatmap ({target_res/1000:.0f}kb)")
            ax.set_xlabel("Genome position")
            ax.set_ylabel("Genome position")

            plt.tight_layout()
            plt.savefig(self.heatmap_file, dpi=self.config.dpi, format=self.config.output_format)
            plt.close()

            self.logger.info(f"热图已生成|Heatmap created: {self.heatmap_file}")
            return True

        except ImportError:
            self.logger.error("缺少依赖包|Missing dependencies: hicstraw, matplotlib, numpy")
            return False
        except Exception as e:
            self.logger.error(f"Python绘图失败|Python plotting failed: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _generate_chr_info(self, output_file: Path):
        """生成PlotHiC所需的染色体信息文件|Generate chromosome info file for PlotHiC

        Args:
            output_file: 输出文件路径|Output file path
        """
        # 读取染色体大小|Read chromosome sizes
        chrom_sizes = {}
        with open(self.config.chrom_sizes_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    chrom_sizes[fields[0]] = int(fields[1])

        # 写入chr_info.txt格式|Write chr_info.txt format
        # 格式: name cumulative_position index
        # 注意：第二列是染色体在.hic文件中的累积位置（cumulative position），不是真实长度
        # Note: Column 2 is cumulative position in .hic file, not actual length
        with open(output_file, 'w') as f:
            f.write("# name\tlength\tindex\n")
            cumulative_pos = 0
            for i, (chrom, length) in enumerate(sorted(chrom_sizes.items()), 1):
                cumulative_pos += length
                f.write(f"{chrom}\t{cumulative_pos}\t{i}\n")

        self.logger.info(f"染色体信息文件已生成|Chromosome info file created: {output_file}")

    def _log_command(self, cmd_list, description=""):
        """记录执行的命令|Log executed command

        Args:
            cmd_list: 命令列表|Command list
            description: 命令描述|Command description
        """
        if isinstance(cmd_list, list):
            cmd_str = ' '.join(
                f'"{arg}"' if ' ' in str(arg) else str(arg)
                for arg in cmd_list
            )
        else:
            cmd_str = str(cmd_list)

        if description:
            self.logger.info(f"{description}命令|{description} command:")
        self.logger.info(f"  {cmd_str}")
