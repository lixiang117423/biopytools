"""
SignalP 6.0信号肽预测工具函数模块|SignalP 6.0 Signal Peptide Prediction Utility Functions Module
"""

import logging
import os
import subprocess
import sys
import shutil
import re
from pathlib import Path
from typing import Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先尝试从命令路径检测|First try to detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command_string(command: str, args: str) -> str:
    """
    构建conda run命令字符串（用于需要shell特性的命令）|Build conda run command string (for commands needing shell features)

    Args:
        command: 命令名称|Command name
        args: 命令参数字符串|Command arguments string

    Returns:
        完整命令字符串|Complete command string
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return f"conda run -n {conda_env} {command} {args}"
    else:
        return f"{command} {args}"


class SignalPLogger:
    """SignalP信号肽预测日志管理器|SignalP Signal Peptide Prediction Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "signalp_prediction.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除旧日志文件|Delete old log file
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class SignalPRunner:
    """SignalP命令执行器|SignalP Command Runner"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def build_command(self):
        """构建SignalP命令|Build SignalP command"""
        cmd_parts = [
            self.config.signalp_path,
            "--fastafile", self.config.fasta_file,
            "--output_dir", self.config.output_dir,
            "--format", self.config.format,
            "--organism", self.config.organism,
            "--mode", self.config.mode,
            "--bsize", str(self.config.bsize),
            "--write_procs", str(self.config.write_procs),
            "--torch_num_threads", str(self.config.torch_num_threads)
        ]

        # 可选参数|Optional parameters
        if self.config.skip_resolve:
            cmd_parts.append("--skip_resolve")

        if self.config.model_dir:
            cmd_parts.extend(["--model_dir", self.config.model_dir])

        return cmd_parts

    def run(self):
        """运行SignalP预测（自动检测conda环境）|Run SignalP prediction (auto-detect conda environment)"""
        self.logger.info("=" * 80)
        self.logger.info("开始SignalP 6.0信号肽预测|Starting SignalP 6.0 Signal Peptide Prediction")
        self.logger.info("=" * 80)

        # 显示配置信息|Show configuration
        self.logger.info(f"输入文件|Input file: {self.config.fasta_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"生物类型|Organism: {self.config.organism}")
        self.logger.info(f"预测模式|Mode: {self.config.mode}")
        self.logger.info(f"输出格式|Format: {self.config.format}")
        self.logger.info(f"批处理大小|Batch size: {self.config.bsize}")
        self.logger.info(f"写入进程数|Write processes: {self.config.write_procs}")
        self.logger.info(f"PyTorch线程数|PyTorch threads: {self.config.torch_num_threads}")

        if self.config.model_dir:
            self.logger.info(f"模型目录|Model directory: {self.config.model_dir}")

        self.logger.info("=" * 80)

        # 构建命令|Build command
        cmd_parts = self.build_command()
        cmd = " ".join(cmd_parts)

        # 自动包装conda环境的命令|Auto-wrap conda environment commands
        signalp_cmd = os.path.basename(self.config.signalp_path)
        cmd_args = " ".join(cmd_parts[1:])  # 跳过命令本身，只保留参数
        wrapped_cmd = build_conda_command_string(signalp_cmd, cmd_args)

        self.logger.info(f"执行命令|Executing command: {wrapped_cmd}")

        try:
            # 执行命令|Execute command
            result = subprocess.run(
                wrapped_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )

            # 输出结果|Output results
            if result.stdout:
                self.logger.info(result.stdout)

            self.logger.info("=" * 80)
            self.logger.info("SignalP预测完成|SignalP prediction completed")
            self.logger.info("=" * 80)

            # 检查输出文件|Check output files
            self._check_output_files()

            # 格式化结果|Format results
            self._format_results()

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"SignalP执行失败|SignalP execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False

        except Exception as e:
            self.logger.error(f"预测出错|Prediction error: {e}")
            return False

    def _check_output_files(self):
        """检查输出文件|Check output files"""
        self.logger.info("检查输出文件|Checking output files:")

        expected_files = [
            "prediction_results.txt",
            "processed_entries.fasta",
            "output.gff3",
            "region_output.gff3",
            "output.json"
        ]

        for filename in expected_files:
            filepath = os.path.join(self.config.output_dir, filename)
            if os.path.exists(filepath):
                size = os.path.getsize(filepath)
                self.logger.info(f"  ✓ {filename} ({size} bytes)")
            else:
                self.logger.warning(f"  ✗ {filename} (未生成|not generated)")

    def _format_results(self):
        """格式化预测结果|Format prediction results"""
        prediction_results = os.path.join(self.config.output_dir, "prediction_results.txt")

        if not os.path.exists(prediction_results):
            self.logger.warning("未找到prediction_results.txt，跳过格式化|prediction_results.txt not found, skipping formatting")
            return

        # 创建易读格式文件|Create readable format file
        readable_output = os.path.join(self.config.output_dir, "prediction_results_readable.txt")

        formatter = SignalPResultFormatter(
            input_file=prediction_results,
            output_file=readable_output,
            logger=self.logger
        )

        # 生成易读格式（带注释的中文版本）|Generate readable format (with comments, Chinese version)
        success = formatter.format()
        if success:
            self.logger.info(f"易读格式结果已保存|Readable format saved: {readable_output}")
        else:
            self.logger.error("格式化失败|Formatting failed")

        # 生成R友好格式（纯TSV，英文列名）|Generate R-friendly format (pure TSV, English columns)
        success_r = formatter.format_for_r()
        if success_r:
            self.logger.info(f"R友好格式已生成|R-friendly format generated: signalp_summary.tsv")
        else:
            self.logger.warning("R格式生成失败|R format generation failed")

        # 清理plot文件|Cleanup plot files
        if self.config.cleanup_plots:
            self._cleanup_plot_files()

    def _cleanup_plot_files(self):
        """清理plot文件|Cleanup plot files

        删除output_*_plot.txt文件以节省空间
        Delete output_*_plot.txt files to save space
        """
        import glob

        plot_files = glob.glob(os.path.join(self.config.output_dir, "output_*_plot.txt"))

        if not plot_files:
            self.logger.info("没有找到plot文件|No plot files found")
            return

        self.logger.info(f"开始清理plot文件|Cleaning up plot files: {len(plot_files)} files")

        # 计算总大小|Calculate total size
        total_size = sum(os.path.getsize(f) for f in plot_files)

        # 删除文件|Delete files
        deleted_count = 0
        for plot_file in plot_files:
            try:
                os.remove(plot_file)
                deleted_count += 1
            except Exception as e:
                self.logger.warning(f"删除失败|Failed to delete {plot_file}: {e}")

        self.logger.info(f"已删除|Deleted {deleted_count}/{len(plot_files)} plot files")
        self.logger.info(f"释放空间|Space freed: {format_number(total_size)} bytes")


def format_number(num: int) -> str:
    """格式化数字|Format number

    大于1百万的数字使用M单位|Numbers larger than 1M use M unit
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


class SignalPResultFormatter:
    """SignalP结果格式化器|SignalP Result Formatter

    将原始prediction_results.txt转换为用户友好的格式
    Convert raw prediction_results.txt to user-friendly format
    """

    # 预测类型映射|Prediction type mapping
    PREDICTION_MAP = {
        'SP': {
            'zh': '分泌蛋白-经典信号肽',
            'en': 'Secreted protein - Classical signal peptide (Sec/SPI)',
            'has_sp': True
        },
        'OTHER': {
            'zh': '非信号肽',
            'en': 'Non-signal peptide',
            'has_sp': False
        },
        'SPI': {
            'zh': '分泌蛋白-经典信号肽',
            'en': 'Secreted protein - Classical signal peptide (Sec/SPI)',
            'has_sp': True
        },
        'SPII': {
            'zh': '脂蛋白-脂蛋白信号肽',
            'en': 'Lipoprotein - Lipoprotein signal peptide (Sec/SPII)',
            'has_sp': True
        },
        'SPIII': {
            'zh': '外膜蛋白-细菌信号肽',
            'en': 'Outer membrane protein - Bacterial signal peptide (Sec/SPIII)',
            'has_sp': True
        },
        'Tat': {
            'zh': '双精氨酸转运蛋白-Tat信号肽',
            'en': 'Tat protein - Twin-arginine signal peptide (Tat/SPI)',
            'has_sp': True
        }
    }

    def __init__(self, input_file: str, output_file: str, logger=None):
        """初始化格式化器|Initialize formatter

        Args:
            input_file: 输入文件路径|Input file path (prediction_results.txt)
            output_file: 输出文件路径|Output file path (prediction_results_readable.txt)
            logger: 日志器|Logger instance
        """
        self.input_file = input_file
        self.output_file = output_file
        self.logger = logger

        # 解析文件头信息|Parse file header information
        self.organism = "Unknown"
        self.timestamp = ""

    def parse_header(self, line: str):
        """解析文件头|Parse file header"""
        # SignalP 6.0 header格式: # SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20260204185141
        parts = line.split('\t')
        for part in parts:
            if 'Organism:' in part:
                self.organism = part.split('Organism:')[1].strip()
            elif 'Timestamp:' in part:
                timestamp_str = part.split('Timestamp:')[1].strip()
                # 格式化时间戳|Format timestamp: 20260204185141 -> 2026-02-04 18:51:41
                if len(timestamp_str) == 14:
                    try:
                        year = timestamp_str[:4]
                        month = timestamp_str[4:6]
                        day = timestamp_str[6:8]
                        hour = timestamp_str[8:10]
                        minute = timestamp_str[10:12]
                        second = timestamp_str[12:14]
                        self.timestamp = f"{year}-{month}-{day} {hour}:{minute}:{second}"
                    except Exception:
                        self.timestamp = timestamp_str

    def parse_cs_position(self, cs_str: str) -> tuple:
        """解析切割位点信息|Parse cleavage site information

        Args:
            cs_str: 切割位点字符串|Cleavage site string (e.g., "CS pos: 20-21. Pr: 0.9766")

        Returns:
            (sp_end, cs_start, cs_end, confidence)
            信号肽结束位置, 切割位点起始, 切割位点结束, 置信度
        """
        if not cs_str or cs_str.strip() == "":
            return ("-", "-", "-", "-")

        try:
            # 提取位置|Extract position
            if "CS pos:" in cs_str:
                pos_part = cs_str.split("CS pos:")[1].split(".")[0].strip()
                if "-" in pos_part:
                    parts = pos_part.split("-")
                    sp_end = parts[0]
                    cs_start = parts[0]
                    cs_end = parts[1]
                else:
                    return ("-", "-", "-", "-")

            # 提取置信度|Extract confidence
            if "Pr:" in cs_str:
                conf_str = cs_str.split("Pr:")[1].strip()
                confidence = f"{float(conf_str) * 100:.2f}%"
            else:
                confidence = "-"

            return (f"1-{sp_end}", f"{cs_start}-{cs_end}", confidence)
        except Exception:
            return ("-", "-", "-", "-")

    def format_probability(self, prob_str: str) -> str:
        """格式化概率|Format probability"""
        try:
            prob = float(prob_str)
            return f"{prob * 100:.2f}%"
        except ValueError:
            return "N/A"

    def format_line(self, line: str) -> str:
        """格式化单行数据|Format single line of data

        Args:
            line: 原始行|Raw line (e.g., "G35\tSP\t0.000266\t0.999684\tCS pos: 20-21. Pr: 0.9766")

        Returns:
            格式化后的行|Formatted line
        """
        parts = line.strip().split('\t')
        if len(parts) < 4:
            return ""

        protein_id = parts[0]
        prediction = parts[1]
        other_prob = parts[2]
        sp_prob = parts[3]
        cs_info = parts[4] if len(parts) > 4 else ""

        # 获取预测类型信息|Get prediction type info
        pred_info = self.PREDICTION_MAP.get(
            prediction,
            {'zh': prediction, 'en': prediction, 'has_sp': False}
        )

        # 格式化概率|Format probability
        sp_prob_formatted = self.format_probability(sp_prob)

        # 解析切割位点|Parse cleavage site
        if pred_info['has_sp'] and cs_info:
            sp_position, cs_position, confidence = self.parse_cs_position(cs_info)
        else:
            sp_position = "-"
            cs_position = "-"
            confidence = "-"

        # 返回格式化的行|Return formatted line
        return (
            f"{protein_id}\t"
            f"{prediction}\t"
            f"{pred_info['zh']}\t"
            f"{sp_prob_formatted}\t"
            f"{sp_position}\t"
            f"{cs_position}\t"
            f"{confidence}"
        )

    def format(self):
        """执行格式化|Execute formatting"""
        if self.logger:
            self.logger.info("开始格式化SignalP预测结果|Starting SignalP prediction result formatting")
            self.logger.info(f"输入文件|Input: {self.input_file}")
            self.logger.info(f"输出文件|Output: {self.output_file}")

        # 统计信息|Statistics
        total_count = 0
        sp_count = 0
        other_count = 0
        formatted_lines = []

        try:
            # 第一步：读取输入文件，解析header，格式化数据行|Step 1: Read input, parse header, format data
            with open(self.input_file, 'r', encoding='utf-8') as f_in:
                for line in f_in:
                    line = line.strip()
                    if not line:
                        continue

                    # 解析文件头|Parse header
                    if line.startswith('#'):
                        self.parse_header(line)
                        continue

                    # 格式化数据行|Format data line
                    formatted_line = self.format_line(line)
                    if formatted_line:
                        formatted_lines.append(formatted_line)

                        # 统计|Count
                        total_count += 1
                        parts = line.split('\t')
                        if len(parts) > 1 and parts[1] != 'OTHER':
                            sp_count += 1
                        else:
                            other_count += 1

            # 第二步：写入输出文件|Step 2: Write output file
            with open(self.output_file, 'w', encoding='utf-8') as f_out:
                # 写入文件头|Write file header
                f_out.write("# " + "=" * 76 + "\n")
                f_out.write("# SignalP 6.0 信号肽预测结果 (易读格式) | SignalP 6.0 Signal Peptide Prediction Results (Readable Format)\n")
                f_out.write("# " + "=" * 76 + "\n")
                f_out.write(f"# 预测时间 | Prediction Time: {self.timestamp}\n")
                f_out.write(f"# 生物类型 | Organism: {self.organism}\n")
                f_out.write("#\n")
                f_out.write("# 列说明 | Column Description:\n")
                f_out.write("# 1. 蛋白质ID | Protein ID\n")
                f_out.write("# 2. 预测类型代码 | Prediction Code (SP/OTHER/SPI/SPII/SPIII/Tat)\n")
                f_out.write("# 3. 预测结果说明 | Prediction Description (Chinese)\n")
                f_out.write("# 4. 信号肽概率 | Signal Peptide Probability\n")
                f_out.write("# 5. 信号肽位置 | Signal Peptide Position (N-terminal to cleavage site)\n")
                f_out.write("# 6. 切割位点 | Cleavage Site Position (predicted signal peptidase cleavage)\n")
                f_out.write("# 7. 切割置信度 | Cleavage Site Confidence Score\n")
                f_out.write("#\n")
                f_out.write("# 预测类型说明 | Prediction Type Legend:\n")
                for code, info in self.PREDICTION_MAP.items():
                    f_out.write(f"#   {code:6s} = {info['zh']} | {info['en']}\n")
                f_out.write("# " + "=" * 76 + "\n")
                f_out.write(
                    "#蛋白质ID\t预测代码\t预测结果\tSP概率\t信号肽位置\t切割位点\t切割置信度\n"
                )
                f_out.write("# " + "-" * 76 + "\n")

                # 写入数据行|Write data lines
                for formatted_line in formatted_lines:
                    f_out.write(formatted_line + "\n")

                # 写入统计信息|Write statistics
                f_out.write("# " + "-" * 76 + "\n")
                f_out.write(f"# 统计信息 | Statistics:\n")
                f_out.write(f"#   总蛋白质数 | Total proteins: {total_count}\n")
                f_out.write(f"#   信号肽蛋白 | Signal peptide proteins: {sp_count} ({sp_count/total_count*100:.2f}%)\n")
                f_out.write(f"#   非信号肽蛋白 | Non-signal peptide proteins: {other_count} ({other_count/total_count*100:.2f}%)\n")
                f_out.write("# " + "=" * 76 + "\n")

            if self.logger:
                self.logger.info("格式化完成|Formatting completed")
                self.logger.info(f"  总蛋白质数|Total: {total_count}")
                self.logger.info(f"  信号肽数量|Signal peptides: {sp_count} ({sp_count/total_count*100:.2f}%)")
                self.logger.info(f"  非信号肽数量|Non-SP: {other_count} ({other_count/total_count*100:.2f}%)")

            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"格式化失败|Formatting failed: {e}")
            return False

    def format_for_r(self):
        """生成R友好的纯TSV格式|Generate R-friendly pure TSV format

        无注释行，英文列名，适合R read.delim()和Excel导入
        No comment lines, English column names, suitable for R read.delim() and Excel import
        """
        if self.logger:
            self.logger.info("生成R友好格式|Generating R-friendly format")

        total_count = 0
        sp_count = 0
        other_count = 0
        data_rows = []

        try:
            # 读取并格式化数据|Read and format data
            with open(self.input_file, 'r', encoding='utf-8') as f_in:
                for line in f_in:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) < 4:
                        continue

                    protein_id = parts[0]
                    prediction = parts[1]
                    other_prob = parts[2]
                    sp_prob = parts[3]
                    cs_info = parts[4] if len(parts) > 4 else ""

                    # 获取预测类型信息|Get prediction type info
                    pred_info = self.PREDICTION_MAP.get(
                        prediction,
                        {'zh': prediction, 'en': prediction, 'has_sp': False}
                    )

                    # 解析概率为小数|Parse probability as decimal (0-1)
                    try:
                        sp_prob_decimal = float(sp_prob)
                    except ValueError:
                        sp_prob_decimal = 0.0

                    # 解析切割位点|Parse cleavage site
                    sp_start = "NA"
                    sp_end = "NA"
                    cs_start = "NA"
                    cs_end = "NA"
                    confidence = "NA"

                    if pred_info['has_sp'] and cs_info:
                        try:
                            if "CS pos:" in cs_info:
                                pos_part = cs_info.split("CS pos:")[1].split(".")[0].strip()
                                if "-" in pos_part:
                                    pos_parts = pos_part.split("-")
                                    sp_start = "1"
                                    sp_end = pos_parts[0]
                                    cs_start = pos_parts[0]
                                    cs_end = pos_parts[1]

                            if "Pr:" in cs_info:
                                conf_str = cs_info.split("Pr:")[1].strip()
                                confidence = str(float(conf_str))
                        except Exception:
                            pass

                    # 构建数据行|Build data row
                    data_rows.append([
                        protein_id,
                        prediction,
                        "TRUE" if pred_info['has_sp'] else "FALSE",
                        str(sp_prob_decimal),
                        sp_start, sp_end,
                        cs_start, cs_end,
                        confidence
                    ])

                    # 统计|Count
                    total_count += 1
                    if prediction != 'OTHER':
                        sp_count += 1
                    else:
                        other_count += 1

            # 写入TSV文件|Write TSV file
            # 文件名|Filename: signalp_summary.tsv
            import os
            output_dir = os.path.dirname(self.output_file)
            r_friendly_file = os.path.join(output_dir, "signalp_summary.tsv")

            with open(r_friendly_file, 'w', encoding='utf-8') as f_out:
                # 写入表头|Write header (中英文对照，方便汇报|Bilingual headers for reporting)
                f_out.write("蛋白质ID\t预测代码\tHas_SP\tSP概率\tSP起始\tSP终止\t切割位点起始\t切割位点终止\t切割置信度\n")

                # 写入数据|Write data
                for row in data_rows:
                    f_out.write("\t".join(row) + "\n")

            if self.logger:
                self.logger.info(f"R友好格式已保存|R-friendly format saved: {r_friendly_file}")
                self.logger.info(f"  可在R中使用|In R use: data <- read.delim('{r_friendly_file}')")

            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"生成R格式失败|R format generation failed: {e}")
            return False
