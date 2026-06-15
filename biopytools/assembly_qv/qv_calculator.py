"""
Merqury QV计算核心模块|Merqury QV Calculation Core Module
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

from .utils import run_command, parse_qv_output, find_r1_r2_pairs, detect_data_type, build_conda_command_string


class MerquryQVCalculator:
    """Merqury QV计算器|Merqury QV Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def build_meryl_db(self) -> Tuple[bool, str]:
        """
        构建Meryl k-mer数据库|Build Meryl k-mer database

        Returns:
            tuple: (成功状态, meryl数据库路径)|(success, meryl_db_path)
        """
        self.logger.info("开始构建Meryl k-mer数据库|Starting Meryl k-mer database construction")

        # 检测数据类型|Detect data type
        data_type = self.config.data_type
        if data_type == "auto":
            data_type = detect_data_type(self.config.fastq_files)

        self.logger.info(f"检测到数据类型|Detected data type: {data_type}")

        # 准备输出路径|Prepare output path
        meryl_db = self.config.output_path / "reads.meryl"

        # 根据数据类型选择k值|Select k value based on data type
        k = self.config.kmer_size
        if k is None:
            if data_type == "hifi":
                k = 21  # HiFi推荐使用较大k值|HiFi recommends larger k
            else:
                k = 21  # Illumina默认k值|Illumina default k value
            self.logger.info(f"自动选择k值|Auto-selected k value: {k}")

        # 检查是否已存在|Check if already exists
        if meryl_db.exists():
            self.logger.info(f"Meryl数据库已存在，跳过构建|Meryl database already exists, skipping build: {meryl_db}")
            return True, str(meryl_db)

        # 构建meryl命令|Build meryl command
        if data_type == "illumina":
            # 配对数据|Paired-end data
            pairs, single_end = find_r1_r2_pairs(self.config.fastq_files)

            # 优先使用配对数据|Prefer paired-end data
            if pairs:
                r1, r2 = pairs[0]
                command = [
                    f"meryl count",
                    f"k={k}",
                    f"output {meryl_db}",
                    f"{r1}",
                    f"{r2}"
                ]
            elif single_end:
                # 如果没有配对数据，使用单端数据|If no pairs, use single-end
                fastq_files = " ".join(single_end[:5])
                command = [
                    f"meryl count",
                    f"k={k}",
                    f"output {meryl_db}",
                    fastq_files
                ]
            else:
                self.logger.error("未找到FASTQ文件|No FASTQ files found")
                return False, ""
        else:
            # HiFi或单端数据|HiFi or single-end data
            fastq_files = " ".join(self.config.fastq_files[:5])  # 限制文件数量避免命令过长
            command = [
                f"meryl count",
                f"k={k}",
                f"output {meryl_db}",
                fastq_files
            ]

        # 设置线程数|Set threads
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(self.config.threads)

        # 确保输出目录存在|Ensure output directory exists
        self.config.output_path.mkdir(parents=True, exist_ok=True)

        # 构建完整命令字符串|Build full command string
        full_command = " ".join(command)

        # 使用conda wrapper包装命令|Use conda wrapper to wrap command
        wrapped_command = build_conda_command_string("meryl", full_command)

        # 执行命令|Execute command
        try:
            result = subprocess.run(
                wrapped_command,
                shell=True,
                env=env,
                capture_output=True,
                text=True,
                check=True
            )

            self.logger.info(f"Meryl数据库构建完成|Meryl database built successfully: {meryl_db}")
            return True, str(meryl_db)

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Meryl数据库构建失败|Meryl database build failed: {e.stderr}")
            return False, ""

    def calculate_qv(self, meryl_db: str) -> Dict:
        """
        计算QV值|Calculate QV value

        Args:
            meryl_db: Meryl数据库路径|Meryl database path

        Returns:
            dict: QV统计结果|QV statistics results
        """
        self.logger.info("开始计算QV值|Starting QV calculation")

        # 准备输出路径|Prepare output path
        output_prefix = self.config.output_path / "qv_result"

        # 构建Merqury QV命令|Build Merqury QV command
        qv_script = os.path.join(self.config.merqury_path, "eval", "qv.sh")

        if not os.path.exists(qv_script):
            self.logger.error(f"Merqury QV脚本不存在|Merqury QV script not found: {qv_script}")
            return {}

        command = [
            "bash",
            qv_script,
            meryl_db,
            self.config.genome_file,
            str(output_prefix)
        ]

        # 执行命令|Execute command
        success, stdout = run_command(
            command,
            self.config.conda_env,
            self.logger,
            "计算QV值|Calculate QV"
        )

        if not success:
            self.logger.error("QV计算失败|QV calculation failed")
            return {}

        # 解析结果|Parse results
        qv_file = str(output_prefix) + ".qv"
        if os.path.exists(qv_file):
            # 添加表头到QV文件|Add header to QV file
            qv_file_with_header = self._add_header_to_qv_file(qv_file)

            stats = parse_qv_output(qv_file_with_header)
            self.logger.info(f"QV值计算完成|QV calculation completed: {qv_file_with_header}")
            return stats
        else:
            self.logger.warning(f"QV文件未生成|QV file not generated: {qv_file}")
            # 尝试从stdout解析|Try to parse from stdout
            return {"raw_output": stdout}

    def _add_header_to_qv_file(self, qv_file: str) -> str:
        """
        为QV文件添加表头|Add header to QV file

        Args:
            qv_file: 原始QV文件路径|Original QV file path

        Returns:
            str: 带表头的QV文件路径|QV file path with header
        """
        # 读取原始数据|Read original data
        with open(qv_file, 'r') as f:
            original_data = f.read()

        # 创建带表头的新文件|Create new file with header
        header = "# 序列名称|Sequence_name\tN50/Contig数|N50_Contigs\t总碱基数|Total_bases\tQV值|QV_value\t错误率|Error_rate\n"
        qv_file_with_header = qv_file.replace('.qv', '_with_header.qv')

        with open(qv_file_with_header, 'w') as f:
            f.write(header)
            f.write(original_data)

        # 删除原文件，重命名新文件|Delete original, rename new file
        import os
        os.remove(qv_file)
        os.rename(qv_file_with_header, qv_file)

        self.logger.info(f"已添加表头到QV文件|Header added to QV file: {qv_file}")
        return qv_file

    def run_full_analysis(self) -> Dict:
        """
        运行完整的QV分析流程|Run full QV analysis pipeline

        Returns:
            dict: QV统计结果|QV statistics results
        """
        # 步骤1: 构建Meryl数据库|Step 1: Build Meryl database
        success, meryl_db = self.build_meryl_db()
        if not success:
            self.logger.error("Meryl数据库构建失败，终止分析|Meryl database build failed, aborting analysis")
            return {}

        # 步骤2: 计算QV值|Step 2: Calculate QV
        qv_stats = self.calculate_qv(meryl_db)

        return qv_stats

    def get_summary(self) -> str:
        """
        获取结果摘要|Get results summary

        Returns:
            str: 格式化的摘要文本|Formatted summary text
        """
        output_prefix = self.config.output_path / "qv_result"
        qv_file = str(output_prefix) + ".qv"

        if not os.path.exists(qv_file):
            return "QV结果文件不存在|QV result file does not exist"

        # 解析QV文件|Parse QV file
        from .utils import parse_qv_output
        stats = parse_qv_output(qv_file)

        summary = []
        summary.append("="*80)
        summary.append("Merqury QV分析结果摘要|Merqury QV Analysis Results Summary")
        summary.append("="*80)
        summary.append(f"基因组文件|Genome file: {self.config.genome_file}")
        summary.append(f"FASTQ目录|FASTQ directory: {self.config.fastq_dir}")
        summary.append(f"输出目录|Output directory: {self.config.output_dir}")
        summary.append("")

        if 'qv_value' in stats:
            # 显示解析后的数据|Display parsed data
            summary.append("-"*80)
            summary.append("QV计算结果|QV Calculation Results:")
            summary.append("-"*80)
            summary.append(f"  序列名称|Sequence name: {stats.get('sequence_name', 'N/A')}")
            summary.append(f"  N50/Contig数|N50/Contigs: {stats.get('n50_or_contigs', 'N/A')}")
            summary.append(f"  总碱基数|Total bases: {stats.get('total_bases', 'N/A')}")
            summary.append(f"  QV值|QV value: {stats.get('qv_value', 'N/A')}")
            summary.append(f"  错误率|Error rate: {stats.get('error_rate', 'N/A')}")
            summary.append("-"*80)
            summary.append("")
            summary.append("列说明|Column descriptions:")
            summary.append("  列1|Col1: 序列名称|Sequence name")
            summary.append("  列2|Col2: N50或contig数量|N50 or contig count")
            summary.append("  列3|Col3: 总碱基数量|Total base count")
            summary.append("  列4|Col4: QV值 (质量值)|QV value (Quality Value)")
            summary.append("  列5|Col5: 错误率|Error rate")
        else:
            # 显示原始数据|Display raw data
            summary.append("-"*80)
            summary.append("原始数据|Raw data:")
            summary.append("-"*80)
            with open(qv_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        summary.append(line)

        summary.append("="*80)

        return "\n".join(summary)
