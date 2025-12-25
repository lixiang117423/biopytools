"""
EGAPx批量工具类 | EGAPx Batch Utility Classes
包含日志和辅助函数 | Contains logging and helper functions
"""

import logging
import sys
import os
import subprocess
from typing import List, Tuple
import re


class EGAPxBatchLogger:
    """EGAPx批量日志管理器 | EGAPx Batch Logger Manager"""

    def __init__(self, output_dir: str):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            output_dir: 输出目录 | Output directory
        """
        self.output_dir = output_dir
        log_file = os.path.join(output_dir, 'egapx_batch.log')
        self.setup_logging(log_file)

    def setup_logging(self, log_file: str):
        """设置日志 | Setup logging"""
        self.logger = logging.getLogger("EGAPxBatch")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式 | Log format
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler | File handler
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class GenomeSplitter:
    """基因组拆分器 | Genome Splitter"""

    def __init__(self, logger):
        """
        初始化拆分器 | Initialize splitter

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def split_genome(self, genome_file: str, output_dir: str, chr_prefix: str = None) -> List[str]:
        """
        拆分基因组文件 | Split genome file

        Args:
            genome_file: 基因组文件路径 | Genome file path
            output_dir: 输出目录 | Output directory
            chr_prefix: 染色体前缀过滤 | Chromosome prefix filter

        Returns:
            拆分后的文件列表 | List of split files
        """
        self.logger.info("开始拆分基因组文件... | Starting to split genome file...")

        os.makedirs(output_dir, exist_ok=True)

        # 使用awk拆分 | Use awk to split
        split_dir = os.path.join(output_dir, "temp_split")
        os.makedirs(split_dir, exist_ok=True)

        # 构建awk命令 | Build awk command
        filter_pattern = chr_prefix if chr_prefix else ""

        awk_script = f'''
/^>/ {{
    id = substr($1, 2)
    if ("{filter_pattern}" != "" && index(id, "{filter_pattern}") != 1) {{
        filename = ""
        next
    }}
    gsub(/[:|\\/\\\\]/, "_", id)
    filename = "{split_dir}/" id ".fa"
    print > filename
    next
}}
filename != "" {{ print >> filename }}
'''

        try:
            with open(genome_file, 'r') as f:
                subprocess.run(
                    ['awk', awk_script],
                    stdin=f,
                    capture_output=True,
                    text=True,
                    check=True
                )

            # 获取拆分后的文件列表 | Get split file list
            split_files = sorted([
                f for f in os.listdir(split_dir)
                if f.endswith('.fa')
            ])

            self.logger.info(f"成功拆分 {len(split_files)} 个序列 | Successfully split {len(split_files)} sequences")

            if len(split_files) == 0:
                self.logger.error("未拆分出任何序列 | No sequences were split")

            return [os.path.join(split_dir, f) for f in split_files]

        except subprocess.CalledProcessError as e:
            self.logger.error(f"拆分基因组失败 | Failed to split genome: {e.stderr}")
            return []

    def extract_sequences(self, genome_file: str, chr_prefix: str = None) -> List[Tuple[str, str]]:
        """
        提取序列信息 | Extract sequence information

        Args:
            genome_file: 基因组文件 | Genome file
            chr_prefix: 染色体前缀过滤 | Chromosome prefix filter

        Returns:
            (序列ID, 文件路径)列表 | List of (sequence ID, file path)
        """
        self.logger.info("正在提取序列信息... | Extracting sequence information...")

        sequences = []
        try:
            with open(genome_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        current_id = line[1:].split()[0]
                        # 清理ID中的特殊字符 | Clean special characters in ID
                        current_id = re.sub(r'[:|/\\]', '_', current_id)

                        # 应用前缀过滤 | Apply prefix filter
                        if chr_prefix and not current_id.startswith(chr_prefix):
                            continue

                        # 只在遇到序列头时添加一次 | Add only when encountering sequence header
                        sequences.append((current_id, None))

            self.logger.info(f"找到 {len(sequences)} 个序列 | Found {len(sequences)} sequences")
            return sequences

        except Exception as e:
            self.logger.error(f"读取序列信息失败 | Failed to read sequence info: {e}")
            return []


class TemplateProcessor:
    """模板处理器 | Template Processor"""

    # 内置YAML模板 | Built-in YAML template
    DEFAULT_YAML_TEMPLATE = """genome: {genome_path}
taxid: {taxid}
short_reads: {short_reads}
long_reads: {long_reads}
locus_tag_prefix: {locus_tag_prefix}"""

    # 内置Shell脚本模板 | Built-in Shell script template
    DEFAULT_SCRIPT_TEMPLATE = """#!/bin/bash
# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx
python3 \\
    ui/egapx.py \\
    {yaml_path} \\
    -e singularity \\
    -w {work_dir} \\
    -o {output_dir} \\
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \\
    -r {report_name}
"""

    def __init__(self, logger):
        """
        初始化处理器 | Initialize processor

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def extract_from_yaml(self, yaml_file: str) -> dict:
        """
        从YAML模板提取信息 | Extract information from YAML template

        Args:
            yaml_file: YAML文件路径 | YAML file path

        Returns:
            提取的信息字典 | Dictionary of extracted information
        """
        self.logger.info("正在解析YAML模板... | Parsing YAML template...")

        info = {
            'genome_path': None,
            'locus_tag_prefix': None,
            'taxid': '71234',
            'short_reads': '',
            'long_reads': ''
        }

        try:
            with open(yaml_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('genome:'):
                        info['genome_path'] = line.split(':', 1)[1].strip()
                    elif line.startswith('locus_tag_prefix:'):
                        info['locus_tag_prefix'] = line.split(':', 1)[1].strip().strip('"')
                    elif line.startswith('taxid:'):
                        info['taxid'] = line.split(':', 1)[1].strip()
                    elif line.startswith('short_reads:'):
                        info['short_reads'] = line.split(':', 1)[1].strip()
                    elif line.startswith('long_reads:'):
                        info['long_reads'] = line.split(':', 1)[1].strip()

            self.logger.info(f"YAML信息提取完成 | YAML info extracted:")
            for key, value in info.items():
                if value:
                    self.logger.info(f"  - {key}: {value}")

            return info

        except Exception as e:
            self.logger.error(f"读取YAML模板失败 | Failed to read YAML template: {e}")
            return info

    def extract_from_script(self, script_file: str) -> dict:
        """
        从脚本模板提取信息 | Extract information from script template

        Args:
            script_file: 脚本文件路径 | Script file path

        Returns:
            提取的信息字典 | Dictionary of extracted information
        """
        self.logger.info("正在解析脚本模板... | Parsing script template...")

        info = {
            'yaml_path': None,
            'work_dir': None,
            'output_dir': None,
            'report_name': 'EGAPx'
        }

        try:
            with open(script_file, 'r') as f:
                content = f.read()

                # 提取YAML路径 | Extract YAML path
                yaml_match = re.search(r'ui/egapx\.py\s+(\S+\.yaml)', content)
                if yaml_match:
                    info['yaml_path'] = yaml_match.group(1)
                else:
                    # 尝试其他模式 | Try other patterns
                    yaml_match = re.search(r'\S+\.yaml', content)
                    if yaml_match:
                        info['yaml_path'] = yaml_match.group(0)

                # 提取work目录 | Extract work directory
                work_match = re.search(r'-w\s+(\S+)', content)
                if work_match:
                    info['work_dir'] = work_match.group(1).strip('"')

                # 提取output目录 | Extract output directory
                output_match = re.search(r'-o\s+(\S+)', content)
                if output_match:
                    info['output_dir'] = output_match.group(1).strip('"')

                # 提取report名称 | Extract report name
                report_match = re.search(r'-r\s+(\S+)', content)
                if report_match:
                    info['report_name'] = report_match.group(1).strip('"')

            self.logger.info(f"脚本信息提取完成 | Script info extracted:")
            for key, value in info.items():
                if value:
                    self.logger.info(f"  - {key}: {value}")

            return info

        except Exception as e:
            self.logger.error(f"读取脚本模板失败 | Failed to read script template: {e}")
            return info


class JobGenerator:
    """作业生成器 | Job Generator"""

    def __init__(self, logger):
        """
        初始化生成器 | Initialize generator

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def generate_yaml(self, genome: str, locus_tag_prefix: str,
                      taxid: str = "71234", short_reads: str = "",
                      long_reads: str = "") -> str:
        """
        生成YAML配置文件 | Generate YAML config file

        Args:
            genome: 基因组路径 | Genome path
            locus_tag_prefix: locus标签前缀 | Locus tag prefix
            taxid: 物种ID | Taxonomy ID
            short_reads: 短读长文件路径 | Short reads file path
            long_reads: 长读长文件路径 | Long reads file path

        Returns:
            生成的YAML内容 | Generated YAML content
        """
        return TemplateProcessor.DEFAULT_YAML_TEMPLATE.format(
            genome_path=genome,
            taxid=taxid,
            short_reads=short_reads,
            long_reads=long_reads,
            locus_tag_prefix=locus_tag_prefix
        )

    def generate_script(self, yaml_path: str, work_dir: str,
                        output_dir: str, report_name: str, chr_dir: str) -> str:
        """
        生成运行脚本 | Generate run script

        Args:
            yaml_path: YAML文件路径 | YAML file path
            work_dir: 工作目录 | Work directory
            output_dir: 输出目录 | Output directory
            report_name: 报告名称 | Report name
            chr_dir: 染色体目录 | Chromosome directory

        Returns:
            生成的脚本内容 | Generated script content
        """
        return TemplateProcessor.DEFAULT_SCRIPT_TEMPLATE.format(
            yaml_path=yaml_path,
            work_dir=work_dir,
            output_dir=output_dir,
            report_name=report_name
        )
