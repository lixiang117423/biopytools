"""
DeepBSA工具模块|DeepBSA Utilities Module
"""

import logging
import sys
import subprocess
from pathlib import Path


class DeepBSALogger:
    """DeepBSA日志管理类|DeepBSA Logger Management Class"""

    def __init__(self, log_file: Path = None):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
        """
        self.log_file = log_file
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """设置日志器|Setup logger

        Returns:
            配置好的日志器|Configured logger
        """
        logger = logging.getLogger('DeepBSA')
        logger.setLevel(logging.DEBUG)

        # 清除已有的处理器|Clear existing handlers
        logger.handlers.clear()

        # 文件处理器|File handler
        if self.log_file:
            fh = logging.FileHandler(self.log_file, encoding='utf-8')
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)

        # 控制台处理器|Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        logger.addHandler(ch)

        # 格式化器|Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        if self.log_file:
            fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger

        Returns:
            日志器对象|Logger object
        """
        return self.logger


def clean_vcf_comments(input_file: Path, output_file: Path, logger) -> bool:
    """清理VCF文件中的注释行（保留VCF格式必需的##和#CHROM行）|Clean comment lines from VCF file (keep required ## and #CHROM lines)

    Args:
        input_file: 输入VCF文件|Input VCF file
        output_file: 输出文件|Output file
        logger: 日志器|Logger

    Returns:
        bool: 是否成功|Whether successful
    """
    logger.info("=" * 60)
    logger.info("预处理|Preprocessing: 清理VCF注释行|Cleaning VCF comment lines")
    logger.info("=" * 60)

    try:
        # 检查是否有注释行|Check if there are comment lines
        result = subprocess.run(
            f"grep -c '^#' {input_file}",
            shell=True,
            capture_output=True,
            text=True
        )

        comment_count = 0
        if result.returncode == 0:
            comment_count = int(result.stdout.strip())

        if comment_count > 0:
            logger.info(f"发现|Found {comment_count} 行注释（#开头）|comment lines (starting with #)")

            # 对于DeepBSA，保留VCF必需的header行（##和#CHROM）
            # DeepBSA需要这些行来解析FORMAT字段
            # 只删除其他元数据注释行，保留VCF规范的行
            logger.info(f"保留VCF格式必需的header行（##和#CHROM）|Keeping required VCF header lines (## and #CHROM)")
            logger.info(f"生成清理后的文件|Generating cleaned file: {output_file}")

            # 保留##开头（VCF元数据）和#CHROM行，删除其他注释行
            result = subprocess.run(
                f"grep '^##\\|^#CHROM\\|^[^#]' {input_file} > {output_file}",
                shell=True,
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                # 统计清理后的行数|Count lines after cleaning
                line_result = subprocess.run(
                    f"wc -l < {output_file}",
                    shell=True,
                    capture_output=True,
                    text=True
                )
                clean_line_count = line_result.stdout.strip()

                logger.info(f"清理完成|Cleaning completed，保留|retained {clean_line_count} 行数据|data lines")
                return True
            else:
                logger.error("清理失败|Cleaning failed，使用原始文件|using original file")
                return False
        else:
            logger.info("未发现注释行|No comment lines found，直接使用原始文件|using original file directly")
            return False

    except Exception as e:
        logger.error(f"清理过程出错|Error during cleaning: {e}")
        return False


def run_deepbsa_method(method: str, input_file: Path, deepbsa_script: Path,
                       output_dir: Path, conda_env_name: str, logger, parallel: bool = False) -> bool:
    """运行单个DeepBSA方法|Run single DeepBSA method

    Args:
        method: 方法名称|Method name
        input_file: 输入文件|Input file
        deepbsa_script: DeepBSA脚本路径|DeepBSA script path
        output_dir: 输出目录|Output directory
        conda_env_name: Conda环境名称|Conda environment name
        logger: 日志器|Logger
        parallel: 是否并行模式|Whether parallel mode

    Returns:
        bool: 是否成功|Whether successful
    """
    log_file = output_dir / f"{method}.log"

    if parallel:
        logger.info(f"启动|Starting {method} (工作目录|work_dir: {output_dir})")

    cmd = [
        "conda", "run", "-n", conda_env_name, "--no-capture-output",
        "python3",
        str(deepbsa_script),
        "--i", str(input_file),
        "--m", method,
        "--p", "0"  # 禁用pretreatment，使用已过滤的数据|Disable pretreatment, use pre-filtered data
    ]

    if not parallel:
        logger.info(f"命令|Command: {' '.join(cmd)}")

    if parallel:
        # 并行模式：后台运行，输出到日志文件|Parallel mode: run in background, output to log file
        with open(log_file, 'w') as f:
            process = subprocess.Popen(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=str(output_dir)  # 设置工作目录|Set working directory
            )
        return True  # 并行模式不等待结果|Parallel mode doesn't wait for result
    else:
        # 串行模式：实时输出|Serial mode: real-time output
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                cwd=str(output_dir)  # 设置工作目录|Set working directory
            )

            # 保存日志|Save log
            with open(log_file, 'w') as f:
                f.write(result.stdout)

            if result.returncode == 0:
                logger.info(f"{method} 完成|completed")
                return True
            else:
                logger.error(f"{method} 失败|failed")
                # 输出最后20行日志|Output last 20 lines of log
                lines = result.stdout.split('\n')
                if len(lines) > 20:
                    logger.error("错误信息（最后20行）|Error message (last 20 lines):")
                    logger.error('\n'.join(lines[-20:]))
                return False

        except Exception as e:
            logger.error(f"{method} 异常|error: {e}")
            return False


def check_parallel_results(output_dir: Path, methods: list, logger):
    """检查并行运行结果|Check parallel run results

    Args:
        output_dir: 输出目录|Output directory
        methods: 方法列表|Method list
        logger: 日志器|Logger
    """
    logger.info("")
    logger.info("运行结果|Run results:")

    # 各方法结果存放在each/目录下|Individual method results are in 'each/' directory
    each_dir = output_dir / "each"

    for method in methods:
        method_dir = each_dir / method
        log_file = method_dir / f"{method}.log"

        if log_file.exists():
            # 检查日志中是否有错误|Check if log has errors
            with open(log_file, 'r') as f:
                content = f.read()

            # 检查是否有真实错误（Traceback、FileNotFoundError等）|Check for real errors
            has_real_error = any(keyword in content for keyword in [
                'Traceback', 'FileNotFoundError', 'FileExistsError',
                'PermissionError', 'RuntimeError', 'ValueError'
            ])

            # 检查conda错误|Check conda errors
            has_conda_error = 'ERROR conda.cli.main_run' in content

            # 检查进程返回码|Check process return code
            if has_real_error or has_conda_error:
                logger.warning(f"  {method} - 失败|failed，请查看|please check {log_file}")
            else:
                logger.info(f"  {method} - 完成|completed")
        else:
            logger.error(f"  {method} - 日志文件不存在|log file not found: {log_file}")
