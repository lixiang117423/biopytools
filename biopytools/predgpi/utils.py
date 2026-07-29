"""
PredGPI工具函数模块|PredGPI utility functions module
"""

import logging
import os
import sys
from pathlib import Path
from typing import List, Tuple


class PredGPILogger:
    """PredGPI日志管理器|PredGPI logger manager"""

    def __init__(self, output_path: Path, log_name: str = "predgpi.log"):
        self.output_path = Path(str(output_path))
        logs_dir = self.output_path / "99_logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = logs_dir / log_name
        self._setup_logging()

    def _setup_logging(self):
        """设置日志(stdout INFO + stderr WARNING+ + 文件 DEBUG)|Setup logging"""
        self.logger = logging.getLogger("PredGPI")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        file_handler = logging.FileHandler(self.log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger|Get logger"""
        return self.logger


def setup_predgpi_env(predgpi_home: str) -> None:
    """
    设置predgpi运行环境|Setup predgpi runtime environment

    将predgpi_home加入sys.path并设PREDGPI_HOME环境变量,
    使后续 import predgpi 可用|Add predgpi_home to sys.path and
    set PREDGPI_HOME so that `import predgpi` works.

    Args:
        predgpi_home: predgpi安装目录(已展开~)|predgpi install dir (~ expanded)
    """
    # 设环境变量|Set environment variable
    os.environ["PREDGPI_HOME"] = predgpi_home

    # 加入sys.path|Add to sys.path
    if predgpi_home not in sys.path:
        sys.path.insert(0, predgpi_home)


def classify_gpi(fpr: float) -> Tuple[str, float]:
    """
    根据FPR分类GPI锚定蛋白|Classify GPI-anchored protein by FPR

    Args:
        fpr: False Positive Rate 估值|FPR estimate

    Returns:
        (分类标签, 概率值)|(classification label, probability score)
    """
    if fpr <= 0.0015:
        return ("Highly Probable", 1.00)
    elif fpr <= 0.005:
        return ("Probable", 0.70)
    elif fpr <= 0.01:
        return ("Weakly Probable", 0.55)
    else:
        return ("Improbable", 1.00)


# TSV列头|TSV header
TSV_HEADER = [
    "ID",
    "Length",
    "GPI_Anchored",
    "Cleavage_Site",
    "FPR",
    "HMM_LogProb",
    "SVM_Score",
    "Probability",
    "Classification",
]


def write_tsv(records: List[dict], output_file: str) -> None:
    """
    写出预测结果TSV|Write prediction results as TSV

    Args:
        records: 记录列表,每条含 TSV_HEADER 对应键|list of record dicts
        output_file: 输出文件路径|output file path
    """
    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write("\t".join(TSV_HEADER) + "\n")
        for r in records:
            fh.write("\t".join([
                r["ID"],
                str(r["Length"]),
                str(r["GPI_Anchored"]),
                str(r["Cleavage_Site"]),
                f"{r['FPR']:.6f}",
                f"{r['HMM_LogProb']:.6f}",
                f"{r['SVM_Score']:.6f}",
                f"{r['Probability']:.2f}",
                r["Classification"],
            ]) + "\n")


def format_number(num: int) -> str:
    """
    格式化大数字(>1M用M单位)|Format large numbers (>1M as M unit)

    Args:
        num: 数字|number

    Returns:
        格式化字符串|formatted string (e.g. 10.00M)
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


def generate_software_version_yml(config, execution_info: dict = None) -> str:
    """
    生成software_versions.yml到00_pipeline_info(§12.5)|Generate software_versions.yml

    PredGPI无--version命令,版本以发表文献标注|PredGPI has no --version flag;
    version is annotated by publication.

    Args:
        config: PredGPIConfig实例|PredGPIConfig instance
        execution_info: 可选执行信息(start_time/end_time/runtime_seconds)|optional exec info

    Returns:
        写入的文件路径|written file path
    """
    import yaml

    info_dir = os.path.join(config.output_dir, "00_pipeline_info")
    os.makedirs(info_dir, exist_ok=True)
    output_file = os.path.join(info_dir, "software_versions.yml")

    hmm_file = "PHMM.TOT.ss.mod_CSDGN" if config.conservative else "PHMM.TOT.ss.mod"
    info = {
        "pipeline": {"name": "biopytools predgpi", "version": "1.0.0"},
        "tools": {
            "predgpi": {
                "version": "PredGPI (Pierleoni et al., BMC Bioinformatics 2008, 9:392)",
                "path": config.predgpi_home,
                "model": hmm_file,
            },
        },
        "parameters": {
            "input_file": config.input_file,
            "conservative": config.conservative,
            "output_prefix": config.output_prefix,
        },
    }
    if execution_info:
        info["execution"] = execution_info

    with open(output_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(info, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    return output_file
