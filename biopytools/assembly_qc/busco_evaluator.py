"""
BUSCO完整性评估模块（完全独立实现）|BUSCO Completeness Evaluation Module (Fully Independent)

直接调用BUSCO命令行工具，不依赖其他模块
Directly calls BUSCO command-line tool, no dependency on other modules
"""

import os
import json
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional, List
from .utils import get_conda_env_from_path, build_conda_command


class BUSCOEvaluator:
    """BUSCO完整性评估器|BUSCO Completeness Evaluator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        # 使用自己的命令执行器|Use own command runner
        self.working_dir = Path(self.config.busco_output_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

        # 从conda环境路径提取环境名|Extract env name from conda env path
        self.busco_env_name = get_conda_env_from_path(self.config.conda_env_busco)

    def evaluate(self) -> Optional[Dict[str, Any]]:
        """运行BUSCO评估|Run BUSCO evaluation"""
        if self.config.skip_busco:
            self.logger.info("跳过BUSCO评估|Skipping BUSCO evaluation")
            return None

        self.logger.info("开始BUSCO完整性评估|Starting BUSCO completeness evaluation")

        # 检查是否已完成|Check if already completed
        # BUSCO v6输出结构：busco_output/short_summary.specific.{lineage}.{output_name}.json
        json_patterns = [
            self.working_dir / "busco_output" / f"short_summary.specific.{self.config.lineage}.busco_output.json",
            self.working_dir / "busco_output" / "short_summary.specific.busco_output.json",
        ]

        already_completed = False
        for json_file in json_patterns:
            if self.config.resume and json_file.exists():
                self.logger.info("BUSCO评估已完成，跳过|BUSCO evaluation already completed, skipping")
                return self._parse_json_results(json_file)

        # 构建并运行BUSCO命令|Build and run BUSCO command
        cmd = self._build_busco_command()

        try:
            self.logger.info(f"运行BUSCO|Running BUSCO: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                cwd=self.working_dir,
                capture_output=True,
                text=True,
                check=False,
                env=os.environ.copy()
            )

            if result.returncode != 0:
                self.logger.error(f"BUSCO运行失败|BUSCO run failed")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                self.logger.error(f"STDERR:\n{result.stderr}")
                self.logger.error(f"STDOUT:\n{result.stdout}")
                return None

            # 解析结果|Parse results
            return self._parse_results()

        except Exception as e:
            self.logger.error(f"BUSCO评估异常|BUSCO evaluation error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _build_busco_command(self) -> list:
        """构建BUSCO命令列表|Build BUSCO command list"""
        # 构建BUSCO数据集的完整路径|Build full path to BUSCO dataset
        # 使用完整路径模式，配合--offline参数|Use full path mode with --offline
        lineage_full_path = os.path.join(self.config.busco_dataset_path, self.config.lineage)

        # BUSCO v6不需要analyze子命令|BUSCO v6 doesn't need analyze subcommand
        args = [
            "-i", self.config.genome,
            "-l", lineage_full_path,  # 使用完整路径|Use full path to dataset
            "-o", "busco_output",  # 输出名称，BUSCO会创建run_busco_output目录
            "-m", self.config.busco_mode,
            "-c", str(self.config.busco_threads),
            "--offline"  # 离线模式，不使用--download_path|Offline mode without --download_path
        ]

        # 使用conda run调用busco|Use conda run to call busco
        return build_conda_command(self.busco_env_name, "busco", args)

    def _parse_results(self) -> Optional[Dict[str, Any]]:
        """解析BUSCO结果|Parse BUSCO results"""
        # BUSCO v6输出结构：busco_output/short_summary.specific.{lineage}.{output_name}.json
        # 优先查找JSON文件|Prefer JSON file
        json_patterns = [
            self.working_dir / "busco_output" / f"short_summary.specific.{self.config.lineage}.busco_output.json",
            self.working_dir / "busco_output" / "short_summary.specific.*.busco_output.json",
            self.working_dir / "busco_output" / "short_summary.generic.*.busco_output.json"
        ]

        import glob
        for pattern in json_patterns:
            matches = glob.glob(str(pattern))
            if matches:
                return self._parse_json_results(Path(matches[0]))

        # 回退到TXT文件|Fallback to TXT file
        txt_patterns = [
            self.working_dir / "busco_output" / f"short_summary.specific.{self.config.lineage}.busco_output.txt",
            self.working_dir / "busco_output" / "short_summary.specific.*.busco_output.txt",
            self.working_dir / "busco_output" / "short_summary.txt"
        ]
        for txt_pattern in txt_patterns:
            matches = glob.glob(str(txt_pattern))
            if matches:
                return self._parse_txt_results(Path(matches[0]))

        self.logger.error("未找到BUSCO结果文件|BUSCO result files not found")
        self.logger.debug(f"查找目录|Search directory: {self.working_dir / 'busco_output'}")
        return None

    def _parse_json_results(self, json_file: Path) -> Optional[Dict[str, Any]]:
        """解析BUSCO JSON结果|Parse BUSCO JSON results"""
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            results = data.get('results', {})
            lineage = data.get('lineage_dataset', {})

            return {
                'lineage': lineage.get('name', self.config.lineage),
                'complete': results.get('Complete percentage', 0),
                'single': results.get('Single copy percentage', 0),
                'duplicated': results.get('Multi copy percentage', 0),
                'fragmented': results.get('Fragmented percentage', 0),
                'missing': results.get('Missing percentage', 0),
                'complete_count': results.get('Complete BUSCOs', 0),
                'single_count': results.get('Single copy BUSCOs', 0),
                'duplicated_count': results.get('Multi copy BUSCOs', 0),
                'fragmented_count': results.get('Fragmented BUSCOs', 0),
                'missing_count': results.get('Missing BUSCOs', 0),
            }

        except Exception as e:
            self.logger.error(f"解析BUSCO JSON结果失败|Failed to parse BUSCO JSON results: {e}")
            return None

    def _parse_txt_results(self, txt_file: Path) -> Optional[Dict[str, Any]]:
        """解析BUSCO TXT结果|Parse BUSCO TXT results"""
        try:
            with open(txt_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # 解析BUSCO输出的总结行|Parse BUSCO summary line
            # 格式: Complete: 95.5%[S:94.2%,D:1.3%], Fragmented: 2.5%, Missing: 2.0%
            import re

            complete_match = re.search(r'Complete:\s*([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]', content)
            if not complete_match:
                # 尝试另一种格式|Try another format
                complete_match = re.search(r'Complete:\s*([\d.]+)%', content)

            fragmented_match = re.search(r'Fragmented:\s*([\d.]+)%', content)
            missing_match = re.search(r'Missing:\s*([\d.]+)%', content)

            if complete_match:
                complete = float(complete_match.group(1))
                single = float(complete_match.group(2)) if len(complete_match.groups()) >= 2 else complete
                duplicated = float(complete_match.group(3)) if len(complete_match.groups()) >= 3 else 0
            else:
                complete = 0.0
                single = 0.0
                duplicated = 0.0

            fragmented = float(fragmented_match.group(1)) if fragmented_match else 0.0
            missing = float(missing_match.group(1)) if missing_match else 0.0

            return {
                'lineage': self.config.lineage,
                'complete': complete,
                'single': single,
                'duplicated': duplicated,
                'fragmented': fragmented,
                'missing': missing,
                'complete_count': 0,
                'single_count': 0,
                'duplicated_count': 0,
                'fragmented_count': 0,
                'missing_count': 0,
            }

        except Exception as e:
            self.logger.error(f"解析BUSCO TXT结果失败|Failed to parse BUSCO TXT results: {e}")
            return None
