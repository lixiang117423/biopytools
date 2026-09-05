"""
NeedleIdentity核心计算模块|Needle Identity Core Calculator
"""

import csv
import itertools
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO

from .utils import build_conda_command


class NeedleIdentityCalculator:
    """needle两两identity计算器|Needle Pairwise Identity Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    ############################################################
    # Public Method
    ############################################################

    def run(self) -> str:
        """运行完整分析|Run full analysis"""
        start_time = datetime.now()
        self.logger.info("开始分析|Starting analysis: pairwise sequence identity")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        # 解析序列|Parse sequences
        records = self._parse_fasta(self.config.input_file)
        self.logger.info(f"序列数量|Sequence count: {len(records)}")

        # 准备输出目录|Prepare output directories
        output_path = Path(self.config.output_dir)
        tmp_dir = output_path / "tmp"
        logs_dir = output_path / "99_logs"
        info_dir = output_path / "00_pipeline_info"
        for d in (tmp_dir, logs_dir, info_dir):
            d.mkdir(parents=True, exist_ok=True)

        # 每条序列写入临时文件|Write per-sequence temp files
        seq_files = self._write_seq_files(records, tmp_dir)

        try:
            # 计算所有两两identity|Compute all pairwise identities
            pair_results = self._compute_all_pairs(records, seq_files, tmp_dir)
            self.logger.info(f"成功计算|Successfully computed: {len(pair_results)} pairs")

            # 输出表格|Write result table
            result_file = self._write_result_table(pair_results, output_path)

            # 版本信息|Write software versions
            self._write_versions(start_time, datetime.now(), info_dir)
        finally:
            # 清理临时目录|Cleanup temp dir
            shutil.rmtree(tmp_dir, ignore_errors=True)

        self.logger.info(f"结果文件|Result file: {result_file}")
        return str(result_file)

    ############################################################
    # Private Method
    ############################################################

    def _parse_fasta(self, fasta_file: str) -> List[dict]:
        """解析FASTA文件|Parse FASTA file"""
        records: List[dict] = []
        seen = set()
        for rec in SeqIO.parse(fasta_file, "fasta"):
            seq_id = rec.id
            if seq_id in seen:
                raise ValueError(f"序列ID重复|Duplicate sequence ID: {seq_id}")
            seen.add(seq_id)
            records.append({"id": seq_id, "seq": str(rec.seq)})
        if len(records) < 2:
            raise ValueError("至少需要2条序列|At least 2 sequences are required")
        return records

    def _write_seq_files(self, records: List[dict], tmp_dir: Path) -> List[Path]:
        """将每条序列写入临时FASTA文件|Write each sequence to a temp FASTA file"""
        seq_files = []
        for idx, rec in enumerate(records):
            seq_file = tmp_dir / f"seq{idx}.fa"
            seq_file.write_text(f">{rec['id']}\n{rec['seq']}\n", encoding='utf-8')
            seq_files.append(seq_file)
        return seq_files

    def _compute_all_pairs(
        self, records: List[dict], seq_files: List[Path], tmp_dir: Path
    ) -> List[dict]:
        """并行计算所有序列对identity|Compute all pairwise identities in parallel"""
        pair_index = list(itertools.combinations(range(len(records)), 2))
        self.logger.info(f"序列对数量|Pair count: {len(pair_index)}")

        result_map: Dict[int, dict] = {}
        with ThreadPoolExecutor(max_workers=self.config.threads) as executor:
            future_to_idx = {
                executor.submit(
                    self._run_pair,
                    records[i], records[j], seq_files[i], seq_files[j], i, j, tmp_dir,
                ): idx
                for idx, (i, j) in enumerate(pair_index)
            }
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                result_map[idx] = future.result()

        # 按原始顺序返回|Return in original pair order
        return [result_map[idx] for idx in range(len(pair_index))]

    def _run_pair(
        self,
        rec_a: dict,
        rec_b: dict,
        file_a: Path,
        file_b: Path,
        idx_a: int,
        idx_b: int,
        tmp_dir: Path,
    ) -> dict:
        """运行单对needle比对|Run needle for one pair"""
        out_file = tmp_dir / f"pair_{idx_a}_{idx_b}.needle"
        args = [
            '-asequence', str(file_a),
            '-bsequence', str(file_b),
            '-outfile', str(out_file),
            '-gapopen', str(self.config.gapopen),
            '-gapextend', str(self.config.gapextend),
            '-aformat', 'pair',
            '-auto',
        ]
        cmd = build_conda_command(self.config.needle_path, args)
        self.logger.info(f"执行|Executing: needle比对|needle alignment: {rec_a['id']} vs {rec_b['id']}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, shell=False, capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"needle执行超时|needle timeout: {rec_a['id']} vs {rec_b['id']}")

        if result.returncode != 0:
            raise RuntimeError(
                f"needle执行失败|needle failed: {rec_a['id']} vs {rec_b['id']}: {result.stderr.strip()}"
            )

        report_text = out_file.read_text(encoding='utf-8') if out_file.exists() else ""
        stats = self._parse_report(report_text)
        stats["seq1"] = rec_a["id"]
        stats["seq2"] = rec_b["id"]
        return stats

    def _parse_report(self, report_text: str) -> dict:
        """解析needle报告|Parse needle report"""
        length_m = re.search(r'^#\s*Length:\s*(\d+)', report_text, re.MULTILINE)
        identity_m = re.search(
            r'^#\s*Identity:\s*(\d+)/(\d+)\s+\(\s*([\d.]+)%\)', report_text, re.MULTILINE
        )
        similarity_m = re.search(
            r'^#\s*Similarity:\s*(\d+)/(\d+)\s+\(\s*([\d.]+)%\)', report_text, re.MULTILINE
        )
        gaps_m = re.search(
            r'^#\s*Gaps:\s*(\d+)/(\d+)\s+\(\s*([\d.]+)%\)', report_text, re.MULTILINE
        )
        score_m = re.search(r'^#\s*Score:\s*([-\d.]+)', report_text, re.MULTILINE)

        if not length_m or not identity_m:
            raise ValueError("needle报告解析失败|Failed to parse needle report")

        return {
            "aligned_length": int(length_m.group(1)),
            "matches": int(identity_m.group(1)),
            "identity_percent": float(identity_m.group(3)),
            "similarity_percent": float(similarity_m.group(3)) if similarity_m else None,
            "gaps": int(gaps_m.group(1)) if gaps_m else None,
            "score": float(score_m.group(1)) if score_m else None,
        }

    def _write_result_table(self, pair_results: List[dict], output_path: Path) -> Path:
        """写入结果表格(TSV)|Write result table (TSV)"""
        result_file = output_path / f"{self.config.output_prefix}_needle_identity.tsv"
        header = [
            "seq1", "seq2", "identity_percent", "matches",
            "aligned_length", "gaps", "similarity_percent", "score",
        ]
        with open(result_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in pair_results:
                writer.writerow({key: row.get(key) for key in header})
        self.logger.info(f"写入表格|Write table: {result_file}")
        return result_file

    def _write_versions(self, start_time: datetime, end_time: datetime, info_dir: Path) -> Path:
        """写入软件版本信息|Write software versions"""
        import yaml

        version_info = self._get_emboss_version()
        runtime_seconds = int((end_time - start_time).total_seconds())
        info = {
            'pipeline': {
                'name': 'biopytools needle_identity',
                'version': '1.0.0',
            },
            'tools': {
                'needle': {
                    'version': version_info,
                    'path': self.config.needle_path,
                    'command': 'embossversion',
                },
            },
            'parameters': {
                'threads': self.config.threads,
                'gapopen': self.config.gapopen,
                'gapextend': self.config.gapextend,
            },
            'execution': {
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime_seconds,
            },
        }
        out_file = info_dir / 'software_versions.yml'
        with open(out_file, 'w', encoding='utf-8') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
        self.logger.info(f"写入版本信息|Write versions: {out_file}")
        return out_file

    def _get_emboss_version(self) -> str:
        """获取EMBOSS版本|Get EMBOSS version"""
        embossversion_path = os.path.join(os.path.dirname(self.config.needle_path), 'embossversion')
        try:
            cmd = build_conda_command(embossversion_path, [])
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, shell=False, capture_output=True, text=True, timeout=60)
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip().splitlines()[0]
        except Exception:
            pass
        return "unknown"
