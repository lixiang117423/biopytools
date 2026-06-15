"""
等位基因鉴定运行器|Allelic Gene Identification Runner

基于JcviPipeline共享管道, 添加等位基因对提取和汇总功能
Based on JcviPipeline shared pipeline, adds allelic pair extraction and summary
"""

import os
import time
from pathlib import Path
from datetime import datetime
from typing import Optional, List

from ..config import JcviBaseConfig
from ..pipeline import JcviPipeline, PipelineResult
from ..utils import JcviLogger, get_sample_name


class AllelicRunner:
    """等位基因鉴定运行器|Allelic Gene Identification Runner"""

    VERSION = "1.0.0"

    def __init__(self, config: JcviBaseConfig, logger: Optional[JcviLogger] = None):
        self.config = config
        self.logger_obj = logger
        self.logger = logger.get_logger() if logger else None
        self.start_time = None

    def run(self) -> bool:
        try:
            self.start_time = time.time()

            if self.logger_obj is None:
                log_dir = Path(self.config.output_dir) / "99_logs"
                log_dir.mkdir(parents=True, exist_ok=True)
                self.logger_obj = JcviLogger(log_dir / "allelic_genes.log", "AllelicGenes")
                self.logger = self.logger_obj.get_logger()

            self._print_header()
            self.config.validate()

            pipeline = JcviPipeline(self.config, self.logger_obj)
            result = pipeline.run()

            if not result.pair_dirs:
                self._print_footer(0, len(result.pair_list), len(result.pair_list))
                return False

            # 步骤5: 提取等位基因对
            self._log_section("步骤5/6: 提取等位基因对|Step 5/6: Extract allelic gene pairs")

            all_pairs_output = []
            success_count = 0
            fail_count = 0

            for (name_a, name_b), pair_dir in result.pair_dirs.items():
                pairs_file = str(pair_dir / f"{name_a}_{name_b}.allelic_pairs.txt")
                lifted_anchors = str(pair_dir / f"{name_a}.{name_b}.lifted.anchors")
                bed_a = str(result.bed_dir / f"{name_a}.uniq.bed")
                bed_b = str(result.bed_dir / f"{name_b}.uniq.bed")

                self._extract_pairs(lifted_anchors, pairs_file, bed_a, bed_b)

                if Path(pairs_file).exists() and Path(pairs_file).stat().st_size > 0:
                    gene_count = self._count_genes(pairs_file)
                    all_pairs_output.append(pairs_file)
                    success_count += 1
                    self.logger.info(f"  完成: 找到 {gene_count} 对等位基因|"
                                    f"Done: {gene_count} allelic gene pairs found")
                else:
                    fail_count += 1
                    self.logger.warning(
                        f"  {name_a} vs {name_b} 未找到等位基因|"
                        f"No allelic genes found for {name_a} vs {name_b}")

            # 步骤6: 汇总
            self._log_section("步骤6/6: 汇总结果|Step 6/6: Summarize results")
            self._summarize_results(all_pairs_output)

            self._print_footer(success_count, fail_count, len(result.pair_list))
            return fail_count == 0

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Program execution error: {e}", exc_info=True)
            return False

    def _extract_pairs(self, anchors_file: str, pairs_file: str,
                       bed_file_a: str, bed_file_b: str):
        self.logger.info("  提取等位基因对|Extracting allelic gene pairs...")

        bed_a = self._load_bed_index(bed_file_a)
        bed_b = self._load_bed_index(bed_file_b)

        count = 0
        with open(anchors_file, 'r') as fin, open(pairs_file, 'w') as fout:
            fout.write("#gene_a\tgene_b\tchrom_a\tstart_a\tend_a\tstrand_a"
                       "\tchrom_b\tstart_b\tend_b\tstrand_b"
                       "\tscore\tpair_type\n")
            for line in fin:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("###"):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    gene_a, gene_b, score_raw = parts[0], parts[1], parts[2]

                    if score_raw.endswith('L'):
                        pair_type = "lifted"
                        score = score_raw[:-1]
                    else:
                        pair_type = "direct"
                        score = score_raw

                    info_a = bed_a.get(gene_a, (".", ".", ".", "."))
                    info_b = bed_b.get(gene_b, (".", ".", ".", "."))

                    fout.write(f"{gene_a}\t{gene_b}"
                               f"\t{info_a[0]}\t{info_a[1]}\t{info_a[2]}\t{info_a[3]}"
                               f"\t{info_b[0]}\t{info_b[1]}\t{info_b[2]}\t{info_b[3]}"
                               f"\t{score}\t{pair_type}\n")
                    count += 1

        self.logger.info(f"  提取 {count} 对等位基因|Extracted {count} allelic gene pairs")

    @staticmethod
    def _load_bed_index(bed_file: str) -> dict:
        index = {}
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split('\t')
                if len(parts) >= 6:
                    index[parts[3]] = (parts[0], parts[1], parts[2], parts[5])
        return index

    def _count_genes(self, pairs_file: str) -> int:
        count = 0
        with open(pairs_file, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                count += 1
        return count

    def _summarize_results(self, pair_files: List[str]):
        summary_file = Path(self.config.output_dir) / "all_allelic_pairs.txt"

        total_pairs = 0
        with open(summary_file, 'w') as fout:
            fout.write("#sample_a\tsample_b\tgene_a\tgene_b"
                       "\tchrom_a\tstart_a\tend_a\tstrand_a"
                       "\tchrom_b\tstart_b\tend_b\tstrand_b"
                       "\tscore\tpair_type\n")

            for pf in sorted(pair_files):
                stem = Path(pf).parent.name
                names = stem.split('_vs_')
                if len(names) == 2:
                    sample_a, sample_b = names
                else:
                    parts = Path(pf).stem.split('_')
                    sample_a, sample_b = parts[0], parts[-1]

                with open(pf, 'r') as fin:
                    for line in fin:
                        if line.startswith("#") or not line.strip():
                            continue
                        cols = line.strip().split('\t')
                        fout.write(f"{sample_a}\t{sample_b}\t"
                                   + "\t".join(cols) + "\n")
                        total_pairs += 1

        self.logger.info(f"  汇总文件|Summary file: {summary_file}")
        self.logger.info(f"  总等位基因对|Total allelic gene pairs: {total_pairs}")

    def _log_section(self, title: str):
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(title)
        self.logger.info("=" * 60)

    def _print_header(self):
        if not self.logger:
            return
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("JCVI 等位基因鉴定|JCVI Allelic Gene Identification")
        self.logger.info("=" * 60)
        self.logger.info(f"版本|Version: {self.VERSION}")
        self.logger.info(f"时间|Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  conda环境|Conda env: {self.config.conda_env}")
        self.logger.info(f"  dbtype: {self.config.dbtype}")
        self.logger.info(f"  C-score: {self.config.cscore}")
        self.logger.info(f"  比对软件|Align soft: {self.config.align_soft}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info("")

    def _print_footer(self, success: int, fail: int, total: int):
        if not self.logger:
            return
        elapsed = time.time() - self.start_time
        hours, remainder = divmod(int(elapsed), 3600)
        minutes, seconds = divmod(remainder, 60)

        self.logger.info("")
        self.logger.info("=" * 60)
        if fail == 0:
            self.logger.info("等位基因鉴定完成|Allelic Gene Identification Completed")
        else:
            self.logger.info(f"等位基因鉴定完成 (部分失败)|Completed with {fail} failures")
        self.logger.info(f"  成功|Success: {success}/{total}")
        self.logger.info(f"  失败|Failed: {fail}/{total}")
        self.logger.info(f"  总耗时|Total time: {hours}h {minutes}m {seconds}s")
        self.logger.info("=" * 60)
