"""
MCscan共线性分析运行器|MCscan Collinearity Analysis Runner

基于JcviPipeline共享管道, 添加共线性统计汇总功能
Based on JcviPipeline shared pipeline, adds collinearity summary
"""

import time
from pathlib import Path
from datetime import datetime
from typing import Optional, List

from ..config import JcviBaseConfig
from ..pipeline import JcviPipeline
from ..utils import JcviLogger


class McscanRunner:
    """MCscan共线性分析运行器|MCscan Collinearity Analysis Runner"""

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
                self.logger_obj = JcviLogger(log_dir / "mcscan.log", "MCscan")
                self.logger = self.logger_obj.get_logger()

            self._print_header()
            self.config.validate()

            pipeline = JcviPipeline(self.config, self.logger_obj)
            result = pipeline.run()

            # 汇总共线性结果
            self._log_section("汇总结果|Summarize results")
            self._summarize_collinearity(result.pair_dirs, result.bed_dir)

            self._print_footer(len(result.pair_dirs),
                              len(result.pair_list) - len(result.pair_dirs),
                              len(result.pair_list))
            return len(result.pair_dirs) > 0

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Program execution error: {e}", exc_info=True)
            return False

    def _summarize_collinearity(self, pair_dirs: dict, bed_dir: Path):
        """汇总所有两两比较的共线性统计|Summarize pairwise collinearity stats"""
        summary_file = Path(self.config.output_dir) / "all_collinearity_summary.txt"

        total_blocks = 0
        total_pairs = 0

        with open(summary_file, 'w') as fout:
            fout.write("#sample_a\tsample_b\tblocks\tgene_pairs\tavg_block_size\tmax_block_size\n")

            for (name_a, name_b), pair_dir in sorted(pair_dirs.items()):
                anchors_file = pair_dir / f"{name_a}.{name_b}.anchors"
                stats = self._parse_anchors(str(anchors_file))

                if stats:
                    fout.write(f"{name_a}\t{name_b}\t"
                               f"{stats['blocks']}\t{stats['pairs']}\t"
                               f"{stats['avg_size']:.1f}\t{stats['max_size']}\n")
                    total_blocks += stats['blocks']
                    total_pairs += stats['pairs']
                    self.logger.info(f"  {name_a} vs {name_b}: "
                                    f"{stats['blocks']} blocks, "
                                    f"{stats['pairs']} gene pairs")
                else:
                    fout.write(f"{name_a}\t{name_b}\t0\t0\t0.0\t0\n")
                    self.logger.warning(f"  {name_a} vs {name_b}: 无共线性结果|No collinearity")

            # 总计行
            fout.write(f"#TOTAL\t\t{total_blocks}\t{total_pairs}\t-\t-\n")

        self.logger.info(f"\n  汇总文件|Summary file: {summary_file}")
        self.logger.info(f"  总共线性区块|Total blocks: {total_blocks}")
        self.logger.info(f"  总基因对|Total gene pairs: {total_pairs}")

    @staticmethod
    def _parse_anchors(anchors_file: str) -> Optional[dict]:
        """解析.anchors文件统计block和pair信息"""
        blocks = 0
        pairs_in_block = []
        current_block_size = 0

        try:
            with open(anchors_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("###"):
                        if current_block_size > 0:
                            pairs_in_block.append(current_block_size)
                        blocks += 1
                        current_block_size = 0
                    elif line and not line.startswith("#"):
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            current_block_size += 1
            # 最后一个block
            if current_block_size > 0:
                pairs_in_block.append(current_block_size)

            if blocks == 0:
                return None

            return {
                'blocks': blocks,
                'pairs': sum(pairs_in_block),
                'avg_size': sum(pairs_in_block) / len(pairs_in_block) if pairs_in_block else 0,
                'max_size': max(pairs_in_block) if pairs_in_block else 0,
            }
        except FileNotFoundError:
            return None

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
        self.logger.info("JCVI MCscan 共线性分析|JCVI MCscan Collinearity Analysis")
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
            self.logger.info("MCscan共线性分析完成|MCscan Collinearity Analysis Completed")
        else:
            self.logger.info(f"MCscan共线性分析完成 (部分失败)|Completed with {fail} failures")
        self.logger.info(f"  成功|Success: {success}/{total}")
        self.logger.info(f"  失败|Failed: {fail}/{total}")
        self.logger.info(f"  总耗时|Total time: {hours}h {minutes}m {seconds}s")
        self.logger.info("=" * 60)
