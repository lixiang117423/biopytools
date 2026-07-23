"""泛基因组Block构建 - 核心算法模块|Pan-Blocks Construction - Core Algorithm Module"""

import os
import shutil
import tempfile
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Tuple, Dict
import logging

from .config import PanBlocksConfig
from .utils import (
    parse_coords_file, build_conda_command, write_bed_file,
    read_bed_file, CommandRunner, format_number
)


class PanBlockBuilder:
    """泛基因组Block构建器|Pan-Genome Block Builder"""

    def __init__(self, config: PanBlocksConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger)
        self.coords_dir = Path(config.coords_dir)
        self.blocks_dir = Path(config.blocks_dir)
        # 临时目录落到 <output>/tmp 下,避免系统 /tmp 爆满|
        # Temp dir under <output>/tmp to avoid system /tmp overflow
        tmp_root = self.config.output_path / "tmp"
        tmp_root.mkdir(parents=True, exist_ok=True)
        self.tmp_dir = tempfile.mkdtemp(prefix="pan_blocks_", dir=str(tmp_root))

    def build_all_chromosomes(self):
        """构建所有染色体的Pan-Blocks|Build pan-blocks for all chromosomes"""
        try:
            return self._build_all_chromosomes_impl()
        finally:
            # 清理临时目录,避免 <output>/tmp/pan_blocks_* 残留|
            # Clean up temp dir to avoid leftover <output>/tmp/pan_blocks_*
            # 只删 mkdtemp 创建的子目录,不动 tmp_root 共享目录|
            # Only remove the mkdtemp-created subdir, never the shared tmp_root
            shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _build_all_chromosomes_impl(self):
        """build_all_chromosomes 主体|build_all_chromosomes body"""
        chromosomes = self.config.get_target_chromosomes()
        self.logger.info(f"开始构建Pan-Blocks|Starting pan-block construction: {len(chromosomes)} 条染色体|chromosomes")

        success_count = 0
        fail_count = 0

        with ThreadPoolExecutor(max_workers=self.config.parallel_alignments) as executor:
            futures = {
                executor.submit(self._run_with_checkpoint, chrom): chrom
                for chrom in chromosomes
            }
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    success = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    self.logger.error(f"染色体构建异常|Chromosome build exception: {chrom}: {e}")
                    fail_count += 1

        self.logger.info(f"Pan-Blocks构建完成|Pan-blocks construction completed: "
                         f"成功|success={success_count}, 失败|fail={fail_count}")
        return fail_count == 0

    def _run_with_checkpoint(self, chrom: str) -> bool:
        """带断点续传的染色体构建|Chromosome build with checkpoint resume"""
        output_file = self.blocks_dir / f"{chrom}.pan_blocks.bed"
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.info(f"跳过已完成染色体|Skipping completed chromosome: {chrom}")
            return True
        return self.build_chromosome(chrom)

    def build_chromosome(self, chrom: str) -> bool:
        """构建单条染色体的Pan-Blocks|Build pan-blocks for a single chromosome

        算法|Algorithm:
        1. 收集所有两两比对中的共线性区域（按来源基因组分组合并）
        2. 按来源基因组对做 bedtools merge
        3. 按 genome_order 迭代去重（bedtools subtract）
        4. 输出 Pan-Blocks
        """
        genome_order = self.config.genome_order_list
        self.logger.info(f"构建染色体|Building chromosome: {chrom} ({len(genome_order)} 个基因组|genomes)")

        # Phase 1: 收集共线性区域
        genome_regions = self._collect_regions(chrom, genome_order)
        total_regions = sum(sum(len(r) for r in partners.values()) for partners in genome_regions.values())
        self.logger.info(f"收集到共线性区域|Collected syntenic regions: {format_number(total_regions)}")

        # Phase 2: 按来源基因组合并重叠区域 (bedtools merge)
        merged_regions = self._merge_regions_per_partner(genome_regions)
        self.logger.info(f"合并重叠区域完成|Region merging completed")

        # Phase 3: 迭代去重 (bedtools subtract)
        filtered_results = self._iterative_subtraction(chrom, genome_order, merged_regions)
        self.logger.info(f"迭代去重完成|Iterative subtraction completed")

        # Phase 4: 合并输出
        pan_blocks_file = self.blocks_dir / f"{chrom}.pan_blocks.bed"
        self._write_pan_blocks(pan_blocks_file, filtered_results, genome_order)
        block_count = sum(len(v) for v in filtered_results.values())
        self.logger.info(f"染色体 {chrom} 完成|Chromosome {chrom} done: {format_number(block_count)} Pan-Blocks")

        return True

    def _collect_regions(self, chrom: str, genome_order: List[str]) -> Dict:
        """收集所有两两比对中的共线性区域|Collect syntenic regions from all pairwise alignments

        从每个 ref.vs.query.filtered.coords 中提取同染色体比对，
        记录 ref 和 query 上对应的区域，标记 contributor。

        Returns:
            genome_regions[genome_name][partner_name] = [(chr, bed_start, bed_end), ...]
            坐标已转为 0-based half-open BED 格式
        """
        genome_regions = {g: {} for g in genome_order}

        for ref_genome in genome_order:
            for query_genome in genome_order:
                if ref_genome == query_genome:
                    continue

                coords_file = self.coords_dir / f"{ref_genome}.vs.{query_genome}.filtered.coords"
                if not coords_file.exists():
                    continue

                alignments = parse_coords_file(str(coords_file))

                for aln in alignments:
                    if aln['ref_chr'] != chrom or aln['qry_chr'] != chrom:
                        continue

                    # 确保 start < end，转 0-based BED
                    ref_start = min(aln['ref_start'], aln['ref_end']) - 1
                    ref_end = max(aln['ref_start'], aln['ref_end'])
                    qry_start = min(aln['qry_start'], aln['qry_end']) - 1
                    qry_end = max(aln['qry_start'], aln['qry_end'])

                    # ref 获得区域，contributor = query
                    if query_genome not in genome_regions[ref_genome]:
                        genome_regions[ref_genome][query_genome] = []
                    genome_regions[ref_genome][query_genome].append((chrom, ref_start, ref_end))

                    # query 获得区域，contributor = ref
                    if ref_genome not in genome_regions[query_genome]:
                        genome_regions[query_genome][ref_genome] = []
                    genome_regions[query_genome][ref_genome].append((chrom, qry_start, qry_end))

        return genome_regions

    def _merge_regions_per_partner(self, genome_regions: Dict) -> Dict:
        """按来源基因组合并重叠区域|Merge overlapping regions per partner using bedtools merge

        Args:
            genome_regions[genome][partner] = [(chr, start, end), ...]

        Returns:
            merged[genome][partner] = [(chr, start, end), ...]
        """
        merged = {}
        for genome, partners in genome_regions.items():
            merged[genome] = {}
            for partner, regions in partners.items():
                if not regions:
                    continue

                # 写入 BED，运行 bedtools merge
                input_bed = os.path.join(self.tmp_dir, f"merge_{genome}_{partner}.bed")
                output_bed = os.path.join(self.tmp_dir, f"merged_{genome}_{partner}.bed")

                write_bed_file(regions, input_bed)

                cmd = [self.config.bedtools_path, 'merge', '-i', input_bed]
                success = self._run_bedtools_to_file(cmd, output_bed)
                if success:
                    merged[genome][partner] = read_bed_file(output_bed)
                else:
                    self.logger.warning(f"bedtools merge失败|bedtools merge failed: {genome} - {partner}")
                    merged[genome][partner] = regions

        return merged

    def _iterative_subtraction(self, chrom: str, genome_order: List[str],
                               merged_regions: Dict) -> Dict[str, List[Tuple]]:
        """迭代去重：按 genome_order 顺序，逐个基因组减去已分配区域|Iterative subtraction

        Args:
            genome_order: 优先级从高到低的基因组列表
            merged_regions[genome][partner] = [(chr, start, end), ...]

        Returns:
            filtered_results[genome] = [(chr, start, end), ...]  # 该基因组贡献的 Pan-Blocks
        """
        filtered_results = {}
        accumulated_bed = os.path.join(self.tmp_dir, f"accumulated_{chrom}.bed")

        # 创建空 accumulated 文件
        Path(accumulated_bed).touch()

        for m, genome in enumerate(genome_order):
            # 收集该基因组所有来源的 merged regions
            all_regions = []
            for partner, regions in merged_regions.get(genome, {}).items():
                all_regions.extend(regions)

            if not all_regions:
                filtered_results[genome] = []
                continue

            # 写入该基因组所有区域的 BED
            genome_all_bed = os.path.join(self.tmp_dir, f"{genome}_{chrom}_all.bed")
            write_bed_file(all_regions, genome_all_bed)

            if m == 0:
                # 第一个基因组：保留全部
                filtered_bed = os.path.join(self.tmp_dir, f"{genome}_{chrom}.filtered")
                cmd = ['sort', '-k1,1', '-k2,2n', '-k3,3n', genome_all_bed]
                self._run_shell_to_file(cmd, filtered_bed)
            else:
                # 后续基因组：减去 accumulated
                filtered_bed = os.path.join(self.tmp_dir, f"{genome}_{chrom}.filtered")
                cmd = [self.config.bedtools_path, 'subtract', '-a', genome_all_bed, '-b', accumulated_bed]
                self._run_bedtools_to_file(cmd, filtered_bed)

            # 读取过滤后的区域
            filtered = read_bed_file(filtered_bed)
            filtered_results[genome] = filtered

            # 将过滤后的区域追加到 accumulated
            self._append_to_file(filtered_bed, accumulated_bed)

        return filtered_results

    def _write_pan_blocks(self, output_file: Path, filtered_results: Dict,
                          genome_order: List[str]):
        """合并所有基因组的过滤结果为 Pan-Blocks 文件|Merge all filtered results into pan-blocks file"""
        with open(output_file, 'w') as out:
            out.write("Chr\tStart\tEnd\tGenome\n")
            for genome in genome_order:
                regions = filtered_results.get(genome, [])
                for chr_name, start, end in regions:
                    # 输出 1-based inclusive 坐标（与原文一致）
                    out.write(f"{chr_name}\t{start + 1}\t{end}\t{genome}\n")

    def _run_bedtools_to_file(self, cmd: List[str], output_file: str) -> bool:
        """运行 bedtools 命令并将 stdout 写入文件|Run bedtools command and write stdout to file"""
        import subprocess
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, shell=False)
            if result.returncode != 0:
                self.logger.error(f"bedtools命令失败|bedtools command failed: {result.stderr.decode('utf-8', errors='ignore')}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"bedtools执行异常|bedtools execution error: {e}")
            return False

    def _run_shell_to_file(self, cmd: List[str], output_file: str) -> bool:
        """运行 shell 命令并将 stdout 写入文件|Run shell command and write stdout to file"""
        import subprocess
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, shell=False)
            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {result.stderr.decode('utf-8', errors='ignore')}")
                return False
            return True
        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {e}")
            return False

    def _append_to_file(self, source: str, target: str):
        """追加文件内容|Append file content"""
        with open(source) as fin, open(target, 'a') as fout:
            for line in fin:
                fout.write(line)
