"""插入检测核心算法|Insert detection core algorithm"""

import re
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


# 正则表达式用于解析CIGAR字符串|Regex for parsing CIGAR strings
_SOFT_START = re.compile(r"^(\d+)S")
_SOFT_END = re.compile(r"(\d+)S$")


@dataclass
class InsertionSite:
    """插入位点数据类|Insertion site dataclass"""
    sample_id: str
    chromosome: str
    start: int
    end: int
    orientation: str
    tdna_position: str
    junction: str
    support_reads: int
    score: int
    group_counts: str  # "group1,group2,group3,group4,group5"


class InsertDetector:
    """插入检测器|Insert detector"""

    def __init__(self, config, logger):
        """初始化检测器|Initialize detector

        Args:
            config: InsertDetectionConfig配置对象|InsertDetectionConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger
        self.tdna_size_dict = {}

    def run(self, samples: Dict[str, Dict]) -> List[InsertionSite]:
        """运行检测流程|Run detection pipeline

        Args:
            samples: 样品信息字典|Sample information dictionary
                {'sample_name': {'r1': 'path1', 'r2': 'path2'}}

        Returns:
            List[InsertionSite]: 插入位点列表|List of insertion sites
        """
        all_results = []

        for sample_id, files in samples.items():
            self.logger.info(f"处理样品|Processing sample: {sample_id}")

            # 步骤1: 比对到基因组|Step 1: Align to genome
            bam_file = self._align_to_genome(sample_id, files)

            # 步骤2: 提取soft-clipped序列|Step 2: Extract soft-clipped sequences
            softclip_fq = self._extract_soft_clipped(sample_id, bam_file)

            # 步骤3: 比对soft-clipped序列到插入序列|Step 3: Align soft-clipped to insert
            tdna_sam = self._map_softclip_to_insert(sample_id, softclip_fq)

            # 步骤4: 过滤和评分|Step 4: Filter and score
            sites = self._filter_and_score(sample_id, tdna_sam)

            # 输出候选位点信息|Output candidate site information
            for site in sites:
                self.logger.info(
                    f"  候选位点|Candidate: {site.chromosome}:{site.start}-{site.end} | "
                    f"方向|Orientation: {site.orientation} | "
                    f"T-DNA位置|T-DNA pos: {site.tdna_position} | "
                    f"支持reads|Support reads: {site.support_reads} | "
                    f"得分|Score: {site.score}"
                )

            all_results.extend(sites)

        return all_results

    def _align_to_genome(self, sample_id: str, files: Dict) -> Path:
        """比对reads到基因组|Align reads to genome

        Args:
            sample_id: 样品ID|Sample ID
            files: 文件字典|File dictionary {'r1': 'path1', 'r2': 'path2'}

        Returns:
            Path: BAM文件路径|BAM file path
        """
        # 按样本分组的目录结构|By Sample directory structure
        # 创建项目子目录，避免与输出目录中的其他文件混杂
        sample_dir = Path(self.config.output_dir) / f"{sample_id}_insert_detection"
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 步骤子目录|Step subdirectory
        step_dir = sample_dir / "01_genome_alignment"
        step_dir.mkdir(parents=True, exist_ok=True)

        bam_path = step_dir / f"{sample_id}.sorted.bam"

        # 检查是否已存在|Check if already exists
        if self.config.skip_existing and bam_path.exists() and Path(f"{bam_path}.bai").exists():
            self.logger.info(f"跳过|Skipping: 基因组比对已存在|Genome alignment exists: {sample_id}")
            return bam_path

        self.logger.info(f"比对到基因组|Aligning to genome: {sample_id}")

        # 检查并构建bowtie2索引|Check and build bowtie2 index
        index_prefix = self.config.genome
        if index_prefix.endswith('.fa'):
            index_prefix = index_prefix[:-3]

            # 检查索引文件是否存在|Check if index files exist
            index_files = [
                Path(f"{index_prefix}.1.bt2"),
                Path(f"{index_prefix}.2.bt2"),
                Path(f"{index_prefix}.3.bt2"),
                Path(f"{index_prefix}.4.bt2"),
                Path(f"{index_prefix}.rev.1.bt2"),
                Path(f"{index_prefix}.rev.2.bt2")
            ]

            if not all(f.exists() for f in index_files):
                self.logger.info(f"构建bowtie2索引|Building bowtie2 index: {index_prefix}")
                try:
                    # 不捕获输出，让进度信息能够显示|Don't capture output to show progress
                    subprocess.run(
                        ["bowtie2-build", self.config.genome, index_prefix],
                        check=True
                    )
                    self.logger.info(f"索引构建完成|Index build completed")
                except subprocess.CalledProcessError as e:
                    self.logger.error(f"索引构建失败|Index build failed: {e}")
                    raise

        # 构建bowtie2命令|Build bowtie2 command
        cmd = [
            self.config.bowtie2_path,
            "-x", index_prefix,
            "-U", files['r1'], files['r2'],
            "--very-sensitive-local",
            "-p", str(self.config.threads),
            "--un", str(step_dir / f"{sample_id}.unmapped.fastq"),
            "|", self.config.samtools_path, "sort", "-O", "bam", "-o", str(bam_path)
        ]

        cmd_str = " ".join(cmd)
        self.logger.debug(f"运行|Running: {cmd_str}")
        self.logger.info(f"开始比对（可能需要较长时间）|Starting alignment (may take a while)...")

        try:
            # 不捕获输出，让进度信息能够显示|Don't capture output to show progress
            subprocess.run(cmd_str, shell=True, check=True)

            # 索引BAM文件|Index BAM file
            self.logger.info(f"建立BAM索引|Building BAM index...")
            subprocess.run(
                [self.config.samtools_path, "index", str(bam_path)],
                check=True
            )

            self.logger.info(f"比对完成|Alignment completed: {sample_id}")
            return bam_path

        except subprocess.CalledProcessError as e:
            self.logger.error(f"比对失败|Alignment failed: {e}")
            raise

    def _extract_soft_clipped(self, sample_id: str, bam_file: Path) -> Path:
        """提取soft-clipped reads|Extract soft-clipped reads

        Args:
            sample_id: 样品ID|Sample ID
            bam_file: BAM文件路径|BAM file path

        Returns:
            Path: 输出FASTQ文件路径|Output FASTQ file path
        """
        # 按样本分组的目录结构|By Sample directory structure
        # 创建项目子目录，避免与输出目录中的其他文件混杂
        sample_dir = Path(self.config.output_dir) / f"{sample_id}_insert_detection"
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 步骤子目录|Step subdirectory
        output_dir = sample_dir / "02_softclip_extraction"
        output_dir.mkdir(parents=True, exist_ok=True)

        output_fq = output_dir / f"{sample_id}.softclip.fastq"

        # 检查是否已存在|Check if already exists
        if self.config.skip_existing and output_fq.exists() and output_fq.stat().st_size > 0:
            self.logger.info(f"跳过|Skipping: soft-clipped序列已存在|Soft-clipped sequences exist: {sample_id}")
            return output_fq

        self.logger.info(f"提取soft-clipped序列|Extracting soft-clipped sequences: {sample_id}")

        # 使用临时目录处理SAM文件|Use temp directory for SAM processing
        # 临时目录落到 <output>/tmp 下,避免系统 /tmp 爆满|
        # Temp dir under <output>/tmp to avoid system /tmp overflow
        import tempfile
        tmp_root = self.config.output_path / "tmp"
        tmp_root.mkdir(parents=True, exist_ok=True)
        tmp_dir = Path(tempfile.mkdtemp(prefix=f"softclip_{sample_id}_", dir=str(tmp_root)))

        # 提取soft-clipped reads到SAM|Extract soft-clipped reads to SAM
        sam_file = tmp_dir / f"{sample_id}.softclip.sam"
        cmd1 = f"{self.config.samtools_path} view {bam_file} | awk '$6 ~ /S/' > {sam_file}"

        # 从SAM提取soft-clipped序列到FASTQ|Extract soft-clipped sequences to FASTQ
        # 预先获取min_clip值以避免格式化问题|Pre-get min_clip value to avoid formatting issues
        min_clip_val = self.config.min_clip
        cmd2 = f"""
        {self.config.samtools_path} view {bam_file} | \
        awk '$6 ~ /S/' | \
        python3 -c "
import sys
import re

_SOFT_START = re.compile(r'^(\\\d+)S')
_SOFT_END = re.compile(r'(\d+)S\$')

for line in sys.stdin:
    if line.startswith('@'):
        continue
    cols = line.rstrip().split('\\\t')
    if len(cols) < 11:
        continue

    readid, flag, chrom, pos, _, cigar, _, _, _, seq_full, qual_full = \
        cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10]

    for match, side in [(_SOFT_START.match(cigar), 'L'), (_SOFT_END.search(cigar), 'R')]:
        if not match:
            continue
        length = int(match.group(1))
        if length < {min_clip_val}:
            continue

        seq = seq_full[:length] if side == 'L' else seq_full[-length:]
        remain = seq_full[length:] if side == 'L' else seq_full[:-length]

        # 检查剩余序列质量|Check remaining sequence quality
        if max(remain.count(b) for b in ('A','T','C','G')) > 0.8 * len(remain):
            continue

        qual = qual_full[:length] if side == 'L' else qual_full[-length:]
        header = f'@{{readid}}::{{chrom}}::{{pos}}::{{flag}}::{{cigar}}::{{side}}'
        print(f'{{header}}\\\n{{seq}}\\\n+{{header}}\\\n{{qual}}')
" > {output_fq}
        """

        try:
            # 使用samtools提取soft-clipped reads|Extract soft-clipped reads using samtools
            # 不捕获输出，让进度信息能够显示|Don't capture output to show progress
            subprocess.run(cmd1, shell=True, check=True)

            # 提取soft-clipped序列|Extract soft-clipped sequences
            self._extract_soft_clipped_from_sam(sam_file, output_fq)

            # 将soft-clipped序列转换为BAM并索引，供IGV查看|Convert soft-clipped to BAM for IGV
            output_bam = output_dir / f"{sample_id}.softclip.bam"
            self.logger.info(f"转换soft-clipped BAM用于IGV查看|Converting soft-clipped BAM for IGV...")
            subprocess.run(
                f"{self.config.samtools_path} view -h {bam_file} -S {sam_file} -b | {self.config.samtools_path} sort -o {output_bam} -",
                shell=True,
                check=True
            )
            subprocess.run(
                [self.config.samtools_path, "index", str(output_bam)],
                check=True
            )

            self.logger.info(f"soft-clipped序列提取完成|Soft-clipped extraction completed: {sample_id}")
            return output_fq

        except subprocess.CalledProcessError as e:
            self.logger.error(f"提取soft-clipped序列失败|Soft-clipped extraction failed: {e}")
            raise
        finally:
            # 清理临时文件|Clean up temp files
            import shutil
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)

    def _extract_soft_clipped_from_sam(self, sam_file: Path, output_fq: Path):
        """从SAM文件提取soft-clipped序列|Extract soft-clipped sequences from SAM file"""
        with open(sam_file) as sam, open(output_fq, 'w') as fq:
            for line in sam:
                if line.startswith('@'):
                    continue

                cols = line.rstrip().split('\t')
                if len(cols) < 11:
                    continue

                readid, flag, chrom, pos, _, cigar, _, _, _, seq_full, qual_full = \
                    cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10]

                for match, side in [(_SOFT_START.match(cigar), 'L'), (_SOFT_END.search(cigar), 'R')]:
                    if not match:
                        continue

                    length = int(match.group(1))
                    if length < self.config.min_clip:
                        continue

                    seq = seq_full[:length] if side == 'L' else seq_full[-length:]
                    remain = seq_full[length:] if side == 'L' else seq_full[:-length]

                    # 检查剩余序列质量|Check remaining sequence quality
                    if max(remain.count(b) for b in ('A', 'T', 'C', 'G')) > 0.8 * len(remain):
                        continue

                    qual = qual_full[:length] if side == 'L' else qual_full[-length:]
                    header = f"@{readid}::{chrom}::{pos}::{flag}::{cigar}::{side}"
                    fq.write(f"{header}\n{seq}\n+{header}\n{qual}\n")

    def _map_softclip_to_insert(self, sample_id: str, softclip_fq: Path) -> Path:
        """比对soft-clipped序列到插入序列|Map soft-clipped sequences to insert

        Args:
            sample_id: 样品ID|Sample ID
            softclip_fq: soft-clipped FASTQ文件|soft-clipped FASTQ file

        Returns:
            Path: 输出SAM文件路径|Output SAM file path
        """
        # 按样本分组的目录结构|By Sample directory structure
        # 创建项目子目录，避免与输出目录中的其他文件混杂
        sample_dir = Path(self.config.output_dir) / f"{sample_id}_insert_detection"
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 步骤子目录|Step subdirectory
        output_dir = sample_dir / "03_insert_alignment"
        output_dir.mkdir(parents=True, exist_ok=True)

        # 使用临时目录处理SAM文件|Use temp directory for SAM processing
        # 临时目录落到 <output>/tmp 下,避免系统 /tmp 爆满|
        # Temp dir under <output>/tmp to avoid system /tmp overflow
        import tempfile
        tmp_root = self.config.output_path / "tmp"
        tmp_root.mkdir(parents=True, exist_ok=True)
        tmp_dir = Path(tempfile.mkdtemp(prefix=f"tdna_{sample_id}_", dir=str(tmp_root)))
        output_sam = tmp_dir / f"{sample_id}.tdna.sam"

        # 检查是否已存在|Check if already exists
        if self.config.skip_existing and output_sam.exists() and output_sam.stat().st_size > 0:
            self.logger.info(f"跳过|Skipping: 插入序列比对已存在|Insert alignment exists: {sample_id}")
            return output_sam

        # 检查并构建bowtie2索引|Check and build bowtie2 index
        index_prefix = self.config.insert_sequence
        if index_prefix.endswith('.fa'):
            index_prefix = index_prefix[:-3]
        elif index_prefix.endswith('.fasta'):
            index_prefix = index_prefix[:-5]
        elif index_prefix.endswith('.fna'):
            index_prefix = index_prefix[:-4]

        index_files = [
            Path(f"{index_prefix}.1.bt2"),
            Path(f"{index_prefix}.2.bt2"),
            Path(f"{index_prefix}.3.bt2"),
            Path(f"{index_prefix}.4.bt2"),
            Path(f"{index_prefix}.rev.1.bt2"),
            Path(f"{index_prefix}.rev.2.bt2")
        ]

        if not all(f.exists() for f in index_files):
            self.logger.info(f"构建插入序列bowtie2索引|Building insert sequence bowtie2 index: {index_prefix}")
            try:
                subprocess.run(
                    ["bowtie2-build", self.config.insert_sequence, index_prefix],
                    check=True
                )
                self.logger.info(f"索引构建完成|Index build completed")
            except subprocess.CalledProcessError as e:
                self.logger.error(f"索引构建失败|Index build failed: {e}")
                raise

        self.logger.info(f"比对到插入序列|Aligning to insert sequence: {sample_id}")

        # 第一步：使用--local模式比对，同时生成SAM和BAM|Step 1: Align with --local mode, generate both SAM and BAM
        output_bam = output_dir / f"{sample_id}.tdna.bam"
        cmd1 = f"{self.config.bowtie2_path} -x {index_prefix} -U {softclip_fq} -L 10 --local -p {self.config.threads} | awk '$3!=\"*\"' | tee {output_sam} | {self.config.samtools_path} view -S - -b | {self.config.samtools_path} sort -o {output_bam} -"

        try:
            # 不捕获输出，让进度信息能够显示|Don't capture output to show progress
            subprocess.run(cmd1, shell=True, check=True)

            # 索引BAM文件|Index BAM file for IGV
            self.logger.info(f"建立T-DNA BAM索引|Building T-DNA BAM index for IGV...")
            subprocess.run(
                [self.config.samtools_path, "index", str(output_bam)],
                check=True
            )

            self.logger.info(f"插入序列比对完成|Insert alignment completed: {sample_id}")

            # 在返回之前读取T-DNA大小信息到内存|Read T-DNA size info into memory before returning
            self._read_tdna_sizes(output_sam)

            # 返回BAM文件路径而不是临时SAM路径|Return BAM path instead of temp SAM path
            return output_bam

        except subprocess.CalledProcessError as e:
            self.logger.error(f"插入序列比对失败|Insert alignment failed: {e}")
            raise
        finally:
            # 清理临时SAM文件|Clean up temp SAM file
            import shutil
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)

    def _get_matched_len(self, md: str) -> int:
        """从MD标签计算匹配长度|Calculate matched length from MD tag"""
        if md.startswith("MD:Z:"):
            md = md.split(':', 2)[2]
        nums = re.findall(r'\d+', md)
        return sum(int(n) for n in nums)

    def _get_mapped_position(self, start: int, cigar: str) -> int:
        """根据CIGAR计算映射结束位置|Calculate mapped end position from CIGAR"""
        cigar_clean = re.sub(r'\d+S', '', cigar)
        cigar_clean = cigar_clean.replace('N', 'M').replace('D', 'M')
        total_m = sum(int(n) for n in re.findall(r'(\d+)M', cigar_clean))
        return start + total_m - 1

    def _filter_and_score(self, sample_id: str, tdna_bam: Path) -> List[InsertionSite]:
        """过滤和评分|Filter and score insertion sites

        Args:
            sample_id: 样品ID|Sample ID
            tdna_bam: T-DNA比对BAM文件|T-DNA alignment BAM file

        Returns:
            List[InsertionSite]: 插入位点列表|List of insertion sites
        """
        self.logger.info(f"过滤和评分|Filtering and scoring: {sample_id}")

        # 读取T-DNA大小（已在_map_softclip_to_insert中读取）|Read T-DNA sizes (already read in _map_softclip_to_insert)
        # 不需要再读取|No need to read again

        # 使用samtools读取比对结果|Use samtools to read alignments
        import subprocess
        try:
            result = subprocess.run(
                [self.config.samtools_path, "view", str(tdna_bam)],
                capture_output=True, text=True, check=True
            )
            lines = result.stdout.split('\n')
        except subprocess.CalledProcessError as e:
            self.logger.error(f"读取BAM文件失败|Failed to read BAM file: {e}")
            return []

        filtered_reads = []
        for line in lines:
            if not line or line.startswith('@'):
                continue

            cols = line.rstrip().split('\t')
            if len(cols) < 12:
                continue

            tdna = cols[2]
            tdna_pos = int(cols[3])
            tdna_cigar = cols[5]
            tdna_seq = cols[9]

            # 解析genomic信息（从QNAME中）|Parse genomic info (from QNAME)
            genomic_info = cols[0]
            s = genomic_info.split('::')
            if len(s) < 6:
                continue
            chrom, pos, flag, cigar, clip_pos = s[1], s[2], s[3], s[4], s[5]

            # 计算错配数（从NM标签）|Count mismatches from NM tag
            try:
                # 查找NM标签|Find NM tag
                nm_tag = None
                for col in cols[11:]:
                    if col.startswith('NM:i:'):
                        nm_tag = col.split(':')[2]
                        break
                if nm_tag is None:
                    continue
                mismatch = int(nm_tag)
            except (ValueError, IndexError):
                continue

            if mismatch > 1:  # 最多允许1个错配|Max 1 mismatch
                continue

            # 获取MD标签|Get MD tag
            try:
                md_tag = None
                for col in cols[11:]:
                    if col.startswith('MD:Z:'):
                        md_tag = col
                        break
                if md_tag is None:
                    continue
            except (ValueError, IndexError):
                continue

            # 计算匹配长度|Calculate matched length
            matched_len = self._get_matched_len(md_tag)

            # 判断是否靠近边界|Check if near border
            tdna_size = self.tdna_size_dict.get(tdna, 0)
            near_border = tdna_pos < self.config.border_window or \
                         tdna_pos > tdna_size - self.config.border_window

            # 应用不同的阈值|Apply different thresholds
            if near_border and matched_len >= 18:  # 靠近边界用l1=18
                pos_out = str(tdna_pos)
            elif not near_border and matched_len >= 30:  # 远离边界用l2=30
                pos_out = str(tdna_pos)
            else:
                continue

            filtered_reads.append({
                'chrom': chrom,
                'start': int(pos),
                'flag': flag,
                'cigar': cigar,
                'clip_pos': clip_pos,
                'tdna': tdna,
                'tdna_pos': pos_out,
                'tdna_cigar': tdna_cigar,
                'matched_len': matched_len,
                'tdna_seq': tdna_seq
            })

        if not filtered_reads:
            self.logger.warning(f"没有找到有效reads|No valid reads found: {sample_id}")
            return []

        # 聚类和评分|Cluster and score
        sites = self._cluster_and_score(sample_id, filtered_reads)

        return sites

    def _read_tdna_sizes(self, tdna_sam: Path):
        """从SAM头读取T-DNA大小|Read T-DNA sizes from SAM header"""
        with open(tdna_sam) as f:
            for line in f:
                if line.startswith('@SQ'):
                    parts = line.rstrip().split('\t')
                    for part in parts:
                        if part.startswith('SN:'):
                            tdna = part.split(':')[1]
                        elif part.startswith('LN:'):
                            size = int(part.split(':')[1])
                            self.tdna_size_dict[tdna] = size

    def _cluster_and_score(self, sample_id: str, reads: List[Dict]) -> List[InsertionSite]:
        """聚类和评分候选位点|Cluster and score candidate sites

        基于TDNAreader的5组reads分类系统|Based on TDNAreader's 5-group classification
        """
        # 权重系统|Weight system: (700, 300, 300, 1, 1)
        weights = (700, 300, 300, 1, 1)

        # T-DNA位置字典和评分字典|T-DNA position dict and score dict
        tdna_pos_dict = {}
        tdna_score_dict = {}
        tis_dict = {}
        tis_cigar_dict = {}
        tis_raw_dict = {}

        # 第一遍：初始化TIS|First pass: Initialize TISs
        for read in reads:
            chrom, start = read['chrom'], read['start']
            cigar, clip_pos = read['cigar'], read['clip_pos']
            tdna, tdna_pos = read['tdna'], read['tdna_pos']

            # 计算末端位置|Calculate end position
            end = self._get_mapped_position(start, cigar)

            # 确定方向|Determine orientation
            ori = '+' if read['flag'] == '0' else '-'

            # 创建TIS ID|Create TIS ID
            if clip_pos == "R":
                tis = f'{chrom}::{end}::{clip_pos}::{ori}'
                cigar_len = int(cigar.split('M')[-1].split('S')[0])
            else:
                tis = f'{chrom}::{start}::{clip_pos}::{ori}'
                cigar_len = int(cigar.split('S')[0])

            tdna_junc = f'{tdna}::{tdna_pos}'

            # 检查是否已存在|Check if already exists
            matched = False
            adjacent = False  # 添加adjacent标志，与TDNAreader保持一致
            for key, members in tis_dict.items():
                base_tdna = tdna_pos_dict[key].split('::')[0]
                base_pos = tdna_pos_dict[key].split('::')[1]

                # 必须是同一个T-DNA|Must be same T-DNA
                if tdna != base_tdna:
                    continue

                # 检查T-DNA位置距离|Check T-DNA position distance
                if self._calculate_distance(tdna_pos, base_pos) >= 150:
                    continue

                # 完全匹配|Exact match
                if tis in members:
                    original_cigar = tis_cigar_dict[tis]
                    diff = abs(cigar_len - original_cigar)
                    if diff == 0:
                        tdna_score_dict[key][2] += 1  # Group 3
                    elif diff == 1:
                        tdna_score_dict[key][1] += 1  # Group 2
                    else:
                        tdna_score_dict[key][0] += 1  # Group 1
                    matched = True
                    tis_raw_dict.setdefault(tis, []).append(read)
                    break  # 完全匹配后立即break

                # 相邻位点|Adjacent sites
                genomic_dist = abs(int(tis.split('::')[1]) - int(key.split('::')[1]))
                if chrom == key.split('::')[0] and \
                   clip_pos == key.split('::')[2] and \
                   genomic_dist < self.config.border_window:
                    members.append(tis)
                    tis_cigar_dict[tis] = cigar_len
                    tdna_score_dict[key][3] += 1  # Group 4
                    adjacent = True  # 只设置adjacent，不设置matched
                    tis_raw_dict.setdefault(tis, []).append(read)

            if not matched and not adjacent:
                # 创建新TIS|Create new TIS
                tis_dict[tis] = [tis]
                tdna_pos_dict[tis] = tdna_junc
                tdna_score_dict[tis] = np.zeros(5, dtype=int)  # 5个组
                tis_cigar_dict[tis] = cigar_len
                tis_raw_dict.setdefault(tis, []).append(read)

        # 计算得分并过滤|Calculate scores and filter
        results = []
        for key, members in tis_dict.items():
            chrom = key.split('::')[0]
            start = min([int(i.split('::')[1]) for i in members])
            end = max([int(i.split('::')[1]) for i in members])
            clip_pos = key.split('::')[2]
            orientation = key.split('::')[3]

            # 计算加权得分|Calculate weighted score
            scores = tdna_score_dict[key]
            total_score = np.dot(scores, weights)

            if total_score < self.config.score_threshold:
                continue

            # 计算总支持reads数|Calculate total supporting reads
            total_supports = int(np.sum(scores))

            # 确定junction类型|Determine junction type
            junction = 'T-DNA:REF' if clip_pos == 'L' else 'REF:T-DNA'

            site = InsertionSite(
                sample_id=sample_id,
                chromosome=chrom,
                start=start,
                end=end if end > start else start + 1,
                orientation=orientation,
                tdna_position=tdna_pos_dict[key],
                junction=junction,
                support_reads=total_supports,
                score=int(total_score),
                group_counts=','.join(map(str, scores))
            )
            results.append(site)

        # 按得分排序|Sort by score
        results.sort(key=lambda x: -x.score)

        self.logger.info(f"找到{len(results)}个候选位点|Found {len(results)} candidate sites: {sample_id}")

        return results

    def _calculate_distance(self, pos1: str, pos2: str) -> int:
        """计算两个T-DNA位置之间的距离|Calculate distance between two T-DNA positions"""
        def parse(p: str):
            if "(" in p:
                a, b = p.split("(", 1)
                return [int(a), int(b.rstrip(")"))]
            return [int(p)]

        vals1, vals2 = parse(pos1), parse(pos2)
        return min(abs(a - b) for a in vals1 for b in vals2)
