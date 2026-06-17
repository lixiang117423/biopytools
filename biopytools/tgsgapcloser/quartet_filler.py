"""
quarTeT GapFiller封装模块|quarTeT GapFiller Wrapper Module

基于quarTeT gapfiller的Gap填充实现，提取核心逻辑并集成到biopytools
Gap filling implementation based on quarTeT gapfiller, integrated into biopytools
"""

import os
import re
import sys
import subprocess
import logging
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import tempfile
import shutil


class QuartetGapFiller:
    """quarTeT风格的Gap填充工具|quarTeT-style Gap Filling Tool"""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """初始化|Initialize"""
        self.logger = logger or logging.getLogger(__name__)

    @staticmethod
    def read_fasta(fasta_file: str) -> Dict[str, str]:
        """读取FASTA文件|Read FASTA file as dictionary"""
        fasta_dict = {}
        with open(fasta_file, 'r') as f:
            all_lines = f.read()
        for id_seq in all_lines.split('>'):
            if not id_seq:
                continue
            try:
                sid_raw, seq_raw = id_seq.split('\n', 1)
                sid = sid_raw.split()[0].strip()
                seq = seq_raw.replace('\n', '').upper()
                fasta_dict[sid] = seq
            except ValueError:
                continue
        return fasta_dict

    @staticmethod
    def read_gfa(gfa_file: str) -> Dict[str, str]:
        """读取GFA文件并提取序列|Read GFA file and extract sequences"""
        fasta_dict = {}
        with open(gfa_file, 'r') as f:
            for line in f:
                if line.startswith('S\t'):  # Sequence line
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        sid = parts[1]
                        seq = parts[2].upper()
                        fasta_dict[sid] = seq
        return fasta_dict

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """生成反向互补序列|Generate reverse complement sequence"""
        trans_table = str.maketrans('ATCG', 'TAGC')
        return seq[::-1].translate(trans_table)

    def count_gaps(self, fasta_file: str, min_gap_length: int = 100) -> Tuple[int, Dict[str, List[Tuple[int, int]]]]:
        """
        统计FASTA文件中的gap数量|Count gaps in FASTA file

        返回|Returns:
            (总gap数量, {序列ID: [(起始位置, 结束位置), ...]})
        """
        fasta_dict = self.read_fasta(fasta_file)
        total_gaps = 0
        gap_dict = {}

        for sid, seq in fasta_dict.items():
            gaps = [(m.start(), m.end()) for m in re.finditer(r'N{' + str(min_gap_length) + ',}', seq)]
            if gaps:
                gap_dict[sid] = gaps
                total_gaps += len(gaps)

        return total_gaps, gap_dict

    def has_gaps(self, fasta_file: str, min_gap_length: int = 100) -> bool:
        """检查是否还有gap|Check if file still has gaps"""
        total_gaps, _ = self.count_gaps(fasta_file, min_gap_length)
        return total_gaps > 0

    def extract_flanking_seqs(self, draft_genome: str, flanking_len: int = 5000,
                              min_gap_length: int = 100) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        提取gap两侧的flanking序列|Extract flanking sequences around gaps

        返回|Returns:
            (flanking_dict, gap_dict)
        """
        draft_dict = self.read_fasta(draft_genome)
        flanking_dict = {}
        gap_dict = {}

        for sid, seq in draft_dict.items():
            if 'N' * min_gap_length not in seq:
                continue

            gap_idx = 1
            gap_sites = [m.span() for m in re.finditer(r'N+', seq)]

            for gap_site in gap_sites:
                start = max(gap_site[0] - flanking_len, 0)
                end = min(gap_site[1] + flanking_len, len(seq))

                left_seq = seq[start:gap_site[0]]
                right_seq = seq[gap_site[1]:end]

                # 检查flanking序列本身是否包含gap
                if 'N' * min_gap_length in left_seq or 'N' * min_gap_length in right_seq:
                    self.logger.warning(
                        f"Gap {sid}.{gap_idx} 的flanking序列包含其他gap，"
                        f"表示两个gap距离太近|Flanking sequence of gap {sid}.{gap_idx} "
                        f"contains another gap. This indicates two gaps are too close."
                    )
                else:
                    flanking_dict[f'{sid}.{gap_idx}.L'] = left_seq
                    flanking_dict[f'{sid}.{gap_idx}.R'] = right_seq

                gap_dict[f'{sid}.{gap_idx}'] = seq[gap_site[0]:gap_site[1]]
                gap_idx += 1

        return flanking_dict, gap_dict

    def run_minimap(self, query_fa: str, ref_fa: str, output_prefix: str,
                    minimap_option: str = '-x asm5', threads: int = 1,
                    overwrite: bool = False) -> str:
        """运行minimap2比对|Run minimap2 alignment"""
        paf_file = f'{output_prefix}.paf'

        # 检查是否已存在
        if not overwrite and os.path.exists(paf_file) and os.path.getsize(paf_file) > 0:
            self.logger.info(f"使用现有PAF文件|Using existing PAF file: {paf_file}")
            return paf_file

        self.logger.info(f"运行minimap2...|Running minimap2...")
        cmd = f'minimap2 {minimap_option} -t {threads} {ref_fa} {query_fa} > {paf_file}'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"minimap2运行失败|minimap2 execution failed: {result.stderr}")

        return paf_file

    def parse_paf(self, paf_file: str, min_align_len: int = 1000,
                  min_identity: float = 0.4) -> List[Dict]:
        """解析PAF文件并过滤|Parse PAF file and filter alignments"""
        alignments = []

        with open(paf_file, 'r') as f:
            for line in f:
                fields = line.split()
                if len(fields) < 11:
                    continue

                aln = {
                    'qryid': fields[0],
                    'qrylen': int(fields[1]),
                    'qrystart': int(fields[2]) + 1,
                    'qryend': int(fields[3]),
                    'strand': fields[4],
                    'refid': fields[5],
                    'reflen': int(fields[6]),
                    'refstart': int(fields[7]) + 1,
                    'refend': int(fields[8]),
                    'match': int(fields[9]),
                    'alignlen': int(fields[10]),
                }

                aln['identity'] = aln['match'] / aln['alignlen']

                # 过滤条件
                if aln['alignlen'] >= min_align_len and aln['identity'] >= min_identity:
                    # 解析gap ID
                    parts = aln['qryid'].split('.')
                    aln['gapid'] = '.'.join(parts[:-1])
                    aln['LR'] = parts[-1]  # L或R
                    alignments.append(aln)

        return alignments

    def parse_paf_debug(self, paf_file: str, min_align_len: int = 1000,
                       min_identity: float = 0.4) -> Tuple[List[Dict], Dict]:
        """解析PAF文件并过滤（带调试信息）|Parse PAF file and filter (with debug info)"""
        alignments = []
        debug_stats = {
            'total_lines': 0,
            'filtered_by_align_len': 0,
            'filtered_by_identity': 0,
            'passed': 0
        }

        with open(paf_file, 'r') as f:
            for line in f:
                debug_stats['total_lines'] += 1
                fields = line.split()
                if len(fields) < 11:
                    continue

                aln = {
                    'qryid': fields[0],
                    'qrylen': int(fields[1]),
                    'qrystart': int(fields[2]) + 1,
                    'qryend': int(fields[3]),
                    'strand': fields[4],
                    'refid': fields[5],
                    'reflen': int(fields[6]),
                    'refstart': int(fields[7]) + 1,
                    'refend': int(fields[8]),
                    'match': int(fields[9]),
                    'alignlen': int(fields[10]),
                }

                aln['identity'] = aln['match'] / aln['alignlen']

                # 统计过滤原因
                if aln['alignlen'] < min_align_len:
                    debug_stats['filtered_by_align_len'] += 1
                    continue
                if aln['identity'] < min_identity:
                    debug_stats['filtered_by_identity'] += 1
                    continue

                debug_stats['passed'] += 1
                # 解析gap ID
                parts = aln['qryid'].split('.')
                aln['gapid'] = '.'.join(parts[:-1])
                aln['LR'] = parts[-1]  # L或R
                alignments.append(aln)

        return alignments, debug_stats

    def fill_gaps(self, draft_genome: str, unitig_file: str,
                  output_prefix: str, flanking_len: int = 5000,
                  min_align_len: int = 1000, min_identity: float = 0.4,
                  max_filling_len: int = 1000000, threads: int = 1,
                  minimap_option: str = '-x asm5', overwrite: bool = False) -> str:
        """
        执行gap填充|Perform gap filling

        返回|Returns:
            填充后的FASTA文件路径|Path to filled FASTA file
        """
        self.logger.info("=" * 80)
        self.logger.info("开始quarTeT Gap填充|Starting quarTeT Gap Filling")
        self.logger.info("=" * 80)

        # 创建临时目录
        tmp_dir = f'{output_prefix}.tmp'
        os.makedirs(tmp_dir, exist_ok=True)

        try:
            # 步骤1: 提取flanking序列
            self.logger.info(f"步骤1: 提取gap flanking序列|Step 1: Extracting gap flanking sequences")
            flanking_dict, gap_dict = self.extract_flanking_seqs(draft_genome, flanking_len)

            if not flanking_dict:
                self.logger.warning("未找到有效的gap|No valid gaps found")
                return draft_genome

            self.logger.info(f"  找到|Found {len(gap_dict)} 个gap")

            # 保存flanking序列
            flanking_fa = f'{tmp_dir}/flanking.fasta'
            with open(flanking_fa, 'w') as f:
                for sid, seq in flanking_dict.items():
                    f.write(f'>{sid}\n{seq}\n')

            # 步骤2: 比对unitig到flanking序列|Step 2: Aligning unitigs to flanking sequences
            self.logger.info(f"步骤2: 比对unitig到flanking序列|Step 2: Aligning unitigs to flanking sequences")

            # 如果unitig文件是GFA格式，需要转换为FASTA
            if unitig_file.endswith('.gfa'):
                unitig_fa = f'{tmp_dir}/unitig.fasta'
                self.logger.info(f"  转换GFA到FASTA|Converting GFA to FASTA: {unitig_fa}")
                with open(unitig_fa, 'w') as out_f, open(unitig_file) as in_f:
                    for line in in_f:
                        if line.startswith('S\t'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                out_f.write(f'>{parts[1]}\n{parts[2]}\n')
                unitig_for_align = unitig_fa
            else:
                unitig_for_align = unitig_file

            paf_file = self.run_minimap(unitig_for_align, flanking_fa,
                                       f'{tmp_dir}/align',
                                       minimap_option, threads, overwrite)

            # 步骤3: 解析比对结果|Step 3: Parsing alignments
            self.logger.info(f"步骤3: 解析比对结果|Step 3: Parsing alignments")
            alignments, debug_stats = self.parse_paf_debug(paf_file, min_align_len, min_identity)
            self.logger.info(f"  总比对数|Total alignments: {debug_stats['total_lines']}")
            self.logger.info(f"  因长度过滤|Filtered by length: {debug_stats['filtered_by_align_len']}")
            self.logger.info(f"  因同一性过滤|Filtered by identity: {debug_stats['filtered_by_identity']}")
            self.logger.info(f"  有效比对数量|Valid alignments: {debug_stats['passed']}")

            if not alignments:
                self.logger.warning("没有有效的比对，无法填充|No valid alignments, cannot fill gaps")
                return draft_genome

            # 步骤4: 为每个gap寻找最佳填充序列|Step 4: Finding best fill sequence for each gap
            self.logger.info(f"步骤4: 为每个gap寻找最佳填充序列|Step 4: Finding best fill sequence for each gap")

            # 读取unitig序列（自动检测格式）|Read unitig sequences (auto-detect format)
            if unitig_file.endswith('.gfa'):
                unitig_dict = self.read_gfa(unitig_file)
                self.logger.info(f"  检测到GFA格式|Detected GFA format")
            else:
                unitig_dict = self.read_fasta(unitig_file)
                self.logger.info(f"  检测到FASTA格式|Detected FASTA format")

            gapcloser_dict = {}  # {gapid: fill_info}

            for gapid in gap_dict:
                left_anchors = [aln for aln in alignments if aln['gapid'] == gapid and aln['LR'] == 'L']
                right_anchors = [aln for aln in alignments if aln['gapid'] == gapid and aln['LR'] == 'R']

                self.logger.info(f"  Gap {gapid}: {len(left_anchors)} left anchors, {len(right_anchors)} right anchors")

                if not left_anchors or not right_anchors:
                    self.logger.info(f"    跳过|Skipped: 缺少一侧的比对|missing alignments on one side")
                    continue

                # 寻找能同时比对到两侧的unitig
                valid_pairs = 0
                for l_aln in left_anchors:
                    for r_aln in right_anchors:
                        if l_aln['refid'] != r_aln['refid'] or l_aln['strand'] != r_aln['strand']:
                            continue

                        valid_pairs += 1

                        # 如果左侧比对在右侧之前
                        if l_aln['refend'] < r_aln['refstart']:
                            score = (l_aln['identity'] + r_aln['identity']) / 2

                            # 检查是否已有更好的填充
                            if gapid not in gapcloser_dict or score > gapcloser_dict[gapid]['score']:
                                fill_start = l_aln['refend'] + l_aln['qrylen'] - l_aln['qryend'] + 1
                                fill_end = r_aln['refstart'] - r_aln['qrystart']
                                fill_seq = unitig_dict[l_aln['refid']][fill_start-1:fill_end]

                                if len(fill_seq) > max_filling_len or not fill_seq:
                                    continue

                                gapcloser_dict[gapid] = {
                                    'seqid': l_aln['refid'],
                                    'range': f'{fill_start}-{fill_end}',
                                    'seq': fill_seq if l_aln['strand'] == '+' else self.reverse_complement(fill_seq),
                                    'strand': l_aln['strand'],
                                    'score': score
                                }

                self.logger.info(f"    找到|Found {valid_pairs} 个有效的左右比对对|valid left-right pairs")

            self.logger.info(f"  成功填充|Successfully filled: {len(gapcloser_dict)}/{len(gap_dict)} gaps")

            if not gapcloser_dict:
                self.logger.warning("没有gap被填充|No gaps were filled")
                return draft_genome

            # 步骤5: 生成填充后的基因组
            self.logger.info(f"步骤5: 生成填充后的基因组|Step 5: Generating filled genome")
            output_fa = f'{output_prefix}.filled.fasta'
            draft_dict = self.read_fasta(draft_genome)

            with open(output_fa, 'w') as w:
                for sid, seq in draft_dict.items():
                    if sid not in ['.'.join(gapid.split('.')[:-1]) for gapid in gapcloser_dict]:
                        # 没有被填充的gap，直接输出原序列
                        w.write(f'>{sid}\n{seq}\n')
                    else:
                        # 有被填充的gap，需要替换
                        matches = list(re.finditer(r'N{100,}', seq))
                        seq_list = []
                        start = 0

                        # 添加序列和gap的交替列表
                        for match in matches:
                            end = match.start()
                            seq_list.append((seq[start:end], 'sequence', sid, start + 1, end, '+'))
                            start = match.end()

                        seq_list.append((seq[start:], 'sequence', sid, start + 1, len(seq), '+'))

                        # 插入填充序列
                        filled_seq = ''
                        for i in range(len(matches)):
                            gapid = f'{sid}.{i+1}'
                            if gapid in gapcloser_dict:
                                # 填充gap
                                filled_seq += seq_list[2*i][0] + gapcloser_dict[gapid]['seq']
                            else:
                                # 保留原gap
                                filled_seq += seq_list[2*i][0] + gap_dict[gapid]

                        # 添加最后一段序列
                        filled_seq += seq_list[-1][0]

                        w.write(f'>{sid}\n{filled_seq}\n')

            self.logger.info(f"输出文件|Output file: {output_fa}")

            # 统计
            final_gap_count, _ = self.count_gaps(output_fa)
            original_gap_count = len(gap_dict)
            filled_count = original_gap_count - final_gap_count

            self.logger.info("=" * 80)
            self.logger.info("填充统计|Filling Statistics:")
            self.logger.info(f"  原始gap数量|Original gaps: {original_gap_count}")
            self.logger.info(f"  填充gap数量|Filled gaps: {filled_count}")
            self.logger.info(f"  剩余gap数量|Remaining gaps: {final_gap_count}")
            self.logger.info(f"  填充率|Fill rate: {filled_count/original_gap_count*100:.1f}%")
            self.logger.info("=" * 80)

            return output_fa

        finally:
            # 清理临时文件
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
