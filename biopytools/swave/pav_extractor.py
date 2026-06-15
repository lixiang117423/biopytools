"""
PAV矩阵提取模块|PAV Matrix Extraction Module
从swave converted VCF中提取Presence/Absence Variation矩阵
Extract Presence/Absence Variation matrix from swave converted VCF
"""

import os
import gzip
import logging
from typing import List, Dict, Tuple, Optional


class PAVExtractor:
    """PAV矩阵提取器|PAV Matrix Extractor"""

    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger

    @staticmethod
    def _infer_svtype(ref_len: int, alt_len: int) -> str:
        """根据REF/ALT长度推断SV类型|Infer SV type from REF/ALT lengths"""
        diff = alt_len - ref_len
        if diff > 0:
            return 'INS'
        elif diff < 0:
            return 'DEL'
        else:
            return 'DUP'

    @staticmethod
    def _parse_gt_alleles(gt: str) -> set:
        """解析GT字段，返回包含的等位基因编号集合|Parse GT field, return set of allele indices"""
        if gt in ('.', './.', '.|.'):
            return set()
        return set(gt.replace('|', '/').split('/'))

    def extract(
        self,
        vcf_file: str,
        output_file: str,
        min_ac: int = 1,
        strip_prefix: bool = True,
        svtype_only: Optional[List[str]] = None
    ) -> str:
        """
        从swave converted VCF提取PAV矩阵

        Args:
            vcf_file: VCF文件路径（支持.gz）|VCF file path (supports .gz)
            output_file: 输出TSV文件路径|Output TSV file path
            min_ac: 最小等位基因数，过滤AC < min_ac的位点|Minimum allele count filter
            strip_prefix: 是否去除CHROM中的样本前缀（如T2T#0#）|Strip sample prefix from CHROM
            svtype_only: 仅保留指定SV类型|Only keep specified SV types (e.g., ['DUP', 'INS', 'DEL'])

        Returns:
            输出文件路径|Output file path
        """
        opener = gzip.open if vcf_file.endswith('.gz') else open

        # 第一遍：解析header获取样本名|Pass 1: parse header to get sample names
        samples = []
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    samples = parts[9:]
                    break

        if not samples:
            raise ValueError("无法解析VCF header中的样本名|Failed to parse sample names from VCF header")

        if self.logger:
            self.logger.info(f"样本数|Samples: {len(samples)}")

        # 第二遍：提取PAV数据|Pass 2: extract PAV data
        records = []
        total_variants = 0
        filtered_by_ac = 0
        filtered_by_type = 0

        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue

                total_variants += 1
                chrom = parts[0]
                pos = parts[1]
                ref_seq = parts[3]
                alt_seqs = parts[4].split(',')
                info_str = parts[7]

                # 解析INFO字段|Parse INFO field
                info = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        info[item.split('=')[0]] = item.split('=', 1)[1]

                # AC过滤|AC filter
                ac = int(info.get('AC', '0'))
                if ac < min_ac:
                    filtered_by_ac += 1
                    continue

                # 去除CHROM前缀|Strip CHROM prefix
                if strip_prefix and '#' in chrom:
                    chrom = chrom.split('#')[-1]

                # 确定SV类型：优先从INFO取，否则从样本FORMAT的TYPE字段取
                # Determine SV type: prefer INFO field, then sample FORMAT TYPE field
                svtype = info.get('SVTYPE', '.')
                if svtype == '.':
                    for i in range(9, len(parts)):
                        fmt_fields = parts[i].split(':')
                        if len(fmt_fields) >= 2 and fmt_fields[1] != '.':
                            svtype = fmt_fields[1]
                            break

                # SV类型过滤（整体位点级别）|SV type filter (site-level)
                if svtype_only and svtype != '.':
                    type_match = False
                    for st in svtype_only:
                        if st in svtype:
                            type_match = True
                            break
                    if not type_match:
                        filtered_by_type += 1
                        continue

                # 对每个ALT等位基因分别生成一条PAV记录|Generate one PAV record per ALT allele
                for alt_idx, alt_seq in enumerate(alt_seqs):
                    allele_idx = alt_idx + 1  # VCF中ALT allele编号从1开始

                    # 该等位基因的SVTYPE：如果父记录有则复用，否则按REF/ALT长度推断
                    allele_svtype = svtype
                    if allele_svtype == '.':
                        allele_svtype = self._infer_svtype(len(ref_seq), len(alt_seq))

                    # 该等位基因的SVLEN：按REF/ALT长度差计算
                    allele_svlen = abs(len(alt_seq) - len(ref_seq))

                    # 构建PAV记录：GT中包含该等位基因编号则为1
                    # GT containing this allele index = 1 (present)
                    pav = {}
                    for i, sample in enumerate(samples):
                        col_idx = i + 9
                        if col_idx < len(parts):
                            gt = parts[col_idx].split(':')[0]
                            alleles = self._parse_gt_alleles(gt)
                            pav[sample] = 1 if str(allele_idx) in alleles else 0
                        else:
                            pav[sample] = 0

                    records.append((
                        chrom, pos, ref_seq, alt_seq,
                        len(ref_seq), len(alt_seq),
                        allele_svtype, allele_svlen, ac, pav
                    ))

        # 写入TSV|Write TSV
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        with open(output_file, 'w') as f:
            # Header
            header = ['CHROM', 'POS', 'REF', 'ALT', 'REF_LEN', 'ALT_LEN', 'SVTYPE', 'SVLEN', 'AC'] + samples
            f.write('\t'.join(header) + '\n')

            # Data rows
            for chrom, pos, ref_seq, alt_seq, ref_len, alt_len, svtype, svlen, ac, pav in records:
                row = [chrom, pos, ref_seq, alt_seq, ref_len, alt_len, svtype, svlen, ac]
                row.extend([pav[s] for s in samples])
                f.write('\t'.join(str(v) for v in row) + '\n')

        if self.logger:
            n_present = sum(1 for r in records if any(v == 1 for v in r[9].values()))
            self.logger.info(f"总变异数|Total variants: {total_variants}")
            self.logger.info(f"AC<{min_ac}过滤|Filtered by AC<{min_ac}: {filtered_by_ac}")
            if svtype_only:
                self.logger.info(f"SVTYPE过滤|Filtered by SVTYPE: {filtered_by_type}")
            self.logger.info(f"最终PAV位点数|Final PAV sites: {len(records)}")
            self.logger.info(f"有变异的位点|Sites with variation: {n_present}")
            self.logger.info(f"输出文件|Output: {output_file}")

        return output_file
