"""
ANNOVAR结果处理模块|ANNOVAR Results Processing Module
功能：处理ANNOVAR生成的注释结果文件，包括外显子注释和所有注释结果
Functions: Process ANNOVAR generated annotation result files, including exonic and all annotation results
作者|Author: Xiang LI
版本|Version: v10 - 集成结果处理版本|Integrated results processing version
日期|Date: 2025-08-26
"""

import os
import re
import subprocess
import sys
from typing import List, Dict, Optional, Union, Any
from pathlib import Path


class ProteinSeqModifier:
    """蛋白序列修改器|Protein Sequence Modifier"""

    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def __init__(self, logger, pep_file: str, cds_file: str = None):
        self.logger = logger
        self.pep_sequences = {}
        self.cds_sequences = {}
        self._load_pep_sequences(pep_file)
        if cds_file:
            self._load_cds_sequences(cds_file)

    def _load_fasta(self, fasta_file: str) -> dict:
        """加载fasta文件到字典|Load fasta file into dictionary"""
        sequences = {}
        if not os.path.exists(fasta_file):
            return sequences
        current_id = None
        seq_parts = []
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(seq_parts)
                    current_id = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
        if current_id:
            sequences[current_id] = ''.join(seq_parts)
        return sequences

    def _load_pep_sequences(self, pep_file: str):
        """使用seqkit fx2tab加载蛋白序列|Load protein sequences using seqkit fx2tab"""
        seqkit_path = '/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/seqkit'
        try:
            result = subprocess.run(
                [seqkit_path, 'fx2tab', '-i', pep_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        self.pep_sequences[parts[0]] = parts[1]
                self.logger.info(f"已加载|Loaded {len(self.pep_sequences)} 条蛋白序列 (via seqkit fx2tab)")
            else:
                self.logger.warning(f"seqkit fx2tab执行失败，回退到普通加载|seqkit fx2tab failed, falling back to normal load: {result.stderr[:200]}")
                self.pep_sequences = self._load_fasta(pep_file)
                if self.pep_sequences:
                    self.logger.info(f"已加载|Loaded {len(self.pep_sequences)} 条蛋白序列 (fallback)")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.warning(f"seqkit不可用，回退到普通加载|seqkit unavailable, falling back to normal load: {e}")
            self.pep_sequences = self._load_fasta(pep_file)
            if self.pep_sequences:
                self.logger.info(f"已加载|Loaded {len(self.pep_sequences)} 条蛋白序列 (fallback)")

    def _load_cds_sequences(self, cds_file: str):
        """使用seqkit fx2tab加载CDS序列，避免大文件读取截断问题|Load CDS sequences using seqkit fx2tab to avoid large file truncation"""
        seqkit_path = '/share/org/YZWL/yzwl_lixg/miniforge3/envs/BioinfTools/bin/seqkit'
        try:
            result = subprocess.run(
                [seqkit_path, 'fx2tab', '-i', cds_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq = parts[1].replace('-', '')  # 移除gap字符|Remove gap characters
                        self.cds_sequences[seq_id] = seq
                self.logger.info(f"已加载|Loaded {len(self.cds_sequences)} 条CDS序列 (via seqkit fx2tab)")
            else:
                self.logger.warning(f"seqkit fx2tab执行失败，回退到普通加载|seqkit fx2tab failed, falling back to normal load: {result.stderr[:200]}")
                self.cds_sequences = self._load_fasta(cds_file)
                if self.cds_sequences:
                    self.logger.info(f"已加载|Loaded {len(self.cds_sequences)} 条CDS序列 (fallback)")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.warning(f"seqkit不可用，回退到普通加载|seqkit unavailable, falling back to normal load: {e}")
            self.cds_sequences = self._load_fasta(cds_file)
            if self.cds_sequences:
                self.logger.info(f"已加载|Loaded {len(self.cds_sequences)} 条CDS序列 (fallback)")

    def _translate(self, cds_seq: str) -> str:
        """翻译CDS序列为蛋白序列|Translate CDS sequence to protein sequence"""
        protein = []
        for i in range(0, len(cds_seq) - 2, 3):
            codon = cds_seq[i:i+3].upper()
            aa = self.CODON_TABLE.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        return ''.join(protein)

    def _apply_cdna_to_cds(self, cds_seq: str, cdna_change: str) -> str:
        """在CDS上引入cDNA变异|Apply cDNA change to CDS sequence"""
        if 'delins' in cdna_change:
            m = re.match(r'(\d+)_(\d+)delins([A-Z]+)', cdna_change)
            if m:
                start = int(m.group(1)) - 1
                end = int(m.group(2))
                alt = m.group(3)
                if start < len(cds_seq):
                    return cds_seq[:start] + alt + cds_seq[end:]
            m = re.match(r'(\d+)delins([A-Z]+)', cdna_change)
            if m:
                pos = int(m.group(1)) - 1
                alt = m.group(2)
                if pos < len(cds_seq):
                    return cds_seq[:pos] + alt + cds_seq[pos + 1:]
        elif 'del' in cdna_change:
            m = re.match(r'(\d+)_(\d+)del([A-Z]*)', cdna_change)
            if m:
                start = int(m.group(1)) - 1
                end = int(m.group(2))
                return cds_seq[:start] + cds_seq[end:]
            else:
                m = re.match(r'(\d+)del([A-Z]+)', cdna_change)
                if m:
                    pos = int(m.group(1)) - 1
                    return cds_seq[:pos] + cds_seq[pos + 1:]
        elif 'ins' in cdna_change:
            m = re.match(r'(\d+)_(\d+)ins([A-Z]+)', cdna_change)
            if m:
                pos = int(m.group(1))
                alt = m.group(3)
                if pos <= len(cds_seq):
                    return cds_seq[:pos] + alt + cds_seq[pos:]
        else:
            m = re.match(r'([A-Z])(\d+)([A-Z])', cdna_change)
            if m:
                pos = int(m.group(2)) - 1
                if pos < len(cds_seq):
                    return cds_seq[:pos] + m.group(3) + cds_seq[pos + 1:]
        return cds_seq

    def get_ref_protein(self, transcript_id: str) -> str:
        return self.pep_sequences.get(transcript_id, 'NA')

    def get_mutant_protein(self, transcript_id: str, cdna_change: str, protein_change: str, ref_protein: str) -> str:
        """生成变异蛋白序列|Generate mutant protein sequence"""
        if ref_protein == 'NA' or not ref_protein:
            return 'NA'

        if not protein_change and not cdna_change:
            return ref_protein

        # 有p.XXX注释时优先走蛋白层面处理|Use protein-level approach when p.XXX annotation exists
        handled = False
        if protein_change:
            if 'delins' in protein_change:
                m = re.match(r'([A-Z]+)(\d+)delins([A-Z]+)', protein_change)
                if m:
                    pos = int(m.group(2))
                    ref_len = len(m.group(1))
                    alt = m.group(3)
                    if pos - 1 + ref_len <= len(ref_protein):
                        return ref_protein[:pos - 1] + alt + ref_protein[pos - 1 + ref_len:]
                    handled = True
            elif 'del' in protein_change:
                m = re.match(r'(\d+)_(\d+)del', protein_change)
                if m:
                    start = int(m.group(1))
                    end = int(m.group(2))
                    if end <= len(ref_protein):
                        return ref_protein[:start - 1] + ref_protein[end:]
                    handled = True
            elif 'X' == protein_change[-1] and re.match(r'[A-Z]\d+X', protein_change):
                m = re.match(r'([A-Z])(\d+)X', protein_change)
                if m:
                    pos = int(m.group(2))
                    if pos <= len(ref_protein):
                        return ref_protein[:pos - 1] + '*'
                    handled = True
            elif 'fs' not in protein_change:
                m = re.match(r'([A-Z])(\d+)([A-Z])', protein_change)
                if m:
                    pos = int(m.group(2))
                    alt = m.group(3)
                    ref_aa = m.group(1)
                    if pos <= len(ref_protein):
                        return ref_protein[:pos - 1] + alt + ref_protein[pos:]
                    elif ref_aa == 'X':
                        # 终止丢失(stoploss): pos超过蛋白长度，走CDS重翻译|Stoploss: pos exceeds protein length, fall through to CDS re-translation
                        pass
                    else:
                        handled = True

        # 蛋白层面无法处理时走CDS重新翻译|Fall back to CDS-level re-translation
        if not handled:
            cds_seq = self.cds_sequences.get(transcript_id, '')
            if cds_seq and cdna_change:
                try:
                    mutant_cds = self._apply_cdna_to_cds(cds_seq, cdna_change)
                    mutant_protein = self._translate(mutant_cds)
                    if mutant_protein:
                        return mutant_protein
                except Exception:
                    pass

        return ref_protein


class ExonicVariantProcessor:
    """外显子变异结果处理器|Exonic Variant Result Processor"""

    def __init__(self, logger):
        self.logger = logger
        self.protein_modifier = None

    def set_protein_modifier(self, modifier: ProteinSeqModifier):
        """设置蛋白序列修改器|Set protein sequence modifier"""
        self.protein_modifier = modifier

    @staticmethod
    def _extract_protein_change(gene_info: str, transcript_id: str) -> str:
        """从原始注释中提取指定转录本的蛋白变化|Extract protein change for a specific transcript from raw annotation"""
        for entry in gene_info.split(','):
            entry = entry.strip()
            if transcript_id in entry:
                m = re.match(r'[^:]+:[^:]+:[^:]+:c\.([^:]+)(?::p\.([^,]+))?', entry)
                if m:
                    return m.group(2) if m.group(2) else ''
        return ''

    @staticmethod
    def _extract_cdna_change(gene_info: str, transcript_id: str) -> str:
        """从原始注释中提取指定转录本的cDNA变化|Extract cDNA change for a specific transcript from raw annotation"""
        for entry in gene_info.split(','):
            entry = entry.strip()
            if transcript_id in entry:
                m = re.match(r'[^:]+:[^:]+:[^:]+:c\.([^:]+)', entry)
                if m:
                    return m.group(1)
        return ''

    def _parse_cdna_change(self, cdna_change: str) -> Dict[str, str]:
        """解析cDNA变化|Parse cDNA change"""
        dna_pos_from = dna_pos_to = dna_ref = dna_alt = 'NA'

        if 'delins' in cdna_change:
            match = re.match(r'(\d+)delins([A-Z]+)', cdna_change)
            if match:
                dna_pos_from = dna_pos_to = match.group(1)
                dna_ref = 'NA'
                dna_alt = match.group(2)
        elif 'del' in cdna_change:
            match = re.match(r'(\d+)_(\d+)del([A-Z]*)', cdna_change)
            if match:
                dna_pos_from = match.group(1)
                dna_pos_to = match.group(2)
                dna_ref = match.group(3) if match.group(3) else 'NA'
                dna_alt = '-'
            else:
                match = re.match(r'(\d+)del([A-Z]+)', cdna_change)
                if match:
                    dna_pos_from = dna_pos_to = match.group(1)
                    dna_ref = match.group(2)
                    dna_alt = '-'
        elif 'ins' in cdna_change:
            match = re.match(r'(\d+)_(\d+)ins([A-Z]+)', cdna_change)
            if match:
                dna_pos_from = match.group(1)
                dna_pos_to = match.group(2)
                dna_alt = match.group(3)
                dna_ref = '-'
        else:
            match = re.match(r'([A-Z])(\d+)([A-Z])', cdna_change)
            if match:
                dna_ref = match.group(1)
                dna_pos_from = dna_pos_to = match.group(2)
                dna_alt = match.group(3)

        return {
            'DNA位置起': dna_pos_from, 'DNA位置止': dna_pos_to,
            'DNA参考': dna_ref, 'DNA变异': dna_alt
        }

    def _parse_protein_change(self, protein_change: str) -> Dict[str, str]:
        """解析蛋白变化|Parse protein change"""
        protein_pos = protein_ref = protein_alt = 'NA'

        if 'fs' in protein_change:
            match = re.match(r'([A-Z])(\d+)fs', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = 'frameshift'
        elif 'delins' in protein_change:
            match = re.match(r'([A-Z]+)(\d+)delins([A-Z]+)', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = match.group(3)
        elif 'del' in protein_change:
            match = re.match(r'(\d+)_(\d+)del', protein_change)
            if match:
                protein_pos = match.group(1)
                protein_ref = 'NA'
                protein_alt = 'NA'
        else:
            match = re.match(r'([A-Z])(\d+)([A-Z])', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = match.group(3)

        return {
            '蛋白位置': protein_pos, '蛋白参考': protein_ref, '蛋白变异': protein_alt
        }

    def parse_annovar_line(self, line: str) -> List[Dict[str, str]]:
        """解析ANNOVAR外显子注释输出的每一行，每个转录本输出一行|Parse each line, output one row per transcript"""
        fields = line.strip().split('\t')

        if len(fields) < 9:
            return []

        line_id = fields[0]
        variant_type = fields[1]  # 变异类型描述|Variant type description
        gene_info = fields[2]     # 基因注释信息|Gene annotation info
        chrom = fields[3]         # 染色体|Chromosome
        start = fields[4]         # 起始位置|Start position
        end = fields[5]           # 终止位置|End position
        ref = fields[6]           # 参考碱基|Reference base
        alt = fields[7]           # 变异碱基|Alternative base

        # 判断突变类型|Determine mutation type (SNP vs INDEL)
        if ref == '-' or alt == '-':
            mutation_type = 'INDEL'
        elif len(ref) == 1 and len(alt) == 1:
            mutation_type = 'SNP'
        elif len(ref) != len(alt):
            mutation_type = 'INDEL'
        else:
            mutation_type = 'COMPLEX'

        # 确定变异结果类型|Determine variant effect type
        # 注意：nonframeshift包含frameshift子串，需先匹配|nonframeshift contains 'frameshift' substring
        if 'nonframeshift' in variant_type:
            effect = '非移码突变'
        elif 'frameshift' in variant_type:
            effect = '移码突变'
        elif 'nonsynonymous' in variant_type:
            effect = '错义突变'
        elif 'synonymous' in variant_type:
            effect = '同义突变'
        elif 'stopgain' in variant_type:
            effect = '无义突变'
        elif 'stoploss' in variant_type:
            effect = '终止丢失'
        else:
            effect = variant_type

        # 解析基因信息，每个转录本生成一行|Parse gene info, one row per transcript
        rows = []
        for entry in gene_info.split(','):
            entry = entry.strip()
            if not entry:
                continue
            m = re.match(r'([^:]+):([^:]+):([^:]+):c\.([^:]+)(?::p\.([^,]+))?', entry)
            if not m:
                continue
            gene = m.group(1)
            transcript = m.group(2)
            cdna_change = m.group(4)
            protein_change = m.group(5) if m.group(5) else ''

            dna = self._parse_cdna_change(cdna_change)
            prot = self._parse_protein_change(protein_change)

            rows.append({
                'Line_ID': line_id,
                '染色体': chrom,
                '变异起始': start,
                '变异终止': end,
                '突变类型': mutation_type,
                '基因': gene,
                '转录本': transcript,
                '变异结果': effect,
                '原始注释': gene_info,
                **dna,
                **prot
            })

        return rows

    def process_exonic_variants(self, input_file: str, output_file: Optional[str] = None) -> List[Dict[str, str]]:
        """处理外显子变异注释结果|Process exonic variant annotation results"""
        self.logger.info(f" 处理外显子变异文件|Processing exonic variant file: {input_file}")

        results = []

        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('line'):  # 跳过表头|Skip header
                    rows = self.parse_annovar_line(line)
                    for row in rows:
                        if self.protein_modifier:
                            transcript = row['转录本']
                            cdna_change_raw = self._extract_cdna_change(row['原始注释'], transcript)
                            protein_change_raw = self._extract_protein_change(row['原始注释'], transcript)
                            ref_pep = self.protein_modifier.get_ref_protein(transcript)
                            mut_pep = self.protein_modifier.get_mutant_protein(
                                transcript, cdna_change_raw, protein_change_raw, ref_pep
                            )
                            row['参考蛋白序列'] = ref_pep
                            row['变异蛋白序列'] = mut_pep
                        else:
                            row['参考蛋白序列'] = 'NA'
                            row['变异蛋白序列'] = 'NA'
                        results.append(row)

        # 输出结果|Output results
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as out:
                self._write_exonic_results(results, out)
            self.logger.info(f" 外显子变异结果已保存|Exonic variant results saved: {output_file}")
        else:
            self._write_exonic_results(results, sys.stdout)

        self.logger.info(f" 共处理|Total processed: {len(results)} 个外显子变异")
        return results

    def _write_exonic_results(self, results: List[Dict[str, str]], output):
        """写入外显子变异结果|Write exonic variant results"""
        # 写入表头|Write headers
        headers = ['Line_ID', '染色体', '变异起始', '变异终止', '突变类型', '转录本',
                   '基因', '变异结果', '原始注释', '参考蛋白序列', '变异蛋白序列',
                   'DNA位置起', 'DNA位置止', 'DNA参考', 'DNA变异', '蛋白位置', '蛋白参考', '蛋白变异']
        output.write('\t'.join(headers) + '\n')

        # 写入数据|Write data
        for result in results:
            output.write('\t'.join(str(result[h]) for h in headers) + '\n')


class AllVariantProcessor:
    """所有变异结果处理器|All Variant Result Processor"""

    def __init__(self, logger):
        self.logger = logger

    def parse_variant_function_line(self, line: str) -> Optional[Dict[str, Union[str, int]]]:
        """解析ANNOVAR variant_function文件的每一行|Parse each line of ANNOVAR variant_function file"""
        fields = line.strip().split('\t')

        if len(fields) < 7:
            return None

        # 基本字段|Basic fields
        region_type = fields[0]      # 区域类型|Region type (intergenic, exonic, intronic等)
        gene_info = fields[1]        # 基因信息|Gene info
        chrom = fields[2]            # 染色体|Chromosome
        start = fields[3]            # 起始位置|Start position
        end = fields[4]              # 终止位置|End position
        ref = fields[5]              # 参考碱基|Reference base
        alt = fields[6]              # 变异碱基|Alternative base

        # 额外信息|Additional info
        if len(fields) >= 10:
            freq = fields[7]         # 频率|Frequency
            qual_score = fields[8]   # 质量分数|Quality score
            depth = fields[9]        # 测序深度|Sequencing depth
        else:
            freq = qual_score = depth = 'NA'

        # 判断突变类型|Determine mutation type
        if ref == '-':
            mutation_type = 'INSERTION'
            var_length = len(alt)
        elif alt == '-':
            mutation_type = 'DELETION'
            var_length = len(ref)
        elif len(ref) == 1 and len(alt) == 1:
            mutation_type = 'SNP'
            var_length = 1
        elif len(ref) == len(alt):
            mutation_type = 'MNP'  # 多核苷酸多态性|Polymorphism
            var_length = len(ref)
        else:
            mutation_type = 'COMPLEX'
            var_length = max(len(ref), len(alt))

        # 解析基因信息|Parse gene information
        # 格式|Format: NONE(dist=NONE),NONE(dist=NONE) 或 GeneA(dist=100),GeneB(dist=200)
        genes = []
        distances = []

        if gene_info != 'NONE(dist=NONE),NONE(dist=NONE)':
            # 提取基因名和距离|Extract gene names and distances
            gene_matches = re.findall(r'([^(,]+)\(dist=([^)]+)\)', gene_info)
            for gene_match in gene_matches:
                gene_name = gene_match[0]
                distance = gene_match[1]
                if gene_name != 'NONE':
                    genes.append(gene_name)
                    distances.append(distance)

        gene_str = ','.join(genes) if genes else 'NONE'
        distance_str = ','.join(distances) if distances else 'NONE'

        # 计算变异大小（对于INDEL）| Calculate variant size (for INDELs)
        if mutation_type in ['INSERTION', 'DELETION']:
            indel_size = abs(int(end) - int(start) + 1)
            if mutation_type == 'INSERTION':
                indel_size = len(alt)
        else:
            indel_size = 0

        return {
            '染色体': chrom,
            '起始位置': start,
            '终止位置': end,
            '区域类型': region_type,
            '基因': gene_str,
            '距离': distance_str,
            '突变类型': mutation_type,
            '参考序列': ref,
            '变异序列': alt,
            '变异长度': var_length,
            'INDEL大小': indel_size if mutation_type in ['INSERTION', 'DELETION'] else 'NA',
            '频率': freq,
            '质量分数': qual_score,
            '测序深度': depth
        }

    def filter_variants(self, results: List[Dict], filters: Optional[Dict[str, Any]] = None) -> List[Dict]:
        """
        过滤变异|Filter variants
        filters: 字典，包含过滤条件|Dictionary containing filter criteria
        例如|e.g: {'region_type': ['exonic', 'splicing'], 'mutation_type': ['SNP']}
        """
        if not filters:
            return results

        filtered = []
        for result in results:
            keep = True

            # 按区域类型过滤|Filter by region type
            if 'region_type' in filters:
                if result['区域类型'] not in filters['region_type']:
                    keep = False

            # 按突变类型过滤|Filter by mutation type
            if 'mutation_type' in filters:
                if result['突变类型'] not in filters['mutation_type']:
                    keep = False

            # 按基因过滤（排除intergenic）| Filter by genes (exclude intergenic)
            if filters.get('exclude_intergenic', False):
                if result['区域类型'] == 'intergenic' and result['基因'] == 'NONE':
                    keep = False

            # 按频率过滤|Filter by frequency
            if 'min_freq' in filters:
                try:
                    if float(result['频率']) < filters['min_freq']:
                        keep = False
                except (ValueError, TypeError):
                    pass

            if keep:
                filtered.append(result)

        return filtered

    def process_all_variants(self, input_file: str, output_file: Optional[str] = None,
                           apply_filters: bool = False, filters: Optional[Dict] = None) -> List[Dict]:
        """处理所有变异注释结果|Process all variant annotation results"""
        self.logger.info(f" 处理所有变异文件|Processing all variants file: {input_file}")

        results = []

        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                parsed = self.parse_variant_function_line(line)
                if parsed:
                    results.append(parsed)

        # 应用过滤（可选）| Apply filters (optional)
        if apply_filters and filters:
            original_count = len(results)
            results = self.filter_variants(results, filters)
            self.logger.info(f" 过滤从 {original_count} 个变异保留到 {len(results)} 个|"
                           f"Filtered from {original_count} to {len(results)} variants")

        # 输出结果|Output results
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as out:
                self._write_all_results(results, out)
            self.logger.info(f" 所有变异结果已保存|All variant results saved: {output_file}")
        else:
            self._write_all_results(results, sys.stdout)

        self.logger.info(f" 共处理|Total processed: {len(results)} 个变异")

        # 输出统计信息|Output statistics
        self._print_statistics(results)

        return results

    def _write_all_results(self, results: List[Dict], output):
        """写入所有变异结果|Write all variant results"""
        # 写入表头|Write headers
        headers = ['染色体', '起始位置', '终止位置', '区域类型', '基因', '距离',
                   '突变类型', '参考序列', '变异序列', '变异长度', 'INDEL大小',
                   '频率', '质量分数', '测序深度']
        output.write('\t'.join(headers) + '\n')

        # 写入数据|Write data
        for result in results:
            output.write('\t'.join(str(result[h]) for h in headers) + '\n')

    def _print_statistics(self, results: List[Dict]):
        """打印统计信息|Print statistics"""
        self.logger.info(" 变异统计信息|Variant Statistics:")

        # 按区域类型统计|Statistics by region type
        region_counts = {}
        for r in results:
            region = r['区域类型']
            region_counts[region] = region_counts.get(region, 0) + 1

        self.logger.info(" 区域类型分布|Region Type Distribution:")
        for region, count in sorted(region_counts.items(), key=lambda x: x[1], reverse=True):
            self.logger.info(f"  {region}: {count}")

        # 按突变类型统计|Statistics by mutation type
        mut_counts = {}
        for r in results:
            mut_type = r['突变类型']
            mut_counts[mut_type] = mut_counts.get(mut_type, 0) + 1

        self.logger.info(" 突变类型分布|Mutation Type Distribution:")
        for mut_type, count in sorted(mut_counts.items(), key=lambda x: x[1], reverse=True):
            self.logger.info(f"  {mut_type}: {count}")

        # 基因相关变异统计|Gene-related variant statistics
        gene_variants = sum(1 for r in results if r['基因'] != 'NONE')
        self.logger.info(f" 与基因相关的变异|Gene-related variants: {gene_variants}")
        self.logger.info(f" 基因间区变异|Intergenic variants: {len(results) - gene_variants}")


class ANNOVARResultsProcessor:
    """ANNOVAR结果处理器主类|Main ANNOVAR Results Processor Class"""

    def __init__(self, logger, output_dir: str, pep_file: str = None, cds_file: str = None):
        self.logger = logger
        self.output_dir = output_dir
        self.pep_file = pep_file
        self.cds_file = cds_file
        self.exonic_processor = ExonicVariantProcessor(logger)
        self.all_processor = AllVariantProcessor(logger)

    def process_exonic_results(self, exonic_file: str, output_prefix: Optional[str] = None) -> Optional[str]:
        """处理外显子注释结果|Process exonic annotation results"""
        if self.pep_file and not self.exonic_processor.protein_modifier:
            if os.path.exists(self.pep_file):
                self.exonic_processor.set_protein_modifier(
                    ProteinSeqModifier(self.logger, self.pep_file, self.cds_file)
                )
            else:
                self.logger.warning(f" 蛋白序列文件不存在，跳过蛋白序列|Protein file does not exist, skipping protein sequences: {self.pep_file}")

        if not os.path.exists(exonic_file):
            self.logger.warning(f" 外显子注释文件不存在|Exonic annotation file does not exist: {exonic_file}")
            return None

        if output_prefix is None:
            output_prefix = os.path.splitext(os.path.basename(exonic_file))[0]

        output_file = os.path.join(self.output_dir, f"{output_prefix}_processed_exonic.tsv")

        try:
            results = self.exonic_processor.process_exonic_variants(exonic_file, output_file)
            self.logger.info(f" 外显子注释结果处理完成|Exonic annotation results processed: {len(results)} variants")
            return output_file
        except Exception as e:
            self.logger.error(f" 处理外显子注释结果失败|Failed to process exonic annotation results: {str(e)}")
            return None

    def process_all_results(self, variant_function_file: str, output_prefix: Optional[str] = None,
                          apply_filters: bool = False, filters: Optional[Dict] = None) -> Optional[str]:
        """处理所有注释结果|Process all annotation results"""
        if not os.path.exists(variant_function_file):
            self.logger.warning(f" 变异功能注释文件不存在|Variant function annotation file does not exist: {variant_function_file}")
            return None

        if output_prefix is None:
            output_prefix = os.path.splitext(os.path.basename(variant_function_file))[0]

        output_file = os.path.join(self.output_dir, f"{output_prefix}_processed_all.tsv")

        # 默认过滤器（如果启用）| Default filters (if enabled)
        if apply_filters and filters is None:
            filters = {'exclude_intergenic': True}

        try:
            results = self.all_processor.process_all_variants(
                variant_function_file, output_file, apply_filters, filters
            )
            self.logger.info(f" 所有注释结果处理完成|All annotation results processed: {len(results)} variants")
            return output_file
        except Exception as e:
            self.logger.error(f" 处理所有注释结果失败|Failed to process all annotation results: {str(e)}")
            return None

    def process_available_results(self, vcf_basename: str, apply_filters: bool = False) -> Dict[str, str]:
        """处理可用的注释结果文件|Process available annotation result files"""
        processed_files = {}

        # 检查外显子注释文件|Check exonic annotation file
        exonic_file = os.path.join(self.output_dir, f"{vcf_basename}.exonic_variant_function")
        if os.path.exists(exonic_file):
            processed_file = self.process_exonic_results(exonic_file, vcf_basename)
            if processed_file:
                processed_files['exonic'] = processed_file

        # 检查所有变异注释文件|Check all variant annotation file
        all_variants_file = os.path.join(self.output_dir, f"{vcf_basename}.variant_function")
        if os.path.exists(all_variants_file):
            processed_file = self.process_all_results(all_variants_file, vcf_basename, apply_filters)
            if processed_file:
                processed_files['all'] = processed_file

        self.logger.info(f" 结果处理完成|Results processing completed: {len(processed_files)} 个文件已处理")
        return processed_files