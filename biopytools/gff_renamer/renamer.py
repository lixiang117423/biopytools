"""
GFF文件重命名核心模块|GFF File Renamer Core Module
"""

import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
from multiprocessing import Pool, cpu_count
from .utils import parse_chromosome, generate_gene_id, load_chr_mapping, build_conda_command, CommandRunner


def _process_line_chunk(args):
    """
    处理行块的辅助函数（用于多进程）|Helper function to process line chunks (for multiprocessing)

    Args:
        args: (lines_chunk, id_mapping, start_line_index, gene_name_suffixes) 元组|tuple

    Returns:
        list: 处理后的行列表|List of processed lines
    """
    lines_chunk, id_mapping, start_line_index, gene_name_suffixes = args
    result = []

    for i, line in enumerate(lines_chunk):
        # 保留注释和空行|Keep comments and empty lines
        if line.startswith('#') or not line:
            result.append(line)
            continue

        fields = line.split('\t')
        if len(fields) < 9:
            result.append(line)
            continue

        # 计算原始文件的line_index|Calculate original file line index
        line_index = start_line_index + i

        # 获取特征类型|Get feature type
        feature_type = fields[2]

        # 更新第9列的属性|Update attributes in column 9
        attributes = fields[8]
        new_attributes = _update_and_clean_attributes(attributes, id_mapping, line_index, feature_type, gene_name_suffixes)
        fields[8] = new_attributes

        result.append('\t'.join(fields))

    return result


def _update_and_clean_attributes(attributes: str, id_mapping: Dict[str, str], line_index: int = None,
                               original_feature_type: str = None, gene_name_suffixes: Dict[str, str] = None) -> str:
    """
    更新ID并清理多余属性|Update IDs and clean redundant attributes

    规则|Rules:
        gene: ID=xxx;Name=xxx (如果是pseudogene则Name=xxx.pseudogene，如果是lncRNA基因则Name=xxx.lncRNA)
        转录本(mRNA/lncRNA等): ID=xxx;Name=xxx;Parent=xxx
        子特征(exon/CDS等): ID=xxx;Parent=xxx

    Args:
        attributes: 属性字符串|Attributes string
        id_mapping: ID映射|ID mapping
        line_index: 行索引（用于CDS映射）|Line index (for CDS mapping)
        original_feature_type: 原始特征类型|Original feature type
        gene_name_suffixes: 基因Name后缀|Gene Name suffixes

    Returns:
        str: 清理后的属性字符串|Cleaned attributes string
    """
    if not attributes:
        return attributes

    # 提取ID|Extract ID
    id_match = re.search(r'\bID=([^;]+)', attributes)
    if not id_match:
        return attributes

    old_id = id_match.group(1)

    # 尝试查找映射|Try to find mapping
    new_id = None

    # 首先尝试直接映射|First try direct mapping
    if old_id in id_mapping:
        new_id = id_mapping[old_id]
    # 如果是CDS，尝试用(line_index, old_id)查找|If CDS, try (line_index, old_id)
    elif line_index is not None and (line_index, old_id) in id_mapping:
        new_id = id_mapping[(line_index, old_id)]

    # 如果找不到映射，仍需检查Parent是否需要替换
    # AGAT清洗可能新增UTR等feature，其ID不在mapping中但Parent需要替换
    # If mapping not found, still check if Parent needs to be replaced
    # AGAT cleaning may add new UTR features whose IDs are not in mapping but Parents need replacement
    if new_id is None:
        parent_match = re.search(r'\bParent=([^;]+)', attributes)
        if parent_match:
            old_parent = parent_match.group(1)
            if old_parent in id_mapping:
                return attributes.replace(f'Parent={old_parent}', f'Parent={id_mapping[old_parent]}')
        return attributes

    # 提取Parent|Extract Parent
    parent_match = re.search(r'\bParent=([^;]+)', attributes)

    # 根据ID格式判断特征类型|Determine feature type based on ID format
    if '.' not in new_id:
        # gene: 没有后缀，是基因|No suffix, it's a gene
        name_suffix = ""

        # 如果是pseudogene，在Name中添加.pseudogene后缀
        # If pseudogene, add .pseudogene suffix to Name
        if original_feature_type == 'pseudogene':
            name_suffix = ".pseudogene"
        # 检查是否有lncRNA后缀|Check if has lncRNA suffix
        elif gene_name_suffixes and new_id in gene_name_suffixes:
            name_suffix = gene_name_suffixes[new_id]

        return f"ID={new_id};Name={new_id}{name_suffix}"
    else:
        # 有后缀，判断是转录本还是子特征|Has suffix, determine if transcript or child feature
        id_parts = new_id.split('.')
        if len(id_parts) >= 2:
            suffix = id_parts[1]
            # 检查是否是转录本类型|Check if it's a transcript type
            transcript_types = ['mRNA', 'lncRNA', 'miRNA', 'ncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA', 'transcript']
            utr_types = ['five_prime_UTR', 'three_prime_UTR']
            is_transcript = any(suffix.startswith(t) for t in transcript_types)
            is_utr = any(suffix.startswith(t) for t in utr_types)

            if is_transcript:
                # 转录本: ID;Name;Parent|Transcript: ID;Name;Parent
                new_parent = parent_match.group(1) if parent_match else ''
                if new_parent and new_parent in id_mapping:
                    new_parent = id_mapping[new_parent]
                return f"ID={new_id};Name={new_id};Parent={new_parent}"
            else:
                # 子特征(exon/CDS/UTR等): ID;Parent|Child feature: ID;Parent
                new_parent = parent_match.group(1) if parent_match else ''
                if new_parent and new_parent in id_mapping:
                    new_parent = id_mapping[new_parent]
                # 如果Parent不存在，也只返回ID（某些边缘情况）|If Parent doesn't exist, return only ID
                if new_parent:
                    return f"ID={new_id};Parent={new_parent}"
                else:
                    return f"ID={new_id}"

    return attributes


class GFFRenamer:
    """GFF文件重命名器|GFF File Renamer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def rename_gff(self) -> bool:
        """
        执行GFF文件重命名|Execute GFF file renaming

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("开始GFF文件重命名|Starting GFF file renaming")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出文件|Output file: {self.config.output_file}")
        self.logger.info(f"前缀|Prefix: {self.config.prefix}")
        self.logger.info(f"物种缩写|Species abbreviation: {self.config.species}")
        self.logger.info(f"命名格式|Naming format: {self.config.naming_format}")
        self.logger.info(f"包含UTR特征|Include UTR features: {self.config.include_utr}")

        # 步骤0：AGAT清洗GFF文件（默认启用）|Step 0: AGAT clean GFF file (enabled by default)
        actual_input = self.config.input_file
        if not self.config.skip_gff_clean:
            self.logger.info("步骤0: 使用AGAT清洗GFF文件|Step 0: Cleaning GFF file with AGAT")
            cmd_runner = CommandRunner(self.logger)

            fd, tmp_clean = tempfile.mkstemp(suffix='.gff', dir=str(Path(self.config.output_file).parent))
            os.close(fd)
            os.unlink(tmp_clean)  # 删除存根，让AGAT创建新文件|Remove stub so AGAT creates the file

            cmd_list = build_conda_command(
                self.config.agat_path,
                ['-g', self.config.input_file, '-o', tmp_clean]
            )
            cmd_str = ' '.join(cmd_list)

            if not cmd_runner.run(cmd_str, "AGAT清洗GFF文件|AGAT clean GFF file"):
                self.logger.error("AGAT清洗失败|AGAT cleaning failed")
                return False

            actual_input = tmp_clean
            self.logger.info(f"AGAT清洗完成，临时文件|AGAT cleaning completed, temp file: {tmp_clean}")
        else:
            self.logger.info("跳过GFF清洗（用户指定）|Skipping GFF cleaning (user specified)")

        # 加载染色体映射文件|Load chromosome mapping file
        chr_mapping = None
        if self.config.chr_mapping_file:
            self.logger.info(f"加载染色体映射文件|Loading chromosome mapping file: {self.config.chr_mapping_file}")
            chr_mapping = load_chr_mapping(self.config.chr_mapping_file)
            self.logger.info(f"加载了|Loaded {len(chr_mapping)} 个映射条目|mapping entries")

        # 步骤1：读取并解析GFF文件|Step 1: Read and parse GFF file
        self.logger.info("步骤1: 读取GFF文件|Step 1: Reading GFF file")
        lines = self._read_gff_file(actual_input)
        if not lines:
            self.logger.error("无法读取GFF文件或文件为空|Cannot read GFF file or file is empty")
            return False

        # 步骤1.5：prefer_mrna 去冗余(默认开启)|Step 1.5: prefer_mrna dedup (on by default)
        # 对含mRNA的基因,丢弃其冗余transcript(misc_RNA)变体及这些转录本的子特征
        # Drop redundant transcript (misc_RNA) variants and their children from mRNA genes
        if self.config.prefer_mrna:
            self.logger.info("步骤1.5: prefer_mrna去冗余|Step 1.5: prefer_mrna dedup")
            lines = self._filter_redundant_transcripts(lines)

        # 步骤2：解析所有特征并建立映射关系|Step 2: Parse all features and build mapping
        self.logger.info("步骤2: 解析特征并建立ID映射|Step 2: Parsing features and building ID mapping")
        id_mapping, gene_name_suffixes = self._build_id_mapping(lines, chr_mapping)
        self.logger.info(f"共映射|Total mapped: {len(id_mapping)} 个特征|features")

        # 步骤3：应用ID映射|Step 3: Apply ID mapping
        self.logger.info("步骤3: 应用ID映射到所有特征|Step 3: Applying ID mapping to all features")
        renamed_lines = self._apply_id_mapping(lines, id_mapping, gene_name_suffixes)

        # 步骤4：写入输出文件|Step 4: Write output file
        self.logger.info("步骤4: 写入输出文件|Step 4: Writing output file")
        self._write_output_file(renamed_lines)

        # 步骤5：生成mRNA映射文件（如果需要）|Step 5: Generate mRNA mapping file (if needed)
        if self.config.output_mrna_mapping:
            self.logger.info("步骤5: 生成mRNA映射文件|Step 5: Generating mRNA mapping file")
            self._write_mrna_mapping_file(lines, id_mapping)

        self.logger.info("GFF文件重命名完成|GFF file renaming completed")
        return True

    def _read_gff_file(self, input_file: str = None) -> List[str]:
        """读取GFF文件|Read GFF file"""
        if input_file is None:
            input_file = self.config.input_file
        lines = []
        with open(input_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                lines.append(line)
        return lines

    def _filter_redundant_transcripts(self, lines: List[str]) -> List[str]:
        """
        丢弃冗余转录本及其子特征|Drop redundant transcripts and their children

        对含有 mRNA 的基因,移除其 transcript(典型为 NCBI/Gnomon 的 misc_RNA 变体)特征,
        以及这些被移除转录本的子特征(exon/CDS/UTR/intron/codon)。
        仅含 transcript、不含 mRNA 的基因保留(可能是真非编码基因)。
        For genes that have mRNA, remove their `transcript` features (typically
        NCBI/Gnomon misc_RNA variants) and the children of those removed transcripts.
        Genes that have only `transcript` (no mRNA) are kept.
        """
        # 1) 含有 mRNA 的基因 parent 集合 | set of gene parents that have mRNA
        mrna_parents = set()
        for line in lines:
            fields = line.split('\t')
            if len(fields) >= 9 and fields[2] == 'mRNA':
                m = re.search(r'Parent=([^;]+)', fields[8])
                if m:
                    mrna_parents.add(m.group(1))
        if not mrna_parents:
            self.logger.info("prefer_mrna: 未发现mRNA特征,跳过去冗余|no mRNA features, skip dedup")
            return lines

        # 2) 冗余 transcript 的 old_id(parent基因含mRNA) | redundant transcript IDs
        redundant = set()
        for line in lines:
            fields = line.split('\t')
            if len(fields) >= 9 and fields[2] == 'transcript':
                pm = re.search(r'Parent=([^;]+)', fields[8])
                im = re.search(r'ID=([^;]+)', fields[8])
                if pm and im and pm.group(1) in mrna_parents:
                    redundant.add(im.group(1))
        if not redundant:
            self.logger.info("prefer_mrna: 含mRNA的基因无冗余transcript变体|no redundant transcript variants under mRNA genes")
            return lines

        # 3) 过滤:丢弃冗余transcript本身 + 其子特征(Parent∈redundant) | filter
        kept = []
        dropped_transcript = 0
        dropped_child = 0
        for line in lines:
            fields = line.split('\t')
            drop = False
            if len(fields) >= 9:
                im = re.search(r'ID=([^;]+)', fields[8])
                pm = re.search(r'Parent=([^;]+)', fields[8])
                old_id = im.group(1) if im else None
                parent = pm.group(1) if pm else None
                # 冗余transcript本身 | the redundant transcript row
                if fields[2] == 'transcript' and old_id in redundant:
                    drop = True
                    dropped_transcript += 1
                # 被丢弃transcript的子特征 | children of a dropped transcript
                elif parent in redundant:
                    drop = True
                    dropped_child += 1
            if not drop:
                kept.append(line)

        self.logger.info(
            f"prefer_mrna去冗余完成|prefer_mrna dedup done: "
            f"丢弃 {dropped_transcript} 个冗余transcript + {dropped_child} 个子特征 "
            f"(涉及含mRNA基因 {len(mrna_parents)} 个)"
        )
        return kept

    def _build_id_mapping(self, lines: List[str], chr_mapping: Dict[str, str] = None) -> Dict[str, str]:
        """
        构建完整的ID映射|Build complete ID mapping

        优化版本：单次遍历收集所有信息，然后在内存中处理
        Optimized version: Single pass to collect all info, then process in memory

        Args:
            lines: GFF文件行|GFF file lines
            chr_mapping: 染色体映射字典|Chromosome mapping dictionary (optional)

        Returns:
            dict: {旧ID: 新ID}|Dictionary of {old_ID: new_ID}
        """
        # 步骤1：单次遍历收集所有特征信息|Step 1: Single pass to collect all feature information
        features = {
            'genes': [],        # [(line_index, start, chr_num, old_gene_id), ...]
            'transcripts': [],  # [(line_index, old_id, parent_id, transcript_num, feature_type), ...]
            'exons': [],        # [(line_index, old_id, parent_id, start, end), ...]  # 添加位置信息
            'introns': [],      # [(line_index, old_id, parent_id, start, end), ...]  # intron特征
            'cds': [],          # [(line_index, old_id, parent_id, start, end), ...]  # 添加位置信息
            'utrs': [],         # [(line_index, old_id, parent_id, start, end, utr_type), ...]  # UTR特征
            'codons': []        # [(line_index, old_id, parent_id, start, end, codon_type), ...]  # start_codon/stop_codon
        }

        self.logger.info("正在收集特征信息|Collecting feature information...")

        for i, line in enumerate(lines):
            # 跳过注释和空行|Skip comments and empty lines
            if line.startswith('#') or not line:
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]

            # 收集gene特征（包括pseudogene假基因）|Collect gene features (including pseudogene)
            if feature_type in ['gene', 'pseudogene']:
                start = int(fields[3])
                seq_id = fields[0]
                std_chr, chr_num = parse_chromosome(seq_id, chr_mapping)

                # 提取gene ID，支持多种格式|Extract gene ID, support multiple formats
                # EGAPx格式: ID=gene-XXX|EGAPx format: ID=gene-XXX
                # EviAnn格式: ID=LOC_00000001 或其他任意ID|EviAnn format: ID=LOC_00000001 or any other ID
                match = re.search(r'ID=gene-([^;]+)', fields[8])
                if not match:
                    # 如果不是EGAPx格式，尝试提取任意ID|If not EGAPx format, try to extract any ID
                    match = re.search(r'ID=([^;]+)', fields[8])

                if match:
                    old_gene_id = match.group(1)
                    features['genes'].append((i, start, chr_num, old_gene_id))

            # 收集转录本特征(mRNA, lncRNA, miRNA, transcript等)|Collect transcript features (mRNA, lncRNA, miRNA, transcript, etc.)
            elif feature_type in ['mRNA', 'lncRNA', 'miRNA', 'ncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA', 'transcript']:
                id_match = re.search(r'ID=([^;]+)', fields[8])
                parent_match = re.search(r'Parent=([^;]+)', fields[8])

                if id_match and parent_match:
                    old_id = id_match.group(1)
                    old_parent = parent_match.group(1)

                    # 提取转录本编号，支持多种格式|Extract transcript number, support multiple formats
                    # BRAKER格式: gXXX.t1, gXXX.t2|BRAKER format: gXXX.t1, gXXX.t2
                    # EGAPx格式: mRNA-XXX-R1|EGAPx format: mRNA-XXX-R1
                    # EviAnn格式: XXX-mRNA-1 或 XXX-transcript-1|EviAnn format: XXX-mRNA-1 or XXX-transcript-1
                    mrna_match = re.search(r'\.t(\d+)$', old_id)
                    if not mrna_match:
                        mrna_match = re.search(r'-R(\d+)$', old_id)
                    if not mrna_match:
                        # 尝试EviAnn格式: -mRNA-N 或 -transcript-N|Try EviAnn format: -mRNA-N or -transcript-N
                        mrna_match = re.search(r'-(?:mRNA|transcript)-(\d+)$', old_id)

                    if mrna_match:
                        transcript_num = int(mrna_match.group(1))
                    else:
                        # 没有编号后缀的RNA，默认为1|RNA without number suffix, default to 1
                        transcript_num = 1

                    # 保持原始feature_type，transcript也当作独立类型处理
                    # Keep original feature_type, treat 'transcript' as independent type
                    features['transcripts'].append((i, old_id, old_parent, transcript_num, feature_type))

            # 收集start_codon和stop_codon特征|Collect start_codon and stop_codon features
            elif feature_type in ['start_codon', 'stop_codon']:
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)
                parent_match = re.search(r'Parent=([^;]+)', fields[8])
                if parent_match:
                    old_parent = parent_match.group(1)
                    features['codons'].append((i, old_id, old_parent, int(fields[3]), int(fields[4]), feature_type))

            # 收集exon特征|Collect exon features
            elif feature_type == 'exon':
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)

                # exon编号由位置顺序决定，不需要从旧ID中提取编号
                # Exon numbering determined by position order, no need to extract from old ID

                parent_match = re.search(r'Parent=([^;]+)', fields[8])
                if parent_match:
                    old_parent = parent_match.group(1)
                    start = int(fields[3])  # 起始位置|Start position
                    end = int(fields[4])    # 结束位置|End position
                    features['exons'].append((i, old_id, old_parent, start, end))

            # 收集intron特征|Collect intron features
            elif feature_type == 'intron':
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)

                parent_match = re.search(r'Parent=([^;]+)', fields[8])
                if parent_match:
                    old_parent = parent_match.group(1)
                    start = int(fields[3])
                    end = int(fields[4])
                    features['introns'].append((i, old_id, old_parent, start, end))

            # 收集CDS特征|Collect CDS features
            elif feature_type == 'CDS':
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)

                # CDS按出现顺序编号，不需要提取原始编号|CDS numbered by appearance, no need to extract original number
                parent_match = re.search(r'Parent=([^;]+)', fields[8])
                if parent_match:
                    old_parent = parent_match.group(1)
                    start = int(fields[3])  # 起始位置|Start position
                    end = int(fields[4])    # 结束位置|End position
                    features['cds'].append((i, old_id, old_parent, start, end))

            # 收集UTR特征(始终收集并重编号,保证输出ID唯一;
            # include_utr 现仅影响UTR是否作为子特征重排序/清理Name)
            # Always collect UTRs so their IDs are renumbered uniquely; include_utr now
            # only controls cosmetic child-reordering / Name stripping
            elif feature_type in ['five_prime_UTR', 'three_prime_UTR']:
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)

                parent_match = re.search(r'Parent=([^;]+)', fields[8])
                if parent_match:
                    old_parent = parent_match.group(1)
                    start = int(fields[3])  # 起始位置|Start position
                    end = int(fields[4])    # 结束位置|End position
                    features['utrs'].append((i, old_id, old_parent, start, end, feature_type))

        self.logger.info(f"收集完成|Collection complete: "
                        f"{len(features['genes'])} genes, "
                        f"{len(features['transcripts'])} transcripts, "
                        f"{len(features['exons'])} exons, "
                        f"{len(features['introns'])} introns, "
                        f"{len(features['cds'])} CDS, "
                        f"{len(features['codons'])} codons, "
                        f"{len(features['utrs'])} UTRs")

        # 步骤2：在内存中处理，建立ID映射|Step 2: Process in memory, build ID mapping
        self.logger.info("正在生成ID映射|Generating ID mapping...")

        id_mapping = {}
        gene_name_suffixes = {}  # 存储每个基因的Name后缀|Store Name suffix for each gene

        # 2.1 对genes按染色体和起始位置排序|Sort genes by chromosome and start position
        genes_sorted = sorted(features['genes'], key=lambda x: (x[2], x[1]))

        # 2.2 统计每条染色体的基因数量，用于生成编号|Count genes per chromosome for numbering
        chr_gene_counts = {}
        for _, _, chr_num, _ in genes_sorted:
            chr_gene_counts[chr_num] = chr_gene_counts.get(chr_num, 0) + 1

        # 2.3 为每个gene分配新ID|Assign new IDs to each gene
        chr_gene_indices = {}  # 当前染色体的基因索引计数器|Current chromosome gene index counter

        for _, start, chr_num, old_gene_id in genes_sorted:
            # 更新基因索引|Update gene index
            chr_gene_indices[chr_num] = chr_gene_indices.get(chr_num, 0) + 1
            gene_index = chr_gene_indices[chr_num]

            # 生成新gene ID|Generate new gene ID
            new_gene_id = generate_gene_id(
                self.config.prefix,
                self.config.species,
                chr_num,
                gene_index,
                self.config.naming_format
            )

            # 映射gene ID（带前缀和不带前缀）|Map gene ID (with and without prefix)
            id_mapping[f"gene-{old_gene_id}"] = new_gene_id
            id_mapping[old_gene_id] = new_gene_id  # 用于Name等属性|For Name and other attributes

        # 2.4 处理转录本|Process transcripts
        # 按parent分组|Group by parent
        transcripts_by_parent = defaultdict(list)
        for _, old_id, old_parent, transcript_num, feature_type in features['transcripts']:
            transcripts_by_parent[old_parent].append((old_id, transcript_num, feature_type))

        # 收集每个基因的转录本类型|Collect transcript types for each gene
        gene_transcript_types = defaultdict(set)
        for old_parent, transcripts in transcripts_by_parent.items():
            for _, _, feature_type in transcripts:
                gene_transcript_types[old_parent].add(feature_type)

        # 为每个基因的转录本分配新ID|Assign new IDs to transcripts of each gene
        for old_parent, transcripts in transcripts_by_parent.items():
            # 查找基因的新ID|Find gene's new ID
            gene_key = f"gene-{old_parent}" if not old_parent.startswith('gene-') else old_parent
            if gene_key not in id_mapping:
                # 如果找不到带前缀的，尝试不带前缀|If not found with prefix, try without prefix
                if old_parent in id_mapping:
                    new_gene_id = id_mapping[old_parent]
                else:
                    continue
            else:
                new_gene_id = id_mapping[gene_key]

            # 检查该基因的转录本类型|Check transcript types for this gene
            transcript_types = gene_transcript_types[old_parent]

            # 如果只有lncRNA（没有mRNA），记录Name后缀
            # If only lncRNA (no mRNA), record Name suffix
            if 'lncRNA' in transcript_types and 'mRNA' not in transcript_types:
                gene_name_suffixes[new_gene_id] = '.lncRNA'

            # 按转录本编号排序|Sort by transcript number
            transcripts_sorted = sorted(transcripts, key=lambda x: x[1])

            for idx, (old_mrna_id, _, feature_type) in enumerate(transcripts_sorted, 1):
                # 使用组内枚举序号，避免AGAT清洗后ID格式变化导致编号丢失
                # Use enumerate index to avoid lost transcript numbers after AGAT cleaning
                new_mrna_id = f"{new_gene_id}.{feature_type}{idx}"
                id_mapping[old_mrna_id] = new_mrna_id

        # 2.5 处理exon、CDS和UTR|Process exons, CDS and UTR
        # 按parent分组|Group by parent
        child_features_by_parent = defaultdict(list)

        # 收集所有exon|Collect all exons
        for line_index, old_id, old_parent, start, end in features['exons']:
            child_features_by_parent[old_parent].append(('exon', line_index, old_id, start, end))

        # 收集所有intron|Collect all introns
        for line_index, old_id, old_parent, start, end in features['introns']:
            child_features_by_parent[old_parent].append(('intron', line_index, old_id, start, end))

        # 收集所有CDS|Collect all CDS
        for line_index, old_id, old_parent, start, end in features['cds']:
            child_features_by_parent[old_parent].append(('cds', line_index, old_id, start, end))

        # 收集所有UTR(始终处理,保证重编号)|Always process UTRs to guarantee renumbering
        for line_index, old_id, old_parent, start, end, utr_type in features['utrs']:
            # 转换UTR类型为简短形式|Convert UTR type to short form
            if utr_type == 'five_prime_UTR':
                short_utr_type = 'five_prime_UTR'
            elif utr_type == 'three_prime_UTR':
                short_utr_type = 'three_prime_UTR'
            else:
                short_utr_type = utr_type
            child_features_by_parent[old_parent].append(('utr', line_index, old_id, start, end, short_utr_type))

        # 收集所有start_codon和stop_codon|Collect all start_codon and stop_codon
        for line_index, old_id, old_parent, start, end, codon_type in features['codons']:
            child_features_by_parent[old_parent].append(('codon', line_index, old_id, start, end, codon_type))

        # 为每个parent的exon、CDS和UTR按基因组位置统一排序后编号
        # Number exons, CDS and UTR sorted by genomic position (mixed together)
        for old_parent, feature_list in child_features_by_parent.items():
            # 查找Parent的新ID|Find Parent's new ID
            if old_parent not in id_mapping:
                continue

            new_parent = id_mapping[old_parent]

            # 按第4列（start位置）升序排序|Sort by column 4 (start position) ascending
            feature_sorted = sorted(feature_list, key=lambda x: x[3])

            # 分别为exon、CDS和UTR分配编号|Assign numbers to exons, CDS and UTR separately
            exon_counter = 0
            intron_counter = 0
            cds_counter = 0
            utr5_counter = 0
            utr3_counter = 0
            start_codon_counter = 0
            stop_codon_counter = 0

            for item in feature_sorted:
                feature_type = item[0]
                line_index = item[1]
                old_id = item[2]

                if feature_type == 'exon':
                    exon_counter += 1
                    new_id = f"{new_parent}.exon{exon_counter}"
                    # 用(line_index, old_id)作key:AGAT/合并可能产生共享old_id的行,
                    # 用纯old_id会被后续覆盖导致输出重复ID
                    # Key by (line_index, old_id): AGAT/merge can emit rows sharing an
                    # old_id; keying by old_id alone gets overwritten -> duplicate output IDs
                    id_mapping[(line_index, old_id)] = new_id
                elif feature_type == 'intron':
                    intron_counter += 1
                    new_id = f"{new_parent}.intron{intron_counter}"
                    id_mapping[(line_index, old_id)] = new_id
                elif feature_type == 'cds':
                    cds_counter += 1
                    new_id = f"{new_parent}.cds{cds_counter}"
                    # CDS使用(line_index, old_id)作为key，因为CDS可能有相同ID
                    # Use (line_index, old_id) as key for CDS, as CDS may have duplicate IDs
                    id_mapping[(line_index, old_id)] = new_id
                elif feature_type == 'utr':
                    # UTR需要区分5'和3'UTR|Need to distinguish 5' and 3' UTR
                    utr_type = item[5]  # five_prime_UTR or three_prime_UTR
                    if utr_type == 'five_prime_UTR':
                        utr5_counter += 1
                        new_id = f"{new_parent}.five_prime_UTR{utr5_counter}"
                    else:  # three_prime_UTR
                        utr3_counter += 1
                        new_id = f"{new_parent}.three_prime_UTR{utr3_counter}"
                    # UTR同样用(line_index, old_id)作key,避免共享old_id覆盖
                    # UTR also keyed by (line_index, old_id) to avoid shared-old_id overwrite
                    id_mapping[(line_index, old_id)] = new_id
                elif feature_type == 'codon':
                    codon_type = item[5]  # start_codon or stop_codon
                    if codon_type == 'start_codon':
                        start_codon_counter += 1
                        new_id = f"{new_parent}.start{start_codon_counter}"
                    else:  # stop_codon
                        stop_codon_counter += 1
                        new_id = f"{new_parent}.stop{stop_codon_counter}"
                    id_mapping[(line_index, old_id)] = new_id

        return id_mapping, gene_name_suffixes

    def _apply_id_mapping(self, lines: List[str], id_mapping: Dict[str, str],
                          gene_name_suffixes: Dict[str, str] = None) -> List[str]:
        """
        应用ID映射到所有行（并行处理）|Apply ID mapping to all lines (parallel processing)

        Args:
            lines: GFF文件行|GFF file lines
            id_mapping: ID映射|ID mapping
            gene_name_suffixes: 基因Name后缀|Gene Name suffixes

        Returns:
            list: 重命名后的行列表|List of renamed lines
        """
        # 获取线程数|Get number of threads
        num_threads = self.config.threads

        # 如果只有1个线程或文件很小，使用单进程|Use single process if only 1 thread or small file
        if num_threads <= 1 or len(lines) < 1000:
            self.logger.info(f"使用单进程处理|Using single process processing")
            return self._apply_id_mapping_single(lines, id_mapping, gene_name_suffixes)

        self.logger.info(f"使用并行处理|Using parallel processing: {num_threads} 线程|threads")

        try:
            # 计算每个块的大小|Calculate chunk size
            total_lines = len(lines)
            chunk_size = (total_lines + num_threads - 1) // num_threads  # 向上取整|Round up

            # 分割行列表|Split lines into chunks
            chunks = []
            for i in range(0, total_lines, chunk_size):
                chunk = lines[i:i + chunk_size]
                chunks.append((chunk, id_mapping, i, gene_name_suffixes))  # 添加起始行索引和gene_name_suffixes|Add start line index and gene_name_suffixes

            self.logger.info(f"分割为|Split into {len(chunks)} 个块|chunks, 每块|each chunk ~{chunk_size} 行|lines")

            # 使用多进程池处理|Use multiprocessing pool
            with Pool(processes=num_threads) as pool:
                # 并行处理所有块|Process all chunks in parallel
                results = pool.map(_process_line_chunk, chunks)

            # 合并结果|Merge results
            renamed_lines = []
            for result in results:
                renamed_lines.extend(result)

            self.logger.info(f"并行处理完成|Parallel processing completed")
            return renamed_lines

        except Exception as e:
            self.logger.warning(f"并行处理失败，回退到单进程模式|Parallel processing failed, fallback to single process: {e}")
            return self._apply_id_mapping_single(lines, id_mapping, gene_name_suffixes)

    def _apply_id_mapping_single(self, lines: List[str], id_mapping: Dict[str, str],
                                  gene_name_suffixes: Dict[str, str] = None) -> List[str]:
        """
        应用ID映射到所有行（单进程模式）|Apply ID mapping to all lines (single process mode)

        Args:
            lines: GFF文件行|GFF file lines
            id_mapping: ID映射|ID mapping
            gene_name_suffixes: 基因Name后缀|Gene Name suffixes

        Returns:
            list: 重命名后的行列表|List of renamed lines
        """
        renamed_lines = []

        for line_index, line in enumerate(lines):
            # 保留注释和空行|Keep comments and empty lines
            if line.startswith('#') or not line:
                renamed_lines.append(line)
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                renamed_lines.append(line)
                continue

            # 获取特征类型|Get feature type
            feature_type = fields[2]

            # 更新第9列的属性|Update attributes in column 9
            attributes = fields[8]
            new_attributes = _update_and_clean_attributes(attributes, id_mapping, line_index, feature_type, gene_name_suffixes)
            fields[8] = new_attributes

            renamed_lines.append('\t'.join(fields))

        return renamed_lines

    def _write_output_file(self, lines: List[str]):
        """写入输出文件|Write output file"""
        self.logger.info("正在清理exon和CDS的Name属性并重新排序|Cleaning Name attributes and reordering exon/CDS...")

        with open(self.config.output_file, 'w') as f:
            i = 0
            while i < len(lines):
                line = lines[i]

                # 跳过注释和空行|Skip comments and empty lines
                if line.startswith('#') or not line:
                    f.write(line + '\n')
                    i += 1
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    f.write(line + '\n')
                    i += 1
                    continue

                feature_type = fields[2]

                # 如果是mRNA等转录本，收集后续的所有exon和CDS进行重新排序
                # If mRNA or other transcript, collect all exons and CDS for reordering
                if feature_type in ['mRNA', 'lncRNA', 'miRNA', 'ncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA', 'transcript']:
                    # 直接写入mRNA行（已经在_apply_id_mapping中处理过）|Write mRNA line directly (already processed)
                    f.write(line + '\n')
                    i += 1

                    # 收集该转录本的所有exon和CDS|Collect all exons and CDS for this transcript
                    child_features = []
                    while i < len(lines):
                        child_line = lines[i]

                        # 跳过注释和空行|Skip comments and empty lines
                        if child_line.startswith('#') or not child_line:
                            i += 1
                            continue

                        child_fields = child_line.split('\t')
                        if len(child_fields) < 9:
                            break

                        child_type = child_fields[2]

                        # 只收集exon、intron、CDS、UTR和codon|Only collect exon, intron, CDS, UTR and codon
                        if child_type in ['exon', 'intron', 'CDS', 'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR']:
                            start = int(child_fields[3])
                            child_features.append((start, i, child_line))
                            i += 1
                        else:
                            # 遇到其他特征类型，停止收集|Stop collecting when hitting other feature types
                            break

                    # 按start位置排序|Sort by start position
                    child_features.sort(key=lambda x: x[0])

                    # 交替写入exon和CDS|Write exons and CDS alternatively
                    for start, line_idx, child_line in child_features:
                        child_fields = child_line.split('\t')
                        if len(child_fields) >= 9:
                            # 删除Name属性|Remove Name attribute
                            attributes = child_fields[8]
                            attributes = re.sub(r'Name=[^;]+;?', '', attributes)
                            attributes = attributes.replace(';;', ';')
                            attributes = attributes.rstrip(';')
                            child_fields[8] = attributes
                            child_line = '\t'.join(child_fields)
                        f.write(child_line + '\n')
                else:
                    # 其他特征类型，删除Name属性（如果需要）并写入
                    # Other feature types, remove Name attribute if needed and write
                    if feature_type in ['exon', 'intron', 'CDS', 'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR']:
                        attributes = fields[8]
                        attributes = re.sub(r'Name=[^;]+;?', '', attributes)
                        attributes = attributes.replace(';;', ';')
                        attributes = attributes.rstrip(';')
                        fields[8] = attributes
                        line = '\t'.join(fields)

                    f.write(line + '\n')
                    i += 1

        self.logger.info(f"输出文件已保存|Output file saved: {self.config.output_file}")

    def _write_mrna_mapping_file(self, lines: List[str], id_mapping: Dict[str, str]):
        """
        生成mRNA映射文件|Generate mRNA mapping file

        文件格式|File format:
        - 第一列：重命名后的ID|Column 1: New ID
        - 第二列：重命名前的ID|Column 2: Old ID
        - 第三列开始：原始第9列的属性（按;分割）|Column 3+: Original attributes from column 9 (split by ;)

        Args:
            lines: GFF文件行|GFF file lines
            id_mapping: ID映射字典|ID mapping dictionary
        """
        mrna_count = 0

        with open(self.config.mrna_mapping_file, 'w') as f:
            for line in lines:
                # 跳过注释和空行|Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    continue

                feature_type = fields[2]

                # 只处理mRNA|Only process mRNA
                if feature_type != 'mRNA':
                    continue

                # 提取旧ID|Extract old ID
                id_match = re.search(r'ID=([^;]+)', fields[8])
                if not id_match:
                    continue

                old_id = id_match.group(1)

                # 查找新ID|Find new ID
                if old_id not in id_mapping:
                    continue

                new_id = id_mapping[old_id]

                # 解析第9列的属性|Parse attributes in column 9
                attributes = fields[8]
                attr_parts = attributes.split(';')

                # 构建输出行：新ID\t旧ID\t属性1\t属性2\t...|Build output line: new_id\told_id\tattr1\tattr2\t...
                output_line = [new_id, old_id] + attr_parts
                f.write('\t'.join(output_line) + '\n')

                mrna_count += 1

        self.logger.info(f"mRNA映射文件已保存|mRNA mapping file saved: {self.config.mrna_mapping_file}")
        self.logger.info(f"共|Total {mrna_count} 个mRNA记录|mRNA records")
