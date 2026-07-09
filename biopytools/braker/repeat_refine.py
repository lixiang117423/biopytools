"""
BRAKER 重复序列库过滤 + 证据驱动还原|Repeat library filtering + evidence-based rescue

背景|Background:
    RepeatModeler(de novo 建库)会把多拷贝基因家族(如效应子)误判为重复家族,
    RepeatMasker 据此将其 soft-mask,导致 BRAKER 注释丢失这些基因。
    本模块提供两层修复:
      1. filter: 从 repeat 库剔除"实为基因家族"的 consensus(TE domain 排除 + RxLR motif)
      2. rescue: 用蛋白/转录组证据把被误 mask 的真基因区还原(unmask)

设计文档|Design doc:
    docs/superpowers/specs/2026-07-09-braker-repeat-family-rescue-design.md
"""

import os
import re
from typing import List, Tuple, Dict, Optional

from Bio.Seq import Seq
from Bio import SeqIO


# ===== 常量|Constants =====

# TE 特征 Pfam 家族关键词(用于判定真 TE,命中即保留)
# 宁窄勿宽:漏判 TE 会保守保留(无害),误判基因家族为 TE 会导致该基因家族仍被 mask
# TE-specific Pfam family keywords (match => keep as real TE)
TE_DOMAIN_KEYWORDS = (
    'transposase', 'integrase', 'reverse transcriptase', 'retrotransposon',
    'retropepsin', 'gag', 'rnase h', 'rnaseh', 'rve', 'rvt',
    'pif', 'harbinger', 'mutator', 'dde', 'piggybac', 'maverick',
    'helicase-rt', 'line-1', 'copia', 'gypsy', 'bel-pao', 'penelope',
    'group-specific', ' capsid',
)

# RxLR/EER 效应子 motif 正则(疫霉效应子特征)|RxLR/EER effector motifs
# RxLR: R-x-L-R (含 QxLR/GxLR 变体); EER: E-E-R (常紧跟 RxLR)
RXLR_MOTIF_RE = re.compile(r'[RGQ].LR', re.IGNORECASE)
EER_MOTIF_RE = re.compile(r'E[EDQ]R', re.IGNORECASE)

# 六框翻译最小 ORF 长度(氨基酸)|Min ORF length (aa) for six-frame translation
DEFAULT_MIN_ORF_LEN = 30
# 标准遗传密码|Standard genetic code
STANDARD_TABLE = 1


# ===== 翻译与 ORF 提取|Translation & ORF extraction =====

def _extract_orfs_from_peptide(pep: str, min_len: int) -> List[str]:
    """
    从单框翻译肽段(含 *)提取所有 ≥ min_len 的连续无 * 段|Extract ORFs (stop-free segments)

    对 repeat consensus,ORF 不强制从 M 开始(TE/基因蛋白可能有变异),
    取所有连续无终止子的长段作为 Pfam 扫描输入。
    """
    orfs = []
    for segment in pep.split('*'):
        seg = segment.replace('X', '').strip()
        if len(seg) >= min_len:
            orfs.append(seg)
    return orfs


def six_frame_orfs(seq_str: str, min_orf_len: int = DEFAULT_MIN_ORF_LEN) -> List[str]:
    """
    六框翻译 DNA 序列,返回所有 ≥ min_orf_len 的 ORF 蛋白|Six-frame translate, return ORFs

    Args:
        seq_str: DNA 序列(可能含小写/N)|DNA sequence (may contain lowercase/N)
        min_orf_len: 最小 ORF 长度(氨基酸)|Min ORF length (aa)

    Returns:
        ORF 蛋白序列列表|List of ORF protein sequences
    """
    seq = seq_str.upper()
    # N 转为 A 以避免翻译异常(Biopython 遇 N 可能报错)|N -> A to avoid translate errors
    seq = seq.replace('N', 'A')
    rc = str(Seq(seq).reverse_complement())
    all_orfs = []
    for frame in range(3):
        pep_fwd = str(Seq(seq[frame:]).translate(table=STANDARD_TABLE))
        pep_rev = str(Seq(rc[frame:]).translate(table=STANDARD_TABLE))
        all_orfs.extend(_extract_orfs_from_peptide(pep_fwd, min_orf_len))
        all_orfs.extend(_extract_orfs_from_peptide(pep_rev, min_orf_len))
    return all_orfs


# ===== Pfam / TE domain 判定|Pfam TE-domain detection =====

def _is_te_domain(family_name: str) -> bool:
    """判断 Pfam 家族名是否为 TE 特征域|Check if Pfam family name is a TE domain"""
    name_lower = family_name.lower()
    return any(kw in name_lower for kw in TE_DOMAIN_KEYWORDS)


def parse_hmmscan_domtblout(domtblout: str, evalue_cutoff: float) -> Dict[str, List[Tuple[str, str, float]]]:
    """
    解析 hmmscan --domtblout 输出|Parse hmmscan domtblout

    Args:
        domtblout: domtblout 文件路径|Path to domtblout file
        evalue_cutoff: domain 独立 E-value 阈值|Domain i-Evalue cutoff

    Returns:
        {orf_id: [(pfam_acc, pfam_name, ievalue)]}
    """
    hits: Dict[str, List[Tuple[str, str, float]]] = {}
    if not os.path.exists(domtblout):
        return hits
    with open(domtblout) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            # domtblout 列: target_name(0) acc(1) tlen(2) query_name(3) query_acc(4)
            # qlen(5) E-value(6) ... '#'(9) of(10) c-Evalue(11) i-Evalue(12) ...
            if len(cols) < 13:
                continue
            orf_id = cols[0]
            pfam_name = cols[3]
            pfam_acc = cols[4]
            try:
                ievalue = float(cols[12])
            except (ValueError, IndexError):
                continue
            if ievalue <= evalue_cutoff:
                hits.setdefault(orf_id, []).append((pfam_acc, pfam_name, ievalue))
    return hits


# ===== RxLR/EER motif 扫描|RxLR/EER motif scan =====

def scan_effector_motif(proteins: List[str]) -> bool:
    """
    扫描蛋白列表是否含 RxLR 或 EER 效应子 motif|Scan for RxLR/EER effector motifs

    RxLR/EER 是疫霉效应子特征 motif,命中强烈提示"效应子基因家族"。
    """
    for prot in proteins:
        if RXLR_MOTIF_RE.search(prot) or EER_MOTIF_RE.search(prot):
            return True
    return False


# ===== masked 区提取与还原|Masked-region extraction & unmask =====

def extract_masked_regions(masked_fa: str, min_len: int) -> List[Tuple[str, int, int]]:
    """
    从 masked genome 提取所有连续小写(soft-mask)段|Extract continuous lowercase regions

    Returns:
        [(chrom, start_1based, end_1based)] 坐标为 1-based inclusive
        长度 < min_len 的短段忽略(避免噪声)
    """
    regions: List[Tuple[str, int, int]] = []
    for rec in SeqIO.parse(masked_fa, 'fasta'):
        seq = str(rec.seq)
        chrom = rec.id
        n = len(seq)
        i = 0
        while i < n:
            if seq[i].islower():
                j = i
                while j < n and seq[j].islower():
                    j += 1
                # i..j-1 为 0-based 小写段 → 1-based [i+1, j]
                if j - i >= min_len:
                    regions.append((chrom, i + 1, j))
                i = j
            else:
                i += 1
    return regions


def write_refined_masked(masked_fa: str, rescue_regions: List[Tuple[str, int, int]],
                         refined_fa: str) -> int:
    """
    把 rescue 区域在 masked genome 中转大写(unmask),写出 refined genome|Unmask rescue regions

    Args:
        masked_fa: 输入 masked genome|Input masked genome
        rescue_regions: [(chrom, start_1based, end_1based)]|Rescue regions (1-based inclusive)
        refined_fa: 输出 refined genome 路径|Output refined genome path

    Returns:
        unmask 的碱基数|Number of unmasked bases
    """
    # 按 chrom 组织 rescue 区,1-based → 0-based 区间集合
    rescue_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, s, e in rescue_regions:
        rescue_by_chrom.setdefault(chrom, []).append((s - 1, e))  # 0-based [s-1, e)

    total_unmasked = 0
    tmp_out = refined_fa + '.tmp'
    with open(tmp_out, 'w') as out:
        for rec in SeqIO.parse(masked_fa, 'fasta'):
            # bytearray 高效原地修改|Efficient in-place edit via bytearray
            ba = bytearray(str(rec.seq).encode('ascii'))
            for start, end in rescue_by_chrom.get(rec.id, []):
                start = max(0, start)
                end = min(len(ba), end)
                if start < end:
                    seg = ba[start:end]
                    ba[start:end] = bytes(seg).upper()
                    total_unmasked += (end - start)
            seq_str = ba.decode('ascii')
            out.write(f">{rec.id}\n")
            for k in range(0, len(seq_str), 60):
                out.write(seq_str[k:k + 60] + '\n')
    os.replace(tmp_out, refined_fa)
    return total_unmasked


# ===== miniprot PAF 解析|miniprot PAF parsing =====

def parse_miniprot_paf(paf_file: str) -> List[Tuple[str, int, int, float]]:
    """
    解析 miniprot 蛋白→基因组比对 PAF|Parse miniprot PAF

    Returns:
        [(chrom, target_start_1based, target_end_1based, identity_percent)]
        target 坐标转 1-based inclusive
    """
    hits: List[Tuple[str, int, int, float]] = []
    if not os.path.exists(paf_file):
        return hits
    with open(paf_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                continue
            chrom = cols[5]
            tstart = int(cols[7])  # 0-based half-open
            tend = int(cols[8])
            identity = 0.0
            # miniprot PAF tag id:f: 为 identity 百分比|identity as percentage
            for tag in cols[12:]:
                if tag.startswith('id:f:'):
                    try:
                        identity = float(tag[5:])  # "id:f:" 长度 5|len("id:f:")==5
                    except ValueError:
                        pass
                    break
            hits.append((chrom, tstart + 1, tend, identity))  # 1-based inclusive
    return hits


# ===== RNA-seq 覆盖度计算|RNA-seq depth =====

def compute_region_mean_depth(bam_files: List[str], regions: List[Tuple[str, int, int]],
                              regions_bed: str, samtools_bin: str, cmd_runner,
                              logger) -> Dict[int, float]:
    """
    用 samtools depth 计算每 region 的平均 RNA-seq 覆盖度(多 bam 取最大)|Mean depth per region

    Args:
        bam_files: BAM 路径列表|List of BAM paths
        regions: [(chrom, start_1based, end_1based)]|Regions (1-based inclusive)
        regions_bed: 临时 bed 文件路径|Temp bed path
        samtools_bin: samtools 路径|samtools path
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        {region_index: mean_depth}
    """
    if not bam_files or not regions:
        return {}

    # 写 BED(0-based half-open)|Write BED (0-based half-open)
    with open(regions_bed, 'w') as f:
        for chrom, s, e in regions:
            f.write(f"{chrom}\t{s - 1}\t{e}\n")

    depth_tsv = regions_bed + '.depth.tsv'
    # samtools depth -a 输出所有 position(含 0)|-a outputs all positions incl. 0
    bam_arg = ' '.join(bam_files)
    cmd = f"{samtools_bin} depth -a -b {regions_bed} {bam_arg} > {depth_tsv}"
    if not cmd_runner.run_command(cmd, "计算RNA-seq覆盖度|Compute RNA-seq depth"):
        logger.warning("RNA-seq覆盖度计算失败,跳过RNA-seq证据|RNA-seq depth failed, skip")
        return {}

    # 解析 per-position depth,按 chrom 组织排序|Parse per-position depth by chrom
    depth_by_chrom: Dict[str, List[Tuple[int, float]]] = {}
    with open(depth_tsv) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            chrom = cols[0]
            pos = int(cols[1])  # 1-based
            depths = [float(x) for x in cols[2:] if x]
            max_d = max(depths) if depths else 0.0
            depth_by_chrom.setdefault(chrom, []).append((pos, max_d))

    for chrom in depth_by_chrom:
        depth_by_chrom[chrom].sort()

    # 对每 region 区间求平均(双指针归并)|Mean over each region interval
    import bisect
    region_mean: Dict[int, float] = {}
    for idx, (chrom, s, e) in enumerate(regions):
        positions = depth_by_chrom.get(chrom, [])
        if not positions:
            region_mean[idx] = 0.0
            continue
        # 用 bisect 定位区间|Locate interval with bisect
        keys = [p[0] for p in positions]
        lo = bisect.bisect_left(keys, s)
        hi = bisect.bisect_right(keys, e)
        chunk = positions[lo:hi]
        if chunk:
            region_mean[idx] = sum(d for _, d in chunk) / (e - s + 1)
        else:
            region_mean[idx] = 0.0
    return region_mean


def _prepare_pfam_db(pfam_db: str, logger) -> Optional[str]:
    """
    确保 Pfam-A.hmm 就绪(必要时从 tar.gz 自动解压)|Ensure Pfam-A.hmm ready

    hmmscan 需要已 hmmpress 的 Pfam-A.hmm(+.h3i/.h3f/.h3m/.h3p 索引)。
    eggnog 打包的 pfam.tar.gz 内已含索引,解压即可用。
    """
    if os.path.exists(pfam_db):
        return pfam_db
    # 候选 tar.gz: 父目录的父目录下 pfam.tar.gz 等|Candidate tar.gz locations
    base_dir = os.path.dirname(os.path.dirname(pfam_db))
    tar_candidates = [
        os.path.join(base_dir, 'pfam.tar.gz'),
        pfam_db + '.tar.gz',
    ]
    tar_gz = next((t for t in tar_candidates if os.path.exists(t)), None)
    if not tar_gz:
        logger.error(f"Pfam 库不存在且找不到 tar.gz|Pfam DB missing, no tar.gz found: {pfam_db}")
        return None
    logger.info(f"首次运行,解压 Pfam 库|First run, extracting Pfam: {tar_gz}")
    import tarfile
    with tarfile.open(tar_gz, 'r:gz') as tar:
        tar.extractall(base_dir)
    if os.path.exists(pfam_db):
        logger.info(f"Pfam 解压完成|Pfam extracted: {pfam_db}")
        return pfam_db
    logger.error(f"解压后仍找不到 Pfam-A.hmm|Pfam-A.hmm still missing after extract: {pfam_db}")
    return None


# ===== 主函数:方案1 repeat 库过滤|Main: repeat library filter =====

def filter_repeat_library(consensi_fa: str, output_dir: str, config, cmd_runner,
                          logger) -> Optional[str]:
    """
    过滤 RepeatModeler 的 repeat 库,剔除被误判的基因家族|Filter repeat library

    判据|Criteria:
        1. 六框翻译 consensus → 蛋白 ORF
        2. hmmscan Pfam-A.hmm:
           - 命中 TE 特征域 → 真 TE,保留
           - 命中非 TE 蛋白域 → 基因家族,剔除
        3. RxLR/EER motif 命中 → 效应子家族,剔除
        4. 无 domain/motif → 保守保留(可能非编码 repeat)

    Args:
        consensi_fa: RepeatModeler 的 consensi.fa.classified|RepeatModeler consensi file
        output_dir: 输出目录(放中间与结果文件)|Output dir
        config: BrakerConfig(提供 hmmscan_bin/pfam_db/threads/te_domain_evalue)|Config
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        过滤后的 consensi 文件路径,失败返回 None|Filtered consensi path or None
    """
    logger.info("=" * 70)
    logger.info("repeat 库过滤(剔除误判的基因家族)|Repeat library filtering")
    logger.info("=" * 70)

    consensi_fa = os.path.abspath(consensi_fa)
    output_dir = os.path.abspath(output_dir)
    hmmscan_bin = config.hmmscan_bin
    pfam_db = config.pfam_db
    evalue_cutoff = config.te_domain_evalue
    min_orf_len = getattr(config, 'filter_min_orf_len', DEFAULT_MIN_ORF_LEN)

    # 确保 Pfam 库就绪(必要时从 tar.gz 自动解压)|Ensure Pfam DB ready
    pfam_db = _prepare_pfam_db(pfam_db, logger)
    if not pfam_db:
        logger.error("Pfam 库不可用,跳过 repeat 库过滤|Pfam unavailable, skip filtering")
        return None

    os.makedirs(output_dir, exist_ok=True)
    orfs_fa = os.path.join(output_dir, "consensi_orfs.fa")
    domtblout = os.path.join(output_dir, "consensi_orfs.domtblout")
    filtered_fa = os.path.join(output_dir, "filtered_consensi.fa.classified")
    report_tsv = os.path.join(output_dir, "filtered_families.tsv")

    # 1. 读 consensi,六框翻译,写 ORF FASTA|Read consensi, six-frame translate
    logger.info("六框翻译 consensus 提取 ORF|Six-frame translating consensi")
    orf_records: List[Tuple[str, str]] = []  # (orf_id, protein)
    family_to_orfs: Dict[str, List[str]] = {}
    total_families = 0
    for rec in SeqIO.parse(consensi_fa, 'fasta'):
        total_families += 1
        # header: >rnd-1_family-619#Unknown (...); family_id 取 '#' 前|family_id before '#'
        family_id = rec.id.split('#')[0]
        orfs = six_frame_orfs(str(rec.seq), min_orf_len)
        for k, prot in enumerate(orfs):
            orf_id = f"{family_id}__orf{k}"
            orf_records.append((orf_id, prot))
            family_to_orfs.setdefault(family_id, []).append(orf_id)

    with open(orfs_fa, 'w') as f:
        for orf_id, prot in orf_records:
            f.write(f">{orf_id}\n{prot}\n")
    logger.info(f"共 {total_families} 个 family,提取 {len(orf_records)} 条 ORF|"
                f"{total_families} families, {len(orf_records)} ORFs")

    if not orf_records:
        logger.warning("未提取到 ORF,跳过过滤(保留原库)|No ORF extracted, skip filtering")
        return None

    # 2. hmmscan 扫 Pfam|Scan Pfam with hmmscan
    logger.info(f"hmmscan 扫描 Pfam domain(线程 {config.threads})|hmmscan Pfam")
    cmd = f"{hmmscan_bin} --cpu {config.threads} --domtblout {domtblout} {pfam_db} {orfs_fa}"
    if not cmd_runner.run_command(cmd, "hmmscan 扫描 Pfam|hmmscan Pfam scan"):
        logger.error("hmmscan 失败,跳过 repeat 库过滤|hmmscan failed, skip filtering")
        return None

    # 3. 解析 domtblout,判定每 family|Parse domtblout, classify each family
    orf_hits = parse_hmmscan_domtblout(domtblout, evalue_cutoff)
    family_decision: Dict[str, str] = {}  # family_id -> 'keep_te'|'remove_gene'|'keep'
    family_evidence: Dict[str, str] = {}

    for family_id, orf_ids in family_to_orfs.items():
        te_hit = False
        gene_hit = False
        evidence_parts = []
        for orf_id in orf_ids:
            for pfam_acc, pfam_name, iev in orf_hits.get(orf_id, []):
                if _is_te_domain(pfam_name):
                    te_hit = True
                    evidence_parts.append(f"TE:{pfam_name}")
                else:
                    gene_hit = True
                    evidence_parts.append(f"gene:{pfam_name}")

        # RxLR/EER motif 检查|Check RxLR/EER motif on this family's ORFs
        family_proteins = [prot for oid, prot in orf_records if oid in orf_ids]
        has_effector_motif = scan_effector_motif(family_proteins)
        if has_effector_motif:
            evidence_parts.append("RxLR/EER")

        if te_hit:
            family_decision[family_id] = 'keep_te'
        elif gene_hit or has_effector_motif:
            family_decision[family_id] = 'remove_gene'
        else:
            family_decision[family_id] = 'keep'
        family_evidence[family_id] = ';'.join(sorted(set(evidence_parts)))

    n_remove = sum(1 for v in family_decision.values() if v == 'remove_gene')
    n_te = sum(1 for v in family_decision.values() if v == 'keep_te')
    logger.info(f"判定结果|Classification: 剔除(基因家族)={n_remove}, "
                f"保留(TE)={n_te}, 保留(其他)={total_families - n_remove - n_te}")

    # 4. 写过滤后的 consensi(剔除 remove_gene 的)|Write filtered consensi
    kept = 0
    with open(filtered_fa, 'w') as out:
        for rec in SeqIO.parse(consensi_fa, 'fasta'):
            family_id = rec.id.split('#')[0]
            if family_decision.get(family_id) == 'remove_gene':
                continue
            # 重写为标准 FASTA(去 description 的括号问题)|Rewrite clean FASTA
            out.write(f">{rec.id}\n{str(rec.seq)}\n")
            kept += 1

    # 5. 写剔除报告|Write removal report
    with open(report_tsv, 'w') as out:
        out.write("family_id\tdesignation\tevidence\n")
        for rec in SeqIO.parse(consensi_fa, 'fasta'):
            family_id = rec.id.split('#')[0]
            decision = family_decision.get(family_id, 'keep')
            label = {'keep_te': 'KEEP_TE', 'remove_gene': 'REMOVE',
                     'keep': 'KEEP'}.get(decision, 'KEEP')
            out.write(f"{family_id}\t{label}\t{family_evidence.get(family_id, '')}\n")

    logger.info(f"过滤完成|Filtering done: 保留 {kept}/{total_families} family → {filtered_fa}")
    logger.info(f"剔除报告|Removal report: {report_tsv}")
    return filtered_fa


# ===== 主函数:方案2 证据驱动还原|Main: evidence-based rescue =====

def rescue_masked_regions(masked_fa: str, prot_seq: Optional[str],
                          bam_files: List[str], output_dir: str, config,
                          cmd_runner, logger) -> Optional[str]:
    """
    用蛋白/转录组证据还原被误 mask 的真基因区|Rescue mis-masked gene regions by evidence

    证据源(自适应)|Evidence sources (adaptive):
        - 蛋白证据(prot_seq 存在时,必用): miniprot 蛋白→基因组,落在 masked 区且高 identity
        - 转录组证据(bam_files 非空时,可选): samtools depth masked 区覆盖度

    满足蛋白或转录组证据的 masked 段 → unmask(小写转大写)|Unmask regions with evidence support

    Args:
        masked_fa: masked genome|Masked genome
        prot_seq: 蛋白序列(可为 None)|Protein sequences (or None)
        bam_files: BAM 路径列表(可为空)|List of BAM paths (may be empty)
        output_dir: 输出目录|Output dir
        config: BrakerConfig|Config
        cmd_runner: 命令执行器|Command runner
        logger: 日志器|Logger

    Returns:
        refined genome 路径,失败/无证据返回 None|Refined genome path or None
    """
    logger.info("=" * 70)
    logger.info("证据驱动还原被误 mask 的真基因区|Evidence-based rescue of mis-masked regions")
    logger.info("=" * 70)

    masked_fa = os.path.abspath(masked_fa)
    output_dir = os.path.abspath(output_dir)
    if prot_seq:
        prot_seq = os.path.abspath(prot_seq)
    bam_files = [os.path.abspath(b) for b in bam_files]
    min_cds_len = config.rescue_min_cds_len
    min_identity = config.rescue_min_identity
    min_depth = config.rescue_min_depth

    if not prot_seq and not bam_files:
        logger.info("无蛋白/转录组证据,跳过 rescue|No protein/transcript evidence, skip rescue")
        return None

    os.makedirs(output_dir, exist_ok=True)
    regions = extract_masked_regions(masked_fa, min_cds_len)
    logger.info(f"提取 masked 区段|Masked regions: {len(regions)} 段(≥{min_cds_len}bp)")
    if not regions:
        logger.info("无 masked 区段,跳过 rescue|No masked regions, skip rescue")
        return None

    rescue_set = set()  # region index 集合|Set of region indices to rescue

    # 蛋白证据|Protein evidence (miniprot)
    if prot_seq:
        logger.info(f"蛋白证据: miniprot 比对 {os.path.basename(prot_seq)} → masked genome|"
                    f"Protein evidence via miniprot")
        paf_file = os.path.join(output_dir, "rescue_protein.paf")
        cmd = f"{config.miniprot_bin} -t {config.threads} --spliceflatten=no " \
              f"{masked_fa} {prot_seq} > {paf_file}"
        if cmd_runner.run_command(cmd, "miniprot 蛋白比对|miniprot protein alignment"):
            hits = parse_miniprot_paf(paf_file)
            logger.info(f"miniprot 比对数|miniprot hits: {len(hits)}")
            prot_rescued = 0
            for idx, (chrom, rs, re_) in enumerate(regions):
                for hchrom, hs, he, hid in hits:
                    if hchrom != chrom:
                        continue
                    overlap = min(re_, he) - max(rs, hs) + 1
                    if overlap >= min_cds_len and hid >= min_identity:
                        rescue_set.add(idx)
                        prot_rescued += 1
                        break
            logger.info(f"蛋白证据支持还原|Protein-supported rescue: {prot_rescued} 段")
        else:
            logger.warning("miniprot 失败,仅用转录组证据(若有)|miniprot failed")

    # 转录组证据|Transcript evidence (samtools depth)
    if bam_files:
        regions_bed = os.path.join(output_dir, "rescue_regions.bed")
        region_depth = compute_region_mean_depth(
            bam_files, regions, regions_bed, config.samtools_bin, cmd_runner, logger)
        rna_rescued = 0
        for idx, (chrom, rs, re_) in enumerate(regions):
            if region_depth.get(idx, 0.0) >= min_depth:
                if idx not in rescue_set:
                    rescue_set.add(idx)
                    rna_rescued += 1
        logger.info(f"转录组证据支持还原|Transcript-supported rescue: {rna_rescued} 段")

    if not rescue_set:
        logger.info("无证据支持还原的区域,refined genome 与 masked 一致|"
                    "No evidence-supported regions; refined = masked")
        return None

    # 写 rescued_regions.bed(0-based half-open)|Write rescued regions BED
    rescued_bed = os.path.join(output_dir, "rescued_regions.bed")
    rescue_regions = [regions[i] for i in sorted(rescue_set)]
    with open(rescued_bed, 'w') as f:
        for chrom, s, e in rescue_regions:
            f.write(f"{chrom}\t{s - 1}\t{e}\n")

    # unmask → refined genome|Unmask to refined genome
    refined_fa = os.path.join(output_dir, "refined_masked.fa")
    total_unmasked = write_refined_masked(masked_fa, rescue_regions, refined_fa)
    logger.info(f"还原完成|Rescue done: {len(rescue_regions)} 段, "
                f"unmask {total_unmasked} bp → {refined_fa}")
    logger.info(f"还原区域|Rescued regions: {rescued_bed}")
    return refined_fa
