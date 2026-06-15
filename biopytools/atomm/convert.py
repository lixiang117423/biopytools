"""ATOMM格式转换模块|ATOMM Format Conversion Module

支持从VCF和交叉感染矩阵表型转换为ATOMM输入格式
Supports conversion from VCF and cross-infection phenotype matrix to ATOMM input format

不引入外部依赖，仅使用Python标准库
No external dependencies, Python standard library only
"""

import csv
import gzip
import logging
from pathlib import Path
import numpy as np

logger = logging.getLogger(__name__)


def _open_vcf(path):
    """根据后缀选择打开方式|Open VCF file with gzip support"""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


# VCF固定列索引
_VCF_CHROM = 0
_VCF_POS = 1
_VCF_ID = 2
_VCF_REF = 3
_VCF_ALT = 4
_VCF_QUAL = 5
_VCF_FILTER = 6
_VCF_FORMAT = 8
_VCF_FIRST_SAMPLE = 9


def _parse_gt_field(gt_field, ref_allele=None):
    """解析单个GT字段|Parse single GT field

    VCF FORMAT字段格式: GT:AD:DP:GQ:PL (如 "0/0:1,0:1:3:0,3,29")
    本函数提取GT部分并返回REF allele计数(0/1编码)

    Args:
        gt_field: 完整FORMAT字段或纯GT字符串|Full FORMAT field or GT string
        ref_allele: 未使用, 保留兼容|Unused, kept for compatibility

    Returns:
        float: REF allele计数 (0=homALT, 1=het, 2=homREF)|REF allele count
    """
    # 提取GT部分(FORMAT字段中冒号前的第一个子字段)
    gt = gt_field.split(':')[0].split('|')[0]

    if gt in ('.', './.', '', '.'):
        return np.nan

    alleles = gt.split('/')
    ref_count = 0

    for allele in alleles:
        if allele == '.':
            return np.nan
        if allele == '0':
            ref_count += 1

    return ref_count


def vcf_to_genotype(vcf_path, output_path, maf_threshold=0.05, encoding='auto',
                     max_snps=None, chrom_prefix_filter=None):
    """将VCF文件转换为ATOMM基因型格式|Convert VCF to ATOMM genotype format

    Args:
        vcf_path: VCF文件路径|VCF file path
        output_path: 输出文件路径|Output file path
        maf_threshold: MAF过滤阈值|MAF filter threshold
        encoding: 编码方式 "haploid"|"dosage"|"auto"
        max_snps: 最大SNP数(0=无限制)|Max SNPs to keep (0=no limit)
        chrom_prefix_filter: 只保留指定染色体前缀|Only keep chromosomes with this prefix
    """
    logger.info(f"读取VCF文件|Reading VCF file: {vcf_path}")

    chrom_ids = []
    positions = []
    snp_records = []
    sample_names = []

    with _open_vcf(vcf_path) as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            if row[0].startswith('#'):
                if row[0].startswith('#CHROM'):
                    sample_names = row[_VCF_FIRST_SAMPLE:]
                continue

            chrom = row[_VCF_CHROM]
            pos = row[_VCF_POS]
            ref = row[_VCF_REF]
            alt = row[_VCF_ALT]
            gt_fields = row[_VCF_FIRST_SAMPLE:] if len(row) > _VCF_FIRST_SAMPLE else []

            if len(alt.split(',')) != 1:
                continue

            if chrom_prefix_filter and not chrom.startswith(chrom_prefix_filter):
                continue

            genotypes = []
            valid = True
            for gt in gt_fields:
                g = _parse_gt_field(gt, ref)
                if np.isnan(g):
                    valid = False
                    break
                genotypes.append(g)

            if not valid:
                continue

            gt_array = np.array(genotypes)
            maf = np.mean(gt_array) / 2.0
            if encoding == 'haploid':
                maf = np.mean(gt_array)

            if maf < maf_threshold or maf > (1 - maf_threshold):
                continue

            if encoding == 'auto':
                is_haploid = np.all((gt_array == 0) | (gt_array == 1))
                if is_haploid:
                    gt_array = gt_array
                else:
                    gt_array = gt_array / 2.0

            chrom_ids.append(chrom)
            positions.append(pos)
            snp_records.append(gt_array)

            if max_snps and len(snp_records) >= max_snps:
                break

    if not snp_records:
        logger.error("VCF中未找到符合条件的SNP|No valid SNPs found in VCF")
        return

    genotype_matrix = np.array(snp_records)
    n_snps, n_ind = genotype_matrix.shape

    logger.info(f"VCF解析完成|VCF parsed: {n_snps} SNPs, {n_ind} samples")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        for i in range(n_snps):
            chrom_str = chrom_ids[i]
            pos = positions[i]
            gt_str = '\t'.join(f"{g:.6f}" for g in genotype_matrix[i])
            f.write(f"{chrom_str}\t{pos}\t{gt_str}\n")

    logger.info(f"基因型文件已保存|Genotype file saved: {output_path} ({n_snps} SNPs x {n_ind} samples)")

    return genotype_matrix, chrom_ids, positions, sample_names


def phenotype_matrix_to_atomm(pheno_path, output_path, missing_value='NA'):
    """将交叉感染矩阵转换为ATOMM表型格式|Convert cross-infection matrix to ATOMM phenotype format

    交叉感染矩阵格式 (CSV/TSV):
        - 第一行第一列为header或空，第一行其余列为病原ID
        - 第一列其余行为宿主ID
        - 单元格为表型值

    输出ATOMM表型格式 (长格式):
        host_id  pathogen_id  intercept  phenotype

    Args:
        pheno_path: 交叉感染矩阵文件路径|Cross-infection matrix file path
        output_path: 输出文件路径|Output file path
        missing_value: 缺失值标记|Missing value marker (e.g., 'NA', 'nan', '.', '-')
    """
    logger.info(f"读取表型矩阵|Reading phenotype matrix: {pheno_path}")

    with open(pheno_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        header = next(reader)

        pathogen_ids_raw = header[1:]
        pathogen_ids = []
        for pid in pathogen_ids_raw:
            pid = pid.strip()
            if pid and pid != '':
                pathogen_ids.append(pid)

        host_ids = []
        phenotype_rows = []

        for row in reader:
            host_id = row[0].strip()
            if not host_id:
                continue

            host_ids.append(host_id)
            pheno_values = []
            for val in row[1:]:
                val = val.strip()
                if val in (missing_value, 'nan', 'NaN', '.', '-', 'NA'):
                    pheno_values.append(np.nan)
                else:
                    try:
                        pheno_values.append(float(val))
                    except ValueError:
                        pheno_values.append(np.nan)
            phenotype_rows.append(pheno_values)

    n_hosts = len(host_ids)
    n_pathogens = len(pathogen_ids)

    if len(phenotype_rows) != n_hosts:
        logger.warning(
            f"表型矩阵行列数不一致|Phenotype matrix row count mismatch: "
            f"宿主数={n_hosts}, 表型行数={len(phenotype_rows)}"
        )
        n_hosts = min(n_hosts, len(phenotype_rows))

    logger.info(f"表型矩阵: {n_hosts} 宿主 x {n_pathogens} 病原|"
                 f"Phenotype matrix: {n_hosts} hosts x {n_pathogens} pathogens")

    host_id_map = {name: i + 1 for i, name in enumerate(host_ids)}
    pathogen_id_map = {name: i + 1 for i, name in enumerate(pathogen_ids)}

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_missing = 0

    with open(output_path, 'w') as f:
        for i in range(n_hosts):
            for j in range(min(n_pathogens, len(phenotype_rows[i]))):
                pheno = phenotype_rows[i][j]

                if np.isnan(pheno):
                    n_missing += 1
                    continue

                h_id = host_id_map[host_ids[i]]
                p_id = pathogen_id_map.get(pathogen_ids[j], j + 1)

                f.write(f"{h_id}\t{p_id}\t1\t{pheno}\n")
                n_written += 1

    logger.info(f"表型文件已保存|Phenotype file saved: {output_path} ({n_written} pairs, {n_missing} missing skipped)")

    return host_ids, pathogen_ids, phenotype_rows


def detect_vcf_encoding(vcf_path, max_records=1000):
    """自动检测VCF的基因型编码方式|Auto-detect VCF genotype encoding

    Args:
        vcf_path: VCF文件路径|VCF file path
        max_records: 最大检测记录数|Max records to check

    Returns:
        str: "haploid" or "dosage"
    """
    is_haploid_count = 0
    total_count = 0

    with _open_vcf(vcf_path) as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            if row[0].startswith('#'):
                continue

            if total_count >= max_records:
                break

            ref = row[_VCF_REF]
            gt_fields = row[_VCF_FIRST_SAMPLE:] if len(row) > _VCF_FIRST_SAMPLE else []

            for gt in gt_fields:
                gt = gt.split('|')[0]
                alleles = gt.split('/')

                if gt in ('.', './.', '', ''):
                    continue

                is_diploid_het = any(a != '.' and a != ref for a in alleles)
                if not is_diploid_het:
                    is_haploid_count += 1

                total_count += 1

    if total_count == 0:
        logger.warning("无法检测编码方式，默认使用dosage|Cannot detect encoding, defaulting to dosage")
        return "dosage"

    haploid_ratio = is_haploid_count / total_count
    if haploid_ratio > 0.95:
        logger.info(f"自动检测编码方式|Auto-detected encoding: haploid (ratio={haploid_ratio:.3f})")
        return "haploid"
    else:
        logger.info(f"自动检测编码方式|Auto-detected encoding: dosage (haploid ratio={haploid_ratio:.3f})")
        return "dosage"
