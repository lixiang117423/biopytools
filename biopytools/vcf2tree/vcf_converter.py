"""
VCF转FASTA转换器|VCF to FASTA Converter

将VCF SNP数据转换为IUPAC编码的FASTA对齐矩阵
|Convert VCF SNP data to IUPAC-encoded FASTA alignment matrix
"""

import gzip
import os
from array import array


# IUPAC核酸模糊代码字典|IUPAC nucleotide ambiguity code dictionary
IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B",
    "ACGT": "N",
}


def get_iupac_code(ref: str, alt: str, genotype: str) -> str:
    """根据基因型返回IUPAC编码|Return IUPAC code based on genotype

    纯合子返回对应碱基，杂合子返回IUPAC模糊码，缺失返回N
    |Homozygous→single base, Heterozygous→IUPAC ambiguity, Missing→N

    Args:
        ref: 参考等位基因|Reference allele
        alt: 替代等位基因|Alternative allele
        genotype: 基因型（如"0/0", "0/1", "1/1", "./."）|Genotype format

    Returns:
        单字符IUPAC编码|Single IUPAC character
    """
    # 入口统一大写: IUPAC_CODES仅含大写键, 小写ref/alt会被误判为N(缺失)。
    # |Normalize to uppercase at entry: IUPAC_CODES has only uppercase keys; lowercase
    # ref/alt would otherwise be misclassified as N (missing).
    ref = ref.upper()
    alt = alt.upper()

    if genotype in ("./.", ".|."):
        return "N"

    sep = "/" if "/" in genotype else "|"
    parts = genotype.split(sep)

    try:
        allele_indices = [int(a) for a in parts]
    except (ValueError, TypeError):
        return "N"

    bases = []
    for idx in allele_indices:
        if idx == 0:
            bases.append(ref)
        elif idx == 1:
            bases.append(alt)
        else:
            return "N"

    # 排序去重以获得一致的IUPAC键|Sort+deduplicate for consistent IUPAC key
    bases_sorted = "".join(sorted(set(bases)))

    # 纯合|Homozygous（所有等位基因相同）
    if len(set(bases)) == 1:
        return bases[0]

    # 杂合→IUPAC编码|Heterozygous→IUPAC code
    return IUPAC_CODES.get(bases_sorted, "N")


class VcfToFastaConverter:
    """VCF转FASTA转换器|VCF to FASTA Converter"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def _open_vcf(self):
        """打开VCF文件（支持gzip）|Open VCF file (supports gzip)"""
        if self.config.input_file.endswith('.gz'):
            return gzip.open(self.config.input_file, 'rt')
        else:
            return open(self.config.input_file, 'r')

    def _extract_sample_names(self):
        """从VCF头部提取样本名称|Extract sample names from VCF header"""
        with self._open_vcf() as f:
            for line in f:
                line = line.strip()
                if line.startswith('#CHROM'):
                    fields = line.split('\t')
                    sample_names = fields[9:]
                    self.logger.info(
                        f"检测到样本数量|Detected sample count: {len(sample_names)}"
                    )
                    self.logger.info(
                        f"样本名称|Sample names: {', '.join(sample_names)}"
                    )
                    return sample_names
        raise ValueError(
            "VCF文件中未找到#CHROM头部行|#CHROM header line not found in VCF"
        )

    def convert(self) -> bool:
        """执行VCF→FASTA转换|Execute VCF→FASTA conversion

        Returns:
            bool: 转换是否成功|Whether conversion succeeded
        """
        self.logger.info("=" * 60)
        self.logger.info(
            "步骤1: VCF转FASTA对齐矩阵|Step 1: VCF to FASTA alignment matrix"
        )
        self.logger.info("=" * 60)

        # 提取样本名|Extract sample names
        sample_names = self._extract_sample_names()
        num_samples = len(sample_names)
        if num_samples == 0:
            raise ValueError("VCF中无样本|No samples in VCF")

        # 调整最小样本数|Adjust min sample count
        if self.config.min_samples_locus > num_samples:
            self.config.min_samples_locus = num_samples
            self.logger.info(
                f"最小样本数已调整为|Min samples adjusted to: {num_samples}"
            )

        # 初始化序列字典: 用array('u')存单字符, 避免list存单字符的指针开销。
        # 50样本×1M SNP时, list≈指针8B×N(且无紧凑布局), array('u')为紧凑4B/字符。
        # |Init sequence dict with array('u') for single chars, avoiding the per-element
        # pointer overhead of a list. At 50 samples x 1M SNPs a list holds ~8B/elem pointers
        # with no compact layout; array('u') stores 4B/char compactly.
        sequences = {name: array('u') for name in sample_names}

        # 统计|Counters
        total_snps = 0
        accepted_snps = 0
        skipped_multiallelic = 0
        skipped_missing = 0
        skipped_indel = 0

        # 处理每个变异位点|Process each variant site
        with self._open_vcf() as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 10:
                    continue

                ref = fields[3]
                alt = fields[4]

                total_snps += 1

                # 仅处理双等位SNP|Only biallelic SNPs
                if ',' in alt:
                    skipped_multiallelic += 1
                    continue

                # 仅处理单碱基SNP|Only single-base SNPs
                if len(ref) != 1 or len(alt) != 1:
                    skipped_indel += 1
                    continue

                # 获取基因型字段|Get genotype fields
                genotypes = fields[9:9 + num_samples]
                if len(genotypes) < num_samples:
                    continue

                # 检查最小样本数|Check min sample count
                called_count = sum(
                    1 for gt in genotypes
                    if gt.split(':')[0] not in ("./.", ".|.")
                )
                if called_count < self.config.min_samples_locus:
                    skipped_missing += 1
                    continue

                # 转换为IUPAC编码|Convert to IUPAC codes
                for i, name in enumerate(sample_names):
                    gt_field = genotypes[i].split(':')[0]
                    iupac = get_iupac_code(ref, alt, gt_field)
                    sequences[name].append(iupac)

                accepted_snps += 1

        # 检查是否有通过筛选的SNP|Check if any SNPs passed
        if accepted_snps == 0:
            self.logger.error(
                "没有SNP通过筛选|No SNPs passed filtering criteria"
            )
            return False

        # 写入FASTA文件|Write FASTA file
        with open(self.config.snps_fa, 'w') as f:
            for name in sample_names:
                seq = sequences[name].tounicode()
                if len(seq) > 0:
                    f.write(f">{name}\n")
                    # 每行最多80字符|Max 80 chars per line
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i + 80] + '\n')

        # 输出统计|Output statistics
        self.logger.info(f"总变异位点|Total variant sites: {total_snps}")
        self.logger.info(f"通过筛选的SNP|Accepted SNPs: {accepted_snps}")
        self.logger.info(
            f"跳过多等位位点|Skipped multi-allelic: {skipped_multiallelic}"
        )
        self.logger.info(
            f"跳过InDel/非单碱基|Skipped indel / non-single-base: {skipped_indel}"
        )
        self.logger.info(
            f"跳过样本不足位点|Skipped insufficient samples: {skipped_missing}"
        )
        self.logger.info(f"FASTA对齐长度(bp)|FASTA alignment length (bp): {accepted_snps}")
        self.logger.info(f"FASTA文件|FASTA file: {self.config.snps_fa}")

        return True
