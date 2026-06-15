"""
单倍型提取与构建模块|Haplotype Extraction and Builder Module

基于geneHapR算法，从VCF文件提取基因单倍型|Extract gene haplotypes from VCF based on geneHapR algorithm
"""

from collections import OrderedDict
from typing import Dict, List, Optional


class AlleleEncoder:
    """等位基因编码器，实现geneHapR的allS编码系统|Allele encoder implementing geneHapR's allS encoding system

    将基因型编码为单字母(纯合)或A|T形式(杂合)|Encode genotypes as single letter (homo) or A|T form (hetero)
    """

    def __init__(self):
        bases = ["A", "C", "G", "T"]
        self._homo: Dict[str, str] = OrderedDict()
        for b in bases:
            self._homo[f"{b}/{b}"] = b

        self._hetero: Dict[str, str] = OrderedDict()
        for a in bases:
            for b in bases:
                if a != b:
                    self._hetero[f"{a}/{b}"] = f"{a}|{b}"

        self._unknown = {"./.": "N", "N/N": "N"}
        self._all: Dict[str, str] = {}
        self._rebuild_all()

    def _rebuild_all(self):
        """重建组合查找表|Rebuild combined lookup table"""
        self._all = {}
        self._all.update(self._homo)
        self._all.update(self._hetero)
        self._all.update(self._unknown)

    def update(self, ref: str, alt: str):
        """用REF/ALT更新编码系统|Update encoding system with REF/ALT alleles

        Args:
            ref: 参考等位基因|Reference allele
            alt: 变异等位基因(逗号分隔)|Alternative alleles (comma-separated)
        """
        alleles = [ref] + [a.strip() for a in alt.split(",")]
        for allele in alleles:
            key = f"{allele}/{allele}"
            if key not in self._homo:
                self._homo[key] = allele

        for a in alleles:
            for b in alleles:
                if a != b:
                    key = f"{a}/{b}"
                    if key not in self._hetero:
                        self._hetero[key] = f"{a}|{b}"
        self._rebuild_all()

    def encode(self, genotype: str) -> str:
        """编码基因型|Encode genotype

        Args:
            genotype: 等位基因格式的基因型|Allele-based genotype (e.g., "A/T", "N/N")

        Returns:
            str: 编码后的基因型|Encoded genotype (e.g., "A", "A|T", "N")
        """
        genotype = genotype.upper().replace("|", "/")
        return self._all.get(genotype, "N")


class HaplotypeBuilder:
    """单倍型构建器，实现geneHapR的vcf2hap算法|Haplotype builder implementing geneHapR's vcf2hap algorithm

    核心逻辑：将每个样本在所有变异位点的基因型拼接为字符串，
    基因型模式完全相同的样本归为同一单倍型，按频率排序分配H001、H002...编号
    """

    def __init__(self, logger):
        self.logger = logger

    def build_from_vcf(
        self,
        vcf_file: str,
        chrom: str,
        start: int,
        end: int,
        hetero_remove: bool = True,
        na_drop: bool = True,
        hap_prefix: str = "H",
        pad: int = 3,
    ) -> dict:
        """从VCF文件构建单倍型|Build haplotypes from VCF file

        Args:
            vcf_file: VCF文件路径|VCF file path
            chrom: 染色体|Chromosome
            start: 起始位置|Start position (1-based)
            end: 终止位置|End position (1-based)
            hetero_remove: 去除杂合位点|Remove heterozygous sites
            na_drop: 去除缺失位点|Remove missing sites
            hap_prefix: 单倍型ID前缀|Haplotype ID prefix
            pad: 单倍型ID位数|Haplotype ID padding digits

        Returns:
            dict: 包含result(hapResult格式)和summary(hapSummary格式)|Contains result and summary
        """
        import cyvcf2

        self.logger.info(f"开始提取单倍型|Starting haplotype extraction")
        self.logger.info(f"区域|Region: {chrom}:{start}-{end}")

        vcf = cyvcf2.VCF(vcf_file)
        sample_names = list(vcf.samples)
        n_samples = len(sample_names)
        self.logger.info(f"样本数|Samples: {n_samples}")

        # 收集区间内变异位点|Collect variants in region
        variants_info = []
        try:
            for variant in vcf(f"{chrom}:{start}-{end}"):
                gt_bases = variant.gt_bases.tolist()
                variants_info.append({
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": ",".join(variant.ALT),
                    "gt_bases": gt_bases,
                })
        except (AssertionError, OSError):
            # 无索引时回退到全文件遍历|Fallback to full scan when no index
            self.logger.info("无tabix索引，回退到全文件遍历|No tabix index, falling back to full scan")
            for variant in vcf:
                if variant.CHROM != chrom or variant.POS < start or variant.POS > end:
                    continue
                gt_bases = variant.gt_bases.tolist()
                variants_info.append({
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": ",".join(variant.ALT),
                    "gt_bases": gt_bases,
                })
        vcf.close()

        n_variants = len(variants_info)
        self.logger.info(f"变异位点数|Variants in region: {n_variants}")

        if not variants_info:
            self.logger.warning("指定区域未找到变异|No variants found in specified region")
            return {}

        # 初始化并更新编码器|Initialize and update encoder
        encoder = AlleleEncoder()
        for vi in variants_info:
            encoder.update(vi["ref"], vi["alt"])

        # 构建基因型矩阵: 样本 x 位点|Build genotype matrix: samples x positions
        geno_matrix = []
        for si in range(n_samples):
            geno_row = []
            for vi in variants_info:
                gt = vi["gt_bases"][si]
                if gt is None or gt in (".", "./."):
                    gt = "N/N"
                gt = gt.replace("|", "/").upper()
                geno_row.append(encoder.encode(gt))
            geno_matrix.append(geno_row)

        # 去除杂合位点(含"|"的列)|Remove heterozygous columns (columns containing "|")
        if hetero_remove:
            cols_to_remove = set()
            for col in range(n_variants):
                if any("|" in geno_matrix[row][col] for row in range(n_samples)):
                    cols_to_remove.add(col)
            if cols_to_remove:
                self.logger.info(f"去除杂合位点|Removing heterozygous sites: {len(cols_to_remove)}")
                geno_matrix = [
                    [g for c, g in enumerate(row) if c not in cols_to_remove]
                    for row in geno_matrix
                ]
                variants_info = [
                    vi for c, vi in enumerate(variants_info) if c not in cols_to_remove
                ]
                n_variants = len(variants_info)

        # 去除缺失位点(含"N"的列)|Remove missing columns (columns containing "N")
        if na_drop:
            cols_to_remove = set()
            for col in range(n_variants):
                if any(geno_matrix[row][col] == "N" for row in range(n_samples)):
                    cols_to_remove.add(col)
            if cols_to_remove:
                self.logger.info(f"去除缺失位点|Removing missing sites: {len(cols_to_remove)}")
                geno_matrix = [
                    [g for c, g in enumerate(row) if c not in cols_to_remove]
                    for row in geno_matrix
                ]
                variants_info = [
                    vi for c, vi in enumerate(variants_info) if c not in cols_to_remove
                ]
                n_variants = len(variants_info)

        self.logger.info(f"过滤后变异位点数|Variants after filtering: {n_variants}")

        if n_variants == 0:
            self.logger.warning("过滤后无变异位点|No variants remaining after filtering")
            return {}

        # 按基因型模式分组|Group by concatenated genotype pattern
        patterns: Dict[str, List[int]] = {}
        for si in range(n_samples):
            pattern = "".join(geno_matrix[si])
            if pattern not in patterns:
                patterns[pattern] = []
            patterns[pattern].append(si)

        # 按频率排序(最多的在前)，频率相同时按R locale字典序|Sort by freq desc, R-locale alpha for ties
        _r_sort_table = str.maketrans("/|", "\x01\x02")

        def _r_sort_key(pattern: str) -> str:
            return pattern.translate(_r_sort_table)

        sorted_patterns = sorted(patterns.items(), key=lambda x: (-len(x[1]), _r_sort_key(x[0])))

        n_haps = len(sorted_patterns)
        self.logger.info(f"单倍型数|Unique haplotypes: {n_haps}")

        # 分配单倍型ID|Assign haplotype IDs
        hap_id_map: Dict[str, str] = {}
        for i, (pattern, _) in enumerate(sorted_patterns):
            hap_id_map[pattern] = f"{hap_prefix}{str(i + 1).zfill(pad)}"

        # 构建hapResult: 每行=一个样本|Build hapResult: each row = one sample
        result_rows = []
        for si in range(n_samples):
            pattern = "".join(geno_matrix[si])
            hap_id = hap_id_map[pattern]
            row = [hap_id] + geno_matrix[si] + [sample_names[si]]
            result_rows.append(row)
        result_rows.sort(key=lambda x: x[0])

        # 构建hapSummary: 每行=一个唯一单倍型|Build hapSummary: each row = one unique haplotype
        summary_rows = []
        for pattern, sample_indices in sorted_patterns:
            hap_id = hap_id_map[pattern]
            geno = geno_matrix[sample_indices[0]]
            acc_str = ";".join(sample_names[si] for si in sample_indices)
            freq = len(sample_indices)
            row = [hap_id] + geno + [acc_str, freq]
            summary_rows.append(row)

        # 元数据行|Metadata rows (geneHapR格式: 含Accession列)|Metadata rows with Accession column
        chr_row = ["CHR"] + [vi["chrom"] for vi in variants_info] + [f"Haplotypes: {n_haps}"]
        pos_row = ["POS"] + [str(vi["pos"]) for vi in variants_info] + [f"Individuals: {n_samples}"]
        info_row = ["INFO"] + ["."] * n_variants + [f"Variants: {n_variants}"]
        allele_row = ["ALLELE"] + [f"{vi['ref']}/{vi['alt']}" for vi in variants_info] + ["Accession"]

        final_result = [chr_row, pos_row, info_row, allele_row] + result_rows
        # hapSummary metadata has extra freq column|hapSummary metadata has extra freq column
        chr_row_s = chr_row + [""]
        pos_row_s = pos_row + [str(n_haps)]
        info_row_s = info_row + ["."]
        allele_row_s = allele_row + ["freq"]
        final_summary = [chr_row_s, pos_row_s, info_row_s, allele_row_s] + summary_rows

        self.logger.info(f"单倍型提取完成|Haplotype extraction completed")

        # 样本-单倍型映射: [sample_name, hap_id]|Sample-haplotype mapping
        sample_hap = [[sample_names[si], hap_id_map["".join(geno_matrix[si])]] for si in range(n_samples)]

        return {
            "result": final_result,
            "summary": final_summary,
            "sample_hap": sample_hap,
            "n_variants": n_variants,
            "n_haplotypes": n_haps,
            "n_samples": n_samples,
        }

    @staticmethod
    def export_table(data: list, output_file: str):
        """导出表格到文件|Export table to file"""
        with open(output_file, "w") as f:
            for row in data:
                f.write("\t".join(str(x) for x in row) + "\n")

    @staticmethod
    def export_xlsx(data: list, output_file: str):
        """导出表格到Excel文件|Export table to Excel file"""
        import openpyxl

        wb = openpyxl.Workbook()
        ws = wb.active
        for row in data:
            ws.append([str(x) for x in row])
        wb.save(output_file)
