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
import sys
from typing import List, Dict, Optional, Union, Any
from pathlib import Path


class ExonicVariantProcessor:
    """外显子变异结果处理器|Exonic Variant Result Processor"""

    def __init__(self, logger):
        self.logger = logger

    def parse_annovar_line(self, line: str) -> Optional[Dict[str, str]]:
        """解析ANNOVAR外显子注释输出的每一行|Parse each line of ANNOVAR exonic annotation output"""
        fields = line.strip().split('\t')

        if len(fields) < 9:
            return None

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

        # 解析基因信息|Parse gene information
        # 格式|Format: nbisL1-mrna-1:G1:exon1:c.T5C:p.V2A,
        gene_match = re.match(r'([^:]+):([^:]+):([^:]+):c\.([^:]+):p\.([^,]+)', gene_info)

        if gene_match:
            transcript = gene_match.group(1)
            gene = gene_match.group(2)
            exon = gene_match.group(3)
            cdna_change = gene_match.group(4)  # DNA变化|DNA change
            protein_change = gene_match.group(5)  # 蛋白变化|Protein change
        else:
            transcript = gene = exon = cdna_change = protein_change = 'NA'

        # 解析DNA变化的位置和碱基|Parse DNA change positions and bases
        dna_pos_from = dna_pos_to = dna_ref = dna_alt = 'NA'

        # 处理不同类型的DNA变化|Handle different types of DNA changes
        if 'ins' in cdna_change:
            # 插入|Insertion: c.17_18insGAAG
            match = re.match(r'(\d+)_(\d+)ins([A-Z]+)', cdna_change)
            if match:
                dna_pos_from = match.group(1)
                dna_pos_to = match.group(2)
                dna_alt = match.group(3)
                dna_ref = '-'
        elif 'del' in cdna_change:
            # 删除|Deletion: c.20_23del
            match = re.match(r'(\d+)_(\d+)del', cdna_change)
            if match:
                dna_pos_from = match.group(1)
                dna_pos_to = match.group(2)
                dna_ref = 'deleted'
                dna_alt = '-'
        else:
            # SNV|SNV: c.T5C 或 c.C28T
            match = re.match(r'([A-Z])(\d+)([A-Z])', cdna_change)
            if match:
                dna_ref = match.group(1)
                dna_pos_from = dna_pos_to = match.group(2)
                dna_alt = match.group(3)

        # 解析蛋白变化|Parse protein change
        protein_pos = protein_ref = protein_alt = 'NA'

        if 'fs' in protein_change:
            # 移码突变|Frameshift: p.V6fs 或 p.L7fs
            match = re.match(r'([A-Z])(\d+)fs', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = 'frameshift'
        elif 'delins' in protein_change:
            # 复杂变化|Complex change: p.L17delinsLRKEVSRRLREST
            match = re.match(r'([A-Z]+)(\d+)delins([A-Z]+)', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = match.group(3)
        else:
            # 简单替换|Simple substitution: p.V2A 或 p.H10Y
            match = re.match(r'([A-Z])(\d+)([A-Z])', protein_change)
            if match:
                protein_ref = match.group(1)
                protein_pos = match.group(2)
                protein_alt = match.group(3)

        # 确定变异结果类型|Determine variant effect type
        if 'frameshift' in variant_type:
            effect = '移码突变'  # Frameshift
        elif 'nonsynonymous' in variant_type:
            effect = '错义突变'  # Nonsynonymous
        elif 'synonymous' in variant_type:
            effect = '同义突变'  # Synonymous
        elif 'stopgain' in variant_type:
            effect = '无义突变'  # Stop gain
        elif 'stoploss' in variant_type:
            effect = '终止丢失'  # Stop loss
        else:
            effect = variant_type

        return {
            'Line_ID': line_id,
            '染色体': chrom,
            '变异起始': start,
            '变异终止': end,
            '突变类型': mutation_type,
            '基因': gene,
            '转录本': transcript,
            '变异结果': effect,
            'DNA位置起': dna_pos_from,
            'DNA位置止': dna_pos_to,
            'DNA参考': dna_ref,
            'DNA变异': dna_alt,
            '蛋白位置': protein_pos,
            '蛋白参考': protein_ref,
            '蛋白变异': protein_alt
        }

    def process_exonic_variants(self, input_file: str, output_file: Optional[str] = None) -> List[Dict[str, str]]:
        """处理外显子变异注释结果|Process exonic variant annotation results"""
        self.logger.info(f" 处理外显子变异文件|Processing exonic variant file: {input_file}")

        results = []

        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('line'):  # 跳过表头|Skip header
                    parsed = self.parse_annovar_line(line)
                    if parsed:
                        results.append(parsed)

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
        headers = ['Line_ID', '染色体', '变异起始', '变异终止', '突变类型', '基因',
                   '转录本', '变异结果', 'DNA位置起', 'DNA位置止', 'DNA参考',
                   'DNA变异', '蛋白位置', '蛋白参考', '蛋白变异']
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
                except:
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

    def __init__(self, logger, output_dir: str):
        self.logger = logger
        self.output_dir = output_dir
        self.exonic_processor = ExonicVariantProcessor(logger)
        self.all_processor = AllVariantProcessor(logger)

    def process_exonic_results(self, exonic_file: str, output_prefix: Optional[str] = None) -> Optional[str]:
        """处理外显子注释结果|Process exonic annotation results"""
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