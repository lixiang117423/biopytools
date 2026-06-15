"""候选基因StringTie组装与gffcompare模块|Candidate Gene StringTie Assembly and GFFCompare Module"""

import os
import logging
from typing import Dict, List, Optional

from .config import GeneRnaseqCheckConfig
from .utils import CommandRunner, build_conda_command


class StringTieRunner:
    """StringTie转录本组装器|StringTie Transcript Assembler"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def assemble_sample(self, sample_name: str, bam_file: str) -> Optional[str]:
        """用StringTie组装单个样本|Assemble a single sample with StringTie

        使用 -G 参数引导组装（非 -e 定量模式），允许发现新的转录本
        """
        out_dir = os.path.join(self.config.output_dir, '04_stringtie')
        os.makedirs(out_dir, exist_ok=True)

        gtf_file = os.path.join(out_dir, f"{sample_name}.gtf")

        if os.path.exists(gtf_file) and not self.config.force:
            self.logger.info(f"StringTie结果已存在，跳过|StringTie result exists, skipping: {gtf_file}")
            return gtf_file

        args = [
            bam_file,
            '-G', self.config.annotation_gff,
            '-o', gtf_file,
            '-p', str(self.config.threads),
            '-l', sample_name,
            '-A', os.path.join(out_dir, f"{sample_name}.abundance.tsv"),
        ]

        # 链特异性参数|Strandness parameter
        if self.config.strandness == 'RF':
            args.append('--rf')
        elif self.config.strandness == 'FR':
            args.append('--fr')

        cmd = build_conda_command(self.config.stringtie_path, args)
        success, stdout, stderr = self.cmd_runner.run(
            cmd, f"StringTie组装|StringTie assembly: {sample_name}",
            timeout=self.config.sample_timeout,
        )
        return gtf_file if success else None

    def assemble_all(self, samples: list, bam_files: List[str]) -> Optional[str]:
        """组装所有样本并合并|Assemble all samples and merge

        Returns:
            合并后的GTF文件路径|Merged GTF file path
        """
        out_dir = os.path.join(self.config.output_dir, '04_stringtie')
        os.makedirs(out_dir, exist_ok=True)

        gtf_files = []
        for sample, bam_file in zip(samples, bam_files):
            gtf = self.assemble_sample(sample['name'], bam_file)
            if gtf and os.path.exists(gtf):
                gtf_files.append(gtf)
            else:
                self.logger.warning(
                    f"样本 {sample['name']} StringTie组装失败|StringTie failed: {sample['name']}"
                )

        if not gtf_files:
            self.logger.error("所有样本StringTie组装失败|StringTie failed for all samples")
            return None

        # 多样本合并|Merge if multiple samples
        if len(gtf_files) == 1:
            return gtf_files[0]

        merged_gtf = os.path.join(out_dir, 'merged.gtf')
        merge_list = os.path.join(out_dir, 'merge_list.txt')

        if os.path.exists(merged_gtf) and not self.config.force:
            self.logger.info(f"StringTie合并结果已存在|StringTie merge exists: {merged_gtf}")
            return merged_gtf

        with open(merge_list, 'w') as f:
            for gtf in gtf_files:
                f.write(gtf + '\n')

        args = [
            '--merge', '-G', self.config.annotation_gff,
            '-o', merged_gtf, '-p', str(self.config.threads),
        ]
        args.extend(gtf_files)

        if self.config.strandness == 'RF':
            args.append('--rf')
        elif self.config.strandness == 'FR':
            args.append('--fr')

        cmd = build_conda_command(self.config.stringtie_path, args)
        success, stdout, stderr = self.cmd_runner.run(
            cmd, "StringTie合并|StringTie merge",
            timeout=self.config.sample_timeout,
        )
        return merged_gtf if success else None


class GFFCompareRunner:
    """GFFCompare运行器|GFFCompare Runner"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self, query_gtf: str) -> Optional[str]:
        """运行gffcompare|Run gffcompare

        Args:
            query_gtf: StringTie组装的GTF（或merged GTF）|Query GTF

        Returns:
            gffcompare输出前缀|gffcompare output prefix
        """
        out_dir = os.path.join(self.config.output_dir, '04_stringtie', 'gffcompare')
        os.makedirs(out_dir, exist_ok=True)

        prefix = os.path.join(out_dir, 'gffcmp')

        if os.path.exists(f"{prefix}.annotated.gtf") and not self.config.force:
            self.logger.info(f"gffcompare结果已存在|gffcompare result exists: {prefix}")
            return prefix

        args = [
            '-r', self.config.annotation_gff,
            '-o', prefix,
            query_gtf,
        ]

        cmd = build_conda_command(self.config.gffcompare_path, args)
        success, stdout, stderr = self.cmd_runner.run(cmd, "gffcompare比较|gffcompare comparison")

        if success:
            return prefix

        self.logger.error(f"gffcompare失败|gffcompare failed: {stderr[:300]}")
        return None


def extract_effector_class_codes(
    gffcompare_prefix: str,
    gene_ids: set,
    logger: logging.Logger,
) -> Dict[str, str]:
    """从gffcompare结果提取每个目标基因的class_code|Extract class_code per target gene

    解析 .tracking 文件和 .annotated.gtf

    Args:
        gffcompare_prefix: gffcompare输出前缀|gffcompare output prefix
        gene_ids: 目标基因ID集合|Target gene ID set
        logger: 日志器|Logger

    Returns:
        基因ID到class_code的映射|Gene ID to class_code mapping
    """
    result = {gid: 'u' for gid in gene_ids}

    # 从.tracking文件提取映射|.tracking maps query_id -> ref_id -> class_code
    tracking_file = f"{gffcompare_prefix}.tracking"
    annotated_file = f"{gffcompare_prefix}.annotated.gtf"

    if not os.path.exists(tracking_file):
        logger.warning(f"gffcompare tracking文件不存在|tracking file not found: {tracking_file}")
        return result

    # 读取tracking文件|Read tracking file
    # 格式: query_id  class_code  ref_id  ...  (列以 tab 分隔)
    query_to_ref = {}
    query_to_class = {}

    with open(tracking_file) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue
            query_id = cols[0]
            class_code = cols[1] if len(cols) > 1 else 'u'

            # ref_id 在后面的列中|ref_id in later columns
            ref_id = ''
            for col in cols[2:]:
                if col and col != '-' and col != '.':
                    ref_id = col.split('|')[0] if '|' in col else col
                    break

            query_to_ref[query_id] = ref_id
            query_to_class[query_id] = class_code

    # 从annotated GTF建立transcript -> gene映射|Build transcript -> gene from annotated GTF
    transcript_to_gene = {}
    if os.path.exists(annotated_file):
        with open(annotated_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 9:
                    continue
                # GFF3 attributes or GTF attributes
                attrs_str = cols[8]
                gene_id = ''
                transcript_id = ''
                for item in attrs_str.split(';'):
                    item = item.strip()
                    if not item:
                        continue
                    if '=' in item:
                        key, val = item.split('=', 1)
                        if key == 'gene_id':
                            gene_id = val.strip().strip('"')
                        elif key == 'transcript_id':
                            transcript_id = val.strip().strip('"')
                        elif key == 'oId':
                            transcript_id = val.strip()

                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id

    # 为每个目标基因确定class_code|Determine class_code for each target gene
    for query_id, class_code in query_to_class.items():
        ref_id = query_to_ref.get(query_id, '')
        gene_id = ref_id or transcript_to_gene.get(query_id, '')

        if gene_id in gene_ids:
            # 如果已经有非'u'的code，保留|Keep existing non-'u' code
            if result.get(gene_id) == 'u' or class_code in ('=', 'c', 'j'):
                result[gene_id] = class_code

    # 统计|Statistics
    from collections import Counter
    code_counts = Counter(result.values())
    logger.info(
        f"gffcompare class_code 分布|class_code distribution: "
        + ' '.join(f"{code}={cnt}" for code, cnt in sorted(code_counts.items()))
    )

    return result
