"""候选基因RNA-seq比对模块|Candidate Gene RNA-seq Alignment Module"""

import os
import logging
import re
from typing import List, Optional, Tuple

from .config import GeneRnaseqCheckConfig
from .utils import CommandRunner, build_conda_command
from .parse_gff import convert_gff3_to_gtf, build_gene_bed


class HISAT2Indexer:
    """HISAT2基因组索引构建器|HISAT2 Genome Index Builder"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_index(self) -> Optional[str]:
        """构建HISAT2索引|Build HISAT2 genome index

        Returns:
            索引前缀路径|Index prefix path, None on failure
        """
        output_dir = os.path.join(self.config.output_dir, '01_align')
        os.makedirs(output_dir, exist_ok=True)

        genome_name = os.path.splitext(os.path.basename(self.config.genome_fa))[0]
        index_prefix = os.path.join(output_dir, f"{genome_name}.hisat2")

        # 断点续传|Checkpoint resume
        if os.path.exists(f"{index_prefix}.1.ht2") and not self.config.force:
            self.logger.info(f"HISAT2索引已存在，跳过|HISAT2 index exists, skipping: {index_prefix}")
            return index_prefix

        # GFF3转GTF供extract脚本使用|Convert GFF3 to GTF for extract scripts
        gtf_file = os.path.join(output_dir, f"{genome_name}.gtf")
        if not os.path.exists(gtf_file):
            convert_gff3_to_gtf(self.config.annotation_gff, gtf_file, self.logger)

        # 提取splice sites和exons|Extract splice sites and exons
        ss_file = os.path.join(output_dir, f"{genome_name}.ss")
        exon_file = os.path.join(output_dir, f"{genome_name}.exon")

        ss_ok = self._extract_splice_sites(gtf_file, ss_file)
        exon_ok = self._extract_exons(gtf_file, exon_file)

        # 构建索引命令|Build index command
        args = ['-p', str(self.config.threads), self.config.genome_fa, index_prefix]
        if ss_ok:
            args.extend(['--ss', ss_file])
        if exon_ok:
            args.extend(['--exon', exon_file])

        cmd = build_conda_command(self.config.hisat2_build_path, args)
        success, stdout, stderr = self.cmd_runner.run(
            cmd, f"HISAT2索引构建|HISAT2 index building"
        )
        if success:
            return index_prefix
        return None

    def _extract_splice_sites(self, gtf_file: str, output_file: str) -> bool:
        """提取剪接位点|Extract splice sites"""
        if os.path.exists(output_file):
            return True
        cmd = build_conda_command(
            self.config.extract_splice_sites_path, [gtf_file]
        )
        success, stdout, stderr = self.cmd_runner.run(
            cmd, f"提取剪接位点|Extracting splice sites"
        )
        if success and stdout:
            with open(output_file, 'w') as f:
                f.write(stdout)
            return True
        self.logger.warning("剪接位点提取失败，不使用splice sites|Splice site extraction failed, skipping")
        return False

    def _extract_exons(self, gtf_file: str, output_file: str) -> bool:
        """提取外显子|Extract exons"""
        if os.path.exists(output_file):
            return True
        cmd = build_conda_command(
            self.config.extract_exons_path, [gtf_file]
        )
        success, stdout, stderr = self.cmd_runner.run(
            cmd, f"提取外显子|Extracting exons"
        )
        if success and stdout:
            with open(output_file, 'w') as f:
                f.write(stdout)
            return True
        self.logger.warning("外显子提取失败，不使用exon信息|Exon extraction failed, skipping")
        return False


class HISAT2Aligner:
    """HISAT2比对器|HISAT2 Aligner"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def align_sample(self, sample: dict, index_prefix: str,
                     strandness: str = '') -> Optional[str]:
        """比对单个样本|Align a single sample

        Args:
            sample: {'name': str, 'fastq1': str, 'fastq2': str}
            index_prefix: HISAT2索引前缀|HISAT2 index prefix
            strandness: 链特异性参数|Strandness parameter ('', 'FR', 'RF')

        Returns:
            BAM文件路径|BAM file path, None on failure
        """
        output_dir = os.path.join(self.config.output_dir, '01_align')
        bam_file = os.path.join(output_dir, f"{sample['name']}.sorted.bam")
        sam_file = os.path.join(output_dir, f"{sample['name']}.sam")

        # 断点续传|Checkpoint resume
        if os.path.exists(f"{bam_file}.bai") and not self.config.force:
            self.logger.info(f"BAM已存在，跳过|BAM exists, skipping: {bam_file}")
            return bam_file

        # Step 1: HISAT2 -> SAM (hisat2 in RNA_Seq env)
        hisat2_args = [
            '-x', index_prefix,
            '-1', sample['fastq1'],
            '-2', sample['fastq2'],
            '-p', str(self.config.threads),
            '-S', sam_file,
        ]
        if strandness == 'RF':
            hisat2_args.extend(['--rna-strandness', 'RF'])
        elif strandness == 'FR':
            hisat2_args.extend(['--rna-strandness', 'FR'])

        cmd_hisat2 = build_conda_command(self.config.hisat2_path, hisat2_args)
        success, stdout, stderr = self.cmd_runner.run(
            cmd_hisat2, f"HISAT2比对|HISAT2 alignment: {sample['name']}",
            timeout=self.config.sample_timeout,
        )
        if not success:
            return None

        # Step 2: samtools sort SAM -> BAM (samtools in GATK env)
        cmd_sort = build_conda_command(self.config.samtools_path, [
            'sort', '-@', str(self.config.threads),
            '-O', 'BAM', '-o', bam_file, sam_file,
        ])
        success, stdout, stderr = self.cmd_runner.run(
            cmd_sort, f"排序BAM|Sorting BAM: {sample['name']}",
            timeout=self.config.sample_timeout,
        )
        if not success:
            return None

        # Step 3: samtools index
        cmd_idx = build_conda_command(self.config.samtools_path, ['index', bam_file])
        success, stdout, stderr = self.cmd_runner.run(
            cmd_idx, f"BAM索引|BAM indexing: {sample['name']}",
            timeout=self.config.sample_timeout,
        )

        # Step 4: samtools flagstat
        flagstat_file = os.path.join(output_dir, f"{sample['name']}.flagstat")
        cmd_flag = build_conda_command(self.config.samtools_path, ['flagstat', bam_file])
        success_f, stdout_f, stderr_f = self.cmd_runner.run(
            cmd_flag, f"BAM flagstat|BAM flagstat: {sample['name']}",
        )
        if success_f and stdout_f:
            with open(flagstat_file, 'w') as f:
                f.write(stdout_f)

        # 清理SAM|Clean up SAM
        if os.path.exists(sam_file):
            os.remove(sam_file)

        return bam_file if os.path.exists(bam_file) else None

    def align_all_samples(self, samples: list, index_prefix: str,
                          strandness: str = '') -> List[str]:
        """比对所有样本|Align all samples"""
        bam_files = []
        for sample in samples:
            bam = self.align_sample(sample, index_prefix, strandness)
            if bam:
                bam_files.append(bam)
            else:
                self.logger.warning(f"样本比对失败|Sample alignment failed: {sample['name']}")
        return bam_files


class StrandnessDetector:
    """建库链特异性检测器|Library Strandness Detector"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def detect(self, bam_file: str,
               chrom_lengths: dict = None) -> Tuple[str, float]:
        """检测建库链特异性|Detect library strandness

        使用RSeQC infer_experiment.py

        Args:
            bam_file: BAM文件路径|BAM file path
            chrom_lengths: 染色体长度字典|Chromosome lengths dict

        Returns:
            (strandness, confidence): 'FR'/'RF'/'unstranded' 和置信度百分比
        """
        output_dir = os.path.join(self.config.output_dir, '01_align')
        os.makedirs(output_dir, exist_ok=True)

        decision_file = os.path.join(output_dir, 'strandness_decision.txt')

        if os.path.exists(decision_file) and not self.config.force:
            self.logger.info(f"链特异性检测结果已存在，跳过|Strandness result exists, skipping")
            with open(decision_file) as f:
                for line in f:
                    if line.startswith('strandness:'):
                        return line.strip().split(':\t')[1], 100.0
            return 'unstranded', 0.0

        # 生成gene BED（RSeQC需要BED格式输入）|Generate gene BED for RSeQC
        gene_bed = os.path.join(output_dir, 'rseqc_gene_regions.bed')
        if not os.path.exists(gene_bed):
            from .parse_gff import parse_gff3_effector_genes
            # 为RSeQC生成全基因组的gene BED（不仅仅是候选基因）
            self._build_genome_gene_bed(gene_bed)

        # 运行infer_experiment.py|Run infer_experiment.py
        cmd = build_conda_command(self.config.infer_experiment_path, [
            '-r', gene_bed,
            '-i', bam_file,
        ])
        success, stdout, stderr = self.cmd_runner.run(
            cmd, "检测链特异性|Detecting strandness (RSeQC)",
            timeout=600,
        )

        if not success:
            self.logger.warning("RSeQC strandness检测失败，使用unstranded|RSeQC failed, using unstranded")
            self._save_decision(decision_file, 'unstranded', 0.0)
            return 'unstranded', 0.0

        # 解析输出|Parse output
        strandness, confidence = self._parse_infer_output(stdout)
        self.logger.info(
            f"链特异性检测结果|Strandness result: {strandness} "
            f"(置信度|confidence: {confidence:.1f}%)"
        )
        self._save_decision(decision_file, strandness, confidence)
        return strandness, confidence

    def _build_genome_gene_bed(self, output_path: str):
        """从GFF3生成全基因组gene BED（用于RSeQC）|Build genome-wide gene BED from GFF3"""
        from .parse_gff import _parse_gff3_attributes

        with open(self.config.annotation_gff) as fin, open(output_path, 'w') as fout:
            for line in fin:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 9 or cols[2] != 'gene':
                    continue
                attrs = _parse_gff3_attributes(cols[8])
                gene_id = attrs.get('ID', '')
                if not gene_id:
                    continue
                chrom = cols[0]
                start = int(cols[3]) - 1  # GFF3 1-based -> BED 0-based
                end = int(cols[4])          # BED half-open, GFF3 1-based end is same
                strand = cols[6]
                fout.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")

    def _parse_infer_output(self, output: str) -> Tuple[str, float]:
        """解析RSeQC输出|Parse RSeQC infer_experiment output"""
        strandness = 'unstranded'
        confidence = 0.0

        for line in output.strip().split('\n'):
            line = line.strip()
            # RSeQC v5 输出格式|Output format varies by version
            if 'this is' in line.lower() and 'data' in line.lower():
                if 'fr' in line.lower() and 'unstranded' not in line.lower():
                    strandness = 'FR'
                elif 'rf' in line.lower() and 'unstranded' not in line.lower():
                    strandness = 'RF'
                elif 'unstranded' in line.lower():
                    strandness = 'unstranded'

            # 尝试提取置信度|Try to extract confidence
            if 'Fraction of reads' in line or 'non-' in line.lower():
                num_match = re.search(r'([\d.]+)%', line)
                if num_match:
                    confidence = float(num_match.group(1))
                    # 如果这行说的是"failed to determine"，反转置信度
                    if 'fail' in line.lower() or 'undetermined' in line.lower():
                        confidence = 100.0 - confidence

        # 如果无法从格式解析，尝试另一种模式|Alternative parsing
        if confidence == 0.0:
            for line in output.strip().split('\n'):
                match = re.search(r'([\d.]+)%\s*(?:of reads|fr|rf|reverse|forward)', line.lower())
                if match:
                    confidence = float(match.group(1))

        # 置信度太低则判定为unstranded|Low confidence = unstranded
        if confidence < self.config.strandness_confidence:
            strandness = 'unstranded'
            confidence = 100.0 - confidence

        return strandness, confidence

    def _save_decision(self, decision_file: str, strandness: str, confidence: float):
        """保存链特异性决定|Save strandness decision"""
        with open(decision_file, 'w') as f:
            f.write(f"strandness:\t{strandness}\n")
            f.write(f"confidence:\t{confidence:.1f}%\n")
            f.write(f"threshold:\t{self.config.strandness_confidence}%\n")
