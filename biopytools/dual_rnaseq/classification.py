"""
互作转录组Reads物种分类模块|Dual RNA-seq Reads Species Classification Module

这是互作转录组分析的核心模块，负责将reads分类到不同物种。
This is the core module for dual RNA-seq analysis, responsible for classifying reads to different species.

分类策略|Classification Strategy:
  - 每个read同时比对到两个参考基因组|Each read is aligned to both reference genomes
  - 选择唯一比对且MAPQ最高的结果|Select uniquely mapped read with highest MAPQ
  - 如果两个基因组都有高质量比对，保留到ambiguous类别|If both genomes have high-quality alignments, keep as ambiguous
"""

import os
import pysam
from .utils import CommandRunner, FileValidator


class ReadsClassifier:
    """Reads物种分类器|Reads Species Classifier"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)

    def align_to_reference(self, fastq1: str, fastq2: str, index_prefix: str,
                          output_bam: str, species_name: str) -> bool:
        """将reads比对到指定参考基因组|Align reads to specified reference genome"""
        threads = self.config.threads

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(output_bam)
        os.makedirs(output_dir, exist_ok=True)

        # 检查BAM文件是否已存在|Check if BAM file already exists
        if self.file_validator.check_file_exists(output_bam, f"BAM文件|BAM file"):
            return True

        # HISAT2比对命令|HISAT2 alignment command
        cmd = (
            f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} "
            f"-p {threads} | "
            f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
        )

        self.logger.info(f"比对reads到{species_name}|Aligning reads to {species_name}")
        success = self.cmd_runner.run(cmd, f"HISAT2比对到{species_name}|HISAT2 alignment to {species_name}")

        if success:
            # 建立索引|Build index
            self.logger.info(f"建立BAM索引|Building BAM index: {output_bam}")
            samtools_cmd = f"samtools index {output_bam}"
            self.cmd_runner.run(samtools_cmd, f"建立{species_name} BAM索引|Building {species_name} BAM index")

        return success

    def classify_reads_by_mapq(self, species1_bam: str, species2_bam: str,
                               sample_name: str, output_dir: str) -> dict:
        """根据MAPQ对reads进行分类|Classify reads based on MAPQ

        分类规则|Classification rules:
          1. 仅在物种1比对上 -> species1
          2. 仅在物种2比对上 -> species2
          3. 两个都比对上，选择MAPQ更高的 -> 相应物种
          4. 两个都比对上，MAPQ相同 -> ambiguous
          5. 两个都没比对上 -> unassigned
        """
        self.logger.info(f"开始分类样本reads|Starting to classify sample reads: {sample_name}")

        # 输出文件路径|Output file paths
        species1_out = os.path.join(output_dir, f"{sample_name}.{self.config.species1_name}.bam")
        species2_out = os.path.join(output_dir, f"{sample_name}.{self.config.species2_name}.bam")
        ambiguous_out = os.path.join(output_dir, f"{sample_name}.ambiguous.bam")
        unassigned_out = os.path.join(output_dir, f"{sample_name}.unassigned.bam")

        # 统计信息|Statistics
        stats = {
            "species1": 0,
            "species2": 0,
            "ambiguous": 0,
            "unassigned": 0,
            "total": 0
        }

        # 打开BAM文件|Open BAM files
        try:
            bam1 = pysam.AlignmentFile(species1_bam, "rb")
            bam2 = pysam.AlignmentFile(species2_bam, "rb")

            # 创建输出BAM文件|Create output BAM files
            out1 = pysam.AlignmentFile(species1_out, "wb", template=bam1)
            out2 = pysam.AlignmentFile(species2_out, "wb", template=bam2)
            out_ambiguous = pysam.AlignmentFile(ambiguous_out, "wb", template=bam1)
            out_unassigned = pysam.AlignmentFile(unassigned_out, "wb", template=bam1)

            # 读取所有reads并按read name分组|Read all reads and group by read name
            self.logger.info("读取比对结果|Reading alignment results")

            # 使用字典存储reads|Use dictionary to store reads
            reads1 = {}
            reads2 = {}

            for read in bam1.fetch():
                if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    read_name = read.query_name
                    mapq = read.mapping_quality
                    reads1[read_name] = (read, mapq)

            for read in bam2.fetch():
                if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    read_name = read.query_name
                    mapq = read.mapping_quality
                    reads2[read_name] = (read, mapq)

            stats["total"] = len(set(list(reads1.keys()) + list(reads2.keys())))

            # 分类reads|Classify reads
            self.logger.info(f"总共|Total {stats['total']}个reads需要分类|reads to classify")

            all_read_names = set(reads1.keys()) | set(reads2.keys())

            for read_name in all_read_names:
                in_species1 = read_name in reads1
                in_species2 = read_name in reads2

                if not in_species1 and not in_species2:
                    # 两个都没比对上|Not mapped to either
                    stats["unassigned"] += 1
                    continue

                if in_species1 and not in_species2:
                    # 仅在物种1|Only in species 1
                    read, mapq = reads1[read_name]
                    if self._is_high_quality(mapq):
                        out1.write(read)
                        stats["species1"] += 1

                elif in_species2 and not in_species1:
                    # 仅在物种2|Only in species 2
                    read, mapq = reads2[read_name]
                    if self._is_high_quality(mapq):
                        out2.write(read)
                        stats["species2"] += 1

                else:
                    # 两个都比对上了|Mapped to both
                    read1, mapq1 = reads1[read_name]
                    read2, mapq2 = reads2[read_name]

                    if not self._is_high_quality(mapq1) and not self._is_high_quality(mapq2):
                        # 两个都是低质量|Both low quality
                        stats["unassigned"] += 1
                    elif mapq1 > mapq2:
                        # 物种1质量更高|Species 1 has higher quality
                        out1.write(read1)
                        stats["species1"] += 1
                    elif mapq2 > mapq1:
                        # 物种2质量更高|Species 2 has higher quality
                        out2.write(read2)
                        stats["species2"] += 1
                    else:
                        # 质量相同，ambiguous|Same quality, ambiguous
                        if self._is_high_quality(mapq1):
                            out_ambiguous.write(read1)
                            stats["ambiguous"] += 1
                        else:
                            stats["unassigned"] += 1

            # 关闭文件|Close files
            bam1.close()
            bam2.close()
            out1.close()
            out2.close()
            out_ambiguous.close()
            out_unassigned.close()

            # 排序并建立索引|Sort and build indexes
            self.logger.info("排序分类结果BAM文件|Sorting classified BAM files")
            for bam_file in [species1_out, species2_out, ambiguous_out]:
                if os.path.exists(bam_file) and os.path.getsize(bam_file) > 0:
                    # 排序BAM文件|Sort BAM file
                    sorted_bam = bam_file.replace('.bam', '.sorted.bam')
                    sort_cmd = f"samtools sort -@ {self.config.threads} -O BAM -o {sorted_bam} {bam_file}"
                    sort_success = self.cmd_runner.run(sort_cmd, f"排序BAM文件|Sorting BAM file: {bam_file}")

                    if sort_success:
                        # 用排序后的文件替换原文件|Replace original file with sorted file
                        import shutil
                        shutil.move(sorted_bam, bam_file)

                        # 建立索引|Build index
                        self.logger.info(f"建立BAM索引|Building BAM index: {bam_file}")
                        samtools_cmd = f"samtools index {bam_file}"
                        self.cmd_runner.run(samtools_cmd, f"建立索引|Building index: {bam_file}")
                    else:
                        self.logger.warning(f"BAM文件排序失败，跳过索引|BAM sorting failed, skipping index: {bam_file}")

            # 输出统计信息|Output statistics
            self.logger.info(f"分类完成|Classification completed for sample: {sample_name}")
            self.logger.info(f"  总reads数|Total reads: {stats['total']}")
            self.logger.info(f"  {self.config.species1_name}|Species 1: {stats['species1']} ({stats['species1']/stats['total']*100:.2f}%)")
            self.logger.info(f"  {self.config.species2_name}|Species 2: {stats['species2']} ({stats['species2']/stats['total']*100:.2f}%)")
            self.logger.info(f"  ambiguous|Ambiguous: {stats['ambiguous']} ({stats['ambiguous']/stats['total']*100:.2f}%)")
            self.logger.info(f"  unassigned|Unassigned: {stats['unassigned']} ({stats['unassigned']/stats['total']*100:.2f}%)")

            return stats

        except Exception as e:
            self.logger.error(f"分类reads时出错|Error classifying reads: {e}")
            return None

    def _is_high_quality(self, mapq: int) -> bool:
        """判断MAPQ是否高质量|Check if MAPQ is high quality"""
        if self.config.unique_only:
            return mapq >= self.config.min_mapq
        else:
            return mapq >= 0

    def process_sample_classification(self, fastq1: str, fastq2: str,
                                      species1_index: str, species2_index: str,
                                      sample_name: str) -> bool:
        """处理单个样本的分类流程|Process classification workflow for a single sample"""
        self.logger.info(f"处理样本分类|Processing sample classification: {sample_name}")

        # 创建临时比对结果目录|Create temporary alignment directory
        temp_dir = os.path.join(self.config.output_dir, "02.classification", ".temp", sample_name)
        os.makedirs(temp_dir, exist_ok=True)

        # 临时BAM文件路径|Temporary BAM file paths
        species1_raw_bam = os.path.join(temp_dir, f"{sample_name}.{self.config.species1_name}.raw.bam")
        species2_raw_bam = os.path.join(temp_dir, f"{sample_name}.{self.config.species2_name}.raw.bam")

        # 输出分类结果目录|Output classification directory
        classification_dir = os.path.join(self.config.output_dir, "02.classification", sample_name)
        os.makedirs(classification_dir, exist_ok=True)

        # 步骤1: 比对到两个参考基因组|Step 1: Align to both reference genomes
        self.logger.info("步骤1: 双基因组比对|Step 1: Dual-genome alignment")

        if not self.align_to_reference(fastq1, fastq2, species1_index, species1_raw_bam, self.config.species1_name):
            self.logger.error(f"{self.config.species1_name}比对失败|{self.config.species1_name} alignment failed")
            return False

        if not self.align_to_reference(fastq1, fastq2, species2_index, species2_raw_bam, self.config.species2_name):
            self.logger.error(f"{self.config.species2_name}比对失败|{self.config.species2_name} alignment failed")
            return False

        # 步骤2: 根据MAPQ分类reads|Step 2: Classify reads based on MAPQ
        self.logger.info("步骤2: Reads物种分类|Step 2: Reads species classification")

        stats = self.classify_reads_by_mapq(species1_raw_bam, species2_raw_bam, sample_name, classification_dir)

        if stats is None:
            self.logger.error(f"样本分类失败|Sample classification failed: {sample_name}")
            return False

        self.logger.info(f"样本分类成功|Sample classification completed: {sample_name}")

        return True
