"""Post-GWAS分析处理器|Post-GWAS Analysis Processor"""

import os
import random
import subprocess
from glob import glob
from ..utils import CommandRunner, format_number, build_conda_command


class PostGwasProcessor:
    """Post-GWAS分析处理器 - 将k-mer映射到参考基因组|Post-GWAS Analysis Processor - Map k-mers to reference genome"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行Post-GWAS分析|Run Post-GWAS analysis"""
        self.logger.info("开始Post-GWAS分析|Starting Post-GWAS analysis")

        # 检查必需参数|Check required parameters
        if not self.config.genome_file:
            self.logger.error("未指定参考基因组文件|Reference genome file not specified")
            return False

        if not os.path.exists(self.config.genome_file):
            self.logger.error(f"参考基因组文件不存在|Reference genome file does not exist: {self.config.genome_file}")
            return False

        # 检查输入目录|Check input directory
        asso_dir = self.config.dirs['association']
        ps_files = glob(os.path.join(asso_dir, '*.ps'))

        if not ps_files:
            self.logger.error(f"未找到关联分析结果文件|No association result files found in {asso_dir}")
            return False

        self.logger.info(f"找到|Found {len(ps_files)} 个关联结果文件|association result files")

        # 创建输出目录（根据比对工具命名）|Create output directory (named by alignment tool)
        if self.config.alignment_tool == 'bwa':
            output_dir_name = '07_post_gwas_bwa'
        else:  # blast
            output_dir_name = '07_post_gwas_blast'
        output_dir = os.path.join(self.config.output_dir, output_dir_name)
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"输出目录|Output directory: {output_dir}")
        self.logger.info(f"比对工具|Alignment tool: {self.config.alignment_tool}")

        # 检查是否需要跳过（断点续传）|Check if should skip (resume)
        final_output = os.path.join(output_dir, 'kmer_gwas_pval_1e-05.txt')
        if not self.config.force and os.path.exists(final_output):
            file_size = os.path.getsize(final_output)
            if file_size > 0:
                self.logger.info(f"输出文件已存在，跳过Post-GWAS分析|Output file exists, skipping Post-GWAS analysis: {final_output}")
                self.logger.info(f"如需重新运行，请使用 --force 参数|Use --force to re-run")
                return True

        # 合并所有.ps文件（带断点续传）|Merge all .ps files (with resume)
        merged_ps = os.path.join(output_dir, 'all_association.ps')
        if not self.config.force and os.path.exists(merged_ps):
            file_size = os.path.getsize(merged_ps)
            if file_size > 0:
                self.logger.info(f"合并的PS文件已存在，跳过合并|Merged PS file exists, skipping: {merged_ps}")
            else:
                self._merge_ps_files(ps_files, merged_ps)
        else:
            self._merge_ps_files(ps_files, merged_ps)

        # 统计总k-mer数|Count total k-mers
        total_kmers = self._count_kmers(merged_ps)
        self.logger.info(f"总k-mer数|Total k-mers: {format_number(total_kmers)}")

        # 计算阈值|Calculate thresholds
        threshold_strict = 1e-5  # 默认严格阈值|Default strict threshold
        threshold_bonferroni = 0.05 / total_kmers  # Bonferroni校正阈值|Bonferroni corrected threshold

        self.logger.info(f"阈值设置|Threshold settings:")
        self.logger.info(f"  严格阈值|Strict threshold: {threshold_strict}")
        self.logger.info(f"  Bonferroni校正阈值|Bonferroni threshold: {threshold_bonferroni}")

        # 创建所有k-mer的FASTA文件（带断点续传）|Create FASTA file for ALL k-mers (with resume)
        kmer_fasta = os.path.join(output_dir, 'all_kmers.fasta')
        kmer_mapping_file = kmer_fasta.replace('.fasta', '.mapping.txt')

        # 检查FASTA和映射文件是否都已存在|Check if both FASTA and mapping files exist
        if not self.config.force and os.path.exists(kmer_fasta) and os.path.exists(kmer_mapping_file):
            fasta_size = os.path.getsize(kmer_fasta)
            mapping_size = os.path.getsize(kmer_mapping_file)
            if fasta_size > 0 and mapping_size > 0:
                self.logger.info(f"k-mer FASTA和映射文件已存在，跳过生成|K-mer FASTA and mapping files exist, skipping: {kmer_fasta}")
                kmer_mapping = kmer_mapping_file
            else:
                kmer_mapping = self._create_kmer_fasta(merged_ps, kmer_fasta)
        else:
            kmer_mapping = self._create_kmer_fasta(merged_ps, kmer_fasta)

        # 选择比对工具并运行比对|Select alignment tool and run alignment
        alignment_result = os.path.join(output_dir, 'all_kmers_alignment_result.txt')

        if self.config.alignment_tool == 'bwa':
            self.logger.info(f"使用BWA进行k-mer比对|Using BWA for k-mer alignment")
            self.logger.info(f"BWA参数|BWA parameters: -k {self.config.bwa_k} -T {self.config.bwa_T} -a")
            self.logger.info(f"智能过滤阈值|Smart filtering threshold: AS_ratio = {self.config.as_ratio}")

            if not self._run_bwa_alignment(kmer_fasta, alignment_result):
                return False
        elif self.config.alignment_tool == 'blast':
            self.logger.info(f"使用BLASTN进行k-mer比对|Using BLASTN for k-mer alignment")

            # 创建BLAST数据库|Create BLAST database
            blast_db = self._create_blast_db()
            if not blast_db:
                return False

            # 运行BLAST比对|Run BLAST alignment
            if not self._run_blast(kmer_fasta, blast_db, alignment_result):
                return False
        else:
            self.logger.error(f"未知的比对工具|Unknown alignment tool: {self.config.alignment_tool}")
            return False

        # 生成GWAS格式文件：显著位点 + 所有位点（抽样）+ 所有位点（不抽样）|Generate GWAS format files: significant + all (sampled) + all (unsampled)
        self._generate_gwas_files(alignment_result, kmer_mapping, output_dir, threshold_strict, self.config.sample_ratio)

        self.logger.info("Post-GWAS分析完成|Post-GWAS analysis completed")

        # 如果提供了GFF文件，进行基因注释|If GFF file provided, perform gene annotation
        if self.config.gff_file:
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("开始基因注释|Starting gene annotation")
            self.logger.info("=" * 60)

            if not os.path.exists(self.config.gff_file):
                self.logger.error(f"GFF文件不存在|GFF file does not exist: {self.config.gff_file}")
                return False

            # 生成基因注释文件（针对严格阈值 p < 1e-5）|Generate gene annotation (for strict threshold)
            gene_annotation_file = os.path.join(output_dir, f'kmer_gene_annotation_pval_{threshold_strict}_window_{self.config.window_size/1000:.0f}kb.txt')

            # 检查是否需要跳过基因注释|Check if should skip gene annotation
            if not self.config.force and os.path.exists(gene_annotation_file):
                file_size = os.path.getsize(gene_annotation_file)
                if file_size > 0:
                    self.logger.info(f"基因注释文件已存在，跳过|Gene annotation file exists, skipping: {gene_annotation_file}")
                else:
                    self._generate_gene_annotation(
                        self.config.gff_file,
                        alignment_result,
                        kmer_mapping,
                        gene_annotation_file,
                        self.config.window_size,
                        threshold_strict
                    )
            else:
                self._generate_gene_annotation(
                    self.config.gff_file,
                    alignment_result,
                    kmer_mapping,
                    gene_annotation_file,
                    self.config.window_size,
                    threshold_strict
                )

            self.logger.info("")
            self.logger.info(f"基因注释文件|Gene annotation file: {gene_annotation_file}")

        return True

    def _merge_ps_files(self, ps_files, output_file):
        """合并所有.ps文件|Merge all .ps files"""
        self.logger.info(f"合并|Merging {len(ps_files)} 个关联结果文件|association result files")

        with open(output_file, 'w') as out_f:
            for ps_file in ps_files:
                with open(ps_file, 'r') as in_f:
                    for line in in_f:
                        if line.strip():
                            out_f.write(line)

        self.logger.info(f"合并完成|Merge completed")

    def _count_kmers(self, ps_file):
        """统计k-mer数量|Count k-mers"""
        count = 0
        with open(ps_file, 'r') as f:
            for line in f:
                if line.strip():
                    count += 1
        return count

    def _create_kmer_fasta(self, ps_file, fasta_file):
        """从.ps文件创建k-mer FASTA文件（包含所有k-mer）|Create k-mer FASTA file from .ps file (all k-mers)"""
        self.logger.info("创建所有k-mer的FASTA文件|Creating FASTA file for all k-mers")

        # 同时创建映射文件|Also create mapping file
        mapping_file = fasta_file.replace('.fasta', '.mapping.txt')

        kmer_count = 0
        with open(fasta_file, 'w') as out_f, open(mapping_file, 'w') as map_f:
            # 映射文件表头：KMER_ID, KMER_SEQ, EFFECT, PVAL
            map_f.write("KMER_ID\tKMER_SEQ\tEFFECT\tPVAL\n")

            with open(ps_file, 'r') as in_f:
                for line in in_f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        kmer = parts[0]
                        effect = float(parts[1])
                        pval = float(parts[2])

                        kmer_count += 1
                        kmer_id = f"kmer_{kmer_count}"
                        # FASTA格式：>kmer_ID\nkmer序列
                        out_f.write(f">{kmer_id}\n{kmer}\n")
                        # 映射文件：kmer_ID kmer_seq effect pval
                        map_f.write(f"{kmer_id}\t{kmer}\t{effect}\t{pval}\n")

        self.logger.info(f"创建了|Created {format_number(kmer_count)} 个k-mer记录|k-mer records")
        self.logger.info(f"k-mer映射文件|K-mer mapping file: {mapping_file}")

        return mapping_file

    def _create_blast_db(self):
        """创建BLAST数据库（如果不存在）|Create BLAST database if not exists"""
        genome_dir = os.path.dirname(self.config.genome_file)
        genome_base = os.path.basename(self.config.genome_file)
        db_name = os.path.join(genome_dir, os.path.splitext(genome_base)[0])

        # 检查数据库是否已存在|Check if database already exists
        if os.path.exists(db_name + '.nin') or os.path.exists(db_name + '.nhr'):
            self.logger.info(f"BLAST数据库已存在|BLAST database already exists: {db_name}")
            return db_name

        self.logger.info(f"创建BLAST数据库|Creating BLAST database for {self.config.genome_file}")

        cmd = [
            'makeblastdb',
            '-in', self.config.genome_file,
            '-dbtype', 'nucl',
            '-out', db_name
        ]

        success, _ = self.cmd_runner.run_command(
            cmd,
            description="创建BLAST数据库|Creating BLAST database"
        )

        if not success:
            self.logger.error("BLAST数据库创建失败|Failed to create BLAST database")
            return None

        return db_name

    def _run_blast(self, query_file, db_name, output_file):
        """运行BLAST比对|Run BLAST alignment"""
        self.logger.info("运行BLAST比对|Running BLAST alignment")
        self.logger.info("这可能需要较长时间...|This may take a while...")

        cmd = [
            'blastn',
            '-query', query_file,
            '-db', db_name,
            '-out', output_file,
            '-evalue', '1e-3',
            '-word_size', '9',
            '-outfmt', '6',  # Tabular format
            '-num_threads', str(self.config.threads)
        ]

        success, _ = self.cmd_runner.run_command(
            cmd,
            description="运行BLAST比对|Running BLAST alignment"
        )

        if not success:
            self.logger.error("BLAST比对失败|BLAST alignment failed")
            return False

        # 统计比对结果|Count alignment results
        hit_count = 0
        with open(output_file, 'r') as f:
            for line in f:
                if line.strip():
                    hit_count += 1

        self.logger.info(f"BLAST比对完成|BLAST alignment completed: {format_number(hit_count)} hits")

        return True

    def _check_bwa_index(self):
        """检查BWA索引是否存在|Check if BWA index exists"""
        # BWA索引文件后缀|BWA index file suffixes
        # 注意：BWA索引是直接在原文件名后加后缀，不是去掉扩展名
        # Note: BWA adds suffix to original filename, not removing extension
        index_suffixes = ['.bwt', '.pac', '.ann', '.amb', '.sa']

        # 检查所有索引文件是否存在|Check if all index files exist
        all_exist = True
        for suffix in index_suffixes:
            if not os.path.exists(self.config.genome_file + suffix):
                all_exist = False
                break

        if all_exist:
            self.logger.info(f"BWA索引已存在，跳过创建|BWA index already exists, skipping: {self.config.genome_file}")
            return self.config.genome_file

        return None

    def _create_bwa_index(self):
        """创建BWA索引|Create BWA index"""
        self.logger.info(f"创建BWA索引|Creating BWA index for {self.config.genome_file}")
        self.logger.info("这可能需要几分钟...|This may take a few minutes...")

        # 使用os.system直接执行，无缓冲|Use os.system for direct execution, no buffering
        cmd = f"bwa index {self.config.genome_file}"

        ret = os.system(cmd)
        if ret == 0:
            self.logger.info(f"BWA索引创建完成|BWA index created: {self.config.genome_file}")
            return self.config.genome_file
        else:
            self.logger.error(f"BWA索引创建失败|Failed to create BWA index (返回码|return code: {ret})")
            return None

    def _split_fasta_batches(self, fasta_file, output_dir, batch_size=1000000):
        """拆分大型FASTA文件为多个批次|Split large FASTA file into multiple batches

        Args:
            fasta_file: 输入FASTA文件|Input FASTA file
            output_dir: 输出目录|Output directory
            batch_size: 每个批次的序列数|Number of sequences per batch (default: 1M)

        Returns:
            批次文件列表|List of batch files
        """
        self.logger.info(f"拆分FASTA文件为批次（每批次|Splitting FASTA into batches ({batch_size} sequences/batch）")

        batch_files = []
        current_batch = 0
        seq_count = 0
        current_batch_file = None
        out_handle = None

        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # 新序列开始|New sequence starts
                        if seq_count % batch_size == 0:
                            # 需要创建新批次|Need to create new batch
                            if out_handle:
                                out_handle.close()

                            current_batch += 1
                            current_batch_file = os.path.join(output_dir, f'batch_{current_batch:03d}.fa')
                            batch_files.append(current_batch_file)
                            out_handle = open(current_batch_file, 'w')

                            if current_batch == 1 or (current_batch % 10) == 0:
                                self.logger.info(f"创建批次|Creating batch {current_batch}: {current_batch_file}")

                        seq_count += 1

                    # 写入行|Write line
                    if out_handle:
                        out_handle.write(line)

            if out_handle:
                out_handle.close()

            self.logger.info(f"FASTA拆分完成|FASTA splitting completed: {len(batch_files)} 个批次|batches, 总计|total {seq_count} 条序列|sequences")
            return batch_files

        except Exception as e:
            self.logger.error(f"FASTA拆分失败|FASTA splitting failed: {e}")
            if out_handle:
                out_handle.close()
            return []

    def _run_bwa_alignment(self, query_file, output_file):
        """运行BWA比对并转换为BLAST格式（分批处理，直接处理SAM）|Run BWA alignment and convert to BLAST format (batch processing, process SAM directly)

        优化版本|Optimized version:
        - 不再转换SAM到BAM|No longer convert SAM to BAM
        - 直接流式处理SAM文件|Stream process SAM files directly
        - 按批次过滤并输出文本结果|Filter and output text results by batch
        - 最后合并文本结果|Finally merge text results

        Args:
            query_file: 查询FASTA文件|Query FASTA file
            output_file: 输出文件（BLAST格式）|Output file (BLAST format)
        """
        self.logger.info("运行BWA比对|Running BWA alignment")
        self.logger.info("这可能需要较长时间...|This may take a while...")

        output_dir = os.path.dirname(output_file)

        # 检查/创建BWA索引|Check/Create BWA index
        genome_base = self._check_bwa_index()
        if not genome_base:
            genome_base = self._create_bwa_index()
            if not genome_base:
                return False

        # 统计序列数量|Count sequences
        self.logger.info("统计查询序列数量|Counting query sequences...")
        seq_count = 0
        with open(query_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seq_count += 1

        self.logger.info(f"查询序列数量|Query sequence count: {format_number(seq_count)}")

        # 创建批次目录（在拆分FASTA之前创建）|Create batch directory (before splitting FASTA)
        batch_dir = os.path.join(output_dir, 'batches')
        os.makedirs(batch_dir, exist_ok=True)

        # 判断是否需要分批处理|Check if batch processing is needed
        BATCH_THRESHOLD = 2000000
        batch_files = []

        if seq_count > BATCH_THRESHOLD:
            self.logger.warning(f"序列数量过大|Large sequence count ({format_number(seq_count)}), 启用分批处理|enabling batch processing")
            batch_files = self._split_fasta_batches(query_file, batch_dir, batch_size=BATCH_THRESHOLD)

            if not batch_files:
                self.logger.error("FASTA拆分失败|FASTA splitting failed")
                return False
        else:
            # 不需要分批，直接使用原文件|No batching needed, use original file directly
            batch_files = [query_file]

        # 处理每个批次|Process each batch

        # 存储每个批次的过滤结果文件|Store filtered result file for each batch
        batch_filtered_files = []

        for batch_idx, batch_file in enumerate(batch_files, 1):
            self.logger.info(f"=" * 60)
            self.logger.info(f"处理批次|Processing batch {batch_idx}/{len(batch_files)}")
            self.logger.info(f"批次文件|Batch file: {batch_file}")
            self.logger.info(f"=" * 60)

            batch_sam = os.path.join(batch_dir, f'batch_{batch_idx:03d}.sam')
            batch_filtered = os.path.join(batch_dir, f'batch_{batch_idx:03d}.filtered.txt')

            # 检查是否已存在过滤结果（断点续传）|Check if filtered result already exists (resume)
            if not self.config.force and os.path.exists(batch_filtered):
                file_size = os.path.getsize(batch_filtered)
                if file_size > 0:
                    self.logger.info(f"批次{batch_idx}过滤结果已存在，跳过|Batch {batch_idx} filtered result exists, skipping: {batch_filtered}")
                    batch_filtered_files.append(batch_filtered)
                    continue

            # 使用os.system执行BWA mem，使用-o参数直接指定输出文件
            # Use os.system to run BWA mem, use -o option to specify output file directly
            bwa_cmd = (f"bwa mem -k {self.config.bwa_k} -T {self.config.bwa_T} -a "
                       f"-t {self.config.threads} -o {batch_sam} {self.config.genome_file} {batch_file}")

            self.logger.info(f"运行BWA mem|Running BWA mem")
            self.logger.debug(f"命令|Command: {bwa_cmd}")

            ret = os.system(bwa_cmd)
            if ret != 0:
                self.logger.error(f"BWA mem比对失败|BWA mem alignment failed (返回码|return code: {ret})")
                return False

            # 检查SAM文件|Check SAM file
            if not os.path.exists(batch_sam) or os.path.getsize(batch_sam) == 0:
                self.logger.warning(f"批次{batch_idx}输出为空，跳过|Batch {batch_idx} output is empty, skipping")
                if os.path.exists(batch_sam):
                    os.remove(batch_sam)
                continue

            self.logger.info(f"BWA mem比对完成|BWA mem alignment completed (批次|batch {batch_idx})")

            # 直接流式处理SAM文件，应用智能过滤并输出|Stream process SAM file directly, apply smart filtering and output
            self._filter_sam_to_blast_format(batch_sam, batch_filtered, self.config.as_ratio)

            # 删除SAM文件以节省空间|Delete SAM file to save space
            if os.path.exists(batch_sam):
                os.remove(batch_sam)

            batch_filtered_files.append(batch_filtered)
            self.logger.info(f"批次{batch_idx}过滤完成|Batch {batch_idx} filtering completed")

        # 检查是否有成功的批次|Check if any batch succeeded
        if not batch_filtered_files:
            self.logger.error("没有成功比对的批次|No successful batches")
            return False

        # 合并所有批次的过滤结果|Merge all batch filtered results
        self.logger.info(f"=" * 60)
        self.logger.info(f"合并批次结果|Merging batch results: {len(batch_filtered_files)} 个批次文件|batch files")
        self.logger.info(f"=" * 60)

        with open(output_file, 'w') as out_f:
            for filtered_file in batch_filtered_files:
                if os.path.exists(filtered_file):
                    with open(filtered_file, 'r') as in_f:
                        out_f.write(in_f.read())

        self.logger.info(f"合并完成|Merge completed: {output_file}")

        # 清理批次文件|Clean batch files
        for filtered_file in batch_filtered_files:
            if os.path.exists(filtered_file):
                os.remove(filtered_file)

        # 清理批次目录|Clean batch directory
        if os.path.exists(batch_dir):
            try:
                import shutil
                shutil.rmtree(batch_dir)
                self.logger.info(f"清理批次目录|Cleaned batch directory: {batch_dir}")
            except Exception as e:
                self.logger.warning(f"清理批次目录失败|Failed to clean batch directory: {e}")

        # 如果是分批处理，删除临时拆分的FASTA文件|If batch processing, clean temporary split FASTA files
        if len(batch_files) > 1:
            for batch_file in batch_files:
                if os.path.exists(batch_file):
                    os.remove(batch_file)

        return True

    def _filter_sam_to_blast_format(self, sam_file, output_file, as_ratio=0.95):
        """流式处理SAM文件，应用智能过滤并转换为BLAST格式|Stream process SAM file, apply smart filtering and convert to BLAST format

        流式处理优化|Stream processing optimization:
        - 逐行读取SAM文件，避免内存溢出|Read SAM file line by line to avoid OOM
        - 按k-mer分组，应用智能AS过滤|Group by k-mer, apply smart AS filtering
        - 直接输出BLAST格式文本|Output BLAST format text directly

        Args:
            sam_file: 输入SAM文件|Input SAM file
            output_file: 输出文件（BLAST格式）|Output file (BLAST format)
            as_ratio: AS过滤阈值|AS filtering threshold (default: 0.95)
        """
        self.logger.info(f"应用智能过滤并转换输出格式（流式处理）|Applying smart filtering and converting output format (stream processing)")
        self.logger.info(f"AS阈值|AS threshold: AS_ratio = {as_ratio}")

        from collections import defaultdict

        # 第一遍扫描：收集每个k-mer的所有hits（流式处理）|First pass: collect all hits for each k-mer (stream processing)
        # 使用defaultdict避免预先分配内存|Use defaultdict to avoid pre-allocating memory
        kmer_hits = defaultdict(list)  # {kmer_id: [(chr, pos, end, mapq, as_score, cigar, pident, nm, gapopen, evalue, bitscore), ...]}

        line_count = 0
        with open(sam_file, 'r') as f:
            for line in f:
                line_count += 1

                # 每1000万行报告一次进度|Report progress every 10M lines
                if line_count % 10000000 == 0:
                    self.logger.info(f"已处理|Processed: {format_number(line_count)} 行|lines")

                if not line or line.startswith('@'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                qname = parts[0]  # k-mer ID
                flag = int(parts[1])
                chr_name = parts[2]
                pos = int(parts[3])  # 1-based
                mapq = int(parts[4])
                cigar = parts[5]
                seq = parts[9]

                # 跳过未比对的reads|Skip unmapped reads
                if flag & 0x4:  # BAM_FUNMAP
                    continue

                # 解析AS标签|Parse AS tag
                as_score = 0
                for i in range(11, len(parts)):
                    if parts[i].startswith('AS:i:'):
                        as_score = int(parts[i].split(':')[2])
                        break

                # 计算比对长度|Calculate alignment length
                aln_length = 0
                i = 0
                while i < len(cigar):
                    num_str = ''
                    while i < len(cigar) and cigar[i].isdigit():
                        num_str += cigar[i]
                        i += 1
                    if i < len(cigar):
                        op = cigar[i]
                        i += 1
                        if num_str:
                            num = int(num_str)
                            if op in 'MND=':  # Match/mismatch/deletion
                                aln_length += num

                # 计算end位置|Calculate end position
                end = pos + aln_length - 1

                # 计算pidentity（简化）|Calculate pidentity (simplified)
                pident = 100.0

                # 从NM标签获取mismatch数|Get mismatch count from NM tag
                nm = 0
                for i in range(11, len(parts)):
                    if parts[i].startswith('NM:i:'):
                        nm = int(parts[i].split(':')[2])
                        break

                # gapopen: 简化处理|gapopen: simplified
                gapopen = 0

                # 计算evalue（使用AS score的近似）|Calculate evalue (approximation using AS score)
                evalue = 10 ** (-as_score / 10.0)

                # bitscore: 使用AS score|bitscore: use AS score
                bitscore = as_score

                # 收集hits|Collect hits
                kmer_hits[qname].append((chr_name, pos, end, mapq, as_score, cigar, pident, nm, gapopen, evalue, bitscore))

        self.logger.info(f"SAM文件读取完成|SAM file reading completed: {format_number(line_count)} 行|lines, {len(kmer_hits)} k-mers")

        # 应用智能过滤并写入输出文件|Apply smart filtering and write to output file
        total_hits = 0
        filtered_hits = 0

        with open(output_file, 'w') as out_f:
            for kmer_id, hits in kmer_hits.items():
                if not hits:
                    continue

                # 找到最高AS score|Find maximum AS score
                max_as = max(hit[4] for hit in hits)

                # 保留AS >= max_AS * as_ratio的所有hits|Keep all hits with AS >= max_AS * as_ratio
                threshold = max_as * as_ratio
                good_hits = [hit for hit in hits if hit[4] >= threshold]

                for hit in good_hits:
                    chr_name, start, end, mapq, as_score, cigar, pident, nm, gapopen, evalue, bitscore = hit

                    # 跳过contig/scaffold|Skip contig/scaffold
                    if any(x in chr_name.lower() for x in ['chr00', 'tig', 'scaf', 'un']):
                        continue

                    # BLAST outfmt 6格式|BLAST outfmt 6 format
                    aln_len = end - start + 1
                    out_f.write(f"{kmer_id}\t{chr_name}\t{pident:.1f}\t{aln_len}\t{nm}\t{gapopen}\t")
                    out_f.write(f"1\t31\t{start}\t{end}\t{evalue:.2e}\t{bitscore:.1f}\n")

                    total_hits += 1
                    if hit in good_hits:
                        filtered_hits += 1

        self.logger.info(f"BWA比对及过滤完成|BWA alignment and filtering completed: {format_number(total_hits)} filtered hits (保留率|retention: {filtered_hits/total_hits*100:.1f}%)")

    def _filter_and_convert_bwa_to_blast_format(self, bam_file, output_file, as_ratio=0.95):
        """过滤BWA结果并转换为BLAST格式（智能全hit模式）|Filter BWA results and convert to BLAST format (smart all-hits mode)

        智能过滤策略|Smart filtering strategy:
        - 保留AS (alignment score) >= 最高AS * as_ratio的所有比对
        - Keep alignments with AS >= max_AS * as_ratio
        - 输出格式兼容BLAST outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        - Output format compatible with BLAST outfmt 6

        Args:
            bam_file: BWA BAM文件|BWA BAM file
            output_file: 输出文件（BLAST格式）|Output file (BLAST format)
            as_ratio: AS过滤阈值|AS filtering threshold (default: 0.95)
        """
        self.logger.info("应用智能过滤并转换输出格式|Applying smart filtering and converting output format")
        self.logger.info(f"AS阈值|AS threshold: AS_ratio = {as_ratio}")

        import subprocess

        # 使用samtools读取BAM并应用智能过滤|Use samtools to read BAM and apply smart filtering
        samtools_cmd = ['samtools', 'view', bam_file]
        wrapped_cmd = build_conda_command('samtools', ['view', bam_file])

        try:
            # 运行samtools view|Run samtools view
            result = subprocess.run(wrapped_cmd, capture_output=True, text=True, check=True)
            sam_lines = result.stdout.split('\n')
        except subprocess.CalledProcessError as e:
            self.logger.error(f"samtools view失败|samtools view failed: {e}")
            return False

        # 收集每个k-mer的所有hits|Collect all hits for each k-mer
        kmer_hits = {}  # {kmer_id: [(chr, pos, mapq, as_score, cigar, seq), ...]}

        for line in sam_lines:
            if not line or line.startswith('@'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            qname = parts[0]  # k-mer ID
            flag = int(parts[1])
            chr_name = parts[2]
            pos = int(parts[3])  # 1-based
            mapq = int(parts[4])
            cigar = parts[5]
            seq = parts[9]

            # 跳过未比对的reads|Skip unmapped reads
            if flag & 0x4:  # BAM_FUNMAP
                continue

            # 解析AS标签|Parse AS tag
            as_score = 0
            for i in range(11, len(parts)):
                if parts[i].startswith('AS:i:'):
                    as_score = int(parts[i].split(':')[2])
                    break

            # 计算比对长度和位置|Calculate alignment length and position
            # 从CIGAR计算比对长度|Calculate alignment length from CIGAR
            aln_length = 0
            cig_parts = []
            i = 0
            while i < len(cigar):
                num_str = ''
                while i < len(cigar) and cigar[i].isdigit():
                    num_str += cigar[i]
                    i += 1
                if i < len(cigar):
                    op = cigar[i]
                    i += 1
                    if num_str:
                        num = int(num_str)
                        cig_parts.append((num, op))
                        if op in 'MND=':  # Match/mismatch/deletion
                            aln_length += num

            # 计算end位置|Calculate end position
            end = pos + aln_length - 1

            # 计算pidentity|Calculate pidentity
            # 简化：假设没有mismatch，实际应该从CIGAR和NM标签计算
            # Simplified: assume no mismatches, should calculate from CIGAR and NM tag in reality
            pident = 100.0

            # 计算mismatch|Calculate mismatch
            # 从NM标签获取mismatch数|Get mismatch count from NM tag
            nm = 0
            for i in range(11, len(parts)):
                if parts[i].startswith('NM:i:'):
                    nm = int(parts[i].split(':')[2])
                    break

            # gapopen: 简化处理，假设为0|gapopen: simplified, assume 0
            gapopen = 0

            # 计算evalue（使用AS score的近似）|Calculate evalue (approximation using AS score)
            # BWA没有evalue，使用AS score的负指数作为近似|BWA doesn't have evalue, use negative exponential of AS score as approximation
            evalue = 10 ** (-as_score / 10.0)

            # bitscore: 使用AS score|bitscore: use AS score
            bitscore = as_score

            # 收集hits|Collect hits
            if qname not in kmer_hits:
                kmer_hits[qname] = []
            kmer_hits[qname].append((chr_name, pos, end, mapq, as_score, cigar, pident, nm, gapopen, evalue, bitscore))

        # 应用智能过滤|Apply smart filtering
        filtered_hits = {}
        for kmer_id, hits in kmer_hits.items():
            if not hits:
                continue

            # 找到最高AS score|Find maximum AS score
            max_as = max(hit[4] for hit in hits)

            # 保留AS >= max_AS * as_ratio的所有hits|Keep all hits with AS >= max_AS * as_ratio
            threshold = max_as * as_ratio
            good_hits = [hit for hit in hits if hit[4] >= threshold]
            filtered_hits[kmer_id] = good_hits

        # 写入BLAST格式文件|Write BLAST format file
        total_hits = 0
        with open(output_file, 'w') as out_f:
            for kmer_id, hits in filtered_hits.items():
                for hit in hits:
                    chr_name, start, end, mapq, as_score, cigar, pident, nm, gapopen, evalue, bitscore = hit

                    # 跳过contig/scaffold|Skip contig/scaffold
                    if any(x in chr_name.lower() for x in ['chr00', 'tig', 'scaf', 'un']):
                        continue

                    # BLAST outfmt 6格式:
                    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                    # 对于k-mer，qstart=1, qend=seq_len|For k-mers, qstart=1, qend=seq_len
                    aln_len = end - start + 1
                    out_f.write(f"{kmer_id}\t{chr_name}\t{pident:.1f}\t{aln_len}\t{nm}\t{gapopen}\t")
                    out_f.write(f"1\t31\t{start}\t{end}\t{evalue:.2e}\t{bitscore:.1f}\n")

                    total_hits += 1

        self.logger.info(f"BWA比对及过滤完成|BWA alignment and filtering completed: {format_number(total_hits)} filtered hits")

    def _extract_chrom_number(self, chr_name):
        """从染色体名称中提取数字编号|Extract chromosome number from name

        Args:
            chr_name: 染色体名称|Chromosome name (e.g., 'Chr1', 'Chr11', 'chr01', 'chrX')

        Returns:
            int: 染色体编号（如果是字母则返回大数字）|Chromosome number (returns large number for non-numeric)
        """
        import re

        # 转小写处理|Convert to lowercase for processing
        chr_lower = chr_name.lower()

        # 提取数字|Extract number
        match = re.search(r'chr(\d+)', chr_lower)
        if match:
            return int(match.group(1))

        # 处理性染色体|Handle sex chromosomes
        if 'x' in chr_lower:
            return 1000
        elif 'y' in chr_lower:
            return 1001
        elif 'm' in chr_lower or 'mt' in chr_lower:
            return 1002  # 线粒体|Mitochondrial

        # 其他情况返回大数字|Return large number for other cases
        return 9999

    def _generate_gwas_files(self, blast_file, mapping_file, output_dir, p_threshold, sample_ratio=0.1):
        """生成GWAS格式文件：显著位点文件 + 所有位点文件（含抽样）|Generate GWAS format files: significant + all (with sampling)

        智能全hit模式策略|Smart all-hits mode strategy:
        - 保留bitscore相等的所有hits（不丢弃重复位点）|Keep all hits with equal bitscore (don't discard duplicate loci)
        - 过滤掉bitscore明显较低的hits（去除低质量比对）|Filter out hits with significantly lower bitscore (remove low-quality alignments)

        抽样策略|Sampling strategy:
        - p < threshold: 全部保留到显著位点文件|Keep all in significant file
        - p >= threshold: 随机抽样指定比例到所有位点文件|Random sampling ratio to all file

        排序策略|Sorting strategy:
        - 按染色体编号 + 物理位置排序|Sort by chromosome number + physical position

        Args:
            blast_file: 比对结果文件（BWA/BLAST）|Alignment result file (BWA/BLAST)
            mapping_file: k-mer映射文件|k-mer mapping file
            output_dir: 输出目录|Output directory
            p_threshold: p值阈值（显著阈值）|p-value threshold (significance threshold)
            sample_ratio: 高p值位点抽样比例|Sampling ratio for high p-value loci (default: 0.1 = 10%)
        """
        self.logger.info(f"生成GWAS格式文件|Generating GWAS format files")
        self.logger.info("采用智能全hit模式（保留bitscore相等的所有重复位点）|Using smart all-hits mode (keep all duplicate loci with equal bitscore)")
        self.logger.info(f"显著阈值|Significance threshold: p < {p_threshold}")
        self.logger.info(f"抽样策略：p >= {p_threshold} 每条染色体随机抽样{sample_ratio*100:.0f}%")
        self.logger.info(f"排序策略：按染色体编号 + 物理位置排序|Sorting: by chromosome number + physical position")

        # 设置随机种子以确保可重复性|Set random seed for reproducibility
        random.seed(42)

        # 读取k-mer映射|Read k-mer mapping
        kmer_map = {}
        with open(mapping_file, 'r') as f:
            next(f)  # 跳过表头|Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    kmer_id = parts[0]
                    kmer_seq = parts[1]
                    effect = float(parts[2])
                    pval = float(parts[3])
                    kmer_map[kmer_id] = (kmer_seq, effect, pval)

        # 过滤比对结果（智能全hit模式）|Filter alignment results (smart all-hits mode)
        # BLAST outfmt 6列: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        kmer_hits = {}  # 每个k-mer的所有hits|All hits for each k-mer

        with open(blast_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue

                kmer_id = parts[0]  # 格式：kmer_序号|Format: kmer_number
                chr_name = parts[1]
                start = int(parts[8])  # sstart
                end = int(parts[9])    # send
                evalue = float(parts[10])
                bitscore = float(parts[11])

                # 确保start < end|Ensure start < end
                if start > end:
                    start, end = end, start

                # 去除contig/scaffold|Remove contig/scaffold
                if any(x in chr_name.lower() for x in ['chr00', 'tig', 'scaf', 'un']):
                    continue

                # 收集每个k-mer的所有hits|Collect all hits for each k-mer
                if kmer_id not in kmer_hits:
                    kmer_hits[kmer_id] = []
                kmer_hits[kmer_id].append((chr_name, start, end, evalue, bitscore))

        # 对每个k-mer的hits应用智能过滤|Apply smart filtering for each k-mer's hits
        filtered_hits = {}
        for kmer_id, hits in kmer_hits.items():
            if not hits:
                continue

            # 找到最高bitscore|Find maximum bitscore
            max_bitscore = max(hit[4] for hit in hits)

            # 保留bitscore >= 最高bitscore * 0.95 的所有hits（容许5%的差异）
            # Keep all hits with bitscore >= 95% of maximum (allow 5% difference)
            threshold = max_bitscore * 0.95
            good_hits = [hit for hit in hits if hit[4] >= threshold]
            filtered_hits[kmer_id] = good_hits

        # 收集数据：显著位点列表 + 所有位点列表（含抽样）+ 所有位点（不抽样）|Collect data: significant + all (sampled) + all (unsampled)
        significant_records = []  # 显著位点|Significant loci
        all_records_sampled = []  # 所有位点（显著 + 抽样）|All loci (significant + sampled)
        all_records_unsampled = []  # 所有位点（不抽样，完整数据）|All loci (unsampled, complete)

        for kmer_id, hits in filtered_hits.items():
            if kmer_id not in kmer_map:
                continue

            kmer_seq, effect, pval = kmer_map[kmer_id]

            # 处理每个hit|Process each hit
            for chr_name, start, end, evalue, bitscore in hits:
                # 创建记录|Create record
                record = {
                    'kmer_id': kmer_id,
                    'kmer_seq': kmer_seq,
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'effect': effect,
                    'pval': pval,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'chr_num': self._extract_chrom_number(chr_name)
                }

                # 判断是否是显著位点|Check if significant
                if pval < p_threshold:
                    # 显著位点：加入所有三个文件|Significant: add to all three files
                    significant_records.append(record)
                    all_records_sampled.append(record)
                    all_records_unsampled.append(record)
                else:
                    # 非显著位点：无抽样文件全部加入，抽样文件按比例抽样|Non-significant: add all to unsampled, sample for sampled file
                    all_records_unsampled.append(record)

                    # 使用k-mer ID作为随机种子进行抽样|Sample using k-mer ID as seed
                    kmer_seed = hash(kmer_id) % (2**32)
                    random.seed(kmer_seed)
                    if random.random() < sample_ratio:
                        all_records_sampled.append(record)

        # 重置随机种子|Reset random seed
        random.seed(42)

        # 排序函数：按染色体编号 + 物理位置排序|Sort function: by chromosome number + physical position
        sort_key = lambda x: (x['chr_num'], x['start'], x['end'], x['kmer_id'])

        significant_records_sorted = sorted(significant_records, key=sort_key)
        all_records_sampled_sorted = sorted(all_records_sampled, key=sort_key)
        all_records_unsampled_sorted = sorted(all_records_unsampled, key=sort_key)

        # 写入显著位点文件|Write significant loci file
        significant_file = os.path.join(output_dir, 'kmer_gwas_significant.txt')
        with open(significant_file, 'w') as out_f:
            # 写入表头|Write header
            out_f.write("KMER_ID\tKMER_SEQ\tCHR\tSTART\tEND\tEFFECT\tP\tALIGN_E_VALUE\tALIGN_BITSCORE\n")

            for rec in significant_records_sorted:
                out_f.write(f"{rec['kmer_id']}\t{rec['kmer_seq']}\t{rec['chr']}\t{rec['start']}\t{rec['end']}\t"
                           f"{rec['effect']}\t{rec['pval']}\t{rec['evalue']}\t{rec['bitscore']}\n")

        # 写入所有位点文件（含抽样）|Write all loci file (with sampling)
        all_sampled_file = os.path.join(output_dir, 'kmer_gwas_all_sampled.txt')
        with open(all_sampled_file, 'w') as out_f:
            # 写入表头|Write header
            out_f.write("KMER_ID\tKMER_SEQ\tCHR\tSTART\tEND\tEFFECT\tP\tALIGN_E_VALUE\tALIGN_BITSCORE\n")

            for rec in all_records_sampled_sorted:
                out_f.write(f"{rec['kmer_id']}\t{rec['kmer_seq']}\t{rec['chr']}\t{rec['start']}\t{rec['end']}\t"
                           f"{rec['effect']}\t{rec['pval']}\t{rec['evalue']}\t{rec['bitscore']}\n")

        # 写入所有位点文件（不抽样，完整数据）|Write all loci file (unsampled, complete)
        all_unsampled_file = os.path.join(output_dir, 'kmer_gwas_all.txt')
        with open(all_unsampled_file, 'w') as out_f:
            # 写入表头|Write header
            out_f.write("KMER_ID\tKMER_SEQ\tCHR\tSTART\tEND\tEFFECT\tP\tALIGN_E_VALUE\tALIGN_BITSCORE\n")

            for rec in all_records_unsampled_sorted:
                out_f.write(f"{rec['kmer_id']}\t{rec['kmer_seq']}\t{rec['chr']}\t{rec['start']}\t{rec['end']}\t"
                           f"{rec['effect']}\t{rec['pval']}\t{rec['evalue']}\t{rec['bitscore']}\n")

        # 统计信息|Statistics
        sig_kmers = len(set(r['kmer_id'] for r in significant_records))
        sampled_kmers = len(set(r['kmer_id'] for r in all_records_sampled)) - sig_kmers
        unsampled_kmers = len(set(r['kmer_id'] for r in all_records_unsampled)) - sig_kmers

        self.logger.info(f"文件生成完成|Files generated:")
        self.logger.info(f"")

        self.logger.info(f"1. 显著位点文件|Significant loci file:")
        self.logger.info(f"   文件|File: {significant_file}")
        self.logger.info(f"   位点数|Loci count: {len(significant_records_sorted)}")
        self.logger.info(f"   k-mer数|k-mer count: {sig_kmers}")
        self.logger.info(f"")

        self.logger.info(f"2. 所有位点文件（含抽样）|All loci file (with sampling):")
        self.logger.info(f"   文件|File: {all_sampled_file}")
        self.logger.info(f"   位点数|Loci count: {len(all_records_sampled_sorted)}")
        self.logger.info(f"   k-mer数|k-mer count: {len(set(r['kmer_id'] for r in all_records_sampled))}")
        self.logger.info(f"     其中显著|Including significant: {sig_kmers}")
        self.logger.info(f"     抽样|Sampled from non-significant: {sampled_kmers}")
        self.logger.info(f"")

        self.logger.info(f"3. 所有位点文件（不抽样，完整数据）|All loci file (unsampled, complete):")
        self.logger.info(f"   文件|File: {all_unsampled_file}")
        self.logger.info(f"   位点数|Loci count: {len(all_records_unsampled_sorted)}")
        self.logger.info(f"   k-mer数|k-mer count: {len(set(r['kmer_id'] for r in all_records_unsampled))}")
        self.logger.info(f"     其中显著|Including significant: {sig_kmers}")
        self.logger.info(f"     非显著（全部）|Non-significant (all): {unsampled_kmers}")

    def _parse_gff(self, gff_file):
        """解析GFF文件，建立基因位置索引|Parse GFF file and build gene position index

        Returns:
            dict: {chr_name: [(start, end, gene_id, gene_name, attributes), ...]}
        """
        self.logger.info(f"解析GFF文件|Parsing GFF file: {gff_file}")

        genes_by_chr = {}

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                seqtype = parts[2]
                if seqtype.lower() not in ['gene', 'mrna', 'cds', 'exon']:
                    # 只提取基因相关特征，或者可以提取所有特征
                    # 这里我们提取 gene 类型的特征
                    if seqtype.lower() != 'gene':
                        continue

                chr_name = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # 解析属性字段|Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                # 提取基因ID和名称|Extract gene ID and name
                gene_id = attr_dict.get('ID', attr_dict.get('gene_id', ''))
                gene_name = attr_dict.get('Name', attr_dict.get('gene_name', attr_dict.get('gene', '')))
                product = attr_dict.get('product', attr_dict.get('description', ''))

                if chr_name not in genes_by_chr:
                    genes_by_chr[chr_name] = []

                # 确保start < end
                if start > end:
                    start, end = end, start

                genes_by_chr[chr_name].append({
                    'start': start,
                    'end': end,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'product': product,
                    'strand': strand
                })

        # 统计基因数量|Count genes
        total_genes = sum(len(genes) for genes in genes_by_chr.values())
        self.logger.info(f"解析完成|Parsed: {format_number(total_genes)} 个基因特征|gene features on {len(genes_by_chr)} 条染色体|chromosomes")

        return genes_by_chr

    def _find_nearby_genes(self, chr_name, start, end, genes_by_chr, window_size):
        """查找指定区间上下游范围内的基因|Find genes within upstream/downstream window

        Args:
            chr_name: 染色体名称|Chromosome name
            start: 区间起始|Interval start
            end: 区间终止|Interval end
            genes_by_chr: 基因索引|Gene index
            window_size: 窗口大小|Window size (bp)

        Returns:
            list: 附近的基因列表|List of nearby genes
        """
        if chr_name not in genes_by_chr:
            return []

        # 扩展的查询范围|Extended query range
        query_start = max(0, start - window_size)
        query_end = end + window_size

        nearby_genes = []

        for gene in genes_by_chr[chr_name]:
            gene_start = gene['start']
            gene_end = gene['end']

            # 检查是否在窗口内|Check if within window
            # 基因的终止位置 > 查询起始 AND 基因的起始位置 < 查询终止
            if gene_end >= query_start and gene_start <= query_end:
                # 计算距离|Calculate distance
                if gene_end < start:
                    # 基因在k-mer上游|Gene is upstream
                    distance = start - gene_end
                elif gene_start > end:
                    # 基因在k-mer下游|Gene is downstream
                    distance = gene_start - end
                else:
                    # 基因与k-mer重叠|Gene overlaps with k-mer
                    distance = 0

                nearby_genes.append({
                    **gene,
                    'distance': distance
                })

        # 按距离排序|Sort by distance
        nearby_genes.sort(key=lambda x: x['distance'])

        return nearby_genes

    def _generate_gene_annotation(self, gff_file, blast_result, mapping_file,
                                   output_file, window_size, p_threshold=1e-5):
        """生成基因注释文件（仅针对显著位点）|Generate gene annotation file (significant loci only)

        Args:
            gff_file: GFF注释文件|GFF annotation file
            blast_result: BLAST结果文件|BLAST result file
            mapping_file: k-mer映射文件|k-mer mapping file
            output_file: 输出文件|Output file
            window_size: 窗口大小|Window size (bp)
            p_threshold: p值阈值（仅注释低于此阈值的位点）|p-value threshold (only annotate loci below this)
        """
        self.logger.info(f"生成基因注释文件 (p < {p_threshold}, window={window_size/1000:.0f}kb)|Generating gene annotation file")
        self.logger.info(f"解析GFF文件|Parsing GFF file: {gff_file}")

        # 解析GFF文件|Parse GFF file
        genes_by_chr = self._parse_gff(gff_file)

        # 读取k-mer映射|Read k-mer mapping
        self.logger.info("读取k-mer映射|Reading k-mer mapping")
        kmer_map = {}
        with open(mapping_file, 'r') as f:
            next(f)  # 跳过表头|Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    kmer_id = parts[0]
                    kmer_seq = parts[1]
                    effect = float(parts[2])
                    pval = float(parts[3])
                    kmer_map[kmer_id] = (kmer_seq, effect, pval)

        # 读取BLAST结果并生成注释|Read BLAST results and generate annotations
        self.logger.info("查找显著位点附近的基因|Finding genes near significant loci")

        with open(output_file, 'w') as out_f:
            # 写入表头|Write header
            out_f.write("KMER_ID\tKMER_SEQ\tCHR\tSTART\tEND\tP\tGENE_ID\tGENE_NAME\tDISTANCE(bp)\tSTRAND\tGENE_START\tGENE_END\tPRODUCT\n")

            total_kmers = 0
            significant_kmers = 0
            annotated_kmers = 0

            with open(blast_result, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue

                    kmer_id = parts[0]
                    chr_name = parts[1]
                    start = int(parts[8])
                    end = int(parts[9])

                    if start > end:
                        start, end = end, start

                    # 去除contig/scaffold|Remove contig/scaffold
                    if any(x in chr_name.lower() for x in ['chr00', 'tig', 'scaf', 'un']):
                        continue

                    total_kmers += 1

                    if kmer_id not in kmer_map:
                        continue

                    kmer_seq, effect, pval = kmer_map[kmer_id]

                    # 只处理显著位点|Only process significant loci
                    if pval >= p_threshold:
                        continue

                    significant_kmers += 1

                    # 查找附近的基因|Find nearby genes
                    nearby_genes = self._find_nearby_genes(chr_name, start, end, genes_by_chr, window_size)

                    if nearby_genes:
                        annotated_kmers += 1

                    # 为每个找到的基因写一行|Write one line for each found gene
                    if nearby_genes:
                        for gene in nearby_genes:
                            out_f.write(f"{kmer_id}\t{kmer_seq}\t{chr_name}\t{start}\t{end}\t{pval}\t"
                                       f"{gene['gene_id']}\t{gene['gene_name']}\t{gene['distance']}\t"
                                       f"{gene['strand']}\t{gene['start']}\t{gene['end']}\t{gene['product']}\n")
                    else:
                        # 即使没有找到基因，也输出一行|Output a line even if no genes found
                        out_f.write(f"{kmer_id}\t{kmer_seq}\t{chr_name}\t{start}\t{end}\t{pval}\t"
                                   f"NA\tNA\tNA\tNA\tNA\tNA\n")

        self.logger.info(f"基因注释完成|Gene annotation completed")
        self.logger.info(f"  总k-mer位点|Total k-mer loci: {format_number(total_kmers)}")
        self.logger.info(f"  显著位点（p < {p_threshold}）|Significant loci (p < {p_threshold}): {format_number(significant_kmers)}")
        self.logger.info(f"  有注释的位点|Annotated loci: {format_number(annotated_kmers)}")
        self.logger.info(f"文件已保存|File saved: {output_file}")
