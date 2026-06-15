"""
转换器模块|Converter Module
"""

from collections import defaultdict
from .utils import print_interval_progress, format_number, open_input


class KmerToVcfConverter:
    """Kmer转VCF转换器|Kmer to VCF Converter"""

    def __init__(self, config, logger):
        """
        初始化转换器|Initialize converter

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

    def process_sorted_file_and_convert(self):
        """
        Pass 3: 扫描排序文件，统计聚合频次并转换为VCF|Process sorted file, calculate aggregated frequencies and convert to VCF

        Returns:
            int: 保留的行数|Number of kept lines
        """
        self.logger.info("=" * 100)
        self.logger.info(" Step 3: 统计聚合频次并转换为VCF|Step 3: Calculating aggregated frequencies and converting to VCF")
        self.logger.info("=" * 100)

        # 判断是否使用分桶模式|Check if use bucket mode
        use_bucket_mode = self.config.should_use_bucket_mode()

        if use_bucket_mode:
            self.logger.info("使用分桶模式处理|Using bucket processing mode")
            return self._process_buckets_and_convert()
        else:
            self.logger.info("使用标准模式处理|Using standard processing mode")
            return self._process_standard_mode()

    def _process_standard_mode(self):
        """
        标准模式：先统计频次，再写VCF|Standard mode: count frequencies first, then write VCF

        Returns:
            int: 保留的行数|Number of kept lines
        """
        sorted_file = self.config.get_temp_file_path('sorted_binary_matrix.txt')

        # 第一遍扫描：统计每个kmer的聚合频次|First scan: calculate aggregated frequencies
        kmer_agg_freq = self._calculate_aggregated_frequencies(sorted_file)

        # 第二遍扫描：从原始文件读取丰度并输出|Second scan: read abundances from original file and output
        kept_lines = self._write_vcf_and_txt(kmer_agg_freq)

        return kept_lines

    def _process_buckets_and_convert(self):
        """
        分桶模式：并行处理所有桶|Bucket mode: process all buckets in parallel

        Returns:
            int: 保留的行数|Number of kept lines
        """
        from .parallel_optimized import BucketProcessor

        bucket_dir = self.config.get_temp_file_path('buckets')
        output_vcf = self.config.output_vcf
        output_txt = self.config.get_output_txt_path()
        min_agg_count = self.config.min_agg_count

        try:
            processor = BucketProcessor(self.config, self.logger)

            # 并行处理所有桶：统计 + 过滤 + 排序 + 写VCF
            kept_lines = processor.process_buckets_in_parallel(
                bucket_dir=bucket_dir,
                output_vcf=output_vcf,
                output_txt=output_txt,
                min_agg_count=min_agg_count
            )

            return kept_lines

        except Exception as e:
            self.logger.error(f"分桶处理失败|Bucket processing failed: {e}")
            raise

    def _calculate_aggregated_frequencies(self, sorted_file):
        """
        计算每个kmer的聚合频次|Calculate aggregated frequencies for each kmer

        Args:
            sorted_file: 排序后的文件|Sorted file path

        Returns:
            dict: kmer -> 聚合频次|Kmer to aggregated frequency mapping
        """
        self.logger.info("统计聚合频次|Calculating aggregated frequencies...")

        kmer_agg_freq = {}
        current_kmer = None
        current_count = 0
        line_count = 0

        with open(sorted_file, 'r') as f:
            header = f.readline()  # 跳过表头|Skip header

            for line in f:
                line_count += 1
                if line_count % 1000000 == 0:
                    print_interval_progress(line_count, 1000000)

                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                kmer_id = parts[0]
                binary_str = ' '.join(parts[1:])  # 重新组合二进制字符串
                count = binary_str.count('1')  # 统计该行中'1'的个数

                # 如果是新kmer，保存上一个kmer的统计|If new kmer, save previous kmer's stats
                if kmer_id != current_kmer:
                    if current_kmer is not None:
                        kmer_agg_freq[current_kmer] = current_count
                    current_kmer = kmer_id
                    current_count = 0

                # 累加当前kmer的频次|Accumulate current kmer's frequency
                current_count += count

            # 保存最后一个kmer|Save last kmer
            if current_kmer is not None:
                kmer_agg_freq[current_kmer] = current_count

        print()  # 换行|New line

        # 统计信息|Statistics
        total_unique_kmers = len(kmer_agg_freq)
        min_freq = min(kmer_agg_freq.values())
        max_freq = max(kmer_agg_freq.values())

        self.logger.info(f"唯一kmer数|Unique kmers: {format_number(total_unique_kmers)}")
        self.logger.info(f"聚合频次范围|Aggregated frequency range: [{min_freq}, {max_freq}]")

        # 统计过滤情况|Filtering statistics
        min_cut = self.config.min_agg_count
        max_cut = self.config.num_samples

        kept = sum(1 for f in kmer_agg_freq.values() if min_cut <= f <= max_cut)
        filtered_low = sum(1 for f in kmer_agg_freq.values() if f < min_cut)
        filtered_high = sum(1 for f in kmer_agg_freq.values() if f > max_cut)

        self.logger.info(f"保留kmer数|Kmers kept: {format_number(kept)}")
        self.logger.info(f"过滤(过低)|Filtered (too rare): {format_number(filtered_low)}")
        self.logger.info(f"过滤(过高)|Filtered (multi-copy): {format_number(filtered_high)}")

        return kmer_agg_freq

    def _write_vcf_and_txt(self, kmer_agg_freq):
        """
        从原始文件读取丰度，过滤并写入VCF和TXT|Read abundances from original file, filter and write VCF and TXT

        Args:
            kmer_agg_freq: kmer聚合频次字典|Kmer aggregated frequency dictionary

        Returns:
            int: 保留的行数|Number of kept lines
        """
        self.logger.info("写入VCF和TXT文件|Writing VCF and TXT files...")

        # 决定输出路径|Determine output paths
        if self.config.is_compressed_output():
            # 如果输出压缩，先写临时文件，最后压缩|If compressed output, write temp file first
            vcf_path = self.config.get_temp_file_path('temp_output.vcf')
        else:
            vcf_path = self.config.output_vcf

        txt_path = self.config.get_output_txt_path()

        kept_lines = 0
        min_cut = self.config.min_agg_count
        max_cut = self.config.num_samples

        # 同时读取原始文件和排序文件|Read original file and sorted file simultaneously
        with open_input(self.config.input_matrix) as f_orig, \
             open(vcf_path, 'w') as f_vcf, \
             open(txt_path, 'w') as f_txt:

            # 写入表头|Write headers
            self._write_vcf_header(f_vcf)
            txt_header = ["KmerID"] + self.config.samples
            f_txt.write("\t".join(txt_header) + "\n")

            line_idx = 0
            header = f_orig.readline()  # 跳过表头|Skip header

            for line in f_orig:
                line_idx += 1
                if line_idx % 100000 == 0:
                    print_interval_progress(line_idx, 100000)

                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                kmer_id = parts[0]
                abundances = parts[1:]

                # 检查是否保留|Check if should keep
                agg_freq = kmer_agg_freq.get(kmer_id, 0)

                if min_cut <= agg_freq <= max_cut:
                    kept_lines += 1

                    # 写入TXT|Write TXT
                    f_txt.write(line)

                    # 写入VCF|Write VCF
                    vcf_row = self._convert_to_vcf_row(kmer_id, abundances, line_idx)
                    f_vcf.write(vcf_row)

        print()  # 换行|New line

        # 如果需要压缩|If compression needed
        if self.config.is_compressed_output():
            from .utils import compress_vcf_with_bgzip
            compress_vcf_with_bgzip(vcf_path, self.config.output_vcf, self.logger)

        self.logger.info(f"保留行数|Kept lines: {format_number(kept_lines)}")
        self.logger.info(f"VCF文件|VCF file: {self.config.output_vcf}")
        self.logger.info(f"TXT文件|TXT file: {txt_path}")

        return kept_lines

    def _write_vcf_header(self, f_vcf):
        """
        写入VCF文件头|Write VCF file header

        Args:
            f_vcf: VCF文件对象|VCF file object
        """
        min_cut = self.config.min_agg_count
        max_cut = self.config.num_samples

        # VCF元信息|VCF meta information
        f_vcf.write("##fileformat=VCFv4.2\n")
        f_vcf.write("##source=kmer2vcf\n")
        f_vcf.write(f"##filter_agg_range={min_cut}_to_{max_cut}\n")
        f_vcf.write("##INFO=<ID=KMER,Number=1,Type=String,Description=\"Kmer sequence\">\n")
        f_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f_vcf.write(f"##contig=<ID=Chr_kmer,length=2147483647>\n")

        # VCF列头|VCF column header
        vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + self.config.samples
        f_vcf.write("\t".join(vcf_header) + "\n")

    def _convert_to_vcf_row(self, kmer_id, abundances, pos):
        """
        转换单行为VCF格式|Convert single line to VCF format

        Args:
            kmer_id: kmer序列|Kmer sequence
            abundances: 丰度列表|Abundance list
            pos: 位置|Position

        Returns:
            str: VCF格式行|VCF format line
        """
        # 转换基因型|Convert genotypes
        # 丰度>0 -> 0/0 (有kmer，参考型)|Abundance>0 -> 0/0 (has kmer, reference)
        # 丰度=0 -> 1/1 (无kmer，缺失)|Abundance=0 -> 1/1 (no kmer, deletion)
        genotypes = []
        for abd in abundances:
            if int(abd) > 0:
                genotypes.append("0/0")
            else:
                genotypes.append("1/1")

        # 构建VCF行|Build VCF row
        # CHROM=Chr_kmer, POS=pos, ID=kmer_id
        # REF=kmer_id, ALT=., QUAL=., FILTER=PASS, INFO=KMER=kmer_id
        vcf_row = [
            "Chr_kmer",          # CHROM
            str(pos),            # POS
            ".",                 # ID
            kmer_id,             # REF
            ".",                 # ALT (deletion)
            ".",                 # QUAL
            "PASS",              # FILTER
            f"KMER={kmer_id}",   # INFO
            "GT"                 # FORMAT
        ] + genotypes

        return "\t".join(vcf_row) + "\n"
