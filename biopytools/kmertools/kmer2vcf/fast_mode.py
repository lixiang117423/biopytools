"""
快速模式转换器|Fast Mode Converter
基于kmer_matrix_to_vcf.py，单次遍历快速转换|Single-pass fast conversion based on kmer_matrix_to_vcf.py
"""

import io
import sys
from pathlib import Path
from typing import Optional, TextIO
from .utils import open_input


class FastModeConverter:
    """快速模式转换器|Fast Mode Converter - Single pass conversion"""

    def __init__(self, config, logger):
        """
        初始化快速模式转换器|Initialize fast mode converter

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

        # 配置参数|Configuration parameters
        self.input_matrix = str(config.input_matrix)
        self.output_vcf = str(config.output_vcf)
        self.chr_length = getattr(config, 'chr_length', 100000000)  # 默认100M|Default 100M
        self.chr_number = getattr(config, 'chr_number', 0)  # 染色体数量（0表示自动）|Chromosome number (0=auto)
        self.min_freq = getattr(config, 'min_freq', 0)  # 最小出现频次|Minimum frequency
        self.kmer_length = getattr(config, 'kmer_length', 51)  # kmer长度|Kmer length
        self.no_header = getattr(config, 'no_header', False)  # 无header行|No header line

    def convert(self) -> int:
        """
        执行快速转换|Execute fast conversion

        Returns:
            int: 输出的kmer数量|Number of kmers output
        """
        self.logger.info("=" * 100)
        self.logger.info(" 快速模式转换（单次遍历）|Fast Mode Conversion (Single Pass)")
        self.logger.info("=" * 100)
        self.logger.info(f"输入文件|Input file: {self.input_matrix}")
        self.logger.info(f"输出VCF|Output VCF: {self.output_vcf}")
        if self.chr_number > 0:
            self.logger.info(f"用户指定染色体数量|User-specified chromosome count: {self.chr_number:,}")
        else:
            self.logger.info(f"每条染色体长度|Chromosome length: {self.chr_length:,}")
        if self.min_freq > 0:
            self.logger.info(f"最小频次过滤|Minimum frequency filter: {self.min_freq}")
        self.logger.info("")

        # 决定输出路径|Determine output paths
        if self.output_vcf.endswith('.vcf.gz'):
            # 如果输出压缩，先写临时文件，最后压缩|If compressed output, write temp file first
            import tempfile
            temp_vcf = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False, dir=self.config.temp_dir)
            actual_output_path = temp_vcf.name
            need_compress = True
        else:
            actual_output_path = self.output_vcf
            need_compress = False

        try:
            # 打开输入文件|Open input file
            input_file = self._open_input_file()

            # 打开输出文件|Open output file
            if need_compress:
                output_file = temp_vcf
            else:
                output_file = open(actual_output_path, 'w')

            # 执行转换|Perform conversion
            output_count = self._convert_fast(input_file, output_file)

            # 关闭文件|Close files
            input_file.close()
            output_file.close()

            # 如果需要压缩|If compression needed
            if need_compress:
                self.logger.info(f"压缩VCF文件|Compressing VCF file with bgzip...")
                from .utils import compress_vcf_with_bgzip
                compress_vcf_with_bgzip(actual_output_path, self.output_vcf, self.logger)
                # 删除临时文件|Delete temp file
                import os
                os.unlink(actual_output_path)

            self.logger.info("")
            self.logger.info(f"转换完成|Conversion completed")
            self.logger.info(f"输出kmer数量|Kmers output: {output_count:,}")

            return output_count

        except Exception as e:
            self.logger.error(f"转换失败|Conversion failed: {e}")
            import traceback
            self.logger.error(f"错误详情|Error details:\n{traceback.format_exc()}")
            raise

    def _open_input_file(self) -> TextIO:
        """
        打开输入文件（支持gzip）|Open input file (supports gzip)

        Returns:
            TextIO: 文件对象|File object
        """
        return open_input(self.input_matrix, self.logger)

    def _convert_fast(self, input_file: TextIO, output_file: TextIO) -> int:
        """
        快速转换核心逻辑|Fast conversion core logic

        Args:
            input_file: 输入文件对象|Input file object
            output_file: 输出文件对象|Output file object

        Returns:
            int: 输出的kmer数量|Number of kmers output
        """
        # 读取第一行|Read first line
        first_line = input_file.readline()

        if not first_line:
            raise ValueError("输入文件为空|Input file is empty")

        # 解析样本名|Parse sample names
        samples = self._parse_header(first_line)
        num_samples = len(samples)

        self.logger.info(f"样本数量|Number of samples: {num_samples}")

        # 计算染色体数量和长度|Calculate chromosome count and length
        # 在no_header模式下，需要减1（因为第一行是样品名，不是数据）
        # In no_header mode, subtract 1 (first line is sample names, not data)
        total_kmers = self._estimate_total_kmers(input_file)
        if self.no_header:
            total_kmers = max(total_kmers - 1, 1)  # 至少1个kmer|At least 1 kmer

        # 根据用户指定的参数计算染色体分配|Calculate chromosome distribution based on user parameters
        if self.chr_number > 0:
            # 用户指定染色体数量|User specified chromosome count
            chr_count = self.chr_number
            chr_length = total_kmers // chr_count if total_kmers >= chr_count else 1
            remainder = total_kmers % chr_count
            last_chrom_length = chr_length + remainder if remainder > 0 else chr_length

            # 使用动态分配模式|Use dynamic allocation mode
            # 不预先计算边界，而是在转换过程中按顺序分配|Don't pre-calculate boundaries, allocate sequentially during conversion
            self.chr_boundaries = []  # 用于记录实际每条染色体的kmer数量|Record actual kmer count per chromosome
            self.chr_allocation_mode = 'dynamic'

            self.logger.info(f"预估kmer总数|Estimated total kmers: {total_kmers:,}")
            self.logger.info(f"用户指定染色体数量|User-specified chromosome count: {chr_count}")
            self.logger.info(f"每条染色体平均长度|Average length per chromosome: {chr_length:,}")
        else:
            # 用户指定每条染色体长度（默认模式）|User specified chromosome length (default mode)
            chr_length = self.chr_length
            chr_count = (total_kmers // chr_length) + 1
            last_chrom_length = total_kmers % chr_length
            self.chr_allocation_mode = 'fixed'  # 固定长度模式|Fixed length mode
            self.logger.info(f"预估kmer总数|Estimated total kmers: {total_kmers:,}")
            self.logger.info(f"每条染色体长度|Chromosome length: {chr_length:,}")
            self.logger.info(f"染色体数量|Number of chromosomes: {chr_count}")

        self.logger.info("")

        # 写入VCF header|Write VCF header
        self._write_vcf_header(output_file, samples, chr_count, last_chrom_length, chr_length)

        # 不需要重置文件位置，继续从当前位置读取
        # No need to reset file position, continue reading from current position
        # 因为第一行已经被读取并解析为样本名了
        # Because first line has been read and parsed as sample names

        # 转换数据行|Convert data lines
        output_count = 0
        count = 0
        current_chr_pos = 0  # 当前染色体上的位置|Position on current chromosome
        current_chr_num = 1  # 当前染色体编号|Current chromosome number

        self.logger.info("开始转换数据行|Starting data line conversion...")

        for line in input_file:
            count += 1

            # 进度显示|Progress display
            if count % 1000000 == 0:
                self.logger.info(f"已处理|Processed: {count:,} kmers, 输出|output: {output_count:,}")

            # 解析数据行|Parse data line
            parts = line.strip().split()
            if len(parts) < 2:
                continue

            kmer_id = parts[0]
            abundances = parts[1:]

            # 可选：简单过滤|Optional: simple filtering
            if self.min_freq > 0:
                freq = sum(1 for x in abundances if x == '1')
                if freq < self.min_freq:
                    continue

            # 分配染色体位置|Assign chromosome position
            if self.chr_allocation_mode == 'dynamic':
                # 动态分配模式：按顺序平均分配到所有染色体|Dynamic allocation: distribute evenly across chromosomes
                # 计算当前kmer应该分配到哪条染色体|Calculate which chromosome current kmer should be assigned to
                target_chr_num = ((count - 1) % chr_count) + 1

                # 如果切换到新染色体，更新位置计数器|If switching to new chromosome, update position counter
                if target_chr_num != current_chr_num:
                    current_chr_num = target_chr_num
                    current_chr_pos = 0

                current_chr_pos += 1
                chromosome = f"chr{current_chr_num}"
                pos = current_chr_pos

                # 记录每条染色体的实际kmer数量|Record actual kmer count for each chromosome
                if len(self.chr_boundaries) < chr_count:
                    self.chr_boundaries.extend([0] * (chr_count - len(self.chr_boundaries)))
                self.chr_boundaries[current_chr_num - 1] = current_chr_pos
            else:
                # 使用固定长度模式|Use fixed length mode (for chr_length parameter)
                chrom_num = (count // chr_length) + 1
                pos = (count % chr_length) + 1
                chromosome = f"chr{chrom_num}"

            # 转换为VCF行|Convert to VCF line
            vcf_line = self._format_vcf_line(
                chromosome=chromosome,
                position=pos,
                kmer_id=kmer_id,
                abundances=abundances,
                num_samples=num_samples
            )

            output_file.write(vcf_line)
            output_count += 1

        self.logger.info(f"转换完成|Conversion completed")
        self.logger.info(f"总处理|Total processed: {count:,} kmers")
        self.logger.info(f"总输出|Total output: {output_count:,} kmers")

        return output_count

    def _parse_header(self, header_line: str) -> list:
        """
        解析header行，提取样本名|Parse header line, extract sample names

        Args:
            header_line: header行|Header line

        Returns:
            list: 样本名列表|List of sample names
        """
        parts = header_line.strip().split()

        # 如果是no_header模式，第一行就是样本名（不是kmer数据）
        # In no_header mode, first line is sample names (not kmer data)
        if self.no_header:
            return parts

        # 检查第一列是否是KMER标识符|Check if first column is KMER identifier
        # 跳过第一列（KMER/ID）|Skip first column (KMER/ID)
        if parts and parts[0].upper() in ['KMER', 'ID', 'KMERID', 'KMER_ID']:
            return parts[1:]

        # 如果没有KMER标识符，假设第一列就是kmer（不是header）
        # If no KMER identifier, assume first column is kmer (not header)
        # 这种情况下，所有列都是样本（包括第一列的kmer）
        # This is problematic - we need to distinguish between:
        #   1. "KMER Sample1 Sample2..." (header line)
        #   2. "AAAAAAAAAAAA 0 0 1..." (data line without header)

        # 尝试检测：如果第一列是纯ATCGN序列且很长（>20bp），可能是数据行
        # Try to detect: if first column is pure ATCGN and long (>20bp), might be data line
        if parts and len(parts[0]) > 20 and all(c in 'ATCGN' for c in parts[0].upper()):
            # 这是数据行，不是header - 这种情况用户应该使用--no-header
            # This is a data line, not header - user should use --no-header
            self.logger.warning("")
            self.logger.warning("=" * 100)
            self.logger.warning(" 警告|WARNING")
            self.logger.warning("=" * 100)
            self.logger.warning("检测到输入文件第一行可能是数据行（而非header行）|Detected that first line might be data (not header)")
            self.logger.warning("第一列是kmer序列|First column appears to be kmer sequence")
            self.logger.warning("")
            self.logger.warning("请检查输入文件格式|Please check input file format:")
            self.logger.warning("  - 如果文件有header行（以KMER开头），请确保第一行是 'KMER Sample1 Sample2...'|If file has header (starts with KMER), ensure first line is 'KMER Sample1 Sample2...'")
            self.logger.warning("  - 如果文件没有header行，请使用 --no-header 参数|If file has no header, please use --no-header parameter")
            self.logger.warning("  - If file has header (starts with KMER), ensure first line is 'KMER Sample1 Sample2...'")
            self.logger.warning("  - If file has no header, please use --no-header parameter")
            self.logger.warning("")
            self.logger.warning("尝试将第一行当作样本名处理...|Attempting to treat first line as sample names...")
            self.logger.warning("=" * 100)
            self.logger.warning("")

        # 返回所有列作为样本名（包括第一列）|Return all columns as sample names (including first)
        return parts

    def _write_vcf_header(self, output_file: TextIO, samples: list,
                         chr_count: int, last_chrom_length: int, chr_length: int):
        """
        写入VCF header|Write VCF header

        Args:
            output_file: 输出文件对象|Output file object
            samples: 样本名列表|List of sample names
            chr_count: 染色体数量|Number of chromosomes
            last_chrom_length: 最后一条染色体长度|Last chromosome length
            chr_length: 每条染色体长度|Length of each chromosome
        """
        # VCF格式标准header|Standard VCF format header
        output_file.write('##fileformat=VCFv4.1\n')
        output_file.write(f'##INFO=<ID=KL,Number=1,Type=Integer,Description="Kmer length ({self.kmer_length})">\n')
        output_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # 写入contig信息|Write contig information
        if self.chr_allocation_mode == 'dynamic':
            # 动态分配模式：使用平均长度|Dynamic mode: use average length
            for i in range(1, chr_count + 1):
                output_file.write(f'##contig=<ID=chr{i},length={chr_length}>\n')
        else:
            # 固定长度模式|Fixed length mode
            for i in range(1, chr_count):
                output_file.write(f'##contig=<ID=chr{i},length={chr_length}>\n')
            output_file.write(f'##contig=<ID=chr{chr_count},length={last_chrom_length}>\n')

        # 写入column header|Write column header
        output_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')

    def _format_vcf_line(self, chromosome: str, position: int, kmer_id: str,
                        abundances: list, num_samples: int) -> str:
        """
        格式化VCF数据行|Format VCF data line

        Args:
            chromosome: 染色体|Chromosome
            position: 位置|Position
            kmer_id: Kmer ID
            abundances: 丰度列表|Abundance list
            num_samples: 样本数量|Number of samples

        Returns:
            str: VCF格式行|VCF format line
        """
        # VCF固定字段|VCF fixed fields
        ref = 'A'
        alt = 'G'
        qual = '100'
        filter_str = 'PASS'
        info = f'KL={self.kmer_length}'
        format_gt = 'GT'

        # 转换丰度为基因型|Convert abundance to genotype
        # 0 → 0/0, 1 → 1/1
        sample_gt = [f'{char}/{char}' for char in abundances[:num_samples]]

        # 组装VCF行|Assemble VCF line
        line = '\t'.join([
            chromosome,
            str(position),
            kmer_id,
            ref,
            alt,
            qual,
            filter_str,
            info,
            format_gt
        ]) + '\t' + '\t'.join(sample_gt) + '\n'

        return line

    def _estimate_total_kmers(self, input_file: TextIO) -> int:
        """
        估算总kmer数量（通过快速扫描）|Estimate total kmers (quick scan)

        Args:
            input_file: 输入文件对象|Input file object

        Returns:
            int: 估算的kmer数量|Estimated number of kmers
        """
        import os

        # 读取前几行获取平均行宽|Read first few lines to get average line width
        first_line = input_file.readline()
        if not first_line:
            return 100000000

        line_length = len(first_line.encode('utf-8'))

        # 基于压缩文件大小估算（对压缩文件估算值偏大，但不影响进度显示）|Estimate based on file size
        if hasattr(input_file, 'buffer') and hasattr(input_file.buffer, 'raw'):
            # zstd流式：用原始文件大小估算|zstd streaming: use raw file size
            raw_file = input_file.buffer.raw
            if hasattr(raw_file, 'name'):
                file_size = os.path.getsize(raw_file.name)
                estimated_total = int(file_size / max(line_length, 1))
            else:
                estimated_total = 100000000
        else:
            try:
                current_pos = input_file.tell()
                input_file.seek(0, 2)
                file_size = input_file.tell()
                input_file.seek(current_pos)
                estimated_total = int(file_size / max(line_length, 1))
            except (io.UnsupportedOperation, OSError):
                estimated_total = 100000000

        return max(estimated_total, 100000000)
