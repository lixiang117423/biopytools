"""
K-mer查询Pipeline|K-mer Query Pipeline
"""

import sys
from pathlib import Path
from typing import Dict, List

from ..config import QueryConfig
from ..rocksdb.querier import RocksDBQuerier
from ..utils import run_command


class KmerQueryPipeline:
    """K-mer查询Pipeline|K-mer Query Pipeline"""

    def __init__(self, config: QueryConfig, logger):
        """
        初始化Pipeline|Initialize Pipeline

        Args:
            config: 查询配置|Query configuration
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def run(self) -> Dict:
        """
        运行完整的查询流程|Run complete query pipeline

        Returns:
            dict: 查询结果统计|Query result statistics
        """
        self.logger.info("=" * 60)
        self.logger.info("开始K-mer查询|Starting K-mer query")
        self.logger.info("=" * 60)

        if self.config.use_kmtricks:
            return self._run_rocksdb_pipeline()
        else:
            return self._run_kmindex_pipeline()

    def _run_rocksdb_pipeline(self) -> Dict:
        """
        运行RocksDB查询流程|Run RocksDB query pipeline

        Returns:
            dict: 查询结果统计|Query result statistics
        """
        try:
            # 步骤1: 从FASTA提取k-mer|Step 1: Extract kmers from FASTA
            self.logger.info("[步骤1/4] 从FASTA提取k-mer|[Step 1/4] Extracting kmers from FASTA")
            kmer_file, pos_file = self._extract_kmers_from_fasta()
            if not kmer_file:
                raise Exception("从FASTA提取k-mer失败|Failed to extract kmers from FASTA")

            # 步骤2: 查询RocksDB|Step 2: Query RocksDB
            self.logger.info("[步骤2/4] 查询RocksDB|[Step 2/4] Querying RocksDB")
            query_result_file = self.config.output_path / "query_results.txt"
            if not self._query_rocksdb(kmer_file, query_result_file):
                raise Exception("查询RocksDB失败|Failed to query RocksDB")

            # 步骤3: 生成矩阵文件|Step 3: Generate matrix file
            self.logger.info("[步骤3/4] 生成矩阵文件|[Step 3/4] Generating matrix file")
            matrix_file = self.config.output_path / "kmer_matrix.txt"
            if not self._generate_matrix(query_result_file, pos_file, matrix_file):
                raise Exception("生成矩阵文件失败|Failed to generate matrix file")

            # 步骤4: 生成位置丰度文件|Step 4: Generate position-abundance file
            self.logger.info("[步骤4/4] 生成位置丰度文件|[Step 4/4] Generating position-abundance file")
            position_abundance_file = self.config.output_path / "kmer_position_abundance.txt"
            if not self._generate_position_abundance(query_result_file, position_abundance_file):
                raise Exception("生成位置丰度文件失败|Failed to generate position-abundance file")

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("K-mer查询完成 (RocksDB模式)|K-mer query completed (RocksDB mode)")
            self.logger.info("=" * 60)
            self.logger.info(f"k-mer文件|K-mer file: {kmer_file}")
            self.logger.info(f"查询结果|Query results: {query_result_file}")
            self.logger.info(f"矩阵文件|Matrix file: {matrix_file}")
            self.logger.info(f"位置丰度文件|Position-abundance file: {position_abundance_file}")

            return {
                'success': True,
                'mode': 'kmtricks',
                'kmer_file': str(kmer_file),
                'query_result': str(query_result_file),
                'matrix_file': str(matrix_file),
                'position_abundance_file': str(position_abundance_file)
            }

        except Exception as e:
            self.logger.error(f"查询失败|Query failed: {e}")
            import traceback
            traceback.print_exc()
            return {
                'success': False,
                'error': str(e)
            }

    def _run_kmindex_pipeline(self) -> Dict:
        """
        运行kmindex查询流程|Run kmindex query pipeline

        Returns:
            dict: 查询结果统计|Query result statistics
        """
        try:
            # 步骤1: 将FASTA转换为kmindex查询格式|Step 1: Convert FASTA to kmindex query format
            self.logger.info("[步骤1/2] 转换FASTA格式|[Step 1/2] Converting FASTA format")
            query_fasta = self._convert_fasta_for_kmindex()
            if not query_fasta:
                raise Exception("转换FASTA格式失败|Failed to convert FASTA format")

            # 步骤2: 运行kmindex query|Step 2: Run kmindex query
            self.logger.info("[步骤2/2] 运行kmindex query|[Step 2/2] Running kmindex query")
            query_output_dir = self.config.output_path / "kmindex_results"
            if not self._run_kmindex_query(query_fasta, query_output_dir):
                raise Exception("kmindex query运行失败|kmindex query failed")

            # 步骤3: 转换输出为热图格式|Step 3: Convert output to heatmap format
            self.logger.info("[步骤3/3] 转换输出格式|[Step 3/3] Converting output format")
            tsv_file = list(query_output_dir.glob("*.tsv"))[0] if query_output_dir.glob("*.tsv") else None
            if not tsv_file:
                raise Exception("未找到kmindex查询结果|kmindex query results not found")

            heatmap_file = self.config.output_path / "kmer_matrix.heatmap.txt"
            if not self._convert_kmindex_to_heatmap(tsv_file, heatmap_file):
                raise Exception("转换热图格式失败|Failed to convert to heatmap format")

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("K-mer查询完成 (kmindex模式)|K-mer query completed (kmindex mode)")
            self.logger.info("=" * 60)
            self.logger.info(f"查询FASTA|Query FASTA: {query_fasta}")
            self.logger.info(f"查询结果目录|Query results directory: {query_output_dir}")
            self.logger.info(f"热图文件|Heatmap file: {heatmap_file}")

            return {
                'success': True,
                'mode': 'kmindex',
                'query_fasta': str(query_fasta),
                'query_output_dir': str(query_output_dir),
                'heatmap_file': str(heatmap_file)
            }

        except Exception as e:
            self.logger.error(f"\n查询失败|Query failed: {e}")
            import traceback
            traceback.print_exc()
            return {
                'success': False,
                'error': str(e)
            }

    def _extract_kmers_from_fasta(self) -> tuple:
        """
        从FASTA文件提取k-mer|Extract kmers from FASTA file

        Returns:
            tuple: (kmer文件路径, 位置文件路径)|(kmer file path, position file path)
        """
        try:
            import pyfastx
        except ImportError:
            self.logger.error("未安装pyfastx|pyfastx not installed. 请先安装|Please install: pip install pyfastx")
            return None, None

        kmer_file = self.config.output_path / "kmers.txt"
        pos_file = self.config.output_path / "kmer_positions.txt"

        try:
            kmer_size = self.config.kmer_size

            with open(kmer_file, 'w') as outfile, open(pos_file, 'w') as posfile:
                outfile.write('kmer\n')
                posfile.write('kmer\tsample\tPos\n')

                kmers_seen = set()

                for name, seq in pyfastx.Fasta(str(self.config.query_fasta_path), build_index=False):
                    self.logger.debug(f"处理序列|Processing sequence: {name}")

                    # 提取k-mer|Extract kmers
                    for i in range(len(seq) - kmer_size + 1):
                        kmer = seq[i:i + kmer_size].upper()

                        # 跳过包含N的k-mer|Skip kmers containing N
                        if 'N' in kmer:
                            continue

                        # 计算canonical k-mer|Calculate canonical k-mer
                        canonical_kmer = self._get_canonical_kmer(kmer)

                        # 记录位置|Record position
                        posfile.write(f"{canonical_kmer}\t{name}\t{i + 1}\n")

                        # 记录k-mer（去重）|Record k-mer (unique)
                        if canonical_kmer not in kmers_seen:
                            outfile.write(f"{canonical_kmer}\n")
                            kmers_seen.add(canonical_kmer)

            self.logger.info(f"k-mer提取完成|K-mer extraction completed")
            self.logger.info(f"唯一k-mer数量|Unique k-mer count: {len(kmers_seen)}")

            return kmer_file, pos_file

        except Exception as e:
            self.logger.error(f"提取k-mer时出错|Error extracting kmers: {e}")
            import traceback
            traceback.print_exc()
            return None, None

    def _get_canonical_kmer(self, kmer: str) -> str:
        """
        计算canonical k-mer|Calculate canonical k-mer

        Args:
            kmer: k-mer序列|K-mer sequence

        Returns:
            str: canonical k-mer（字典序较小的）|Canonical k-mer (lexicographically smaller)
        """
        rc_kmer = self._reverse_complement(kmer)
        return min(kmer, rc_kmer)

    def _reverse_complement(self, sequence: str) -> str:
        """
        计算反向互补序列|Calculate reverse complement

        Args:
            sequence: DNA序列|DNA sequence

        Returns:
            str: 反向互补序列|Reverse complement sequence
        """
        complement_map = str.maketrans("ATCG", "TAGC")
        complement_seq = sequence.translate(complement_map)
        return complement_seq[::-1]

    def _query_rocksdb(self, kmer_file: Path, output_file: Path) -> bool:
        """
        查询RocksDB|Query RocksDB

        Args:
            kmer_file: k-mer文件路径|K-mer file path
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 创建查询器|Create querier
            querier = RocksDBQuerier(
                db_path=str(self.config.rocksdb_path),
                bloom_bits=15,
                header_db_key=self.config.header_db_key,
                logger=self.logger
            )

            # 打开数据库|Open database
            if not querier.open_db():
                return False

            # 读取k-mer列表|Read k-mer list
            kmer_list = []
            with open(kmer_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and line != 'kmer':  # 跳过表头|Skip header
                        kmer_list.append(line)

            self.logger.info(f"共|Total {len(kmer_list)} 个k-mer待查询|kmers to query")

            # 执行查询|Execute query
            results = querier.query_kmers(kmer_list, include_rc=True)

            # 写入结果|Write results
            querier.write_results(results, str(output_file))

            # 关闭数据库|Close database
            querier.close()

            return True

        except Exception as e:
            self.logger.error(f"查询RocksDB时出错|Error querying RocksDB: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _generate_matrix(self, query_result: Path, pos_file: Path, output_file: Path) -> bool:
        """
        生成矩阵文件|Generate matrix file

        Args:
            query_result: 查询结果文件|Query result file
            pos_file: 位置文件|Position file
            output_file: 输出文件|Output file

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 读取位置信息|Read position info
            pos_dict = {}
            with open(pos_file, 'r') as f:
                for line in f:
                    if line.startswith('kmer'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        kmer, sample, pos = parts[0], parts[1], parts[2]
                        key = f"{kmer}_{sample}_{pos}"
                        pos_dict[key] = line.strip()

            # 处理查询结果并添加位置信息|Process query results and add position info
            processed_lines = []
            with open(query_result, 'r') as f:
                for line in f:
                    if line.startswith('kmer') or line.startswith('ID'):
                        # 更改表头|Change header
                        line = line.replace('ID', 'kmer')
                        processed_lines.append(line)
                    else:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            kmer = parts[0]
                            abundance = parts[1]

                            # 查找匹配的位置信息|Find matching position info
                            matched_key = None
                            for key in pos_dict:
                                if key.startswith(kmer):
                                    matched_key = key
                                    break

                            if matched_key:
                                pos_line = pos_dict[matched_key]
                                combined = f"{line.strip()}\t{pos_line}"
                                processed_lines.append(combined)

            # 写入处理后的文件|Write processed file
            temp_file = output_file.parent / "temp_combined.txt"
            with open(temp_file, 'w') as f:
                f.write('\n'.join(processed_lines))

            # 转换为最终矩阵格式|Convert to final matrix format
            with open(output_file, 'w') as f:
                f.write('Gene_Pos\t')  # 新表头|New header

                # 写入第一个有效行的列名|Write column names from first valid line
                first_line = None
                for line in processed_lines:
                    if line.startswith('kmer'):
                        parts = line.split('\t')
                        if len(parts) > 4:
                            f.write('\t'.join(parts[4:]) + '\n')
                            first_line = line
                            break

                # 写入数据行|Write data rows
                for line in processed_lines:
                    if line.startswith('kmer'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        gene = parts[1]
                        pos = parts[2]
                        abundance = parts[4:]
                        f.write(f"{gene}_{pos}\t{' '.join(abundance)}\n")

            # 删除临时文件|Delete temp file
            temp_file.unlink()

            self.logger.info(f"矩阵文件生成完成|Matrix file generated: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成矩阵文件时出错|Error generating matrix file: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _generate_position_abundance(self, query_result_file: Path, output_file: Path) -> bool:
        """
        生成位置丰度文件（BED格式+样本丰度）|Generate position-abundance file (BED format + sample abundances)

        使用left join逻辑：以BED文件为主表，query_results.txt为从表
        Use left join logic: BED file as primary table, query_results.txt as secondary table

        Args:
            query_result_file: 查询结果文件|Query result file
            output_file: 输出文件|Output file

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 检查是否提供了BED文件|Check if BED file is provided
            if not hasattr(self.config, 'bed_file') or not self.config.bed_file:
                self.logger.warning("未提供BED文件，跳过位置丰度文件生成|No BED file provided, skipping position-abundance file generation")
                # 仍然创建一个空文件以保持兼容性|Still create an empty file for compatibility
                output_file.touch()
                return True

            bed_path = self.config.bed_file_path if hasattr(self.config, 'bed_file_path') else Path(self.config.bed_file)
            if not bed_path.exists():
                self.logger.error(f"BED文件不存在|BED file does not exist: {bed_path}")
                return False

            # 步骤1: 读取query_results.txt，建立canonical kmer到丰度的映射|Step 1: Read query_results.txt, build canonical kmer to abundance mapping
            kmer_to_abundance = {}
            sample_names = []

            with open(query_result_file, 'r') as f:
                # 读取表头|Read header
                header = f.readline().strip()
                if header.startswith('kmer') or header.startswith('ID'):
                    header_parts = header.split('\t')
                    sample_names = header_parts[1:] if len(header_parts) > 1 else []

                # 读取数据行|Read data rows
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 2:
                        continue

                    kmer_seq = parts[0]
                    abundances = parts[1:] if len(parts) > 1 else []

                    # 将"-"替换为"0"|Replace "-" with "0"
                    abundances = [ab.replace('-', '0') for ab in abundances]

                    kmer_to_abundance[kmer_seq] = abundances

            self.logger.info(f"从query_results.txt读取了|Read from query_results.txt: {len(kmer_to_abundance)} 个kmer")

            # 步骤2: 读取BED文件，将kmer转换为canonical形式，建立原始kmer到位置信息的映射|Step 2: Read BED file, convert kmer to canonical form
            bed_records = []
            with open(bed_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        chrom = parts[0]
                        start = parts[1]
                        end = parts[2]
                        kmer_seq = parts[3]
                        # 转换为canonical kmer|Convert to canonical kmer
                        canonical_kmer = self._get_canonical_kmer(kmer_seq)
                        bed_records.append({
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'original_kmer': kmer_seq,
                            'canonical_kmer': canonical_kmer
                        })

            self.logger.info(f"从BED文件读取了|Read from BED file: {len(bed_records)} 条记录")

            # 步骤3: Left join - 以BED文件为主，生成输出文件|Step 3: Left join - BED file as primary, generate output file
            with open(output_file, 'w') as f_out:
                # 写入表头|Write header
                f_out.write('Chr\tStart\tEnd\tKmer\t' + '\t'.join(sample_names) + '\n')

                # 统计|Statistics
                matched_count = 0
                unmatched_count = 0

                # 遍历BED文件中的每一条记录|Iterate through each record in BED file
                for record in bed_records:
                    chrom = record['chrom']
                    start = record['start']
                    end = record['end']
                    original_kmer = record['original_kmer']
                    canonical_kmer = record['canonical_kmer']

                    # 在query_results中查找canonical kmer|Find canonical kmer in query_results
                    if canonical_kmer in kmer_to_abundance:
                        abundances = kmer_to_abundance[canonical_kmer]
                        f_out.write(f'{chrom}\t{start}\t{end}\t{original_kmer}\t' + '\t'.join(abundances) + '\n')
                        matched_count += 1
                    else:
                        # 未找到匹配，填充0|Not found, fill with 0
                        zero_abundances = ['0'] * len(sample_names)
                        f_out.write(f'{chrom}\t{start}\t{end}\t{original_kmer}\t' + '\t'.join(zero_abundances) + '\n')
                        unmatched_count += 1

            self.logger.info(f"位置丰度文件生成完成|Position-abundance file generated: {output_file}")
            self.logger.info(f"匹配的k-mer数量|Matched k-mer count: {matched_count}")
            self.logger.info(f"未匹配的k-mer数量|Unmatched k-mer count: {unmatched_count}")
            self.logger.info(f"匹配率|Match rate: {matched_count}/{len(bed_records)} ({100*matched_count/len(bed_records):.1f}%)")

            return True

        except Exception as e:
            self.logger.error(f"生成位置丰度文件时出错|Error generating position-abundance file: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _convert_fasta_for_kmindex(self) -> Path:
        """
        转换FASTA为kmindex查询格式|Convert FASTA to kmindex query format

        Returns:
            Path: 转换后的FASTA文件路径|Converted FASTA file path
        """
        try:
            import pyfastx
        except ImportError:
            self.logger.error("未安装pyfastx|pyfastx not installed. 请先安装|Please install: pip install pyfastx")
            return None

        output_fasta = self.config.output_path / "query_kmindex.fasta"

        try:
            kmer_size = self.config.kmer_size

            with open(output_fasta, 'w') as outfile:
                for name, seq in pyfastx.Fasta(str(self.config.query_fasta_path), build_index=False):
                    self.logger.debug(f"处理序列|Processing sequence: {name}")

                    # 提取k-mer|Extract kmers
                    for i in range(len(seq) - kmer_size + 1):
                        kmer = seq[i:i + kmer_size].upper()

                        # 跳过包含N的k-mer|Skip kmers containing N
                        if 'N' in kmer:
                            continue

                        # 计算canonical k-mer|Calculate canonical k-mer
                        canonical_kmer = self._get_canonical_kmer(kmer)

                        # kmindex格式: >name_pos kmer
                        outfile.write(f">{name}_{i+1}\n{canonical_kmer}\n")

            self.logger.info(f"FASTA格式转换完成|FASTA conversion completed: {output_fasta}")
            return output_fasta

        except Exception as e:
            self.logger.error(f"转换FASTA格式时出错|Error converting FASTA format: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _run_kmindex_query(self, query_fasta: Path, output_dir: Path) -> bool:
        """
        运行kmindex query|Run kmindex query

        Args:
            query_fasta: 查询FASTA文件|Query FASTA file
            output_dir: 输出目录|Output directory

        Returns:
            bool: 是否成功|Success
        """
        try:
            import os

            # 创建输出目录|Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # 构建命令|Build command
            cmd = [
                self.config.kmindex_path,
                'query',
                '-i', str(self.config.index_path),
                '-q', str(query_fasta),
                '-n', self.config.index_name,
                '-z', str(self.config.zvalue),
                '-r', str(self.config.threshold),
                '-f', self.config.output_format,
                '-o', str(output_dir),
                '-t', str(self.config.threads)
            ]

            result = run_command(cmd, self.logger, check=True, capture_output=False)

            self.logger.info("kmindex query完成|kmindex query completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"kmindex query失败|kmindex query failed: {e}")
            return False
        except Exception as e:
            self.logger.error(f"运行kmindex query时出错|Error running kmindex query: {e}")
            return False

    def _convert_kmindex_to_heatmap(self, tsv_file: Path, output_file: Path) -> bool:
        """
        转换kmindex TSV输出为热图格式|Convert kmindex TSV output to heatmap format

        Args:
            tsv_file: kmindex输出的TSV文件|kmindex output TSV file
            output_file: 输出热图文件|Output heatmap file

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 用Python转置矩阵(替代awk+sed,避免shell引号/管道问题,§13)|Transpose matrix in Python (replaces awk+sed, avoids shell quoting/pipe issues)
            self.logger.info(f"转置矩阵生成热图|Transposing matrix for heatmap: {tsv_file} -> {output_file}")
            with open(tsv_file, 'r') as f_in:
                # 按空白分隔字段(等价awk默认FS)|split fields on whitespace (equivalent to awk default FS)
                rows = [line.split() for line in f_in if line.strip()]

            # 转置: 行<->列; 空矩阵则生成空输出|Transpose rows<->cols; empty matrix yields empty output
            transposed = list(zip(*rows))
            with open(output_file, 'w') as f_out:
                for col in transposed:
                    f_out.write('\t'.join(col) + '\n')

            self.logger.info(f"热图文件生成完成|Heatmap file generated: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"转换热图格式时出错|Error converting heatmap format: {e}")
            return False
