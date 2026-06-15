"""
VCF抽样核心逻辑模块|VCF Sampling Core Logic Module
"""

import gzip
import os
import random
import subprocess
from collections import defaultdict
from typing import Dict, Set, Tuple, Optional

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


class VCFSamplerCore:
    """VCF抽样核心类|VCF Sampling Core Class"""

    def __init__(self, config, logger):
        """
        初始化抽样核心类|Initialize sampling core class

        Args:
            config: VCFSamplerConfig配置对象|VCFSamplerConfig object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

        # 设置随机种子|Set random seed
        if config.random_seed is not None:
            random.seed(config.random_seed)

        # 追踪是否自动创建了索引|Track if index was auto-created
        self.auto_created_index = False

    def _check_index_exists(self, vcf_file: str) -> bool:
        """
        检查VCF文件是否存在索引|Check if VCF file has index

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            是否存在索引|Whether index exists
        """
        index_file = f"{vcf_file}.tbi"
        if os.path.exists(index_file):
            return True

        # 也检查csi索引|Also check for csi index
        index_file_csi = f"{vcf_file}.csi"
        if os.path.exists(index_file_csi):
            return True

        return False

    def _create_index_with_tabix(self, vcf_file: str) -> bool:
        """
        使用tabix创建VCF索引|Create VCF index using tabix

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            是否成功创建索引|Whether index was successfully created
        """
        self.logger.info(f"VCF索引不存在，开始创建索引|VCF index not found, creating index...")
        self.logger.info(f"文件|File: {vcf_file}")

        try:
            # 使用tabix创建索引|Use tabix to create index
            cmd = ['tabix', '-p', 'vcf', vcf_file]

            self.logger.debug(f"执行命令|Running command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10分钟超时|10 minutes timeout
            )

            if result.returncode == 0:
                self.logger.info("索引创建成功|Index created successfully")
                return True
            else:
                self.logger.error(f"创建索引失败|Failed to create index: {result.stderr}")
                return False

        except FileNotFoundError:
            self.logger.warning(
                "tabix命令未找到，请安装htslib：| "
                "tabix command not found, please install htslib: "
                "conda install -c bioconda htslib"
            )
            return False
        except subprocess.TimeoutExpired:
            self.logger.error("创建索引超时（超过10分钟）| Index creation timed out (exceeded 10 minutes)")
            return False
        except Exception as e:
            self.logger.error(f"创建索引时出错|Error creating index: {str(e)}")
            return False

    def _ensure_index_exists(self, vcf_file: str) -> bool:
        """
        确保VCF文件有索引，如果没有则自动创建|Ensure VCF has index, create if missing

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            是否有可用的索引|Whether index is available
        """
        if self._check_index_exists(vcf_file):
            self.logger.info("检测到VCF索引|VCF index detected")
            return True

        # 尝试自动创建索引|Try to auto-create index
        self.auto_created_index = self._create_index_with_tabix(vcf_file)

        if not self.auto_created_index:
            self.logger.warning(
                "无法创建VCF索引，将使用Python模式（较慢）| "
                "Cannot create VCF index, will use Python mode (slow)"
            )

        return self.auto_created_index

    def _cleanup_auto_created_index(self):
        """
        清理自动创建的索引|Clean up auto-created index
        """
        if self.auto_created_index:
            index_file = f"{self.config.input_vcf}.tbi"
            try:
                if os.path.exists(index_file):
                    os.remove(index_file)
                    self.logger.info(f"已删除自动创建的索引|Removed auto-created index: {index_file}")
            except Exception as e:
                self.logger.warning(f"删除索引失败|Failed to remove index: {str(e)}")

    def count_snp_by_chromosome(self) -> Dict[str, int]:
        """
        统计每条染色体的SNP数量|Count SNP number by chromosome

        优先使用pysam加速，如果不可用则使用Python逐行读取|Use pysam for speed if available,
        otherwise fall back to Python line-by-line reading

        Returns:
            染色体到SNP数量的映射|Mapping from chromosome to SNP count
        """
        self.logger.info("开始统计每条染色体的SNP数量|Start counting SNPs by chromosome")

        if PYSAM_AVAILABLE:
            return self._count_snp_with_pysam()
        else:
            self.logger.warning(
                "pysam未安装，使用Python逐行读取（建议安装pysam以加速）| "
                "pysam not installed, using Python line-by-line reading "
                "(recommend installing pysam for speed)"
            )
            return self._count_snp_with_python()

    def _count_snp_with_pysam(self) -> Dict[str, int]:
        """
        使用pysam快速统计SNP数量|Count SNPs quickly using pysam

        Returns:
            染色体到SNP数量的映射|Mapping from chromosome to SNP count
        """
        # 确保索引存在|Ensure index exists
        if not self._ensure_index_exists(self.config.input_vcf):
            # 无法创建索引，回退到Python方法|Cannot create index, fall back to Python
            return self._count_snp_with_python()

        chr_counts = defaultdict(int)

        try:
            with pysam.VariantFile(self.config.input_vcf) as vcf:
                # 获取所有染色体|Get all chromosomes
                chromosomes = list(vcf.header.contigs)

                self.logger.info(f"检测到|Detected {len(chromosomes)} 条染色体|chromosomes")

                # 为每条染色体统计记录数|Count records for each chromosome
                for chrom in chromosomes:
                    try:
                        # 使用fetch遍历并计数|Use fetch to iterate and count
                        count = 0
                        for record in vcf.fetch(contig=chrom):
                            count += 1
                        chr_counts[chrom] = count
                        self.logger.debug(f"染色体|Chromosome {chrom}: {count:,} 个SNP|SNPs")
                    except Exception as e:
                        self.logger.warning(
                            f"无法统计染色体|Failed to count chromosome {chrom}: {str(e)}"
                        )

        except Exception as e:
            self.logger.error(
                f"使用pysam统计失败，回退到Python方法|"
                f"pysam counting failed, falling back to Python method: {str(e)}"
            )
            return self._count_snp_with_python()

        total_snp = sum(chr_counts.values())
        self.logger.info(f"总SNP数量|Total SNP count: {total_snp:,}")
        self.logger.info(f"染色体数量|Chromosome count: {len(chr_counts)}")

        return chr_counts

    def _count_snp_with_python(self) -> Dict[str, int]:
        """
        使用Python逐行统计SNP数量（较慢）| Count SNPs using Python line-by-line (slow)

        Returns:
            染色体到SNP数量的映射|Mapping from chromosome to SNP count
        """
        chr_counts = defaultdict(int)

        with gzip.open(self.config.input_vcf, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                if not line.startswith('#'):
                    chrom = line.split('\t')[0]
                    chr_counts[chrom] += 1

                # 每处理100万行输出一次进度|Output progress every 1M lines
                if line_num % 1000000 == 0:
                    self.logger.debug(f"已处理|Processed: {line_num:,} 行|lines")

        total_snp = sum(chr_counts.values())
        self.logger.info(f"总SNP数量|Total SNP count: {total_snp:,}")
        self.logger.info(f"染色体数量|Chromosome count: {len(chr_counts)}")

        return chr_counts

    def select_snp_indices(self, chr_counts: Dict[str, int]) -> Dict[str, Set[int]]:
        """
        选择要抽取的SNP索引|Select SNP indices to sample

        Args:
            chr_counts: 染色体到SNP数量的映射|Mapping from chromosome to SNP count

        Returns:
            染色体到选中索引集合的映射|Mapping from chromosome to selected index set
        """
        self.logger.info(
            f"开始随机抽取SNP (抽样比例|sampling rate: {self.config.sample_rate:.1%})"
        )

        selected_indices = {}
        total_selected = 0

        for chrom, count in chr_counts.items():
            n_select = int(count * self.config.sample_rate)

            # 确保至少选择1个SNP|Ensure at least 1 SNP is selected
            if n_select == 0 and count > 0:
                n_select = 1

            selected_indices[chrom] = set(random.sample(range(count), n_select))
            total_selected += n_select

        self.logger.info(f"总共选择|Total selected: {total_selected:,} 个SNP|SNPs")

        return selected_indices

    def write_sampled_vcf(self, selected_indices: Dict[str, Set[int]]) -> Tuple[int, int]:
        """
        写入抽样的VCF文件|Write sampled VCF file

        优先使用pysam加速，如果不可用则使用Python逐行写入 |
        Use pysam for speed if available, otherwise fall back to Python

        Args:
            selected_indices: 染色体到选中索引集合的映射|Mapping from chromosome to selected index set

        Returns:
            (总SNP数, 写入的SNP数)|(total SNP count, written SNP count)
        """
        self.logger.info("开始写入抽样的VCF文件|Start writing sampled VCF file")

        if PYSAM_AVAILABLE:
            return self._write_sampled_vcf_with_pysam(selected_indices)
        else:
            return self._write_sampled_vcf_with_python(selected_indices)

    def _write_sampled_vcf_with_pysam(
        self, selected_indices: Dict[str, Set[int]]
    ) -> Tuple[int, int]:
        """
        使用pysam写入抽样的VCF文件|Write sampled VCF using pysam

        Args:
            selected_indices: 染色体到选中索引集合的映射|Mapping from chromosome to selected index set

        Returns:
            (总SNP数, 写入的SNP数)|(total SNP count, written SNP count)
        """
        total_snp = 0
        written_snp = 0
        chr_current_idx = defaultdict(int)

        try:
            with pysam.VariantFile(self.config.input_vcf) as invcf:
                # 创建输出VCF文件|Create output VCF file
                with pysam.VariantFile(
                    self.config.output_vcf, 'wz', header=invcf.header
                ) as outvcf:

                    for chrom in selected_indices.keys():
                        chr_current_idx[chrom] = 0

                        # 使用pysam直接访问特定染色体的记录|Use pysam to access specific chromosome records
                        try:
                            for record in invcf.fetch(contig=chrom):
                                current_idx = chr_current_idx[chrom]

                                if current_idx in selected_indices[chrom]:
                                    outvcf.write(record)
                                    written_snp += 1

                                chr_current_idx[chrom] += 1
                                total_snp += 1

                        except ValueError as e:
                            # 染色体不存在于VCF中|Chromosome not in VCF
                            self.logger.warning(
                                f"染色体|Chromosome {chrom} 不存在或无法访问|"
                                f"does not exist or cannot be accessed: {str(e)}"
                            )
                            continue

        except Exception as e:
            self.logger.error(
                f"使用pysam写入失败，回退到Python方法|"
                f"pysam writing failed, falling back to Python method: {str(e)}"
            )
            return self._write_sampled_vcf_with_python(selected_indices)

        self.logger.info(f"总SNP数|Total SNP count: {total_snp:,}")
        self.logger.info(f"写入的SNP数|Written SNP count: {written_snp:,}")
        self.logger.info(
            f"实际抽样比例|Actual sampling rate: {written_snp/total_snp:.2%} "
            f"(目标|target: {self.config.sample_rate:.1%})"
        )
        self.logger.info(f"输出文件|Output file: {self.config.output_vcf}")

        # 创建索引|Create index
        try:
            self.logger.info("创建VCF索引|Creating VCF index...")
            pysam.tabix_index(self.config.output_vcf, preset='vcf')
            self.logger.info("索引创建成功|Index created successfully")
        except Exception as e:
            self.logger.warning(f"创建索引失败|Failed to create index: {str(e)}")

        # 清理自动创建的输入索引|Clean up auto-created input index
        self._cleanup_auto_created_index()

        return total_snp, written_snp

    def _write_sampled_vcf_with_python(
        self, selected_indices: Dict[str, Set[int]]
    ) -> Tuple[int, int]:
        """
        使用Python逐行写入抽样的VCF文件（较慢）| Write sampled VCF using Python (slow)

        Args:
            selected_indices: 染色体到选中索引集合的映射|Mapping from chromosome to selected index set

        Returns:
            (总SNP数, 写入的SNP数)|(total SNP count, written SNP count)
        """
        chr_current_idx = defaultdict(int)
        header_lines = []
        total_snp = 0
        written_snp = 0

        with gzip.open(self.config.input_vcf, 'rt') as fin, \
             gzip.open(self.config.output_vcf, 'wt') as fout:

            for line_num, line in enumerate(fin, 1):
                if line.startswith('#'):
                    # 保存头信息|Save header lines
                    if self.config.keep_header:
                        header_lines.append(line)
                        fout.write(line)
                else:
                    # 处理SNP行|Process SNP line
                    chrom = line.split('\t')[0]
                    current_idx = chr_current_idx[chrom]

                    if current_idx in selected_indices.get(chrom, set()):
                        fout.write(line)
                        written_snp += 1

                    chr_current_idx[chrom] += 1
                    total_snp += 1

                # 每处理100万行输出一次进度|Output progress every 1M lines
                if line_num % 1000000 == 0:
                    self.logger.debug(f"已处理|Processed: {line_num:,} 行|lines, "
                                    f"已写入|written: {written_snp:,} 个SNP|SNPs")

        self.logger.info(f"总SNP数|Total SNP count: {total_snp:,}")
        self.logger.info(f"写入的SNP数|Written SNP count: {written_snp:,}")
        self.logger.info(
            f"实际抽样比例|Actual sampling rate: {written_snp/total_snp:.2%} "
            f"(目标|target: {self.config.sample_rate:.1%})"
        )
        self.logger.info(f"输出文件|Output file: {self.config.output_vcf}")

        # 清理自动创建的输入索引|Clean up auto-created input index
        self._cleanup_auto_created_index()

        return total_snp, written_snp
