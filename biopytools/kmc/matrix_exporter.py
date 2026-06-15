"""
KMC矩阵导出模块|KMC Matrix Exporter Module
"""

import h5py
import numpy as np
from pathlib import Path
from typing import Optional

from .config import KMCConfig
from .utils import KMCLogger


class KMCMatrixExporter:
    """KMC矩阵导出器|KMC Matrix Exporter"""

    def __init__(self, config: KMCConfig, logger: Optional[KMCLogger] = None):
        """初始化导出器|Initialize exporter"""
        self.config = config
        self.config.validate()

        if logger is None:
            self.logger_manager = KMCLogger(self.config.output_path)
            self.logger = self.logger_manager.get_logger()
        else:
            self.logger_manager = logger
            self.logger = logger.get_logger()

        # HDF5文件路径|HDF5 file path
        self.matrix_file = self.config.output_path / 'abundance_matrix.h5'
        self.kmer_dict_file = self.config.output_path / 'kmer_dictionary.h5'

    def export_to_tsv(self, output_file: str, format: str = 'sparse',
                      min_abundance: int = 1) -> bool:
        """导出矩阵为TSV格式|Export matrix to TSV format

        Args:
            output_file: 输出TSV文件路径|Output TSV file path
            format: 输出格式 ('full' 或 'sparse')|Output format ('full' or 'sparse')
            min_abundance: 最小丰度阈值|Minimum abundance threshold

        Returns:
            是否成功|Success or not
        """
        self.logger.info(f"开始导出矩阵到TSV|Starting to export matrix to TSV: {output_file}")
        self.logger.info(f"输出格式|Output format: {format}")
        self.logger.info(f"最小丰度阈值|Minimum abundance threshold: {min_abundance}")

        # 检查文件是否存在|Check if files exist
        if not self.matrix_file.exists():
            self.logger.error(f"矩阵文件不存在|Matrix file not found: {self.matrix_file}")
            return False

        if not self.kmer_dict_file.exists():
            self.logger.error(f"k-mer字典文件不存在|K-mer dictionary file not found: {self.kmer_dict_file}")
            return False

        try:
            # 读取HDF5文件|Read HDF5 files
            with h5py.File(self.matrix_file, 'r') as matrix_f:
                with h5py.File(self.kmer_dict_file, 'r') as dict_f:
                    storage_type = matrix_f.attrs['storage_type']
                    # 安全地decode bytes对象|Safely decode bytes objects
                    sample_names = [s.decode('utf-8') if isinstance(s, bytes) else s for s in matrix_f['sample_names'][:]]
                    kmers = [k.decode('utf-8') if isinstance(k, bytes) else k for k in dict_f['kmer'][:]]
                    storage_type_str = storage_type.decode('utf-8') if isinstance(storage_type, bytes) else storage_type

                    self.logger.info(f"矩阵信息|Matrix info:")
                    self.logger.info(f"  样本数|Samples: {len(sample_names)}")
                    self.logger.info(f"  k-mer数|K-mers: {len(kmers)}")
                    self.logger.info(f"  存储类型|Storage type: {storage_type_str}")

                    if format == 'sparse':
                        self._export_sparse(matrix_f, dict_f, kmers, sample_names,
                                          output_file, min_abundance)
                    else:  # format == 'full'
                        self._export_full(matrix_f, dict_f, kmers, sample_names,
                                        output_file, min_abundance)

            self.logger.info(f"导出完成|Export completed: {output_file}")
            return True

        except Exception as e:
            self.logger.error(f"导出失败|Export failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return False

    def _export_sparse(self, matrix_f, dict_f, kmers: list, sample_names: list,
                       output_file: str, min_abundance: int):
        """导出为稀疏格式|Export to sparse format"""
        self.logger.info("使用稀疏格式导出|Exporting in sparse format")

        kmer_ids = matrix_f['kmer_id'][:]
        sample_ids = matrix_f['sample_id'][:]
        abundances = matrix_f['abundance'][:]

        # 过滤低丰度k-mer|Filter low abundance k-mers
        total_entries = len(abundances)
        if min_abundance > 1:
            mask = abundances >= min_abundance
            kmer_ids = kmer_ids[mask]
            sample_ids = sample_ids[mask]
            abundances = abundances[mask]
            filtered = total_entries - len(abundances)
            self.logger.info(f"过滤掉 {filtered} 个丰度<{min_abundance} 的条目|"
                           f"Filtered {filtered} entries with abundance<{min_abundance}")

        # 写入TSV|Write to TSV
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w') as out:
            # 写入表头|Write header
            header = ['kmer_id', 'kmer', 'sample_id', 'sample_name', 'abundance']
            out.write('\t'.join(header) + '\n')

            # 写入数据|Write data
            for kid, sid, ab in zip(kmer_ids, sample_ids, abundances):
                kmer = kmers[kid]
                sample = sample_names[sid]
                out.write(f'{kid}\t{kmer}\t{sid}\t{sample}\t{ab}\n')

        self.logger.info(f"写入 {len(abundances)} 条记录|Wrote {len(abundances)} entries")

    def _export_full(self, matrix_f, dict_f, kmers: list, sample_names: list,
                     output_file: str, min_abundance: int):
        """导出为完整矩阵格式|Export to full matrix format"""
        self.logger.info("使用完整矩阵格式导出|Exporting in full matrix format")

        storage_type = matrix_f.attrs['storage_type']
        # 兼容字符串和bytes|Compatible with both str and bytes
        if isinstance(storage_type, bytes):
            storage_type_str = storage_type.decode('utf-8')
        else:
            storage_type_str = storage_type

        if storage_type_str == 'dense':
            # 密集矩阵|Dense matrix
            matrix = matrix_f['abundance'][:]

            # 过滤低丰度k-mer（过滤掉整行都低于阈值的k-mer）
            if min_abundance > 1:
                row_max = matrix.max(axis=1)
                mask = row_max >= min_abundance
                matrix = matrix[mask]
                kmers = [kmers[i] for i in range(len(kmers)) if mask[i]]
                self.logger.info(f"过滤掉 {sum(~mask)} 个低丰度k-mer|"
                               f"Filtered {sum(~mask)} low abundance k-mers")

            # 写入TSV|Write to TSV
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, 'w') as out:
                # 写入表头|Write header
                out.write('kmer\t' + '\t'.join(sample_names) + '\n')

                # 写入数据|Write data
                for i, kmer in enumerate(kmers):
                    row = [kmer] + [str(matrix[i, j]) for j in range(len(sample_names))]
                    out.write('\t'.join(row) + '\n')

            self.logger.info(f"写入 {len(kmers)} 个k-mer x {len(sample_names)} 个样本|"
                           f"Wrote {len(kmers)} k-mers x {len(sample_names)} samples")

        else:  # sparse matrix
            # 稀疏矩阵转密集矩阵（仅适合小规模数据）|Sparse to dense (for small data only)
            self.logger.warning("从稀疏矩阵转换为完整矩阵可能需要大量内存和时间|"
                              "Converting sparse to full matrix may require lots of memory and time")

            n_kmers = len(kmers)
            n_samples = len(sample_names)

            # 检查规模|Check size
            estimated_size = n_kmers * n_samples * 2  # 每个元素2字节（uint16）
            estimated_size_gb = estimated_size / (1024**3)

            if estimated_size_gb > 4:  # 超过4GB警告|Warn if >4GB
                self.logger.warning(f"预估矩阵大小: {estimated_size_gb:.2f}GB，可能导致内存不足|"
                                  f"Estimated matrix size: {estimated_size_gb:.2f}GB, may cause OOM")

            # 构建密集矩阵|Build dense matrix
            matrix = np.zeros((n_kmers, n_samples), dtype=np.uint16)

            kmer_ids = matrix_f['kmer_id'][:]
            sample_ids = matrix_f['sample_id'][:]
            abundances = matrix_f['abundance'][:]

            for kid, sid, ab in zip(kmer_ids, sample_ids, abundances):
                if ab >= min_abundance:
                    matrix[kid, sid] = ab

            # 不过滤全零行，导出所有k-mer|Don't filter all-zero rows, export all k-mers
            # 过滤低丰度k-mer（可选）|Filter low abundance k-mers (optional)
            if min_abundance > 1:
                row_max = matrix.max(axis=1)
                mask = row_max >= min_abundance
                # 保留符合条件的k-mer|Keep matching k-mers
                matrix_filtered = matrix[mask]
                kmers_filtered = [kmers[i] for i in range(len(kmers)) if mask[i]]
                self.logger.info(f"过滤掉 {sum(~mask)} 个低丰度k-mer|"
                               f"Filtered {sum(~mask)} low abundance k-mers (<{min_abundance})")
                matrix = matrix_filtered
                kmers = kmers_filtered

            # 写入TSV|Write to TSV
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, 'w') as out:
                # 写入表头|Write header
                out.write('kmer\t' + '\t'.join(sample_names) + '\n')

                # 写入数据|Write data
                for i, kmer in enumerate(kmers):
                    row = [kmer] + [str(matrix[i, j]) for j in range(len(sample_names))]
                    out.write('\t'.join(row) + '\n')

            self.logger.info(f"写入 {len(kmers)} 个k-mer x {len(sample_names)} 个样本|"
                           f"Wrote {len(kmers)} k-mers x {len(sample_names)} samples")
