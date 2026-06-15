"""
K-mer库构建Pipeline|K-mer Database Build Pipeline
"""

import subprocess
import sys
from pathlib import Path
from typing import Dict

from .config import BuildConfig
from .utils import generate_fof_file, run_command, build_conda_command


class KmerBuildPipeline:
    """K-mer库构建Pipeline|K-mer Database Build Pipeline"""

    def __init__(self, config: BuildConfig, logger):
        """
        初始化Pipeline|Initialize Pipeline

        Args:
            config: 构建配置|Build configuration
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def run(self) -> Dict:
        """
        运行完整的建库流程|Run complete build pipeline

        Returns:
            dict: 构建结果统计|Build result statistics
        """
        self.logger.info("=" * 60)
        self.logger.info("开始K-mer库构建|Starting K-mer database build")
        self.logger.info("=" * 60)

        if self.config.use_kmtricks:
            return self._run_kmtricks_pipeline_full()
        else:
            return self._run_kmindex_pipeline()

    def _run_kmtricks_pipeline_full(self) -> Dict:
        """
        运行完整的kmtricks+RocksDB流程|Run complete kmtricks+RocksDB pipeline

        使用PipelineProcessor来处理，包含：
        - 样品命名规范化（处理特殊字符如①②）
        - 单端和双端fastq支持
        - 自动计算nb_partitions

        Returns:
            dict: 构建结果统计|Build result statistics
        """
        try:
            from .pipeline_processor import PipelineProcessor
            from .utils import CommandRunner

            # 创建命令运行器|Create command runner
            cmd_runner = CommandRunner(self.logger)

            # 创建PipelineProcessor|Create PipelineProcessor
            processor = PipelineProcessor(self.config, self.logger, cmd_runner)

            # 运行完整流程|Run full pipeline
            success = processor.run_full_pipeline()

            if not success:
                raise Exception("PipelineProcessor执行失败|PipelineProcessor execution failed")

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("K-mer库构建完成 (kmtricks+RocksDB模式)|K-mer database build completed (kmtricks+RocksDB mode)")
            self.logger.info("=" * 60)
            self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            return {
                'success': True,
                'mode': 'kmtricks',
                'input_dir': str(self.config.input_path),
                'output_dir': str(self.config.output_path)
            }

        except Exception as e:
            self.logger.error(f"建库失败|Build failed: {e}")
            return {
                'success': False,
                'error': str(e)
            }

    def _run_kmindex_pipeline(self) -> Dict:
        """
        运行kmindex流程|Run kmindex pipeline

        Returns:
            dict: 构建结果统计|Build result statistics
        """
        try:
            # 步骤1: 生成FOF文件|Step 1: Generate FOF file
            self.logger.info("[步骤1/2] 生成FOF文件|[Step 1/2] Generating FOF file")
            fof_file = self.config.output_path / "fof.txt"
            if not generate_fof_file(
                self.config.input_path,
                fof_file,
                self.config.fof_suffix_1,
                self.config.fof_suffix_2,
                self.logger
            ):
                raise Exception("生成FOF文件失败|Failed to generate FOF file")

            # 步骤2: 运行kmindex build|Step 2: Run kmindex build
            self.logger.info("[步骤2/2] 运行kmindex build|[Step 2/2] Running kmindex build")
            index_dir = self.config.output_path / "kmindex_index"
            index_file = self.config.output_path / self.config.index_name

            if not self._run_kmindex_build(fof_file, index_dir, index_file):
                raise Exception("kmindex build运行失败|kmindex build failed")

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("K-mer库构建完成 (kmindex模式)|K-mer database build completed (kmindex mode)")
            self.logger.info("=" * 60)
            self.logger.info(f"FOF文件|FOF file: {fof_file}")
            self.logger.info(f"kmindex索引目录|kmindex index directory: {index_dir}")
            self.logger.info(f"kmindex全局索引|kmindex global index: {index_file}")

            return {
                'success': True,
                'mode': 'kmindex',
                'fof_file': str(fof_file),
                'index_dir': str(index_dir),
                'index_file': str(index_file)
            }

        except Exception as e:
            self.logger.error(f"建库失败|Build failed: {e}")
            return {
                'success': False,
                'error': str(e)
            }

    def _run_kmtricks_pipeline(self, fof_file: Path, run_dir: Path) -> bool:
        """
        运行kmtricks pipeline|Run kmtricks pipeline

        Args:
            fof_file: FOF文件路径|FOF file path
            run_dir: 运行目录|Run directory

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 删除已存在的运行目录|Remove existing run directory
            if run_dir.exists():
                self.logger.info(f"删除已存在的运行目录|Removing existing run directory: {run_dir}")
                import shutil
                shutil.rmtree(run_dir)

            # 构建命令|Build command
            # 减少线程数以避免并发问题（181个样本时使用较少线程）
            threads = min(self.config.threads, 32)

            # 构建参数列表|Build argument list
            args = [
                'pipeline',
                '--file', str(fof_file),
                '--run-dir', str(run_dir),
                '--kmer-size', str(self.config.kmer_size),
                '--hard-min', str(self.config.hard_min),
                '--recurrence-min', str(self.config.recurrence_min),
                '--mode', self.config.mode,
                '--minimizer-size', str(self.config.minimizer_size),
                '--nb-partitions', '1000',  # 限制分区数以减少并发文件创建
                '-t', str(threads),
                '--cpr'  # 启用压缩|Enable compression
            ]

            # 使用conda命令包装器|Use conda command wrapper
            cmd = build_conda_command(self.config.kmtricks_path, args)

            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            result = run_command(cmd, self.logger, check=True)

            self.logger.info("kmtricks pipeline完成|kmtricks pipeline completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"kmtricks pipeline失败|kmtricks pipeline failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行kmtricks pipeline时出错|Error running kmtricks pipeline: {e}")
            return False

    def _run_kmtricks_aggregate(self, run_dir: Path, output_file: Path) -> bool:
        """
        运行kmtricks aggregate|Run kmtricks aggregate

        Args:
            run_dir: 运行目录|Run directory
            output_file: 输出文件|Output file

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 解析mode|Parse mode
            mode_parts = self.config.mode.split(':')
            matrix_type = mode_parts[0] if len(mode_parts) > 0 else 'kmer'

            # 构建参数列表|Build argument list
            args = [
                'aggregate',
                '--run-dir', str(run_dir),
                '--cpr-in',  # 输入文件是压缩的|Input files are compressed
                '--pa-matrix', matrix_type,
                '--format', 'text',
                '-t', str(self.config.threads),
                '--output', str(output_file)
            ]

            # 使用conda命令包装器|Use conda command wrapper
            cmd = build_conda_command(self.config.kmtricks_path, args)

            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            result = run_command(cmd, self.logger, check=True)

            self.logger.info("kmtricks aggregate完成|kmtricks aggregate completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"kmtricks aggregate失败|kmtricks aggregate failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行kmtricks aggregate时出错|Error running kmtricks aggregate: {e}")
            return False

    def _compress_matrix(self, input_file: Path, output_file: Path) -> bool:
        """
        使用bgzip压缩矩阵文件|Compress matrix file using bgzip

        Args:
            input_file: 输入文件|Input file
            output_file: 输出文件|Output file

        Returns:
            bool: 是否成功|Success
        """
        try:
            # 删除已存在的输出文件|Remove existing output file
            if output_file.exists():
                output_file.unlink()

            # 构建参数列表|Build argument list
            args = [
                '-c',  # 写入stdout|Write to stdout
                '-k',  # 保留输入文件|Keep input file
                '-@', str(self.config.threads),
                str(input_file)
            ]

            # 使用conda命令包装器|Use conda command wrapper
            cmd = build_conda_command(self.config.bgzip_path, args)

            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            # 运行命令并保存输出|Run command and save output
            result = run_command(cmd, self.logger, check=True)

            with open(output_file, 'w') as f:
                f.write(result.stdout)

            self.logger.info(f"矩阵文件压缩完成|Matrix file compressed: {output_file}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"bgzip压缩失败|bgzip compression failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"压缩矩阵文件时出错|Error compressing matrix file: {e}")
            return False

    def _extract_header(self, fof_file: Path) -> Path:
        """
        从FOF文件提取样本名作为header|Extract sample names from FOF file as header

        Args:
            fof_file: FOF文件路径|FOF file path

        Returns:
            Path: header文件路径|Header file path
        """
        header_file = fof_file.parent / "header.txt"

        try:
            sample_names = []
            with open(fof_file, 'r') as f:
                for line in f:
                    if line.strip():
                        # FOF格式: sample\tR1\tR2
                        sample = line.split('\t')[0].strip()
                        sample_names.append(sample)

            with open(header_file, 'w') as f:
                for sample in sample_names:
                    f.write(f"{sample}\n")

            self.logger.info(f"Header文件生成成功|Header file generated: {header_file}")
            self.logger.info(f"共|Total {len(sample_names)} 个样本|samples")

            return header_file

        except Exception as e:
            self.logger.error(f"生成header文件失败|Failed to generate header file: {e}")
            return None

    def _import_rocksdb(self, matrix_gz: Path, rocksdb_dir: Path,
                        header_file: Path = None) -> bool:
        """
        导入数据到RocksDB|Import data to RocksDB

        Args:
            matrix_gz: 压缩的矩阵文件|Compressed matrix file
            rocksdb_dir: RocksDB输出目录|RocksDB output directory
            header_file: Header文件路径|Header file path

        Returns:
            bool: 是否成功|Success
        """
        try:
            from .rocksdb.importer import RocksDBImporter
        except ImportError:
            self.logger.error("未找到RocksDB导入模块|RocksDB importer module not found")
            return False

        try:
            # 创建导入器|Create importer
            importer = RocksDBImporter(
                matrix_file=str(matrix_gz),
                db_path=str(rocksdb_dir),
                input_delimiter='\t',
                batch_size=self.config.batch_size,
                bloom_bits=self.config.bloom_bits,
                force_overwrite=True,
                header_file=str(header_file) if header_file else None,
                header_db_key="kmer_header",
                logger=self.logger
            )

            # 执行导入|Execute import
            success = importer.import_to_rocksdb()

            if success:
                self.logger.info(f"RocksDB导入成功|RocksDB import completed: {rocksdb_dir}")
            else:
                self.logger.error("RocksDB导入失败|RocksDB import failed")

            return success

        except Exception as e:
            self.logger.error(f"导入RocksDB时出错|Error importing to RocksDB: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _run_kmindex_build(self, fof_file: Path, index_dir: Path, index_file: Path) -> bool:
        """
        运行kmindex build|Run kmindex build

        Args:
            fof_file: FOF文件路径|FOF file path
            index_dir: 索引输出目录|Index output directory
            index_file: 全局索引文件路径|Global index file path

        Returns:
            bool: 是否成功|Success
        """
        try:
            import shutil

            # 删除已存在的运行目录|Remove existing run directory
            if index_dir.exists():
                self.logger.info(f"删除已存在的索引目录|Removing existing index directory: {index_dir}")
                shutil.rmtree(index_dir)

            # 删除已存在的全局索引|Remove existing global index
            if index_file.exists():
                self.logger.info(f"删除已存在的全局索引|Removing existing global index: {index_file}")
                shutil.rmtree(index_file)

            # 构建参数列表|Build argument list
            args = [
                'build',
                '-f', str(fof_file),
                '-d', str(index_dir),
                '-i', str(index_file),
                '-r', self.config.index_name,
                '-k', str(self.config.kmer_size),
                '--hard-min', str(self.config.hard_min),
                '--bloom-size', str(self.config.bloom_size),
                '-t', str(self.config.threads),
                '--cpr'  # 启用压缩|Enable compression
            ]

            # 使用conda命令包装器|Use conda command wrapper
            cmd = build_conda_command(self.config.kmindex_path, args)

            self.logger.info(f"运行命令|Running command: {' '.join(cmd)}")

            result = run_command(cmd, self.logger, check=True)

            self.logger.info("kmindex build完成|kmindex build completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"kmindex build失败|kmindex build failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行kmindex build时出错|Error running kmindex build: {e}")
            return False
