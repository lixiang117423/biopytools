"""
主程序模块|Main Program Module
"""

import os
import sys
import argparse
import shutil
import pandas as pd
import subprocess
from typing import List
from pathlib import Path
from .utils import setup_logger, check_dependencies, CommandRunner
from .file_processor import FileProcessor
from .jellyfish_processor import JellyfishProcessor
from .data_processor import DataProcessor
from .window_analyzer import WindowAnalyzer
from .config import KmerCountConfig
from .progress_manager import ProgressManager

# 自定义异常类
class FileIntegrityError(Exception):
    """文件完整性错误"""
    pass

class KmerCountAnalyzer:
    """K-mer计数分析器|K-mer Count Analyzer"""
    
    def __init__(self, config):
        self.config = config
        self.logger = setup_logger(config.output_dir, config.verbose)
        self.cmd_runner = CommandRunner(self.logger)
        
        # 初始化处理器|Initialize processors
        self.file_processor = FileProcessor(config, self.logger, self.cmd_runner)
        self.jellyfish_processor = JellyfishProcessor(config, self.logger, self.cmd_runner)
        self.data_processor = DataProcessor(config, self.logger)
        self.window_analyzer = WindowAnalyzer(config, self.logger)




    def _merge_sample_files(self, sample_files: List[Path]) -> pd.DataFrame:
        """合并样品文件（SQLite内存优化版本）

        优化策略：
        1. 使用SQLite数据库代替内存DataFrame
        2. 所有BED的k-mer初始化为0
        3. 逐样本更新检测到的k-mer
        4. 分块导出为TSV，避免内存爆炸

        内存占用：<1GB（恒定，不受数据规模影响）
        """
        import sqlite3

        self.logger.info("合并所有样本结果|Merge all sample results（SQLite优化版本|SQLite optimized version）")

        # SQLite数据库路径
        db_path = self.config.output_dir / "abundance_matrix_temp.db"
        output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"

        # ========== 步骤1：初始化数据库 ==========
        if not db_path.exists():
            self.logger.info("初始化SQLite数据库|Initializing SQLite database")

            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()

            # 性能优化
            cursor.execute('PRAGMA journal_mode=OFF')
            cursor.execute('PRAGMA synchronous=OFF')
            cursor.execute('PRAGMA cache_size=-500000')  # 500MB缓存

            # 创建表结构
            if self.config.bed_file:
                # 有BED文件：表结构包含位置信息
                cursor.execute('''
                    CREATE TABLE abundance_matrix (
                        unique_key TEXT PRIMARY KEY,
                        chr TEXT,
                        start INTEGER,
                        end INTEGER,
                        kmer TEXT,
                        score TEXT,
                        strand TEXT
                    )
                ''')

                # 获取样本名列表
                sample_names = [f.stem.replace('_kmer_abundance', '') for f in sample_files]

                # 为每个样本添加列（使用ALTER TABLE）
                for sample_name in sample_names:
                    try:
                        cursor.execute(f'ALTER TABLE abundance_matrix ADD COLUMN "{sample_name}" INTEGER DEFAULT 0')
                    except:
                        pass  # 列可能已存在

                # 批量插入BED数据（只插入正向序列）
                self.logger.info("插入BED数据到数据库|Inserting BED data into database")

                base_bed = self.data_processor.bed_info[self.data_processor.bed_info['strand'] == '+']
                bed_size = len(base_bed)

                # 分批插入，每批10万条
                batch_size = 100000
                for i in range(0, bed_size, batch_size):
                    chunk = base_bed.iloc[i:i+batch_size]

                    # 生成unique_key
                    chunk_data = []
                    for _, row in chunk.iterrows():
                        unique_key = f"{row['chr']}_{row['start']}_{row['end']}_{row['kmer']}"
                        chunk_data.append((
                            unique_key,
                            row['chr'],
                            int(row['start']),
                            int(row['end']),
                            row['kmer'],
                            row['score'],
                            row['strand']
                        ))

                    cursor.executemany('''
                        INSERT INTO abundance_matrix (unique_key, chr, start, end, kmer, score, strand)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', chunk_data)

                    if (i + batch_size) % 1000000 == 0 or i + batch_size >= bed_size:
                        conn.commit()
                        self.logger.info(f"已插入|Inserted {min(i + batch_size, bed_size):,} / {bed_size:,} 条BED记录|BED records")

                self.logger.info(f"BED数据插入完成|BED data insertion complete，共 {bed_size:,} 条记录|total records")

            else:
                # 无BED文件的情况
                raise NotImplementedError("无BED文件的SQLite合并暂未实现")

            conn.close()

        else:
            self.logger.info("使用已有数据库|Using existing database")

        # ========== 步骤2：逐样本更新丰度 ==========
        self.logger.info(f"更新样本丰度数据|Updating sample abundance data")

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        for idx, sample_file in enumerate(sample_files, 1):
            sample_name = sample_file.stem.replace('_kmer_abundance', '')

            try:
                # 检查文件是否为空
                if sample_file.stat().st_size < 10:
                    self.logger.warning(f"样本文件为空|Sample file is empty: {sample_name}")
                    continue

                # 读取样本文件
                sample_df = pd.read_csv(sample_file, sep='\t')

                if len(sample_df) == 0:
                    self.logger.warning(f"样本文件无数据|Sample file has no data: {sample_name}")
                    continue

                if 'unique_key' not in sample_df.columns:
                    self.logger.warning(f"样本文件缺少unique_key列|Sample file missing unique_key column: {sample_name}")
                    continue

                # 批量更新数据库
                self.logger.info(f"[{idx}/{len(sample_files)}] 更新样本|Updating sample: {sample_name} ({len(sample_df):,} 条记录|records)")

                # 分批更新，每批1万条
                update_batch_size = 10000
                for i in range(0, len(sample_df), update_batch_size):
                    chunk = sample_df.iloc[i:i+update_batch_size]

                    # 构建批量更新语句
                    for _, row in chunk.iterrows():
                        try:
                            cursor.execute(f'''
                                UPDATE abundance_matrix
                                SET "{sample_name}" = {row[sample_name]}
                                WHERE unique_key = ?
                            ''', (row['unique_key'],))
                        except:
                            pass  # k-mer可能不在BED中

                    if i % 50000 == 0:
                        conn.commit()

                conn.commit()

            except Exception as e:
                self.logger.error(f"处理样本失败|Sample processing failed: {sample_name}, 错误|error: {e}")

        conn.close()

        # ========== 步骤3：分块导出为TSV ==========
        self.logger.info("导出最终矩阵为TSV文件|Exporting final matrix to TSV")

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 获取总行数
        cursor.execute('SELECT COUNT(*) FROM abundance_matrix')
        total_rows = cursor.fetchone()[0]

        # 获取所有列名
        cursor.execute('PRAGMA table_info(abundance_matrix)')
        columns_info = cursor.fetchall()
        columns = [col[1] for col in columns_info]

        self.logger.info(f"导出|Exporting \1 行|rows × \2 列|columns")

        # 分块导出，每块10万行
        chunk_size = 100000
        offset = 0

        with open(output_file, 'w') as f:
            # 写入表头
            f.write('\t'.join(columns) + '\n')

            while offset < total_rows:
                # 为列名添加双引号（处理带特殊字符的列名如OV8-105）
                quoted_columns = [f'"{col}"' for col in columns]
                cursor.execute(f'''
                    SELECT {','.join(quoted_columns)}
                    FROM abundance_matrix
                    ORDER BY chr, start
                    LIMIT {chunk_size} OFFSET {offset}
                ''')

                chunk_rows = cursor.fetchall()

                # 写入数据
                for row in chunk_rows:
                    f.write('\t'.join(map(str, row)) + '\n')

                offset += chunk_size

                if offset % 1000000 == 0 or offset >= total_rows:
                    self.logger.info(f"已导出|Exported {min(offset, total_rows):,} / {total_rows:,} 行|rows")

        conn.close()

        # 删除临时数据库
        if db_path.exists():
            db_path.unlink()
            self.logger.info("临时数据库已删除|Temporary database deleted")

        self.logger.info(f"最终矩阵已保存|Final matrix saved: {output_file}")

        # 返回None（后续步骤会从文件读取）
        return None


    def run_analysis(self):
        """运行完整分析（支持断点续传）|Run complete analysis (with checkpoint resumption)"""
        try:
            self.logger.info("="*60)
            self.logger.info("开始K-mer丰度分析|Starting K-mer abundance analysis")
            self.logger.info("="*60)

            # 创建进度文件路径|Create progress file path
            progress_file = self.config.output_dir / ".kmer_count_progress.json"
            progress_mgr = ProgressManager(progress_file, self.logger)

            # 显示进度摘要|Display progress summary
            progress_mgr.print_progress_summary()
            
            # ========== 步骤1: 检查依赖软件 ==========
            if not progress_mgr.is_step_completed('check_dependencies'):
                progress_mgr.mark_step_in_progress('check_dependencies')
                self.logger.info("步骤1: 检查依赖软件|Step 1: Checking dependencies")
                if not check_dependencies(self.config.jellyfish_path, self.logger):
                    raise RuntimeError("依赖检查失败|Dependency check failed")
                progress_mgr.mark_step_completed('check_dependencies')
            else:
                self.logger.info("步骤1已完成，跳过|Step 1 completed, skipping")

            # ========== 步骤2: 设置临时目录 ==========
            # 临时目录每次运行都需要重新创建（不使用断点续传）
            # 因为临时文件（如generator文件）只在本次运行中有效
            progress_mgr.mark_step_in_progress('setup_temp_dir')
            self.logger.info("步骤2: 设置临时目录|Step 2: Setting up temporary directory")
            self.config.setup_temp_dir()
            self.logger.info(f"临时目录|Temporary directory: {self.config.temp_dir}")
            progress_mgr.mark_step_completed('setup_temp_dir')

            # ========== 步骤3: 查找输入文件 ==========
            if not progress_mgr.is_step_completed('find_input_files'):
                progress_mgr.mark_step_in_progress('find_input_files')
                self.logger.info("步骤3: 查找输入文件|Step 3: Finding input files")
                samples = self.file_processor.find_fastq_files()
                progress_mgr.mark_step_completed('find_input_files',
                                                metadata={'sample_count': len(samples)})
            else:
                self.logger.info("步骤3已完成，跳过|Step 3 completed, skipping")
                # 重新查找样本（因为需要samples变量）
                samples = self.file_processor.find_fastq_files()
                self.logger.info(f"找到{len(samples)}个样本|Found {len(samples)} samples")
            
            # ========== 步骤4: 准备k-mer库 ==========
            # 检查kmer_lib文件是否存在
            kmer_lib = self.config.output_dir / "intermediate_files" / "kmers_expanded.fasta"
            if not progress_mgr.is_step_completed('prepare_kmer_lib') or not kmer_lib.exists():
                progress_mgr.mark_step_in_progress('prepare_kmer_lib')
                self.logger.info("步骤4: 准备k-mer库|Step 4: Preparing k-mer library")
                kmer_lib = self.file_processor.prepare_kmer_library()
                progress_mgr.mark_step_completed('prepare_kmer_lib')
            else:
                self.logger.info("步骤4已完成，跳过|Step 4 completed, skipping")
                self.logger.info(f"使用已处理的k-mer库|Using processed k-mer library: {kmer_lib}")

            # 更新配置：步骤5及后续步骤使用扩展后的k-mer库
            self.config.kmer_lib = kmer_lib

            # ========== 步骤5: 解析注释信息 ==========
            intermediate_dir = self.config.output_dir / "intermediate_files"

            if self.config.bed_file:
                # 有BED文件：只处理BED文件，不需要SQLite数据库
                processed_bed = intermediate_dir / "processed_kmers.bed"

                need_parse = (not progress_mgr.is_step_completed('parse_annotations') or
                             not processed_bed.exists())

                if need_parse:
                    progress_mgr.mark_step_in_progress('parse_annotations')
                    self.logger.info("步骤5: 解析注释信息|Step 5: Parsing annotation information")
                    self.logger.info("检测到BED文件，跳过SQLite数据库构建|BED file detected, skipping SQLite database construction")
                    self.data_processor.parse_bed_file()
                    progress_mgr.mark_step_completed('parse_annotations')
                else:
                    self.logger.info("步骤5已完成，加载已有数据|Step 5 completed, loading existing data")
                    self.data_processor.load_existing_bed_file(processed_bed)
            else:
                # 无BED文件：需要SQLite数据库
                kmer_db = intermediate_dir / "kmer_mapping.db"

                need_parse = (not progress_mgr.is_step_completed('parse_annotations') or
                             not kmer_db.exists())

                if need_parse:
                    progress_mgr.mark_step_in_progress('parse_annotations')
                    self.logger.info("步骤5: 解析注释信息|Step 5: Parsing annotation information")
                    self.data_processor.parse_kmer_library()
                    progress_mgr.mark_step_completed('parse_annotations')
                else:
                    self.logger.info("步骤5已完成，加载已有数据|Step 5 completed, loading existing data")
                    self.data_processor.load_existing_database(kmer_db)
            
            # ========== 步骤6: 处理样本 ==========
            sample_files = []
            processed_samples = 0
            skipped_samples = 0
            
            # 检查是否已完成样本处理
            if progress_mgr.is_step_completed('process_samples'):
                self.logger.info("步骤6已完成，加载已有结果|Step 6 completed, loading existing results")
                # 加载已有的样本文件
                each_sample_dir = self.config.output_dir / "each_sample"
                if each_sample_dir.exists():
                    sample_files = list(each_sample_dir.glob("*_kmer_abundance.tsv"))
                    processed_samples = len(sample_files)
                    self.logger.info(f"找到{processed_samples}个已有的样本结果文件|Found {processed_samples} existing sample result files")
            else:
                progress_mgr.mark_step_in_progress('process_samples')
                self.logger.info("步骤6: 处理样本|Step 6: Processing samples")
                
                # 将file_processor引用传给jellyfish_processor
                self.jellyfish_processor.file_processor = self.file_processor

                total_samples = len(samples)

                for sample_name, r1_file, r2_file in samples:
                    try:
                        self.logger.info(f"处理样本|Processing sample: {sample_name}")
                        
                        # 准备文件
                        fastq_files = [r1_file, r2_file]
                        
                        # k-mer计数
                        jf_file = self.jellyfish_processor.count_kmers(sample_name, fastq_files)
                        
                        # 查询k-mer丰度
                        count_file = self.jellyfish_processor.query_kmers(sample_name, jf_file, kmer_lib)
                        
                        # 解析结果
                        sample_df = self.data_processor.parse_jellyfish_output(count_file, sample_name)
                        
                        # 与BED合并并保存
                        merged_sample_df = self.data_processor.merge_single_sample_with_bed(sample_df, sample_name)

                        # 创建each_sample子目录
                        each_sample_dir = self.config.output_dir / "each_sample"
                        each_sample_dir.mkdir(exist_ok=True)

                        sample_file = each_sample_dir / f"{sample_name}_kmer_abundance.tsv"
                        merged_sample_df.to_csv(sample_file, sep='	', index=False)
                        sample_files.append(sample_file)
                        processed_samples += 1
                        self.logger.info(f"样品{sample_name}结果已保存|Sample {sample_name} results saved: {sample_file}")
                        
                        # 清理临时文件
                        for f in fastq_files:
                            if f and os.path.exists(f):
                                if str(self.config.temp_dir) in f:
                                    os.remove(f)

                    except Exception as e:
                        skipped_samples += 1
                        error_msg = str(e)
                        
                        if ("killed by signal 13" in error_msg or 
                            "Some generator commands failed" in error_msg or
                            "SIGPIPE" in error_msg.upper()):
                            self.logger.error(f"样本 {sample_name} 遇到SIGPIPE错误|encountered SIGPIPE error")
                        elif "FileIntegrityError" in error_msg or "文件完整性" in error_msg:
                            self.logger.error(f"样本 {sample_name} 文件完整性检查失败|file integrity check failed")
                        else:
                            self.logger.error(f"样本 {sample_name} 处理失败|processing failed: {error_msg}")
                        
                        self.logger.warning(f"跳过该样本，继续处理其他样本|Skipping this sample, continuing with others")
                        continue

                # 输出处理总结
                self.logger.info("="*50)
                self.logger.info(f"样本处理总结|Sample Processing Summary:")
                self.logger.info(f"  总样本数|Total samples: \1
                self.logger.info(f"  成功处理|Successfully processed: {processed_samples}")
                self.logger.info(f"   跳过样本|Skipped samples: {skipped_samples}")
                if skipped_samples > 0:
                    self.logger.warning(f"   跳过比例|Skip ratio: {skipped_samples/total_samples:.1%}")
                self.logger.info("="*50)

                # 标记步骤完成
                progress_mgr.mark_step_completed('process_samples',
                                                metadata={
                                                    'total_samples': total_samples,
                                                    'processed_samples': processed_samples,
                                                    'skipped_samples': skipped_samples
                                                })

            # 检查是否有成功处理的样本
            if not sample_files:
                raise RuntimeError("没有样本成功处理，分析无法继续|No samples processed successfully, analysis cannot continue")

            # ========== 步骤7: 合并所有样本结果 ==========
            if not progress_mgr.is_step_completed('merge_results'):
                progress_mgr.mark_step_in_progress('merge_results')
                self.logger.info("步骤7: 合并所有样本结果|Step 7: Merging all sample results")
                final_df = self._merge_sample_files(sample_files)
                progress_mgr.mark_step_completed('merge_results')
            else:
                self.logger.info("步骤7已完成，跳过|Step 7 completed, skipping")
                # 需要重新加载或生成final_df用于后续步骤
                output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                if output_file.exists():
                    final_df = pd.read_csv(output_file, sep='	')
                    self.logger.info(f"加载已有结果矩阵|Loading existing result matrix: {output_file}")
                else:
                    # 如果结果文件不存在，重新合并
                    final_df = self._merge_sample_files(sample_files)

            # ========== 步骤8: 保存最终结果 ==========
            if not progress_mgr.is_step_completed('save_results'):
                progress_mgr.mark_step_in_progress('save_results')
                self.logger.info("步骤8: 保存最终结果|Step 8: Saving final results")

                # SQLite版本已经在步骤7保存了文件
                if final_df is None:
                    self.logger.info("结果矩阵已在步骤7保存|Result matrix already saved in step 7")
                    output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                else:
                    # 旧版本：需要保存DataFrame
                    output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                    final_df.to_csv(output_file, sep='\t', index=False)
                    self.logger.info(f"最终结果已保存|Final results saved: {output_file}")

                progress_mgr.mark_step_completed('save_results',
                                                metadata={'output_file': str(output_file)})
            else:
                self.logger.info("步骤8已完成，跳过|Step 8 completed, skipping")

            # ========== 步骤9: 滑动窗口分析（可选）==========
            if self.config.bed_file and self.config.window_size:
                if not progress_mgr.is_step_completed('sliding_window'):
                    progress_mgr.mark_step_in_progress('sliding_window')
                    self.logger.info("步骤9: 滑动窗口分析|Step 9: Sliding window analysis")

                    # 如果内存中没有DataFrame，从文件读取
                    if final_df is None:
                        output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                        self.logger.info(f"从文件加载数据用于滑动窗口分析|Loading data from file for sliding window analysis: {output_file}")
                        final_df = pd.read_csv(output_file, sep='\t')

                    window_df = self.window_analyzer.sliding_window_analysis(final_df)
                    if window_df is not None:
                        window_file = self.config.output_dir / "sliding_window_analysis.tsv"
                        window_df.to_csv(window_file, sep='\t', index=False)
                        self.logger.info(f"滑动窗口结果已保存|Sliding window results saved: {window_file}")
                    progress_mgr.mark_step_completed('sliding_window')
                else:
                    self.logger.info("步骤9已完成，跳过|Step 9 completed, skipping")
            else:
                self.logger.info("步骤9跳过（未指定窗口大小或无BED文件）|Step 9 skipped (no window size or no BED file)")

            # ========== 步骤10: 创建0/1矩阵（默认生成）==========
            if not progress_mgr.is_step_completed('binary_matrix'):
                progress_mgr.mark_step_in_progress('binary_matrix')
                self.logger.info("步骤10: 创建0/1矩阵|Step 10: Creating 0/1 matrix")

                # 如果内存中没有DataFrame，从文件读取
                if final_df is None:
                    output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                    self.logger.info(f"从文件加载数据用于0/1矩阵创建|Loading data from file for binary matrix creation: {output_file}")
                    final_df = pd.read_csv(output_file, sep='\t')

                binary_df = self.window_analyzer.create_binary_matrix(final_df)
                binary_file = self.config.output_dir / "kmer_binary_matrix.tsv"
                binary_df.to_csv(binary_file, sep='\t', index=False)
                self.logger.info(f"0/1矩阵已保存|0/1 matrix saved: {binary_file}")
                progress_mgr.mark_step_completed('binary_matrix')
            else:
                self.logger.info("步骤10已完成，跳过|Step 10 completed, skipping")

            # ========== 生成统计报告 ==========
            # 如果内存中没有DataFrame，从文件读取
            if final_df is None:
                output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
                if output_file.exists():
                    self.logger.info(f"从文件加载数据用于统计报告|Loading data from file for summary report: {output_file}")
                    final_df = pd.read_csv(output_file, sep='\t')
                else:
                    self.logger.warning("无法找到结果文件，跳过统计报告|Cannot find result file, skipping summary report")
                    final_df = None

            if final_df is not None:
                self._generate_summary_report(final_df, processed_samples, skipped_samples)
            
            # 所有步骤完成，删除进度文件
            progress_mgr.reset_progress()
            
            self.logger.info("="*60)
            if skipped_samples > 0:
                self.logger.warning(f"K-mer丰度分析完成，但有{skipped_samples}个样本被跳过|K-mer abundance analysis completed, but {skipped_samples} samples were skipped")
            else:
                self.logger.info("K-mer丰度分析完成|K-mer abundance analysis completed!")
            self.logger.info("="*60)
            
        except Exception as e:
            self.logger.error(f"分析过程中出现错误|Error during analysis: {e}")
            # 保存错误信息到进度文件
            if 'progress_mgr' in locals():
                progress_mgr.progress_data['last_error'] = str(e)
                progress_mgr._save_progress()
            raise
        finally:
            # 清理资源
            if hasattr(self, 'data_processor'):
                self.data_processor.cleanup()
            
            # 清理临时目录
            if self.config.temp_dir and not self.config.keep_temp:
                try:
                    shutil.rmtree(self.config.temp_dir)
                    self.logger.info("临时目录已清理|Temporary directory cleaned")
                except Exception as e:
                    self.logger.warning(f"清理临时目录失败|Failed to clean temporary directory: {e}")


    def _generate_summary_report(self, final_df: pd.DataFrame, processed_samples: int = 0, skipped_samples: int = 0):
        """生成统计摘要报告|Generate summary statistics report"""
        self.logger.info("生成统计摘要报告|Generating summary statistics report")
        
        # 获取样本列|Get sample columns
        sample_columns = [col for col in final_df.columns
                         if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id', 'score', 'strand', 'unique_key']
                         and not col.startswith('tmp_')]
        
        report_lines = []
        report_lines.append("K-mer丰度分析统计报告|K-mer Abundance Analysis Statistics Report")
        report_lines.append("="*80)
        report_lines.append("")
        
        # 处理统计|Processing statistics
        if processed_samples > 0 or skipped_samples > 0:
            total_samples = processed_samples + skipped_samples
            report_lines.append("处理统计信息|Processing Statistics:")
            report_lines.append(f"  总样本数|Total samples: {total_samples}")
            report_lines.append(f"  成功处理|Successfully processed: {processed_samples}")
            report_lines.append(f"   跳过样本|Skipped samples: {skipped_samples}")
            if total_samples > 0:
                report_lines.append(f"  成功率|Success rate: {processed_samples/total_samples:.1%}")
            report_lines.append("")
        
        # 基本统计|Basic statistics
        report_lines.append("基本统计信息|Basic Statistics:")
        report_lines.append(f"  总k-mer数量|Total k-mers: {len(final_df):,}")
        report_lines.append(f"  最终样本数量|Final sample count: {len(sample_columns)}")
        report_lines.append(f"  K-mer长度|K-mer size: {self.config.kmer_size}")
        report_lines.append("")

        # 样本统计|Sample statistics
        if sample_columns:
            report_lines.append("样本统计信息|Sample Statistics:")
            for sample in sample_columns:
                try:
                    # 添加这一行：确保列是数值类型
                    sample_data = pd.to_numeric(final_df[sample], errors='coerce').fillna(0)
                    
                    total_abundance = sample_data.sum()
                    present_kmers = (sample_data > 0).sum()
                    avg_abundance = sample_data.mean()
                    
                    report_lines.append(f"  {sample}:")
                    report_lines.append(f"    总丰度|Total abundance: {total_abundance:,}")
                    report_lines.append(f"    存在k-mer数|Present k-mers: {present_kmers:,}")
                    report_lines.append(f"    平均丰度|Average abundance: {avg_abundance:.2f}")
                    report_lines.append(f"    存在比例|Presence ratio: {present_kmers/len(final_df):.2%}")
                    
                except Exception as e:
                    self.logger.warning(f" 样本 {sample} 统计计算失败|statistics calculation failed: {e}")
                    report_lines.append(f"  {sample}: 统计计算失败")
            
            report_lines.append("")
        
        # 保存报告|Save report
        report_file = self.config.output_dir / "analysis_summary.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_lines))
        
        self.logger.info(f"统计报告已保存|Summary report saved: {report_file}")

def create_parser():
    """创建命令行参数解析器|Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='K-mer丰度分析工具|K-mer abundance analysis tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i /data/fastq -p "*_1.fq.gz" -k kmers.fasta -o results/
        """
    )
    
    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input', required=True,
                        help='输入文件目录|Input files directory')
    required.add_argument('-p', '--pattern', required=True,
                        help='文件模式，支持FASTQ和FASTA格式，如*_1.fq.gz、*.fa|File pattern, support FASTQ and FASTA formats, e.g. *_1.fq.gz, *.fa')
    required.add_argument('-k', '--kmer-lib', required=True,
                        help='K-mer库文件(FASTA格式)|K-mer library file (FASTA format)')
    required.add_argument('-o', '--output', required=True,
                        help='输出目录|Output directory')
    
    # 可选参数|Optional arguments
    optional = parser.add_argument_group(' 可选参数|Optional arguments')
    optional.add_argument('-b', '--bed-file',
                        help='BED文件路径|BED file path')
    optional.add_argument('-m', '--kmer-size', type=int, default=51,
                        help='K-mer长度 |K-mer size (default: %(default)s)')
    optional.add_argument('-s', '--hash-size', default='1000M',
                        help=' 哈希表大小 |Hash table size (default: %(default)s)')
    optional.add_argument('-t', '--threads', type=int, default=8,
                        help='线程数 |Number of threads (default: %(default)s)')
    optional.add_argument('-w', '--window-size', type=int, default=500000,
                        help='滑动窗口大小bp |Sliding window size in bp (default: %(default)s)')
    optional.add_argument('--step-size', type=int,
                        help='滑动窗口步长bp (默认: window-size/5)|Sliding window step size in bp (default: window-size/5)')
    optional.add_argument('-C', '--canonical', action='store_true',
                        help='统计正向和反向互补链|Count both forward and reverse complement')
    optional.add_argument('--keep-temp', action='store_true',
                        help='保留临时文件|Keep temporary files')
    optional.add_argument('--keep-binary', action='store_true',
                        help='保留0/1存在缺失矩阵|Keep 0/1 presence/absence matrix')
    optional.add_argument('--jellyfish-path', default='jellyfish',
                        help='Jellyfish程序路径 |Jellyfish program path (default: %(default)s)')
    optional.add_argument('-v', '--verbose', action='store_true',
                        help='详细输出|Verbose output')
    
    return parser


def main():
    """主函数 - 这是run_kmer_count脚本的入口点|Main function - Entry point for run_kmer_count script"""
    parser = create_parser()
    args = parser.parse_args()

    try:
        # 创建配置|Create configuration
        config = KmerCountConfig.from_args(args)

        # 验证配置|Validate configuration
        config.validate()

        # 创建分析器并运行|Create analyzer and run
        analyzer = KmerCountAnalyzer(config)
        analyzer.run_analysis()

    except KeyboardInterrupt:
        print("分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()