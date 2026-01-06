"""
BLAST分析主程序模块|BLAST Analysis Main Module
标准化BLAST分析器，遵循开发规范
"""

import os
import sys
import subprocess
import glob
import re
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from ..core.base_analyzer import BaseAnalyzer
from .config import BLASTConfig


class BLASTAnalyzer(BaseAnalyzer):
    """BLAST分析主类|Main BLAST Analysis Class"""

    def __init__(self, config: BLASTConfig):
        """
        初始化BLAST分析器|Initialize BLAST analyzer

        Args:
            config: BLAST配置对象|BLAST configuration object
        """
        super().__init__(config)
        self.sample_mapping = {}
        self.database_path = None

    def validate_inputs(self) -> bool:
        """
        验证输入|Validate inputs

        Returns:
            bool: 验证是否通过|Whether validation passed
        """
        try:
            # 调用父类验证
            super().validate_inputs()

            # 检查依赖工具|Check dependencies
            self._check_dependencies()

            # 生成或验证样品映射|Generate or validate sample mapping
            self._generate_sample_mapping()

            return True

        except Exception as e:
            self.logger.error(f"输入验证失败|Input validation failed: {e}")
            return False

    def _check_dependencies(self):
        """检查BLAST依赖工具|Check BLAST dependencies"""
        # 获取正确的BLAST程序路径|Get correct BLAST program path
        blast_program = getattr(self.config, f'{self.config.blast_type}_path', self.config.blast_type)

        required_tools = [self.config.makeblastdb_path, blast_program]

        for tool in required_tools:
            if not shutil.which(tool):
                raise RuntimeError(f"未找到依赖工具|Required tool not found: {tool}")

    def _generate_sample_mapping(self):
        """生成或验证样品映射|Generate or validate sample mapping"""
        if self.config.sample_map_file:
            # 使用现有的样品映射文件|Use existing sample mapping file
            self._load_sample_mapping()
        elif self.config.input:
            # 自动生成样品映射文件|Auto-generate sample mapping
            self._auto_generate_sample_mapping()
        else:
            raise ValueError("必须指定输入文件或样品映射文件|Must specify input file or sample mapping file")

    def _load_sample_mapping(self):
        """加载样品映射文件|Load sample mapping file"""
        self.logger.info(f"加载样品映射文件|Loading sample mapping file: {self.config.sample_map_file}")

        try:
            with open(self.config.sample_map_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) >= 2:
                        file_path = parts[0].strip()
                        sample_name = parts[1].strip()

                        if os.path.exists(file_path):
                            self.sample_mapping[sample_name] = file_path
                        else:
                            self.logger.warning(f"文件不存在|File not found: {file_path}")
                    else:
                        self.logger.warning(f"第{line_num}行格式错误|Line {line_num} format error: {line}")

            if not self.sample_mapping:
                raise ValueError("样品映射文件为空或格式错误|Sample mapping file is empty or incorrectly formatted")

            self.logger.info(f"加载了{len(self.sample_mapping)}个样品|Loaded {len(self.sample_mapping)} samples")

        except Exception as e:
            raise RuntimeError(f"加载样品映射文件失败|Failed to load sample mapping file: {e}")

    def _auto_generate_sample_mapping(self):
        """自动生成样品映射|Auto-generate sample mapping"""
        self.logger.info(f"自动扫描输入文件|Auto-scanning input files: {self.config.input}")

        input_path = Path(self.config.input)

        if input_path.is_file():
            # 单个文件|Single file
            sample_name = self.config.sample_name or input_path.stem
            self.sample_mapping[sample_name] = str(input_path)
            self.config.auto_generated_map = True

        elif input_path.is_dir():
            # 目录，扫描匹配的文件|Directory, scan matching files
            pattern = self.config.input_suffix
            files = glob.glob(os.path.join(str(input_path), pattern))

            for file_path in files:
                file_name = os.path.basename(file_path)
                # 使用正则表达式提取样品名称
                match = re.search(self.config.sample_name_pattern, file_name)
                if match:
                    sample_name = match.group(1)
                    self.sample_mapping[sample_name] = file_path

            self.config.auto_generated_map = True

        else:
            raise FileNotFoundError(f"输入路径不存在|Input path not found: {self.config.input}")

        if not self.sample_mapping:
            raise ValueError(f"在{self.config.input}中未找到匹配的文件|No matching files found in {self.config.input}")

        self.logger.info(f"发现了{len(self.sample_mapping)}个样品|Found {len(self.sample_mapping)} samples")

    def run_analysis(self) -> bool:
        """
        运行BLAST分析流程|Run BLAST analysis pipeline

        Returns:
            bool: 分析是否成功|Whether analysis was successful
        """
        try:
            # 步骤1: 创建输出目录|Step 1: Create output directory
            if not self._create_output_directory():
                return False

            # 步骤2: 创建BLAST数据库|Step 2: Create BLAST database
            if not self._create_blast_database():
                return False

            # 步骤3: 运行BLAST比对|Step 3: Run BLAST alignment
            if not self._run_blast_alignment():
                return False

            # 步骤4: 处理结果|Step 4: Process results
            if not self._process_results():
                return False

            # 步骤5: 生成统计报告|Step 5: Generate statistics report
            if not self._generate_statistics_report():
                return False

            # 步骤6: 生成比对可视化|Step 6: Generate alignment visualization
            if self.config.alignment_output != 'none':
                if not self._generate_alignment_visualization():
                    return False

            return True

        except Exception as e:
            self.logger.error(f"BLAST分析失败|BLAST analysis failed: {e}", exc_info=True)
            return False

    def _create_output_directory(self) -> bool:
        """创建输出目录|Create output directory"""
        try:
            os.makedirs(self.config.output, exist_ok=True)
            self.logger.info(f"输出目录已创建|Output directory created: {self.config.output}")
            return True
        except Exception as e:
            self.logger.error(f"创建输出目录失败|Failed to create output directory: {e}")
            return False

    def _create_blast_database(self) -> bool:
        """创建BLAST数据库|Create BLAST database"""
        try:
            self.log_step_start("创建BLAST数据库|Creating BLAST Database", 1)

            db_path = self.config.get_target_db_path()

            # 检查数据库是否已存在|Check if database already exists
            if os.path.exists(db_path) and not self.config.force:
                self.logger.info(f"数据库已存在，跳过创建|Database already exists, skipping creation: {db_path}")
                self.database_path = db_path
                self.log_step_end("创建BLAST数据库|Creating BLAST Database", True)
                return True

            # 构建数据库|Build database
            cmd = [
                self.config.makeblastdb_path,
                '-in', self.config.reference,
                '-dbtype', self.config.target_db_type,
                '-out', db_path
            ]

            self.logger.info(f"执行命令|Executing command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            self.database_path = db_path
            self.logger.info(f"数据库创建成功|Database created successfully: {db_path}")

            self.log_step_end("创建BLAST数据库|Creating BLAST Database", True)
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"创建数据库失败|Database creation failed: {e}")
            if e.stdout:
                self.logger.error(f"stdout: {e.stdout}")
            if e.stderr:
                self.logger.error(f"stderr: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"创建数据库时发生错误|Error creating database: {e}")
            return False

    def _run_blast_alignment(self) -> bool:
        """运行BLAST比对|Run BLAST alignment"""
        try:
            self.log_step_start("运行BLAST比对|Running BLAST Alignment", 2)

            all_results = []
            total_samples = len(self.sample_mapping)

            for i, (sample_name, sample_file) in enumerate(self.sample_mapping.items(), 1):
                self.logger.info(f"处理样品 {i}/{total_samples}: {sample_name}|Processing sample {i}/{total_samples}: {sample_name}")

                result_file = self._run_single_blast(sample_name, sample_file)
                if result_file:
                    all_results.append(result_file)
                else:
                    self.logger.warning(f"样品 {sample_name} 比对失败|Sample {sample_name} alignment failed")

            if not all_results:
                raise RuntimeError("所有样品比对都失败了|All sample alignments failed")

            # 合并结果|Merge results
            merged_file = self._merge_results(all_results)

            self.logger.info(f"BLAST比对完成|BLAST alignment completed. Results: {merged_file}")

            self.log_step_end("运行BLAST比对|Running BLAST Alignment", True)
            return True

        except Exception as e:
            self.logger.error(f"BLAST比对失败|BLAST alignment failed: {e}")
            self.log_step_end("运行BLAST比对|Running BLAST Alignment", False)
            return False

    def _run_single_blast(self, sample_name: str, sample_file: str) -> Optional[str]:
        """运行单个样品的BLAST比对|Run BLAST alignment for single sample"""
        try:
            output_file = os.path.join(self.config.output, f"{sample_name}_{self.config.blast_type}_results.tsv")

            # 获取正确的BLAST程序路径|Get correct BLAST program path
            blast_program = getattr(self.config, f'{self.config.blast_type}_path', self.config.blast_type)

            cmd = [
                blast_program,
                '-query', sample_file,
                '-db', self.database_path,
                '-out', output_file,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq',
                '-evalue', str(self.config.evalue),
                '-max_target_seqs', str(self.config.max_target_seqs),
                '-word_size', str(self.config.word_size),
                '-num_threads', str(self.config.threads)
            ]

            self.logger.debug(f"执行BLAST命令|Executing BLAST command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # 检查输出文件|Check output file
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                self.logger.info(f"样品 {sample_name} 比对完成|Sample {sample_name} alignment completed")
                return output_file
            else:
                self.logger.warning(f"样品 {sample_name} 比对输出为空|Sample {sample_name} alignment output is empty")
                return None

        except subprocess.CalledProcessError as e:
            self.logger.error(f"样品 {sample_name} BLAST比对失败|Sample {sample_name} BLAST alignment failed: {e}")
            if e.stderr:
                self.logger.error(f"stderr: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"样品 {sample_name} 比对时发生错误|Error during sample {sample_name} alignment: {e}")
            return None

    def _merge_results(self, result_files: List[str]) -> str:
        """合并所有结果文件并计算覆盖度|Merge all result files and calculate coverage"""
        try:
            merged_file = self.config.get_summary_output_path()

            with open(merged_file, 'w', encoding='utf-8') as out_f:
                # 写入中文表头|Write Chinese header (16列)
                out_f.write("样品名称\t查询序列ID\t目标序列ID\t序列相似度(%)\t比对长度\t错配数\tGap数\t")
                out_f.write("查询起始位置\t查询结束位置\t目标起始位置\t目标结束位置\tE-value\tBit_Score\t")
                out_f.write("目标序列长度\t目标序列覆盖度(%)\t查询序列\t目标序列\n")

                for result_file in result_files:
                    # Extract sample name from result file
                    # Format: {sample_name}_{blast_type}_results.tsv
                    basename = os.path.basename(result_file)
                    # Remove blast_type and suffix
                    parts = basename.replace('_results.tsv', '').rsplit('_', 1)
                    sample_name = parts[0] if len(parts) > 1 else basename

                    try:
                        with open(result_file, 'r', encoding='utf-8') as in_f:
                            for line in in_f:
                                line = line.strip()
                                if line:
                                    parts = line.split('\t')
                                    if len(parts) >= 14:
                                        try:
                                            # 计算覆盖度|Calculate coverage
                                            # BLAST output format (15 columns): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq
                                            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, slen, qseq, sseq = parts[:15]

                                            # 计算目标序列覆盖度|Calculate subject sequence coverage
                                            try:
                                                sstart_int = int(sstart)
                                                send_int = int(send)
                                                slen_int = int(slen)
                                                coverage = abs(send_int - sstart_int + 1) / slen_int * 100
                                                coverage = min(coverage, 100.0)
                                            except (ValueError, ZeroDivisionError):
                                                coverage = 0.0

                                            # 写入合并文件，添加覆盖度列|Write to merged file with coverage column
                                            out_line = '\t'.join([
                                                sample_name, qseqid, sseqid, pident, length, mismatch, gapopen,
                                                qstart, qend, sstart, send, evalue, bitscore, slen,
                                                f"{coverage:.2f}", qseq, sseq
                                            ])
                                            out_f.write(out_line + '\n')
                                        except (ValueError, IndexError) as e:
                                            self.logger.debug(f"跳过格式错误的行|Skipping malformed line: {e}")
                                            continue
                    except Exception as e:
                        self.logger.warning(f"读取结果文件失败|Failed to read result file: {result_file}, error: {e}")

            self.logger.info(f"合并结果文件|Merged results file: {merged_file}")
            return merged_file

        except Exception as e:
            raise RuntimeError(f"合并结果文件失败|Failed to merge result files: {e}")

    def _process_results(self) -> bool:
        """处理结果|Process results"""
        try:
            self.log_step_start("处理结果|Processing Results", 3)

            summary_file = self.config.get_summary_output_path()

            if not os.path.exists(summary_file):
                raise FileNotFoundError(f"汇总结果文件不存在|Summary result file not found: {summary_file}")

            # 统计结果|Generate statistics
            stats = self._generate_statistics(summary_file)
            self.log_statistics(stats)

            # 对结果排序|Sort results
            sorted_file = self._sort_results(summary_file)

            # 创建高质量结果文件|Create high quality results file
            self._create_high_quality_results(sorted_file)

            self.log_step_end("处理结果|Processing Results", True)
            return True

        except Exception as e:
            self.logger.error(f"处理结果失败|Failed to process results: {e}")
            self.log_step_end("处理结果|Processing Results", False)
            return False

    def _generate_statistics(self, result_file: str) -> Dict:
        """生成统计信息|Generate statistics"""
        try:
            stats = {
                'total_alignments': 0,
                'samples_count': len(self.sample_mapping),
                'unique_queries': set(),
                'unique_subjects': set()
            }

            with open(result_file, 'r') as f:
                next(f)  # 跳过表头

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 13:  # 13列: Sample + 12列BLAST输出
                        stats['total_alignments'] += 1
                        stats['unique_queries'].add(parts[1])  # qseqid
                        stats['unique_subjects'].add(parts[2])  # sseqid

            stats['unique_queries'] = len(stats['unique_queries'])
            stats['unique_subjects'] = len(stats['unique_subjects'])

            return stats

        except Exception as e:
            self.logger.error(f"生成统计信息失败|Failed to generate statistics: {e}")
            return {'total_alignments': 0, 'samples_count': 0, 'unique_queries': 0, 'unique_subjects': 0}

    def _sort_results(self, summary_file: str) -> str:
        """对结果按覆盖度、相似度和E-value排序|Sort results by coverage, identity and E-value"""
        try:
            self.logger.info("按覆盖度、相似度和E-value排序结果|Sorting results by coverage, identity and E-value")

            sorted_file = summary_file.replace('.tsv', '_sorted.tsv')

            # 读取所有数据行
            header = None
            data_lines = []

            with open(summary_file, 'r', encoding='utf-8') as f:
                header = f.readline().strip()
                for line in f:
                    line = line.strip()
                    if line:
                        data_lines.append(line)

            if not data_lines:
                self.logger.warning("没有数据需要排序|No data to sort")
                return summary_file

            # 解析并排序
            def parse_line(line):
                parts = line.split('\t')
                # 新格式：16列 - Sample qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen coverage qseq sseq
                if len(parts) < 16:
                    return None

                try:
                    # 解析列索引
                    # 0: Sample, 1: qseqid, 2: sseqid, 3: pident, 4: length,
                    # 5: mismatch, 6: gapopen, 7: qstart, 8: qend, 9: sstart, 10: send,
                    # 11: evalue, 12: bitscore, 13: slen, 14: coverage, 15: qseq, 16: sseq
                    identity = float(parts[3])
                    coverage = float(parts[14])
                    evalue_str = parts[11]

                    # Parse e-value
                    try:
                        if evalue_str.startswith('e') or evalue_str.startswith('E'):
                            evalue = float(f"1{evalue_str}")
                        else:
                            evalue = float(evalue_str)
                    except:
                        evalue = 1.0

                    return {
                        'line': line,
                        'coverage': coverage,
                        'identity': identity,
                        'evalue': evalue
                    }
                except (ValueError, IndexError):
                    return None

            # 解析所有行
            parsed_data = [parse_line(line) for line in data_lines]
            parsed_data = [d for d in parsed_data if d is not None]

            # 排序：覆盖度降序、相似度降序、E-value升序
            parsed_data.sort(key=lambda x: (-x['coverage'], -x['identity'], x['evalue']))

            # 写入排序后的文件
            with open(sorted_file, 'w', encoding='utf-8') as f:
                f.write(header + '\n')
                for item in parsed_data:
                    f.write(item['line'] + '\n')

            self.logger.info(f"排序完成|Sorting completed: {sorted_file}")
            self.logger.info(f"排序结果数|Sorted results: {len(parsed_data)}")

            return sorted_file

        except Exception as e:
            self.logger.error(f"排序失败|Sorting failed: {e}")
            return summary_file

    def _create_high_quality_results(self, sorted_file: str) -> bool:
        """从排序文件创建高质量比对结果文件|Create high quality alignment results from sorted file"""
        try:
            self.logger.info("从排序文件筛选高质量比对结果|Filtering high quality alignments from sorted file")

            # 从 sorted.tsv 生成 sorted_high_quality.tsv
            if '_sorted.tsv' in sorted_file:
                high_quality_file = sorted_file.replace('_sorted.tsv', '_sorted_high_quality.tsv')
            else:
                high_quality_file = sorted_file.replace('.tsv', '_high_quality.tsv')

            # 读取排序后的结果
            high_quality_alignments = []
            total_count = 0

            with open(sorted_file, 'r', encoding='utf-8') as f:
                # 读取表头
                header = f.readline().strip()

                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    total_count += 1
                    parts = line.split('\t')

                    # 新格式：16列
                    if len(parts) < 16:
                        continue

                    try:
                        # Sample, qseqid, sseqid, pident, length, mismatch, gapopen,
                        # qstart, qend, sstart, send, evalue, bitscore, slen, coverage, qseq, sseq
                        sample_name = parts[0]
                        query_id = parts[1]
                        subject_id = parts[2]
                        identity = float(parts[3])
                        length = int(parts[4])
                        evalue = parts[11]
                        coverage = float(parts[14])

                        # 解析e-value
                        try:
                            if evalue.startswith('e') or evalue.startswith('E'):
                                evalue_float = float(f"1{evalue}")
                            else:
                                evalue_float = float(evalue)
                        except:
                            evalue_float = 1.0

                        # 筛选条件：E-value、相似度和覆盖度
                        if (evalue_float <= self.config.high_quality_evalue and
                            identity >= self.config.min_identity and
                            coverage >= self.config.min_coverage):
                            high_quality_alignments.append(line)

                    except (ValueError, IndexError) as e:
                        self.logger.debug(f"跳过格式错误的行|Skipping malformed line: {e}")
                        continue

            # 写入高质量结果文件
            if high_quality_alignments:
                with open(high_quality_file, 'w', encoding='utf-8') as f:
                    f.write(header + '\n')
                    for alignment in high_quality_alignments:
                        f.write(alignment + '\n')

                self.logger.info(f"高质量结果文件已创建|High quality results file created: {high_quality_file}")
                self.logger.info(f"高质量比对数量|High quality alignments: {len(high_quality_alignments)} / {total_count}")
            else:
                self.logger.warning(f"未找到符合条件的高质量比对结果|No high quality alignments found")
                self.logger.info(f"筛选条件|Filter criteria: E-value <= {self.config.high_quality_evalue}, "
                               f"Identity >= {self.config.min_identity}%, "
                               f"Coverage >= {self.config.min_coverage}%")

            return True

        except Exception as e:
            self.logger.error(f"创建高质量结果文件失败|Failed to create high quality results file: {e}")
            return False

    def _generate_statistics_report(self) -> bool:
        """生成统计报告|Generate statistics report"""
        try:
            from .statistics import StatisticsGenerator

            self.log_step_start("生成统计报告|Generating Statistics Report", 5)

            summary_file = self.config.get_summary_output_path()

            if not os.path.exists(summary_file) or os.path.getsize(summary_file) == 0:
                self.logger.warning("汇总结果文件为空，跳过统计报告生成|Summary file is empty, skipping statistics report")
                self.log_step_end("生成统计报告|Generating Statistics Report", True)
                return True

            stats_generator = StatisticsGenerator(self.config, self.logger)
            stats_file = stats_generator.generate_statistics_report(summary_file)

            if stats_file:
                self.logger.info(f"统计报告已生成|Statistics report generated: {stats_file}")
                self.log_step_end("生成统计报告|Generating Statistics Report", True)
                return True
            else:
                self.logger.warning("统计报告生成失败，但继续执行|Statistics report generation failed, but continuing")
                self.log_step_end("生成统计报告|Generating Statistics Report", True)
                return True

        except Exception as e:
            self.logger.error(f"生成统计报告失败|Failed to generate statistics report: {e}")
            self.log_step_end("生成统计报告|Generating Statistics Report", False)
            return False

    def _generate_alignment_visualization(self) -> bool:
        """生成比对可视化|Generate alignment visualization"""
        try:
            from .alignment_visualizer import AlignmentVisualizer

            self.log_step_start("生成比对可视化|Generating Alignment Visualization", 6)

            # 准备BLAST结果数据|Prepare BLAST result data
            blast_results = []
            for sample_name, sample_file in self.sample_mapping.items():
                result_file = os.path.join(self.config.output, f"{sample_name}_{self.config.blast_type}_results.tsv")
                if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                    blast_results.append((os.path.basename(sample_file), sample_name, result_file))

            if not blast_results:
                self.logger.warning("没有可用于可视化的比对结果|No alignments available for visualization")
                self.log_step_end("生成比对可视化|Generating Alignment Visualization", True)
                return True

            # 使用AlignmentVisualizer生成可视化|Use AlignmentVisualizer
            visualizer = AlignmentVisualizer(self.config, self.logger)
            output_files = visualizer.generate_visualizations(blast_results)

            if output_files:
                self.logger.info(f"比对可视化生成完成|Alignment visualization generated successfully")
                if 'html' in output_files:
                    self.logger.info(f"  HTML索引文件|HTML index file: {output_files['html']['index']}")
                if 'text' in output_files:
                    self.logger.info(f"  文本摘要文件|Text summary file: {output_files['text']['summary']}")
            else:
                self.logger.warning("比对可视化生成失败|Failed to generate alignment visualization")

            self.log_step_end("生成比对可视化|Generating Alignment Visualization", True)
            return True

        except Exception as e:
            self.logger.error(f"生成比对可视化失败|Failed to generate alignment visualization: {e}")
            self.log_step_end("生成比对可视化|Generating Alignment Visualization", False)
            return False

    def _load_alignments_for_visualization(self, result_file: str) -> List[Dict]:
        """加载用于可视化的比对结果|Load alignments for visualization"""
        alignments = []

        try:
            with open(result_file, 'r') as f:
                next(f)  # 跳过表头

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 13:  # 13列: Sample + 12列BLAST输出
                        alignment = {
                            'sample': parts[0],
                            'qseqid': parts[1],
                            'sseqid': parts[2],
                            'pident': float(parts[3]),
                            'length': int(parts[4]),
                            'mismatch': int(parts[5]),
                            'gapopen': int(parts[6]),
                            'qstart': int(parts[7]),
                            'qend': int(parts[8]),
                            'sstart': int(parts[9]),
                            'send': int(parts[10]),
                            'evalue': float(parts[11]),
                            'bitscore': float(parts[12])
                        }

                        # 应用过滤条件|Apply filters
                        if (alignment['pident'] >= self.config.alignment_min_identity and
                            (alignment['length'] / alignment['qend'] * 100) >= self.config.alignment_min_coverage):
                            alignments.append(alignment)

            # 按样品分组并限制每个样品的最大比对数|Group by sample and limit max alignments per sample
            sample_alignments = {}
            for alignment in alignments:
                sample = alignment['sample']
                if sample not in sample_alignments:
                    sample_alignments[sample] = []
                if len(sample_alignments[sample]) < self.config.alignment_max_per_sample:
                    sample_alignments[sample].append(alignment)

            # 展平并按score排序|Flatten and sort by score
            filtered_alignments = []
            for sample_aligns in sample_alignments.values():
                filtered_alignments.extend(sample_aligns)

            filtered_alignments.sort(key=lambda x: x['bitscore'], reverse=True)
            return filtered_alignments

        except Exception as e:
            self.logger.error(f"加载比对结果失败|Failed to load alignments: {e}")
            return []

    def _generate_text_visualization(self, alignments: List[Dict], output_file: str):
        """生成文本格式可视化|Generate text format visualization"""
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write(" BLAST序列比对可视化结果|BLAST Alignment Visualization Results\n")
                f.write("=" * 80 + "\n\n")

                f.write(f"总比对数|Total alignments: {len(alignments)}\n")
                f.write(f"最小相似度|Minimum identity: {self.config.alignment_min_identity}%\n")
                f.write(f"最小覆盖度|Minimum coverage: {self.config.alignment_min_coverage}%\n")
                f.write(f"每样品最大比对数|Max alignments per sample: {self.config.alignment_max_per_sample}\n\n")

                # 按样品分组显示|Display by sample groups
                current_sample = None
                for i, alignment in enumerate(alignments, 1):
                    if current_sample != alignment['sample']:
                        current_sample = alignment['sample']
                        f.write(f"\n{'-' * 60}\n")
                        f.write(f" 样品|Sample: {current_sample}\n")
                        f.write(f"{'-' * 60}\n")

                    f.write(f"\n{i:3d}. 查询序列|Query: {alignment['qseqid']}\n")
                    f.write(f"    目标序列|Subject: {alignment['sseqid']}\n")
                    f.write(f"    相似度|Identity: {alignment['pident']:.1f}%\n")
                    f.write(f"    比对长度|Alignment length: {alignment['length']}\n")
                    f.write(f"    错配数|Mismatches: {alignment['mismatch']}\n")
                    f.write(f"    Gap数|Gaps: {alignment['gapopen']}\n")
                    f.write(f"    查询起止|Query range: {alignment['qstart']}-{alignment['qend']}\n")
                    f.write(f"    目标起止|Subject range: {alignment['sstart']}-{alignment['send']}\n")
                    f.write(f"    E-value: {alignment['evalue']:.2e}\n")
                    f.write(f"    Bit score: {alignment['bitscore']:.1f}\n")

            self.logger.info(f"文本格式比对可视化已生成|Text alignment visualization generated: {output_file}")

        except Exception as e:
            raise RuntimeError(f"生成文本可视化失败|Failed to generate text visualization: {e}")

    def _generate_html_visualization(self, alignments: List[Dict], output_file: str):
        """生成HTML格式可视化|Generate HTML format visualization"""
        try:
            # 按样品分组|Group by samples
            sample_groups = {}
            for alignment in alignments:
                if alignment['sample'] not in sample_groups:
                    sample_groups[alignment['sample']] = []
                sample_groups[alignment['sample']].append(alignment)

            html_content = self._create_html_template(sample_groups, len(alignments))

            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)

            self.logger.info(f"HTML格式比对可视化已生成|HTML alignment visualization generated: {output_file}")

        except Exception as e:
            raise RuntimeError(f"生成HTML可视化失败|Failed to generate HTML visualization: {e}")

    def _create_html_template(self, sample_groups: Dict[str, List[Dict]], total_alignments: int) -> str:
        """创建HTML模板|Create HTML template"""
        theme_styles = {
            'modern': {
                'bg': '#ffffff',
                'header': '#2c3e50',
                'sample': '#3498db',
                'border': '#ecf0f1'
            },
            'classic': {
                'bg': '#f8f9fa',
                'header': '#495057',
                'sample': '#007bff',
                'border': '#dee2e6'
            },
            'dark': {
                'bg': '#2d3748',
                'header': '#1a202c',
                'sample': '#63b3ed',
                'border': '#4a5568'
            }
        }

        theme = theme_styles.get(self.config.html_theme, theme_styles['modern'])

        html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title> BLAST序列比对可视化结果|BLAST Alignment Visualization Results</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: {theme['bg']}; color: #333; }}
        .header {{ background-color: {theme['header']}; color: white; padding: 20px; border-radius: 8px; margin-bottom: 20px; text-align: center; }}
        .summary {{ background-color: {theme['border']}; padding: 15px; border-radius: 8px; margin-bottom: 20px; }}
        .sample-group {{ margin-bottom: 30px; border: 1px solid {theme['border']}; border-radius: 8px; overflow: hidden; }}
        .sample-header {{ background-color: {theme['sample']}; color: white; padding: 15px; font-weight: bold; }}
        .alignment {{ border-bottom: 1px solid {theme['border']}; padding: 15px; }}
        .alignment:last-child {{ border-bottom: none; }}
        .alignment-title {{ font-weight: bold; color: #2c3e50; margin-bottom: 8px; }}
        .alignment-details {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; }}
        .detail-item {{ display: flex; justify-content: space-between; }}
        .detail-label {{ font-weight: 600; }}
        .detail-value {{ font-family: monospace; }}
        .high-quality {{ background-color: #d4edda; }}
        .medium-quality {{ background-color: #fff3cd; }}
        .low-quality {{ background-color: #f8d7da; }}
    </style>
</head>
<body>
    <div class="header">
        <h1> BLAST序列比对可视化结果</h1>
        <p>BLAST Alignment Visualization Results</p>
    </div>

    <div class="summary">
        <h3> 分析摘要|Analysis Summary</h3>
        <p><strong>总比对数|Total alignments:</strong> {total_alignments}</p>
        <p><strong>样品数量|Sample count:</strong> {len(sample_groups)}</p>
        <p><strong>最小相似度|Minimum identity:</strong> {self.config.alignment_min_identity}%</p>
        <p><strong>最小覆盖度|Minimum coverage:</strong> {self.config.alignment_min_coverage}%</p>
        <p><strong>每样品最大比对数|Max alignments per sample:</strong> {self.config.alignment_max_per_sample}</p>
    </div>
"""

        for sample_name, sample_alignments in sample_groups.items():
            html += f"""
    <div class="sample-group">
        <div class="sample-header">
             样品|Sample: {sample_name} ({len(sample_alignments)} 比对|alignments)
        </div>
"""

            for alignment in sample_alignments:
                # 根据相似度设置质量等级|Set quality level based on identity
                quality_class = ""
                if alignment['pident'] >= 90:
                    quality_class = "high-quality"
                elif alignment['pident'] >= 70:
                    quality_class = "medium-quality"
                else:
                    quality_class = "low-quality"

                html += f"""
        <div class="alignment {quality_class}">
            <div class="alignment-title">
                查询序列|Query: {alignment['qseqid']} → 目标序列|Subject: {alignment['sseqid']}
            </div>
            <div class="alignment-details">
                <div class="detail-item">
                    <span class="detail-label">相似度|Identity:</span>
                    <span class="detail-value">{alignment['pident']:.1f}%</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">比对长度|Length:</span>
                    <span class="detail-value">{alignment['length']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">错配数|Mismatches:</span>
                    <span class="detail-value">{alignment['mismatch']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">Gap数|Gaps:</span>
                    <span class="detail-value">{alignment['gapopen']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">查询范围|Query range:</span>
                    <span class="detail-value">{alignment['qstart']}-{alignment['qend']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">目标范围|Subject range:</span>
                    <span class="detail-value">{alignment['sstart']}-{alignment['send']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">E-value:</span>
                    <span class="detail-value">{alignment['evalue']:.2e}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">Bit score:</span>
                    <span class="detail-value">{alignment['bitscore']:.1f}</span>
                </div>
            </div>
        </div>
"""

            html += "\n    </div>\n"

        html += """
</body>
</html>"""

        return html

    def log_additional_summary(self):
        """记录额外的摘要信息|Log additional summary information"""
        self.logger.info(f"样品数量|Sample count: {len(self.sample_mapping)}")
        self.logger.info(f"BLAST类型|BLAST type: {self.config.blast_type}")
        self.logger.info(f"目标数据库|Target database: {self.config.reference}")
        self.logger.info(f"输出目录|Output directory: {self.config.output}")

        # 记录比对可视化信息|Log alignment visualization information
        if self.config.alignment_output != 'none':
            self.logger.info(f"比对可视化|Alignment visualization: {self.config.alignment_output}")
            if self.config.alignment_output in ['html', 'both']:
                html_path = self.config.get_alignment_output_path('html')
                self.logger.info(f"HTML可视化路径|HTML visualization path: {html_path}")
            if self.config.alignment_output in ['text', 'both']:
                text_path = self.config.get_alignment_output_path('text')
                self.logger.info(f"文本可视化路径|Text visualization path: {text_path}")

        if self.config.auto_generated_map:
            self.logger.info("使用自动生成的样品映射|Used auto-generated sample mapping")


def main():
    """命令行入口函数|Command line entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="BLAST序列比对分析工具|BLAST Sequence Alignment Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入文件或目录路径|Input file or directory path')
    parser.add_argument('-r', '--reference', required=True,
                       help='目标基因序列文件|Target gene sequence file')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output', default='./blast_output',
                       help='输出目录路径|Output directory path')
    parser.add_argument('-p', '--prefix', default='blast_output',
                       help='输出文件前缀|Output prefix')
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Thread count')
    parser.add_argument('-q', '--quality', type=float, default=1e-5,
                       help='E-value阈值|E-value threshold')
    parser.add_argument('-m', '--memory', default='8G',
                       help='内存限制|Memory limit')

    # 样本信息参数|Sample information parameters
    parser.add_argument('--sample-id', default=None,
                       help='样本ID|Sample ID')
    parser.add_argument('--sample-name', default=None,
                       help='样本名称|Sample name')

    # 质控参数|Quality control parameters
    parser.add_argument('--min-quality', type=float, default=20,
                       help='最小质量值|Minimum quality value')
    parser.add_argument('--min-length', type=int, default=50,
                       help='最小序列长度|Minimum sequence length')
    parser.add_argument('--min-depth', type=int, default=10,
                       help='最小测序深度|Minimum sequencing depth')
    parser.add_argument('--max-depth', type=int, default=1000,
                       help='最大测序深度|Maximum sequencing depth')
    parser.add_argument('--mapping-quality', type=int, default=20,
                       help='最小mapping质量|Minimum mapping quality')

    # 日志控制参数|Logging parameters
    parser.add_argument('-v', '--verbose', action='count', default=0,
                       help='详细输出模式|Verbose output mode')
    parser.add_argument('--quiet', action='store_true',
                       help='静默模式|Quiet mode')
    parser.add_argument('--log-file', default=None,
                       help='日志文件路径|Log file path')

    # 执行控制参数|Execution control parameters
    parser.add_argument('-f', '--force', action='store_true',
                       help='强制覆盖已存在文件|Force overwrite existing files')
    parser.add_argument('--dry-run', action='store_true',
                       help='模拟运行不执行|Dry run without execution')
    parser.add_argument('--keep-intermediate', action='store_true',
                       help='保留中间文件|Keep intermediate files')

    # BLAST特定参数|BLAST-specific parameters
    parser.add_argument('--blast-type', default='blastn',
                       choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'],
                       help='BLAST程序类型|BLAST program type')
    parser.add_argument('--max-target-seqs', type=int, default=10,
                       help='最大目标序列数|Maximum target sequences')
    parser.add_argument('--min-identity', type=float, default=70.0,
                       help='最小序列相似度|Minimum sequence identity')
    parser.add_argument('--min-coverage', type=float, default=50.0,
                       help='最小覆盖度|Minimum coverage')
    parser.add_argument('--target-db-type', default='nucl',
                       choices=['nucl', 'prot'],
                       help='目标数据库类型|Target database type')
    parser.add_argument('--high-quality-evalue', type=float, default=1e-10,
                       help='高质量比对E-value阈值|High quality alignment E-value threshold')

    # 工具路径参数|Tool path parameters
    parser.add_argument('--makeblastdb-path', default='makeblastdb',
                       help='makeblastdb程序路径|makeblastdb program path')
    parser.add_argument('--blastn-path', default='blastn',
                       help='blastn程序路径|blastn program path')
    parser.add_argument('--blastp-path', default='blastp',
                       help='blastp程序路径|blastp program path')
    parser.add_argument('--blastx-path', default='blastx',
                       help='blastx程序路径|blastx program path')
    parser.add_argument('--tblastn-path', default='tblastn',
                       help='tblastn程序路径|tblastn program path')
    parser.add_argument('--tblastx-path', default='tblastx',
                       help='tblastx程序路径|tblastx program path')

    # 比对可视化参数|Alignment visualization parameters
    parser.add_argument('--alignment-output', default='both',
                       choices=['none', 'text', 'html', 'both'],
                       help='比对可视化输出格式|Alignment visualization output format')
    parser.add_argument('--alignment-width', type=int, default=80,
                       help='比对每行显示的字符数|Characters per line in alignment display')
    parser.add_argument('--alignment-min-identity', type=float, default=0.0,
                       help='比对可视化最小相似度过滤|Minimum identity for alignment visualization')
    parser.add_argument('--alignment-min-coverage', type=float, default=0.0,
                       help='比对可视化最小覆盖度过滤|Minimum coverage for alignment visualization')
    parser.add_argument('--alignment-max-per-sample', type=int, default=100,
                       help='每个样品最多显示的比对数|Maximum alignments to display per sample')
    parser.add_argument('--html-theme', default='modern',
                       choices=['modern', 'classic', 'dark'],
                       help='HTML主题样式|HTML theme style')

    args = parser.parse_args()

    # 创建配置|Create configuration
    config = BLASTConfig(
        input=args.input,
        reference=args.reference,
        output=args.output,
        prefix=args.prefix,
        threads=args.threads,
        quality=args.quality,
        memory=args.memory,
        sample_id=args.sample_id,
        sample_name=args.sample_name,
        min_quality=args.min_quality,
        min_length=args.min_length,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        mapping_quality=args.mapping_quality,
        log_level="DEBUG" if args.verbose >= 2 else "INFO" if args.verbose == 1 else "WARNING",
        log_file=args.log_file,
        force=args.force,
        dry_run=args.dry_run,
        keep_intermediate=args.keep_intermediate,
        blast_type=args.blast_type,
        max_target_seqs=args.max_target_seqs,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        target_db_type=args.target_db_type,
        high_quality_evalue=args.high_quality_evalue,
        makeblastdb_path=args.makeblastdb_path,
        blastn_path=args.blastn_path,
        blastp_path=args.blastp_path,
        blastx_path=args.blastx_path,
        tblastn_path=args.tblastn_path,
        tblastx_path=args.tblastx_path,
        alignment_output=args.alignment_output,
        alignment_width=args.alignment_width,
        alignment_min_identity=args.alignment_min_identity,
        alignment_min_coverage=args.alignment_min_coverage,
        alignment_max_per_sample=args.alignment_max_per_sample,
        html_theme=args.html_theme
    )

    # 创建分析器并运行|Create analyzer and run
    analyzer = BLASTAnalyzer(config)
    success = analyzer.run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()