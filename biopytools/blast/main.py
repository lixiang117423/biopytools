"""
BLAST分析主程序模块 | BLAST Analysis Main Module
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
    """BLAST分析主类 | Main BLAST Analysis Class"""

    def __init__(self, config: BLASTConfig):
        """
        初始化BLAST分析器 | Initialize BLAST analyzer

        Args:
            config: BLAST配置对象 | BLAST configuration object
        """
        super().__init__(config)
        self.sample_mapping = {}
        self.database_path = None

    def validate_inputs(self) -> bool:
        """
        验证输入 | Validate inputs

        Returns:
            bool: 验证是否通过 | Whether validation passed
        """
        try:
            # 调用父类验证
            super().validate_inputs()

            # 检查依赖工具 | Check dependencies
            self._check_dependencies()

            # 生成或验证样品映射 | Generate or validate sample mapping
            self._generate_sample_mapping()

            return True

        except Exception as e:
            self.logger.error(f"输入验证失败 | Input validation failed: {e}")
            return False

    def _check_dependencies(self):
        """检查BLAST依赖工具 | Check BLAST dependencies"""
        # 获取正确的BLAST程序路径 | Get correct BLAST program path
        blast_program = getattr(self.config, f'{self.config.blast_type}_path', self.config.blast_type)

        required_tools = [self.config.makeblastdb_path, blast_program]

        for tool in required_tools:
            if not shutil.which(tool):
                raise RuntimeError(f"未找到依赖工具 | Required tool not found: {tool}")

    def _generate_sample_mapping(self):
        """生成或验证样品映射 | Generate or validate sample mapping"""
        if self.config.sample_map_file:
            # 使用现有的样品映射文件 | Use existing sample mapping file
            self._load_sample_mapping()
        elif self.config.input:
            # 自动生成样品映射文件 | Auto-generate sample mapping
            self._auto_generate_sample_mapping()
        else:
            raise ValueError("必须指定输入文件或样品映射文件 | Must specify input file or sample mapping file")

    def _load_sample_mapping(self):
        """加载样品映射文件 | Load sample mapping file"""
        self.logger.info(f"加载样品映射文件 | Loading sample mapping file: {self.config.sample_map_file}")

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
                            self.logger.warning(f"文件不存在 | File not found: {file_path}")
                    else:
                        self.logger.warning(f"第{line_num}行格式错误 | Line {line_num} format error: {line}")

            if not self.sample_mapping:
                raise ValueError("样品映射文件为空或格式错误 | Sample mapping file is empty or incorrectly formatted")

            self.logger.info(f"加载了{len(self.sample_mapping)}个样品 | Loaded {len(self.sample_mapping)} samples")

        except Exception as e:
            raise RuntimeError(f"加载样品映射文件失败 | Failed to load sample mapping file: {e}")

    def _auto_generate_sample_mapping(self):
        """自动生成样品映射 | Auto-generate sample mapping"""
        self.logger.info(f"自动扫描输入文件 | Auto-scanning input files: {self.config.input}")

        input_path = Path(self.config.input)

        if input_path.is_file():
            # 单个文件 | Single file
            sample_name = self.config.sample_name or input_path.stem
            self.sample_mapping[sample_name] = str(input_path)
            self.config.auto_generated_map = True

        elif input_path.is_dir():
            # 目录，扫描匹配的文件 | Directory, scan matching files
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
            raise FileNotFoundError(f"输入路径不存在 | Input path not found: {self.config.input}")

        if not self.sample_mapping:
            raise ValueError(f"在{self.config.input}中未找到匹配的文件 | No matching files found in {self.config.input}")

        self.logger.info(f"发现了{len(self.sample_mapping)}个样品 | Found {len(self.sample_mapping)} samples")

    def run_analysis(self) -> bool:
        """
        运行BLAST分析流程 | Run BLAST analysis pipeline

        Returns:
            bool: 分析是否成功 | Whether analysis was successful
        """
        try:
            # 步骤1: 创建输出目录 | Step 1: Create output directory
            if not self._create_output_directory():
                return False

            # 步骤2: 创建BLAST数据库 | Step 2: Create BLAST database
            if not self._create_blast_database():
                return False

            # 步骤3: 运行BLAST比对 | Step 3: Run BLAST alignment
            if not self._run_blast_alignment():
                return False

            # 步骤4: 处理结果 | Step 4: Process results
            if not self._process_results():
                return False

            # 步骤5: 生成比对可视化 | Step 5: Generate alignment visualization
            if self.config.alignment_output != 'none':
                if not self._generate_alignment_visualization():
                    return False

            return True

        except Exception as e:
            self.logger.error(f"BLAST分析失败 | BLAST analysis failed: {e}", exc_info=True)
            return False

    def _create_output_directory(self) -> bool:
        """创建输出目录 | Create output directory"""
        try:
            os.makedirs(self.config.output, exist_ok=True)
            self.logger.info(f"输出目录已创建 | Output directory created: {self.config.output}")
            return True
        except Exception as e:
            self.logger.error(f"创建输出目录失败 | Failed to create output directory: {e}")
            return False

    def _create_blast_database(self) -> bool:
        """创建BLAST数据库 | Create BLAST database"""
        try:
            self.log_step_start("创建BLAST数据库 | Creating BLAST Database", 1)

            db_path = self.config.get_target_db_path()

            # 检查数据库是否已存在 | Check if database already exists
            if os.path.exists(db_path) and not self.config.force:
                self.logger.info(f"数据库已存在，跳过创建 | Database already exists, skipping creation: {db_path}")
                self.database_path = db_path
                self.log_step_end("创建BLAST数据库 | Creating BLAST Database", True)
                return True

            # 构建数据库 | Build database
            cmd = [
                self.config.makeblastdb_path,
                '-in', self.config.reference,
                '-dbtype', self.config.target_db_type,
                '-out', db_path
            ]

            self.logger.info(f"执行命令 | Executing command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            self.database_path = db_path
            self.logger.info(f"数据库创建成功 | Database created successfully: {db_path}")

            self.log_step_end("创建BLAST数据库 | Creating BLAST Database", True)
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"创建数据库失败 | Database creation failed: {e}")
            if e.stdout:
                self.logger.error(f"stdout: {e.stdout}")
            if e.stderr:
                self.logger.error(f"stderr: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"创建数据库时发生错误 | Error creating database: {e}")
            return False

    def _run_blast_alignment(self) -> bool:
        """运行BLAST比对 | Run BLAST alignment"""
        try:
            self.log_step_start("运行BLAST比对 | Running BLAST Alignment", 2)

            all_results = []
            total_samples = len(self.sample_mapping)

            for i, (sample_name, sample_file) in enumerate(self.sample_mapping.items(), 1):
                self.logger.info(f"处理样品 {i}/{total_samples}: {sample_name} | Processing sample {i}/{total_samples}: {sample_name}")

                result_file = self._run_single_blast(sample_name, sample_file)
                if result_file:
                    all_results.append(result_file)
                else:
                    self.logger.warning(f"样品 {sample_name} 比对失败 | Sample {sample_name} alignment failed")

            if not all_results:
                raise RuntimeError("所有样品比对都失败了 | All sample alignments failed")

            # 合并结果 | Merge results
            merged_file = self._merge_results(all_results)

            self.logger.info(f"BLAST比对完成 | BLAST alignment completed. Results: {merged_file}")

            self.log_step_end("运行BLAST比对 | Running BLAST Alignment", True)
            return True

        except Exception as e:
            self.logger.error(f"BLAST比对失败 | BLAST alignment failed: {e}")
            self.log_step_end("运行BLAST比对 | Running BLAST Alignment", False)
            return False

    def _run_single_blast(self, sample_name: str, sample_file: str) -> Optional[str]:
        """运行单个样品的BLAST比对 | Run BLAST alignment for single sample"""
        try:
            output_file = os.path.join(self.config.output, f"{sample_name}_{self.config.blast_type}_results.tsv")

            # 获取正确的BLAST程序路径 | Get correct BLAST program path
            blast_program = getattr(self.config, f'{self.config.blast_type}_path', self.config.blast_type)

            cmd = [
                blast_program,
                '-query', sample_file,
                '-db', self.database_path,
                '-out', output_file,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                '-evalue', str(self.config.evalue),
                '-max_target_seqs', str(self.config.max_target_seqs),
                '-word_size', str(self.config.word_size),
                '-num_threads', str(self.config.threads)
            ]

            self.logger.debug(f"执行BLAST命令 | Executing BLAST command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # 检查输出文件 | Check output file
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                self.logger.info(f"样品 {sample_name} 比对完成 | Sample {sample_name} alignment completed")
                return output_file
            else:
                self.logger.warning(f"样品 {sample_name} 比对输出为空 | Sample {sample_name} alignment output is empty")
                return None

        except subprocess.CalledProcessError as e:
            self.logger.error(f"样品 {sample_name} BLAST比对失败 | Sample {sample_name} BLAST alignment failed: {e}")
            if e.stderr:
                self.logger.error(f"stderr: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"样品 {sample_name} 比对时发生错误 | Error during sample {sample_name} alignment: {e}")
            return None

    def _merge_results(self, result_files: List[str]) -> str:
        """合并所有结果文件 | Merge all result files"""
        try:
            merged_file = self.config.get_summary_output_path()

            with open(merged_file, 'w') as out_f:
                # 写入表头 | Write header
                out_f.write("Sample\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n")

                for result_file in result_files:
                    sample_name = os.path.basename(result_file).split('_')[0]

                    try:
                        with open(result_file, 'r') as in_f:
                            for line in in_f:
                                out_f.write(f"{sample_name}\t{line}")
                    except Exception as e:
                        self.logger.warning(f"读取结果文件失败 | Failed to read result file: {result_file}, error: {e}")

            return merged_file

        except Exception as e:
            raise RuntimeError(f"合并结果文件失败 | Failed to merge result files: {e}")

    def _process_results(self) -> bool:
        """处理结果 | Process results"""
        try:
            self.log_step_start("处理结果 | Processing Results", 3)

            summary_file = self.config.get_summary_output_path()

            if not os.path.exists(summary_file):
                raise FileNotFoundError(f"汇总结果文件不存在 | Summary result file not found: {summary_file}")

            # 统计结果 | Generate statistics
            stats = self._generate_statistics(summary_file)
            self.log_statistics(stats)

            self.log_step_end("处理结果 | Processing Results", True)
            return True

        except Exception as e:
            self.logger.error(f"处理结果失败 | Failed to process results: {e}")
            self.log_step_end("处理结果 | Processing Results", False)
            return False

    def _generate_statistics(self, result_file: str) -> Dict:
        """生成统计信息 | Generate statistics"""
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
                    if len(parts) >= 14:
                        stats['total_alignments'] += 1
                        stats['unique_queries'].add(parts[1])  # qseqid
                        stats['unique_subjects'].add(parts[2])  # sseqid

            stats['unique_queries'] = len(stats['unique_queries'])
            stats['unique_subjects'] = len(stats['unique_subjects'])

            return stats

        except Exception as e:
            self.logger.error(f"生成统计信息失败 | Failed to generate statistics: {e}")
            return {'total_alignments': 0, 'samples_count': 0, 'unique_queries': 0, 'unique_subjects': 0}

    def _generate_alignment_visualization(self) -> bool:
        """生成比对可视化 | Generate alignment visualization"""
        try:
            self.log_step_start("生成比对可视化 | Generating Alignment Visualization", 4)

            summary_file = self.config.get_summary_output_path()

            if not os.path.exists(summary_file):
                raise FileNotFoundError(f"汇总结果文件不存在 | Summary result file not found: {summary_file}")

            # 加载BLAST结果 | Load BLAST results
            alignments = self._load_alignments_for_visualization(summary_file)

            if not alignments:
                self.logger.warning("没有可用于可视化的比对结果 | No alignments available for visualization")
                self.log_step_end("生成比对可视化 | Generating Alignment Visualization", True)
                return True

            # 生成文本格式可视化 | Generate text format visualization
            if self.config.alignment_output in ['text', 'both']:
                text_output = self.config.get_alignment_output_path("text")
                self._generate_text_visualization(alignments, text_output)

            # 生成HTML格式可视化 | Generate HTML format visualization
            if self.config.alignment_output in ['html', 'both']:
                html_output = self.config.get_alignment_output_path("html")
                self._generate_html_visualization(alignments, html_output)

            self.log_step_end("生成比对可视化 | Generating Alignment Visualization", True)
            return True

        except Exception as e:
            self.logger.error(f"生成比对可视化失败 | Failed to generate alignment visualization: {e}")
            self.log_step_end("生成比对可视化 | Generating Alignment Visualization", False)
            return False

    def _load_alignments_for_visualization(self, result_file: str) -> List[Dict]:
        """加载用于可视化的比对结果 | Load alignments for visualization"""
        alignments = []

        try:
            with open(result_file, 'r') as f:
                next(f)  # 跳过表头

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 14:
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

                        # 应用过滤条件 | Apply filters
                        if (alignment['pident'] >= self.config.alignment_min_identity and
                            (alignment['length'] / alignment['qend'] * 100) >= self.config.alignment_min_coverage):
                            alignments.append(alignment)

            # 按样品分组并限制每个样品的最大比对数 | Group by sample and limit max alignments per sample
            sample_alignments = {}
            for alignment in alignments:
                sample = alignment['sample']
                if sample not in sample_alignments:
                    sample_alignments[sample] = []
                if len(sample_alignments[sample]) < self.config.alignment_max_per_sample:
                    sample_alignments[sample].append(alignment)

            # 展平并按score排序 | Flatten and sort by score
            filtered_alignments = []
            for sample_aligns in sample_alignments.values():
                filtered_alignments.extend(sample_aligns)

            filtered_alignments.sort(key=lambda x: x['bitscore'], reverse=True)
            return filtered_alignments

        except Exception as e:
            self.logger.error(f"加载比对结果失败 | Failed to load alignments: {e}")
            return []

    def _generate_text_visualization(self, alignments: List[Dict], output_file: str):
        """生成文本格式可视化 | Generate text format visualization"""
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("🧬 BLAST序列比对可视化结果 | BLAST Alignment Visualization Results\n")
                f.write("=" * 80 + "\n\n")

                f.write(f"总比对数 | Total alignments: {len(alignments)}\n")
                f.write(f"最小相似度 | Minimum identity: {self.config.alignment_min_identity}%\n")
                f.write(f"最小覆盖度 | Minimum coverage: {self.config.alignment_min_coverage}%\n")
                f.write(f"每样品最大比对数 | Max alignments per sample: {self.config.alignment_max_per_sample}\n\n")

                # 按样品分组显示 | Display by sample groups
                current_sample = None
                for i, alignment in enumerate(alignments, 1):
                    if current_sample != alignment['sample']:
                        current_sample = alignment['sample']
                        f.write(f"\n{'-' * 60}\n")
                        f.write(f"📂 样品 | Sample: {current_sample}\n")
                        f.write(f"{'-' * 60}\n")

                    f.write(f"\n{i:3d}. 查询序列 | Query: {alignment['qseqid']}\n")
                    f.write(f"    目标序列 | Subject: {alignment['sseqid']}\n")
                    f.write(f"    相似度 | Identity: {alignment['pident']:.1f}%\n")
                    f.write(f"    比对长度 | Alignment length: {alignment['length']}\n")
                    f.write(f"    错配数 | Mismatches: {alignment['mismatch']}\n")
                    f.write(f"    Gap数 | Gaps: {alignment['gapopen']}\n")
                    f.write(f"    查询起止 | Query range: {alignment['qstart']}-{alignment['qend']}\n")
                    f.write(f"    目标起止 | Subject range: {alignment['sstart']}-{alignment['send']}\n")
                    f.write(f"    E-value: {alignment['evalue']:.2e}\n")
                    f.write(f"    Bit score: {alignment['bitscore']:.1f}\n")

            self.logger.info(f"文本格式比对可视化已生成 | Text alignment visualization generated: {output_file}")

        except Exception as e:
            raise RuntimeError(f"生成文本可视化失败 | Failed to generate text visualization: {e}")

    def _generate_html_visualization(self, alignments: List[Dict], output_file: str):
        """生成HTML格式可视化 | Generate HTML format visualization"""
        try:
            # 按样品分组 | Group by samples
            sample_groups = {}
            for alignment in alignments:
                if alignment['sample'] not in sample_groups:
                    sample_groups[alignment['sample']] = []
                sample_groups[alignment['sample']].append(alignment)

            html_content = self._create_html_template(sample_groups, len(alignments))

            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)

            self.logger.info(f"HTML格式比对可视化已生成 | HTML alignment visualization generated: {output_file}")

        except Exception as e:
            raise RuntimeError(f"生成HTML可视化失败 | Failed to generate HTML visualization: {e}")

    def _create_html_template(self, sample_groups: Dict[str, List[Dict]], total_alignments: int) -> str:
        """创建HTML模板 | Create HTML template"""
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
    <title>🧬 BLAST序列比对可视化结果 | BLAST Alignment Visualization Results</title>
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
        <h1>🧬 BLAST序列比对可视化结果</h1>
        <p>BLAST Alignment Visualization Results</p>
    </div>

    <div class="summary">
        <h3>📊 分析摘要 | Analysis Summary</h3>
        <p><strong>总比对数 | Total alignments:</strong> {total_alignments}</p>
        <p><strong>样品数量 | Sample count:</strong> {len(sample_groups)}</p>
        <p><strong>最小相似度 | Minimum identity:</strong> {self.config.alignment_min_identity}%</p>
        <p><strong>最小覆盖度 | Minimum coverage:</strong> {self.config.alignment_min_coverage}%</p>
        <p><strong>每样品最大比对数 | Max alignments per sample:</strong> {self.config.alignment_max_per_sample}</p>
    </div>
"""

        for sample_name, sample_alignments in sample_groups.items():
            html += f"""
    <div class="sample-group">
        <div class="sample-header">
            📂 样品 | Sample: {sample_name} ({len(sample_alignments)} 比对 | alignments)
        </div>
"""

            for alignment in sample_alignments:
                # 根据相似度设置质量等级 | Set quality level based on identity
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
                查询序列 | Query: {alignment['qseqid']} → 目标序列 | Subject: {alignment['sseqid']}
            </div>
            <div class="alignment-details">
                <div class="detail-item">
                    <span class="detail-label">相似度 | Identity:</span>
                    <span class="detail-value">{alignment['pident']:.1f}%</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">比对长度 | Length:</span>
                    <span class="detail-value">{alignment['length']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">错配数 | Mismatches:</span>
                    <span class="detail-value">{alignment['mismatch']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">Gap数 | Gaps:</span>
                    <span class="detail-value">{alignment['gapopen']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">查询范围 | Query range:</span>
                    <span class="detail-value">{alignment['qstart']}-{alignment['qend']}</span>
                </div>
                <div class="detail-item">
                    <span class="detail-label">目标范围 | Subject range:</span>
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
        """记录额外的摘要信息 | Log additional summary information"""
        self.logger.info(f"样品数量 | Sample count: {len(self.sample_mapping)}")
        self.logger.info(f"BLAST类型 | BLAST type: {self.config.blast_type}")
        self.logger.info(f"目标数据库 | Target database: {self.config.reference}")
        self.logger.info(f"输出目录 | Output directory: {self.config.output}")

        # 记录比对可视化信息 | Log alignment visualization information
        if self.config.alignment_output != 'none':
            self.logger.info(f"比对可视化 | Alignment visualization: {self.config.alignment_output}")
            if self.config.alignment_output in ['html', 'both']:
                html_path = self.config.get_alignment_output_path('html')
                self.logger.info(f"HTML可视化路径 | HTML visualization path: {html_path}")
            if self.config.alignment_output in ['text', 'both']:
                text_path = self.config.get_alignment_output_path('text')
                self.logger.info(f"文本可视化路径 | Text visualization path: {text_path}")

        if self.config.auto_generated_map:
            self.logger.info("使用自动生成的样品映射 | Used auto-generated sample mapping")


def main():
    """命令行入口函数 | Command line entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description="BLAST序列比对分析工具 | BLAST Sequence Alignment Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入文件或目录路径 | Input file or directory path')
    parser.add_argument('-r', '--reference', required=True,
                       help='目标基因序列文件 | Target gene sequence file')

    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./blast_output',
                       help='输出目录路径 | Output directory path (default: ./blast_output)')
    parser.add_argument('-p', '--prefix', default='blast_output',
                       help='输出文件前缀 | Output prefix (default: blast_output)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='线程数 | Thread count (default: 4)')
    parser.add_argument('-q', '--quality', type=float, default=1e-5,
                       help='E-value阈值 | E-value threshold (default: 1e-5)')
    parser.add_argument('-m', '--memory', default='8G',
                       help='内存限制 | Memory limit (default: 8G)')

    # 样本信息参数 | Sample information parameters
    parser.add_argument('--sample-id', default=None,
                       help='样本ID | Sample ID')
    parser.add_argument('--sample-name', default=None,
                       help='样本名称 | Sample name')

    # 质控参数 | Quality control parameters
    parser.add_argument('--min-quality', type=float, default=20,
                       help='最小质量值 | Minimum quality value (default: 20)')
    parser.add_argument('--min-length', type=int, default=50,
                       help='最小序列长度 | Minimum sequence length (default: 50)')
    parser.add_argument('--min-depth', type=int, default=10,
                       help='最小测序深度 | Minimum sequencing depth (default: 10)')
    parser.add_argument('--max-depth', type=int, default=1000,
                       help='最大测序深度 | Maximum sequencing depth (default: 1000)')
    parser.add_argument('--mapping-quality', type=int, default=20,
                       help='最小mapping质量 | Minimum mapping quality (default: 20)')

    # 日志控制参数 | Logging parameters
    parser.add_argument('-v', '--verbose', action='count', default=0,
                       help='详细输出模式 | Verbose output mode')
    parser.add_argument('--quiet', action='store_true',
                       help='静默模式 | Quiet mode')
    parser.add_argument('--log-file', default=None,
                       help='日志文件路径 | Log file path')

    # 执行控制参数 | Execution control parameters
    parser.add_argument('-f', '--force', action='store_true',
                       help='强制覆盖已存在文件 | Force overwrite existing files')
    parser.add_argument('--dry-run', action='store_true',
                       help='模拟运行不执行 | Dry run without execution')
    parser.add_argument('--keep-intermediate', action='store_true',
                       help='保留中间文件 | Keep intermediate files')

    # BLAST特定参数 | BLAST-specific parameters
    parser.add_argument('--blast-type', default='blastn',
                       choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'],
                       help='BLAST程序类型 | BLAST program type (default: blastn)')
    parser.add_argument('--max-target-seqs', type=int, default=10,
                       help='最大目标序列数 | Maximum target sequences (default: 10)')
    parser.add_argument('--min-identity', type=float, default=70.0,
                       help='最小序列相似度 | Minimum sequence identity (default: 70.0)')
    parser.add_argument('--min-coverage', type=float, default=50.0,
                       help='最小覆盖度 | Minimum coverage (default: 50.0)')
    parser.add_argument('--target-db-type', default='nucl',
                       choices=['nucl', 'prot'],
                       help='目标数据库类型 | Target database type (default: nucl)')
    parser.add_argument('--high-quality-evalue', type=float, default=1e-10,
                       help='高质量比对E-value阈值 | High quality alignment E-value threshold (default: 1e-10)')

    # 工具路径参数 | Tool path parameters
    parser.add_argument('--makeblastdb-path', default='makeblastdb',
                       help='makeblastdb程序路径 | makeblastdb program path (default: makeblastdb)')
    parser.add_argument('--blastn-path', default='blastn',
                       help='blastn程序路径 | blastn program path (default: blastn)')
    parser.add_argument('--blastp-path', default='blastp',
                       help='blastp程序路径 | blastp program path (default: blastp)')
    parser.add_argument('--blastx-path', default='blastx',
                       help='blastx程序路径 | blastx program path (default: blastx)')
    parser.add_argument('--tblastn-path', default='tblastn',
                       help='tblastn程序路径 | tblastn program path (default: tblastn)')
    parser.add_argument('--tblastx-path', default='tblastx',
                       help='tblastx程序路径 | tblastx program path (default: tblastx)')

    # 比对可视化参数 | Alignment visualization parameters
    parser.add_argument('--alignment-output', default='both',
                       choices=['none', 'text', 'html', 'both'],
                       help='比对可视化输出格式 | Alignment visualization output format (default: both)')
    parser.add_argument('--alignment-width', type=int, default=80,
                       help='比对每行显示的字符数 | Characters per line in alignment display (default: 80)')
    parser.add_argument('--alignment-min-identity', type=float, default=0.0,
                       help='比对可视化最小相似度过滤 | Minimum identity for alignment visualization (default: 0.0)')
    parser.add_argument('--alignment-min-coverage', type=float, default=0.0,
                       help='比对可视化最小覆盖度过滤 | Minimum coverage for alignment visualization (default: 0.0)')
    parser.add_argument('--alignment-max-per-sample', type=int, default=100,
                       help='每个样品最多显示的比对数 | Maximum alignments to display per sample (default: 100)')
    parser.add_argument('--html-theme', default='modern',
                       choices=['modern', 'classic', 'dark'],
                       help='HTML主题样式 | HTML theme style (default: modern)')

    args = parser.parse_args()

    # 创建配置 | Create configuration
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

    # 创建分析器并运行 | Create analyzer and run
    analyzer = BLASTAnalyzer(config)
    success = analyzer.run_pipeline()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()