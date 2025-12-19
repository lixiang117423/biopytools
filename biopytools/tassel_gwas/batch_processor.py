"""
TASSEL GWAS批量处理器 | TASSEL GWAS Batch Processor
"""

import os
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

from .core import TASSELGWASAnalyzer, TASSELGWASConfig
from .utils import TASSELLogger


class TASSELBatchProcessor:
    """TASSEL GWAS批量处理器 | TASSEL GWAS Batch Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.results = []
        self.failed_traits = []

    def extract_single_trait(self, pheno_source: Path, trait_column: int, output_dir: Path) -> Tuple[Optional[Path], Optional[str]]:
        """提取单个表型 | Extract single trait"""
        try:
            # 读取表头获取表型名称 | Read header to get trait name
            with open(pheno_source, 'r') as f:
                header_line = f.readline().strip()
                # 先尝试制表符分割，如果只有一列则尝试空格分割
                columns = header_line.split('\t')
                use_space_delim = False
                if len(columns) <= 1:
                    # 尝试空格分割（支持多个连续空格）
                    columns = header_line.split()
                    use_space_delim = True

            if trait_column >= len(columns):
                return None, f"列索引超出范围 | Column index out of range: {trait_column} >= {len(columns)}"

            trait_name = columns[trait_column]
            # 清理文件名中的特殊字符 | Clean special characters in filename
            trait_name_clean = trait_name.replace('/', '_').replace('\\', '_').replace(' ', '_')

            # 创建表型文件 | Create phenotype file
            single_pheno_file = output_dir / f"{trait_name_clean}.pheno.txt"

            self.logger.info(f"正在处理表型 | Processing trait: {trait_name} (列索引 | column index: {trait_column + 1})")

            # 提取第1列和第trait_column+1列 | Extract column 1 and trait_column+1
            # 根据分隔符类型选择合适的awk命令
            if use_space_delim:
                # 对于空格分隔的文件，使用默认的空格分隔
                awk_command = f"awk -v col={trait_column + 1} 'BEGIN{{OFS=\"\\t\"}} {{print $1, $col}}' '{pheno_source}' > '{single_pheno_file}'"
            else:
                # 对于制表符分隔的文件，明确指定制表符
                awk_command = f"awk -v col={trait_column + 1} 'BEGIN{{FS=\"\\t\"; OFS=\"\\t\"}} {{print $1, $col}}' '{pheno_source}' > '{single_pheno_file}'"

            result = os.system(awk_command)

            if result != 0:
                return None, "表型数据提取失败 | Phenotype data extraction failed"

            if not single_pheno_file.exists() or single_pheno_file.stat().st_size == 0:
                return None, "生成的表型文件为空 | Generated phenotype file is empty"

            self.logger.debug(f"✅ 表型文件提取完成 | Phenotype file extracted: {single_pheno_file}")
            return single_pheno_file, None

        except Exception as e:
            return None, f"表型提取异常 | Phenotype extraction error: {str(e)}"

    def process_single_trait(self, pheno_source: Path, vcf_file: Path, trait_column: int,
                            output_base_dir: Path, global_config: Dict) -> Dict:
        """处理单个表型的完整GWAS分析 | Process complete GWAS analysis for single trait"""
        try:
            # 读取表头获取表型名称 | Read header to get trait name
            with open(pheno_source, 'r') as f:
                header_line = f.readline().strip()
                # 先尝试制表符分割，如果只有一列则尝试空格分割
                columns = header_line.split('\t')
                if len(columns) <= 1:
                    # 尝试空格分割（支持多个连续空格）
                    columns = header_line.split()

            if trait_column >= len(columns):
                return {
                    'trait_name': f'Column_{trait_column + 1}',
                    'status': 'failed',
                    'error': f"列索引超出范围 | Column index out of range: {trait_column} >= {len(columns)}"
                }

            trait_name = columns[trait_column]
            # 清理文件名中的特殊字符 | Clean special characters in filename
            trait_name_clean = trait_name.replace('/', '_').replace('\\', '_').replace(' ', '_')

            # 创建表型专用目录 | Create trait-specific directory
            trait_dir = output_base_dir / trait_name_clean
            trait_dir.mkdir(parents=True, exist_ok=True)

            self.logger.info(f"------------------------------------------------------")
            self.logger.info(f"开始处理表型 | Starting trait processing: {trait_name} (列索引 | column index: {trait_column + 1})")

            # 提取表型数据 | Extract phenotype data
            single_pheno_file, error = self.extract_single_trait(pheno_source, trait_column, trait_dir)
            if error:
                return {
                    'trait_name': trait_name,
                    'status': 'failed',
                    'error': error
                }

            # 配置GWAS分析器 | Configure GWAS analyzer
            analyzer_config = {
                'vcf_file': str(vcf_file),
                'pheno_file': str(single_pheno_file),
                'output_prefix': f'{trait_name_clean}_GWAS',
                'output_dir': str(trait_dir),
                'model': global_config.get('model', 'MLM'),
                'memory_max': global_config.get('memory_max', '100g'),
                'threads': global_config.get('threads', 4),
                'maf_filter': global_config.get('maf_filter'),
                'miss_filter': global_config.get('miss_filter'),
                'skip_sort': True,  # 跳过排序，因为批量处理时已经处理过
                'pca_components': global_config.get('pca_components', 5),  # 添加PCA数量参数
                'q_matrix': global_config.get('precomputed_q_matrix') or global_config.get('q_matrix'),
                'kinship': global_config.get('precomputed_kinship') or global_config.get('kinship'),
                'keep_temp': global_config.get('keep_temp', True),  # 批量处理时默认保留临时文件
                'log_level': global_config.get('log_level', 'INFO')
            }

            # 创建分析器 | Create analyzer
            analyzer = TASSELGWASAnalyzer(**analyzer_config)

            self.logger.info(f"   使用预计算矩阵 | Using pre-computed matrices:")
            if analyzer_config['kinship']:
                self.logger.info(f"   • K矩阵 | K matrix: {analyzer_config['kinship']}")
            if analyzer_config['q_matrix']:
                self.logger.info(f"   • PCA协变量 | PCA covariate: {analyzer_config['q_matrix']}")
            else:
                self.logger.info(f"   • PCA协变量: 未使用 | PCA covariate: Not used")

            # 运行GWAS分析 | Run GWAS analysis
            success = analyzer.run_analysis()

            if success:
                self.logger.info(f"✅ 表型 {trait_name} 运行成功 | Trait {trait_name} completed successfully")
                return {
                    'trait_name': trait_name,
                    'status': 'success',
                    'output_dir': str(trait_dir),
                    'output_prefix': f'{trait_name_clean}_GWAS'
                }
            else:
                self.logger.error(f"❌ 表型 {trait_name} 运行失败 | Trait {trait_name} failed")
                return {
                    'trait_name': trait_name,
                    'status': 'failed',
                    'error': 'GWAS分析失败 | GWAS analysis failed'
                }

        except Exception as e:
            self.logger.error(f"❌ 处理表型时发生异常 | Exception occurred while processing trait: {str(e)}")
            trait_name = f'Trait_{trait_column + 1}'
            return {
                'trait_name': trait_name,
                'status': 'failed',
                'error': f'处理异常 | Processing error: {str(e)}'
            }

    def process_batch(self, pheno_source: Path, vcf_file: Path, output_base_dir: Path,
                      global_config: Dict, parallel: bool = False, max_workers: int = 4):
        """批量处理所有表型 | Batch process all traits"""
        try:
            # 获取表型文件的表头信息 | Get phenotype file header information
            with open(pheno_source, 'r') as f:
                header_line = f.readline().strip()
                # 先尝试制表符分割，如果只有一列则尝试空格分割
                columns = header_line.split('\t')
                if len(columns) <= 1:
                    # 尝试空格分割（支持多个连续空格）
                    columns = header_line.split()
                total_cols = len(columns)

            self.logger.info(f"检测到总列数 | Detected total columns: {total_cols}")
            self.logger.info(f"表头信息 | Header info: {columns}")

            # 检查第一列是否为<Trait> | Check if first column is <Trait>
            if total_cols <= 1:
                self.logger.error(f"❌ 表型文件至少需要2列 | Phenotype file needs at least 2 columns")
                return False

            if columns[0] != "<Trait>":
                self.logger.warning(f"⚠️ 第一列不是<Trait>，将尝试自动识别表型 | First column is not <Trait>, will try auto-detection")

            # 确定要处理的表型列范围 | Determine trait columns to process
            trait_columns = list(range(1, total_cols))  # 从第2列开始

            if not trait_columns:
                self.logger.error(f"❌ 没有找到表型数据列 | No trait data columns found")
                return False

            self.logger.info(f"开始处理表型列 | Starting to process trait columns: {trait_columns}")
            self.logger.info(f"表型名称 | Trait names: {[columns[i] for i in trait_columns]}")

            # 预计算K矩阵和PCA协变量（如果需要MLM模型）| Pre-compute K matrix and PCA covariates (if MLM model needed)
            model = global_config.get('model', 'MLM')
            kinship_file = None
            pca_covariate_file = None

            if model in ['MLM', 'BOTH']:
                self.logger.info("🔬 检测到MLM模型，预计算共享矩阵 | Detected MLM model, pre-computing shared matrices")

                # 创建临时分析器用于矩阵计算 | Create temporary analyzer for matrix calculation
                temp_pheno_file = output_base_dir / "temp_pheno.txt"

                # 创建临时表型文件（使用第一个表型）| Create temporary phenotype file (using first trait)
                try:
                    temp_pheno, _ = self.extract_single_trait(pheno_source, 1, output_base_dir)
                    if temp_pheno:
                        temp_pheno_file = temp_pheno
                    else:
                        # 如果提取失败，创建一个简单的临时表型文件
                        with open(temp_pheno_file, 'w') as f:
                            f.write("<Trait>\t<Trait>\n")
                            f.write("<Data>\t<Data>\n")
                            f.write("temp_sample\t1.0\n")
                except Exception as e:
                    self.logger.warning(f"临时表型文件创建失败，继续使用默认文件 | Temporary phenotype file creation failed: {e}")

                # 配置临时分析器 | Configure temporary analyzer
                temp_config = {
                    'vcf_file': str(vcf_file),
                    'pheno_file': str(temp_pheno_file),
                    'output_prefix': 'temp_matrix_calculation',
                    'output_dir': str(output_base_dir),
                    'model': 'MLM',
                    'memory_max': global_config.get('memory_max', '100g'),
                    'threads': global_config.get('threads', 4),
                    'maf_filter': global_config.get('maf_filter'),
                    'miss_filter': global_config.get('miss_filter'),
                    'skip_sort': False,  # 矩阵预计算时需要完整的处理流程
                    'pca_components': global_config.get('pca_components', 5),  # 添加PCA数量参数
                    'q_matrix': None,  # 不使用外部Q矩阵，需要自己计算
                    'kinship': None,   # 不使用外部K矩阵，需要自己计算
                    'keep_temp': True,
                    'log_level': global_config.get('log_level', 'INFO')
                }

                temp_analyzer = TASSELGWASAnalyzer(**temp_config)

                # 计算K矩阵 | Calculate K matrix
                if global_config.get('kinship') and Path(global_config.get('kinship')).exists():
                    kinship_file = Path(global_config.get('kinship'))
                    self.logger.info(f"使用现有的K矩阵 | Using existing K matrix: {kinship_file}")
                else:
                    try:
                        self.logger.info("🔄 计算K矩阵（Kinship）| Calculating K matrix (Kinship)...")
                        vcf_sorted = temp_analyzer.sort_and_filter_vcf()
                        kinship_file = temp_analyzer.calculate_kinship(vcf_sorted)
                        self.logger.info(f"✅ K矩阵计算完成 | K matrix calculated: {kinship_file}")
                    except Exception as e:
                        self.logger.error(f"❌ K矩阵计算失败 | K matrix calculation failed: {e}")
                        kinship_file = None

                # 计算PCA协变量 | Calculate PCA covariates
                if global_config.get('q_matrix') and Path(global_config.get('q_matrix')).exists():
                    pca_covariate_file = Path(global_config.get('q_matrix'))
                    self.logger.info(f"使用现有的PCA协变量 | Using existing PCA covariate: {pca_covariate_file}")
                else:
                    try:
                        self.logger.info("🔄 计算PCA协变量（用于MLM）| Calculating PCA covariates (for MLM)...")
                        pca_covariate_file = temp_analyzer.calculate_pca(vcf_file)
                        self.logger.info(f"✅ PCA协变量计算完成 | PCA covariate calculated: {pca_covariate_file}")
                    except Exception as e:
                        self.logger.warning(f"⚠️ PCA协变量计算失败，将仅使用K矩阵 | PCA covariate calculation failed, will use K matrix only: {e}")
                        pca_covariate_file = None

                # 清理临时文件 | Clean up temporary files
                try:
                    if temp_pheno_file.exists() and temp_pheno_file.name.startswith("temp"):
                        temp_pheno_file.unlink()
                except Exception:
                    pass

                # 将预计算的矩阵路径添加到全局配置中 | Add pre-computed matrix paths to global config
                global_config['precomputed_kinship'] = str(kinship_file) if kinship_file else None
                global_config['precomputed_q_matrix'] = str(pca_covariate_file) if pca_covariate_file else None

                self.logger.info("🎯 矩阵预计算完成，所有表型将共享这些矩阵 | Matrix pre-computation completed, all traits will share these matrices")

            if parallel and max_workers > 1:
                # 并行处理 | Parallel processing
                self.logger.info(f"🚀 使用并行处理，工作线程数 | Using parallel processing with workers: {max_workers}")

                # 创建失败日志文件 | Create failed log file
                failed_log_file = output_base_dir / "failed_traits.log"
                failed_traits = []

                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    # 提交所有任务 | Submit all tasks
                    future_to_trait = {}
                    for i in trait_columns:  # 使用trait_columns而不是range
                        future = executor.submit(self.process_single_trait, pheno_source, vcf_file, i,
                                                output_base_dir, global_config)
                        future_to_trait[future] = i  # 列索引

                    # 收集结果 | Collect results
                    for future in as_completed(future_to_trait):
                        col_idx = future_to_trait[future]
                        try:
                            result = future.result()
                            self.results.append(result)
                            if result['status'] == 'failed':
                                failed_traits.append(result)
                                with open(failed_log_file, 'a') as f:
                                    f.write(f"{result['trait_name']}: {result['error']}\n")
                        except Exception as e:
                            self.logger.error(f"处理第{col_idx+1}列时发生异常 | Exception occurred while processing column {col_idx+1}: {e}")
                            failed_trait = {
                                'trait_name': f'Column_{col_idx+1}',
                                'status': 'failed',
                                'error': f'处理异常 | Processing error: {str(e)}'
                            }
                            self.results.append(failed_trait)
                            failed_traits.append(failed_trait)
                            with open(failed_log_file, 'a') as f:
                                f.write(f"{failed_trait['trait_name']}: {failed_trait['error']}\n")
            else:
                # 串行处理 | Serial processing
                self.logger.info("📝 使用串行处理 | Using serial processing")

                # 创建失败日志文件 | Create failed log file
                failed_log_file = output_base_dir / "failed_traits.log"
                failed_traits = []

                for i in trait_columns:  # 使用trait_columns而不是range
                    result = self.process_single_trait(pheno_source, vcf_file, i,
                                                       output_base_dir, global_config)
                    self.results.append(result)
                    if result['status'] == 'failed':
                        failed_traits.append(result)
                        with open(failed_log_file, 'a') as f:
                            f.write(f"{result['trait_name']}: {result['error']}\n")

            # 生成批量处理报告 | Generate batch processing report
            self._generate_batch_report(output_base_dir, len(trait_columns), failed_traits)

            return True

        except Exception as e:
            self.logger.error(f"❌ 批量处理失败 | Batch processing failed: {str(e)}")
            return False

    def _generate_batch_report(self, output_dir: Path, total_traits: int, failed_traits: List[Dict]):
        """生成批量处理报告 | Generate batch processing report"""
        try:
            report_file = output_dir / "batch_processing_report.txt"

            successful_traits = total_traits - len(failed_traits)
            success_rate = (successful_traits / total_traits) * 100 if total_traits > 0 else 0

            with open(report_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("                    TASSEL GWAS Batch Processing Report\n")
                f.write("=" * 80 + "\n")
                f.write(f"Processing Date:    {Path.cwd()}\n")
                f.write(f"Output Directory:  {output_dir}\n\n")

                f.write("PROCESSING SUMMARY\n")
                f.write("-" * 50 + "\n")
                f.write(f"Total Traits:      {total_traits}\n")
                f.write(f"Successful:         {successful_traits} ({success_rate:.1f}%)\n")
                f.write(f"Failed:            {len(failed_traits)} ({100-success_rate:.1f}%)\n\n")

                if failed_traits:
                    f.write("FAILED TRAITS\n")
                    f.write("-" * 50 + "\n")
                    for trait in failed_traits:
                        f.write(f"• {trait['trait_name']}: {trait['error']}\n")

                f.write("\nSUCCESSFUL TRAITS\n")
                f.write("-" * 50 + "\n")
                for result in self.results:
                    if result['status'] == 'success':
                        f.write(f"✓ {result['trait_name']}\n")
                        f.write(f"  Output Directory: {result['output_dir']}\n")

            self.logger.info(f"✅ 批量处理报告已生成 | Batch processing report generated: {report_file}")
            self.logger.info(f"   • 成功处理 | Successful: {successful_traits}/{total_traits} ({success_rate:.1f}%)")
            self.logger.info(f"   • 处理失败 | Failed: {len(failed_traits)}/{total_traits}")

        except Exception as e:
            self.logger.warning(f"⚠️ 批量处理报告生成失败 | Batch processing report generation failed: {e}")

    def get_summary(self) -> Dict:
        """获取处理摘要 | Get processing summary"""
        if not self.results:
            return {'total': 0, 'successful': 0, 'failed': 0, 'success_rate': 0.0}

        total = len(self.results)
        successful = len([r for r in self.results if r['status'] == 'success'])
        failed = total - successful
        success_rate = (successful / total) * 100 if total > 0 else 0.0

        return {
            'total': total,
            'successful': successful,
            'failed': failed,
            'success_rate': success_rate,
            'failed_traits': [r for r in self.results if r['status'] == 'failed']
        }