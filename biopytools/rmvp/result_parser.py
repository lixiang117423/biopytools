"""
rMVP结果解析和整合|rMVP Result Parsing and Integration
"""

import os
import re
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd


class RMVPResultParser:
    """rMVP结果解析器|rMVP Result Parser"""

    def __init__(self, output_dir: Path, output_prefix: str, logger):
        """
        初始化|Initialize

        Args:
            output_dir: 输出目录|Output directory
            output_prefix: 输出前缀|Output prefix
            logger: 日志对象|Logger object
        """
        self.output_dir = Path(output_dir)
        self.output_prefix = output_prefix
        self.logger = logger

    def parse_and_integrate_results(self, trait_names: List[str],
                                    models: List[str]) -> Dict[str, pd.DataFrame]:
        """
        解析并整合所有结果|Parse and integrate all results

        Args:
            trait_names: 表型名称列表|List of trait names
            models: 模型列表|List of models

        Returns:
            字典，包含整合后的结果|Dictionary containing integrated results
            {
                "glm": DataFrame,  # GLM模型所有表型的结果
                "mlm": DataFrame,  # MLM模型所有表型的结果
                "farmcpu": DataFrame,  # FarmCPU模型所有表型的结果
                "all": DataFrame  # 所有模型所有表型的结果
            }
        """
        self.logger.info(" 开始解析和整合结果|Starting to parse and integrate results")

        integrated_results = {}

        # 为每个模型整合结果|Integrate results for each model
        for model in models:
            model_lower = model.lower()
            self.logger.info(f"   处理{model}模型|Processing {model} model")

            # 查找该模型所有表型的结果文件|Find result files for all traits of this model
            result_dfs = []

            for trait in trait_names:
                # rMVP输出文件命名格式|rMVP output file naming format
                # {output_prefix}_{trait}.{model}.pmap
                pmap_file = self.output_dir / f"{self.output_prefix}_{trait}.{model_lower}.pmap"

                if pmap_file.exists():
                    try:
                        # 读取结果文件|Read result file
                        df = pd.read_csv(pmap_file, sep='\t')
                        # 重命名P值列|Rename P-value column
                        if 'p_value' in df.columns:
                            df = df.rename(columns={'p_value': f"{trait}_{model}"})
                            # 只保留SNP、Chromosome、Position和P值|Keep only SNP, Chr, Pos, and P-value
                            df_result = df[['SNP', 'Chr', 'Pos', f"{trait}_{model}"]]
                            result_dfs.append(df_result)
                        else:
                            self.logger.warning(f"     文件格式异常|File format abnormal: {pmap_file}")
                    except Exception as e:
                        self.logger.warning(f"     读取文件失败|Failed to read file {pmap_file}: {e}")
                else:
                    self.logger.warning(f"     文件不存在|File does not exist: {pmap_file}")

            # 合并该模型所有表型的结果|Merge all traits results for this model
            if result_dfs:
                # 按SNP合并|Merge by SNP
                merged_df = result_dfs[0]
                for df in result_dfs[1:]:
                    merged_df = pd.merge(merged_df, df, on=['SNP', 'Chr', 'Pos'], how='outer')

                integrated_results[model_lower] = merged_df
                self.logger.info(f"     {model}模型结果整合完成|{model} model results integrated: {len(merged_df)} SNPs")
            else:
                self.logger.warning(f"     {model}模型没有找到结果文件|No result files found for {model} model")

        # 整合所有模型到一个文件|Integrate all models into one file
        if integrated_results:
            all_models_df = integrated_results[list(integrated_results.keys())[0]]

            for model in list(integrated_results.keys())[1:]:
                all_models_df = pd.merge(
                    all_models_df,
                    integrated_results[model],
                    on=['SNP', 'Chr', 'Pos'],
                    how='outer',
                    suffixes=('', f'_{model}')
                )

            integrated_results['all'] = all_models_df
            self.logger.info(f" 所有模型整合完成|All models integrated: {len(all_models_df)} SNPs")

        return integrated_results

    def save_integrated_results(self, integrated_results: Dict[str, pd.DataFrame]) -> List[Path]:
        """
        保存整合后的结果|Save integrated results

        Args:
            integrated_results: 整合后的结果|Integrated results

        Returns:
            保存的文件列表|List of saved files
        """
        self.logger.info(" 保存整合结果|Saving integrated results")

        saved_files = []

        # 保存每个模型的结果|Save results for each model
        for model_name, df in integrated_results.items():
            if model_name == 'all':
                output_file = self.output_dir / f"{self.output_prefix}_all_models_integrated.txt"
            else:
                output_file = self.output_dir / f"{self.output_prefix}_{model_name}_integrated.txt"

            try:
                # 保存为制表符分隔的文本文件|Save as tab-separated text file
                df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
                self.logger.info(f"   已保存|Saved: {output_file} ({len(df)} SNPs)")
                saved_files.append(output_file)
            except Exception as e:
                self.logger.error(f"   保存失败|Failed to save {output_file}: {e}")

        return saved_files

    def collect_output_files(self, trait_names: List[str],
                            models: List[str]) -> Dict[str, List[Path]]:
        """
        收集所有输出文件|Collect all output files

        Args:
            trait_names: 表型名称列表|List of trait names
            models: 模型列表|List of models

        Returns:
            文件分类字典|File classification dictionary
            {
                "tables": [路径列表],  # 表格文件
                "figures": [路径列表],  # 图片文件
                "logs": [路径列表]      # 日志文件
            }
        """
        self.logger.info(" 收集输出文件|Collecting output files")

        files = {
            "tables": [],
            "figures": [],
            "logs": []
        }

        # 收集每个表型每个模型的文件|Collect files for each trait and model
        for trait in trait_names:
            for model in models:
                model_lower = model.lower()
                trait_prefix = f"{self.output_prefix}_{trait}"

                # 表格文件|Table files
                table_patterns = [
                    f"{trait_prefix}.{model_lower}.pmap",  # 所有SNP的P值
                    f"{trait_prefix}.{model_lower}.pmap.signal",  # 显著SNP
                    f"{trait_prefix}.{model_lower}.eff",  # SNP效应
                    f"{trait_prefix}.RData"  # R结果对象
                ]

                for pattern in table_patterns:
                    matching_files = list(self.output_dir.glob(pattern))
                    files["tables"].extend(matching_files)

                # 图片文件|Figure files
                figure_patterns = [
                    f"{trait_prefix}.{model_lower}.Manhattan*.{self._get_file_extension()}",
                    f"{trait_prefix}.{model_lower}.QQ*.{self._get_file_extension()}",
                    f"{trait_prefix}.Phe_Distribution.{self._get_file_extension()}",
                    f"{trait_prefix}.PCA*.{self._get_file_extension()}"
                ]

                for pattern in figure_patterns:
                    matching_files = list(self.output_dir.glob(pattern))
                    files["figures"].extend(matching_files)

        # 日志文件|Log files
        log_patterns = [
            f"{self.output_prefix}*.log",
            "*.RData"
        ]

        for pattern in log_patterns:
            matching_files = list(self.output_dir.glob(pattern))
            files["logs"].extend(matching_files)

        # 去重|Remove duplicates
        for key in files:
            files[key] = list(set(files[key]))

        # 输出统计|Output statistics
        self.logger.info("   文件统计|File statistics:")
        self.logger.info(f"     表格文件|Table files: {len(files['tables'])}")
        self.logger.info(f"     图片文件|Figure files: {len(files['figures'])}")
        self.logger.info(f"     日志文件|Log files: {len(files['logs'])}")

        return files

    def generate_summary_report(self, trait_names: List[str],
                               models: List[str]) -> Path:
        """
        生成汇总报告|Generate summary report

        Args:
            trait_names: 表型名称列表|List of trait names
            models: 模型列表|List of models

        Returns:
            报告文件路径|Report file path
        """
        self.logger.info(" 生成汇总报告|Generating summary report")

        report_file = self.output_dir / f"{self.output_prefix}_summary_report.txt"

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                # 报告头部|Report header
                f.write("=" * 80 + "\n")
                f.write("                    rMVP GWAS分析汇总报告|rMVP GWAS Summary Report\n")
                f.write("=" * 80 + "\n\n")

                # 基本信息|Basic information
                f.write("分析配置|Analysis Configuration\n")
                f.write("-" * 50 + "\n")
                f.write(f"输出目录|Output directory: {self.output_dir}\n")
                f.write(f"输出前缀|Output prefix: {self.output_prefix}\n")
                f.write(f"表型数量|Number of traits: {len(trait_names)}\n")
                f.write(f"分析模型|Models: {', '.join(models)}\n\n")

                # 表型列表|Trait list
                f.write("表型列表|Trait List\n")
                f.write("-" * 50 + "\n")
                for i, trait in enumerate(trait_names, 1):
                    f.write(f"  {i}. {trait}\n")
                f.write("\n")

                # 文件统计|File statistics
                files = self.collect_output_files(trait_names, models)
                f.write("输出文件统计|Output File Statistics\n")
                f.write("-" * 50 + "\n")
                f.write(f"表格文件|Table files: {len(files['tables'])}\n")
                f.write(f"图片文件|Figure files: {len(files['figures'])}\n")
                f.write(f"日志文件|Log files: {len(files['logs'])}\n\n")

                # 整合结果文件|Integrated result files
                f.write("整合结果文件|Integrated Result Files\n")
                f.write("-" * 50 + "\n")

                for model in models:
                    model_lower = model.lower()
                    integrated_file = self.output_dir / f"{self.output_prefix}_{model_lower}_integrated.txt"
                    if integrated_file.exists():
                        file_size = integrated_file.stat().st_size
                        f.write(f"  ✓ {model}模型|{model} model: {integrated_file.name} ({file_size:,} bytes)\n")

                all_integrated_file = self.output_dir / f"{self.output_prefix}_all_models_integrated.txt"
                if all_integrated_file.exists():
                    file_size = all_integrated_file.stat().st_size
                    f.write(f"  ✓ 所有模型|All models: {all_integrated_file.name} ({file_size:,} bytes)\n")

            self.logger.info(f"   汇总报告已生成|Summary report generated: {report_file}")
            return report_file

        except Exception as e:
            self.logger.error(f"   生成汇总报告失败|Failed to generate summary report: {e}")
            return Path("")

    def _get_file_extension(self) -> str:
        """
        获取文件扩展名（从配置中）|Get file extension (from config)

        Returns:
            文件扩展名|File extension
        """
        # 这个方法需要访问config，暂时返回jpg
        # TODO: 从config中读取
        return "jpg"
