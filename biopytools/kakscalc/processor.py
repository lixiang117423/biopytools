"""
📈 Ka/Ks Calculator结果处理模块
功能: 解析和处理Ka/Ks计算结果 | Parse and process Ka/Ks calculation results
"""

import os
import json
import pandas as pd
from typing import Dict, Any
from datetime import datetime
from .config import KaKsConfig
from .logger import Logger

class ResultProcessor:
    """📈 Ka/Ks结果处理器 | Ka/Ks result processor"""
    
    def __init__(self, logger: Logger):
        """
        🏗️ 初始化结果处理器 | Initialize result processor
        
        Args:
            logger: 日志器实例 | Logger instance
        """
        self.logger = logger
        self.config = KaKsConfig()
    
    def parse_results(self, output_file: str) -> pd.DataFrame:
        """
        📊 解析Ka/Ks计算结果 | Parse Ka/Ks calculation results
        
        Args:
            output_file: KaKs_Calculator输出文件路径 | KaKs_Calculator output file path
            
        Returns:
            解析后的结果DataFrame | Parsed results DataFrame
        """
        try:
            self.logger.info("解析计算结果 | Parsing calculation results", "📊")
            
            # 🔍 检查文件内容 | Check file content
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                raise ValueError("输出文件为空 | Output file is empty")
            
            self.logger.debug(f"输出文件总行数 | Total lines: {len(lines)}")
            
            # 检查第一行是否为表头 | Check if first line is header
            first_line = lines[0].strip()
            self.logger.debug(f"第一行内容 | First line: '{first_line}'")
            
            has_header = False
            if 'Sequence' in first_line and 'Method' in first_line and 'Ka' in first_line:
                has_header = True
                self.logger.debug("检测到标准KaKs_Calculator输出格式（有表头） | Detected standard KaKs_Calculator format with header")
            else:
                self.logger.debug("检测到无表头格式 | Detected headerless format")
            
            # 📖 根据格式选择解析方式 | Choose parsing method based on format
            if has_header:
                # 标准格式：有表头的制表符分隔 | Standard format: tab-separated with header
                df = pd.read_csv(output_file, sep='\t')
                self.logger.debug(f"使用有表头模式解析：{df.shape[0]} 行 {df.shape[1]} 列 | Parsed with header: {df.shape[0]} rows {df.shape[1]} columns")
                
                # 标准化列名 | Standardize column names
                df.columns = df.columns.str.strip()
                
                # 确保标准列名存在 | Ensure standard column names exist
                if 'Ka/Ks' not in df.columns and 'Ka_Ks' in df.columns:
                    df['Ka/Ks'] = df['Ka_Ks']
                
                # 创建Pair_ID | Create Pair_ID
                if 'Sequence' in df.columns:
                    df['Pair_ID'] = df['Sequence'].astype(str)
                    # 如果需要，可以从Sequence列提取更多信息
                    df['Seq1_Name'] = df['Sequence'].astype(str)
                    df['Seq2_Name'] = df['Sequence'].astype(str)
                
            else:
                # 处理无表头格式（我们之前的情况）| Handle headerless format (our previous case)
                self.logger.debug("使用无表头模式解析 | Parsing without header")
                
                # 尝试多种分隔方式 | Try multiple separation methods
                df = None
                for sep_method in [r'\s+', '\t', ' ']:
                    try:
                        if sep_method == r'\s+':
                            df = pd.read_csv(output_file, sep=sep_method, header=None, engine='python')
                        else:
                            df = pd.read_csv(output_file, sep=sep_method, header=None)
                        
                        if df.shape[1] >= 10:  # 至少需要10列
                            self.logger.debug(f"分隔符 '{sep_method}' 成功：{df.shape[1]} 列 | Separator '{sep_method}' successful: {df.shape[1]} columns")
                            break
                    except Exception as e:
                        self.logger.debug(f"分隔符 '{sep_method}' 失败: {e} | Separator '{sep_method}' failed: {e}")
                        continue
                
                if df is None or df.shape[1] < 10:
                    raise ValueError("无法正确解析输出文件 | Cannot parse output file correctly")
                
                # 为无表头格式定义列名 | Define column names for headerless format
                if df.shape[1] >= 13:
                    # AXT格式混合结果：Index Seq1 Seq2 Start1 Length1 Strand1 Start2 Length2 Strand2 Method Ka Ks Ka/Ks ...
                    column_names = [
                        'Index', 'Seq1_Name', 'Seq2_Name', 'Start1', 'Length1', 'Strand1', 
                        'Start2', 'Length2', 'Strand2', 'Method', 'Ka', 'Ks', 'Ka/Ks'
                    ]
                    
                    # 添加额外列名 | Add extra column names
                    extra_cols = df.shape[1] - len(column_names)
                    if extra_cols > 0:
                        extra_names = [
                            'P_Value', 'Length_seq', 'S_Sites', 'N_Sites', 'Fold_Sites',
                            'Substitutions', 'S_Substitutions', 'N_Substitutions', 'Fold_S_Sub', 'Fold_N_Sub'
                        ]
                        for i in range(extra_cols):
                            if i < len(extra_names):
                                column_names.append(extra_names[i])
                            else:
                                column_names.append(f'Extra_Col_{i+1}')
                    
                    df.columns = column_names[:df.shape[1]]
                    
                    # 创建Pair_ID | Create Pair_ID
                    if 'Seq1_Name' in df.columns and 'Seq2_Name' in df.columns:
                        df['Pair_ID'] = df['Seq1_Name'].astype(str) + '_vs_' + df['Seq2_Name'].astype(str)
                    else:
                        df['Pair_ID'] = [f"pair_{i+1}" for i in range(len(df))]
                else:
                    # 列数不足，使用通用列名 | Insufficient columns, use generic names
                    df.columns = [f'Col_{i+1}' for i in range(df.shape[1])]
                    df['Pair_ID'] = [f"pair_{i+1}" for i in range(len(df))]
            
            self.logger.info(f"成功解析：{df.shape[0]} 行 {df.shape[1]} 列 | Successfully parsed: {df.shape[0]} rows {df.shape[1]} columns")
            self.logger.debug(f"列名 | Column names: {list(df.columns)}")
            
            # 🔧 数据类型转换 | Data type conversion
            self._convert_data_types(df)
            
            # 显示解析结果样本 | Show parsing sample
            if len(df) > 0:
                sample_data = df.iloc[0].to_dict()
                # 只显示关键列 | Only show key columns
                key_cols = ['Pair_ID', 'Method', 'Ka', 'Ks', 'Ka/Ks']
                sample_display = {k: v for k, v in sample_data.items() if k in key_cols and k in df.columns}
                self.logger.debug(f"解析样本 | Parsed sample: {sample_display}")
            
            # 🧮 计算附加指标 | Calculate additional metrics
            df = self._add_calculated_metrics(df)
            
            # 🎯 分类选择压力 | Classify selection pressure
            if 'Ka/Ks' in df.columns:
                df['Selection_Type'] = df['Ka/Ks'].apply(self._classify_selection)
                df['Selection_Strength'] = df['Ka/Ks'].apply(self._classify_selection_strength)
            
            # 📏 添加序列信息 | Add sequence information
            length_cols = ['Length', 'Length1', 'Length_seq']
            for col in length_cols:
                if col in df.columns:
                    df['Sequence_Length'] = pd.to_numeric(df[col], errors='coerce')
                    break
            
            # 🔍 质量控制标记 | Quality control flags
            df = self._add_quality_flags(df)
            
            self.logger.success(f"结果解析完成：{len(df)} 个序列对 | Results parsed: {len(df)} sequence pairs")
            return df
            
        except Exception as e:
            self.logger.error(f"结果解析失败 | Failed to parse results: {e}")
            raise
    
    def _convert_data_types(self, df: pd.DataFrame):
        """🔧 数据类型转换 | Convert data types"""
        try:
            # 数值列转换 | Convert numeric columns
            numeric_columns = ['Ka', 'Ks', 'Ka/Ks', 'Ka_Ks', 'P_Value', 'P-Value(Fisher)', 
                             'Length', 'Length1', 'S-Sites', 'N-Sites', 'S_Sites', 'N_Sites']
            
            for col in numeric_columns:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # 字符串列转换 | Convert string columns
            string_columns = ['Sequence', 'Method', 'Seq1_Name', 'Seq2_Name', 'Pair_ID']
            for col in string_columns:
                if col in df.columns:
                    df[col] = df[col].astype(str)
            
            # 确保Ka/Ks列存在 | Ensure Ka/Ks column exists
            if 'Ka/Ks' not in df.columns:
                if 'Ka_Ks' in df.columns:
                    df['Ka/Ks'] = df['Ka_Ks']
                elif 'Ka' in df.columns and 'Ks' in df.columns:
                    # 计算Ka/Ks | Calculate Ka/Ks
                    df['Ka/Ks'] = df['Ka'] / df['Ks']
                    df['Ka/Ks'] = df['Ka/Ks'].replace([float('inf'), float('-inf')], float('nan'))
            
            self.logger.debug("数据类型转换完成 | Data type conversion completed")
            
        except Exception as e:
            self.logger.warning(f"数据类型转换警告 | Data type conversion warning: {e}")
    
    def _add_calculated_metrics(self, df: pd.DataFrame) -> pd.DataFrame:
        """🧮 添加计算指标 | Add calculated metrics"""
        try:
            # 检查必需的列是否存在 | Check if required columns exist
            required_cols = ['Ka', 'Ks']
            omega_col = None
            
            # 寻找Ka/Ks列 | Find Ka/Ks column
            if 'Ka/Ks' in df.columns:
                omega_col = 'Ka/Ks'
                required_cols.append('Ka/Ks')
            elif 'Ka_Ks' in df.columns:
                omega_col = 'Ka_Ks'
                required_cols.append('Ka_Ks')
                # 创建标准的Ka/Ks列 | Create standard Ka/Ks column
                df['Ka/Ks'] = df['Ka_Ks']
            
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                self.logger.warning(f"缺少必需列 | Missing required columns: {missing_cols}")
                return df
            
            if not omega_col:
                self.logger.warning("未找到Ka/Ks列 | Ka/Ks column not found")
                return df
            
            # Ka/Ks比率别名 | Ka/Ks ratio alias
            df['Omega'] = df[omega_col]
            
            # 替换率差异 | Substitution rate difference
            df['Ka_minus_Ks'] = df['Ka'] - df['Ks']
            
            # 相对替换率 | Relative substitution rates
            ka_median = df['Ka'].median()
            ks_median = df['Ks'].median()
            
            df['Ka_relative'] = df['Ka'] / ka_median if ka_median > 0 else 0
            df['Ks_relative'] = df['Ks'] / ks_median if ks_median > 0 else 0
            
            # 标准化得分 | Normalized scores
            ka_std = df['Ka'].std()
            ks_std = df['Ks'].std()
            omega_std = df[omega_col].std()
            
            if ka_std > 0:
                df['Ka_zscore'] = (df['Ka'] - df['Ka'].mean()) / ka_std
            else:
                df['Ka_zscore'] = 0
                
            if ks_std > 0:
                df['Ks_zscore'] = (df['Ks'] - df['Ks'].mean()) / ks_std
            else:
                df['Ks_zscore'] = 0
                
            if omega_std > 0:
                df['Omega_zscore'] = (df[omega_col] - df[omega_col].mean()) / omega_std
            else:
                df['Omega_zscore'] = 0
            
            self.logger.debug("计算指标添加完成 | Calculated metrics added successfully")
            return df
            
        except Exception as e:
            self.logger.warning(f"添加计算指标时出错 | Error adding calculated metrics: {e}")
            return df
    
    def _classify_selection(self, omega: float) -> str:
        """🎯 根据omega值分类选择类型 | Classify selection type based on omega value"""
        if pd.isna(omega):
            return "Unknown"
        
        thresholds = self.config.SELECTION_THRESHOLDS
        
        if omega < thresholds['strong_negative']:
            return "Strong_Negative"
        elif omega < thresholds['moderate_negative']:
            return "Moderate_Negative"
        elif thresholds['neutral_lower'] <= omega <= thresholds['neutral_upper']:
            return "Neutral"
        elif omega <= thresholds['weak_positive']:
            return "Weak_Positive"
        else:
            return "Strong_Positive"
    
    def _classify_selection_strength(self, omega: float) -> str:
        """💪 分类选择强度 | Classify selection strength"""
        if pd.isna(omega):
            return "Unknown"
        
        if omega < 0.1:
            return "Very_Strong_Negative"
        elif omega < 0.5:
            return "Strong_Negative"
        elif omega < 1.0:
            return "Negative"
        elif omega < 1.1:
            return "Nearly_Neutral"
        elif omega < 2.0:
            return "Positive"
        else:
            return "Strong_Positive"
    
    def _add_quality_flags(self, df: pd.DataFrame) -> pd.DataFrame:
        """🚩 添加质量控制标记 | Add quality control flags"""
        try:
            # 初始化质量标记 | Initialize quality flags
            df['Quality_Flag'] = "Good"
            df['Quality_Issues'] = ""
            
            # 检查异常值 | Check outliers
            omega_q99 = df['Ka/Ks'].quantile(0.99)
            ks_q01 = df['Ks'].quantile(0.01)
            
            # 标记异常值 | Flag outliers
            mask_high_omega = df['Ka/Ks'] > omega_q99
            mask_low_ks = df['Ks'] < ks_q01
            mask_na_values = df[['Ka', 'Ks', 'Ka/Ks']].isna().any(axis=1)
            
            df.loc[mask_high_omega, 'Quality_Flag'] = "Warning"
            df.loc[mask_high_omega, 'Quality_Issues'] += "High_Omega;"
            
            df.loc[mask_low_ks, 'Quality_Flag'] = "Warning"
            df.loc[mask_low_ks, 'Quality_Issues'] += "Low_Ks;"
            
            df.loc[mask_na_values, 'Quality_Flag'] = "Poor"
            df.loc[mask_na_values, 'Quality_Issues'] += "Missing_Values;"
            
            # 清理质量问题字符串 | Clean quality issues string
            df['Quality_Issues'] = df['Quality_Issues'].str.rstrip(';')
            
            return df
            
        except Exception as e:
            self.logger.warning(f"添加质量标记时出错 | Error adding quality flags: {e}")
            return df
    
    def generate_summary_stats(self, df: pd.DataFrame) -> Dict[str, Any]:
        """📊 生成汇总统计 | Generate summary statistics"""
        try:
            self.logger.info("生成汇总统计 | Generating summary statistics", "📊")
            
            # 基本统计 | Basic statistics
            total_pairs = len(df)
            successful_calcs = len(df.dropna(subset=['Ka/Ks']))
            failed_calcs = total_pairs - successful_calcs
            
            # 计算统计值 | Calculate statistics
            stats = {
                'analysis_info': {
                    'timestamp': datetime.now().isoformat(),
                    'total_pairs': total_pairs,
                    'successful_calculations': successful_calcs,
                    'failed_calculations': failed_calcs,
                    'success_rate': successful_calcs / total_pairs if total_pairs > 0 else 0
                },
                'ka_statistics': self._calculate_variable_stats(df, 'Ka'),
                'ks_statistics': self._calculate_variable_stats(df, 'Ks'),
                'omega_statistics': self._calculate_variable_stats(df, 'Ka/Ks'),
                'selection_distribution': df['Selection_Type'].value_counts().to_dict(),
                'selection_strength_distribution': df['Selection_Strength'].value_counts().to_dict(),
                'quality_distribution': df['Quality_Flag'].value_counts().to_dict()
            }
            
            # 添加生物学解释 | Add biological interpretation
            stats['biological_interpretation'] = self._interpret_results(stats)
            
            self.logger.success("汇总统计生成完成 | Summary statistics generated")
            return stats
            
        except Exception as e:
            self.logger.error(f"生成汇总统计失败 | Failed to generate summary statistics: {e}")
            return {}
    
    def _calculate_variable_stats(self, df: pd.DataFrame, column: str) -> Dict[str, float]:
        """📏 计算变量统计值 | Calculate variable statistics"""
        try:
            data = df[column].dropna()
            if len(data) == 0:
                return {'count': 0}
            
            return {
                'count': len(data),
                'mean': float(data.mean()),
                'median': float(data.median()),
                'std': float(data.std()),
                'min': float(data.min()),
                'max': float(data.max()),
                'q25': float(data.quantile(0.25)),
                'q75': float(data.quantile(0.75)),
                'skewness': float(data.skew()) if len(data) > 1 else 0,
                'kurtosis': float(data.kurtosis()) if len(data) > 1 else 0
            }
        except Exception:
            return {'count': 0, 'error': 'calculation_failed'}
    
    def _interpret_results(self, stats: Dict[str, Any]) -> Dict[str, str]:
        """🧠 生物学结果解释 | Biological result interpretation"""
        interpretation = {}
        
        try:
            # 整体选择压力解释 | Overall selection pressure interpretation
            omega_mean = stats.get('omega_statistics', {}).get('mean', 0)
            
            if omega_mean < 0.5:
                interpretation['overall_selection'] = "🔒 强负选择：基因受到强烈功能约束 | Strong negative selection: genes under strong functional constraints"
            elif omega_mean < 1.0:
                interpretation['overall_selection'] = "⚖️ 净化选择：有害突变被清除 | Purifying selection: deleterious mutations being removed"
            elif omega_mean > 1.5:
                interpretation['overall_selection'] = "🚀 正选择：适应性进化 | Positive selection: adaptive evolution"
            else:
                interpretation['overall_selection'] = "🎯 近中性进化：接近中性突变 | Nearly neutral evolution: close to neutral mutation"
            
            # 选择类型分布解释 | Selection type distribution interpretation
            selection_dist = stats.get('selection_distribution', {})
            negative_ratio = (selection_dist.get('Strong_Negative', 0) + selection_dist.get('Moderate_Negative', 0)) / stats['analysis_info']['successful_calculations']
            
            if negative_ratio > 0.8:
                interpretation['selection_pattern'] = "🛡️ 高度保守：大多数基因受到强选择约束 | Highly conserved: most genes under strong selective constraints"
            elif negative_ratio > 0.5:
                interpretation['selection_pattern'] = "⚖️ 平衡选择：负选择为主，部分基因快速进化 | Balanced selection: predominantly negative with some rapidly evolving genes"
            else:
                interpretation['selection_pattern'] = "🌪️ 多样化选择：选择压力多样化 | Diversifying selection: diverse selective pressures"
            
            return interpretation
            
        except Exception as e:
            self.logger.warning(f"结果解释生成失败 | Failed to generate interpretation: {e}")
            return {'error': 'interpretation_failed'}
    
    def save_results(self, df: pd.DataFrame, stats: Dict[str, Any], output_dir: str):
        """💾 保存处理后的结果 | Save processed results"""
        try:
            self.logger.info("保存分析结果 | Saving analysis results", "💾")
            
            # 📊 保存详细结果(TSV格式) | Save detailed results (TSV format)
            detailed_file = os.path.join(output_dir, self.config.OUTPUT_FILES['detailed'])
            df.to_csv(detailed_file, sep='\t', index=False, encoding='utf-8')
            self.logger.debug(f"TSV详细结果已保存 | TSV detailed results saved: {detailed_file}")
            
            # 📊 保存详细结果(CSV格式) | Save detailed results (CSV format)
            csv_file = os.path.join(output_dir, self.config.OUTPUT_FILES['detailed_csv'])
            df.to_csv(csv_file, index=False, encoding='utf-8')
            self.logger.debug(f"CSV详细结果已保存 | CSV detailed results saved: {csv_file}")
            
            # 📈 保存Excel汇总文件 | Save Excel summary file
            summary_file = os.path.join(output_dir, self.config.OUTPUT_FILES['summary'])
            self._save_excel_summary(df, stats, summary_file)
            
            # 💾 保存JSON统计文件 | Save JSON statistics file
            stats_file = os.path.join(output_dir, self.config.OUTPUT_FILES['statistics'])
            with open(stats_file, 'w', encoding='utf-8') as f:
                json.dump(stats, f, indent=2, ensure_ascii=False, default=str)
            self.logger.debug(f"统计数据已保存 | Statistics saved: {stats_file}")
            
            self.logger.success(f"所有结果已保存到 | All results saved to: {output_dir}")
            
            # 📋 列出生成的文件 | List generated files
            generated_files = [
                self.config.OUTPUT_FILES['detailed'],
                self.config.OUTPUT_FILES['detailed_csv'],
                self.config.OUTPUT_FILES['summary'],
                self.config.OUTPUT_FILES['statistics']
            ]
            
            self.logger.info("生成的输出文件 | Generated output files:")
            for file_name in generated_files:
                file_path = os.path.join(output_dir, file_name)
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    self.logger.info(f"  📄 {file_name} ({file_size} bytes)")
                    
        except Exception as e:
            self.logger.error(f"保存结果失败 | Failed to save results: {e}")
            raise
    
    def _save_excel_summary(self, df: pd.DataFrame, stats: Dict[str, Any], summary_file: str):
        """📈 保存Excel汇总文件 | Save Excel summary file"""
        try:
            # 🧹 清理数据框用于Excel输出 | Clean dataframe for Excel output
            df_clean = df.copy()
            
            # 替换可能导致Excel问题的值 | Replace values that might cause Excel issues
            df_clean = df_clean.replace([float('inf'), float('-inf')], ['INF', '-INF'])
            df_clean = df_clean.fillna('')  # 将NaN替换为空字符串
            
            # 确保所有列名都是字符串 | Ensure all column names are strings
            df_clean.columns = [str(col) for col in df_clean.columns]
            
            with pd.ExcelWriter(summary_file, engine='openpyxl') as writer:
                # 主要结果 | Main results
                df_clean.to_excel(writer, sheet_name='Results', index=False)
                
                # 统计摘要 | Statistical summary
                if stats:
                    # 基本信息 | Basic info
                    analysis_info = stats.get('analysis_info', {})
                    if analysis_info:
                        basic_info = pd.DataFrame([analysis_info])
                        basic_info.to_excel(writer, sheet_name='Summary_Info', index=False)
                    
                    # 选择分布 | Selection distribution
                    selection_dist = stats.get('selection_distribution', {})
                    if selection_dist:
                        selection_df = pd.DataFrame(list(selection_dist.items()),
                                                  columns=['Selection_Type', 'Count'])
                        if len(selection_df) > 0:
                            total_count = selection_df['Count'].sum()
                            selection_df['Percentage'] = (selection_df['Count'] / total_count * 100).round(2)
                            selection_df.to_excel(writer, sheet_name='Selection_Distribution', index=False)
                    
                    # 质量分布 | Quality distribution
                    quality_dist = stats.get('quality_distribution', {})
                    if quality_dist:
                        quality_df = pd.DataFrame(list(quality_dist.items()),
                                                columns=['Quality_Flag', 'Count'])
                        if len(quality_df) > 0:
                            total_count = quality_df['Count'].sum()
                            quality_df['Percentage'] = (quality_df['Count'] / total_count * 100).round(2)
                            quality_df.to_excel(writer, sheet_name='Quality_Distribution', index=False)
                    
                    # 生物学解释 | Biological interpretation
                    bio_interp = stats.get('biological_interpretation', {})
                    if bio_interp:
                        interp_df = pd.DataFrame(list(bio_interp.items()),
                                               columns=['Aspect', 'Interpretation'])
                        interp_df.to_excel(writer, sheet_name='Interpretation', index=False)
                    
                    # Ka/Ks统计 | Ka/Ks statistics
                    omega_stats = stats.get('omega_statistics', {})
                    if omega_stats:
                        omega_df = pd.DataFrame([omega_stats])
                        omega_df.to_excel(writer, sheet_name='KaKs_Statistics', index=False)
            
            self.logger.debug(f"Excel汇总文件已保存 | Excel summary saved: {summary_file}")
            
        except Exception as e:
            self.logger.warning(f"Excel文件保存失败，尝试保存为CSV | Failed to save Excel file, trying CSV: {e}")
            # 如果Excel保存失败，保存为CSV备份 | If Excel save fails, save as CSV backup
            csv_backup = summary_file.replace('.xlsx', '_backup.csv')
            df.to_csv(csv_backup, index=False, encoding='utf-8')
            self.logger.info(f"已保存CSV备份文件 | Saved CSV backup: {csv_backup}")
