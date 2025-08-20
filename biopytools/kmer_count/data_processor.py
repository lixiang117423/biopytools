"""
📊 数据处理模块 | Data Processing Module
"""

import os
import pandas as pd
from pathlib import Path
from typing import Dict, List

class DataProcessor:
    """📊 数据处理器 | Data Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.kmer_ids = {}
        self.bed_info = None
    
    # def parse_kmer_library(self) -> Dict[str, str]:
    #     """🧬 解析k-mer库 | Parse k-mer library"""
    #     self.logger.info("🧬 解析k-mer库 | Parsing k-mer library")
        
    #     kmer_dict = {}
    #     current_id = None
        
    #     with open(self.config.kmer_lib, 'r') as f:
    #         for line in f:
    #             line = line.strip()
    #             if line.startswith('>'):
    #                 current_id = line[1:]  # 去掉'>' | Remove '>'
    #             elif line and current_id:
    #                 kmer_seq = line.upper()
    #                 kmer_dict[kmer_seq] = current_id
    #                 self.kmer_ids[kmer_seq] = current_id
        
    #     self.logger.info(f"✅ 解析完成，共{len(kmer_dict)}个k-mer | Parsing completed, {len(kmer_dict)} k-mers total")
    #     return kmer_dict
    def parse_kmer_library(self) -> Dict[str, str]:
        """🧬 解析k-mer库 | Parse k-mer library"""
        self.logger.info("🧬 解析k-mer库 | Parsing k-mer library")
        
        kmer_dict = {}
        self.kmer_pairs = []  # 🔥 新增：保存所有ID-序列对
        current_id = None
        
        with open(self.config.kmer_lib, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]  # 去掉'>' | Remove '>'
                elif line and current_id:
                    kmer_seq = line.upper()
                    kmer_dict[kmer_seq] = current_id
                    # 🔥 保存每个ID-序列对
                    self.kmer_pairs.append({
                        'kmer_id': current_id,
                        'kmer': kmer_seq,
                        'unique_key': f"{current_id}_{kmer_seq}"
                    })
                    self.kmer_ids[kmer_seq] = current_id  # 保持原有逻辑
        
        self.logger.info(f"✅ 解析完成，共{len(self.kmer_pairs)}个k-mer对 | Parsing completed, {len(self.kmer_pairs)} k-mer pairs total")
        return kmer_dict
    
    def parse_bed_file(self) -> pd.DataFrame:
        """📋 解析BED文件 | Parse BED file"""
        if not self.config.bed_file:
            return None
        
        self.logger.info("📋 解析BED文件 | Parsing BED file")
        
        # 🔥 只读取前6列，忽略多余的列 | Only read first 6 columns, ignore extra columns
        try:
            bed_df = pd.read_csv(
                self.config.bed_file, 
                sep='\t', 
                header=None,
                usecols=[0, 1, 2, 3, 4, 5],  # 只读取前6列
                names=['chr', 'start', 'end', 'kmer', 'score', 'strand'],  # 直接指定列名
                dtype={
                    'chr': str,
                    'start': int, 
                    'end': int,
                    'kmer': str,
                    'score': str,  # score可能是数字或'.'
                    'strand': str
                }
            )
            self.logger.info(f"✅ 成功读取BED文件前6列 | Successfully read first 6 columns of BED file")
            
        except Exception as e:
            self.logger.warning(f"⚠️ 读取前6列失败，尝试灵活读取 | Failed to read 6 columns, trying flexible reading: {e}")
            
            # 🔥 备用方案：灵活读取，跳过问题行 | Fallback: flexible reading, skip problematic lines
            try:
                bed_df = pd.read_csv(
                    self.config.bed_file, 
                    sep='\t', 
                    header=None,
                    on_bad_lines='skip'  # 跳过问题行
                )
                
                # 只保留前6列
                if bed_df.shape[1] >= 6:
                    bed_df = bed_df.iloc[:, :6]
                    bed_df.columns = ['chr', 'start', 'end', 'kmer', 'score', 'strand']
                elif bed_df.shape[1] >= 4:
                    # 如果只有4列，添加默认的score和strand
                    bed_df = bed_df.iloc[:, :4]
                    bed_df.columns = ['chr', 'start', 'end', 'kmer']
                    bed_df['score'] = '0'
                    bed_df['strand'] = '.'
                else:
                    raise ValueError(f"BED文件列数不足，至少需要4列，实际有{bed_df.shape[1]}列")
                
                self.logger.info(f"✅ 灵活读取成功，跳过了问题行 | Flexible reading successful, skipped problematic lines")
                
            except Exception as e2:
                raise ValueError(f"❌ BED文件格式错误，无法解析 | BED file format error, cannot parse: {e2}")
        
        # 🔥 确保start和end列为整数类型 | Ensure start and end columns are integer type
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
        
        # 转换k-mer序列为大写 | Convert k-mer sequences to uppercase
        bed_df['kmer'] = bed_df['kmer'].str.upper()
        
        self.bed_info = bed_df
        self.logger.info(f"✅ BED文件解析完成，共{len(bed_df)}条记录 | BED file parsing completed, {len(bed_df)} records total")
        return bed_df
    
    # def parse_bed_file(self) -> pd.DataFrame:
    #     """📋 解析BED文件 | Parse BED file"""
    #     if not self.config.bed_file:
    #         return None
        
    #     self.logger.info("📋 解析BED文件 | Parsing BED file")
        
    #     # 读取BED文件 | Read BED file
    #     bed_df = pd.read_csv(self.config.bed_file, sep='\t', header=None)
        
    #     # 设置列名 | Set column names
    #     # col_names = ['chr', 'start', 'end', 'kmer']
    #     # if bed_df.shape[1] > 4:
    #     #     col_names.extend([f'tmp_{i}' for i in range(1, bed_df.shape[1] - 3)])
    #     # 设置列名 | Set column names
    #     col_names = ['chr', 'start', 'end', 'kmer']
    #     if bed_df.shape[1] > 4:
    #         if bed_df.shape[1] == 6:
    #             col_names.extend(['score', 'strand'])
    #         else:
    #             col_names.extend([f'tmp_{i}' for i in range(1, bed_df.shape[1] - 3)])
        
    #     bed_df.columns = col_names[:bed_df.shape[1]]
        
    #     # 转换k-mer序列为大写 | Convert k-mer sequences to uppercase
    #     bed_df['kmer'] = bed_df['kmer'].str.upper()
        
    #     self.bed_info = bed_df
    #     self.logger.info(f"✅ BED文件解析完成，共{len(bed_df)}条记录 | BED file parsing completed, {len(bed_df)} records total")
    #     return bed_df
    
    # def parse_jellyfish_output(self, count_file: Path, sample_name: str) -> pd.DataFrame:
    #     """📊 解析Jellyfish输出 | Parse Jellyfish output"""
    #     self.logger.info(f"📊 解析Jellyfish输出 | Parsing Jellyfish output for {sample_name}")
        
    #     results = []
    #     with open(count_file, 'r') as f:
    #         for line in f:
    #             line = line.strip()
    #             if line:
    #                 parts = line.split()
    #                 if len(parts) == 2:
    #                     kmer_seq, count = parts
    #                     results.append({
    #                         'kmer': kmer_seq.upper(),
    #                         sample_name: int(count)
    #                     })
        
    #     df = pd.DataFrame(results)
    #     self.logger.info(f"✅ 解析完成，共{len(df)}个k-mer结果 | Parsing completed, {len(df)} k-mer results")
    #     return df

    def parse_jellyfish_output(self, count_file: Path, sample_name: str) -> pd.DataFrame:
        """📊 解析Jellyfish输出 | Parse Jellyfish output"""
        self.logger.info(f"📊 解析Jellyfish输出 | Parsing Jellyfish output for {sample_name}")
        
        kmer_counts = {}  # 使用字典去重，累加丰度
        
        with open(count_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) == 2:
                        kmer_seq, count = parts
                        kmer_seq = kmer_seq.upper()
                        count = int(count)
                        
                        # 累加相同k-mer的丰度
                        if kmer_seq in kmer_counts:
                            kmer_counts[kmer_seq] += count
                        else:
                            kmer_counts[kmer_seq] = count
        
        # 转换为DataFrame
        results = [{'kmer': kmer, sample_name: count} for kmer, count in kmer_counts.items()]
        df = pd.DataFrame(results)
        
        self.logger.info(f"✅ 解析完成，去重后共{len(df)}个k-mer结果 | Parsing completed, {len(df)} unique k-mer results after deduplication")
        return df
    
    # 单个样品合并bed文件
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """🔗 单个样品与BED文件合并 | Merge single sample with BED file"""
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并 | Merging sample {sample_name} with BED file")
        
    #     if self.bed_info is not None:
    #         # 创建unique key：chr+start+end+kmer
    #         self.bed_info['unique_key'] = self.bed_info['chr'] + '_' + self.bed_info['start'].astype(str) + '_' + self.bed_info['end'].astype(str) + '_' + self.bed_info['kmer']
    #         sample_df['unique_key'] = sample_df['kmer']  # 这里需要后续处理
            
    #         # 与BED信息合并，使用kmer合并但保留unique_key
    #         result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
    #         # 填充缺失值为0
    #         result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
    #         # 删除丰度为0的行
    #         result_df = result_df[result_df[sample_name] > 0]
    #         result_df = result_df.sort_values(['chr', 'start'])
            
    #         self.logger.info(f"🗑️ 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f"✅ 样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
        """🔗 单个样品与BED文件合并 | Merge single sample with BED file"""
        self.logger.info(f"🔗 样品{sample_name}与BED文件合并 | Merging sample {sample_name} with BED file")
        
        if self.bed_info is not None:
            # 🔥 不要在这里创建unique_key，直接用kmer合并
            result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
            # 填充缺失值为0
            result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
            # 删除丰度为0的行
            result_df = result_df[result_df[sample_name] > 0]
            
            # 🔥 在这里创建unique_key用于后续合并
            result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            
            result_df = result_df.sort_values(['chr', 'start'])
            
            self.logger.info(f"🗑️ 已删除丰度为0的行，剩余{len(result_df)}行")
        else:
            # 添加k-mer ID信息
            sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
            # 创建unique key：kmer_id+kmer
            sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
            # 🔥 添加：删除丰度为0的行 | Add: Remove rows with abundance = 0
            sample_df = sample_df[sample_df[sample_name] > 0]

            cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
            result_df = sample_df[cols]
        
        self.logger.info(f"✅ 样品{sample_name}合并完成，共{len(result_df)}行")
        return result_df
    
    def merge_sample_results(self, sample_results: List[pd.DataFrame]) -> pd.DataFrame:
        """🔗 合并所有样本结果 | Merge all sample results"""
        self.logger.info("🔗 合并样本结果 | Merging sample results")
        
        # 从第一个样本开始 | Start with first sample
        merged_df = sample_results[0].copy()
        
        # 逐个合并其他样本 | Merge other samples one by one
        for df in sample_results[1:]:
            merged_df = pd.merge(merged_df, df, on='kmer', how='outer')
        
        # 填充缺失值为0 | Fill missing values with 0
        sample_columns = [col for col in merged_df.columns if col != 'kmer']
        merged_df[sample_columns] = merged_df[sample_columns].fillna(0).astype(int)
        
        self.logger.info(f"✅ 合并完成，共{len(merged_df)}个k-mer，{len(sample_columns)}个样本")
        return merged_df
    
    def add_annotation_info(self, merged_df: pd.DataFrame) -> pd.DataFrame:
        """📝 添加注释信息 | Add annotation information"""
        if self.bed_info is not None:
            # 使用BED信息 | Use BED information
            self.logger.info("📋 添加BED注释信息 | Adding BED annotation information")
            
            # 合并BED信息 | Merge BED information
            result_df = pd.merge(self.bed_info, merged_df, on='kmer', how='left')
            
            # 填充缺失的丰度值为0 | Fill missing abundance values with 0
            sample_columns = [col for col in merged_df.columns if col != 'kmer']
            result_df[sample_columns] = result_df[sample_columns].fillna(0).astype(int)
            
            # 按照染色体和位置排序 | Sort by chromosome and position
            result_df = result_df.sort_values(['chr', 'start'])
            
        else:
            # 使用k-mer ID信息 | Use k-mer ID information
            self.logger.info("🧬 添加k-mer ID信息 | Adding k-mer ID information")
            
            # 添加k-mer ID列 | Add k-mer ID column
            merged_df['kmer_id'] = merged_df['kmer'].map(self.kmer_ids)
            
            # 重新排列列顺序 | Rearrange column order
            cols = ['kmer_id', 'kmer'] + [col for col in merged_df.columns if col not in ['kmer_id', 'kmer']]
            result_df = merged_df[cols]
        
        return result_df
