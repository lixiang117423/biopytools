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
        """📋 解析BED文件并进行标准化处理 | Parse BED file and standardize"""
        if not self.config.bed_file:
            return None
        
        self.logger.info("📋 解析BED文件 | Parsing BED file")
        
        try:
            # 读取BED文件前6列
            bed_df = pd.read_csv(
                self.config.bed_file, 
                sep='\t', 
                header=None,
                usecols=[0, 1, 2, 3, 4, 5],
                names=['chr', 'start', 'end', 'kmer', 'score', 'strand'],
                dtype={
                    'chr': str,
                    'start': int, 
                    'end': int,
                    'kmer': str,
                    'score': str,
                    'strand': str
                }
            )
            self.logger.info(f"✅ 成功读取BED文件前6列 | Successfully read first 6 columns of BED file")
            
        except Exception as e:
            self.logger.warning(f"⚠️ 读取前6列失败，尝试灵活读取 | Failed to read 6 columns, trying flexible reading: {e}")
            
            try:
                bed_df = pd.read_csv(
                    self.config.bed_file, 
                    sep='\t', 
                    header=None,
                    on_bad_lines='skip'
                )
                
                if bed_df.shape[1] >= 6:
                    bed_df = bed_df.iloc[:, :6]
                    bed_df.columns = ['chr', 'start', 'end', 'kmer', 'score', 'strand']
                elif bed_df.shape[1] >= 4:
                    bed_df = bed_df.iloc[:, :4]
                    bed_df.columns = ['chr', 'start', 'end', 'kmer']
                    bed_df['score'] = '0'
                    bed_df['strand'] = '.'
                else:
                    raise ValueError(f"BED文件列数不足，至少需要4列，实际有{bed_df.shape[1]}列")
                
                self.logger.info(f"✅ 灵活读取成功 | Flexible reading successful")
                
            except Exception as e2:
                raise ValueError(f"❌ BED文件格式错误，无法解析 | BED file format error: {e2}")
        
        # 确保数据类型正确
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)
        
        # 🔥 关键步骤：转换k-mer序列为大写
        original_kmer_count = len(bed_df)
        bed_df['kmer'] = bed_df['kmer'].str.upper()
        # 🔥 强制扩展，添加反向互补序列
        bed_df = self._expand_bed_with_reverse_complement(bed_df)

        self.logger.info(f"🔤 已将{original_kmer_count}个k-mer序列转换为大写 | Converted {original_kmer_count} k-mer sequences to uppercase")
        
        # 🔥 新增：如果启用扩展功能，添加反向互补序列
        if hasattr(self.config, 'expand_kmer_lib') and self.config.expand_kmer_lib:
            bed_df = self._expand_bed_with_reverse_complement(bed_df)
        
        # 🔥 保存处理后的BED文件到临时目录
        processed_bed_file = self.config.temp_dir / "processed_kmers.bed"
        bed_df.to_csv(processed_bed_file, sep='\t', header=False, index=False)
        self.logger.info(f"💾 已保存处理后的BED文件 | Saved processed BED file: {processed_bed_file}")
        
        # 更新配置以使用处理后的BED文件
        self.config.bed_file = processed_bed_file
        
        self.bed_info = bed_df
        self.logger.info(f"✅ BED文件解析完成，共{len(bed_df)}条记录 | BED file parsing completed, {len(bed_df)} records total")
        return bed_df

    def _expand_bed_with_reverse_complement(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """🔄 扩展BED文件，添加反向互补序列 | Expand BED file with reverse complement sequences"""
        self.logger.info("🔄 扩展BED文件，添加反向互补序列 | Expanding BED file with reverse complement sequences")
        
        def generate_reverse_complement(seq):
            """生成反向互补序列"""
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(complement.get(base, base) for base in reversed(seq))
        
        def is_self_complement(seq):
            """检查序列是否与自身反向互补相同"""
            return seq == generate_reverse_complement(seq)
        
        expanded_rows = []
        original_count = len(bed_df)
        
        for _, row in bed_df.iterrows():
            # 添加原始序列（标记为正向）
            original_row = row.copy()
            original_row['strand'] = '+'
            expanded_rows.append(original_row)
            
            # 生成反向互补序列
            rc_seq = generate_reverse_complement(row['kmer'])
            
            # 检查是否为自互补序列或与原序列相同
            if not is_self_complement(row['kmer']) and rc_seq != row['kmer']:
                rc_row = row.copy()
                rc_row['kmer'] = rc_seq
                rc_row['strand'] = '-'
                expanded_rows.append(rc_row)
        
        expanded_df = pd.DataFrame(expanded_rows)
        added_count = len(expanded_df) - original_count
        
        self.logger.info(f"🔄 扩展完成：原始{original_count}条，新增{added_count}条反向互补序列 | Expansion completed: {original_count} original, {added_count} reverse complement sequences added")
        
        return expanded_df
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
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """🔗 单个样品与BED文件合并 | Merge single sample with BED file"""
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并 | Merging sample {sample_name} with BED file")
        
    #     if self.bed_info is not None:
    #         # 🔥 不要在这里创建unique_key，直接用kmer合并
    #         result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
    #         # 填充缺失值为0
    #         result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
    #         # 删除丰度为0的行
    #         result_df = result_df[result_df[sample_name] > 0]
            
    #         # 🔥 在这里创建unique_key用于后续合并
    #         result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            
    #         result_df = result_df.sort_values(['chr', 'start'])
            
    #         self.logger.info(f"🗑️ 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         # 🔥 添加：删除丰度为0的行 | Add: Remove rows with abundance = 0
    #         sample_df = sample_df[sample_df[sample_name] > 0]

    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f"✅ 样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    
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
    
    # def merge_reverse_complement_abundance(self, sample_df, sample_name):
    #     """合并反向互补序列的丰度到原始序列"""
    #     abundance_map = {}
        
    #     for _, row in sample_df.iterrows():
    #         kmer_seq = row['kmer']
    #         abundance = row[sample_name]
            
    #         # 通过BED信息查找对应关系
    #         bed_matches = self.bed_info[self.bed_info['kmer'] == kmer_seq]
            
    #         for _, bed_row in bed_matches.iterrows():
    #             if bed_row['strand'] == '-':
    #                 # 反向互补序列，找到对应的正向序列
    #                 original_key = f"{bed_row['chr']}_{bed_row['start']}_{bed_row['end']}"
    #             else:
    #                 # 正向序列
    #                 original_key = f"{bed_row['chr']}_{bed_row['start']}_{bed_row['end']}"
                
    #             if original_key in abundance_map:
    #                 abundance_map[original_key] += abundance
    #             else:
    #                 abundance_map[original_key] = abundance
        
    #     return abundance_map
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """🔗 单个样品与BED文件合并（处理反向互补序列）| Merge single sample with BED file (handle reverse complement)"""
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并 | Merging sample {sample_name} with BED file")
        
    #     if self.bed_info is not None:
    #         # 首先进行反向互补序列的丰度合并
    #         merged_abundance = self.merge_reverse_complement_abundance(sample_df, sample_name)
            
    #         # 与BED信息合并，只保留正向序列（strand为'+'或'.'）
    #         forward_bed = self.bed_info[self.bed_info['strand'].isin(['+', '.'])].copy()
    #         result_df = pd.merge(forward_bed, sample_df, on='kmer', how='left')
            
    #         # 应用合并后的丰度值
    #         for index, row in result_df.iterrows():
    #             position_key = f"{row['chr']}_{row['start']}_{row['end']}"
    #             if position_key in merged_abundance:
    #                 result_df.at[index, sample_name] = merged_abundance[position_key]
    #             else:
    #                 result_df.at[index, sample_name] = 0
            
    #         # 确保丰度列为整数类型
    #         result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
    #         # 删除丰度为0的行
    #         result_df = result_df[result_df[sample_name] > 0]
            
    #         # 创建unique_key用于后续合并
    #         result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            
    #         result_df = result_df.sort_values(['chr', 'start'])
            
    #         self.logger.info(f"🗑️ 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         # 删除丰度为0的行
    #         sample_df = sample_df[sample_df[sample_name] > 0]

    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f"✅ 样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """合并单个样品与BED文件（新逻辑：保留原始序列，用strand标识方向）"""
    #     self.logger.info(f"合并样品{sample_name}与BED文件")
        
    #     if self.bed_info is not None:
    #         # 与BED信息合并
    #         result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
    #         result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
    #         # 按基因组位置分组，合并正向和反向序列的结果
    #         final_results = []
            
    #         for (chr_name, start, end), group in result_df.groupby(['chr', 'start', 'end']):
    #             # 获取正向和反向序列的丰度
    #             forward_rows = group[group['strand'] == '+']
    #             reverse_rows = group[group['strand'] == '-']
                
    #             forward_abundance = forward_rows[sample_name].iloc[0] if len(forward_rows) > 0 else 0
    #             reverse_abundance = reverse_rows[sample_name].iloc[0] if len(reverse_rows) > 0 else 0
                
    #             # 确定最终结果
    #             if forward_abundance > 0 and reverse_abundance > 0:
    #                 # 两个方向都检测到，选择丰度高的
    #                 if forward_abundance >= reverse_abundance:
    #                     final_row = forward_rows.iloc[0].copy()
    #                     final_row[sample_name] = forward_abundance
    #                     final_row['strand'] = '+'
    #                 else:
    #                     # 虽然检测到反向，但显示原始序列
    #                     final_row = forward_rows.iloc[0].copy()
    #                     final_row[sample_name] = reverse_abundance
    #                     final_row['strand'] = '-'
    #             elif forward_abundance > 0:
    #                 # 只检测到正向
    #                 final_row = forward_rows.iloc[0].copy()
    #                 final_row[sample_name] = forward_abundance
    #                 final_row['strand'] = '+'
    #             elif reverse_abundance > 0:
    #                 # 只检测到反向，但显示原始序列
    #                 final_row = forward_rows.iloc[0].copy()
    #                 final_row[sample_name] = reverse_abundance
    #                 final_row['strand'] = '-'
    #             else:
    #                 # 两个方向都没检测到，跳过
    #                 continue
                
    #             final_results.append(final_row)
            
    #         if final_results:
    #             result_df = pd.DataFrame(final_results)
    #             result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
    #             result_df = result_df.sort_values(['chr', 'start'])
    #         else:
    #             result_df = pd.DataFrame()
            
    #         self.logger.info(f"处理完成，剩余{len(result_df)}行")
    #     else:
    #         # 没有BED文件的情况保持不变
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         sample_df = sample_df[sample_df[sample_name] > 0]
    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f"样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
        """合并单个样品与BED文件（避免重复计数）"""
        self.logger.info(f"合并样品{sample_name}与BED文件")
        
        if self.bed_info is not None:
            # 与BED信息合并
            result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
            result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
            # 按基因组位置分组，选择最佳检测结果（避免重复计数）
            final_results = []
            
            for (chr_name, start, end), group in result_df.groupby(['chr', 'start', 'end']):
                # 获取正向和反向序列的丰度
                forward_rows = group[group['strand'] == '+']
                reverse_rows = group[group['strand'] == '-']
                
                forward_abundance = forward_rows[sample_name].iloc[0] if len(forward_rows) > 0 else 0
                reverse_abundance = reverse_rows[sample_name].iloc[0] if len(reverse_rows) > 0 else 0
                
                # 选择检测到的结果（不求和）
                if forward_abundance > 0 and reverse_abundance > 0:
                    # 两个方向都检测到，选择丰度高的，但避免重复计数
                    # 实际上这种情况下，丰度应该是相同的，选择正向作为代表
                    final_row = forward_rows.iloc[0].copy()
                    final_row[sample_name] = forward_abundance  # 不要 max()，直接用正向的
                    final_row['strand'] = '+'  # 单个文件中保留检测到的方向信息
                elif forward_abundance > 0:
                    final_row = forward_rows.iloc[0].copy()
                    final_row[sample_name] = forward_abundance
                    final_row['strand'] = '+'
                elif reverse_abundance > 0:
                    # 检测到反向，保留原始序列但标记为反向
                    final_row = forward_rows.iloc[0].copy()
                    final_row[sample_name] = reverse_abundance
                    final_row['strand'] = '-'
                else:
                    # 都没检测到，跳过
                    continue
                
                final_results.append(final_row)
            
            if final_results:
                result_df = pd.DataFrame(final_results)
                result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
                result_df = result_df.sort_values(['chr', 'start'])
            else:
                result_df = pd.DataFrame()
                
        else:
            # 无BED文件的情况
            sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
            sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
            sample_df = sample_df[sample_df[sample_name] > 0]
            cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
            result_df = sample_df[cols]
        
        self.logger.info(f"样品{sample_name}合并完成，共{len(result_df)}行")
        return result_df

    def merge_reverse_complement_abundance(self, sample_df: pd.DataFrame, sample_name: str) -> dict:
        """🔄 合并反向互补序列的丰度到原始序列 | Merge reverse complement abundance to original sequences"""
        self.logger.info(f"🔄 合并反向互补序列丰度 | Merging reverse complement abundance")
        
        abundance_map = {}
        
        for _, row in sample_df.iterrows():
            kmer_seq = row['kmer']
            abundance = row[sample_name]
            
            # 在BED信息中查找匹配的k-mer
            bed_matches = self.bed_info[self.bed_info['kmer'] == kmer_seq]
            
            for _, bed_row in bed_matches.iterrows():
                # 创建位置键，用于标识同一个基因组位置
                position_key = f"{bed_row['chr']}_{bed_row['start']}_{bed_row['end']}"
                
                # 累加丰度值（无论是正向还是反向序列）
                if position_key in abundance_map:
                    abundance_map[position_key] += abundance
                else:
                    abundance_map[position_key] = abundance
        
        self.logger.info(f"🔄 完成反向互补序列合并，共{len(abundance_map)}个位置")
        return abundance_map
