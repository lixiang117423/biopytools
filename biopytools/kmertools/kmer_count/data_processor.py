"""
 数据处理模块|Data Processing Module
"""

import sqlite3
from typing import Dict, List, Optional
import os
import pandas as pd
from pathlib import Path

class DataProcessor:
    """ 数据处理器|Data Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.kmer_ids = {}
        self.bed_info = None

        # SQLite相关属性|SQLite related attributes
        self.kmer_db_conn: Optional[sqlite3.Connection] = None
        self.kmer_db_path: Optional[Path] = None
        self._use_sqlite: bool = False  # 是否使用SQLite|Whether to use SQLite
    
    def parse_kmer_library(self) -> None:
        """解析k-mer库（支持大规模数据）|Parse k-mer library (support large-scale data)

        根据数据规模自动选择存储策略:
        - 小规模(<1000万): 使用字典（更快）
        - 大规模(>=1000万): 使用SQLite（内存友好）

        Returns:
            None
        """
        self.logger.info("解析k-mer库|Parsing k-mer library")

        # 创建中间文件目录（在输出目录中）
        intermediate_dir = self.config.output_dir / "intermediate_files"
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        # 创建SQLite数据库
        db_path = intermediate_dir / "kmer_mapping.db"
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        self.logger.info(f"k-mer映射数据库|K-mer mapping database: {db_path}")

        # 性能优化配置（针对超算环境优化）
        cursor.execute('PRAGMA journal_mode=OFF')  # 完全禁用journal（最快，但有风险）
        cursor.execute('PRAGMA synchronous=OFF')  # 最大写入性能
        cursor.execute('PRAGMA cache_size=-64000')  # 64MB缓存
        cursor.execute('PRAGMA temp_store=MEMORY')  # 临时数据在内存中
        cursor.execute('PRAGMA mmap_size=30000000000')  # 30GB内存映射
        # PRIMARY KEY在kmer_id上，允许同一序列在多个位置出现
        # 注意：不在导入时创建索引，导入完成后批量创建（避免每次INSERT维护索引）
        cursor.execute('CREATE TABLE IF NOT EXISTS kmer_ids (kmer_id TEXT PRIMARY KEY, kmer_seq TEXT)')

        self.logger.info("性能优化配置|Performance optimizations: batch_size=1M, journal=OFF, synchronous=OFF, 无内存字典, 延迟创建索引")

        # 流式读取FASTA文件
        current_id = None
        count = 0
        batch_size = 1000000  # 100万批量写入（平衡内存和性能）
        batch_data = []

        self.logger.info(f"开始读取FASTA文件|Start reading FASTA file: {self.config.kmer_lib}")

        with open(self.config.kmer_lib, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_id = line[1:]
                elif line and current_id:
                    kmer_seq = line.upper()
                    # 修改：先插入kmer_id（PRIMARY KEY），再插入kmer_seq
                    batch_data.append((current_id, kmer_seq))

                    count += 1

                    # 批量写入数据库
                    if len(batch_data) >= batch_size:
                        # kmer_id是PRIMARY KEY，保证唯一性，直接用INSERT
                        cursor.executemany('INSERT INTO kmer_ids VALUES (?, ?)', batch_data)
                        conn.commit()
                        batch_data = []

                    # 每100万打印一次进度
                    if count % 1000000 == 0:
                        self.logger.info(f"已处理{count/1000000:.1f}M个k-mer|Processed {count/1000000:.1f}M k-mers")

        # 写入剩余数据
        if batch_data:
            cursor.executemany('INSERT INTO kmer_ids VALUES (?, ?)', batch_data)
            conn.commit()

        self.logger.info(f"数据导入完成，共{count}个k-mer|Data import completed, {count} k-mers total")

        # 数据导入完成后再创建索引（避免导入过程中维护索引开销）
        self.logger.info("创建索引以加速查询|Creating index for faster queries")
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_kmer_seq ON kmer_ids(kmer_seq)')
        conn.commit()

        self.logger.info(f"解析完成，共{count}个k-mer|Parsing completed, {count} k-mers total")
        self.logger.info(f"k-mer映射已保存到数据库|K-mer mapping saved to: {db_path}")

        # 保存数据库连接
        self.kmer_db_conn = conn
        self.kmer_db_path = db_path
        self._use_sqlite = count >= 10000000  # 大规模数据使用SQLite

        # 如果是小规模数据，关闭数据库节省内存
        if not self._use_sqlite:
            conn.close()
            self.kmer_db_conn = None
            self.logger.info("数据规模较小，使用内存字典|Small dataset, using in-memory dictionary")


    def parse_bed_file(self) -> pd.DataFrame:
        """ 解析BED文件并进行标准化处理|Parse BED file and standardize"""
        if not self.config.bed_file:
            return None
        
        self.logger.info(" 解析BED文件|Parsing BED file")
        
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
            self.logger.info(f" 成功读取BED文件前6列|Successfully read first 6 columns of BED file")
            
        except Exception as e:
            self.logger.warning(f" 读取前6列失败，尝试灵活读取|Failed to read 6 columns, trying flexible reading: {e}")
            
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
                
                self.logger.info(f" 灵活读取成功|Flexible reading successful")
                
            except Exception as e2:
                raise ValueError(f" BED文件格式错误，无法解析|BED file format error: {e2}")
        
        # 确保数据类型正确
        bed_df['start'] = bed_df['start'].astype(int)
        bed_df['end'] = bed_df['end'].astype(int)

        # 转换k-mer序列为大写
        original_kmer_count = len(bed_df)
        bed_df['kmer'] = bed_df['kmer'].str.upper()

        # 扩展k-mer库，添加反向互补序列（无条件执行，确保strand正确设置）
        bed_df = self._expand_bed_with_reverse_complement(bed_df)

        self.logger.info(f"已处理{len(bed_df)}个k-mer序列|Processed {len(bed_df)} k-mer sequences")

        # 保存处理后的BED文件到输出目录的intermediate_files子目录
        intermediate_dir = self.config.output_dir / "intermediate_files"
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        processed_bed_file = intermediate_dir / "processed_kmers.bed"
        bed_df.to_csv(processed_bed_file, sep='\t', header=False, index=False)
        self.logger.info(f"已保存处理后的BED文件|Saved processed BED file: {processed_bed_file}")

        # 更新配置以使用处理后的BED文件
        self.config.bed_file = processed_bed_file
        
        self.bed_info = bed_df
        self.logger.info(f" BED文件解析完成，共{len(bed_df)}条记录|BED file parsing completed, {len(bed_df)} records total")
        return bed_df

    def _expand_bed_with_reverse_complement(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """扩展BED文件，添加反向互补序列|Expand BED file with reverse complement sequences

        使用向量化操作替代逐行迭代，提升性能并降低内存占用
        Use vectorized operations instead of row iteration for better performance
        """
        self.logger.info("扩展BED文件，添加反向互补序列|Expanding BED file with reverse complement sequences")

        def generate_reverse_complement(seq):
            """生成反向互补序列|Generate reverse complement sequence"""
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(complement.get(base, base) for base in reversed(seq))

        def is_self_complement(seq):
            """检查序列是否与自身反向互补相同|Check if sequence is self-complementary"""
            return seq == generate_reverse_complement(seq)

        original_count = len(bed_df)

        # 向量化操作：筛选需要添加反向互补的序列
        # 使用pandas apply进行批量判断，避免逐行迭代
        mask = bed_df['kmer'].apply(
            lambda x: not is_self_complement(x) and generate_reverse_complement(x) != x
        )

        # 创建反向互补DataFrame（只对需要的序列）
        rc_df = bed_df[mask].copy()
        rc_df['kmer'] = rc_df['kmer'].apply(generate_reverse_complement)
        rc_df['strand'] = '-'

        # 确保原始序列标记为正向
        bed_df['strand'] = '+'

        # 使用pd.concat合并（比列表append更高效）
        expanded_df = pd.concat([bed_df, rc_df], ignore_index=True)
        added_count = len(expanded_df) - original_count

        self.logger.info(f"扩展完成：原始{original_count}条，新增{added_count}条反向互补序列|Expansion completed: {original_count} original, {added_count} reverse complement sequences added")

        return expanded_df
    # def parse_bed_file(self) -> pd.DataFrame:
    #     """ 解析BED文件|Parse BED file"""
    #     if not self.config.bed_file:
    #         return None
        
    #     self.logger.info(" 解析BED文件|Parsing BED file")
        
    #     # 读取BED文件|Read BED file
    #     bed_df = pd.read_csv(self.config.bed_file, sep='\t', header=None)
        
    #     # 设置列名|Set column names
    #     # col_names = ['chr', 'start', 'end', 'kmer']
    #     # if bed_df.shape[1] > 4:
    #     #     col_names.extend([f'tmp_{i}' for i in range(1, bed_df.shape[1] - 3)])
    #     # 设置列名|Set column names
    #     col_names = ['chr', 'start', 'end', 'kmer']
    #     if bed_df.shape[1] > 4:
    #         if bed_df.shape[1] == 6:
    #             col_names.extend(['score', 'strand'])
    #         else:
    #             col_names.extend([f'tmp_{i}' for i in range(1, bed_df.shape[1] - 3)])
        
    #     bed_df.columns = col_names[:bed_df.shape[1]]
        
    #     # 转换k-mer序列为大写|Convert k-mer sequences to uppercase
    #     bed_df['kmer'] = bed_df['kmer'].str.upper()
        
    #     self.bed_info = bed_df
    #     self.logger.info(f" BED文件解析完成，共{len(bed_df)}条记录|BED file parsing completed, {len(bed_df)} records total")
    #     return bed_df
    
    # def parse_jellyfish_output(self, count_file: Path, sample_name: str) -> pd.DataFrame:
    #     """ 解析Jellyfish输出|Parse Jellyfish output"""
    #     self.logger.info(f" 解析Jellyfish输出|Parsing Jellyfish output for {sample_name}")
        
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
    #     self.logger.info(f" 解析完成，共{len(df)}个k-mer结果|Parsing completed, {len(df)} k-mer results")
    #     return df

    def parse_jellyfish_output(self, count_file: Path, sample_name: str) -> pd.DataFrame:
        """ 解析Jellyfish输出|Parse Jellyfish output"""
        self.logger.info(f" 解析Jellyfish输出|Parsing Jellyfish output for {sample_name}")
        
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
        
        self.logger.info(f" 解析完成，去重后共{len(df)}个k-mer结果|Parsing completed, {len(df)} unique k-mer results after deduplication")
        return df
    
    # 单个样品合并bed文件
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """🔗 单个样品与BED文件合并|Merge single sample with BED file"""
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并|Merging sample {sample_name} with BED file")
        
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
            
    #         self.logger.info(f" 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f" 样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    # def merge_single_sample_with_bed(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
    #     """🔗 单个样品与BED文件合并|Merge single sample with BED file"""
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并|Merging sample {sample_name} with BED file")
        
    #     if self.bed_info is not None:
    #         #  不要在这里创建unique_key，直接用kmer合并
    #         result_df = pd.merge(self.bed_info, sample_df, on='kmer', how='left')
    #         # 填充缺失值为0
    #         result_df[sample_name] = result_df[sample_name].fillna(0).astype(int)
            
    #         # 删除丰度为0的行
    #         result_df = result_df[result_df[sample_name] > 0]
            
    #         #  在这里创建unique_key用于后续合并
    #         result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            
    #         result_df = result_df.sort_values(['chr', 'start'])
            
    #         self.logger.info(f" 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         #  添加：删除丰度为0的行|Add: Remove rows with abundance = 0
    #         sample_df = sample_df[sample_df[sample_name] > 0]

    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f" 样品{sample_name}合并完成，共{len(result_df)}行")
    #     return result_df
    
    def merge_sample_results(self, sample_results: List[pd.DataFrame]) -> pd.DataFrame:
        """🔗 合并所有样本结果|Merge all sample results"""
        self.logger.info("🔗 合并样本结果|Merging sample results")
        
        # 从第一个样本开始|Start with first sample
        merged_df = sample_results[0].copy()
        
        # 逐个合并其他样本|Merge other samples one by one
        for df in sample_results[1:]:
            merged_df = pd.merge(merged_df, df, on='kmer', how='outer')
        
        # 填充缺失值为0|Fill missing values with 0
        sample_columns = [col for col in merged_df.columns if col != 'kmer']
        merged_df[sample_columns] = merged_df[sample_columns].fillna(0).astype(int)
        
        self.logger.info(f" 合并完成，共{len(merged_df)}个k-mer，{len(sample_columns)}个样本")
        return merged_df
    
    def add_annotation_info(self, merged_df: pd.DataFrame) -> pd.DataFrame:
        """ 添加注释信息|Add annotation information"""
        if self.bed_info is not None:
            # 使用BED信息|Use BED information
            self.logger.info(" 添加BED注释信息|Adding BED annotation information")
            
            # 合并BED信息|Merge BED information
            result_df = pd.merge(self.bed_info, merged_df, on='kmer', how='left')
            
            # 填充缺失的丰度值为0|Fill missing abundance values with 0
            sample_columns = [col for col in merged_df.columns if col != 'kmer']
            result_df[sample_columns] = result_df[sample_columns].fillna(0).astype(int)
            
            # 按照染色体和位置排序|Sort by chromosome and position
            result_df = result_df.sort_values(['chr', 'start'])
            
        else:
            # 使用k-mer ID信息|Use k-mer ID information
            self.logger.info(" 添加k-mer ID信息|Adding k-mer ID information")

            # 添加k-mer ID列（使用批量查询优化）|Add k-mer ID column (using batch query)
            if self._use_sqlite:
                kmer_ids_map = self.get_kmer_ids_batch(merged_df['kmer'].tolist())
                merged_df['kmer_id'] = merged_df['kmer'].map(lambda x: kmer_ids_map.get(x.upper(), ''))
            else:
                merged_df['kmer_id'] = merged_df['kmer'].map(self.kmer_ids)

            # 重新排列列顺序|Rearrange column order
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
    #     self.logger.info(f"🔗 样品{sample_name}与BED文件合并|Merging sample {sample_name} with BED file")
        
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
            
    #         self.logger.info(f" 已删除丰度为0的行，剩余{len(result_df)}行")
    #     else:
    #         # 添加k-mer ID信息
    #         sample_df['kmer_id'] = sample_df['kmer'].map(self.kmer_ids)
    #         # 创建unique key：kmer_id+kmer
    #         sample_df['unique_key'] = sample_df['kmer_id'] + '_' + sample_df['kmer']
    #         # 删除丰度为0的行
    #         sample_df = sample_df[sample_df[sample_name] > 0]

    #         cols = ['unique_key', 'kmer_id', 'kmer', sample_name]
    #         result_df = sample_df[cols]
        
    #     self.logger.info(f" 样品{sample_name}合并完成，共{len(result_df)}行")
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
        """合并单个样品与BED文件（内存优化版本）

        优化策略：
        - 小规模(<100万行): 使用反向merge（快速）
        - 大规模(>=100万行): 使用字典映射（内存友好）
        """
        self.logger.info(f"合并样品{sample_name}与BED文件|Merging sample {sample_name} with BED file")

        if self.bed_info is not None:
            bed_size = len(self.bed_info)
            sample_size = len(sample_df)

            self.logger.info(f"BED记录数: {bed_size:,}, 样本k-mer数: {sample_size:,}|BED records: {bed_size:,}, sample k-mers: {sample_size:,}")

            # 根据数据规模选择策略
            if bed_size < 1000000:
                # 方案A：反向merge（小规模，快速）
                self.logger.info("使用反向merge策略|Using reverse merge strategy")
                result_df = self._merge_by_reverse_merge(sample_df, sample_name)
            else:
                # 方案B：字典映射（大规模，内存友好）
                self.logger.info("使用字典映射策略|Using dictionary mapping strategy")
                result_df = self._merge_by_dict_mapping(sample_df, sample_name)

        else:
            # 无BED文件的情况
            result_df = sample_df.copy()
            result_df['unique_key'] = result_df['kmer']
            result_df = result_df[result_df[sample_name] > 0]
            cols = ['unique_key', 'kmer', sample_name]
            result_df = result_df[cols]

        self.logger.info(f"样品{sample_name}合并完成，共{len(result_df)}行|Sample {sample_name} merge completed, {len(result_df)} rows total")
        return result_df

    def _merge_by_reverse_merge(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
        """方案A：反向merge（小规模数据）

        优化：将小表（样本）与大表（BED）merge，而不是相反
        原来的方式：pd.merge(bed_info, sample_df) → 结果=bed大小（巨大）
        优化后：pd.merge(sample_df, bed_info) → 结果=sample大小（小）
        """
        # 只merge检测到的k-mer（how='inner'）
        merged = pd.merge(sample_df, self.bed_info, on='kmer', how='inner')
        merged[sample_name] = merged[sample_name].fillna(0).astype(int)

        # 按基因组位置分组，选择最佳检测结果（避免重复计数）
        final_results = []

        for (chr_name, start, end), group in merged.groupby(['chr', 'start', 'end']):
            # 获取正向和反向序列的丰度
            forward_rows = group[group['strand'] == '+']
            reverse_rows = group[group['strand'] == '-']

            forward_abundance = forward_rows[sample_name].iloc[0] if len(forward_rows) > 0 else 0
            reverse_abundance = reverse_rows[sample_name].iloc[0] if len(reverse_rows) > 0 else 0

            # 选择检测到的结果（不求和）
            if forward_abundance > 0 and reverse_abundance > 0:
                # 两个方向都检测到，选择正向
                final_row = forward_rows.iloc[0].copy()
                final_row[sample_name] = forward_abundance
                final_row['strand'] = '+'
            elif forward_abundance > 0:
                final_row = forward_rows.iloc[0].copy()
                final_row[sample_name] = forward_abundance
                final_row['strand'] = '+'
            elif reverse_abundance > 0:
                final_row = forward_rows.iloc[0].copy()
                final_row[sample_name] = reverse_abundance
                final_row['strand'] = '-'
            else:
                continue

            final_results.append(final_row)

        if final_results:
            result_df = pd.DataFrame(final_results)
            result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            result_df = result_df.sort_values(['chr', 'start'])
        else:
            result_df = pd.DataFrame()

        return result_df

    def _merge_by_dict_mapping(self, sample_df: pd.DataFrame, sample_name: str) -> pd.DataFrame:
        """方案B：字典映射（大规模数据，内存优化）

        优点：
        - 内存占用极低（只保留样本的k-mer映射）
        - 处理速度快（O(n)复杂度）
        - 抗住上百G的kmer.fasta文件
        """
        # 构建kmer → 位置信息的映射字典
        # 只保留样本中实际检测到的k-mer的位置信息
        sample_kmers = set(sample_df['kmer'].tolist())

        self.logger.info(f"构建k-mer位置映射字典，样本中{len(sample_kmers):,}个k-mer|Building k-mer position mapping dictionary, {len(sample_kmers):,} k-mers in sample")

        # 从bed_info中提取样本k-mer的位置信息
        kmer_to_pos = {}
        for kmer in sample_kmers:
            # 使用布尔索引加速查询
            matches = self.bed_info[self.bed_info['kmer'] == kmer]
            if len(matches) > 0:
                # 如果有多个位置，取第一个
                kmer_to_pos[kmer] = matches.iloc[0].to_dict()

        self.logger.info(f"映射字典构建完成，共{len(kmer_to_pos):,}个k-mer有位置信息|Mapping dictionary built, {len(kmer_to_pos):,} k-mers have position info")

        # 为sample_df添加位置信息
        result_rows = []
        for _, row in sample_df.iterrows():
            kmer = row['kmer']
            abundance = row[sample_name]

            if kmer in kmer_to_pos:
                pos_info = kmer_to_pos[kmer]
                result_row = {
                    'chr': pos_info['chr'],
                    'start': pos_info['start'],
                    'end': pos_info['end'],
                    'kmer': kmer,
                    'score': pos_info['score'],
                    'strand': pos_info['strand'],
                    sample_name: abundance
                }
                result_rows.append(result_row)

        if result_rows:
            result_df = pd.DataFrame(result_rows)
            result_df['unique_key'] = result_df['chr'] + '_' + result_df['start'].astype(str) + '_' + result_df['end'].astype(str) + '_' + result_df['kmer']
            result_df = result_df.sort_values(['chr', 'start'])
        else:
            result_df = pd.DataFrame()

        return result_df

    def merge_reverse_complement_abundance(self, sample_df: pd.DataFrame, sample_name: str) -> dict:
        """ 合并反向互补序列的丰度到原始序列|Merge reverse complement abundance to original sequences"""
        self.logger.info(f" 合并反向互补序列丰度|Merging reverse complement abundance")
        
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
        
        self.logger.info(f" 完成反向互补序列合并，共{len(abundance_map)}个位置|Reverse complement merge completed, {len(abundance_map)} positions total")
        return abundance_map

    def get_kmer_id(self, kmer_seq: str) -> Optional[str]:
        """获取k-mer对应的ID（自动选择查询方式）|Get k-mer ID (auto-select query method)

        Args:
            kmer_seq: k-mer序列|k-mer sequence

        Returns:
            k-mer ID或None|k-mer ID or None
        """
        kmer_seq = kmer_seq.upper()

        # 优先使用内存字典
        if kmer_seq in self.kmer_ids:
            return self.kmer_ids[kmer_seq]

        # 如果有SQLite，查询数据库
        if self._use_sqlite and self.kmer_db_conn:
            cursor = self.kmer_db_conn.cursor()
            cursor.execute('SELECT kmer_id FROM kmer_ids WHERE kmer_seq = ?', (kmer_seq,))
            result = cursor.fetchone()
            if result:
                # 缓存到内存
                self.kmer_ids[kmer_seq] = result[0]
                return result[0]

        return None  # 未找到

    def get_kmer_ids_batch(self, kmer_seqs: list) -> dict:
        """批量获取k-mer ID|Batch get k-mer IDs

        Args:
            kmer_seqs: k-mer序列列表|k-mer sequence list

        Returns:
            序列到ID的映射字典|Sequence to ID mapping dictionary
        """
        result = {}
        missing_seqs = []

        # 先从内存字典查找
        for seq in kmer_seqs:
            seq_upper = seq.upper()
            if seq_upper in self.kmer_ids:
                result[seq_upper] = self.kmer_ids[seq_upper]
            else:
                missing_seqs.append(seq_upper)

        # 批量查询SQLite
        if missing_seqs and self._use_sqlite and self.kmer_db_conn:
            cursor = self.kmer_db_conn.cursor()

            # 使用IN查询批量获取
            placeholders = ','.join(['?' for _ in missing_seqs])
            query = f'SELECT kmer_seq, kmer_id FROM kmer_ids WHERE kmer_seq IN ({placeholders})'
            cursor.execute(query, missing_seqs)

            for row in cursor.fetchall():
                result[row[0]] = row[1]
                # 缓存到内存
                self.kmer_ids[row[0]] = row[1]

        return result



    def load_existing_database(self, db_path: Path):
        """加载已有的k-mer映射数据库|Load existing k-mer mapping database

        用于断点续传时恢复SQLite连接
        Used for checkpoint resumption to restore SQLite connection

        Args:
            db_path: 数据库文件路径|Database file path
        """
        if not db_path.exists():
            raise FileNotFoundError(f"数据库文件不存在|Database file not found: {db_path}")

        self.logger.info(f"加载已有数据库|Loading existing database: {db_path}")

        # 创建SQLite连接
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 检查表是否存在
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='kmer_ids'")
        if not cursor.fetchone():
            raise ValueError(f"数据库格式错误，缺少kmer_ids表|Invalid database format, missing kmer_ids table")

        # 获取k-mer数量
        cursor.execute('SELECT COUNT(*) FROM kmer_ids')
        count = cursor.fetchone()[0]
        self.logger.info(f"数据库包含{count}个k-mer|Database contains {count} k-mers")

        # 保存连接
        self.kmer_db_conn = conn
        self.kmer_db_path = db_path
        self._use_sqlite = count >= 10000000

        self.logger.info(f"数据库连接已恢复|Database connection restored")

    def load_existing_bed_file(self, bed_path: Path):
        """加载已有的处理后BED文件|Load existing processed BED file

        用于断点续传时恢复bed_info
        Used for checkpoint resumption to restore bed_info

        Args:
            bed_path: BED文件路径|BED file path
        """
        if not bed_path.exists():
            raise FileNotFoundError(f"BED文件不存在|BED file not found: {bed_path}")

        self.logger.info(f"加载已有BED文件|Loading existing BED file: {bed_path}")

        # 读取BED文件
        bed_df = pd.read_csv(
            bed_path,
            sep='	',
            header=None,  # 新版pandas要求用None代替False
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

        self.bed_info = bed_df
        self.logger.info(f"BED文件已加载，共{len(bed_df)}条记录|BED file loaded, {len(bed_df)} records")


    def cleanup(self):
        """清理资源|Cleanup resources

        关闭数据库连接，释放资源
        Close database connection and release resources
        """
        if self.kmer_db_conn:
            self.kmer_db_conn.close()
            self.kmer_db_conn = None
            self.logger.info("已关闭数据库连接|Database connection closed")
