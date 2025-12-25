"""
🪟 滑动窗口分析模块 | Sliding Window Analysis Module
"""

import pandas as pd

class WindowAnalyzer:
    """🪟 滑动窗口分析器 | Sliding Window Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    # def sliding_window_analysis(self, data_df: pd.DataFrame) -> pd.DataFrame:
    #     """🪟 滑动窗口分析 | Sliding window analysis"""
    #     if not self.config.bed_file:
    #         self.logger.warning("⚠️ 滑动窗口分析需要BED文件 | Sliding window analysis requires BED file")
    #         return None
        
    #     self.logger.info(f"🪟 开始滑动窗口分析 | Starting sliding window analysis")
    #     self.logger.info(f"📏 窗口大小: {self.config.window_size}bp, 步长: {self.config.step_size}bp")
        
    #     # 获取样本列 | Get sample columns
    #     sample_columns = [col for col in data_df.columns 
    #                      if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id'] 
    #                      and not col.startswith('tmp_')]
        
    #     # 按染色体分组处理 | Process by chromosome
    #     results = []
        
    #     for chrom in data_df['chr'].unique():
    #         chrom_data = data_df[data_df['chr'] == chrom].copy()
    #         chrom_data = chrom_data.sort_values('start')
            
    #         self.logger.info(f"🧬 处理染色体 | Processing chromosome: {chrom}")
            
    #         # 获取染色体范围 | Get chromosome range
    #         min_pos = chrom_data['start'].min()
    #         max_pos = chrom_data['end'].max()
            
    #         # 滑动窗口 | Sliding window
    #         window_start = min_pos
    #         while window_start <= max_pos:
    #             window_end = window_start + self.config.window_size
                
    #             # 获取窗口内的k-mer | Get k-mers in window
    #             window_kmers = chrom_data[
    #                 (chrom_data['start'] >= window_start) & 
    #                 (chrom_data['end'] <= window_end)
    #             ]
                
    #             if len(window_kmers) > 0:
    #                 # 计算每个样本的存在比例 | Calculate presence ratio for each sample
    #                 window_result = {
    #                     'chr': chrom,
    #                     'start': window_start,
    #                     'end': window_end,
    #                     'kmer_count': len(window_kmers)
    #                 }
                    
    #                 for sample in sample_columns:
    #                     present_kmers = (window_kmers[sample] > 0).sum()
    #                     total_kmers = len(window_kmers)
    #                     ratio = present_kmers / total_kmers if total_kmers > 0 else 0
    #                     window_result[sample] = round(ratio, 4)
                    
    #                 results.append(window_result)
                
    #             window_start += self.config.step_size
        
    #     if results:
    #         window_df = pd.DataFrame(results)
    #         window_df = window_df.sort_values(['chr', 'start'])
    #         self.logger.info(f"✅ 滑动窗口分析完成，共{len(window_df)}个窗口")
    #         return window_df
    #     else:
    #         self.logger.warning("⚠️ 未生成滑动窗口结果 | No sliding window results generated")
    #         return None

    def sliding_window_analysis(self, data_df: pd.DataFrame) -> pd.DataFrame:
        """🪟 滑动窗口分析 | Sliding window analysis"""
        if not self.config.bed_file:
            self.logger.warning("⚠️ 滑动窗口分析需要BED文件 | Sliding window analysis requires BED file")
            return None
        
        self.logger.info(f"🪟 开始滑动窗口分析 | Starting sliding window analysis")
        self.logger.info(f"📏 窗口大小: {self.config.window_size}bp, 步长: {self.config.step_size}bp")
        
        # 🔥 确保start和end列为整数类型 | Ensure start and end columns are integer type
        data_df = data_df.copy()
        data_df['start'] = data_df['start'].astype(int)
        data_df['end'] = data_df['end'].astype(int)
        
        # 获取样本列 | Get sample columns
        sample_columns = [col for col in data_df.columns 
                        if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id', 'score', 'strand'] 
                        and not col.startswith('tmp_')]
        
        # 按染色体分组处理 | Process by chromosome
        results = []
        
        for chrom in data_df['chr'].unique():
            chrom_data = data_df[data_df['chr'] == chrom].copy()
            chrom_data = chrom_data.sort_values('start')
            
            self.logger.info(f"🧬 处理染色体 | Processing chromosome: {chrom}")
            
            # 获取染色体范围 | Get chromosome range
            min_pos = chrom_data['start'].min()
            max_pos = chrom_data['end'].max()
            
            # 滑动窗口 | Sliding window
            window_start = min_pos
            while window_start <= max_pos:
                window_end = window_start + self.config.window_size
                
                # 获取窗口内的k-mer | Get k-mers in window
                window_kmers = chrom_data[
                    (chrom_data['start'] >= window_start) & 
                    (chrom_data['end'] <= window_end)
                ]
                
                if len(window_kmers) > 0:
                    # 计算每个样本的存在比例 | Calculate presence ratio for each sample
                    window_result = {
                        'chr': chrom,
                        'start': window_start,
                        'end': window_end,
                        'kmer_count': len(window_kmers)
                    }
                    
                    for sample in sample_columns:
                        present_kmers = (window_kmers[sample] > 0).sum()
                        total_kmers = len(window_kmers)
                        ratio = present_kmers / total_kmers if total_kmers > 0 else 0
                        window_result[sample] = round(ratio, 4)
                    
                    results.append(window_result)
                
                window_start += self.config.step_size
        
        if results:
            window_df = pd.DataFrame(results)
            window_df = window_df.sort_values(['chr', 'start'])
            self.logger.info(f"✅ 滑动窗口分析完成，共{len(window_df)}个窗口")
            return window_df
        else:
            self.logger.warning("⚠️ 未生成滑动窗口结果 | No sliding window results generated")
            return None
    
    def create_binary_matrix(self, data_df: pd.DataFrame) -> pd.DataFrame:
        """🔢 创建0/1存在缺失矩阵 | Create 0/1 presence/absence matrix"""
        self.logger.info("🔢 创建0/1存在缺失矩阵 | Creating 0/1 presence/absence matrix")
        
        # 获取样本列 | Get sample columns
        sample_columns = [col for col in data_df.columns 
                         if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id'] 
                         and not col.startswith('tmp_')]
        
        # 复制数据框 | Copy dataframe
        binary_df = data_df.copy()
        
        # 转换为0/1矩阵 | Convert to 0/1 matrix
        for col in sample_columns:
            binary_df[col] = (binary_df[col] > 0).astype(int)
        
        self.logger.info("✅ 0/1矩阵创建完成 | 0/1 matrix creation completed")
        return binary_df
