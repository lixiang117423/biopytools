"""
Minimap2数据处理模块 | Minimap2 Data Processing Module
"""

import pandas as pd
import os
from .utils import IntervalProcessor

class PAFProcessor:
    """PAF文件处理器 | PAF File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.interval_processor = IntervalProcessor(logger, config)
    
    def load_paf_file(self):
        """加载PAF文件 | Load PAF file"""
        paf_file = self.config.paf_file
        
        self.logger.info(f"加载PAF文件 | Loading PAF file: {paf_file}")
        
        try:
            # 先读取几行检查列数 | Read a few lines to check column count
            with open(paf_file, 'r') as f:
                first_line = f.readline().strip()
                columns_count = len(first_line.split('\t'))
            
            self.logger.info(f"检测到PAF文件有 {columns_count} 列 | Detected {columns_count} columns in PAF file")
            
            # 定义PAF文件的前12列（标准列） | Define first 12 columns (standard columns) of PAF file
            base_column_names = [
                "query_name", "query_len", "query_start", "query_end", "strand",
                "target_name", "target_len", "target_start", "target_end", "number_match",
                "all_len", "map_quality"
            ]
            
            # 为额外的列生成名称 | Generate names for additional columns
            column_names = base_column_names[:]
            for i in range(12, columns_count):
                column_names.append(f"tag_{i-11}")
            
            self.logger.info(f"使用列名 | Using column names: {column_names}")
            
            # 读取PAF文件 | Read PAF file
            df = pd.read_csv(paf_file, sep='\t', header=None, names=column_names)
            self.logger.info(f"成功加载PAF文件，共 {len(df)} 行 | Successfully loaded PAF file with {len(df)} rows")
            
            # 检查关键列是否存在 | Check if key columns exist
            required_columns = ["query_name", "query_len", "query_start", "query_end", "number_match"]
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                self.logger.error(f"缺少必要的列 | Missing required columns: {missing_columns}")
                return None
            
            # 显示数据框信息 | Show dataframe info
            self.logger.info(f"数据框列名 | DataFrame columns: {list(df.columns)}")
            self.logger.info(f"前几行数据 | First few rows:")
            for i, (_, row) in enumerate(df.head(3).iterrows()):
                self.logger.info(f"  行 {i+1} | Row {i+1}: query={row['query_name']}, len={row['query_len']}, match={row['number_match']}")
            
            return df
            
        except Exception as e:
            self.logger.error(f"加载PAF文件失败 | Failed to load PAF file: {e}")
            return None
    
    def filter_paf_data(self, df):
        """筛选PAF数据 | Filter PAF data"""
        self.logger.info("开始筛选PAF数据 | Starting PAF data filtering")
        
        original_count = len(df)
        self.logger.info(f"原始记录数 | Original records: {original_count}")
        
        # 找到包含tp信息的列 | Find column containing tp information
        tp_column = None
        for col in df.columns:
            if col.startswith('tag_'):
                # 检查这一列是否包含tp:A:S或tp:A:P | Check if this column contains tp:A:S or tp:A:P
                sample_values = df[col].dropna().head(10)
                if any('tp:A:' in str(val) for val in sample_values):
                    tp_column = col
                    self.logger.info(f"找到tp标签列 | Found tp tag column: {col}")
                    break
        
        if tp_column is None:
            self.logger.warning("未找到tp标签列，跳过tp过滤 | No tp tag column found, skipping tp filtering")
            df_filtered = df.copy()
        else:
            # 根据tp_type参数构建筛选条件 | Build filter condition based on tp_type parameter
            if self.config.tp_type.upper() == 'S':
                filter_pattern = 'S'
                filter_desc = "tp_tag包含'S'的记录 | Records with 'S' in tp_tag"
            elif self.config.tp_type.upper() == 'P':
                filter_pattern = 'P'
                filter_desc = "tp_tag包含'P'的记录 | Records with 'P' in tp_tag"
            elif self.config.tp_type.upper() == 'SP':
                filter_pattern = 'S|P'
                filter_desc = "tp_tag包含'S'或'P'的记录 | Records with 'S' or 'P' in tp_tag"
            else:
                self.logger.warning(f"无效的tp_type参数 | Invalid tp_type parameter: {self.config.tp_type}，使用默认值'P' | Using default 'P'")
                filter_pattern = 'P'
                filter_desc = "tp_tag包含'P'的记录 | Records with 'P' in tp_tag"
            
            # 步骤2: 筛选指定tp类型的行 | Filter rows with specified tp type
            df_filtered = df[df[tp_column].str.contains(filter_pattern, na=False, regex=True)]
            self.logger.info(f"筛选 {filter_desc}: {len(df_filtered)}")
            
            # 显示tp标签的分布 | Show distribution of tp tags
            if len(df_filtered) > 0:
                tp_counts = df_filtered[tp_column].value_counts()
                self.logger.info(f"tp标签分布 | tp tag distribution:")
                for tp_val, count in tp_counts.head(10).items():
                    self.logger.info(f"  {tp_val}: {count}")
            else:
                self.logger.warning(f"没有找到符合条件的tp类型记录 | No records found with tp type: {self.config.tp_type}")
        # 步骤3: 筛选number_match >= min_match_length的行 | Filter rows with number_match >= min_match_length
        df_filtered = df_filtered[df_filtered['number_match'] >= self.config.min_match_length]
        self.logger.info(f"筛选number_match >= {self.config.min_match_length}的记录 | Records with number_match >= {self.config.min_match_length}: {len(df_filtered)}")
        
        # 显示number_match的分布 | Show distribution of number_match
        if len(df_filtered) > 0:
            self.logger.info(f"number_match统计 | number_match statistics:")
            self.logger.info(f"  最大值 | Max: {df_filtered['number_match'].max():,}")
            self.logger.info(f"  最小值 | Min: {df_filtered['number_match'].min():,}")
            self.logger.info(f"  平均值 | Mean: {df_filtered['number_match'].mean():,.0f}")
        
        # 步骤4: 提取指定列 | Extract specified columns
        selected_columns = [
            'query_name', 'query_len', 'query_start', 'query_end', 'strand',
            'target_name', 'target_len', 'target_start', 'target_end', 'number_match'
        ]
        
        # 确保所有需要的列都存在 | Ensure all required columns exist
        available_columns = [col for col in selected_columns if col in df_filtered.columns]
        missing_columns = [col for col in selected_columns if col not in df_filtered.columns]
        
        if missing_columns:
            self.logger.warning(f"缺少部分列 | Missing some columns: {missing_columns}")
        
        df_result = df_filtered[available_columns].copy()
        self.logger.info(f"提取指定列，最终筛选结果 | Extracted specified columns, final filtered results: {len(df_result)}")
        
        # 显示筛选后的前几条记录 | Show first few filtered records
        if len(df_result) > 0:
            self.logger.info(f"筛选后的记录样本 | Sample of filtered records:")
            for i, (_, row) in enumerate(df_result.head(3).iterrows()):
                self.logger.info(f"  #{i+1}: {row['query_name']} {row['query_start']:,}-{row['query_end']:,} (match: {row['number_match']:,})")
        
        return df_result
    
    def find_unmapped_regions(self, df_filtered):
        """找到未比对的区间 | Find unmapped regions"""
        self.logger.info("开始查找未比对区间 | Starting to find unmapped regions")
        
        # 使用减法策略 | Use subtraction strategy
        unmapped_regions = self.interval_processor.find_unmapped_intervals_by_subtraction(df_filtered)
        
        self.logger.info(f"找到 {len(unmapped_regions)} 个未比对区间（长度 >= {self.config.min_unmapped_length}bp）")
        self.logger.info(f"Found {len(unmapped_regions)} unmapped regions (length >= {self.config.min_unmapped_length}bp)")
        
        return unmapped_regions
    
    def save_bed_file(self, unmapped_regions):
        """保存BED格式文件 | Save BED format file"""
        bed_file = self.config.bed_file
        
        self.logger.info(f"保存BED文件 | Saving BED file: {bed_file}")
        
        try:
            with open(bed_file, 'w') as f:
                for region in unmapped_regions:
                    # # BED格式：query_name\tstart\tend | BED format: query_name\tstart\tend
                    # f.write(f"{region['query_name']}\t{region['start']}\t{region['end']}\n")
                    # 将内部的1-based, closed [start, end] 转换为 BED 格式的 0-based, half-open [start-1, end)
                    # To convert our internal 1-based, closed [start, end] to BED format's 0-based, half-open [start-1, end)
                    bed_start = region['start'] - 1
                    bed_end = region['end']
                    f.write(f"{region['query_name']}\t{bed_start}\t{bed_end}\n")
            
            self.logger.info(f"BED文件保存成功 | BED file saved successfully: {bed_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"保存BED文件失败 | Failed to save BED file: {e}")
            return False
