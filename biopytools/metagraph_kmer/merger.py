"""
🔗 K-mer结果合并模块 | K-mer Results Merger Module
"""

import pandas as pd
from pathlib import Path

class ResultsMerger:
    """结果合并器 | Results Merger"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def merge_results(self, coords_file: Path, counts_file: Path) -> Path:
        """
        合并坐标和计数信息 | Merge coordinates and count information
        
        Args:
            coords_file: 坐标文件路径 | Coordinates file path
            counts_file: 计数文件路径 | Counts file path
        
        Returns:
            最终输出文件路径 | Final output file path
        """
        self.logger.info("🔗 步骤4: 合并k-mer坐标和丰度信息 | Step 4: Merging k-mer coordinates and abundance")
        
        # 加载坐标信息 | Load coordinates
        self.logger.info("   加载坐标信息 | Loading coordinates...")
        coords_df = pd.read_csv(coords_file, sep='\t')
        self.logger.info(f"   ✅ 加载了 {len(coords_df):,} 条坐标记录 | Loaded {len(coords_df):,} coordinate records")
        
        # 加载计数信息 | Load counts
        self.logger.info("   加载计数信息 | Loading counts...")
        counts_df = pd.read_csv(
            counts_file,
            sep='\t',
            header=None,
            names=['kmer', 'count']
        )
        self.logger.info(f"   ✅ 加载了 {len(counts_df):,} 条计数记录 | Loaded {len(counts_df):,} count records")
        
        # 合并数据 | Merge data
        self.logger.info("   合并数据 | Merging data...")
        merged = coords_df.merge(
            counts_df,
            on='kmer',
            how='left'
        )
        
        # 填充缺失值 | Fill missing values
        merged['count'] = merged['count'].fillna(0).astype(int)
        
        # 添加存在性列 | Add presence column
        merged['presence'] = (merged['count'] > 0).astype(int)
        
        # 保存结果 | Save results
        output_file = self.config.output_path / "final_output.tsv"
        merged.to_csv(output_file, sep='\t', index=False)
        
        # 统计信息 | Statistics
        total_kmers = len(merged)
        present_kmers = merged['presence'].sum()
        absent_kmers = total_kmers - present_kmers
        total_abundance = merged['count'].sum()
        
        self.logger.info(f"\n📈 统计摘要 | Statistics Summary:")
        self.logger.info(f"   📊 库中k-mer总数 | Total k-mers in library: {total_kmers:,}")
        self.logger.info(f"   ✅ 查询中存在的k-mer | K-mers present in query: {present_kmers:,} ({present_kmers/total_kmers*100:.2f}%)")
        self.logger.info(f"   ❌ 查询中缺失的k-mer | K-mers absent in query: {absent_kmers:,} ({absent_kmers/total_kmers*100:.2f}%)")
        self.logger.info(f"   🔢 总丰度 | Total abundance: {total_abundance:,}")
        if present_kmers > 0:
            self.logger.info(f"   📊 平均丰度 | Average abundance: {total_abundance/present_kmers:.2f}")
        
        self.logger.info(f"\n✅ 最终结果已保存 | Final results saved: {output_file}")
        
        return output_file
