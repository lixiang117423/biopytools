"""
📊 存在性矩阵生成模块 | Presence Matrix Generation Module
"""

import pandas as pd
import numpy as np
import subprocess
import os
from pathlib import Path

class WindowAnalyzer:
    """🪟 窗口分析器 | Window Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def analyze_windows(self):
        """🔍 进行窗口分析 | Perform window analysis"""
        self.logger.info("🚀 开始窗口分析 | Starting window analysis")
        self.logger.info(f"📏 窗口大小: {self.config.window_size:,} bp | Window size: {self.config.window_size:,} bp")
        
        if self.config.overlapping:
            self.logger.info(f"👣 步长: {self.config.step_size:,} bp (重叠窗口) | Step size: {self.config.step_size:,} bp (overlapping windows)")
        else:
            self.logger.info("📦 使用非重叠窗口 | Using non-overlapping windows")
        
        try:
            # 读取样本列表 | Read sample list
            samples = self._read_sample_list()
            
            # 读取位置信息 | Read position information
            positions_df = self._read_position_info()
            if positions_df is None or len(positions_df) == 0:
                self.logger.error("❌ 无法读取k-mer位置信息，跳过窗口分析 | Cannot read k-mer position information, skipping window analysis")
                return False
            
            # 创建窗口 | Create windows
            windows = self._create_windows(positions_df)
            
            # 分析每个窗口 | Analyze each window
            window_results = self._analyze_each_window(windows, positions_df, samples)
            
            # 保存结果 | Save results
            self._save_window_results(window_results, samples)
            
            self.logger.info("✅ 窗口分析完成 | Window analysis completed")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 窗口分析失败: {e} | Window analysis failed: {e}")
            return False
    
    def _read_sample_list(self):
        """📋 读取样本列表 | Read sample list"""
        with open(self.config.sample_list_file, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]
        return samples
    
    def _read_position_info(self):
        """📍 读取位置信息 | Read position information"""
        self.logger.info("📖 读取k-mer位置信息 | Reading k-mer position information")
        
        if not self.config.genome_positions_file.exists():
            self.logger.error(f"❌ 位置文件不存在: {self.config.genome_positions_file} | Position file does not exist")
            return None
        
        try:
            positions_df = pd.read_csv(self.config.genome_positions_file, 
                                      sep='\t', header=None,
                                      names=['chr', 'start', 'end', 'kmer', 'score', 'strand'])
            
            self.logger.info(f"📊 读取到 {len(positions_df)} 个k-mer位置记录 | Read {len(positions_df)} k-mer position records")
            return positions_df
            
        except Exception as e:
            self.logger.error(f"❌ 读取位置文件失败: {e} | Failed to read position file: {e}")
            return None
    
    def _create_windows(self, positions_df):
        """🏗️ 创建窗口 | Create windows"""
        self.logger.info("🏗️ 创建基因组窗口 | Creating genome windows")
        
        windows = []
        
        # 按染色体分组 | Group by chromosome
        for chr_name in positions_df['chr'].unique():
            chr_positions = positions_df[positions_df['chr'] == chr_name]
            chr_start = chr_positions['start'].min()
            chr_end = chr_positions['end'].max()
            
            self.logger.info(f"🧬 处理染色体 {chr_name}: {chr_start:,} - {chr_end:,} bp | Processing chromosome {chr_name}: {chr_start:,} - {chr_end:,} bp")
            
            # 创建窗口 | Create windows
            current_start = chr_start
            window_id = 1
            
            while current_start < chr_end:
                window_end = min(current_start + self.config.window_size, chr_end)
                
                windows.append({
                    'chr': chr_name,
                    'start': current_start,
                    'end': window_end,
                    'window_id': f"{chr_name}_W{window_id:04d}"
                })
                
                current_start += self.config.step_size
                window_id += 1
                
                # 如果是非重叠窗口且当前窗口已经到达末尾，退出 | If non-overlapping and reached end, break
                if not self.config.overlapping and window_end >= chr_end:
                    break
        
        self.logger.info(f"📦 创建了 {len(windows)} 个窗口 | Created {len(windows)} windows")
        return windows
    
    def _analyze_each_window(self, windows, positions_df, samples):
        """🔍 分析每个窗口 | Analyze each window"""
        self.logger.info("📈 分析每个窗口中的k-mer存在比例 | Analyzing k-mer presence ratio in each window")
        
        window_results = []
        
        for i, window in enumerate(windows, 1):
            if i % 100 == 0 or i == len(windows):
                self.logger.info(f"🏃 处理窗口 {i}/{len(windows)} | Processing window {i}/{len(windows)}")
            
            # 找到窗口内的k-mer | Find k-mers within window
            window_kmers = positions_df[
                (positions_df['chr'] == window['chr']) &
                (positions_df['start'] >= window['start']) &
                (positions_df['end'] <= window['end'])
            ]['kmer'].tolist()
            
            if len(window_kmers) == 0:
                # 窗口内没有k-mer | No k-mers in window
                result = {
                    'chr': window['chr'],
                    'start': window['start'],
                    'end': window['end'],
                    'window_id': window['window_id'],
                    'total_kmers': 0
                }
                # 为每个样本添加0比例 | Add 0 ratio for each sample
                for sample in samples:
                    result[f"{sample}_ratio"] = 0.0
                    result[f"{sample}_present"] = 0
                    result[f"{sample}_total"] = 0
            else:
                # 计算每个样本在此窗口中的k-mer存在比例 | Calculate k-mer presence ratio for each sample in this window
                result = {
                    'chr': window['chr'],
                    'start': window['start'],
                    'end': window['end'],
                    'window_id': window['window_id'],
                    'total_kmers': len(window_kmers)
                }
                
                for sample in samples:
                    present_kmers = self._get_sample_present_kmers(sample, window_kmers)
                    ratio = len(present_kmers) / len(window_kmers) if len(window_kmers) > 0 else 0.0
                    
                    result[f"{sample}_ratio"] = ratio
                    result[f"{sample}_present"] = len(present_kmers)
                    result[f"{sample}_total"] = len(window_kmers)
            
            window_results.append(result)
        
        return window_results
    
    def _get_sample_present_kmers(self, sample_name, window_kmers):
        """🔍 获取样本中存在的窗口k-mer | Get present k-mers in sample for the window"""
        present_file = self.config.presence_results_dir / f"{sample_name}_present.txt"
        
        if not present_file.exists():
            return set()
        
        # 读取样本中存在的k-mer | Read present k-mers in sample
        present_kmers = set()
        try:
            with open(present_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if parts:
                            kmer = parts[0]
                            present_kmers.add(kmer)
        except Exception:
            pass
        
        # 返回窗口内存在的k-mer | Return k-mers present in the window
        return present_kmers.intersection(set(window_kmers))
    
    def _save_window_results(self, window_results, samples):
        """💾 保存窗口分析结果 | Save window analysis results"""
        self.logger.info("💾 保存窗口分析结果 | Saving window analysis results")
        
        # 创建DataFrame | Create DataFrame
        df = pd.DataFrame(window_results)
        
        # 保存完整结果 | Save complete results
        complete_file = self.config.window_results_dir / f"kmer{self.config.kmer_size}_window_analysis_complete.csv"
        df.to_csv(complete_file, index=False)
        
        # 创建类似BED格式的简化结果（只包含比例信息）| Create BED-like simplified results (ratios only)
        bed_columns = ['chr', 'start', 'end'] + [f"{sample}_ratio" for sample in samples]
        bed_df = df[bed_columns].copy()
        
        # 重命名列以符合BED格式 | Rename columns to match BED format
        ratio_columns = {f"{sample}_ratio": sample for sample in samples}
        bed_df = bed_df.rename(columns=ratio_columns)
        
        bed_file = self.config.window_results_dir / f"kmer{self.config.kmer_size}_window_presence_ratios.bed"
        bed_df.to_csv(bed_file, sep='\t', index=False, header=True)
        
        # 生成统计摘要 | Generate statistics summary
        self._generate_window_statistics(df, samples)
        
        self.logger.info(f"📊 窗口分析结果已保存 | Window analysis results saved:")
        self.logger.info(f"  📄 完整结果: {complete_file.name} | Complete results: {complete_file.name}")
        self.logger.info(f"  🛏️ BED格式: {bed_file.name} | BED format: {bed_file.name}")
    
    def _generate_window_statistics(self, df, samples):
        """📈 生成窗口统计信息 | Generate window statistics"""
        self.logger.info("📊 生成窗口分析统计 | Generating window analysis statistics")
        
        # 总体统计 | Overall statistics
        total_windows = len(df)
        windows_with_kmers = len(df[df['total_kmers'] > 0])
        
        self.logger.info(f"\n=== 📊 窗口分析统计 | Window Analysis Statistics ===")
        self.logger.info(f"📦 总窗口数: {total_windows:,} | Total windows: {total_windows:,}")
        self.logger.info(f"🎯 包含k-mer的窗口数: {windows_with_kmers:,} | Windows with k-mers: {windows_with_kmers:,}")
        self.logger.info(f"📏 窗口大小: {self.config.window_size:,} bp | Window size: {self.config.window_size:,} bp")
        
        if self.config.overlapping:
            self.logger.info(f"👣 步长: {self.config.step_size:,} bp | Step size: {self.config.step_size:,} bp")
        
        # 每个样本的统计 | Statistics for each sample
        stats_data = []
        for sample in samples:
            ratio_col = f"{sample}_ratio"
            if ratio_col in df.columns:
                ratios = df[ratio_col]
                stats_data.append({
                    'sample': sample,
                    'mean_ratio': ratios.mean(),
                    'median_ratio': ratios.median(),
                    'min_ratio': ratios.min(),
                    'max_ratio': ratios.max(),
                    'windows_with_100pct': len(ratios[ratios == 1.0]),
                    'windows_with_0pct': len(ratios[ratios == 0.0])
                })
        
        stats_df = pd.DataFrame(stats_data)
        stats_file = self.config.window_results_dir / f"kmer{self.config.kmer_size}_window_sample_stats.csv"
        stats_df.to_csv(stats_file, index=False)
        
        self.logger.info(f"\n📊 样本窗口统计摘要 | Sample window statistics summary:")
        for _, row in stats_df.iterrows():
            self.logger.info(f"  🧪 {row['sample']}: 平均比例={row['mean_ratio']:.3f}, 100%窗口={row['windows_with_100pct']}, 0%窗口={row['windows_with_0pct']}")
        
        self.logger.info(f"💾 详细统计已保存: {stats_file.name} | Detailed statistics saved: {stats_file.name}")


class PresenceMatrixGenerator:
    """📊 存在性矩阵生成器 | Presence Matrix Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_matrix(self):
        """🏗️ 生成0/1存在性矩阵 | Generate 0/1 presence matrix"""
        self.logger.info("🚀 开始生成存在性矩阵 | Starting presence matrix generation")
        
        try:
            # 读取样本列表 | Read sample list
            samples = self._read_sample_list()
            
            # 读取基因组k-mer | Read genome k-mers
            genome_kmers = self._read_genome_kmers()
            
            # 读取位置信息 | Read position information
            kmer_to_position = self._read_position_info()
            
            # 生成存在性矩阵 | Generate presence matrix
            presence_matrix = self._create_presence_matrix(samples, genome_kmers)
            
            # 保存结果 | Save results
            self._save_results(samples, genome_kmers, presence_matrix, kmer_to_position)
            
            # 生成统计摘要 | Generate statistics summary
            self._generate_statistics(samples, genome_kmers, presence_matrix)
            
            self.logger.info("✅ 存在性矩阵生成完成 | Presence matrix generation completed")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 存在性矩阵生成失败: {e} | Presence matrix generation failed: {e}")
            return False
    
    def _read_sample_list(self):
        """📋 读取样本列表 | Read sample list"""
        self.logger.info("📖 读取样本列表 | Reading sample list")
        
        with open(self.config.sample_list_file, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]
        
        self.logger.info(f"🎯 发现 {len(samples)} 个样本 | Found {len(samples)} samples")
        return samples
    
    def _read_genome_kmers(self):
        """🧬 读取基因组k-mer列表 | Read genome k-mer list"""
        self.logger.info("📖 读取基因组k-mer列表 | Reading genome k-mer list")
        
        # 转换为文本格式 | Convert to text format
        genome_kmers_txt = self.config.output_path / "genome_kmers_list.txt"
        
        try:
            with open(genome_kmers_txt, "w") as f:
                subprocess.run([self.config.unikmer_path, "view", str(self.config.genome_kmers_file)], 
                              stdout=f, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"❌ 转换基因组k-mer为文本格式失败: {e}")
        
        # 读取k-mer | Read k-mers
        genome_kmers = set()
        with open(genome_kmers_txt, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if parts:
                        kmer = parts[0]
                        genome_kmers.add(kmer)
        
        genome_kmers = list(genome_kmers)
        self.logger.info(f"🧮 总共有 {len(genome_kmers)} 个基因组k-mer | Total {len(genome_kmers)} genome k-mers")
        
        if len(genome_kmers) == 0:
            raise ValueError("❌ 没有找到基因组k-mer | No genome k-mers found")
        
        return genome_kmers
    
    def _read_position_info(self):
        """📍 读取位置信息 | Read position information"""
        self.logger.info("📖 读取位置信息 | Reading position information")
        
        kmer_to_position = {}
        
        if self.config.genome_positions_file.exists():
            try:
                positions_df = pd.read_csv(self.config.genome_positions_file, 
                                          sep='\t', header=None,
                                          names=['chr', 'start', 'end', 'kmer', 'score', 'strand'])
                
                self.logger.info(f"📊 读取到 {len(positions_df)} 个位置记录 | Read {len(positions_df)} position records")
                
                # 创建k-mer到位置的映射 | Create k-mer to position mapping
                for _, row in positions_df.iterrows():
                    kmer_to_position[row['kmer']] = f"{row['chr']}:{row['start']}-{row['end']}:{row['strand']}"
                    
            except Exception as e:
                self.logger.warning(f"⚠️ 读取位置文件失败: {e} | Failed to read position file: {e}")
        else:
            self.logger.warning(f"⚠️ 位置文件不存在: {self.config.genome_positions_file} | Position file does not exist")
        
        return kmer_to_position
    
    def _read_kmers_from_file(self, filepath):
        """📖 读取k-mer文件，返回k-mer集合 | Read k-mer file and return k-mer set"""
        kmers = set()
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if parts:
                            kmer = parts[0]
                            kmers.add(kmer)
        except FileNotFoundError:
            self.logger.warning(f"⚠️ 文件不存在: {filepath} | File not found: {filepath}")
        except Exception as e:
            self.logger.error(f"❌ 读取文件出错 {filepath}: {e} | Error reading file {filepath}: {e}")
        
        return kmers
    
    def _create_presence_matrix(self, samples, genome_kmers):
        """🏗️ 创建存在性矩阵 | Create presence matrix"""
        self.logger.info("🧮 生成存在性矩阵 | Generating presence matrix")
        
        n_samples = len(samples)
        n_kmers = len(genome_kmers)
        presence_matrix = np.zeros((n_kmers, n_samples), dtype=int)
        
        self.logger.info(f"📊 矩阵大小: {n_kmers} x {n_samples} | Matrix size: {n_kmers} x {n_samples}")
        
        # 处理每个样本 | Process each sample
        for sample_idx, sample_name in enumerate(samples):
            self.logger.info(f"🧪 处理样本 {sample_idx+1}/{n_samples}: {sample_name} | Processing sample {sample_idx+1}/{n_samples}: {sample_name}")
            
            present_file = self.config.presence_results_dir / f"{sample_name}_present.txt"
            present_kmers = self._read_kmers_from_file(present_file)
            
            self.logger.info(f"  🎯 样本 {sample_name} 中有 {len(present_kmers)} 个基因组k-mer存在 | Sample {sample_name} has {len(present_kmers)} genome k-mers present")
            
            # 更新存在性矩阵 | Update presence matrix
            for kmer_idx, kmer in enumerate(genome_kmers):
                if kmer in present_kmers:
                    presence_matrix[kmer_idx, sample_idx] = 1
        
        return presence_matrix
    
    def _save_results(self, samples, genome_kmers, presence_matrix, kmer_to_position):
        """💾 保存结果 | Save results"""
        self.logger.info("💾 保存结果矩阵 | Saving result matrix")
        
        # 添加位置信息到k-mer标识符 | Add position information to k-mer identifiers
        kmer_info = []
        for kmer in genome_kmers:
            position = kmer_to_position.get(kmer, "Unknown")
            kmer_info.append(f"{kmer}|{position}")
        
        # 创建结果DataFrame | Create result DataFrame
        result_df = pd.DataFrame(presence_matrix, 
                               index=kmer_info, 
                               columns=samples)
        
        # 保存结果 | Save results
        output_file = self.config.output_path / f"kmer{self.config.kmer_size}_presence_matrix.csv"
        result_df.to_csv(output_file)
        
        self.logger.info(f"✅ 存在性矩阵已保存: {output_file} | Presence matrix saved: {output_file}")
    
    def _generate_statistics(self, samples, genome_kmers, presence_matrix):
        """📈 生成统计摘要 | Generate statistics summary"""
        self.logger.info("📊 生成分析摘要 | Generating analysis summary")
        
        n_samples = len(samples)
        n_kmers = len(genome_kmers)
        
        self.logger.info(f"\n=== 📊 分析摘要 | Analysis Summary ===")
        self.logger.info(f"🧬 基因组k-mer总数: {n_kmers} | Total genome k-mers: {n_kmers}")
        self.logger.info(f"🧪 样本总数: {n_samples} | Total samples: {n_samples}")
        
        # 每个k-mer在多少个样本中存在 | How many samples each k-mer is present in
        presence_counts = presence_matrix.sum(axis=1)
        self.logger.info(f"\n🎯 k-mer存在性统计 | K-mer presence statistics:")
        self.logger.info(f"  💯 在所有样本中都存在的k-mer: {np.sum(presence_counts == n_samples)} | K-mers present in all samples: {np.sum(presence_counts == n_samples)}")
        self.logger.info(f"  📈 在50%以上样本中存在的k-mer: {np.sum(presence_counts >= n_samples/2)} | K-mers present in >50% samples: {np.sum(presence_counts >= n_samples/2)}")
        self.logger.info(f"  🚫 在任何样本中都不存在的k-mer: {np.sum(presence_counts == 0)} | K-mers absent in all samples: {np.sum(presence_counts == 0)}")
        
        # 每个样本中有多少基因组k-mer存在 | How many genome k-mers are present in each sample
        sample_counts = presence_matrix.sum(axis=0)
        self.logger.info(f"\n🧪 样本k-mer检出率统计 | Sample k-mer detection statistics:")
        self.logger.info(f"  📊 平均每个样本检出k-mer数: {np.mean(sample_counts):.0f} | Average k-mers detected per sample: {np.mean(sample_counts):.0f}")
        self.logger.info(f"  🔝 最大检出数: {np.max(sample_counts)} | Maximum detected: {np.max(sample_counts)}")
        self.logger.info(f"  🔻 最小检出数: {np.min(sample_counts)} | Minimum detected: {np.min(sample_counts)}")
        
        # 保存样本统计信息 | Save sample statistics
        stats_df = pd.DataFrame({
            'sample': samples,
            'detected_kmers': sample_counts,
            'detection_rate': sample_counts / n_kmers * 100
        })
        
        stats_file = self.config.output_path / f"kmer{self.config.kmer_size}_sample_stats.csv"
        stats_df.to_csv(stats_file, index=False)
        
        self.logger.info(f"\n📁 结果文件 | Result files:")
        self.logger.info(f"  📊 存在性矩阵 | Presence matrix: kmer{self.config.kmer_size}_presence_matrix.csv")
        self.logger.info(f"  📈 样本统计 | Sample statistics: kmer{self.config.kmer_size}_sample_stats.csv")
