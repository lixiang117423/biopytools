"""
📈 覆盖度统计分析模块 | Depth Statistics Analysis Module
"""

from pathlib import Path
from collections import defaultdict
import statistics

class DepthStatistics:
    """📈 覆盖度统计分析器 | Depth Statistics Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_statistics(self):
        """📊 生成统计报告 | Generate statistics report"""
        self.logger.info("📈 生成覆盖度统计报告 | Generating depth statistics report")
        
        # 读取结果文件 | Read results file
        sample_stats = defaultdict(list)
        chromosome_stats = defaultdict(list)
        total_positions = 0
        
        try:
            with open(self.config.output_file, 'r') as f:
                header = f.readline().strip()  # 跳过表头 | Skip header
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        sample, chrom, pos, depth = parts[0], parts[1], parts[2], parts[3]
                        
                        try:
                            depth_value = float(depth)
                            sample_stats[sample].append(depth_value)
                            chromosome_stats[chrom].append(depth_value)
                            total_positions += 1
                        except ValueError:
                            continue
            
            # 生成统计报告 | Generate statistics report
            stats_file = Path(self.config.output_file).with_suffix('.stats.txt')
            
            with open(stats_file, 'w') as f:
                f.write("📈 覆盖度统计报告 | Depth Statistics Report\n")
                f.write("=" * 60 + "\n\n")
                
                # 总体统计 | Overall statistics
                f.write("🔢 总体统计 | Overall Statistics:\n")
                f.write(f"  📍 总位点数 | Total positions: {total_positions:,}\n")
                f.write(f"  🏷️ 样品数量 | Number of samples: {len(sample_stats)}\n")
                f.write(f"  🧬 染色体数量 | Number of chromosomes: {len(chromosome_stats)}\n\n")
                
                # 样品统计 | Sample statistics
                f.write("🏷️ 各样品统计 | Sample Statistics:\n")
                for sample, depths in sample_stats.items():
                    if depths:
                        f.write(f"  📄 {sample}:\n")
                        f.write(f"    位点数 | Positions: {len(depths):,}\n")
                        f.write(f"    平均覆盖度 | Mean depth: {statistics.mean(depths):.2f}\n")
                        f.write(f"    中位覆盖度 | Median depth: {statistics.median(depths):.2f}\n")
                        f.write(f"    最小覆盖度 | Min depth: {min(depths):.2f}\n")
                        f.write(f"    最大覆盖度 | Max depth: {max(depths):.2f}\n")
                        if len(depths) > 1:
                            f.write(f"    标准差 | Std dev: {statistics.stdev(depths):.2f}\n")
                        f.write("\n")
                
                # 染色体统计 | Chromosome statistics
                f.write("🧬 各染色体统计 | Chromosome Statistics:\n")
                for chrom, depths in chromosome_stats.items():
                    if depths:
                        f.write(f"  📍 {chrom}:\n")
                        f.write(f"    位点数 | Positions: {len(depths):,}\n")
                        f.write(f"    平均覆盖度 | Mean depth: {statistics.mean(depths):.2f}\n")
                        f.write(f"    中位覆盖度 | Median depth: {statistics.median(depths):.2f}\n")
                        f.write(f"    最小覆盖度 | Min depth: {min(depths):.2f}\n")
                        f.write(f"    最大覆盖度 | Max depth: {max(depths):.2f}\n")
                        f.write("\n")
                
                # 输出文件信息 | Output file information
                f.write("📁 输出文件 | Output Files:\n")
                f.write(f"  📄 主结果文件 | Main results file: {self.config.output_file}\n")
                f.write(f"  📊 统计报告文件 | Statistics report file: {stats_file}\n")
                
                # 检查是否有窗口分析文件 | Check if window analysis file exists
                if self.config.enable_window_analysis:
                    window_file = Path(self.config.output_file).with_name(f"{Path(self.config.output_file).stem}_windows_{self.config.window_size}bp{Path(self.config.output_file).suffix}")
                    if window_file.exists():
                        f.write(f"  📊 窗口分析文件 | Window analysis file: {window_file}\n")
                        f.write(f"    📏 窗口大小 | Window size: {self.config.window_size}bp\n")
                        f.write(f"    👣 步长 | Step size: {self.config.window_step}bp\n")
                        if self.config.window_step < self.config.window_size:
                            overlap = self.config.window_size - self.config.window_step
                            f.write(f"    🔄 重叠区域 | Overlap: {overlap}bp ({overlap/self.config.window_size*100:.1f}%)\n")
                        else:
                            f.write(f"    📍 无重叠 | No overlap\n")
            
            self.logger.info(f"✅ 统计报告已生成: {stats_file} | Statistics report generated: {stats_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 生成统计报告失败: {e} | Failed to generate statistics report: {e}")
            return False
