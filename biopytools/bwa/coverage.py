"""
覆盖度分析模块|Coverage Analysis Module
"""

import subprocess
from pathlib import Path

class CoverageAnalyzer:
    """覆盖度分析器|Coverage Analyzer"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def analyze_sample_coverage(self, sample_name: str, bam_file: str):
        """分析样品覆盖度|Analyze sample coverage"""
        self.logger.info(f"开始覆盖度分析|Starting coverage analysis: {sample_name}")
        
        # Step 1: 计算位点覆盖度|Calculate per-base coverage
        depth_file = self._calculate_depth(sample_name, bam_file)
        if not depth_file:
            return False
        
        # Step 2: 计算滑窗覆盖度|Calculate windowed coverage
        if not self._calculate_windowed_coverage(sample_name, depth_file):
            return False
        
        self.logger.info(f"覆盖度分析完成|Coverage analysis completed: {sample_name}")
        return True
    
    def _calculate_depth(self, sample_name: str, bam_file: str) -> Path:
        """计算位点覆盖度|Calculate per-base depth"""
        output_file = self.config.coverage_dir / f"{sample_name}.depth.txt"
        
        # 构建samtools depth命令|Build samtools depth command
        cmd_parts = [f"{self.config.samtools_path} depth"]
        
        if self.config.min_base_quality > 0:
            cmd_parts.append(f"-q {self.config.min_base_quality}")
        
        if self.config.min_mapping_quality > 0:
            cmd_parts.append(f"-Q {self.config.min_mapping_quality}")
        
        if self.config.max_depth > 0:
            cmd_parts.append(f"-d {self.config.max_depth}")

        cmd_parts.append(f"{Path(bam_file).absolute()} > {output_file.absolute()}")

        cmd = " ".join(cmd_parts)
        
        success = self.cmd_runner.run(cmd, f"计算位点覆盖度|Calculating depth: {sample_name}")
        
        return output_file if success else None
    
    def _calculate_windowed_coverage(self, sample_name: str, depth_file: Path) -> bool:
        """计算滑窗覆盖度|Calculate windowed coverage"""
        output_file = self.config.window_dir / f"{sample_name}.window.bed"
        
        self.logger.info(f"计算滑窗覆盖度|Calculating windowed coverage")
        self.logger.info(f"   窗口大小|Window size: {self.config.window_size:,} bp")
        self.logger.info(f"   步长|Step size: {self.config.step_size:,} bp")
        
        try:
            # 读取深度文件并计算滑窗|Read depth file and calculate windows
            windows = {}
            
            with open(depth_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    
                    chrom = parts[0]
                    pos = int(parts[1])
                    depth = int(parts[2])
                    
                    # 计算窗口索引|Calculate window index
                    window_idx = (pos - 1) // self.config.step_size
                    window_start = window_idx * self.config.step_size + 1
                    window_end = window_start + self.config.window_size - 1
                    
                    window_key = (chrom, window_start, window_end)
                    
                    if window_key not in windows:
                        windows[window_key] = {'sum': 0, 'count': 0}
                    
                    windows[window_key]['sum'] += depth
                    windows[window_key]['count'] += 1
            
            # 写入结果|Write results
            with open(output_file, 'w') as f:
                for (chrom, start, end), data in sorted(windows.items()):
                    avg_depth = data['sum'] / data['count'] if data['count'] > 0 else 0
                    f.write(f"{chrom}\t{start}\t{end}\t{avg_depth:.2f}\n")
            
            self.logger.info(f"滑窗覆盖度计算完成|Windowed coverage completed")
            self.logger.info(f"   输出文件|Output file: {output_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f" 滑窗覆盖度计算失败|Windowed coverage failed: {e}")
            return False
