"""
VCF距离矩阵计算模块 | VCF Distance Matrix Calculation Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class DistanceCalculator:
    """距离矩阵计算器 | Distance Matrix Calculator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def calculate_distance_matrix(self) -> bool:
        """计算距离矩阵 | Calculate distance matrix"""
        if self.config.skip_vcf2dis:
            self.logger.info("跳过VCF2Dis步骤 | Skipping VCF2Dis step")
            return True
        
        self.logger.info("开始计算VCF距离矩阵 | Starting VCF distance matrix calculation")
        
        # 获取VCF统计信息 | Get VCF statistics
        stats = self._get_vcf_stats()
        if stats:
            self.logger.info(f"VCF文件包含 {stats.get('samples', 'unknown')} 个样本和 {stats.get('variants', 'unknown')} 个变异")
        
        # 构建VCF2Dis命令 | Build VCF2Dis command
        cmd = (
            f"{self.config.vcf2dis_path} "
            f"-InPut {self.config.vcf_file} "
            f"-OutPut {self.config.distance_matrix}"
        )
        
        success = self.cmd_runner.run(
            cmd, 
            "使用VCF2Dis计算距离矩阵 | Calculate distance matrix using VCF2Dis"
        )
        
        if success:
            self.logger.info(f"距离矩阵已生成 | Distance matrix generated: {self.config.distance_matrix}")
            self._log_matrix_stats()
        
        return success
    
    def _get_vcf_stats(self):
        """获取VCF文件统计信息 | Get VCF file statistics"""
        try:
            stats = {}
            
            # 统计样本数量 | Count samples
            with open(self.config.vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#CHROM'):
                        headers = line.strip().split('\t')
                        # 前9列是固定的，从第10列开始是样本 | First 9 columns are fixed, samples start from column 10
                        stats['samples'] = len(headers) - 9
                        break
            
            # 统计变异数量 | Count variants
            variant_count = 0
            with open(self.config.vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        variant_count += 1
            stats['variants'] = variant_count
            
            return stats
            
        except Exception as e:
            self.logger.warning(f"无法获取VCF统计信息 | Cannot get VCF statistics: {e}")
            return {}
    
    def _log_matrix_stats(self):
        """记录矩阵统计信息 | Log matrix statistics"""
        try:
            if os.path.exists(self.config.distance_matrix):
                with open(self.config.distance_matrix, 'r') as f:
                    lines = f.readlines()
                    
                # 跳过可能的标题行 | Skip possible header line
                data_lines = [line for line in lines if line.strip() and not line.startswith('#')]
                
                # 第一行可能是矩阵维度 | First line might be matrix dimension
                first_line = data_lines[0].strip() if data_lines else ""
                try:
                    sample_count = int(first_line)
                    self.logger.info(f"距离矩阵包含 {sample_count} 个样本 | Distance matrix contains {sample_count} samples")
                except ValueError:
                    # 如果第一行不是数字，说明是数据行 | If first line is not a number, it's a data line
                    sample_count = len(data_lines)
                    self.logger.info(f"距离矩阵包含约 {sample_count} 个样本 | Distance matrix contains approximately {sample_count} samples")
                    
        except Exception as e:
            self.logger.warning(f"无法读取距离矩阵统计信息 | Cannot read distance matrix statistics: {e}")
