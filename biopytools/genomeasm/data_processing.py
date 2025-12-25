"""
基因组组装数据处理模块 | Genome Assembly Data Processing Module 📊
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
from .utils import CommandRunner, get_file_stats, estimate_genome_coverage

class DataQualityController:
    """数据质量控制器 | Data Quality Controller"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.qc_dir = Path(config.output_dir) / "qc"
        self.qc_dir.mkdir(parents=True, exist_ok=True)
    
    def run_quality_control(self) -> Dict[str, any]:
        """运行完整的质量控制流程 | Run complete quality control pipeline"""
        self.logger.info("🔬 开始数据质量控制 | Starting data quality control")
        
        qc_results = {}
        
        # HiFi数据质控
        qc_results['hifi'] = self._qc_hifi_data()
        
        # Hi-C数据质控
        if 'Hi-C' in self.config.detected_data_types:
            qc_results['hic'] = self._qc_hic_data()
        
        # ONT数据质控
        if 'ONT' in self.config.detected_data_types:
            qc_results['ont'] = self._qc_ont_data()
        
        # NGS数据质控
        if 'NGS' in self.config.detected_data_types:
            qc_results['ngs'] = self._qc_ngs_data()
        
        # 生成质控报告
        self._generate_qc_report(qc_results)
        
        # 评估数据充分性
        self._assess_data_sufficiency(qc_results)
        
        self.logger.info("✅ 数据质量控制完成 | Data quality control completed")
        return qc_results
    
    def _qc_hifi_data(self) -> Dict[str, any]:
        """HiFi数据质量控制 | HiFi data quality control"""
        self.logger.info("🧬 检查HiFi数据质量 | Checking HiFi data quality")
        
        # 基本统计
        stats = get_file_stats(self.config.hifi_reads, self.logger)
        
        # 估算覆盖度
        coverage = estimate_genome_coverage(
            self.config.hifi_reads, 
            self.config.genome_size, 
            self.logger
        )
        
        # 读长分布分析
        read_lengths = self._analyze_read_lengths(self.config.hifi_reads)
        
        results = {
            'file': self.config.hifi_reads,
            'stats': stats,
            'coverage': coverage,
            'read_lengths': read_lengths,
            'pass_qc': coverage >= self.config.min_hifi_coverage
        }
        
        if not results['pass_qc']:
            self.logger.warning(f"⚠️ HiFi覆盖度偏低 | Low HiFi coverage: {coverage:.1f}X < {self.config.min_hifi_coverage}X")
        else:
            self.logger.info(f"✅ HiFi数据质量良好 | HiFi data quality good: {coverage:.1f}X coverage")
        
        return results
    
    # def _qc_hic_data(self) -> Dict[str, any]:
    #     """Hi-C数据质量控制 | Hi-C data quality control"""
    #     self.logger.info("🔗 检查Hi-C数据质量 | Checking Hi-C data quality")
        
    #     results = {}
        
    #     # 检查R1和R2文件
    #     for read_type, file_path in [('R1', self.config.hic_r1), ('R2', self.config.hic_r2)]:
    #         # FastQC质量检查
    #         fastqc_output = self.qc_dir / f"hic_{read_type.lower()}_fastqc"
    #         cmd = f"fastqc -t {min(self.config.threads, 50)} -o {self.qc_dir} {file_path}"
    #         self.cmd_runner.run(cmd, f"Hi-C {read_type} FastQC质量检查")
            
    #         # 基本统计
    #         stats = get_file_stats(file_path, self.logger)
    #         results[read_type.lower()] = {
    #             'file': file_path,
    #             'stats': stats,
    #             'fastqc_dir': fastqc_output
    #         }

    def _qc_hic_data(self) -> Dict[str, any]:
        """Hi-C数据质量控制 | Hi-C data quality control"""
        self.logger.info("🔗 检查Hi-C数据质量 | Checking Hi-C data quality")
        
        results = {}
        
        # 检查R1和R2文件
        for read_type, file_path in [('R1', self.config.hic_r1), ('R2', self.config.hic_r2)]:
            # 条件性运行FastQC质量检查
            if not self.config.skip_fastqc:
                self.logger.info(f"🔬 运行FastQC质量检查: Hi-C {read_type}")
                fastqc_output = self.qc_dir / f"hic_{read_type.lower()}_fastqc"
                cmd = f"fastqc -t {min(self.config.threads, 8)} -o {self.qc_dir} {file_path}"
                self.cmd_runner.run(cmd, f"Hi-C {read_type} FastQC质量检查")
                
                results[read_type.lower()] = {
                    'file': file_path,
                    'fastqc_dir': fastqc_output
                }
            else:
                self.logger.info(f"⏭️ 跳过FastQC检查: Hi-C {read_type} (节省时间)")
                results[read_type.lower()] = {
                    'file': file_path,
                    'fastqc_skipped': True
                }
            
            # 基本统计（总是运行，速度很快）
            stats = get_file_stats(file_path, self.logger)
            results[read_type.lower()]['stats'] = stats
        
        # 估算Hi-C覆盖度
        if results.get('r1', {}).get('stats', {}).get('sum_len'):
            total_bases = int(results['r1']['stats']['sum_len'].replace(',', '')) * 2  # R1+R2
            genome_size_bp = self._parse_genome_size_bp(self.config.genome_size)
            hic_coverage = total_bases / genome_size_bp if genome_size_bp > 0 else 0
            
            results['coverage'] = hic_coverage
            results['pass_qc'] = hic_coverage >= self.config.min_hic_coverage
            
            if not results['pass_qc']:
                self.logger.warning(f"⚠️ Hi-C覆盖度偏低 | Low Hi-C coverage: {hic_coverage:.1f}X")
            else:
                self.logger.info(f"✅ Hi-C数据质量良好 | Hi-C data quality good: {hic_coverage:.1f}X coverage")
        
        return results
    
    def _qc_ont_data(self) -> Dict[str, any]:
        """ONT数据质量控制 | ONT data quality control"""
        self.logger.info("🔬 检查ONT数据质量 | Checking ONT data quality")
        
        # 基本统计
        stats = get_file_stats(self.config.ont_reads, self.logger)
        
        # 估算覆盖度
        coverage = estimate_genome_coverage(
            self.config.ont_reads,
            self.config.genome_size,
            self.logger
        )
        
        # 读长分布分析
        read_lengths = self._analyze_read_lengths(self.config.ont_reads)
        
        results = {
            'file': self.config.ont_reads,
            'stats': stats,
            'coverage': coverage,
            'read_lengths': read_lengths,
            'pass_qc': coverage >= 20  # ONT minimum coverage
        }
        
        if not results['pass_qc']:
            self.logger.warning(f"⚠️ ONT覆盖度偏低 | Low ONT coverage: {coverage:.1f}X")
        else:
            self.logger.info(f"✅ ONT数据质量良好 | ONT data quality good: {coverage:.1f}X coverage")
        
        return results
    
    def _qc_ngs_data(self) -> Dict[str, any]:
        """NGS数据质量控制 | NGS data quality control"""
        self.logger.info("🧪 检查NGS数据质量 | Checking NGS data quality")
        
        results = {}
        
        # 检查R1和R2文件
        for read_type, file_path in [('R1', self.config.ngs_r1), ('R2', self.config.ngs_r2)]:
            if file_path:
                # 条件性运行FastQC质量检查
                if not self.config.skip_fastqc:
                    self.logger.info(f"🔬 运行FastQC质量检查: NGS {read_type}")
                    cmd = f"fastqc -t {min(self.config.threads, 8)} -o {self.qc_dir} {file_path}"
                    self.cmd_runner.run(cmd, f"NGS {read_type} FastQC质量检查")
                else:
                    self.logger.info(f"⏭️ 跳过FastQC检查: NGS {read_type} (节省时间)")
                
                # 基本统计（总是运行）
                stats = get_file_stats(file_path, self.logger)
                results[read_type.lower()] = {
                    'file': file_path,
                    'stats': stats,
                    'fastqc_skipped': self.config.skip_fastqc
                }
        
        return results
    
    def _analyze_read_lengths(self, reads_file: str) -> Dict[str, float]:
        """分析读长分布 | Analyze read length distribution"""
        try:
            cmd = f"seqkit fx2tab -l {reads_file} | cut -f2 | sort -n"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            lengths = [int(x) for x in result.stdout.strip().split('\n') if x.strip()]
            
            if lengths:
                lengths.sort()
                n = len(lengths)
                return {
                    'min': lengths[0],
                    'max': lengths[-1],
                    'mean': sum(lengths) / n,
                    'median': lengths[n//2],
                    'n50': self._calculate_n50(lengths)
                }
            else:
                return {}
        except Exception as e:
            self.logger.warning(f"⚠️ 读长分析失败 | Read length analysis failed: {e}")
            return {}
    
    def _calculate_n50(self, lengths: List[int]) -> int:
        """计算N50 | Calculate N50"""
        lengths.sort(reverse=True)
        total = sum(lengths)
        target = total * 0.5
        cumulative = 0
        
        for length in lengths:
            cumulative += length
            if cumulative >= target:
                return length
        return 0
    
    def _parse_genome_size_bp(self, genome_size: str) -> int:
        """解析基因组大小为bp | Parse genome size to bp"""
        try:
            size = genome_size.lower()
            if size.endswith('k'):
                return int(float(size[:-1]) * 1000)
            elif size.endswith('m'):
                return int(float(size[:-1]) * 1000000)
            elif size.endswith('g'):
                return int(float(size[:-1]) * 1000000000)
            else:
                return int(size)
        except:
            return 0
    
    def _generate_qc_report(self, qc_results: Dict[str, any]):
        """生成质控报告 | Generate QC report"""
        report_file = self.qc_dir / "quality_control_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🔬 基因组数据质量控制报告 | Genome Data Quality Control Report\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"项目名称 | Project: {self.config.project_name}\n")
            f.write(f"检测到的数据类型 | Detected data types: {', '.join(self.config.detected_data_types)}\n")
            f.write(f"组装策略 | Assembly strategy: {self.config.assembly_strategy}\n")
            f.write(f"FastQC检查 | FastQC check: {'跳过 (节省时间)' if self.config.skip_fastqc else '已运行'}\n\n")
            
            # HiFi数据报告
            if 'hifi' in qc_results:
                hifi = qc_results['hifi']
                f.write("🧬 HiFi数据质量 | HiFi Data Quality:\n")
                f.write(f"  文件 | File: {hifi['file']}\n")
                f.write(f"  覆盖度 | Coverage: {hifi['coverage']:.1f}X\n")
                f.write(f"  质量评估 | QC Pass: {'✅ 通过' if hifi['pass_qc'] else '❌ 不通过'}\n\n")
            
            # Hi-C数据报告
            if 'hic' in qc_results:
                hic = qc_results['hic']
                f.write("🔗 Hi-C数据质量 | Hi-C Data Quality:\n")
                if 'coverage' in hic:
                    f.write(f"  覆盖度 | Coverage: {hic['coverage']:.1f}X\n")
                    f.write(f"  质量评估 | QC Pass: {'✅ 通过' if hic['pass_qc'] else '❌ 不通过'}\n")
                f.write(f"  R1文件 | R1 file: {hic.get('r1', {}).get('file', 'N/A')}\n")
                f.write(f"  R2文件 | R2 file: {hic.get('r2', {}).get('file', 'N/A')}\n\n")
            
            # 其他数据类型...
            
        self.logger.info(f"📋 质控报告已生成 | QC report generated: {report_file}")
    
    def _assess_data_sufficiency(self, qc_results: Dict[str, any]):
        """评估数据充分性 | Assess data sufficiency"""
        issues = []
        recommendations = []
        
        # HiFi数据评估
        if 'hifi' in qc_results:
            hifi = qc_results['hifi']
            if hifi['coverage'] < self.config.min_hifi_coverage:
                issues.append(f"HiFi覆盖度不足: {hifi['coverage']:.1f}X < {self.config.min_hifi_coverage}X")
                recommendations.append("考虑增加HiFi测序深度或降低质量要求")
        
        # Hi-C数据评估
        if 'hic' in qc_results:
            hic = qc_results['hic']
            if 'coverage' in hic and hic['coverage'] < self.config.min_hic_coverage:
                issues.append(f"Hi-C覆盖度不足: {hic['coverage']:.1f}X < {self.config.min_hic_coverage}X")
                recommendations.append("考虑增加Hi-C测序深度")
        
        if issues:
            self.logger.warning("⚠️ 数据质量问题 | Data quality issues:")
            for issue in issues:
                self.logger.warning(f"  - {issue}")
            self.logger.info("💡 建议 | Recommendations:")
            for rec in recommendations:
                self.logger.info(f"  - {rec}")
        else:
            self.logger.info("🎉 所有数据质量检查通过 | All data quality checks passed")
