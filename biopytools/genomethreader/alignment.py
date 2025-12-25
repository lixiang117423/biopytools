"""
GenomeThreader 比对分析模块 | GenomeThreader Alignment Analysis Module
"""

import os
from pathlib import Path
from typing import List, Optional
from .utils import CommandRunner

class GTHAligner:
    """GenomeThreader 比对器 | GenomeThreader Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_gth_command(self) -> str:
        """构建GenomeThreader命令 | Build GenomeThreader command"""
        cmd_parts = [self.config.gth_path]
        
        # 必需参数：基因组文件
        cmd_parts.extend(['-genomic', self.config.genomic_file])
        
        # 序列输入文件
        if self.config.cdna_file:
            cmd_parts.extend(['-cdna', self.config.cdna_file])
        if self.config.protein_file:
            cmd_parts.extend(['-protein', self.config.protein_file])
        if self.config.est_file:
            cmd_parts.extend(['-cdna', self.config.est_file])  # EST files use -cdna parameter
        
        # 物种模型
        if self.config.species:
            cmd_parts.extend(['-species', self.config.species])
        
        # 链方向参数
        if self.config.forward_only:
            cmd_parts.append('-f')
        if self.config.reverse_only:
            cmd_parts.append('-r')
        if self.config.cdna_forward:
            cmd_parts.append('-cdnaforward')
        
        # 比对质量参数
        if self.config.min_alignment_score != 0.65:  # 非默认值
            cmd_parts.extend(['-minalignmentscore', str(self.config.min_alignment_score)])
        if self.config.gc_min_coverage != 50:  # 非默认值
            cmd_parts.extend(['-gcmincoverage', str(self.config.gc_min_coverage)])
        
        # 特殊分析参数
        if self.config.paralogs:
            cmd_parts.append('-paralogs')
        
        # 性能优化参数
        if self.config.intron_cutout:
            cmd_parts.append('-introncutout')
        if self.config.auto_intron_cutout > 0:
            cmd_parts.extend(['-autointroncutout', str(self.config.auto_intron_cutout)])
        if self.config.fastdp:
            cmd_parts.append('-fastdp')
        
        # 高级参数
        if self.config.bssm_file:
            cmd_parts.extend(['-bssm', self.config.bssm_file])
        if self.config.score_matrix != 'BLOSUM62':  # 非默认值
            cmd_parts.extend(['-scorematrix', self.config.score_matrix])
        if self.config.translation_table != 1:  # 非默认值
            cmd_parts.extend(['-translationtable', str(self.config.translation_table)])
        
        # 序列范围参数
        if self.config.from_pos != 1:  # 非默认值
            cmd_parts.extend(['-frompos', str(self.config.from_pos)])
        if self.config.to_pos:
            cmd_parts.extend(['-topos', str(self.config.to_pos)])
        if self.config.width != 1000000:  # 非默认值
            cmd_parts.extend(['-width', str(self.config.width)])
        
        # 输出限制参数
        if self.config.first_alignments > 0:
            cmd_parts.extend(['-first', str(self.config.first_alignments)])
        
        # 输出格式参数
        if self.config.gff3_output:
            cmd_parts.append('-gff3out')
            # GFF3输出需要配合-skipalignmentout或-intermediate使用
            if not self.config.intermediate and not self.config.skip_alignment_out:
                cmd_parts.append('-skipalignmentout')  # 自动添加-skipalignmentout
                self.logger.info("🔧 自动启用-skipalignmentout以支持GFF3输出 | Auto-enabled -skipalignmentout for GFF3 output")
        if self.config.xml_output:
            cmd_parts.append('-xmlout')
        if self.config.intermediate:
            cmd_parts.append('-intermediate')
        if self.config.skip_alignment_out:
            cmd_parts.append('-skipalignmentout')
        
        # 输出文件
        output_file = self.config.output_path / f"{self.config.base_name}.out"
        cmd_parts.extend(['-o', str(output_file)])
        
        return ' '.join(cmd_parts)
    
    def run_alignment(self) -> bool:
        """运行GenomeThreader比对分析 | Run GenomeThreader alignment analysis"""
        self.logger.info("🧬 开始GenomeThreader比对分析 | Starting GenomeThreader alignment analysis")
        
        # 构建命令
        cmd = self.build_gth_command()
        
        # 记录分析参数
        self._log_analysis_parameters()
        
        # 执行比对
        success = self.cmd_runner.run(
            cmd, 
            f"GenomeThreader基因结构预测 | GenomeThreader gene structure prediction"
        )
        
        if success:
            self.logger.info("✅ GenomeThreader比对分析完成 | GenomeThreader alignment analysis completed")
            return True
        else:
            self.logger.error("❌ GenomeThreader比对分析失败 | GenomeThreader alignment analysis failed")
            return False
    
    def _log_analysis_parameters(self):
        """记录分析参数 | Log analysis parameters"""
        self.logger.info("📋 分析参数总结 | Analysis Parameters Summary")
        self.logger.info(f"  🧬 基因组文件 | Genomic file: {self.config.genomic_file}")
        
        if self.config.cdna_file:
            self.logger.info(f"  🧵 cDNA文件 | cDNA file: {self.config.cdna_file}")
        if self.config.protein_file:
            self.logger.info(f"  🔤 蛋白质文件 | Protein file: {self.config.protein_file}")
        if self.config.est_file:
            self.logger.info(f"  📝 EST文件 | EST file: {self.config.est_file}")
        
        if self.config.species:
            self.logger.info(f"  🐾 物种模型 | Species model: {self.config.species}")
        
        self.logger.info(f"  ⚡ 线程数 | Threads: {self.config.threads}")
        self.logger.info(f"  🎯 最小比对得分 | Min alignment score: {self.config.min_alignment_score}")
        self.logger.info(f"  📊 全局链最小覆盖度 | GC min coverage: {self.config.gc_min_coverage}%")
        
        if self.config.paralogs:
            self.logger.info("  🔄 启用旁系同源基因分析 | Paralogous gene analysis enabled")
        if self.config.intron_cutout:
            self.logger.info("  ✂️ 启用内含子切除技术 | Intron cutout technique enabled")
        if self.config.fastdp:
            self.logger.info("  🚀 启用快速DP算法 | Fast DP algorithm enabled")

class ConsensusProcessor:
    """共识序列处理器 | Consensus Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def process_intermediate_results(self) -> bool:
        """处理中间结果生成最终预测 | Process intermediate results to generate final predictions"""
        if not self.config.intermediate:
            return True
        
        self.logger.info("🔄 处理中间结果 | Processing intermediate results")
        
        intermediate_file = self.config.output_path / f"{self.config.base_name}.out"
        final_output = self.config.output_path / f"{self.config.base_name}_final.gff3"
        
        # 构建gthconsensus命令
        cmd_parts = [
            self.config.gthconsensus_path,
            '-intermediate',
            '-gff3out',
            str(intermediate_file)
        ]
        
        cmd = f"{' '.join(cmd_parts)} > {final_output}"
        
        success = self.cmd_runner.run(
            cmd,
            "处理中间结果生成最终GFF3 | Process intermediate results to final GFF3"
        )
        
        if success:
            self.logger.info(f"✅ 最终预测结果已生成 | Final predictions generated: {final_output}")
        else:
            self.logger.error("❌ 中间结果处理失败 | Failed to process intermediate results")
        
        return success
