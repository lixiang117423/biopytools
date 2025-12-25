"""
GenomeThreader 输出处理模块 | GenomeThreader Output Processing Module
"""

import os
import json
from pathlib import Path
from typing import Dict, List
from .utils import CommandRunner

class OutputParser:
    """输出结果解析器 | Output Result Parser"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        
    def parse_gff3_output(self) -> Dict[str, any]:
        """解析GFF3输出结果 | Parse GFF3 output results"""
        self.logger.info("📊 解析GFF3输出结果 | Parsing GFF3 output results")
        
        gff3_file = self.config.output_path / f"{self.config.base_name}.out"
        if not gff3_file.exists():
            # 尝试查找其他可能的输出文件
            gff3_file = self.config.output_path / f"{self.config.base_name}_final.gff3"
        
        if not gff3_file.exists():
            self.logger.warning("⚠️ 未找到GFF3输出文件 | GFF3 output file not found")
            return {}
        
        stats = {
            'genes': 0,
            'mrnas': 0,
            'exons': 0,
            'cds': 0,
            'chromosomes': set()
        }
        
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue
                    
                    chromosome = parts[0]
                    feature_type = parts[2]
                    
                    stats['chromosomes'].add(chromosome)
                    
                    if feature_type == 'gene':
                        stats['genes'] += 1
                    elif feature_type in ['mRNA', 'transcript']:
                        stats['mrnas'] += 1
                    elif feature_type == 'exon':
                        stats['exons'] += 1
                    elif feature_type == 'CDS':
                        stats['cds'] += 1
            
            stats['chromosomes'] = len(stats['chromosomes'])
            
            self.logger.info(f"📈 GFF3解析统计 | GFF3 parsing statistics:")
            self.logger.info(f"  🧬 基因数 | Genes: {stats['genes']}")
            self.logger.info(f"  🧵 转录本数 | mRNAs: {stats['mrnas']}")
            self.logger.info(f"  📝 外显子数 | Exons: {stats['exons']}")
            self.logger.info(f"  🔤 CDS数 | CDS: {stats['cds']}")
            self.logger.info(f"  🗺️ 染色体数 | Chromosomes: {stats['chromosomes']}")
            
        except Exception as e:
            self.logger.error(f"❌ 解析GFF3文件时出错 | Error parsing GFF3 file: {e}")
            return {}
        
        return stats
    
    def validate_output_quality(self) -> bool:
        """验证输出质量 | Validate output quality"""
        self.logger.info("🔍 验证输出质量 | Validating output quality")
        
        stats = self.parse_gff3_output()
        
        if not stats:
            self.logger.error("❌ 输出验证失败：无法解析结果文件 | Output validation failed: cannot parse result file")
            return False
        
        # 基本质量检查
        if stats['genes'] == 0:
            self.logger.warning("⚠️ 警告：未预测到任何基因 | Warning: No genes predicted")
            return False
        
        if stats['mrnas'] == 0:
            self.logger.warning("⚠️ 警告：未预测到任何转录本 | Warning: No transcripts predicted")
        
        # 结构完整性检查
        if stats['genes'] > 0 and stats['exons'] == 0:
            self.logger.warning("⚠️ 警告：有基因但无外显子结构 | Warning: Genes found but no exon structures")
        
        self.logger.info("✅ 输出质量验证通过 | Output quality validation passed")
        return True

class SummaryGenerator:
    """结果总结生成器 | Results Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_analysis_summary(self, sequence_stats: Dict, output_stats: Dict):
        """生成分析总结报告 | Generate analysis summary report"""
        self.logger.info("📝 生成分析总结报告 | Generating analysis summary report")
        
        summary_file = self.config.output_path / "analysis_summary.txt"
        
        try:
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write("GenomeThreader 基因预测分析总结报告\n")
                f.write("GenomeThreader Gene Prediction Analysis Summary Report\n")
                f.write("=" * 70 + "\n\n")
                
                # 输入文件信息
                f.write("输入文件信息 | Input File Information:\n")
                f.write("-" * 40 + "\n")
                f.write(f"基因组文件 | Genomic file: {self.config.genomic_file}\n")
                if self.config.cdna_file:
                    f.write(f"cDNA文件 | cDNA file: {self.config.cdna_file}\n")
                if self.config.protein_file:
                    f.write(f"蛋白质文件 | Protein file: {self.config.protein_file}\n")
                if self.config.est_file:
                    f.write(f"EST文件 | EST file: {self.config.est_file}\n")
                f.write(f"物种模型 | Species model: {self.config.species or 'undefined'}\n\n")
                
                # 输入序列统计
                f.write("输入序列统计 | Input Sequence Statistics:\n")
                f.write("-" * 40 + "\n")
                if 'genomic' in sequence_stats:
                    genomic = sequence_stats['genomic']
                    f.write(f"基因组序列数 | Genomic sequences: {genomic['sequences']}\n")
                    f.write(f"基因组总长度 | Genomic total length: {genomic['total_length']:,} bp\n")
                
                for seq_type in ['cdna', 'protein', 'est']:
                    if seq_type in sequence_stats:
                        stats = sequence_stats[seq_type]
                        type_name = {'cdna': 'cDNA', 'protein': '蛋白质', 'est': 'EST'}[seq_type]
                        f.write(f"{type_name}序列数 | {type_name.upper()} sequences: {stats['sequences']}\n")
                        f.write(f"{type_name}总长度 | {type_name.upper()} total length: {stats['total_length']:,} bp\n")
                f.write("\n")
                
                # 分析参数
                f.write("分析参数 | Analysis Parameters:\n")
                f.write("-" * 40 + "\n")
                f.write(f"最小比对得分 | Min alignment score: {self.config.min_alignment_score}\n")
                f.write(f"全局链最小覆盖度 | GC min coverage: {self.config.gc_min_coverage}%\n")
                f.write(f"旁系同源基因分析 | Paralogous gene analysis: {'是' if self.config.paralogs else '否'}\n")
                f.write(f"内含子切除技术 | Intron cutout technique: {'是' if self.config.intron_cutout else '否'}\n")
                f.write(f"快速DP算法 | Fast DP algorithm: {'是' if self.config.fastdp else '否'}\n\n")
                
                # 预测结果统计
                if output_stats:
                    f.write("预测结果统计 | Prediction Result Statistics:\n")
                    f.write("-" * 40 + "\n")
                    f.write(f"预测基因数 | Predicted genes: {output_stats.get('genes', 0)}\n")
                    f.write(f"预测转录本数 | Predicted transcripts: {output_stats.get('mrnas', 0)}\n")
                    f.write(f"预测外显子数 | Predicted exons: {output_stats.get('exons', 0)}\n")
                    f.write(f"预测CDS数 | Predicted CDS: {output_stats.get('cds', 0)}\n")
                    f.write(f"涉及染色体数 | Chromosomes involved: {output_stats.get('chromosomes', 0)}\n")
                    
                    if output_stats.get('genes', 0) > 0:
                        avg_transcripts = output_stats.get('mrnas', 0) / output_stats['genes']
                        avg_exons = output_stats.get('exons', 0) / output_stats['genes']
                        f.write(f"平均每个基因转录本数 | Avg transcripts per gene: {avg_transcripts:.2f}\n")
                        f.write(f"平均每个基因外显子数 | Avg exons per gene: {avg_exons:.2f}\n")
                f.write("\n")
                
                # 输出文件
                f.write("输出文件 | Output Files:\n")
                f.write("-" * 40 + "\n")
                output_files = list(self.config.output_path.glob("*.out")) + list(self.config.output_path.glob("*.gff3"))
                for file_path in output_files:
                    f.write(f"  - {file_path.name}\n")
                f.write(f"\n输出目录 | Output directory: {self.config.output_dir}\n")
            
            self.logger.info(f"✅ 总结报告已生成 | Summary report generated: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 生成总结报告时出错 | Error generating summary report: {e}")
    
    def generate_json_summary(self, sequence_stats: Dict, output_stats: Dict):
        """生成JSON格式的总结 | Generate JSON format summary"""
        summary_data = {
            "analysis_type": "GenomeThreader Gene Prediction",
            "timestamp": str(Path().cwd()),  # 简化时间戳
            "input_files": {
                "genomic_file": self.config.genomic_file,
                "cdna_file": self.config.cdna_file,
                "protein_file": self.config.protein_file,
                "est_file": self.config.est_file,
                "species": self.config.species
            },
            "parameters": {
                "min_alignment_score": self.config.min_alignment_score,
                "gc_min_coverage": self.config.gc_min_coverage,
                "paralogs": self.config.paralogs,
                "intron_cutout": self.config.intron_cutout,
                "fastdp": self.config.fastdp
            },
            "sequence_statistics": sequence_stats,
            "output_statistics": output_stats,
            "output_directory": self.config.output_dir
        }
        
        json_file = self.config.output_path / "analysis_summary.json"
        try:
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(summary_data, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"✅ JSON总结已生成 | JSON summary generated: {json_file}")
        except Exception as e:
            self.logger.error(f"❌ 生成JSON总结时出错 | Error generating JSON summary: {e}")
