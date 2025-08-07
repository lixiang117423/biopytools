"""
基因组比对模块 | Genome Alignment Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class GenomeAligner:
    """基因组比对器 | Genome Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 创建比对结果目录 | Create alignment results directory
        self.alignment_dir = Path(self.config.output_dir) / "alignments"
        self.alignment_dir.mkdir(parents=True, exist_ok=True)
    
    def run_alignments(self) -> dict:
        """运行所有成对基因组比对 | Run all pairwise genome alignments"""
        self.logger.info("🧬 开始基因组比对 | Starting genome alignments")
        
        alignment_files = {}
        
        for i in range(len(self.config.sample_list) - 1):
            ref_sample = self.config.sample_list[i]
            query_sample = self.config.sample_list[i + 1]
            
            ref_genome = self.config.genome_paths[ref_sample]
            query_genome = self.config.genome_paths[query_sample]
            
            # 生成输出文件名 | Generate output filename
            bam_file = self.alignment_dir / f"{ref_sample}_{query_sample}.bam"
            
            # 运行minimap2比对 | Run minimap2 alignment
            success = self._run_minimap2_alignment(
                ref_genome, query_genome, bam_file, ref_sample, query_sample
            )
            
            if success:
                alignment_files[f"{ref_sample}_{query_sample}"] = str(bam_file)
                
                # 创建BAM索引 | Create BAM index
                self._index_bam_file(bam_file)
            else:
                self.logger.error(f"❌ 比对失败 | Alignment failed: {ref_sample} vs {query_sample}")
                return {}
        
        self.logger.info(f"🎉 完成 {len(alignment_files)} 个基因组比对 | Completed {len(alignment_files)} genome alignments")
        return alignment_files
    
    def _run_minimap2_alignment(self, ref_genome: str, query_genome: str, 
                               output_bam: Path, ref_name: str, query_name: str) -> bool:
        """运行minimap2比对单个基因组对 | Run minimap2 alignment for single genome pair"""
        
        self.logger.info(f"🔄 比对 {ref_name} vs {query_name} | Aligning {ref_name} vs {query_name}")
        
        # 构建minimap2命令 | Build minimap2 command
        cmd = (
            f"{self.config.minimap2_path} -ax {self.config.minimap2_preset} "
            f"-t {self.config.threads} --eqx {ref_genome} {query_genome} | "
            f"{self.config.samtools_path} sort -O BAM -@ {self.config.threads} - > {output_bam}"
        )
        
        description = f"🧬 Minimap2 alignment: {ref_name} vs {query_name}"
        return self.cmd_runner.run(cmd, description)
    
    def _index_bam_file(self, bam_file: Path) -> bool:
        """为BAM文件创建索引 | Create index for BAM file"""
        
        cmd = f"{self.config.samtools_path} index {bam_file}"
        description = f"📇 Indexing BAM file: {bam_file.name}"
        
        return self.cmd_runner.run(cmd, description)
