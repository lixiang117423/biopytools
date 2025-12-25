"""
🧩 转录本组装模块 | Transcript Assembly Module
"""

from pathlib import Path
from .utils import CommandRunner, SequencingTypeDetector

class StringTieAssembler:
    """🧩 StringTie转录本重构器 | StringTie Transcript Reconstructor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def reconstruct_transcripts(self):
        """🔧 重构转录本 | Reconstruct transcripts"""
        if not hasattr(self.config, 'bam_files') or not self.config.bam_files:
            self.logger.error("❌ 未找到BAM文件，请先执行比对步骤 | No BAM files found, please run alignment step first")
            return False
        
        output_dir = self.config.output_path / "stringtie_assembly"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.config.gtf_files = []
        
        # 第一轮：重构每个样本的转录本 | First round: reconstruct transcripts for each sample
        for i, bam_file in enumerate(self.config.bam_files):
            sample_name = f"sample_{i+1:03d}"
            gtf_file = output_dir / f"{sample_name}.gtf"
            
            self.logger.info(f"🧩 重构转录本 | Reconstructing transcripts: {sample_name}")
            
            # 构建StringTie命令 | Build StringTie command
            cmd = (
                f"{self.config.stringtie_path} "
                f"{bam_file} "
                f"-o {gtf_file} "
                f"-p {self.config.threads} "
                f"-m {self.config.stringtie_min_length} "
                f"-c {self.config.stringtie_min_coverage} "
                f"-f {self.config.stringtie_min_iso} "
                f"-l {sample_name}"
            )
            
            # 添加可选参数 | Add optional parameters
            if self.config.stringtie_conservative:
                cmd += " --conservative"
            
            if not self.cmd_runner.run(cmd, f"🧩 StringTie转录本重构 | StringTie transcript reconstruction: {sample_name}"):
                return False
            
            self.config.gtf_files.append(str(gtf_file))
        
        # 第二轮：合并转录本 | Second round: merge transcripts
        if len(self.config.gtf_files) > 1:
            if not self._merge_transcripts(output_dir):
                return False
        else:
            self.config.merged_gtf = self.config.gtf_files[0]
        
        self.logger.info(f"🎉 StringTie转录本重构完成 | StringTie transcript reconstruction completed")
        return True
    
    def _merge_transcripts(self, output_dir: Path) -> bool:
        """🔀 合并多个样本的转录本 | Merge transcripts from multiple samples"""
        gtf_list_file = output_dir / "gtf_list.txt"
        merged_gtf = output_dir / "merged.gtf"
        
        # 创建GTF文件列表 | Create GTF file list
        with open(gtf_list_file, 'w') as f:
            for gtf_file in self.config.gtf_files:
                f.write(f"{gtf_file}\n")
        
        self.logger.info("🔀 合并转录本 | Merging transcripts")
        
        # 构建StringTie merge命令 | Build StringTie merge command
        cmd = (
            f"{self.config.stringtie_path} --merge "
            f"{gtf_list_file} "
            f"-o {merged_gtf} "
            f"-m {self.config.stringtie_min_length} "
            f"-c {self.config.stringtie_min_coverage} "
            f"-F {self.config.stringtie_min_fpkm}"
        )
        
        if not self.cmd_runner.run(cmd, "🔀 StringTie转录本合并 | StringTie transcript merging"):
            return False
        
        self.config.merged_gtf = str(merged_gtf)
        return True

class TrinityAssembler:
    """🔗 Trinity de novo组装器 | Trinity de novo Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def de_novo_assembly(self):
        """🔬 de novo转录本组装 | de novo transcript assembly"""
        output_dir = self.config.output_path / "trinity_assembly"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.info("🔬 开始Trinity de novo组装 | Starting Trinity de novo assembly")
        
        # 检查是否有样本文件 | Check if samples file exists
        if self.config.samples_file:
            return self._assembly_with_samples_file(output_dir)
        else:
            return self._assembly_with_fastq_files(output_dir)
    
    def _assembly_with_samples_file(self, output_dir: Path) -> bool:
        """📁 使用样本文件进行组装 | Assembly with samples file"""
        cmd = (
            f"{self.config.trinity_path} "
            f"--seqType fq "
            f"--samples_file {self.config.samples_file} "
            f"--max_memory {self.config.trinity_max_memory} "
            f"--CPU {self.config.trinity_cpu} "
            f"--output {output_dir} "
            f"--min_contig_length {self.config.trinity_min_contig_length}"
        )
        
        # 添加链特异性参数 | Add strand specificity parameter
        if self.config.trinity_ss_lib_type:
            cmd += f" --SS_lib_type {self.config.trinity_ss_lib_type}"
        
        if not self.cmd_runner.run(cmd, "📁 Trinity de novo组装（样本文件模式） | Trinity de novo assembly (samples file mode)"):
            return False
        
        return self._check_trinity_output(output_dir)
    
    def _assembly_with_fastq_files(self, output_dir: Path) -> bool:
        """📊 使用FASTQ文件进行组装 | Assembly with FASTQ files"""
        # 检测测序类型和样本组织 | Detect sequencing type and sample organization
        layout, samples = SequencingTypeDetector.detect_sequencing_layout(self.config.rna_seq_files)
        
        # 构建Trinity命令 | Build Trinity command
        cmd = (
            f"{self.config.trinity_path} "
            f"--seqType fq "
            f"--max_memory {self.config.trinity_max_memory} "
            f"--CPU {self.config.trinity_cpu} "
            f"--output {output_dir} "
            f"--min_contig_length {self.config.trinity_min_contig_length}"
        )
        
        # 处理输入文件 | Process input files
        if layout == 'paired':
            # 配对末端 | Paired-end
            left_files = [sample[0] for sample in samples]
            right_files = [sample[1] for sample in samples]
            cmd += f" --left {','.join(left_files)} --right {','.join(right_files)}"
        elif layout == 'single':
            # 单端 | Single-end
            single_files = [sample[0] for sample in samples]
            cmd += f" --single {','.join(single_files)}"
        elif layout == 'mixed':
            # 混合模式：分别处理配对和单端 | Mixed mode: handle paired and single separately
            paired_samples = [s for s in samples if len(s) == 2]
            single_samples = [s for s in samples if len(s) == 1]
            
            if paired_samples:
                left_files = [sample[0] for sample in paired_samples]
                right_files = [sample[1] for sample in paired_samples]
                cmd += f" --left {','.join(left_files)} --right {','.join(right_files)}"
            
            if single_samples:
                single_files = [sample[0] for sample in single_samples]
                if paired_samples:
                    # 如果已有配对文件，单端文件作为附加 | If paired files exist, single files as additional
                    self.logger.warning("⚠️ 混合模式检测到，建议分别处理配对和单端数据 | Mixed mode detected, recommend processing paired and single data separately")
                cmd += f" --single {','.join(single_files)}"
        
        # 添加链特异性参数 | Add strand specificity parameter
        if self.config.trinity_ss_lib_type:
            cmd += f" --SS_lib_type {self.config.trinity_ss_lib_type}"
        
        if not self.cmd_runner.run(cmd, "📊 Trinity de novo组装 | Trinity de novo assembly"):
            return False
        
        return self._check_trinity_output(output_dir)
    
    # def _check_trinity_output(self, output_dir: Path) -> bool:
    #     """✅ 检查Trinity输出 | Check Trinity output"""
    #     # Trinity输出文件路径 | Trinity output file path
    #     trinity_fasta = output_dir / "Trinity.fasta"
        
    #     if trinity_fasta.exists():
    #         self.config.trinity_assembly = str(trinity_fasta)
    #         self.logger.info(f"🎉 Trinity组装完成 | Trinity assembly completed: {trinity_fasta}")
    #         return True
    #     else:
    #         self.logger.error("❌ Trinity组装失败，未找到输出文件 | Trinity assembly failed, output file not found")
    #         return False
    def _check_trinity_output(self, output_dir: Path) -> bool:
        """✅ 检查Trinity输出 | Check Trinity output"""
        # Trinity输出文件可能的路径 | Possible Trinity output file paths
        possible_trinity_files = [
            # 标准路径：在指定的输出目录内 | Standard path: inside specified output directory  
            output_dir / "Trinity.fasta",
            # 备用路径1：在父目录中，以目录名为前缀 | Alternative path 1: in parent directory with directory name prefix
            output_dir.parent / f"{output_dir.name}.Trinity.fasta",
            # 备用路径2：在输出目录中，以目录名为前缀 | Alternative path 2: in output directory with directory name prefix
            output_dir / f"{output_dir.name}.Trinity.fasta",
            # 备用路径3：在父目录中 | Alternative path 3: in parent directory
            output_dir.parent / "Trinity.fasta"
        ]
        
        for trinity_fasta in possible_trinity_files:
            if trinity_fasta.exists():
                self.config.trinity_assembly = str(trinity_fasta)
                self.logger.info(f"🎉 Trinity组装完成 | Trinity assembly completed: {trinity_fasta}")
                return True
        
        # 如果都找不到，列出实际的输出目录内容以便调试 | If not found, list actual output directory contents for debugging
        self.logger.error("❌ Trinity组装失败，未找到输出文件 | Trinity assembly failed, output file not found")
        try:
            self.logger.error(f"📁 输出目录内容 | Output directory contents ({output_dir}):")
            if output_dir.exists():
                for item in output_dir.iterdir():
                    self.logger.error(f"  - {item.name}")
            
            self.logger.error(f"📁 父目录内容 | Parent directory contents ({output_dir.parent}):")
            for item in output_dir.parent.iterdir():
                if item.name.startswith('trinity') or item.name.startswith('Trinity'):
                    self.logger.error(f"  - {item.name}")
        except Exception as e:
            self.logger.error(f"无法列出目录内容: {e}")
        
        return False
