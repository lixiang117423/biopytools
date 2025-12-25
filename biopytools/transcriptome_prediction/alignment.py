"""
🧬 HISAT2比对模块 | HISAT2 Alignment Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, SequencingTypeDetector

class HISAT2Aligner:
    """🎯 HISAT2比对器 | HISAT2 Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_index(self):
        """🏗️  构建HISAT2索引 | Build HISAT2 index"""
        index_prefix = self.config.output_path / "hisat2_index" / self.config.base_name
        index_prefix.parent.mkdir(parents=True, exist_ok=True)
        
        self.logger.info("🏗️  构建HISAT2索引 | Building HISAT2 index")
        
        cmd = (
            f"{self.config.hisat2_build_path} "
            f"-p {self.config.threads} "
            f"{self.config.genome_file} "
            f"{index_prefix}"
        )
        
        success = self.cmd_runner.run(cmd, "🏗️  构建HISAT2索引 | Build HISAT2 index")
        
        if success:
            self.config.hisat2_index = str(index_prefix)
            self.logger.info(f"✅ HISAT2索引已生成 | HISAT2 index generated: {index_prefix}")
        
        return success
    
    def _parse_samples_file(self, samples_file: str):
        """📋 解析Trinity样本文件 | Parse Trinity samples file"""
        samples_dict = {}
        rna_seq_files = []
        
        try:
            with open(samples_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    
                    # 检测文件格式
                    if len(parts) >= 4:
                        # Trinity 4列格式: condition  sample_name  left.fq  right.fq
                        condition = parts[0].strip()
                        sample_name = parts[1].strip() 
                        left_file = parts[2].strip()
                        right_file = parts[3].strip() if len(parts) > 3 and parts[3].strip() else None
                        
                        # 检查文件是否存在
                        if not os.path.exists(left_file):
                            self.logger.warning(f"⚠️ 文件不存在 | File not found: {left_file}")
                            continue
                        
                        rna_seq_files.append(left_file)
                        
                        if sample_name not in samples_dict:
                            samples_dict[sample_name] = []
                        samples_dict[sample_name].append(left_file)
                        
                        # 如果有右端文件
                        if right_file and os.path.exists(right_file):
                            rna_seq_files.append(right_file)
                            samples_dict[sample_name].append(right_file)
                        elif right_file and not os.path.exists(right_file):
                            self.logger.warning(f"⚠️ 右端文件不存在 | Right file not found: {right_file}")
                    
                    elif len(parts) >= 3:
                        # Trinity 3列格式 (单端): condition  sample_name  reads.fq
                        condition = parts[0].strip()
                        sample_name = parts[1].strip()
                        reads_file = parts[2].strip()
                        
                        if not os.path.exists(reads_file):
                            self.logger.warning(f"⚠️ 文件不存在 | File not found: {reads_file}")
                            continue
                        
                        rna_seq_files.append(reads_file)
                        
                        if sample_name not in samples_dict:
                            samples_dict[sample_name] = []
                        samples_dict[sample_name].append(reads_file)
                    
                    elif len(parts) >= 2:
                        # 兼容原有2列格式: file_path  sample_name
                        file_path = parts[0].strip()
                        sample_name = parts[1].strip()
                        
                        if not os.path.exists(file_path):
                            self.logger.warning(f"⚠️ 文件不存在 | File not found: {file_path}")
                            continue
                        
                        rna_seq_files.append(file_path)
                        
                        if sample_name not in samples_dict:
                            samples_dict[sample_name] = []
                        samples_dict[sample_name].append(file_path)
                    
                    else:
                        self.logger.warning(f"⚠️ 第{line_num}行格式不正确，跳过 | Line {line_num} format incorrect, skipping: {line}")
                        continue
            
            self.logger.info(f"📋 从样本文件解析出 {len(samples_dict)} 个样本，{len(rna_seq_files)} 个文件 | "
                        f"Parsed {len(samples_dict)} samples and {len(rna_seq_files)} files from samples file")
            
            # 显示解析结果
            for sample_name, files in samples_dict.items():
                self.logger.info(f"  - {sample_name}: {len(files)} 个文件 | files")
            
            return rna_seq_files, samples_dict
        
        except Exception as e:
            self.logger.error(f"❌ 解析样本文件失败 | Failed to parse samples file: {e}")
            return [], {}

    def _organize_samples_from_dict(self, samples_dict):
        """🗂️ 从样本字典组织样本文件 | Organize sample files from samples dictionary"""
        samples = []
        all_paired = True
        all_single = True
        
        for sample_name, files in samples_dict.items():
            if len(files) == 2:
                # 配对末端：尝试识别R1和R2
                r1_file = None
                r2_file = None
                
                for file_path in files:
                    filename = Path(file_path).name.lower()
                    # 更准确的模式匹配
                    if any(pattern in filename for pattern in [
                        '_f1.fq', '_r1.fq', '_1.fq', '.1.fq', 
                        '_f1.fastq', '_r1.fastq', '_1.fastq', '.1.fastq',
                        '_forward.fq', '_left.fq',
                        '_f1.clean.fq', '_1.clean.fq'
                    ]):
                        r1_file = file_path
                    elif any(pattern in filename for pattern in [
                        '_r2.fq', '_f2.fq', '_2.fq', '.2.fq',
                        '_r2.fastq', '_f2.fastq', '_2.fastq', '.2.fastq', 
                        '_reverse.fq', '_right.fq',
                        '_r2.clean.fq', '_2.clean.fq'
                    ]):
                        r2_file = file_path
                
                if r1_file and r2_file:
                    samples.append((r1_file, r2_file))
                    all_single = False
                    self.logger.info(f"  📌 配对样本 | Paired sample: {sample_name}")
                    self.logger.info(f"    - R1: {Path(r1_file).name}")
                    self.logger.info(f"    - R2: {Path(r2_file).name}")
                else:
                    # 无法识别R1/R2，按字母顺序排列
                    sorted_files = sorted(files)
                    samples.append(tuple(sorted_files))
                    all_single = False
                    self.logger.warning(f"⚠️ 无法识别R1/R2，按顺序处理 | Cannot identify R1/R2, processing in order: {sample_name}")
                    
            elif len(files) == 1:
                # 单端
                samples.append((files[0],))
                all_paired = False
                self.logger.info(f"  📌 单端样本 | Single-end sample: {sample_name} - {Path(files[0]).name}")
                
            else:
                self.logger.warning(f"⚠️ 样本 {sample_name} 有 {len(files)} 个文件，跳过 | Sample {sample_name} has {len(files)} files, skipping")
                continue
        
        # 确定整体布局
        if all_paired:
            layout = 'paired'
            self.logger.info("🔍 检测结果：全部为配对末端测序 | Detection result: All paired-end sequencing")
        elif all_single:
            layout = 'single' 
            self.logger.info("🔍 检测结果：全部为单端测序 | Detection result: All single-end sequencing")
        else:
            layout = 'mixed'
            self.logger.info("🔍 检测结果：混合测序类型 | Detection result: Mixed sequencing types")
        
        return layout, samples
    
    # def _parse_samples_file(self, samples_file: str):
    #     """📋 解析Trinity样本文件 | Parse Trinity samples file"""
    #     samples_dict = {}
    #     rna_seq_files = []
        
    #     try:
    #         with open(samples_file, 'r') as f:
    #             for line in f:
    #                 line = line.strip()
    #                 if not line or line.startswith('#'):
    #                     continue
                    
    #                 # Trinity samples file format: file_path<tab>sample_name
    #                 # parts = line.split('\t')

    #                 parts = line.split()
                    
    #                 if len(parts) >= 2:
    #                     file_path = parts[0].strip()
    #                     sample_name = parts[1].strip()
                        
    #                     if not os.path.exists(file_path):
    #                         self.logger.warning(f"⚠️ 文件不存在 | File not found: {file_path}")
    #                         continue
                        
    #                     rna_seq_files.append(file_path)
                        
    #                     if sample_name not in samples_dict:
    #                         samples_dict[sample_name] = []
    #                     samples_dict[sample_name].append(file_path)
            
    #         self.logger.info(f"📋 从样本文件解析出 {len(samples_dict)} 个样本，{len(rna_seq_files)} 个文件 | "
    #                        f"Parsed {len(samples_dict)} samples and {len(rna_seq_files)} files from samples file")
            
    #         return rna_seq_files, samples_dict
        
    #     except Exception as e:
    #         self.logger.error(f"❌ 解析样本文件失败 | Failed to parse samples file: {e}")
    #         return [], {}
    
    def align_reads(self):
        """🎯 比对RNA-seq读段 | Align RNA-seq reads"""
        if not hasattr(self.config, 'hisat2_index'):
            if not self.build_index():
                return False
        
        # 获取RNA-seq文件列表
        if self.config.samples_file:
            # 如果使用样本文件，解析它
            self.logger.info("📋 使用样本文件模式 | Using samples file mode")
            rna_seq_files, samples_dict = self._parse_samples_file(self.config.samples_file)
            if not rna_seq_files:
                self.logger.error("❌ 从样本文件中未解析到有效的RNA-seq文件 | No valid RNA-seq files parsed from samples file")
                return False
        else:
            # 使用直接提供的RNA-seq文件
            rna_seq_files = self.config.rna_seq_files
            samples_dict = {}
        
        # 检测测序类型和样本组织 | Detect sequencing type and sample organization
        if self.config.samples_file:
            # 对于samples file，我们需要按样本组织文件
            layout, samples = self._organize_samples_from_dict(samples_dict)
        else:
            layout, samples = SequencingTypeDetector.detect_sequencing_layout(rna_seq_files)
        
        self.logger.info(f"🔍 检测到测序布局 | Detected sequencing layout: {layout}")
        self.logger.info(f"📊 检测到 {len(samples)} 个样本 | Detected {len(samples)} samples")
        
        if len(samples) == 0:
            self.logger.error("❌ 未检测到有效样本 | No valid samples detected")
            return False
        
        output_dir = self.config.output_path / "hisat2_alignment"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.config.bam_files = []
        
        # 处理每个样本 | Process each sample
        for i, sample_files in enumerate(samples):
            if self.config.samples_file:
                # 从样本字典中获取样本名
                sample_name = self._get_sample_name_from_dict(sample_files, samples_dict)
            else:
                sample_name = f"sample_{i+1:03d}"
            
            if layout == 'paired' or (layout == 'mixed' and len(sample_files) == 2):
                # 配对末端测序 | Paired-end sequencing
                success = self._align_paired_end(sample_files, sample_name, output_dir)
            else:
                # 单端测序 | Single-end sequencing
                success = self._align_single_end(sample_files[0], sample_name, output_dir)
            
            if not success:
                return False
        
        self.logger.info(f"🎉 比对完成，生成 {len(self.config.bam_files)} 个BAM文件 | Alignment completed, generated {len(self.config.bam_files)} BAM files")
        return True
    
    # def _organize_samples_from_dict(self, samples_dict):
    #     """🗂️ 从样本字典组织样本文件 | Organize sample files from samples dictionary"""
    #     samples = []
    #     all_paired = True
    #     all_single = True
        
    #     for sample_name, files in samples_dict.items():
    #         if len(files) == 2:
    #             # 配对末端：尝试识别R1和R2
    #             r1_file = None
    #             r2_file = None
                
    #             for file_path in files:
    #                 filename = Path(file_path).name.lower()
    #                 if any(pattern in filename for pattern in ['_1.clean.fq', '_r1', '_f1', '_1.fq', '.1.fq', '_forward']):
    #                     r1_file = file_path
    #                 elif any(pattern in filename for pattern in ['_2.clean.fq', '_r2', '_r2', '_2.fq', '.2.fq', '_reverse']):
    #                     r2_file = file_path
                
    #             if r1_file and r2_file:
    #                 samples.append((r1_file, r2_file))
    #                 all_single = False
    #             else:
    #                 # 无法识别R1/R2，按顺序排列
    #                 samples.append(tuple(sorted(files)))
    #                 all_single = False
        #     elif len(files) == 1:
        #         # 单端
        #         samples.append((files[0],))
        #         all_paired = False
        #     else:
        #         self.logger.warning(f"⚠️ 样本 {sample_name} 有 {len(files)} 个文件，跳过 | Sample {sample_name} has {len(files)} files, skipping")
        #         continue
        
        # # 确定整体布局
        # if all_paired:
        #     layout = 'paired'
        # elif all_single:
        #     layout = 'single'
        # else:
        #     layout = 'mixed'
        
        # return layout, samples
    
    def _get_sample_name_from_dict(self, sample_files, samples_dict):
        """🏷️ 从样本字典获取样本名 | Get sample name from samples dictionary"""
        for sample_name, files in samples_dict.items():
            if set(sample_files).issubset(set(files)):
                return sample_name
        return "unknown_sample"
    
    # def _align_paired_end(self, sample_files, sample_name, output_dir):
    #     """👥 配对末端比对 | Paired-end alignment"""
    #     r1_file, r2_file = sample_files
    #     sam_file = output_dir / f"{sample_name}.sam"
    #     sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        
    #     self.logger.info(f"👥 执行配对末端读段比对 | Performing paired-end read alignment: {sample_name}")
        
    #     # 构建HISAT2命令 | Build HISAT2 command
    #     cmd = (
    #         f"{self.config.hisat2_path} "
    #         f"-x {self.config.hisat2_index} "
    #         f"-1 {r1_file} "
    #         f"-2 {r2_file} "
    #         f"-S {sam_file} "
    #         f"-p {self.config.threads} "
    #         f"--min-intronlen {self.config.hisat2_min_intron} "
    #         f"--max-intronlen {self.config.hisat2_max_intron}"
    #     )
        
    #     # 添加可选参数 | Add optional parameters
    #     if self.config.hisat2_dta:
    #         cmd += " --dta"
        
    #     if self.config.hisat2_novel_splicesite:
    #         novel_ss_file = output_dir / f"{sample_name}_novel_splicesites.txt"
    #         cmd += f" --novel-splicesite-outfile {novel_ss_file}"
        
    #     if not self.cmd_runner.run(cmd, f"👥 HISAT2配对末端比对 | HISAT2 paired-end alignment: {sample_name}"):
    #         return False
        
    #     # 转换SAM为BAM并排序 | Convert SAM to BAM and sort
    #     if not self._sam_to_sorted_bam(sam_file, sorted_bam):
    #         return False
        
    #     self.config.bam_files.append(str(sorted_bam))
    #     return True
    def _align_paired_end(self, sample_files, sample_name, output_dir):
        """配对末端比对 | Paired-end alignment"""
        r1_file, r2_file = sample_files
        sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        
        self.logger.info(f"执行配对末端读段比对 | Performing paired-end read alignment: {sample_name}")
        
        # HISAT2直接管道到排序BAM | HISAT2 pipeline to sorted BAM
        cmd = (
            f"{self.config.hisat2_path} "
            f"-x {self.config.hisat2_index} "
            f"-1 {r1_file} "
            f"-2 {r2_file} "
            f"-p {self.config.threads} "
            f"--min-intronlen {self.config.hisat2_min_intron} "
            f"--max-intronlen {self.config.hisat2_max_intron}"
        )
        
        # 添加可选参数
        if self.config.hisat2_dta:
            cmd += " --dta"
        
        if self.config.hisat2_novel_splicesite:
            novel_ss_file = output_dir / f"{sample_name}_novel_splicesites.txt"
            cmd += f" --novel-splicesite-outfile {novel_ss_file}"
        
        # 添加管道到samtools
        cmd += (
            f" | {self.config.samtools_path} view -@ {self.config.threads} -bS - "
            f"| {self.config.samtools_path} sort -@ {self.config.threads} - -o {sorted_bam}"
        )
        
        if not self.cmd_runner.run(cmd, f"HISAT2配对末端比对到排序BAM | HISAT2 paired-end alignment to sorted BAM: {sample_name}"):
            return False
        
        # 建立索引
        cmd_index = f"{self.config.samtools_path} index -@ {self.config.threads} {sorted_bam}"
        if not self.cmd_runner.run(cmd_index, f"BAM索引 | BAM indexing: {sorted_bam.name}"):
            return False
        
        self.config.bam_files.append(str(sorted_bam))
        return True
    
    # def _align_single_end(self, read_file, sample_name, output_dir):
    #     """🔗 单端比对 | Single-end alignment"""
    #     sam_file = output_dir / f"{sample_name}.sam"
    #     sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        
    #     self.logger.info(f"🔗 执行单端读段比对 | Performing single-end read alignment: {sample_name}")
        
    #     # 构建HISAT2命令 | Build HISAT2 command
    #     cmd = (
    #         f"{self.config.hisat2_path} "
    #         f"-x {self.config.hisat2_index} "
    #         f"-U {read_file} "
    #         f"-S {sam_file} "
    #         f"-p {self.config.threads} "
    #         f"--min-intronlen {self.config.hisat2_min_intron} "
    #         f"--max-intronlen {self.config.hisat2_max_intron}"
    #     )
        
    #     # 添加可选参数 | Add optional parameters
    #     if self.config.hisat2_dta:
    #         cmd += " --dta"
        
    #     if self.config.hisat2_novel_splicesite:
    #         novel_ss_file = output_dir / f"{sample_name}_novel_splicesites.txt"
    #         cmd += f" --novel-splicesite-outfile {novel_ss_file}"
        
    #     if not self.cmd_runner.run(cmd, f"🔗 HISAT2单端比对 | HISAT2 single-end alignment: {sample_name}"):
    #         return False
        
    #     # 转换SAM为BAM并排序 | Convert SAM to BAM and sort
    #     if not self._sam_to_sorted_bam(sam_file, sorted_bam):
    #         return False
        
    #     self.config.bam_files.append(str(sorted_bam))
    #     return True

    def _align_single_end(self, read_file, sample_name, output_dir):
        """单端比对 | Single-end alignment"""
        sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        
        self.logger.info(f"执行单端读段比对 | Performing single-end read alignment: {sample_name}")
        
        # HISAT2直接管道到排序BAM | HISAT2 pipeline to sorted BAM
        cmd = (
            f"{self.config.hisat2_path} "
            f"-x {self.config.hisat2_index} "
            f"-U {read_file} "
            f"-p {self.config.threads} "
            f"--min-intronlen {self.config.hisat2_min_intron} "
            f"--max-intronlen {self.config.hisat2_max_intron}"
        )
        
        # 添加可选参数
        if self.config.hisat2_dta:
            cmd += " --dta"
        
        if self.config.hisat2_novel_splicesite:
            novel_ss_file = output_dir / f"{sample_name}_novel_splicesites.txt"
            cmd += f" --novel-splicesite-outfile {novel_ss_file}"
        
        # 添加管道到samtools
        cmd += (
            f" | {self.config.samtools_path} view -@ {self.config.threads} -bS - "
            f"| {self.config.samtools_path} sort -@ {self.config.threads} - -o {sorted_bam}"
        )
        
        if not self.cmd_runner.run(cmd, f"HISAT2单端比对到排序BAM | HISAT2 single-end alignment to sorted BAM: {sample_name}"):
            return False
        
        # 建立索引
        cmd_index = f"{self.config.samtools_path} index -@ {self.config.threads} {sorted_bam}"
        if not self.cmd_runner.run(cmd_index, f"BAM索引 | BAM indexing: {sorted_bam.name}"):
            return False
        
        self.config.bam_files.append(str(sorted_bam))
        return True
    
    # def _sam_to_sorted_bam(self, sam_file: Path, sorted_bam: Path) -> bool:
    #     """📄 将SAM文件转换为排序的BAM文件 | Convert SAM file to sorted BAM file"""
    #     # SAM转BAM | SAM to BAM
    #     bam_file = sam_file.with_suffix('.bam')
    #     cmd_sam2bam = f"{self.config.samtools_path} view -@ {self.config.threads} -bS {sam_file} > {bam_file}"
        
    #     if not self.cmd_runner.run(cmd_sam2bam, f"📄 SAM转BAM | SAM to BAM: {sam_file.name}"):
    #         return False
        
    #     # 排序BAM | Sort BAM
    #     cmd_sort = f"{self.config.samtools_path} sort -@ {self.config.threads} {bam_file} -o {sorted_bam}"
        
    #     if not self.cmd_runner.run(cmd_sort, f"📊 BAM排序 | BAM sorting: {bam_file.name}"):
    #         return False
        
    #     # 建立索引 | Create index
    #     cmd_index = f"{self.config.samtools_path} index -@ {self.config.threads} {sorted_bam}"
        
    #     if not self.cmd_runner.run(cmd_index, f"📇 BAM索引 | BAM indexing: {sorted_bam.name}"):
    #         return False
        
    #     # 删除临时文件 | Remove temporary files
    #     try:
    #         sam_file.unlink()
    #         bam_file.unlink()
    #     except:
    #         pass
        
    #     return True