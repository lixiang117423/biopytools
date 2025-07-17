"""
Augustus多转录组数据处理模块 | Augustus Multiple RNA-seq Data Processing Module
"""

import os
from pathlib import Path
from .utils import CommandRunner, FileValidator

class HISAT2Manager:
    """HISAT2索引和比对管理器 | HISAT2 Index and Alignment Manager"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)
    
    def build_hisat2_index(self) -> bool:
        """构建HISAT2索引 | Build HISAT2 index"""
        import os
        
        index_file = f"{self.config.hisat2_index}.1.ht2"
        
        if self.file_validator.check_file_exists(index_file, "HISAT2索引 | HISAT2 index"):
            return True
        
        self.logger.info("构建HISAT2索引 | Building HISAT2 index...")
        cmd = f"hisat2-build -p {self.config.threads} {self.config.genome_file} {self.config.hisat2_index}"
        
        success, output = self.cmd_runner.run(cmd, "构建HISAT2索引 | Building HISAT2 index")
        if success:
            self.logger.info("HISAT2索引构建完成 | HISAT2 index building completed")
        
        return success
    
    def align_sample(self, sample_info: dict) -> tuple:
        """对单个样本进行比对 | Align single sample"""
        import os
        
        sample_name = sample_info['name']
        r1_file = sample_info['r1_file']
        r2_file = sample_info['r2_file']
        
        sample_dir = self.config.output_path / sample_name
        sample_dir.mkdir(exist_ok=True)
        
        sam_file = sample_dir / f"{sample_name}.sam"
        bam_file = sample_dir / f"{sample_name}_sorted.bam"
        
        self.logger.info(f"处理样本 | Processing sample: {sample_name}")
        
        # HISAT2比对 | HISAT2 alignment
        if not bam_file.exists():
            self.logger.info(f"  - HISAT2比对 | HISAT2 alignment...")
            hisat2_cmd = (
                f"hisat2 -p {self.config.threads} "
                f"-x {self.config.hisat2_index} "
                f"-1 {r1_file} -2 {r2_file} "
                f"-S {sam_file} "
                f"2> {sample_dir}/{sample_name}_hisat2.log"
            )
            
            success, _ = self.cmd_runner.run(hisat2_cmd, f"HISAT2比对 | HISAT2 alignment - {sample_name}")
            if not success:
                return False, None
            
            # SAM转BAM并排序 | Convert SAM to BAM and sort
            self.logger.info(f"  - 转换和排序BAM | Convert and sort BAM...")
            sort_cmd = (
                f"samtools view -@ {self.config.threads} -bS {sam_file} | "
                f"samtools sort -@ {self.config.threads} -o {bam_file}"
            )
            
            success, _ = self.cmd_runner.run(sort_cmd, f"排序BAM | Sort BAM - {sample_name}")
            if not success:
                return False, None
            
            # 建立索引 | Build index
            index_cmd = f"samtools index -@ {self.config.threads} {bam_file}"
            self.cmd_runner.run(index_cmd, f"建立BAM索引 | Build BAM index - {sample_name}")
            
            # 删除SAM文件节省空间 | Remove SAM file to save space
            if sam_file.exists():
                sam_file.unlink()
        
        return True, str(bam_file)

class BAMProcessor:
    """BAM文件处理器 | BAM File Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner, has_augustus_tools: bool = True):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.has_augustus_tools = has_augustus_tools
    
    # def filter_bam(self, bam_file: str) -> str:
    #     """过滤BAM文件 | Filter BAM file"""
    #     if not self.config.filter_bam or not self.has_augustus_tools:
    #         if not self.has_augustus_tools:
    #             self.logger.info(f"  - 跳过BAM过滤（缺少filterBam工具）| Skipping BAM filtering (missing filterBam tool): {bam_file}")
    #         else:
    #             self.logger.info(f"  - 跳过BAM过滤（用户指定）| Skipping BAM filtering (user specified): {bam_file}")
    #         return bam_file
        
    #     sample_dir = Path(bam_file).parent
    #     sample_name = Path(bam_file).stem.replace('_sorted', '')
    #     filtered_bam = sample_dir / f"{sample_name}_filtered.bam"
        
    #     if filtered_bam.exists():
    #         self.logger.info(f"  - 过滤后的BAM文件已存在 | Filtered BAM file already exists: {filtered_bam}")
    #         return str(filtered_bam)
        
    #     self.logger.info(f"  - 过滤BAM文件 | Filtering BAM file...")
    #     filter_cmd = (
    #         f"filterBam --uniq --paired --pairwiseAlignments "
    #         f"--in {bam_file} --out {filtered_bam}"
    #     )
        
    #     success, _ = self.cmd_runner.run(filter_cmd, "过滤BAM文件 | Filter BAM file")
    #     if success:
    #         # 建立过滤后BAM的索引 | Build index for filtered BAM
    #         index_cmd = f"samtools index -@ {self.config.threads} {filtered_bam}"
    #         self.cmd_runner.run(index_cmd, "建立过滤后BAM索引 | Build filtered BAM index")
            
    #         return str(filtered_bam)
        
    #     return bam_file

    def filter_bam(self, bam_file: str) -> str:
        """过滤BAM文件 | Filter BAM file"""
        if not self.config.filter_bam or not self.has_augustus_tools:
            if not self.has_augustus_tools:
                self.logger.info(f"  - 跳过BAM过滤（缺少filterBam工具）| Skipping BAM filtering (missing filterBam tool): {bam_file}")
            else:
                self.logger.info(f"  - 跳过BAM过滤（用户指定）| Skipping BAM filtering (user specified): {bam_file}")
            return bam_file
        
        sample_dir = Path(bam_file).parent
        sample_name = Path(bam_file).stem.replace('_sorted', '')
        
        # 按查询名称排序的BAM文件 | BAM file sorted by query name
        queryname_sorted_bam = sample_dir / f"{sample_name}_queryname_sorted.bam"
        filtered_bam = sample_dir / f"{sample_name}_filtered.bam"
        
        if filtered_bam.exists():
            self.logger.info(f"  - 过滤后的BAM文件已存在 | Filtered BAM file already exists: {filtered_bam}")
            return str(filtered_bam)
        
        # 步骤1: 按查询名称重新排序BAM文件 | Step 1: Resort BAM file by query name
        if not queryname_sorted_bam.exists():
            self.logger.info(f"  - 按查询名称排序BAM文件 | Sorting BAM file by query name...")
            sort_cmd = f"samtools sort -n -@ {self.config.threads} -o {queryname_sorted_bam} {bam_file}"
            
            success, _ = self.cmd_runner.run(sort_cmd, "按查询名称排序BAM | Sort BAM by query name")
            if not success:
                return bam_file
        
        # 步骤2: 过滤BAM文件 | Step 2: Filter BAM file
        self.logger.info(f"  - 过滤BAM文件 | Filtering BAM file...")
        filter_cmd = (
            f"filterBam --uniq --paired --pairwiseAlignments "
            f"--in {queryname_sorted_bam} --out {filtered_bam}"
        )
        
        success, _ = self.cmd_runner.run(filter_cmd, "过滤BAM文件 | Filter BAM file")
        if success:
            # 步骤3: 将过滤后的BAM文件重新按坐标排序 | Step 3: Resort filtered BAM by coordinates
            final_filtered_bam = sample_dir / f"{sample_name}_filtered_sorted.bam"
            resort_cmd = f"samtools sort -@ {self.config.threads} -o {final_filtered_bam} {filtered_bam}"
            
            resort_success, _ = self.cmd_runner.run(resort_cmd, "重新排序过滤后的BAM | Resort filtered BAM")
            if resort_success:
                # 建立索引 | Build index
                index_cmd = f"samtools index  -@ {self.config.threads} {final_filtered_bam}"
                self.cmd_runner.run(index_cmd, "建立过滤后BAM索引 | Build filtered BAM index")
                
                # 删除中间文件 | Remove intermediate files
                if queryname_sorted_bam.exists():
                    queryname_sorted_bam.unlink()
                if filtered_bam.exists():
                    filtered_bam.unlink()
                
                return str(final_filtered_bam)
            else:
                # 如果重新排序失败，返回原来的过滤文件 | If resort fails, return original filtered file
                index_cmd = f"samtools index  -@ {self.config.threads} {filtered_bam}"
                self.cmd_runner.run(index_cmd, "建立过滤后BAM索引 | Build filtered BAM index")
                return str(filtered_bam)
        
        return bam_file
    
    def get_alignment_stats(self, bam_file: str) -> dict:
        """获取比对统计信息 | Get alignment statistics"""
        sample_dir = Path(bam_file).parent
        sample_name = Path(bam_file).stem.replace('_sorted', '').replace('_filtered', '')
        stats_file = sample_dir / f"{sample_name}_stats.txt"
        
        stats_cmd = f"samtools flagstat {bam_file}"
        success, output = self.cmd_runner.run(stats_cmd, "获取比对统计 | Get alignment statistics")
        
        if success:
            with open(stats_file, 'w') as f:
                f.write(output)
            
            self.logger.info(f"  - 比对统计 | Alignment statistics:")
            for line in output.split('\n'):
                if line.strip():
                    self.logger.info(f"    {line}")
            
            return {'stats_file': str(stats_file), 'stats_content': output}
        
        return {}

class HintsGenerator:
    """Hints文件生成器 | Hints File Generator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner, has_augustus_tools: bool = True):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.has_augustus_tools = has_augustus_tools
    
    def generate_hints_for_sample(self, bam_file: str) -> str:
        """为单个样本生成hints | Generate hints for single sample"""
        sample_dir = Path(bam_file).parent
        sample_name = Path(bam_file).stem.replace('_sorted', '').replace('_filtered', '')
        hints_file = sample_dir / f"{sample_name}_introns.gff"
        
        if hints_file.exists():
            self.logger.info(f"  - Hints文件已存在 | Hints file already exists: {hints_file}")
            return str(hints_file)
        
        if not self.has_augustus_tools:
            self.logger.warning(f"  - 跳过hints生成（缺少bam2hints工具）| Skipping hints generation (missing bam2hints tool)")
            # 创建空的hints文件 | Create empty hints file
            hints_file.touch()
            return str(hints_file)
        
        self.logger.info(f"  - 生成hints文件 | Generating hints file...")
        hints_cmd = f"bam2hints --intronsonly --in={bam_file} --out={hints_file}"
        
        success, _ = self.cmd_runner.run(hints_cmd, f"生成hints | Generate hints - {sample_name}")
        if success:
            return str(hints_file)
        
        # 如果生成失败，创建空文件 | If generation fails, create empty file
        self.logger.warning(f"  - hints生成失败，创建空文件 | Hints generation failed, creating empty file")
        hints_file.touch()
        return str(hints_file)
    
    def merge_and_filter_hints(self, hints_files: list) -> str:
        """合并和过滤hints文件 | Merge and filter hints files"""
        import os
        
        self.logger.info("合并hints文件 | Merging hints files...")
        
        combined_hints = self.config.output_path / "combined_hints.gff"
        sorted_hints = self.config.output_path / "sorted_hints.gff"
        filtered_hints = self.config.output_path / "filtered_hints.gff"
        
        # 检查是否有有效的hints文件 | Check if there are valid hints files
        valid_hints_files = []
        for hints_file in hints_files:
            if os.path.exists(hints_file) and os.path.getsize(hints_file) > 0:
                valid_hints_files.append(hints_file)
        
        if not valid_hints_files:
            self.logger.warning("没有有效的hints文件，创建空的hints文件 | No valid hints files, creating empty hints file")
            filtered_hints.touch()
            return str(filtered_hints)
        
        # 合并所有hints文件 | Merge all hints files
        with open(combined_hints, 'w') as outfile:
            for hints_file in valid_hints_files:
                with open(hints_file, 'r') as infile:
                    outfile.write(infile.read())
        
        # 如果合并后的文件为空，直接创建空的过滤文件 | If merged file is empty, create empty filtered file
        if os.path.getsize(combined_hints) == 0:
            self.logger.warning("合并后的hints文件为空 | Merged hints file is empty")
            filtered_hints.touch()
            return str(filtered_hints)
        
        # 排序hints | Sort hints
        sort_cmd = f"sort -k1,1 -k4,4n {combined_hints} > {sorted_hints}"
        self.cmd_runner.run(sort_cmd, "排序hints | Sort hints")
        
        # 过滤低支持度的hints | Filter low-support hints
        filter_cmd = (
            f"awk '$2==\"b2h\" && $3==\"intron\" && $6>={self.config.min_intron_support}' "
            f"{sorted_hints} > {filtered_hints}"
        )
        self.cmd_runner.run(filter_cmd, "过滤hints | Filter hints")
        
        # 统计hints | Count hints
        self._count_hints(combined_hints, filtered_hints)
        
        return str(filtered_hints)
    
    def _count_hints(self, combined_hints: Path, filtered_hints: Path):
        """统计hints信息 | Count hints information"""
        import os
        
        try:
            with open(combined_hints, 'r') as f:
                total_hints = sum(1 for line in f)
            
            with open(filtered_hints, 'r') as f:
                filtered_count = sum(1 for line in f)
            
            with open(filtered_hints, 'r') as f:
                intron_count = sum(1 for line in f if 'intron' in line)
            
            self.logger.info("Hints统计 | Hints statistics:")
            self.logger.info(f"  总hints数 | Total hints: {total_hints}")
            self.logger.info(f"  过滤后hints数 | Filtered hints: {filtered_count}")
            self.logger.info(f"  内含子hints | Intron hints: {intron_count}")
            
        except Exception as e:
            self.logger.warning(f"统计hints失败 | Failed to count hints: {e}")
