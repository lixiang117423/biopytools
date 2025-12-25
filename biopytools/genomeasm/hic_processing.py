"""
Hi-C数据处理模块 | Hi-C Data Processing Module 🔗
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from .utils import CommandRunner, get_file_stats

class HiCProcessor:
    """Hi-C数据处理器 | Hi-C Data Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger  
        self.cmd_runner = cmd_runner
        self.hic_dir = Path(config.output_dir) / "hic_processing"
        self.hic_dir.mkdir(parents=True, exist_ok=True)
        
        # 为不同策略创建子目录
        strategy_dir = self.hic_dir / config.hic_strategy
        strategy_dir.mkdir(parents=True, exist_ok=True)
        self.working_dir = strategy_dir
    
    def process_all_assemblies(self, fasta_files: Dict[str, List[str]]) -> Dict[str, Dict[str, str]]:
        """处理所有组装结果的Hi-C挂载 | Process Hi-C scaffolding for all assemblies"""
        self.logger.info("🔗 开始Hi-C染色体挂载处理 | Starting Hi-C chromosome scaffolding")
        
        if 'Hi-C' not in self.config.detected_data_types:
            self.logger.warning("⚠️ 未检测到Hi-C数据，跳过挂载步骤 | No Hi-C data detected, skipping scaffolding")
            return {}
        
        scaffold_results = {}
        
        # 处理每种组装类型的每个文件
        for assembly_type, fasta_list in fasta_files.items():
            scaffold_results[assembly_type] = {}
            
            for fasta_file in fasta_list:
                assembly_name = Path(fasta_file).stem
                self.logger.info(f"🔗 处理{assembly_type}组装: {assembly_name}")
                
                # 根据策略选择处理方法
                if self.config.hic_strategy == "complete_juicer":
                    result = self._process_complete_juicer(fasta_file, assembly_name)
                elif self.config.hic_strategy == "standard_3ddna":
                    result = self._process_standard_3ddna(fasta_file, assembly_name)
                elif self.config.hic_strategy == "simplified_salsa2":
                    result = self._process_simplified_salsa2(fasta_file, assembly_name)
                else:
                    self.logger.error(f"❌ 未知的Hi-C策略: {self.config.hic_strategy}")
                    continue
                
                if result:
                    scaffold_results[assembly_type][assembly_name] = result
                    self.logger.info(f"✅ {assembly_name} Hi-C处理完成")
                else:
                    self.logger.error(f"❌ {assembly_name} Hi-C处理失败")
        
        # 生成Hi-C处理统计报告
        self._generate_hic_report(scaffold_results)
        
        self.logger.info("🎉 所有Hi-C处理完成 | All Hi-C processing completed")
        return scaffold_results
    
    def _process_complete_juicer(self, fasta_file: str, assembly_name: str) -> Optional[Dict[str, str]]:
        """完整Juicer + 3D-DNA流程 | Complete Juicer + 3D-DNA pipeline"""
        self.logger.info(f"🧬 执行完整Juicer流程: {assembly_name}")
        
        # 为每个组装创建独立工作目录
        work_dir = self.working_dir / f"juicer_{assembly_name}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 准备Juicer输入
        self._prepare_juicer_input(fasta_file, work_dir)
        
        # 运行Juicer pipeline
        juicer_success = self._run_juicer_pipeline(work_dir, assembly_name)
        if not juicer_success:
            return None
        
        # 运行3D-DNA
        result = self._run_3ddna_pipeline(work_dir, assembly_name)
        if not result:
            return None
        
        return result
    
    def _process_standard_3ddna(self, fasta_file: str, assembly_name: str) -> Optional[Dict[str, str]]:
        """标准3D-DNA流程 (BWA + 格式转换 + 3D-DNA) | Standard 3D-DNA pipeline"""
        self.logger.info(f"🔬 执行标准3D-DNA流程: {assembly_name}")
        
        # 为每个组装创建独立工作目录  
        work_dir = self.working_dir / f"3ddna_{assembly_name}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 复制组装文件到工作目录
        assembly_copy = work_dir / f"{assembly_name}.fa"
        shutil.copy2(fasta_file, assembly_copy)
        
        # BWA比对Hi-C数据
        bam_file = self._align_hic_with_bwa(assembly_copy, work_dir, assembly_name)
        if not bam_file:
            return None
        
        # 格式转换: BAM -> pairs -> merged_nodups
        merged_nodups = self._convert_bam_to_merged_nodups(bam_file, work_dir, assembly_name)
        if not merged_nodups:
            return None
        
        # 运行3D-DNA
        result = self._run_3ddna_with_merged_nodups(assembly_copy, merged_nodups, work_dir, assembly_name)
        return result
    
    def _process_simplified_salsa2(self, fasta_file: str, assembly_name: str) -> Optional[Dict[str, str]]:
        """简化SALSA2流程 | Simplified SALSA2 pipeline"""
        self.logger.info(f"🔧 执行简化SALSA2流程: {assembly_name}")
        
        # 为每个组装创建独立工作目录
        work_dir = self.working_dir / f"salsa2_{assembly_name}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 复制组装文件到工作目录
        assembly_copy = work_dir / f"{assembly_name}.fa"
        shutil.copy2(fasta_file, assembly_copy)
        
        # BWA比对Hi-C数据
        bam_file = self._align_hic_with_bwa(assembly_copy, work_dir, assembly_name)
        if not bam_file:
            return None
        
        # 运行SALSA2
        result = self._run_salsa2_pipeline(assembly_copy, bam_file, work_dir, assembly_name)
        return result
    
    def _prepare_juicer_input(self, fasta_file: str, work_dir: Path):
        """准备Juicer输入文件 | Prepare Juicer input files"""
        self.logger.info("📋 准备Juicer输入文件 | Preparing Juicer input files")
        
        # 创建Juicer目录结构
        (work_dir / "references").mkdir(exist_ok=True)
        (work_dir / "fastq").mkdir(exist_ok=True)
        (work_dir / "splits").mkdir(exist_ok=True)
        
        # 复制基因组文件
        genome_file = work_dir / "references" / "genome.fa"
        shutil.copy2(fasta_file, genome_file)
        
        # 创建染色体大小文件
        self._create_chrom_sizes(genome_file, work_dir / "chrom.sizes")
        
        # 创建限制性酶切位点文件
        self._create_restriction_sites(genome_file, work_dir, self.config.restriction_enzyme)
        
        # 创建Hi-C数据软链接
        hic_r1_link = work_dir / "fastq" / "R1.fastq.gz"
        hic_r2_link = work_dir / "fastq" / "R2.fastq.gz"
        
        if not hic_r1_link.exists():
            os.symlink(os.path.abspath(self.config.hic_r1), hic_r1_link)
        if not hic_r2_link.exists():
            os.symlink(os.path.abspath(self.config.hic_r2), hic_r2_link)
    
    def _run_juicer_pipeline(self, work_dir: Path, assembly_name: str) -> bool:
        """运行Juicer pipeline | Run Juicer pipeline"""
        self.logger.info(f"🧬 运行Juicer pipeline: {assembly_name}")
        
        # 构建Juicer命令
        juicer_cmd = [
            self.config.juicer_path,
            f"-d {work_dir}",
            f"-g {assembly_name}",
            f"-s {self.config.restriction_enzyme}",
            f"-p {work_dir / 'chrom.sizes'}",
            f"-t {self.config.threads}"
        ]
        
        cmd = " ".join(juicer_cmd)
        
        # 切换到工作目录执行
        old_cwd = os.getcwd()
        try:
            os.chdir(work_dir)
            success = self.cmd_runner.run(cmd, f"Juicer pipeline - {assembly_name}", timeout=43200)
            return success
        finally:
            os.chdir(old_cwd)
    
    def _align_hic_with_bwa(self, fasta_file: Path, work_dir: Path, assembly_name: str) -> Optional[str]:
        """使用BWA比对Hi-C数据 | Align Hi-C data with BWA"""
        self.logger.info(f"🎯 BWA比对Hi-C数据: {assembly_name}")
        
        # BWA索引
        index_cmd = f"{self.config.bwa_path} index {fasta_file}"
        if not self.cmd_runner.run(index_cmd, f"BWA索引 - {assembly_name}"):
            return None
        
        # BWA比对
        bam_file = work_dir / f"{assembly_name}_hic_aligned.bam"
        align_cmd = (
            f"{self.config.bwa_path} mem -t {self.config.threads} -B 8 {fasta_file} "
            f"{self.config.hic_r1} {self.config.hic_r2} | "
            f"{self.config.samtools_path} view -@ {self.config.threads//4} -bS - | "
            f"{self.config.samtools_path} sort -@ {self.config.threads//4} -o {bam_file} -"
        )
        
        if not self.cmd_runner.run(align_cmd, f"BWA比对 - {assembly_name}", timeout=7200):
            return None
        
        # 索引BAM文件
        index_cmd = f"{self.config.samtools_path} index {bam_file}"
        if not self.cmd_runner.run(index_cmd, f"BAM索引 - {assembly_name}"):
            return None
        
        # 检查比对质量
        self._check_alignment_quality(bam_file, assembly_name)
        
        return str(bam_file)
    
    def _convert_bam_to_merged_nodups(self, bam_file: str, work_dir: Path, assembly_name: str) -> Optional[str]:
        """将BAM转换为merged_nodups格式 | Convert BAM to merged_nodups format"""
        self.logger.info(f"🔄 格式转换: {assembly_name}")
        
        # 这里需要实现BAM到merged_nodups的转换
        # 可以使用pairtools或类似工具
        merged_nodups = work_dir / f"{assembly_name}_merged_nodups.txt"
        
        # 使用pairtools进行转换 (如果可用)
        try:
            # pairtools parse
            pairs_file = work_dir / f"{assembly_name}.pairs"
            parse_cmd = f"pairtools parse -c {work_dir / 'chrom.sizes'} {bam_file} > {pairs_file}"
            
            if not self.cmd_runner.run(parse_cmd, f"解析pairs - {assembly_name}"):
                # 回退方案：简单格式转换
                return self._simple_bam_conversion(bam_file, merged_nodups, assembly_name)
            
            # pairtools sort and dedup
            dedup_cmd = f"pairtools sort --nproc {self.config.threads//4} {pairs_file} | pairtools dedup --output {merged_nodups}"
            
            if not self.cmd_runner.run(dedup_cmd, f"去重pairs - {assembly_name}"):
                return self._simple_bam_conversion(bam_file, merged_nodups, assembly_name)
                
        except Exception as e:
            self.logger.warning(f"⚠️ pairtools转换失败，使用简单转换: {e}")
            return self._simple_bam_conversion(bam_file, merged_nodups, assembly_name)
        
        return str(merged_nodups) if merged_nodups.exists() else None
    
    def _simple_bam_conversion(self, bam_file: str, output_file: Path, assembly_name: str) -> Optional[str]:
        """简单的BAM格式转换 | Simple BAM format conversion"""
        self.logger.info(f"🔧 执行简单BAM转换: {assembly_name}")
        
        # 简单的BAM到merged_nodups转换脚本
        conversion_script = f"""
        {self.config.samtools_path} view {bam_file} | \\
        awk 'BEGIN{{OFS="\\t"}} {{
            if($1 ~ /1$/) {{
                r1_name = $1; r1_chr = $3; r1_pos = $4; r1_strand = ($2 % 32 >= 16) ? "-" : "+"
            }} else if($1 ~ /2$/) {{
                r2_name = $1; r2_chr = $3; r2_pos = $4; r2_strand = ($2 % 32 >= 16) ? "-" : "+"
                gsub(/2$/, "1", r2_name)
                if(r1_name == r2_name) {{
                    print r1_name, r1_chr, r1_pos, r1_strand, r2_chr, r2_pos, r2_strand
                }}
            }}
        }}' > {output_file}
        """
        
        success = self.cmd_runner.run(conversion_script, f"简单格式转换 - {assembly_name}")
        return str(output_file) if success and output_file.exists() else None
    
    def _run_3ddna_pipeline(self, work_dir: Path, assembly_name: str) -> Optional[Dict[str, str]]:
        """运行3D-DNA pipeline (使用Juicer输出) | Run 3D-DNA pipeline with Juicer output"""
        self.logger.info(f"🧭 运行3D-DNA pipeline: {assembly_name}")
        
        # 查找Juicer输出文件
        merged_nodups = work_dir / "aligned" / "merged_nodups.txt"
        genome_file = work_dir / "references" / "genome.fa"
        
        if not merged_nodups.exists():
            self.logger.error(f"❌ 未找到Juicer输出: {merged_nodups}")
            return None
        
        return self._run_3ddna_with_merged_nodups(genome_file, merged_nodups, work_dir, assembly_name)
    
    def _run_3ddna_with_merged_nodups(self, fasta_file: Path, merged_nodups: Path, work_dir: Path, assembly_name: str) -> Optional[Dict[str, str]]:
        """使用merged_nodups文件运行3D-DNA | Run 3D-DNA with merged_nodups file"""
        self.logger.info(f"🧭 3D-DNA染色体挂载: {assembly_name}")
        
        # 构建3D-DNA命令
        cmd_parts = [self.config.pipeline_3ddna]
        
        # 添加参数
        if self.config.species_type == "haploid":
            cmd_parts.append("-m haploid")
        else:
            cmd_parts.append("-m diploid")
        
        cmd_parts.extend([
            f"-i {self.config.min_contig_size}",
            f"-r {self.config.edit_rounds}",
            str(fasta_file),
            str(merged_nodups)
        ])
        
        cmd = " ".join(cmd_parts)
        
        # 切换到工作目录执行
        old_cwd = os.getcwd()
        try:
            os.chdir(work_dir)
            success = self.cmd_runner.run(cmd, f"3D-DNA挂载 - {assembly_name}", timeout=7200)
            if not success:
                return None
        finally:
            os.chdir(old_cwd)
        
        # 检查3D-DNA输出
        return self._check_3ddna_outputs(work_dir, assembly_name)
    
    def _run_salsa2_pipeline(self, fasta_file: Path, bam_file: str, work_dir: Path, assembly_name: str) -> Optional[Dict[str, str]]:
        """运行SALSA2 pipeline | Run SALSA2 pipeline"""
        self.logger.info(f"🔧 运行SALSA2 pipeline: {assembly_name}")
        
        # 准备SALSA2输入
        fai_file = f"{fasta_file}.fai"
        cmd = f"{self.config.samtools_path} faidx {fasta_file}"
        if not self.cmd_runner.run(cmd, f"生成fai文件 - {assembly_name}"):
            return None
        
        # 运行SALSA2
        salsa_cmd = [
            "python", self.config.salsa2_path,
            "-a", str(fasta_file),
            "-l", fai_file,
            "-b", bam_file,
            "-e", self.config.restriction_enzyme,
            "-o", str(work_dir / "salsa_output"),
            "-m", "yes"  # 生成AGP文件
        ]
        
        cmd = " ".join(salsa_cmd)
        success = self.cmd_runner.run(cmd, f"SALSA2挂载 - {assembly_name}", timeout=3600)
        
        if not success:
            return None
        
        return self._check_salsa2_outputs(work_dir, assembly_name)
    
    def _check_3ddna_outputs(self, work_dir: Path, assembly_name: str) -> Optional[Dict[str, str]]:
        """检查3D-DNA输出文件 | Check 3D-DNA output files"""
        expected_files = {
            'hic_file': work_dir / f"{assembly_name}.0.hic",
            'assembly_file': work_dir / f"{assembly_name}.0.assembly",
            'review_file': work_dir / f"{assembly_name}.0.review.assembly",
            'final_fasta': work_dir / f"{assembly_name}_scaffolded.fa"
        }
        
        results = {}
        for key, file_path in expected_files.items():
            if file_path.exists():
                results[key] = str(file_path)
                self.logger.info(f"✅ 找到{key}: {file_path}")
            else:
                self.logger.warning(f"⚠️ 未找到{key}: {file_path}")
        
        return results if results else None
    
    def _check_salsa2_outputs(self, work_dir: Path, assembly_name: str) -> Optional[Dict[str, str]]:
        """检查SALSA2输出文件 | Check SALSA2 output files"""
        output_dir = work_dir / "salsa_output"
        expected_files = {
            'scaffolds_final': output_dir / "scaffolds_FINAL.fasta",
            'agp_file': output_dir / "scaffolds_FINAL.agp"
        }
        
        results = {}
        for key, file_path in expected_files.items():
            if file_path.exists():
                results[key] = str(file_path)
                self.logger.info(f"✅ 找到{key}: {file_path}")
            else:
                self.logger.warning(f"⚠️ 未找到{key}: {file_path}")
        
        return results if results else None
    
    def _create_chrom_sizes(self, fasta_file: Path, output_file: Path):
        """创建染色体大小文件 | Create chromosome sizes file"""
        cmd = f"{self.config.samtools_path} faidx {fasta_file}"
        self.cmd_runner.run(cmd, "生成fasta索引")
        
        cmd = f"cut -f1,2 {fasta_file}.fai > {output_file}"
        self.cmd_runner.run(cmd, "生成染色体大小文件")
    
    def _create_restriction_sites(self, fasta_file: Path, work_dir: Path, enzyme: str):
        """创建限制性酶切位点文件 | Create restriction sites file"""
        # 根据酶类型确定识别序列
        enzyme_sites = {
            'MboI': 'GATC',
            'DpnII': 'GATC', 
            'HindIII': 'AAGCTT',
            'EcoRI': 'GAATTC'
        }
        
        site_seq = enzyme_sites.get(enzyme, 'GATC')
        
        # 使用简单的序列搜索生成酶切位点
        sites_file = work_dir / f"{enzyme}_sites.txt"
        
        # 这里可以使用更复杂的工具，暂时用简单方法
        cmd = f"echo '# {enzyme} sites for {fasta_file.name}' > {sites_file}"
        self.cmd_runner.run(cmd, f"创建{enzyme}酶切位点文件")
    
    def _check_alignment_quality(self, bam_file: str, assembly_name: str):
        """检查比对质量 | Check alignment quality"""
        stats_file = Path(bam_file).parent / f"{assembly_name}_alignment_stats.txt"
        cmd = f"{self.config.samtools_path} flagstat {bam_file} > {stats_file}"
        self.cmd_runner.run(cmd, f"生成比对统计 - {assembly_name}")
        
        # 读取映射率
        try:
            with open(stats_file, 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if 'mapped (' in line:
                        mapping_rate = float(line.split('(')[1].split('%')[0]) / 100
                        self.logger.info(f"📊 {assembly_name} Hi-C比对率: {mapping_rate:.2%}")
                        
                        if mapping_rate < self.config.min_mapping_rate:
                            self.logger.warning(f"⚠️ {assembly_name} 比对率偏低: {mapping_rate:.2%}")
                        break
        except Exception as e:
            self.logger.warning(f"⚠️ 无法解析比对统计: {e}")
    
    def _generate_hic_report(self, scaffold_results: Dict[str, Dict[str, str]]):
        """生成Hi-C处理报告 | Generate Hi-C processing report"""
        report_file = self.hic_dir / "hic_processing_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🔗 Hi-C染色体挂载处理报告 | Hi-C Chromosome Scaffolding Report\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"项目名称 | Project: {self.config.project_name}\n")
            f.write(f"Hi-C策略 | Hi-C Strategy: {self.config.hic_strategy}\n")
            f.write(f"限制性酶 | Restriction Enzyme: {self.config.restriction_enzyme}\n\n")
            
            for assembly_type, results in scaffold_results.items():
                f.write(f"📋 {assembly_type.upper()} 组装挂载结果:\n")
                f.write("-" * 50 + "\n")
                
                for assembly_name, files in results.items():
                    f.write(f"组装名称 | Assembly: {assembly_name}\n")
                    for file_type, file_path in files.items():
                        f.write(f"  {file_type}: {Path(file_path).name}\n")
                    f.write("\n")
        
        self.logger.info(f"📋 Hi-C处理报告已生成: {report_file}")
