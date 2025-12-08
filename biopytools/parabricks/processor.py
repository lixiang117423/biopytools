# """
# 🔬 parabricks WGS数据处理模块 | parabricks WGS Data Processing Module 🔬
# """

# from pathlib import Path
# from .utils import CommandRunner, FileProcessor
# import shutil

# class parabricksProcessor:
#     """🔬 parabricks数据处理器 | parabricks Data Processor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.file_processor = FileProcessor(config, logger)
    
#     def run_fq2bam(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
#         """
#         步骤1: FASTQ转BAM (比对和排序) 🧬
#         Step 1: FASTQ to BAM (alignment and sorting)
#         """
#         self.logger.info(f"🧬 [fq2bam] 开始处理样品 | Starting fq2bam for sample: {sample_name}")
        
#         output_bam = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
#         # 检查输出文件是否已存在
#         if output_bam.exists():
#             self.logger.info(f"  ⏭️ BAM文件已存在，跳过 | BAM file exists, skipping: {output_bam.name}")
#             return True
        
#         self.logger.info(f"  📄 R1文件 | R1 file: {Path(r1_file).name}")
#         self.logger.info(f"  📄 R2文件 | R2 file: {Path(r2_file).name}")
#         self.logger.info(f"  💾 输出BAM | Output BAM: {output_bam.name}")
        
#         # 构建fq2bam命令
#         fq2bam_cmd = [
#             "pbrun", "fq2bam",
#             "--ref", str(self.config.reference),
#             "--in-fq", str(r1_file), str(r2_file),
#             "--out-bam", str(output_bam),
#             "--read-group-sm", sample_name,  # 🔑 关键：设置样品名
#             "--read-group-pl", "ILLUMINA",   # 测序平台
#         ]
        
#         success = self.cmd_runner.run_container(
#             fq2bam_cmd,
#             f"🧬 fq2bam (样品: {sample_name})"
#         )
        
#         if success:
#             self.logger.info(f"  ✅ fq2bam完成 | fq2bam completed: {sample_name}")
#             if output_bam.exists():
#                 bam_size = self.file_processor.get_file_size(str(output_bam))
#                 self.logger.info(f"  📏 BAM文件大小 | BAM file size: {bam_size}")
#             return True
#         else:
#             self.logger.error(f"  ❌ fq2bam失败 | fq2bam failed: {sample_name}")
#             return False
    
#     def run_haplotypecaller(self, sample_name: str) -> bool:
#         """
#         步骤2: BAM转VCF/GVCF (变异检测) 📜
#         Step 2: BAM to VCF/GVCF (variant calling)
#         """
#         self.logger.info(f"📜 [haplotypecaller] 开始处理样品 | Starting haplotypecaller for sample: {sample_name}")
        
#         input_bam = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
#         # 检查输入BAM文件是否存在
#         if not input_bam.exists():
#             self.logger.error(f"  ❌ BAM文件不存在 | BAM file not found: {input_bam.name}")
#             self.logger.error(f"     请先运行 fq2bam 步骤 | Please run fq2bam step first")
#             return False
        
#         # 根据gvcf参数决定输出文件名
#         if self.config.gvcf:
#             output_vcf = self.config.vcf_output_dir / f"{sample_name}.g.vcf.gz"
#         else:
#             output_vcf = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
        
#         # 检查输出文件是否已存在
#         if output_vcf.exists():
#             self.logger.info(f"  ⏭️ VCF文件已存在，跳过 | VCF file exists, skipping: {output_vcf.name}")
#             return True
        
#         self.logger.info(f"  💾 输入BAM | Input BAM: {input_bam.name}")
#         self.logger.info(f"  📜 输出{'GVCF' if self.config.gvcf else 'VCF'} | Output {'GVCF' if self.config.gvcf else 'VCF'}: {output_vcf.name}")
        
#         # 构建haplotypecaller命令
#         haplotypecaller_cmd = [
#             "pbrun", "haplotypecaller",
#             "--ref", str(self.config.reference),
#             "--in-bam", str(input_bam),
#             "--ploidy", str(self.config.ploidy),
#             "--out-variants", str(output_vcf),
#         ]
        
#         # 如果启用GVCF，添加--gvcf参数
#         if self.config.gvcf:
#             haplotypecaller_cmd.append("--gvcf")
        
#         success = self.cmd_runner.run_container(
#             haplotypecaller_cmd,
#             f"📜 haplotypecaller (样品: {sample_name})"
#         )
        
#         if success:
#             self.logger.info(f"  ✅ haplotypecaller完成 | haplotypecaller completed: {sample_name}")
#             if output_vcf.exists():
#                 vcf_size = self.file_processor.get_file_size(str(output_vcf))
#                 self.logger.info(f"  📏 VCF文件大小 | VCF file size: {vcf_size}")
#             return True
#         else:
#             self.logger.error(f"  ❌ haplotypecaller失败 | haplotypecaller failed: {sample_name}")
#             return False
    
#     def run_genotypegvcf(self, sample_names: list) -> bool:
#         """
#         步骤3: Joint Calling - 使用genotypegvcf进行联合基因分型 🔗
#         Step 3: Joint Calling - Use genotypegvcf for joint genotyping
#         """
#         # 检查是否启用了GVCF
#         if not self.config.gvcf:
#             self.logger.warning("⚠️ genotypegvcf需要GVCF格式 | genotypegvcf requires GVCF format")
#             return False
        
#         self.logger.info("=" * 60)
#         self.logger.info("🔗 [genotypegvcf] 开始Joint Calling | Starting Joint Calling")
        
#         # 收集所有GVCF文件
#         gvcf_files = []
#         missing_samples = []
        
#         for sample_name in sample_names:
#             gvcf_file = self.config.vcf_output_dir / f"{sample_name}.g.vcf.gz"
#             if gvcf_file.exists():
#                 gvcf_files.append(gvcf_file)
#             else:
#                 missing_samples.append(sample_name)
        
#         if missing_samples:
#             self.logger.warning(f"⚠️ 缺少GVCF文件 | Missing GVCF files: {', '.join(missing_samples)}")
        
#         if len(gvcf_files) < 1:
#             self.logger.error(f"❌ 未找到GVCF文件 | No GVCF files found")
#             return False
        
#         self.logger.info(f"📊 找到 {len(gvcf_files)} 个GVCF | Found {len(gvcf_files)} GVCF files")
#         for gvcf in gvcf_files:
#             self.logger.info(f"   - {gvcf.name}")
        
#         # 定义输出文件
#         output_vcf_name = Path(self.config.combined_output_name).with_suffix('').with_suffix('.vcf.gz').name
#         combined_output_vcf = self.config.vcf_output_dir / output_vcf_name
        
#         # 检查输出文件是否已存在
#         if combined_output_vcf.exists():
#             self.logger.info(f"⏭️ 输出VCF已存在 | Output VCF exists: {combined_output_vcf.name}")
#             return True
        
#         self.logger.info("")
#         self.logger.info("🔗 使用genotypegvcf进行Joint Calling | Using genotypegvcf for Joint Calling")
#         self.logger.info(f"📜 输出VCF: {combined_output_vcf.name}")
        
#         # 构建genotypegvcf命令
#         genotypegvcf_cmd = [
#             "pbrun", "genotypegvcf",
#             "--ref", str(self.config.reference),
#             "--out-vcf", str(combined_output_vcf)
#         ]
        
#         # 添加所有GVCF文件
#         for gvcf_file in gvcf_files:
#             genotypegvcf_cmd.extend(["--in-gvcf", str(gvcf_file)])
        
#         self.logger.info(f"🔗 输入 {len(gvcf_files)} 个GVCF文件 | Input {len(gvcf_files)} GVCF files")
        
#         success = self.cmd_runner.run_container(
#             genotypegvcf_cmd,
#             f"🔗 genotypegvcf - Joint Calling ({len(gvcf_files)} samples)"
#         )
        
#         if not success:
#             self.logger.error("❌ genotypegvcf失败 | genotypegvcf failed")
#             return False
        
#         self.logger.info("✅ genotypegvcf完成 | genotypegvcf completed")
        
#         # 报告文件大小和样品信息
#         if combined_output_vcf.exists():
#             vcf_size = self.file_processor.get_file_size(str(combined_output_vcf))
#             self.logger.info(f"📏 输出VCF大小 | Output VCF size: {vcf_size}")
            
#             # 创建索引
#             self.logger.info("  |-- 建立索引 | Creating index")
#             indexing_success = self.cmd_runner.run(
#                 f"bcftools index -t {str(combined_output_vcf)}",
#                 "🧬 索引VCF"
#             )
#             if indexing_success:
#                 self.logger.info("  ✅ 索引完成 | Index created")
#             else:
#                 self.logger.warning("  ⚠️ 索引失败 | Index failed")
            
#             # 显示样品信息
#             self.logger.info("")
#             self.logger.info("📋 最终VCF中的样品 | Samples in final VCF:")
#             list_cmd = f"bcftools query -l {str(combined_output_vcf)}"
#             self.cmd_runner.run(list_cmd, "📊 样品列表")
            
#             return True
#         else:
#             self.logger.error("❌ 输出VCF未生成 | Output VCF not generated")
#             return False
    
#     def process_sample(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
#         """
#         根据workflow配置处理单个样品 🔀
#         Process single sample based on workflow configuration
#         """
#         self.logger.info(f"🚀 开始处理样品 | Processing sample: {sample_name}")
#         self.logger.info(f"🔀 工作流程 | Workflow: {self.config.workflow}")
        
#         success = True
        
#         # 根据workflow执行相应步骤
#         if self.config.workflow in ["fq2bam", "all"]:
#             success = self.run_fq2bam(sample_name, r1_file, r2_file)
#             if not success:
#                 return False
        
#         if self.config.workflow in ["haplotypecaller", "all"]:
#             success = self.run_haplotypecaller(sample_name)
#             if not success:
#                 return False
        
#         # genotypegvcf在所有样本处理完后统一执行，不在这里处理
        
#         if success:
#             self.logger.info(f"✅ 样品 {sample_name} 处理完成 | Sample completed: {sample_name}")
        
#         return success
    
#     def joint_calling(self, sample_names: list) -> bool:
#         """执行Joint Calling的包装方法，保持向后兼容"""
#         return self.run_genotypegvcf(sample_names)


"""
🔬 parabricks WGS数据处理模块 | parabricks WGS Data Processing Module 🔬
"""

from pathlib import Path
from .utils import CommandRunner, FileProcessor

class parabricksProcessor:
    """🔬 parabricks数据处理器 | parabricks Data Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_processor = FileProcessor(config, logger)
    
    def process_sample(self, sample_name: str, r1_file: str, r2_file: str) -> bool:
        """处理单个样品 🧬 | Process single sample"""
        self.logger.info(f"🚀 开始处理样品 | Starting to process sample: {sample_name}")
        
        # 检查输出文件是否已存在 🧐 | Check if output files already exist
        if self.file_processor.check_output_exists(sample_name):
            self.logger.info(f"⏭️ 样品 {sample_name} 已处理完成，跳过 | Sample {sample_name} already processed, skipping")
            return True
        
        # 定义输出文件路径 📂 | Define output file paths
        output_bam = self.config.bam_output_dir / f"{sample_name}.sorted.bam"
        
        # --- MODIFIED: 根据配置决定输出是VCF还是GVCF ---
        if self.config.gvcf:
            output_vcf = self.config.vcf_output_dir / f"{sample_name}.g.vcf.gz"
            vcf_type_str = "GVCF"
        else:
            output_vcf = self.config.vcf_output_dir / f"{sample_name}.vcf.gz"
            vcf_type_str = "VCF"
        
        self.logger.info(f"  📄 R1文件 | R1 file: {Path(r1_file).name}")
        self.logger.info(f"  📄 R2文件 | R2 file: {Path(r2_file).name}")
        self.logger.info(f"  📜 输出 {vcf_type_str} | Output {vcf_type_str}: vcf/{output_vcf.name}")
        self.logger.info(f"  💾 输出 BAM | Output BAM: bam/{output_bam.name}")
        
        # 步骤1: 运行 fq2bam 🧬 | Step 1: Run fq2bam
        self.logger.info(f"  1️⃣  步骤1: FASTQ to BAM | Step 1: FASTQ to BAM")
        fq2bam_cmd = [
            "pbrun", "fq2bam",
            "--ref", str(self.config.reference),
            "--in-fq", str(r1_file), str(r2_file),
            "--out-bam", str(output_bam),
            # --- NEW: 添加Read Group信息，使用正确的样本名 ---
            "--read-group-sm", sample_name,
            "--read-group-pl", "ILLUMINA",  # 平台信息，通常是ILLUMINA
        ]
        
        success = self.cmd_runner.run_container(
            fq2bam_cmd,
            f"🧬 fq2bam (样品: {sample_name}) | fq2bam (sample: {sample_name})"
        )
        
        if not success:
            self.logger.error(f"  ❌ 样品 {sample_name} fq2bam 步骤失败 | Sample {sample_name} fq2bam step failed")
            return False
        
        # 步骤2: 运行 haplotypecaller 📜 | Step 2: Run haplotypecaller
        self.logger.info(f"  2️⃣  步骤2: BAM to VCF/GVCF | Step 2: BAM to VCF/GVCF")
        haplotypecaller_cmd = [
            "pbrun", "haplotypecaller",
            "--ref", str(self.config.reference),
            "--in-bam", str(output_bam),
            "--out-variants", str(output_vcf),
        ]
        
        # --- MODIFIED: 根据配置添加 --gvcf 标志 ---
        if self.config.gvcf:
            haplotypecaller_cmd.append("--gvcf")
        
        success = self.cmd_runner.run_container(
            haplotypecaller_cmd,
            f"📜 haplotypecaller (样品: {sample_name}) | haplotypecaller (sample: {sample_name})"
        )
        
        if success:
            self.logger.info(f"  ✅ ✓ 样品 {sample_name} 处理完成 | Sample {sample_name} processing completed")
            
            # 检查并报告输出文件大小 📊 | Check and report output file sizes
            if output_vcf.exists():
                vcf_size = self.file_processor.get_file_size(str(output_vcf))
                self.logger.info(f"  📏 {vcf_type_str} 文件大小 | {vcf_type_str} file size: {vcf_size}")
            
            if output_bam.exists():
                bam_size = self.file_processor.get_file_size(str(output_bam))
                self.logger.info(f"  📏 BAM文件大小 | BAM file size: {bam_size}")
            
            return True
        else:
            self.logger.error(f"  ❌ ✗ 样品 {sample_name} haplotypecaller 步骤失败 | Sample {sample_name} haplotypecaller step failed")
            return False

    def joint_calling(self, sample_names: list) -> bool:
        """执行Joint Calling，使用 genotypegvcf 合并所有GVCF文件并生成最终VCF 🧬 | Perform Joint Calling using genotypegvcf to combine all GVCFs and generate a final VCF"""
        
        if not self.config.gvcf:
            self.logger.warning("⚠️ Joint calling需要GVCF格式，但当前未启用GVCF输出 | Joint calling requires GVCF format, but GVCF output is not enabled")
            return False
        
        if not self.config.joint_calling:
            self.logger.info("ℹ️ Joint calling未启用，跳过 | Joint calling not enabled, skipping")
            return True
        
        if len(sample_names) < 2:
            self.logger.warning(f"⚠️ 成功处理的样本少于2个 ({len(sample_names)}个)，无法进行joint calling | Less than 2 successfully processed samples ({len(sample_names)}), cannot perform joint calling")
            return False

        self.logger.info("=" * 60)
        self.logger.info("🔗 开始Joint Calling (使用 genotypegvcf) | Starting Joint Calling (using genotypegvcf)")
        
        gvcf_files = []
        missing_samples = []
        
        for sample_name in sample_names:
            gvcf_file = self.config.vcf_output_dir / f"{sample_name}.g.vcf.gz"
            if gvcf_file.exists():
                gvcf_files.append(gvcf_file)
            else:
                missing_samples.append(sample_name)
        
        if missing_samples:
            self.logger.warning(f"⚠️ 以下样本的GVCF文件不存在，将被跳过 | Missing GVCF files for samples: {', '.join(missing_samples)}")
        
        if len(gvcf_files) < 2:
            self.logger.warning(f"⚠️ 找到的GVCF文件少于2个 ({len(gvcf_files)}个)，无法进行joint calling | Found less than 2 GVCF files ({len(gvcf_files)}), cannot perform joint calling")
            return False
        
        self.logger.info(f"📊 找到 {len(gvcf_files)} 个GVCF文件用于合并 | Found {len(gvcf_files)} GVCF files for combining")
        
        output_filename = Path(self.config.combined_output_name).with_suffix('').with_suffix('.vcf.gz').name
        combined_output_vcf = self.config.vcf_output_dir / output_filename
        
        if combined_output_vcf.exists():
            self.logger.info(f"⏭️ Combined VCF文件已存在，跳过 | Combined VCF file already exists: {combined_output_vcf.name}")
            return True
        
        genotype_cmd = [
            "pbrun", "genotypegvcf",
            "--ref", str(self.config.reference),
            "--out-vcf", str(combined_output_vcf)
        ]
        
        for gvcf_file in gvcf_files:
            genotype_cmd.extend(["--in-gvcf", str(gvcf_file)])
        
        self.logger.info(f"📜 输出Combined VCF | Output Combined VCF: {combined_output_vcf.name}")
        self.logger.info(f"🔗 合并并基因分型 {len(gvcf_files)} 个GVCF文件 | Combining and genotyping {len(gvcf_files)} GVCF files")
        
        success = self.cmd_runner.run_container(
            genotype_cmd,
            f"🔗 Joint Calling - genotypegvcf ({len(gvcf_files)} samples)"
        )
        
        if success:
            self.logger.info("✅ Joint Calling (genotypegvcf) 完成 | Joint Calling (genotypegvcf) completed")
            
            if combined_output_vcf.exists():
                vcf_size = self.file_processor.get_file_size(str(combined_output_vcf))
                self.logger.info(f"📏 Combined VCF文件大小 | Combined VCF file size: {vcf_size}")
                
                self.logger.info("  |-- 建立VCF.gz索引 | Indexing VCF.gz file")
                indexing_success = self.cmd_runner.run(
                    f"bcftools index -t {str(combined_output_vcf)}",
                    "🧬 Indexing final VCF.gz"
                )
                if indexing_success:
                    self.logger.info("  ✅ ✓ VCF.gz索引创建成功 | VCF.gz index created successfully")
                else:
                    self.logger.warning("  ⚠️ VCF.gz索引创建失败 | Failed to create VCF.gz index")
            return True
        else:
            self.logger.error("❌ Joint Calling (genotypegvcf) 失败 | Joint Calling (genotypegvcf) failed")
            return False