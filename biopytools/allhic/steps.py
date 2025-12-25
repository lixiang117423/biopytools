"""
ALLHiC流水线步骤模块 🚀 | ALLHiC Pipeline Steps Module
"""

import os
from utils import link_file, check_file_exists, clean_work_directory

class ALLHiCSteps:
    """ALLHiC流水线步骤执行器 | ALLHiC Pipeline Steps Executor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def step0_prepare_data(self):
        """Step 0: 准备数据 | Prepare data"""
        if self.config.skip_steps.get("data", False):
            self.logger.log("⏭️ 跳过数据准备 | Skipping data preparation")
            return

        self.logger.log_section("📂 [Step 0] 准备数据 | Preparing input data")
        data_dir = self.config.directories["data"]

        # 清理工作目录 | Clean work directory
        clean_work_directory(data_dir, self.logger)

        os.chdir(data_dir)

        # 链接输入文件 | Link input files
        link_file(self.config.reference, "draft.asm.fasta", self.logger)
        link_file(self.config.read1, "reads_R1.fastq.gz", self.logger)
        link_file(self.config.read2, "reads_R2.fastq.gz", self.logger)
        
        # 构建索引 | Build indices
        if not check_file_exists("draft.asm.fasta.bwt"):
            self.logger.run_command("bwa index draft.asm.fasta", "构建BWA索引 | Building BWA index")
        
        if not check_file_exists("draft.asm.fasta.fai"):
            self.logger.run_command("samtools faidx draft.asm.fasta", "构建samtools索引 | Building samtools index")
        
        self.logger.log_success("数据准备完成 | Data preparation completed")
        self.logger.log_step_time()
    
    def step1_mapping(self):
        """Step 1: 比对 | Mapping"""
        if self.config.skip_steps["mapping"]:
            self.logger.log("⏭️ 跳过比对 | Skipping mapping")
            return
        
        self.logger.log_section("🔍 [Step 1] 比对Hi-C读段 | Mapping Hi-C reads")
        map_dir = self.config.directories["mapping"]
        os.chdir(map_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        
        for ext in ["amb", "ann", "bwt", "pac", "sa", "fai"]:
            link_file(f"{data_dir}/draft.asm.fasta.{ext}", f"draft.asm.fasta.{ext}", self.logger)
        
        link_file(f"{data_dir}/reads_R1.fastq.gz", "reads_R1.fastq.gz", self.logger)
        link_file(f"{data_dir}/reads_R2.fastq.gz", "reads_R2.fastq.gz", self.logger)
        
        clean_bam = "sample.clean.bam"
        if check_file_exists(clean_bam):
            self.logger.log(f"✅ 找到现有文件 | Found existing: {clean_bam}")
            self.logger.log_step_time()
            return
        
        cmd = (f"bwa mem -t {self.config.threads} -5SPM draft.asm.fasta reads_R1.fastq.gz reads_R2.fastq.gz "
               f"2>bwa_mapping.log | "
               f"samtools view -@ {self.config.threads} -hF 2316 -q {self.config.mapq_step1} - | "
               f"samtools sort -@ {self.config.threads} -n -o {clean_bam} -")
        
        self.logger.run_command(cmd, "执行BWA比对和排序 | Executing BWA mapping and sorting")
        self.logger.log_success("比对完成 | Mapping completed")
        self.logger.log_step_time()
    
    def step1_5_allele_detection(self):
        """Step 1.5: 等位基因检测 | Allele detection"""
        if self.config.skip_steps["allele"]:
            self.logger.log("⏭️ 跳过等位基因检测 | Skipping allele detection")
            return
        
        self.logger.log_section("🧬 [Step 1.5] 等位基因检测 | Allele detection")
        allele_dir = self.config.directories["allele"]
        os.chdir(allele_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        map_dir = self.config.directories["mapping"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        link_file(f"{map_dir}/sample.clean.bam", "sample.clean.bam", self.logger)
        
        if check_file_exists("alleles.table"):
            self.logger.log("✅ 找到现有文件 | Found existing: alleles.table")
            return
        
        # 提取counts | Extract counts
        self.logger.run_command(
            f"allhic extract sample.clean.bam draft.asm.fasta --RE {self.config.enzyme}",
            "提取Hi-C counts | Extracting Hi-C counts"
        )
        
        # 构建paf文件 | Build PAF file
        self.logger.run_command(
            f"minimap2 -DP -k19 -w19 -m200 -t {self.config.threads} "
            f"draft.asm.fasta draft.asm.fasta 2>minimap2.log > draft.asm.paf",
            "构建自比对PAF文件 | Building self-alignment PAF file"
        )
        
        # 生成alleles表 | Generate alleles table
        self.logger.run_command(
            f"allhic alleles draft.asm.paf sample.clean.counts_{self.config.enzyme}.txt > alleles.table",
            "生成等位基因表 | Generating alleles table"
        )
        
        self.logger.log_success("等位基因检测完成 | Allele detection completed")
        self.logger.log_step_time()
    
    def step2_pruning(self):
        """Step 2: 修剪 | Pruning"""
        if self.config.skip_steps["prune"]:
            self.logger.log("⏭️ 跳过修剪 | Skipping pruning")
            return
        
        self.logger.log_section("✂️ [Step 2] 修剪嵌合contigs | Pruning chimeric contigs")
        prune_dir = self.config.directories["prune"]
        os.chdir(prune_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        map_dir = self.config.directories["mapping"]
        allele_dir = self.config.directories["allele"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        link_file(f"{map_dir}/sample.clean.bam", "sample.clean.bam", self.logger)
        link_file(f"{allele_dir}/alleles.table", "alleles.table", self.logger)
        
        if check_file_exists("prunning.bam"):
            self.logger.log("✅ 找到现有文件 | Found existing: prunning.bam")
            return
        
        if check_file_exists("alleles.table"):
            self.logger.run_command(
                "ALLHiC_prune -i alleles.table -b sample.clean.bam -r draft.asm.fasta",
                "执行修剪 | Executing pruning"
            )
        else:
            self.logger.log_warning("未找到等位基因表，使用原始BAM | No allele table found, using original BAM")
            link_file("sample.clean.bam", "prunning.bam", self.logger)
        
        self.logger.log_success("修剪完成 | Pruning completed")
        self.logger.log_step_time()
    
    def step3_partition(self):
        """Step 3: 分区 | Partition"""
        if self.config.skip_steps["partition"]:
            self.logger.log("⏭️ 跳过分区 | Skipping partition")
            return
        
        self.logger.log_section("📦 [Step 3] 分区contigs | Partitioning contigs")
        part_dir = self.config.directories["partition"]
        os.chdir(part_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        prune_dir = self.config.directories["prune"]
        map_dir = self.config.directories["mapping"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        
        if check_file_exists(f"{prune_dir}/prunning.bam"):
            link_file(f"{prune_dir}/prunning.bam", "prunning.bam", self.logger)
        else:
            link_file(f"{map_dir}/sample.clean.bam", "prunning.bam", self.logger)
        
        if check_file_exists("prunning.clusters.txt"):
            self.logger.log("✅ 找到现有文件 | Found existing: prunning.clusters.txt")
            return
        
        self.logger.run_command(
            f"ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e {self.config.enzyme} -k {self.config.chr_num}",
            "执行分区 | Executing partition"
        )
        
        self.logger.log_success("分区完成 | Partition completed")
        self.logger.log_step_time()
    
    def step3_5_extract(self):
        """Step 3.5: 提取矩阵 | Extract matrix"""
        if self.config.skip_steps["extract"]:
            self.logger.log("⏭️ 跳过矩阵提取 | Skipping matrix extraction")
            return
        
        self.logger.log_section("🧬 [Step 3.5] 提取接触矩阵 | Extracting contact matrix")
        extract_dir = self.config.directories["extract"]
        os.chdir(extract_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        map_dir = self.config.directories["mapping"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        link_file(f"{map_dir}/sample.clean.bam", "sample.clean.bam", self.logger)
        
        if check_file_exists("sample.clean.clm"):
            self.logger.log("✅ 找到现有文件 | Found existing: sample.clean.clm")
            return
        
        self.logger.run_command(
            f"allhic extract sample.clean.bam draft.asm.fasta --RE {self.config.enzyme}",
            "提取接触矩阵 | Extracting contact matrix"
        )
        
        self.logger.log_success("矩阵提取完成 | Matrix extraction completed")
        self.logger.log_step_time()
    
    def step4_rescue(self):
        """Step 4: 拯救 | Rescue"""
        if self.config.skip_steps["rescue"]:
            self.logger.log("⏭️ 跳过拯救 | Skipping rescue")
            return
        
        self.logger.log_section("🚑 [Step 4] 拯救未定位contigs | Rescuing unplaced contigs")
        rescue_dir = self.config.directories["rescue"]
        os.chdir(rescue_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        map_dir = self.config.directories["mapping"]
        part_dir = self.config.directories["partition"]
        extract_dir = self.config.directories["extract"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        link_file(f"{map_dir}/sample.clean.bam", "sample.clean.bam", self.logger)
        link_file(f"{part_dir}/prunning.clusters.txt", "prunning.clusters.txt", self.logger)
        link_file(f"{extract_dir}/sample.clean.counts_{self.config.enzyme}.txt", 
                 f"sample.clean.counts_{self.config.enzyme}.txt", self.logger)
        
        rescue_file = f"prunning.counts_{self.config.enzyme}.{self.config.chr_num}g1.txt"
        if check_file_exists(rescue_file):
            self.logger.log("✅ 找到现有拯救结果 | Found existing rescue results")
            return
        
        self.logger.run_command(
            f"ALLHiC_rescue -b sample.clean.bam -r draft.asm.fasta "
            f"-c prunning.clusters.txt -i sample.clean.counts_{self.config.enzyme}.txt",
            "执行拯救操作 | Executing rescue operation"
        )
        
        self.logger.log_success("拯救完成 | Rescue completed")
        self.logger.log_step_time()
    
    def step5_optimize(self):
        """Step 5: 优化 | Optimize"""
        if self.config.skip_steps["optimize"]:
            self.logger.log("⏭️ 跳过优化 | Skipping optimization")
            return
        
        self.logger.log_section("⚙️ [Step 5] 优化contig顺序 | Optimizing contig order")
        opt_dir = self.config.directories["optimize"]
        os.chdir(opt_dir)
        
        # 链接文件 | Link files
        extract_dir = self.config.directories["extract"]
        rescue_dir = self.config.directories["rescue"]
        part_dir = self.config.directories["partition"]
        
        link_file(f"{extract_dir}/sample.clean.clm", "sample.clean.clm", self.logger)
        
        # 确定源目录 | Determine source directory
        source_dir = rescue_dir
        if not os.path.exists(source_dir) or not any(f.startswith("prunning.counts_") for f in os.listdir(source_dir)):
            source_dir = part_dir
        
        # 复制分组文件 | Copy group files
        group_count = 0
        for gfile in os.listdir(source_dir):
            if gfile.startswith(f"prunning.counts_{self.config.enzyme}.{self.config.chr_num}g") and gfile.endswith(".txt"):
                group_num = gfile.split(f".{self.config.chr_num}g")[1].split(".txt")[0]
                dest_file = f"group{group_num}.txt"
                self.logger.log(f"  处理 | Processing: {gfile} -> {dest_file}")
                link_file(os.path.join(source_dir, gfile), dest_file, self.logger)
                group_count += 1
        
        # 确定需要优化的组 | Determine groups to optimize
        to_optimize = []
        for i in range(1, group_count + 1):
            if not check_file_exists(f"group{i}.tour"):
                to_optimize.append(str(i))
        
        if not to_optimize:
            self.logger.log("✅ 所有组已优化 | All groups already optimized")
            return
        
        self.logger.log(f"需要优化 {len(to_optimize)} 个组 | Need to optimize {len(to_optimize)} groups: {', '.join(to_optimize)}")
        
        # 并行优化 | Parallel optimization
        for group_num in to_optimize:
            self.logger.run_command(
                f"allhic optimize group{group_num}.txt sample.clean.clm > optimize_group{group_num}.log 2>&1",
                f"优化组 {group_num} | Optimizing group {group_num}"
            )
        
        # 检查结果 | Check results
        failed = 0
        for i in to_optimize:
            if not check_file_exists(f"group{i}.tour") or os.path.getsize(f"group{i}.tour") == 0:
                self.logger.log_error(f"  ❌ 组 {i}: 失败 | Group {i}: FAILED")
                failed += 1
        
        if failed > 0:
            raise RuntimeError(f"❌ {failed} 个组优化失败 | {failed} groups failed optimization")
        
        self.logger.log_success("优化完成 | Optimization completed")
        self.logger.log_step_time()
    
    def step6_build(self):
        """Step 6: 构建 | Build"""
        if self.config.skip_steps["build"]:
            self.logger.log("⏭️ 跳过构建 | Skipping build")
            return
        
        self.logger.log_section("🏗️ [Step 6] 构建最终组装 | Building final assembly")
        build_dir = self.config.directories["build"]
        os.chdir(build_dir)
        
        # 链接文件 | Link files
        data_dir = self.config.directories["data"]
        opt_dir = self.config.directories["optimize"]
        
        link_file(f"{data_dir}/draft.asm.fasta", "draft.asm.fasta", self.logger)
        
        for tour_file in os.listdir(opt_dir):
            if tour_file.endswith(".tour"):
                link_file(os.path.join(opt_dir, tour_file), tour_file, self.logger)
        
        if check_file_exists("groups.asm.fasta"):
            self.logger.log("✅ 找到现有文件 | Found existing: groups.asm.fasta")
        else:
            self.logger.run_command("ALLHiC_build draft.asm.fasta", "构建组装 | Building assembly")
        
        if not check_file_exists("groups.asm.fasta.fai"):
            self.logger.run_command("samtools faidx groups.asm.fasta", "构建FASTA索引 | Building FASTA index")
        
        self.logger.log_success("构建完成 | Build completed")
        self.logger.log_step_time()
    
    def step7_plot(self):
        """Step 7: 绘图 | Plot"""
        if self.config.skip_steps["plot"]:
            self.logger.log("⏭️ 跳过绘图 | Skipping plot")
            return
        
        self.logger.log_section("📊 [Step 7] 生成接触图谱 | Generating contact maps")
        plot_dir = self.config.directories["plot"]
        os.chdir(plot_dir)
        
        # 链接文件 | Link files
        map_dir = self.config.directories["mapping"]
        build_dir = self.config.directories["build"]
        
        link_file(f"{map_dir}/sample.clean.bam", "sample.clean.bam", self.logger)
        link_file(f"{build_dir}/groups.agp", "groups.agp", self.logger)
        link_file(f"{build_dir}/groups.asm.fasta.fai", "groups.asm.fasta.fai", self.logger)
        
        if not check_file_exists("chrn.list"):
            self.logger.run_command("cut -f1,2 groups.asm.fasta.fai > chrn.list", "生成染色体列表 | Generating chromosome list")
        
        self.logger.run_command(
            f"ALLHiC_plot -b sample.clean.bam -a groups.agp -l chrn.list "
            f"-m {self.config.min_bin_size} -s {self.config.bin_size} -o ./",
            "生成接触图谱 | Generating contact maps"
        )
        
        self.logger.log_success("绘图完成 | Plotting completed")
        self.logger.log_step_time()
