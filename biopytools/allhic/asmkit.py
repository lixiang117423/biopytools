"""
asmkit JBAT生成模块 🎯 | asmkit JBAT Generation Module
"""

import os
from utils import link_file, check_file_exists

class AsmkitProcessor:
    """asmkit处理器 | asmkit Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def run_asmkit_jbat(self):
        """Step 8: 生成asmkit JBAT文件 | Generate asmkit JBAT files"""
        if self.config.skip_steps["asmkit"]:
            self.logger.log("⏭️ 跳过asmkit JBAT生成 | Skipping asmkit JBAT generation")
            return
        
        self.logger.log_section("🎯 [Step 8] 使用asmkit生成JBAT文件 | Generating JBAT files with asmkit")
        self.logger.log("📌 策略: 使用ALLHiC结果 (sample.clean.bam + groups.agp)")
        
        jbat_dir = self.config.directories["jbat"]
        os.chdir(jbat_dir)
        
        # 1. 链接必要文件
        self._link_input_files()
        
        # 2. 检查asmkit工具
        asmkit_cmd = self._check_asmkit_tool()
        
        # 3. 检查是否已有结果
        if check_file_exists("groups.assembly") and check_file_exists("out.links"):
            self.logger.log("✅ JBAT文件已存在 | JBAT files already exist")
            if check_file_exists("groups.hic"):
                self.logger.log("✅ HiC文件也存在 | HiC file also exists")
            self.logger.log_step_time()
            return
        
        # 4. 生成links文件
        self._generate_links_file(asmkit_cmd)
        
        # 5. 生成assembly文件
        self._generate_assembly_file(asmkit_cmd)
        
        # 6. 生成.hic文件
        self._generate_hic_file()
        
        # 7. 输出摘要
        self._output_summary()
        
        self.logger.log_success("JBAT生成完成 | JBAT generation completed")
        self.logger.log_step_time()
    
    def _link_input_files(self):
        """链接输入文件 | Link input files"""
        self.logger.log("链接输入文件 | Linking input files...")
        
        # 链接BAM文件
        source_bam = os.path.join(self.config.directories["mapping"], "sample.clean.bam")
        if not check_file_exists(source_bam):
            raise RuntimeError(f"❌ BAM文件未找到 | BAM file not found: {source_bam}")
        link_file(source_bam, "sample.clean.bam", self.logger)
        
        # 链接AGP文件
        source_agp = os.path.join(self.config.directories["build"], "groups.agp")
        if not check_file_exists(source_agp):
            raise RuntimeError(f"❌ AGP文件未找到 | AGP file not found: {source_agp}")
        link_file(source_agp, "groups.agp", self.logger)
        
        self.logger.log_success("输入文件准备完成 | Input files prepared")
    
    def _check_asmkit_tool(self):
        """检查asmkit工具 | Check asmkit tool"""
        if os.path.exists(self.config.asmkit_path) and os.access(self.config.asmkit_path, os.X_OK):
            asmkit_cmd = self.config.asmkit_path
            self.logger.log_success(f"找到asmkit: {asmkit_cmd}")
        else:
            # 检查PATH中是否有asmkit
            import shutil
            asmkit_path = shutil.which("asmkit")
            if asmkit_path:
                asmkit_cmd = asmkit_path
                self.logger.log_success(f"在PATH中找到asmkit: {asmkit_cmd}")
            else:
                raise RuntimeError(f"❌ asmkit未找到！| asmkit not found! 期望位置 | Expected location: {self.config.asmkit_path}")
        
        return asmkit_cmd
    
    def _generate_links_file(self, asmkit_cmd):
        """生成links文件 | Generate links file"""
        if check_file_exists("out.links"):
            return
        
        self.logger.log_section("Step 8.1: 从BAM生成links文件 | Generating links file from BAM")
        
        import time
        start_time = time.time()
        
        success = self.logger.run_command(
            f"{asmkit_cmd} bam2links sample.clean.bam out.links",
            "生成links文件 | Generating links file"
        )
        
        if success and check_file_exists("out.links"):
            end_time = time.time()
            elapsed = int(end_time - start_time)
            link_count = self._count_lines("out.links")
            file_size = self._get_file_size("out.links")
            self.logger.log_success(f"成功生成out.links ({link_count} 行, {file_size})")
            self.logger.log(f"   耗时 | Time: {self.logger.format_time(elapsed)}")
        else:
            raise RuntimeError("❌ links文件生成失败 | Failed to generate links file")
    
    def _generate_assembly_file(self, asmkit_cmd):
        """生成assembly文件 | Generate assembly file"""
        if check_file_exists("groups.assembly"):
            return
        
        self.logger.log_section("Step 8.2: 从AGP生成assembly文件 | Generating assembly file from AGP")
        
        import time
        start_time = time.time()
        
        success = self.logger.run_command(
            f"{asmkit_cmd} agp2assembly groups.agp groups.assembly",
            "生成assembly文件 | Generating assembly file"
        )
        
        if success and check_file_exists("groups.assembly"):
            end_time = time.time()
            elapsed = int(end_time - start_time)
            self.logger.log_success(f"成功生成groups.assembly")
            self.logger.log(f"   耗时 | Time: {self.logger.format_time(elapsed)}")
        else:
            raise RuntimeError("❌ assembly文件生成失败 | Failed to generate assembly file")
    
    def _generate_hic_file(self):
        """生成.hic文件 | Generate .hic file"""
        if check_file_exists("groups.hic"):
            return
        
        self.logger.log_section("Step 8.3: 使用visualizer生成.hic文件 | Generating .hic file using visualizer")
        
        if not check_file_exists(self.config.assembly_visualizer):
            self.logger.log_warning(f"Visualizer未找到: {self.config.assembly_visualizer}")
            self.logger.log("   跳过.hic生成。您可以手动运行。| Skipping .hic generation. You can run it manually.")
            return
        
        import time
        start_time = time.time()
        
        success = self.logger.run_command(
            f"bash {self.config.assembly_visualizer} groups.assembly out.links",
            "生成.hic文件 | Generating .hic file"
        )
        
        if success and check_file_exists("groups.hic"):
            end_time = time.time()
            elapsed = int(end_time - start_time)
            file_size = self._get_file_size("groups.hic")
            self.logger.log_success(f"成功生成groups.hic ({file_size})")
            self.logger.log(f"   耗时 | Time: {self.logger.format_time(elapsed)}")
        else:
            self.logger.log_warning("⚠️ Visualizer失败。请检查日志。| Visualizer failed. Check logs.")
    
    def _output_summary(self):
        """输出摘要 | Output summary"""
        self.logger.log_section("📦 JBAT生成完成！| JBAT Generation Complete!")
        
        # 列出生成的文件
        import os
        output_files = []
        for filename in ["groups.assembly", "out.links", "groups.hic"]:
            if check_file_exists(filename):
                output_files.append(filename)
        
        if output_files:
            self.logger.log(f"输出目录 | Output directory: {os.getcwd()}")
            for filename in output_files:
                size = self._get_file_size(filename)
                self.logger.log(f"   📄 {filename} ({size})")
        
        # 下一步建议
        self.logger.log("")
        self.logger.log("💡 下一步建议 | Next Steps (Manual Curation):")
        self.logger.log("   1. 下载'groups.assembly'和'groups.hic'(如果生成) | 1. Download 'groups.assembly' and 'groups.hic' (if generated)")
        self.logger.log("   2. 加载到Juicebox Assembly Tools | 2. Load into Juicebox Assembly Tools")
        self.logger.log("   3. 编辑并保存为'reviewed.assembly' | 3. Edit and save as 'reviewed.assembly'")
    
    def _count_lines(self, filename):
        """计算文件行数 | Count file lines"""
        try:
            with open(filename) as f:
                return sum(1 for _ in f)
        except:
            return 0
    
    def _get_file_size(self, filename):
        """获取文件大小 | Get file size"""
        try:
            size = os.path.getsize(filename)
            for unit in ['B', 'KB', 'MB', 'GB']:
                if size < 1024.0:
                    return f"{size:.1f}{unit}"
                size /= 1024.0
            return f"{size:.1f}TB"
        except:
            return "0B"
