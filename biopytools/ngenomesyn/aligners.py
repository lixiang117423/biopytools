"""
比对器模块 | Aligners Module
"""

import os
import subprocess
from typing import Dict, List, Optional
from .data_processing import GenomeProcessor

class BaseAligner:
    """基础比对器类 | Base Aligner Class"""
    
    def __init__(self, logger, threads: int = 16, output_dir: str = "./"):
        self.logger = logger
        self.threads = threads
        self.output_dir = output_dir
        self.command_log = []  # 存储执行的命令
    
    def log_command(self, command: str, working_dir: str = None, description: str = ""):
        """记录执行的命令到日志 | Log executed command"""
        if working_dir is None:
            working_dir = self.output_dir
        
        self.logger.info("🔗" + "=" * 80)
        self.logger.info(f"🔗 执行命令 | Executing command: {description}")
        self.logger.info(f"📁 工作目录 | Working directory: {working_dir}")
        self.logger.info(f"💻 命令内容 | Command: {command}")
        self.logger.info("🔗" + "=" * 80)
        
        # 保存到命令历史
        self.command_log.append({
            "description": description,
            "command": command,
            "working_dir": working_dir
        })
    
    def save_command_script(self, script_name: str = "run_commands.sh"):
        """保存所有执行过的命令到脚本文件 | Save all executed commands to script file"""
        script_path = os.path.join(self.output_dir, script_name)
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 🚀 自动生成的命令脚本 | Auto-generated command script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")
            
            for i, cmd_info in enumerate(self.command_log, 1):
                f.write(f"# 📋 步骤 {i}: {cmd_info['description']}\n")
                f.write(f"echo \"🔄 执行步骤 {i}: {cmd_info['description']}\"\n")
                f.write(f"cd \"{cmd_info['working_dir']}\"\n")
                f.write(f"{cmd_info['command']}\n")
                f.write("if [ $? -ne 0 ]; then\n")
                f.write(f"    echo \"❌ 错误: 步骤 {i} 执行失败\"\n")
                f.write("    exit 1\n")
                f.write("fi\n\n")
        
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 命令脚本已保存 | Command script saved: {script_path}")
        self.logger.info("🚀 您可以使用以下命令手动执行 | You can manually execute with:")
        self.logger.info(f"bash {script_path}")
        
        return script_path
    
    def align(self, genome1: Dict, genome2: Dict, selected_chromosomes: Optional[List[int]] = None, **params) -> str:
        """执行比对 | Perform alignment"""
        raise NotImplementedError

class Minimap2Aligner(BaseAligner):
    """Minimap2比对器 | Minimap2 Aligner"""
    
    def __init__(self, logger, threads: int = 16, output_dir: str = "./"):
        super().__init__(logger, threads, output_dir)
        self.logger.info(f"🔧 Minimap2比对器初始化 | Minimap2 aligner initialized with {threads} threads")
    
    def align(self, genome1: Dict, genome2: Dict, selected_chromosomes: Optional[List[int]] = None, **params) -> str:
        """执行Minimap2比对 | Perform Minimap2 alignment"""
        preset = params.get("preset", "asm5")
        min_length = params.get("min_length", 5000)
        
        # 确定输入文件
        genome1_file = genome1['file']
        genome2_file = genome2['file']
        
        # 如果指定了染色体，先过滤序列文件
        if selected_chromosomes:
            self.logger.info(f"🧬 过滤染色体 | Filtering chromosomes: {selected_chromosomes}")
            
            filtered_genome1 = os.path.join(self.output_dir, f"{genome1['name']}_filtered.fa")
            filtered_genome2 = os.path.join(self.output_dir, f"{genome2['name']}_filtered.fa")
            
            processor = GenomeProcessor(self.logger)
            
            # 过滤第一个基因组
            if not os.path.exists(filtered_genome1):
                processor.filter_sequences_by_chromosome(genome1['file'], filtered_genome1, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome1}")
            
            # 过滤第二个基因组
            if not os.path.exists(filtered_genome2):
                processor.filter_sequences_by_chromosome(genome2['file'], filtered_genome2, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome2}")
            
            genome1_file = filtered_genome1
            genome2_file = filtered_genome2
            
            # 更新输出文件名以反映过滤状态
            chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))
        else:
            chr_suffix = ""
        
        # 生成输出文件的绝对路径
        paf_file = os.path.abspath(os.path.join(self.output_dir, f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}.paf"))
        link_file = os.path.abspath(os.path.join(self.output_dir, f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}.link"))
        
        # 确保输出目录存在
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.logger.info(f"🚀 开始Minimap2比对 | Starting Minimap2 alignment: {genome1['name']} vs {genome2['name']}")
        if selected_chromosomes:
            self.logger.info(f"🧬 分析染色体 | Analyzing chromosomes: {selected_chromosomes}")
        self.logger.info(f"📄 输出PAF文件 | Output PAF file: {paf_file}")
        self.logger.info(f"🔗 输出LINK文件 | Output LINK file: {link_file}")
        
        try:
            # 运行minimap2 - 使用传入的线程数
            cmd_minimap2 = f"minimap2 -x {preset} -t {self.threads} \"{genome2_file}\" \"{genome1_file}\" > \"{paf_file}\""
            
            # 记录命令到日志
            self.log_command(
                command=cmd_minimap2,
                working_dir=self.output_dir,
                description=f"🧬 Minimap2比对: {genome1['name']} vs {genome2['name']}"
            )
            
            result = subprocess.run(
                cmd_minimap2, 
                shell=True, 
                capture_output=True, 
                text=True
            )
            
            if result.returncode != 0:
                self.logger.error(f"❌ Minimap2执行失败 | Minimap2 failed: {result.stderr}")
                # 即使失败也保存命令脚本
                self.save_command_script(f"minimap2_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ Minimap2执行失败 | Minimap2 failed")
            
            # 检查PAF文件是否成功创建
            if not os.path.exists(paf_file):
                self.logger.error(f"❌ PAF文件未创建 | PAF file not created: {paf_file}")
                self.save_command_script(f"minimap2_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ PAF文件未创建 | PAF file not created")
            
            # 检查PAF文件是否为空
            if os.path.getsize(paf_file) == 0:
                self.logger.warning(f"⚠️ PAF文件为空，可能没有找到比对结果 | PAF file is empty, no alignments found")
                # 不抛出异常，因为这可能是正常的生物学现象
            
            self.logger.info(f"✅ Minimap2比对完成 | Minimap2 alignment completed: {paf_file}")
            
            # 转换为link格式
            if self._check_conversion_script():
                cmd_convert = f"GetTwoGenomeSyn.pl Paf2Link \"{os.path.basename(paf_file)}\" {min_length} \"{os.path.basename(link_file)}\""
                
                # 记录转换命令到日志
                self.log_command(
                    command=cmd_convert,
                    working_dir=self.output_dir,
                    description=f"🔄 PAF转LINK格式转换: {os.path.basename(paf_file)} -> {os.path.basename(link_file)}"
                )
                
                result = subprocess.run(
                    cmd_convert, 
                    shell=True, 
                    capture_output=True, 
                    text=True,
                    cwd=self.output_dir
                )
                
                if result.returncode != 0:
                    self.logger.warning(f"⚠️ PAF转换失败，将使用PAF文件 | PAF conversion failed, will use PAF file: {result.stderr}")
                    # 保存命令脚本以便调试
                    self.save_command_script(f"alignment_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                    return paf_file
                else:
                    self.logger.info(f"✅ LINK文件转换完成 | LINK file conversion completed: {link_file}")
                    # 保存成功的命令脚本
                    self.save_command_script(f"alignment_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                    return link_file
            else:
                self.logger.warning(f"⚠️ GetTwoGenomeSyn.pl脚本不可用，返回PAF文件 | GetTwoGenomeSyn.pl script not available, returning PAF file")
                # 保存命令脚本
                self.save_command_script(f"alignment_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                return paf_file
            
        except Exception as e:
            self.logger.error(f"💥 Minimap2比对失败 | Minimap2 alignment failed: {e}")
            # 无论如何都保存命令脚本以便调试
            self.save_command_script(f"alignment_commands_{genome1['name']}_vs_{genome2['name']}_FAILED.sh")
            raise
    
    def _check_conversion_script(self) -> bool:
        """检查转换脚本是否可用 | Check if conversion script is available"""
        try:
            result = subprocess.run(
                ["which", "GetTwoGenomeSyn.pl"],
                capture_output=True,
                text=True,
                timeout=5
            )

            # 记录检查命令
            self.log_command(
                command="which GetTwoGenomeSyn.pl",
                working_dir=self.output_dir,
                description="🔍 检查GetTwoGenomeSyn.pl脚本是否可用"
            )

            return result.returncode == 0
        except:
            return False

    def align_for_syri(self, genome1: Dict, genome2: Dict, selected_chromosomes: Optional[List[int]] = None, **params) -> str:
        """执行Minimap2比对并生成SAM格式用于SyRI分析 | Perform Minimap2 alignment and generate SAM format for SyRI

        SyRI需要SAM格式的输入，使用不同的参数：
        - 使用 -ax asm5 而不是 -x asm5（生成SAM而不是PAF）
        - 添加 --eqx 参数生成扩展的CIGAR字符串

        Args:
            genome1: 参考基因组信息
            genome2: 查询基因组信息
            selected_chromosomes: 选择的染色体列表
            **params: 额外参数

        Returns:
            sam_file: 生成的SAM文件路径
        """
        preset = params.get("preset", "asm5")

        # 确定输入文件
        genome1_file = genome1['file']
        genome2_file = genome2['file']

        # 如果指定了染色体，先过滤序列文件
        if selected_chromosomes:
            self.logger.info(f"🧬 过滤染色体 | Filtering chromosomes: {selected_chromosomes}")

            filtered_genome1 = os.path.join(self.output_dir, f"{genome1['name']}_filtered.fa")
            filtered_genome2 = os.path.join(self.output_dir, f"{genome2['name']}_filtered.fa")

            processor = GenomeProcessor(self.logger)

            # 过滤第一个基因组
            if not os.path.exists(filtered_genome1):
                processor.filter_sequences_by_chromosome(genome1['file'], filtered_genome1, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome1}")

            # 过滤第二个基因组
            if not os.path.exists(filtered_genome2):
                processor.filter_sequences_by_chromosome(genome2['file'], filtered_genome2, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome2}")

            genome1_file = filtered_genome1
            genome2_file = filtered_genome2

            # 更新输出文件名以反映过滤状态
            chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))
        else:
            chr_suffix = ""

        # 生成SAM输出文件的绝对路径
        sam_file = os.path.abspath(os.path.join(self.output_dir, f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}.sam"))

        # 确保输出目录存在
        os.makedirs(self.output_dir, exist_ok=True)

        self.logger.info(f"🚀 开始Minimap2比对 (SAM格式) | Starting Minimap2 alignment (SAM format): {genome1['name']} vs {genome2['name']}")
        if selected_chromosomes:
            self.logger.info(f"🧬 分析染色体 | Analyzing chromosomes: {selected_chromosomes}")
        self.logger.info(f"📄 输出SAM文件 | Output SAM file: {sam_file}")

        try:
            # 运行minimap2生成SAM格式 - 使用 -ax 和 --eqx 参数
            # -ax asm5: 生成SAM格式输出（用于长序列比对）
            # --eqx: 生成扩展的CIGAR字符串（包含=和X操作符）
            cmd_minimap2 = f"minimap2 -ax {preset} --eqx -t {self.threads} \"{genome2_file}\" \"{genome1_file}\" > \"{sam_file}\""

            # 记录命令到日志
            self.log_command(
                command=cmd_minimap2,
                working_dir=self.output_dir,
                description=f"🧬 Minimap2比对 (SAM格式 for SyRI): {genome1['name']} vs {genome2['name']}"
            )

            result = subprocess.run(
                cmd_minimap2,
                shell=True,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"❌ Minimap2执行失败 | Minimap2 failed: {result.stderr}")
                # 即使失败也保存命令脚本
                self.save_command_script(f"minimap2_syri_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ Minimap2执行失败 | Minimap2 failed")

            # 检查SAM文件是否成功创建
            if not os.path.exists(sam_file):
                self.logger.error(f"❌ SAM文件未创建 | SAM file not created: {sam_file}")
                self.save_command_script(f"minimap2_syri_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ SAM文件未创建 | SAM file not created")

            # 检查SAM文件是否为空
            if os.path.getsize(sam_file) == 0:
                self.logger.warning(f"⚠️ SAM文件为空，可能没有找到比对结果 | SAM file is empty, no alignments found")
                # 不抛出异常，因为这可能是正常的生物学现象

            self.logger.info(f"✅ Minimap2比对完成 (SAM格式) | Minimap2 alignment completed (SAM format): {sam_file}")
            # 保存成功的命令脚本
            self.save_command_script(f"minimap2_syri_commands_{genome1['name']}_vs_{genome2['name']}.sh")

            return sam_file

        except Exception as e:
            self.logger.error(f"💥 Minimap2比对失败 | Minimap2 alignment failed: {e}")
            # 无论如何都保存命令脚本以便调试
            self.save_command_script(f"minimap2_syri_commands_{genome1['name']}_vs_{genome2['name']}_FAILED.sh")
            raise

class MummerAligner(BaseAligner):
    """MUMmer比对器 | MUMmer Aligner"""
    
    def __init__(self, logger, threads: int = 16, output_dir: str = "./"):
        super().__init__(logger, threads, output_dir)
        self.logger.info(f"🧮 MUMmer比对器初始化 | MUMmer aligner initialized")
    
    def align(self, genome1: Dict, genome2: Dict, selected_chromosomes: Optional[List[int]] = None, **params) -> str:
        """执行MUMmer比对 | Perform MUMmer alignment"""
        prefix = params.get("prefix", "mummer_out")
        min_match = params.get("min_match", 20)
        min_cluster = params.get("min_cluster", 65)
        max_gap = params.get("max_gap", 90)
        match_type = params.get("match_type", "mumreference")  # mum, mumreference, maxmatch
        
        # 确定输入文件
        genome1_file = genome1['file']
        genome2_file = genome2['file']
        
        # 如果指定了染色体，先过滤序列文件
        if selected_chromosomes:
            self.logger.info(f"🧬 过滤染色体 | Filtering chromosomes: {selected_chromosomes}")
            
            filtered_genome1 = os.path.join(self.output_dir, f"{genome1['name']}_filtered.fa")
            filtered_genome2 = os.path.join(self.output_dir, f"{genome2['name']}_filtered.fa")
            
            processor = GenomeProcessor(self.logger)
            
            # 过滤第一个基因组
            if not os.path.exists(filtered_genome1):
                processor.filter_sequences_by_chromosome(genome1['file'], filtered_genome1, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome1}")
            
            # 过滤第二个基因组
            if not os.path.exists(filtered_genome2):
                processor.filter_sequences_by_chromosome(genome2['file'], filtered_genome2, selected_chromosomes)
            else:
                self.logger.info(f"♻️ 使用已存在的过滤文件 | Using existing filtered file: {filtered_genome2}")
            
            genome1_file = filtered_genome1
            genome2_file = filtered_genome2
            
            # 更新输出文件名以反映过滤状态
            chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))
        else:
            chr_suffix = ""
        
        # 生成输出文件的绝对路径
        output_prefix = f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}_mummer"
        delta_file = os.path.abspath(os.path.join(self.output_dir, f"{output_prefix}.delta"))
        coords_file = os.path.abspath(os.path.join(self.output_dir, f"{output_prefix}.coords"))
        
        # 确保输出目录存在
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.logger.info(f"🚀 开始MUMmer比对 | Starting MUMmer alignment: {genome1['name']} vs {genome2['name']}")
        if selected_chromosomes:
            self.logger.info(f"🧬 分析染色体 | Analyzing chromosomes: {selected_chromosomes}")
        self.logger.info(f"📄 输出Delta文件 | Output Delta file: {delta_file}")
        self.logger.info(f"📊 输出Coords文件 | Output Coords file: {coords_file}")
        
        try:
            # 构建nucmer命令
            cmd_nucmer = f"nucmer --{match_type} -l {min_match} -c {min_cluster} -g {max_gap} -p \"{output_prefix}\" \"{genome2_file}\" \"{genome1_file}\""
            
            # 记录命令到日志
            self.log_command(
                command=cmd_nucmer,
                working_dir=self.output_dir,
                description=f"🧮 MUMmer nucmer比对: {genome1['name']} vs {genome2['name']}"
            )
            
            result = subprocess.run(
                cmd_nucmer, 
                shell=True, 
                capture_output=True, 
                text=True,
                cwd=self.output_dir
            )
            
            if result.returncode != 0:
                self.logger.error(f"❌ MUMmer nucmer执行失败 | MUMmer nucmer failed: {result.stderr}")
                # 即使失败也保存命令脚本
                self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ MUMmer nucmer执行失败 | MUMmer nucmer failed")
            
            # 检查Delta文件是否成功创建
            if not os.path.exists(delta_file):
                self.logger.error(f"❌ Delta文件未创建 | Delta file not created: {delta_file}")
                self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                raise RuntimeError(f"❌ Delta文件未创建 | Delta file not created")
            
            # 检查Delta文件是否为空
            if os.path.getsize(delta_file) == 0:
                self.logger.warning(f"⚠️ Delta文件为空，可能没有找到比对结果 | Delta file is empty, no alignments found")
                # 不抛出异常，因为这可能是正常的生物学现象
            
            self.logger.info(f"✅ MUMmer nucmer比对完成 | MUMmer nucmer alignment completed: {delta_file}")
            
            # 转换为coords格式以便后续处理
            if self._check_show_coords():
                cmd_coords = f"show-coords -r -c -l \"{os.path.basename(delta_file)}\" > \"{os.path.basename(coords_file)}\""
                
                # 记录转换命令到日志
                self.log_command(
                    command=cmd_coords,
                    working_dir=self.output_dir,
                    description=f"🔄 Delta转Coords格式转换: {os.path.basename(delta_file)} -> {os.path.basename(coords_file)}"
                )
                
                result = subprocess.run(
                    cmd_coords, 
                    shell=True, 
                    capture_output=True, 
                    text=True,
                    cwd=self.output_dir
                )
                
                if result.returncode != 0:
                    self.logger.warning(f"⚠️ Delta转换失败，将使用Delta文件 | Delta conversion failed, will use Delta file: {result.stderr}")
                    # 保存命令脚本以便调试
                    self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                    return delta_file
                else:
                    self.logger.info(f"✅ Coords文件转换完成 | Coords file conversion completed: {coords_file}")
                    # 保存成功的命令脚本
                    self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                    return coords_file
            else:
                self.logger.warning(f"⚠️ show-coords工具不可用，返回Delta文件 | show-coords tool not available, returning Delta file")
                # 保存命令脚本
                self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}.sh")
                return delta_file
            
        except Exception as e:
            self.logger.error(f"💥 MUMmer比对失败 | MUMmer alignment failed: {e}")
            # 无论如何都保存命令脚本以便调试
            self.save_command_script(f"mummer_commands_{genome1['name']}_vs_{genome2['name']}_FAILED.sh")
            raise
    
    def _check_show_coords(self) -> bool:
        """检查show-coords工具是否可用 | Check if show-coords tool is available"""
        try:
            result = subprocess.run(
                ["which", "show-coords"],
                capture_output=True,
                text=True,
                timeout=5
            )

            # 记录检查命令
            self.log_command(
                command="which show-coords",
                working_dir=self.output_dir,
                description="🔍 检查show-coords工具是否可用"
            )

            return result.returncode == 0
        except:
            return False

class SyriProcessor:
    """SyRI结构变异分析处理器 | SyRI Structural Variation Processor"""

    def __init__(self, logger, syri_bin: Optional[str] = None, output_dir: str = "./"):
        self.logger = logger
        self.syri_bin = syri_bin
        self.output_dir = output_dir
        self.command_log = []
        self.logger.info("🧬 SyRI处理器初始化 | SyRI processor initialized")

    def log_command(self, command: str, working_dir: str = None, description: str = ""):
        """记录执行的命令到日志 | Log executed command"""
        if working_dir is None:
            working_dir = self.output_dir

        self.logger.info("🧬" + "=" * 80)
        self.logger.info(f"🧬 执行命令 | Executing command: {description}")
        self.logger.info(f"📁 工作目录 | Working directory: {working_dir}")
        self.logger.info(f"💻 命令内容 | Command: {command}")
        self.logger.info("🧬" + "=" * 80)

        self.command_log.append({
            "description": description,
            "command": command,
            "working_dir": working_dir
        })

    def find_syri_binary(self) -> Optional[str]:
        """查找SyRI二进制文件 | Find SyRI binary"""
        import shutil

        # 如果已经指定，直接返回
        if self.syri_bin and os.path.exists(self.syri_bin):
            return self.syri_bin

        # 尝试从PATH查找syri命令
        syri_path = shutil.which("syri")
        if syri_path:
            return syri_path

        # 尝试常见安装路径
        possible_paths = [
            "/usr/local/bin/syri",
            os.path.expanduser("~/bin/syri"),
            os.path.expanduser("~/.local/bin/syri"),
        ]

        for path in possible_paths:
            if os.path.exists(path):
                return path

        return None

    def run_syri_from_minimap2(
        self,
        genome1: Dict,
        genome2: Dict,
        sam_file: str,
        min_length: int = 5000,
        selected_chromosomes: Optional[List[int]] = None
    ) -> str:
        """从Minimap2 SAM结果运行SyRI分析 | Run SyRI analysis from Minimap2 SAM output

        Args:
            genome1: 参考基因组信息
            genome2: 查询基因组信息
            sam_file: Minimap2生成的SAM文件
            min_length: 最小比对长度
            selected_chromosomes: 选择的染色体列表

        Returns:
            link_file: 生成的LINK文件路径
        """
        # 确定输出文件名
        chr_suffix = ""
        if selected_chromosomes:
            chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))

        base_name = f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}"
        syri_out_file = os.path.abspath(os.path.join(self.output_dir, f"{base_name}.syri.out"))
        syri_vcf_file = os.path.abspath(os.path.join(self.output_dir, f"{base_name}.syri.vcf"))
        link_file = os.path.abspath(os.path.join(self.output_dir, f"{base_name}.link"))

        # 查找SyRI二进制文件
        syri_bin = self.find_syri_binary()
        if not syri_bin:
            raise RuntimeError("❌ 找不到SyRI程序，请安装或指定路径 | SyRI binary not found, please install or specify path")

        self.logger.info(f"🧬 开始SyRI分析 | Starting SyRI analysis")
        self.logger.info(f"📄 输入SAM文件 | Input SAM file: {sam_file}")
        self.logger.info(f"🧬 参考基因组 | Reference genome: {genome1['file']}")
        self.logger.info(f"🧬 查询基因组 | Query genome: {genome2['file']}")

        try:
            # 运行SyRI - 从SAM文件
            # 使用syri命令直接调用（不是python3）
            cmd_syri = f"\"{syri_bin}\" " \
                      f"-c \"{sam_file}\" " \
                      f"-r \"{genome1['file']}\" " \
                      f"-q \"{genome2['file']}\" " \
                      f"-F B --prefix {base_name}"

            self.log_command(
                command=cmd_syri,
                working_dir=self.output_dir,
                description=f"🧬 SyRI结构变异分析: {genome1['name']} vs {genome2['name']}"
            )

            result = subprocess.run(
                cmd_syri,
                shell=True,
                capture_output=True,
                text=True,
                cwd=self.output_dir
            )

            if result.returncode != 0:
                self.logger.error(f"❌ SyRI执行失败 | SyRI failed: {result.stderr}")
                raise RuntimeError(f"❌ SyRI执行失败 | SyRI failed")

            # 检查输出文件
            if not os.path.exists(syri_out_file):
                # SyRI有时会生成带不同格式的文件
                possible_files = [
                    f"{base_name}.syri.out",
                    f"{base_name}_syri.out",
                    f"{base_name}syri.out",  # 实际观察到的格式
                    "syri.out"
                ]
                for pf in possible_files:
                    full_path = os.path.join(self.output_dir, pf)
                    if os.path.exists(full_path):
                        syri_out_file = full_path
                        break
                else:
                    raise RuntimeError(f"❌ SyRI输出文件未找到 | SyRI output file not found: {syri_out_file}")

            self.logger.info(f"✅ SyRI分析完成 | SyRI analysis completed: {syri_out_file}")

            # 使用GetTwoGenomeSyn.pl转换为LINK格式
            link_file = self._convert_syri_to_link(syri_out_file, link_file, min_length)

            return link_file

        except Exception as e:
            self.logger.error(f"💥 SyRI分析失败 | SyRI analysis failed: {e}")
            raise

    def run_syri_from_mummer_coords(
        self,
        genome1: Dict,
        genome2: Dict,
        coords_file: str,
        min_length: int = 5000,
        selected_chromosomes: Optional[List[int]] = None
    ) -> str:
        """从MUMmer coords结果运行SyRI分析 | Run SyRI analysis from MUMmer coords output

        Args:
            genome1: 参考基因组信息
            genome2: 查询基因组信息
            coords_file: MUMmer生成的coords文件
            min_length: 最小比对长度
            selected_chromosomes: 选择的染色体列表

        Returns:
            link_file: 生成的LINK文件路径
        """
        # 确定输出文件名
        chr_suffix = ""
        if selected_chromosomes:
            chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))

        base_name = f"{genome1['name']}_vs_{genome2['name']}{chr_suffix}"
        syri_out_file = os.path.abspath(os.path.join(self.output_dir, f"{base_name}.syri.out"))
        link_file = os.path.abspath(os.path.join(self.output_dir, f"{base_name}.link"))

        # 查找SyRI二进制文件
        syri_bin = self.find_syri_binary()
        if not syri_bin:
            raise RuntimeError("❌ 找不到SyRI程序，请安装或指定路径 | SyRI binary not found, please install or specify path")

        self.logger.info(f"🧬 开始SyRI分析 | Starting SyRI analysis")
        self.logger.info(f"📊 输入Coords文件 | Input Coords file: {coords_file}")
        self.logger.info(f"🧬 参考基因组 | Reference genome: {genome1['file']}")
        self.logger.info(f"🧬 查询基因组 | Query genome: {genome2['file']}")

        try:
            # 运行SyRI - 从Coords文件
            # 使用syri命令直接调用（不是python3）
            cmd_syri = f"\"{syri_bin}\" " \
                      f"-c \"{coords_file}\" " \
                      f"-r \"{genome1['file']}\" " \
                      f"-q \"{genome2['file']}\" " \
                      f"-F B --prefix {base_name}"

            self.log_command(
                command=cmd_syri,
                working_dir=self.output_dir,
                description=f"🧬 SyRI结构变异分析 (from coords): {genome1['name']} vs {genome2['name']}"
            )

            result = subprocess.run(
                cmd_syri,
                shell=True,
                capture_output=True,
                text=True,
                cwd=self.output_dir
            )

            if result.returncode != 0:
                self.logger.error(f"❌ SyRI执行失败 | SyRI failed: {result.stderr}")
                raise RuntimeError(f"❌ SyRI执行失败 | SyRI failed")

            # 检查输出文件
            if not os.path.exists(syri_out_file):
                # SyRI有时会生成带不同格式的文件
                possible_files = [
                    f"{base_name}.syri.out",
                    f"{base_name}_syri.out",
                    f"{base_name}syri.out",  # 实际观察到的格式
                    "syri.out"
                ]
                for pf in possible_files:
                    full_path = os.path.join(self.output_dir, pf)
                    if os.path.exists(full_path):
                        syri_out_file = full_path
                        break
                else:
                    raise RuntimeError(f"❌ SyRI输出文件未找到 | SyRI output file not found: {syri_out_file}")

            self.logger.info(f"✅ SyRI分析完成 | SyRI analysis completed: {syri_out_file}")

            # 使用GetTwoGenomeSyn.pl转换为LINK格式
            link_file = self._convert_syri_to_link(syri_out_file, link_file, min_length)

            return link_file

        except Exception as e:
            self.logger.error(f"💥 SyRI分析失败 | SyRI analysis failed: {e}")
            raise

    def _convert_syri_to_link(self, syri_out_file: str, link_file: str, min_length: int) -> str:
        """使用GetTwoGenomeSyn.pl将SyRI输出转换为LINK格式 | Convert SyRI output to LINK format

        Args:
            syri_out_file: SyRI生成的输出文件
            link_file: 输出的LINK文件路径
            min_length: 最小比对长度

        Returns:
            link_file: 生成的LINK文件路径
        """
        self.logger.info(f"🔄 转换SyRI输出为LINK格式 | Converting SyRI output to LINK format")

        try:
            # 使用GetTwoGenomeSyn.pl的Syri2Link模式
            cmd_convert = f"GetTwoGenomeSyn.pl Syri2Link " \
                          f"\"{os.path.basename(syri_out_file)}\" " \
                          f"{min_length} " \
                          f"\"{os.path.basename(link_file)}\""

            self.log_command(
                command=cmd_convert,
                working_dir=self.output_dir,
                description=f"🔄 SyRI转LINK格式转换: {os.path.basename(syri_out_file)} -> {os.path.basename(link_file)}"
            )

            result = subprocess.run(
                cmd_convert,
                shell=True,
                capture_output=True,
                text=True,
                cwd=self.output_dir
            )

            if result.returncode != 0:
                self.logger.warning(f"⚠️ SyRI转换失败 | SyRI conversion failed: {result.stderr}")
                self.logger.warning(f"📄 stderr输出 | stderr output: {result.stderr}")
                raise RuntimeError(f"❌ SyRI转LINK格式失败 | SyRI to LINK conversion failed")

            # 检查LINK文件是否生成
            if not os.path.exists(link_file):
                raise RuntimeError(f"❌ LINK文件未生成 | LINK file not generated: {link_file}")

            self.logger.info(f"✅ LINK文件生成完成 | LINK file generated: {link_file}")

            return link_file

        except Exception as e:
            self.logger.error(f"💥 SyRI转换失败 | SyRI conversion failed: {e}")
            raise

    def save_command_script(self, script_name: str = "syri_commands.sh"):
        """保存所有执行过的命令到脚本文件 | Save all executed commands to script file"""
        script_path = os.path.join(self.output_dir, script_name)

        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 🧬 SyRI分析命令脚本 | SyRI Analysis Command Script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")

            for i, cmd_info in enumerate(self.command_log, 1):
                f.write(f"# 📋 步骤 {i}: {cmd_info['description']}\n")
                f.write(f"echo \"🔄 执行步骤 {i}: {cmd_info['description']}\"\n")
                f.write(f"cd \"{cmd_info['working_dir']}\"\n")
                f.write(f"{cmd_info['command']}\n")
                f.write("if [ $? -ne 0 ]; then\n")
                f.write(f"    echo \"❌ 错误: 步骤 {i} 执行失败\"\n")
                f.write("    exit 1\n")
                f.write("fi\n\n")

        # 设置可执行权限
        os.chmod(script_path, 0o755)

        self.logger.info(f"📝 SyRI命令脚本已保存 | SyRI command script saved: {script_path}")

        return script_path

class AlignerFactory:
    """比对器工厂 | Aligner Factory"""
    
    @staticmethod
    def create_aligner(aligner_type: str, logger, threads: int = 16, output_dir: str = "./") -> BaseAligner:
        """创建比对器实例 | Create aligner instance"""
        aligners = {
            "minimap2": Minimap2Aligner,
            "mummer": MummerAligner
        }
        
        if aligner_type not in aligners:
            raise ValueError(f"❌ 不支持的比对器类型 | Unsupported aligner type: {aligner_type}")
        
        return aligners[aligner_type](logger, threads, output_dir)