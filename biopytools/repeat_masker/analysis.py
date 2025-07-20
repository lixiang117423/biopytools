"""
重复序列分析核心模块 | Repeat Sequence Analysis Core Module
"""

import os
from .utils import CommandRunner

class RepeatModelerRunner:
    """RepeatModeler运行器 | RepeatModeler Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_database(self) -> bool:
        """构建RepeatModeler数据库 | Build RepeatModeler database"""
        genome_file = self.config.genome_file
        species = self.config.species
        database_name = f"{species}_genome"
        
        # 复制基因组文件到工作目录 | Copy genome file to working directory
        genome_copy = os.path.join(self.config.output_dir, os.path.basename(genome_file))
        if not os.path.exists(genome_copy):
            import shutil
            shutil.copy2(genome_file, genome_copy)
        
        command = f"BuildDatabase -name {database_name} {os.path.basename(genome_file)}"
        
        success = self.cmd_runner.run(command, "构建RepeatModeler数据库 | Building RepeatModeler database")
        if success:
            self.config.modeler_database = database_name
            self.logger.info(f"RepeatModeler数据库构建完成 | RepeatModeler database built: {database_name}")
        
        return success
    
    def run_repeat_modeler(self) -> bool:
        """运行RepeatModeler | Run RepeatModeler"""
        if not hasattr(self.config, 'modeler_database'):
            self.logger.error("RepeatModeler数据库未构建 | RepeatModeler database not built")
            return False
        
        database_name = self.config.modeler_database
        threads = self.config.rm_model_threads
        
        command_parts = [
            "RepeatModeler",
            f"-database {database_name}",
            f"-pa {threads}"
        ]
        
        if self.config.rm_use_ltr_struct:
            command_parts.append("-LTRStruct")
        
        command = " ".join(command_parts)
        
        success = self.cmd_runner.run(command, "运行RepeatModeler | Running RepeatModeler")
        if success:
            # 查找生成的重复序列库 | Find generated repeat library
            import glob
            lib_files = glob.glob(os.path.join(self.config.output_dir, "RM_*/consensi.fa.classified"))
            if lib_files:
                self.config.custom_repeat_lib = lib_files[0]
                self.logger.info(f"RepeatModeler重复序列库生成 | RepeatModeler library generated: {self.config.custom_repeat_lib}")
        
        return success

class RepeatMaskerRunner:
    """RepeatMasker运行器 | RepeatMasker Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_repeat_masker(self) -> bool:
        """运行RepeatMasker | Run RepeatMasker"""
        genome_file = self.config.genome_file
        species = self.config.species
        threads = self.config.rm_threads
        output_dir = os.path.join(self.config.output_dir, "repeatmasker_output")
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        command_parts = [
            "RepeatMasker",
            f"-species \"{species}\"",
            f"-dir {output_dir}",
            f"-pa {threads}"
        ]
        
        # 添加可选参数 | Add optional parameters
        if self.config.rm_soft_mask:
            command_parts.append("-xsmall")
        
        if self.config.rm_generate_gff:
            command_parts.append("-gff")
        
        if self.config.rm_exclude_low_complexity:
            command_parts.append("-excln")
        
        # 使用自定义库或数据库 | Use custom library or database
        if hasattr(self.config, 'custom_repeat_lib') and self.config.custom_repeat_lib:
            command_parts.append(f"-lib {self.config.custom_repeat_lib}")
        elif self.config.custom_lib:
            command_parts.append(f"-lib {self.config.custom_lib}")
        
        command_parts.append(genome_file)
        command = " ".join(command_parts)
        
        success = self.cmd_runner.run(command, "运行RepeatMasker | Running RepeatMasker")
        
        if success:
            # 记录输出文件 | Record output files
            genome_basename = os.path.basename(genome_file)
            output_files = [
                os.path.join(output_dir, f"{genome_basename}.masked"),
                os.path.join(output_dir, f"{genome_basename}.out"),
                os.path.join(output_dir, f"{genome_basename}.tbl")
            ]
            
            if self.config.rm_generate_gff:
                output_files.append(os.path.join(output_dir, f"{genome_basename}.out.gff"))
            
            self.config.repeatmasker_outputs = [f for f in output_files if os.path.exists(f)]
            self.logger.info("RepeatMasker分析完成 | RepeatMasker analysis completed")
            
            for file in self.config.repeatmasker_outputs:
                self.logger.info(f"  输出文件 | Output file: {file}")
        
        return success

class TRFRunner:
    """TRF (Tandem Repeats Finder) 运行器 | TRF Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_trf(self) -> bool:
        """运行TRF分析 | Run TRF analysis"""
        genome_file = self.config.genome_file
        output_dir = os.path.join(self.config.output_dir, "trf_output")
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 复制基因组文件到TRF输出目录 | Copy genome file to TRF output directory
        import shutil
        genome_copy = os.path.join(output_dir, os.path.basename(genome_file))
        if not os.path.exists(genome_copy):
            shutil.copy2(genome_file, genome_copy)
        
        # 构建TRF命令 | Build TRF command
        command = (
            f"cd {output_dir} && "
            f"trf {os.path.basename(genome_file)} "
            f"{self.config.trf_match_weight} "
            f"{self.config.trf_mismatch_penalty} "
            f"{self.config.trf_indel_penalty} "
            f"{self.config.trf_min_score} "
            f"{self.config.trf_max_period} "
            f"{self.config.trf_min_copies} "
            f"{self.config.trf_max_length} "
            f"-h -d"
        )
        
        success = self.cmd_runner.run(command, "运行TRF串联重复分析 | Running TRF tandem repeat analysis")
        
        if success:
            # 查找TRF输出文件 | Find TRF output files
            import glob
            trf_files = glob.glob(os.path.join(output_dir, "*.dat"))
            if trf_files:
                self.config.trf_outputs = trf_files
                self.logger.info(f"TRF分析完成，生成 {len(trf_files)} 个输出文件 | TRF analysis completed, generated {len(trf_files)} output files")
            else:
                self.logger.warning("未找到TRF输出文件 | No TRF output files found")
        
        return success

class EDTARunner:
    """EDTA (Extensive de-novo TE Annotator) 运行器 | EDTA Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_edta(self) -> bool:
        """运行EDTA转座元件分析 | Run EDTA transposable element analysis"""
        genome_file = self.config.genome_file
        output_dir = os.path.join(self.config.output_dir, "edta_output")
        
        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # 复制基因组文件到EDTA输出目录 | Copy genome file to EDTA output directory
        import shutil
        genome_copy = os.path.join(output_dir, os.path.basename(genome_file))
        if not os.path.exists(genome_copy):
            shutil.copy2(genome_file, genome_copy)
        
        # 构建EDTA命令 | Build EDTA command
        command_parts = [
            f"cd {output_dir} &&",
            "EDTA.pl",
            f"--genome {os.path.basename(genome_file)}",
            f"--species {self.config.edta_species}",
            f"--step {self.config.edta_step}",
            f"--threads {self.config.edta_threads}"
        ]
        
        # 添加可选参数 | Add optional parameters
        if self.config.edta_sensitive:
            command_parts.append("--sensitive 1")
        
        if self.config.edta_anno:
            command_parts.append("--anno 1")
        
        if self.config.edta_evaluate:
            command_parts.append("--evaluate 1")
        
        if self.config.edta_overwrite:
            command_parts.append("--overwrite 1")
        
        command = " ".join(command_parts)
        
        success = self.cmd_runner.run(command, "运行EDTA转座元件分析 | Running EDTA transposable element analysis")
        
        if success:
            # 查找EDTA输出文件 | Find EDTA output files
            import glob
            
            # 主要输出文件 | Main output files
            genome_basename = os.path.splitext(os.path.basename(genome_file))[0]
            expected_files = [
                f"{genome_basename}.EDTA.TEanno.gff3",  # TE注释文件
                f"{genome_basename}.EDTA.TElib.fa",     # TE库文件
                f"{genome_basename}.EDTA.TEanno.sum",   # 统计文件
                f"{genome_basename}.EDTA.intact.fa",    # 完整TE序列
                f"{genome_basename}.EDTA.intact.gff3"   # 完整TE注释
            ]
            
            edta_outputs = []
            for expected_file in expected_files:
                full_path = os.path.join(output_dir, expected_file)
                if os.path.exists(full_path):
                    edta_outputs.append(full_path)
            
            if edta_outputs:
                self.config.edta_outputs = edta_outputs
                self.config.edta_lib = os.path.join(output_dir, f"{genome_basename}.EDTA.TElib.fa")
                self.logger.info(f"EDTA分析完成，生成 {len(edta_outputs)} 个输出文件 | EDTA analysis completed, generated {len(edta_outputs)} output files")
                
                for file in edta_outputs:
                    self.logger.info(f"  EDTA输出文件 | EDTA output file: {os.path.basename(file)}")
            else:
                self.logger.warning("未找到EDTA输出文件 | No EDTA output files found")
        
        return success
