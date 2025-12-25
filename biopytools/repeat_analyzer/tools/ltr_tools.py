"""
LTR分析工具包装器 | LTR Analysis Tool Wrappers
"""

import os
from pathlib import Path

class LTRFinderRunner:
    """🔍 LTR_FINDER运行器 | LTR_FINDER Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_file = None
    
    def run_ltr_finder(self) -> bool:
        """🚀 运行LTR_FINDER | Run LTR_FINDER"""
        output_file = self.config.output_path / f"{self.config.base_name}_ltrfinder.scn"
        
        cmd = f"{self.config.ltr_finder_path} -C -w 2 {self.config.genome_file} > {output_file}"
        
        success = self.cmd_runner.run(
            cmd,
            "🔍 运行LTR_FINDER识别长末端重复序列 | Run LTR_FINDER to identify long terminal repeats"
        )
        
        if success and output_file.exists():
            self.output_file = output_file
            self.logger.info(f"✅ LTR_FINDER输出文件生成 | LTR_FINDER output file generated: {output_file}")
        
        return success
    
    def get_output_path(self) -> Path:
        """📁 获取输出文件路径 | Get output file path"""
        return self.output_file

class LTRHarvestRunner:
    """🌾 LTRharvest运行器 | LTRharvest Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_file = None
    
    # def run_ltrharvest(self) -> bool:
    #     """🚀 运行LTRharvest | Run LTRharvest"""
    #     output_file = self.config.output_path / f"{self.config.base_name}_ltrharvest.out"
        
    #     cmd = f"{self.config.ltrharvest_path} -out {output_file} -outinner -gff3 {output_file}.gff3 {self.config.genome_file}"
        
    #     success = self.cmd_runner.run(
    #         cmd,
    #         "🌾 运行LTRharvest识别LTR逆转录转座子 | Run LTRharvest to identify LTR retrotransposons"
    #     )
        
    #     if success and output_file.exists():
    #         self.output_file = output_file
    #         self.logger.info(f"✅ LTRharvest输出文件生成 | LTRharvest output file generated: {output_file}")
        
    #     return success
    def run_ltrharvest(self) -> bool:
        """🚀 运行LTRharvest | Run LTRharvest"""
        output_file = self.config.output_path / f"{self.config.base_name}_ltrharvest.out"
        index_name = self.config.output_path / f"{self.config.base_name}_index"
        
        # 第一步：创建索引
        index_cmd = f"gt suffixerator -db {self.config.genome_file} -indexname {index_name} -tis -suf -lcp -des -ssp -sds -dna"
        
        if not self.cmd_runner.run(index_cmd, "创建基因组索引 | Create genome index"):
            return False
        
        # 第二步：运行LTRharvest
        cmd = f"{self.config.ltrharvest_path} -index {index_name} -out {output_file} -outinner {output_file}.inner -gff3 {output_file}.gff3"
        
        success = self.cmd_runner.run(
            cmd,
            "🌾 运行LTRharvest识别LTR逆转录转座子 | Run LTRharvest to identify LTR retrotransposons"
        )
        
        if success and output_file.exists():
            self.output_file = output_file
            self.logger.info(f"✅ LTRharvest输出文件生成 | LTRharvest output file generated: {output_file}")
        
        return success
        
    def get_output_path(self) -> Path:
        """📁 获取输出文件路径 | Get output file path"""
        return self.output_file

class LTRRetrieverRunner:
    """🔄 LTR_retriever运行器 | LTR_retriever Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_lib = None
    
    def run_ltr_retriever(self, ltr_finder_out: Path, ltrharvest_out: Path) -> bool:
        """🚀 运行LTR_retriever | Run LTR_retriever"""
        # 准备输入文件列表 | Prepare input file list
        input_files = []
        if ltr_finder_out and ltr_finder_out.exists():
            input_files.append(str(ltr_finder_out))
        if ltrharvest_out and ltrharvest_out.exists():
            input_files.append(str(ltrharvest_out))
        
        if not input_files:
            self.logger.error("❌ 没有可用的LTR输入文件 | No available LTR input files")
            return False
        
        # LTR_retriever输出前缀 | LTR_retriever output prefix
        output_prefix = self.config.output_path / self.config.base_name
        
        cmd = f"{self.config.ltr_retriever_path} -genome {self.config.genome_file} -inharvest {' '.join(input_files)} -threads {self.config.threads}"
        
        success = self.cmd_runner.run(
            cmd,
            "🔄 运行LTR_retriever整合LTR分析结果 | Run LTR_retriever to integrate LTR analysis results"
        )
        
        if success:
            # 查找生成的库文件 | Find generated library file
            expected_lib = Path(f"{self.config.genome_file}.mod.EDTA.intact.fa")
            alt_lib = self.config.output_path / f"{self.config.base_name}.LTRlib.fa"
            
            if expected_lib.exists():
                self.output_lib = expected_lib
            elif alt_lib.exists():
                self.output_lib = alt_lib
            else:
                # 查找其他可能的输出文件 | Find other possible output files
                possible_libs = list(self.config.output_path.glob("*.LTRlib.fa"))
                if possible_libs:
                    self.output_lib = possible_libs[0]
                    self.logger.warning(f"⚠️ 使用找到的LTR库文件 | Using found LTR library file: {self.output_lib}")
                else:
                    self.logger.warning("⚠️ 未找到LTR_retriever库文件 | LTR_retriever library file not found")
        
        return success
    
    def get_library_path(self) -> Path:
        """📁 获取生成的库文件路径 | Get generated library file path"""
        return self.output_lib
