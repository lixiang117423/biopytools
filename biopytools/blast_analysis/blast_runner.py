# """
# 🔬 BLAST比对执行模块
# """

# import os
# import stat
# import time
# from datetime import datetime
# from pathlib import Path
# from typing import List, Tuple
# from .utils import CommandRunner, SampleManager

# class BLASTRunner:
#     """🔬 BLAST比对执行器"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.sample_manager = SampleManager(config, logger)
        
#     def run_blast_analysis(self, target_db_path: str) -> List[Tuple[str, str, str]]:
#         """运行BLAST分析"""
        
#         self.sample_manager.load_sample_mapping()
#         input_files = self.sample_manager.get_file_paths()
        
#         if not input_files:
#             error_msg = "❌ 样品映射文件中未找到有效的输入文件"
#             self.logger.error(error_msg)
#             raise FileNotFoundError(error_msg)
        
#         blast_results = []
#         failed_commands = []  # 收集失败的命令
#         total_files = len(input_files)
        
#         self.logger.info(f"🚀 开始处理 {total_files} 个输入文件")
        
#         for idx, input_file in enumerate(input_files, 1):
#             file_name = os.path.basename(input_file)
#             sample_name = self.sample_manager.extract_sample_name(file_name)
            
#             self.logger.info(f"🧬 [{idx}/{total_files}] 正在处理: {file_name} (样品: {sample_name})")
            
#             # 设置输出文件名
#             output_filename = f"{sample_name}_vs_{self.config.blast_type}_results.txt"
#             output_file = self.config.output_path / output_filename
            
#             # 构建BLAST命令
#             blast_cmd = self._build_blast_command(input_file, target_db_path, output_file)
            
#             # 运行BLAST（添加重试机制）
#             max_retries = 1  # 减少重试次数，避免浪费时间
#             success = False
            
#             for attempt in range(max_retries + 1):
#                 if attempt > 0:
#                     self.logger.warning(f"  🔄 第{attempt}次重试: {file_name}")
#                     time.sleep(2)
                
#                 success = self.cmd_runner.run(blast_cmd, f"{self.config.blast_type.upper()}比对")
                
#                 if success:
#                     break
#                 else:
#                     if attempt < max_retries:
#                         self.logger.warning(f"  ⏳ 等待后重试...")
            
#             if success:
#                 if output_file.exists() and output_file.stat().st_size > 0:
#                     hit_count = sum(1 for _ in open(output_file, 'r'))
#                     self.logger.info(f"  ✅ BLAST完成，发现 {hit_count} 个比对结果")
#                 else:
#                     self.logger.info(f"  ℹ️  无比对结果")
#                 blast_results.append((file_name, sample_name, str(output_file)))
#             else:
#                 self.logger.error(f"  ❌ BLAST运行失败: {file_name}")
#                 self.logger.warning(f"  📝 将添加到手动执行脚本中")
                
#                 # 收集失败的命令
#                 failed_commands.append({
#                     'command': blast_cmd,
#                     'file_name': file_name,
#                     'sample_name': sample_name,
#                     'output_file': str(output_file)
#                 })
                
#                 # 创建空的结果文件以保持一致性
#                 output_file.touch()
#                 blast_results.append((file_name, sample_name, str(output_file)))
            
#             # 添加短暂延迟
#             time.sleep(0.5)
        
#         # 生成失败命令的shell脚本
#         if failed_commands:
#             self._generate_failed_commands_script(failed_commands)
        
#         self.logger.info(f"🎉 BLAST分析完成！处理了 {len(input_files)} 个文件")
#         return blast_results

#     def _generate_failed_commands_script(self, failed_commands):
#         """生成失败命令的shell脚本"""
#         script_file = self.config.output_path / "failed_blast_commands.sh"
        
#         with open(script_file, 'w', encoding='utf-8') as f:
#             f.write("#!/bin/bash\n")
#             f.write("# BLAST失败命令手动执行脚本\n")
#             f.write("# Script for manually running failed BLAST commands\n")
#             f.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
#             f.write(f"# 失败文件数量: {len(failed_commands)}\n")
#             f.write("#" + "="*60 + "\n\n")
            
#             f.write("echo \"开始执行失败的BLAST命令...\"\n")
#             f.write("echo \"Starting failed BLAST commands...\"\n\n")
            
#             for i, cmd_info in enumerate(failed_commands, 1):
#                 f.write(f"# 任务 {i}: {cmd_info['file_name']} -> {cmd_info['sample_name']}\n")
#                 f.write(f"echo \"[{i}/{len(failed_commands)}] 处理: {cmd_info['file_name']}\"\n")
#                 f.write(f"{cmd_info['command']}\n")
#                 f.write("if [ $? -eq 0 ]; then\n")
#                 f.write(f"    echo \"  ✅ 成功: {cmd_info['file_name']}\"\n")
#                 f.write("else\n")
#                 f.write(f"    echo \"  ❌ 失败: {cmd_info['file_name']}\"\n")
#                 f.write("fi\n")
#                 f.write("echo \"\"\n\n")
            
#             f.write("echo \"所有失败命令执行完成！\"\n")
#             f.write("echo \"All failed commands completed!\"\n")
        
#         # 设置可执行权限
#         script_file.chmod(script_file.stat().st_mode | stat.S_IEXEC)
        
#         self.logger.warning(f"⚠️  生成失败命令脚本: {script_file}")
#         self.logger.warning(f"💡 请手动执行: bash {script_file}")
#         self.logger.warning(f"🔍 或逐个检查失败的文件格式和内容")
    
    
#     def _build_blast_command(self, query_file: str, db_path: str, output_file: Path) -> str:
#         """构建BLAST命令"""
#         blast_executable = self.config.get_blast_executable()
        
#         # 确保使用绝对路径
#         output_file_abs = output_file.resolve()
        
#         cmd_parts = [
#             blast_executable,
#             f"-query {query_file}",
#             f"-db {db_path}",
#             f"-out {output_file_abs}",
#             f"-outfmt \"{self.config.outfmt}\"",
#             f"-num_threads {self.config.threads}",
#             f"-evalue {self.config.evalue}",
#             f"-max_target_seqs {self.config.max_target_seqs}"
#         ]
        
#         # 添加特定于程序的参数
#         if self.config.blast_type in ['blastn', 'tblastx']:
#             cmd_parts.append(f"-word_size {self.config.word_size}")
        
#         return " ".join(cmd_parts)

"""
🔬 BLAST比对执行模块
"""

import os
import stat
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from .utils import CommandRunner, SampleManager

class BLASTRunner:
    """🔬 BLAST比对执行器"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.sample_manager = SampleManager(config, logger)
        
    def run_blast_analysis(self, target_db_path: str) -> List[Tuple[str, str, str]]:
        """运行BLAST分析"""
        
        self.sample_manager.load_sample_mapping()
        input_files = self.sample_manager.get_file_paths()
        
        if not input_files:
            error_msg = "❌ 样品映射文件中未找到有效的输入文件"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        blast_results = []
        failed_commands = []  # 收集失败的命令
        total_files = len(input_files)
        
        self.logger.info(f"🚀 开始处理 {total_files} 个输入文件")
        
        for idx, input_file in enumerate(input_files, 1):
            file_name = os.path.basename(input_file)
            sample_name = self.sample_manager.extract_sample_name(file_name)
            
            self.logger.info(f"🧬 [{idx}/{total_files}] 正在处理: {file_name} (样品: {sample_name})")
            
            # 设置输出文件名 (*** 修改点 ***)
            output_filename = f"{sample_name}_vs_{self.config.target_base_name}_{self.config.blast_type}_results.txt"
            output_file = self.config.output_path / output_filename
            
            # 构建BLAST命令
            blast_cmd = self._build_blast_command(input_file, target_db_path, output_file)
            
            # 运行BLAST（添加重试机制）
            max_retries = 1  # 减少重试次数，避免浪费时间
            success = False
            
            for attempt in range(max_retries + 1):
                if attempt > 0:
                    self.logger.warning(f"  🔄 第{attempt}次重试: {file_name}")
                    time.sleep(2)
                
                success = self.cmd_runner.run(blast_cmd, f"{self.config.blast_type.upper()}比对")
                
                if success:
                    break
                else:
                    if attempt < max_retries:
                        self.logger.warning(f"  ⏳ 等待后重试...")
            
            if success:
                if output_file.exists() and output_file.stat().st_size > 0:
                    hit_count = sum(1 for _ in open(output_file, 'r'))
                    self.logger.info(f"  ✅ BLAST完成，发现 {hit_count} 个比对结果")
                else:
                    self.logger.info(f"  ℹ️  无比对结果")
                blast_results.append((file_name, sample_name, str(output_file)))
            else:
                self.logger.error(f"  ❌ BLAST运行失败: {file_name}")
                self.logger.warning(f"  📝 将添加到手动执行脚本中")
                
                # 收集失败的命令
                failed_commands.append({
                    'command': blast_cmd,
                    'file_name': file_name,
                    'sample_name': sample_name,
                    'output_file': str(output_file)
                })
                
                # 创建空的结果文件以保持一致性
                output_file.touch()
                blast_results.append((file_name, sample_name, str(output_file)))
            
            # 添加短暂延迟
            time.sleep(0.5)
        
        # 生成失败命令的shell脚本
        if failed_commands:
            self._generate_failed_commands_script(failed_commands)
        
        self.logger.info(f"🎉 BLAST分析完成！处理了 {len(input_files)} 个文件")
        return blast_results

    def _generate_failed_commands_script(self, failed_commands):
        """生成失败命令的shell脚本"""
        script_file = self.config.output_path / "failed_blast_commands.sh"
        
        with open(script_file, 'w', encoding='utf-8') as f:
            f.write("#!/bin/bash\n")
            f.write("# BLAST失败命令手动执行脚本\n")
            f.write("# Script for manually running failed BLAST commands\n")
            f.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 失败文件数量: {len(failed_commands)}\n")
            f.write("#" + "="*60 + "\n\n")
            
            f.write("echo \"开始执行失败的BLAST命令...\"\n")
            f.write("echo \"Starting failed BLAST commands...\"\n\n")
            
            for i, cmd_info in enumerate(failed_commands, 1):
                f.write(f"# 任务 {i}: {cmd_info['file_name']} -> {cmd_info['sample_name']}\n")
                f.write(f"echo \"[{i}/{len(failed_commands)}] 处理: {cmd_info['file_name']}\"\n")
                f.write(f"{cmd_info['command']}\n")
                f.write("if [ $? -eq 0 ]; then\n")
                f.write(f"    echo \"  ✅ 成功: {cmd_info['file_name']}\"\n")
                f.write("else\n")
                f.write(f"    echo \"  ❌ 失败: {cmd_info['file_name']}\"\n")
                f.write("fi\n")
                f.write("echo \"\"\n\n")
            
            f.write("echo \"所有失败命令执行完成！\"\n")
            f.write("echo \"All failed commands completed!\"\n")
        
        # 设置可执行权限
        script_file.chmod(script_file.stat().st_mode | stat.S_IEXEC)
        
        self.logger.warning(f"⚠️  生成失败命令脚本: {script_file}")
        self.logger.warning(f"💡 请手动执行: bash {script_file}")
        self.logger.warning(f"🔍 或逐个检查失败的文件格式和内容")
    
    
    def _build_blast_command(self, query_file: str, db_path: str, output_file: Path) -> str:
        """构建BLAST命令"""
        blast_executable = self.config.get_blast_executable()
        
        # 确保使用绝对路径
        output_file_abs = output_file.resolve()
        
        cmd_parts = [
            blast_executable,
            f"-query {query_file}",
            f"-db {db_path}",
            f"-out {output_file_abs}",
            f"-outfmt \"{self.config.outfmt}\"",
            f"-num_threads {self.config.threads}",
            f"-evalue {self.config.evalue}",
            f"-max_target_seqs {self.config.max_target_seqs}"
        ]
        
        # 添加特定于程序的参数
        if self.config.blast_type in ['blastn', 'tblastx']:
            cmd_parts.append(f"-word_size {self.config.word_size}")
        
        return " ".join(cmd_parts)