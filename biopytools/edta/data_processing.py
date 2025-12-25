# """
# 🌾 EDTA数据处理模块 | EDTA Data Processing Module
# """

# import os
# import json
# from pathlib import Path
# from typing import Dict, List
# import subprocess
# from .utils import CommandRunner, save_checkpoint

# class EDTAProcessor:
#     """EDTA数据处理器 | EDTA Data Processor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner, directories: Dict[str, Path]):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
#         self.directories = directories
        
#     def prepare_input_files(self) -> Dict[str, str]:
#         """准备输入文件 | Prepare input files"""
#         self.logger.info("📋 准备输入文件 | Preparing input files")
        
#         input_files = {}
        
#         if self.config.batch_mode:
#             # 批量模式处理
#             self.logger.info(f"📊 批量模式：处理 {len(self.config.genome_list)} 个基因组 | Batch mode: processing {len(self.config.genome_list)} genomes")
            
#             for i, genome in enumerate(self.config.genome_list):
#                 genome_name = Path(genome).stem
#                 target_path = self.directories["input"] / f"{genome_name}_{i:03d}.fa"
                
#                 # 复制到工作目录
#                 if not target_path.exists():
#                     import shutil
#                     shutil.copy2(genome, target_path)
#                     self.logger.info(f"📂 复制基因组文件 | Copied genome file: {genome} -> {target_path}")
                
#                 input_files[f"genome_{i:03d}"] = str(target_path)
#         else:
#             # 单个基因组处理
#             genome_name = Path(self.config.genome).stem
#             target_path = self.directories["input"] / f"{genome_name}.fa"
            
#             if not target_path.exists():
#                 import shutil
#                 shutil.copy2(self.config.genome, target_path)
#                 self.logger.info(f"📂 复制基因组文件 | Copied genome file: {self.config.genome} -> {target_path}")
            
#             input_files["genome"] = str(target_path)
        
#         # 保存输入文件信息
#         input_info_file = self.directories["input"] / "input_files_info.json"
#         with open(input_info_file, 'w', encoding='utf-8') as f:
#             json.dump(input_files, f, indent=2, ensure_ascii=False)
        
#         self.logger.info(f"✅ 输入文件准备完成 | Input files preparation completed")
#         return input_files
    
#     def run_edta_analysis(self, input_files: Dict[str, str]) -> Dict[str, Dict]:
#         """运行EDTA分析 | Run EDTA analysis"""
#         self.logger.info("🚀 开始EDTA分析 | Starting EDTA analysis")
        
#         results = {}
        
#         if self.config.batch_mode:
#             # 批量处理模式
#             for genome_key, genome_file in input_files.items():
#                 if genome_key.startswith("genome_"):
#                     self.logger.info(f"📊 处理基因组 | Processing genome: {genome_key}")
#                     result = self._run_single_edta_analysis(genome_file, genome_key)
#                     results[genome_key] = result
                    
#                     # 保存检查点
#                     save_checkpoint(self.config.output_path, f"edta_completed_{genome_key}", 
#                                   {"genome": genome_key, "result": result}, self.logger)
#         else:
#             # 单基因组处理模式
#             genome_file = input_files["genome"]
#             self.logger.info(f"📊 处理单个基因组 | Processing single genome: {genome_file}")
#             result = self._run_single_edta_analysis(genome_file, "single_genome")
#             results["single_genome"] = result
            
#             # 保存检查点
#             save_checkpoint(self.config.output_path, "edta_completed", 
#                           {"result": result}, self.logger)
        
#         self.logger.info("✅ EDTA分析完成 | EDTA analysis completed")
#         return results
    
#     def _run_single_edta_analysis(self, genome_file: str, analysis_name: str) -> Dict:
#         """运行单个EDTA分析 | Run single EDTA analysis"""
#         self.logger.info(f"🔧 运行EDTA分析 | Running EDTA analysis: {analysis_name}")
        
#         # 创建专门的输出目录
#         analysis_dir = self.directories["edta_raw"] / analysis_name
#         analysis_dir.mkdir(parents=True, exist_ok=True)
        
#         # 构建EDTA命令
#         cmd = self.config.to_edta_command(genome_file)
        
#         self.logger.info("🚀 开始EDTA运行，这可能需要较长时间... | Starting EDTA run, this may take a long time...")
        
#         try:
#             # 运行EDTA
#             result = self.cmd_runner.run(
#                 cmd,
#                 description=f"EDTA分析 | EDTA analysis - {analysis_name}",
#                 check=True,
#                 show_progress=True
#             )
            
#             # 整理输出文件
#             output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
            
#             return {
#                 "status": "success",
#                 "genome_file": genome_file,
#                 "output_dir": str(analysis_dir),
#                 "command": " ".join(cmd),
#                 "output_files": output_info,
#                 "stdout": result.stdout if hasattr(result, 'stdout') else ""
#             }
            
#         # except subprocess.CalledProcessError as e:
#         #     self.logger.error(f"❌ EDTA分析失败 | EDTA analysis failed: {analysis_name}")
#         #     return {
#         #         "status": "failed",
#         #         "genome_file": genome_file,
#         #         "output_dir": str(analysis_dir),
#         #         "command": " ".join(cmd),
#         #         "error": str(e),
#         #         "returncode": e.returncode
#         #     }
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"❌ EDTA分析失败 | EDTA analysis failed: {analysis_name}")
            
#             # 如果错误信息包含LTR相关问题，尝试使用--force参数重试
#             error_msg = str(e.stdout) if hasattr(e, 'stdout') else ""
#             if "Raw LTR results not found" in error_msg or "intact LTRs" in error_msg:
#                 self.logger.warning("⚠️ 检测到LTR问题，尝试使用--force参数重试 | LTR issue detected, retrying with --force")
                
#                 # 临时设置force参数
#                 original_force = self.config.force
#                 self.config.force = 1
                
#                 try:
#                     # 重新构建命令并重试
#                     retry_cmd = self.config.to_edta_command(genome_file)
#                     self.logger.info("🔄 使用--force 1重试EDTA分析 | Retrying EDTA analysis with --force 1")
                    
#                     retry_result = self.cmd_runner.run(
#                         retry_cmd,
#                         description=f"EDTA分析重试 | EDTA analysis retry - {analysis_name}",
#                         check=True,
#                         show_progress=True
#                     )
                    
#                     # 恢复原始force设置
#                     self.config.force = original_force
                    
#                     # 整理输出文件
#                     output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
                    
#                     return {
#                         "status": "success",
#                         "genome_file": genome_file,
#                         "output_dir": str(analysis_dir),
#                         "command": " ".join(retry_cmd),
#                         "output_files": output_info,
#                         "stdout": retry_result.stdout if hasattr(retry_result, 'stdout') else "",
#                         "retry_with_force": True
#                     }
                    
#                 except subprocess.CalledProcessError as retry_e:
#                     self.config.force = original_force  # 恢复原始设置
#                     self.logger.error(f"❌ 使用--force参数重试仍然失败 | Retry with --force also failed: {analysis_name}")
            
#             return {
#                 "status": "failed",
#                 "genome_file": genome_file,
#                 "output_dir": str(analysis_dir),
#                 "command": " ".join(cmd),
#                 "error": str(e),
#                 "returncode": e.returncode
#             }
    
#     def _organize_edta_outputs(self, genome_file: str, analysis_dir: Path, analysis_name: str) -> Dict[str, str]:
#         """整理EDTA输出文件 | Organize EDTA output files"""
#         self.logger.info(f"📂 整理EDTA输出文件 | Organizing EDTA output files: {analysis_name}")
        
#         # EDTA的主要输出文件模式
#         genome_basename = Path(genome_file).stem
        
#         expected_outputs = [
#             f"{genome_basename}.fa.mod.EDTA.TEanno.gff3",  # 主要注释文件
#             f"{genome_basename}.fa.mod.EDTA.TElib.fa",     # TE库文件
#             f"{genome_basename}.fa.mod.EDTA.intact.fa",    # 完整TE序列
#             f"{genome_basename}.fa.mod.EDTA.intact.gff3",  # 完整TE注释
#             f"{genome_basename}.fa.mod.MAKER.masked",      # 屏蔽基因组
#             f"{genome_basename}.fa.mod.out",              # RepeatMasker输出
#             f"{genome_basename}.fa.mod.EDTA.anno.sum"      # 统计摘要
#         ]
        
#         output_files = {}
#         source_dir = Path(genome_file).parent  # EDTA通常在输入文件目录生成输出
        
#         for expected_file in expected_outputs:
#             source_file = source_dir / expected_file
#             target_file = analysis_dir / expected_file
            
#             if source_file.exists():
#                 # 移动到指定输出目录
#                 import shutil
#                 shutil.move(str(source_file), str(target_file))
#                 output_files[expected_file] = str(target_file)
#                 self.logger.info(f"📄 移动输出文件 | Moved output file: {expected_file}")
#             else:
#                 self.logger.warning(f"⚠️ 预期输出文件不存在 | Expected output file not found: {expected_file}")
        
#         return output_files

"""
🌾 EDTA数据处理模块 | EDTA Data Processing Module
"""

import os
import json
from pathlib import Path
from typing import Dict, List
import subprocess
from .utils import CommandRunner, save_checkpoint

class EDTAProcessor:
    """EDTA数据处理器 | EDTA Data Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner, directories: Dict[str, Path]):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.directories = directories
        
    def prepare_input_files(self) -> Dict[str, str]:
        """准备输入文件 | Prepare input files"""
        self.logger.info("📋 准备输入文件 | Preparing input files")
        
        input_files = {}
        
        if self.config.batch_mode:
            # 批量模式处理
            self.logger.info(f"📊 批量模式：处理 {len(self.config.genome_list)} 个基因组 | Batch mode: processing {len(self.config.genome_list)} genomes")
            
            for i, genome in enumerate(self.config.genome_list):
                genome_name = Path(genome).stem
                target_path = self.directories["input"] / f"{genome_name}_{i:03d}.fa"
                
                # 复制到工作目录
                if not target_path.exists():
                    import shutil
                    shutil.copy2(genome, target_path)
                    self.logger.info(f"📂 复制基因组文件 | Copied genome file: {genome} -> {target_path}")
                
                input_files[f"genome_{i:03d}"] = str(target_path)
        else:
            # 单个基因组处理
            genome_name = Path(self.config.genome).stem
            target_path = self.directories["input"] / f"{genome_name}.fa"
            
            if not target_path.exists():
                import shutil
                shutil.copy2(self.config.genome, target_path)
                self.logger.info(f"📂 复制基因组文件 | Copied genome file: {self.config.genome} -> {target_path}")
            
            input_files["genome"] = str(target_path)
        
        # 保存输入文件信息
        input_info_file = self.directories["input"] / "input_files_info.json"
        with open(input_info_file, 'w', encoding='utf-8') as f:
            json.dump(input_files, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"✅ 输入文件准备完成 | Input files preparation completed")
        return input_files
    
    def run_edta_analysis(self, input_files: Dict[str, str]) -> Dict[str, Dict]:
        """运行EDTA分析 | Run EDTA analysis"""
        self.logger.info("🚀 开始EDTA分析 | Starting EDTA analysis")
        
        results = {}
        
        if self.config.batch_mode:
            # 批量处理模式
            for genome_key, genome_file in input_files.items():
                if genome_key.startswith("genome_"):
                    self.logger.info(f"📊 处理基因组 | Processing genome: {genome_key}")
                    # 修改为新的方法名
                    result = self._run_single_edta_analysis_with_output_check(genome_file, genome_key)
                    results[genome_key] = result
                    
                    # 保存检查点
                    save_checkpoint(self.config.output_path, f"edta_completed_{genome_key}", 
                                  {"genome": genome_key, "result": result}, self.logger)
        else:
            # 单基因组处理模式
            genome_file = input_files["genome"]
            self.logger.info(f"📊 处理单个基因组 | Processing single genome: {genome_file}")
            # 修改为新的方法名
            result = self._run_single_edta_analysis_with_output_check(genome_file, "single_genome")
            results["single_genome"] = result
            
            # 保存检查点
            save_checkpoint(self.config.output_path, "edta_completed", 
                          {"result": result}, self.logger)
        
        self.logger.info("✅ EDTA分析完成 | EDTA analysis completed")
        return results
    
    def _run_single_edta_analysis_with_output_check(self, genome_file: str, analysis_name: str) -> Dict:
        """运行单个EDTA分析（带输出检查）| Run single EDTA analysis with output checking"""
        self.logger.info(f"🔧 运行EDTA分析 | Running EDTA analysis: {analysis_name}")
        
        # 创建专门的输出目录
        analysis_dir = self.directories["edta_raw"] / analysis_name
        analysis_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建EDTA命令
        cmd = self.config.to_edta_command(genome_file)
        
        self.logger.info("🚀 开始EDTA运行，这可能需要较长时间... | Starting EDTA run, this may take a long time...")
        
        try:
            # 运行EDTA
            result = self.cmd_runner.run(
                cmd,
                description=f"EDTA分析 | EDTA analysis - {analysis_name}",
                check=True,
                show_progress=True
            )
            
            # 整理输出文件
            output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
            
            return {
                "status": "success",
                "genome_file": genome_file,
                "output_dir": str(analysis_dir),
                "command": " ".join(cmd),
                "output_files": output_info,
                "stdout": result.stdout if hasattr(result, 'stdout') else ""
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ EDTA分析报错 | EDTA analysis reported error: {analysis_name}")
            
            # 在声明失败之前，先检查主要输出文件是否存在
            # Before declaring failure, check if main outputs exist
            genome_basename = Path(genome_file).stem
            source_dir = Path(genome_file).parent
            
            # 定义主要输出文件列表
            main_outputs = [
                f"{genome_basename}.fa.mod.EDTA.TEanno.gff3",  # 主要注释文件
                f"{genome_basename}.fa.mod.EDTA.TElib.fa",     # TE库文件  
                f"{genome_basename}.fa.mod.EDTA.TEanno.sum",   # 统计摘要
                f"{genome_basename}.fa.mod.MAKER.masked"       # 屏蔽基因组
            ]
            
            # 检查文件存在性和大小
            existing_outputs = []
            file_status = {}
            
            for output_file in main_outputs:
                file_path = source_dir / output_file
                exists = file_path.exists()
                size = file_path.stat().st_size if exists else 0
                size_mb = round(size / (1024 * 1024), 2) if exists else 0
                
                file_status[output_file] = {
                    'exists': exists,
                    'size_mb': size_mb,
                    'non_empty': size > 0
                }
                
                if exists and size > 0:
                    existing_outputs.append(output_file)
            
            # 记录文件检查结果
            self.logger.info("🔍 检查主要输出文件状态 | Checking main output files status:")
            for filename, status in file_status.items():
                status_icon = "✅" if status['exists'] and status['non_empty'] else "❌"
                size_info = f"({status['size_mb']} MB)" if status['exists'] else "(不存在)"
                self.logger.info(f"  {status_icon} {filename} {size_info}")
            
            # 判断是否应该视为成功
            critical_files = [
                f"{genome_basename}.fa.mod.EDTA.TEanno.gff3",
                f"{genome_basename}.fa.mod.EDTA.TElib.fa"
            ]
            
            critical_exists = sum(1 for f in critical_files if file_status.get(f, {}).get('non_empty', False))
            total_exists = len(existing_outputs)
            
            # 如果至少有2个主要文件存在且非空，视为成功
            if critical_exists >= 2 or total_exists >= 2:
                self.logger.warning("⚠️ EDTA报错但主要输出文件存在 - 视为成功 | EDTA reported errors but main outputs exist - treating as success")
                self.logger.info(f"📊 发现 {total_exists}/{len(main_outputs)} 个主要输出文件 | Found {total_exists}/{len(main_outputs)} main output files")
                
                # 整理输出文件
                try:
                    output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
                    
                    return {
                        "status": "success_with_warnings",
                        "genome_file": genome_file,
                        "output_dir": str(analysis_dir),
                        "command": " ".join(cmd),
                        "output_files": output_info,
                        "warnings": f"Process reported errors but {total_exists} main outputs found",
                        "error_message": str(e),
                        "file_status": file_status,
                        "stdout": getattr(e, 'stdout', '') or str(e)
                    }
                except Exception as organize_error:
                    self.logger.error(f"❌ 整理输出文件时出错 | Error organizing output files: {organize_error}")
                    
                    # 即使整理失败，也返回文件信息
                    manual_output_info = {}
                    for output_file in main_outputs:
                        file_path = source_dir / output_file  
                        if file_path.exists():
                            manual_output_info[output_file] = str(file_path)
                    
                    return {
                        "status": "success_with_warnings",
                        "genome_file": genome_file,
                        "output_dir": str(analysis_dir),
                        "command": " ".join(cmd),
                        "output_files": manual_output_info,
                        "warnings": f"Process errors and file organization issues, but {total_exists} outputs found",
                        "error_message": str(e),
                        "organize_error": str(organize_error),
                        "file_status": file_status
                    }
            
            # 如果关键文件都不存在，尝试重试逻辑
            else:
                self.logger.error(f"❌ 关键输出文件缺失 ({critical_exists}/2 critical files found) | Critical output files missing")
                
                # 检查是否为已知的LTR相关错误
                error_msg = getattr(e, 'stdout', '') or str(e)
                ltr_error_patterns = [
                    "Raw LTR results not found",
                    "intact LTRs", 
                    "No such file or directory",
                    ".cln.cln",
                    "LTR.intact.fa.ori.dusted.cln.cln",
                    "output_by_list.pl"
                ]
                
                is_ltr_error = any(pattern in error_msg for pattern in ltr_error_patterns)
                
                if is_ltr_error:
                    self.logger.warning("⚠️ 检测到LTR相关错误，尝试使用--force重试 | Detected LTR-related error, trying --force retry")
                    return self._retry_with_force(genome_file, analysis_dir, analysis_name, cmd, e)
                else:
                    # 非LTR错误，直接返回失败
                    return {
                        "status": "failed",
                        "genome_file": genome_file,
                        "output_dir": str(analysis_dir),
                        "command": " ".join(cmd),
                        "error": str(e),
                        "error_message": error_msg,
                        "returncode": e.returncode,
                        "file_status": file_status
                    }

    def _retry_with_force(self, genome_file: str, analysis_dir: Path, analysis_name: str, 
                         original_cmd: List[str], original_error: Exception) -> Dict:
        """使用--force参数重试EDTA | Retry EDTA with --force parameter"""
        
        # 临时设置force参数
        original_force = self.config.force
        self.config.force = 1
        
        try:
            # 重新构建命令并重试
            retry_cmd = self.config.to_edta_command(genome_file)
            self.logger.info("🔄 使用--force 1重试EDTA分析 | Retrying EDTA analysis with --force 1")
            
            retry_result = self.cmd_runner.run(
                retry_cmd,
                description=f"EDTA分析重试 | EDTA analysis retry - {analysis_name}",
                check=True,
                show_progress=True
            )
            
            # 恢复原始force设置
            self.config.force = original_force
            
            # 整理输出文件
            output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
            
            return {
                "status": "success", 
                "genome_file": genome_file,
                "output_dir": str(analysis_dir),
                "command": " ".join(retry_cmd),
                "output_files": output_info,
                "stdout": retry_result.stdout if hasattr(retry_result, 'stdout') else "",
                "retry_with_force": True,
                "original_error": str(original_error)
            }
            
        except subprocess.CalledProcessError as retry_e:
            self.config.force = original_force  # 恢复原始设置
            self.logger.error(f"❌ 使用--force参数重试仍然失败 | Retry with --force also failed: {analysis_name}")
            
            # 最后再检查一次输出文件
            genome_basename = Path(genome_file).stem  
            source_dir = Path(genome_file).parent
            final_check_files = [
                f"{genome_basename}.fa.mod.EDTA.TEanno.gff3",
                f"{genome_basename}.fa.mod.EDTA.TElib.fa"
            ]
            
            final_existing = [f for f in final_check_files 
                             if (source_dir / f).exists() and (source_dir / f).stat().st_size > 0]
            
            if len(final_existing) >= 1:
                self.logger.warning("⚠️ 重试失败但检测到部分输出文件，标记为部分成功 | Retry failed but some outputs detected, marking as partial success")
                
                try:
                    output_info = self._organize_edta_outputs(genome_file, analysis_dir, analysis_name)
                    return {
                        "status": "success_with_warnings",
                        "genome_file": genome_file,
                        "output_dir": str(analysis_dir),
                        "command": " ".join(retry_cmd),
                        "output_files": output_info,
                        "warnings": "Both original and retry failed but partial outputs exist",
                        "original_error": str(original_error),
                        "retry_error": str(retry_e),
                        "final_check_files": final_existing
                    }
                except:
                    pass
            
            # 完全失败
            return {
                "status": "failed",
                "genome_file": genome_file,
                "output_dir": str(analysis_dir),
                "command": " ".join(original_cmd),
                "error": str(original_error),
                "retry_error": str(retry_e),
                "returncode": retry_e.returncode
            }
    
    def _organize_edta_outputs(self, genome_file: str, analysis_dir: Path, analysis_name: str) -> Dict[str, str]:
        """整理EDTA输出文件 | Organize EDTA output files"""
        self.logger.info(f"📂 整理EDTA输出文件 | Organizing EDTA output files: {analysis_name}")
        
        # EDTA的主要输出文件模式
        genome_basename = Path(genome_file).stem
        
        expected_outputs = [
            f"{genome_basename}.fa.mod.EDTA.TEanno.gff3",  # 主要注释文件
            f"{genome_basename}.fa.mod.EDTA.TElib.fa",     # TE库文件
            f"{genome_basename}.fa.mod.EDTA.intact.fa",    # 完整TE序列
            f"{genome_basename}.fa.mod.EDTA.intact.gff3",  # 完整TE注释
            f"{genome_basename}.fa.mod.MAKER.masked",      # 屏蔽基因组
            f"{genome_basename}.fa.mod.out",              # RepeatMasker输出
            f"{genome_basename}.fa.mod.EDTA.anno.sum",     # 统计摘要
            f"{genome_basename}.fa.mod.EDTA.TEanno.sum"    # 详细统计摘要
        ]
        
        output_files = {}
        source_dir = Path(genome_file).parent  # EDTA通常在输入文件目录生成输出
        
        for expected_file in expected_outputs:
            source_file = source_dir / expected_file
            target_file = analysis_dir / expected_file
            
            if source_file.exists():
                # 移动到指定输出目录
                import shutil
                try:
                    shutil.move(str(source_file), str(target_file))
                    output_files[expected_file] = str(target_file)
                    self.logger.info(f"📄 移动输出文件 | Moved output file: {expected_file}")
                except Exception as move_error:
                    # 如果移动失败，尝试复制
                    try:
                        shutil.copy2(str(source_file), str(target_file))
                        output_files[expected_file] = str(target_file)
                        self.logger.warning(f"⚠️ 文件移动失败，已复制 | File move failed, copied instead: {expected_file}")
                    except Exception as copy_error:
                        self.logger.error(f"❌ 文件处理失败 | File processing failed: {expected_file}, {copy_error}")
            else:
                self.logger.warning(f"⚠️ 预期输出文件不存在 | Expected output file not found: {expected_file}")
        
        # 如果没有找到主要输出文件，记录现有文件
        if not output_files:
            self.logger.warning("⚠️ 未找到预期输出文件，列出源目录内容 | No expected outputs found, listing source directory")
            try:
                all_files = list(source_dir.glob(f"{genome_basename}*"))
                for file_path in all_files[:20]:  # 只显示前20个文件
                    self.logger.info(f"📄 发现文件 | Found file: {file_path.name}")
            except Exception as e:
                self.logger.error(f"❌ 无法列出源目录 | Cannot list source directory: {e}")
        
        return output_files