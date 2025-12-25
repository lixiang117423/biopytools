"""
🔍 编码区预测模块 | Coding Region Prediction Module
"""

from pathlib import Path
from .utils import CommandRunner

class TransDecoderPredictor:
    """🔍 TransDecoder编码区预测器 | TransDecoder Coding Region Predictor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def identify_coding_regions(self):
        """🔬 识别候选编码区域 | Identify candidate coding regions"""
        if not hasattr(self.config, 'trinity_assembly'):
            self.logger.error("❌ 未找到Trinity组装文件，请先执行组装步骤 | No Trinity assembly found, please run assembly step first")
            return False
        
        output_dir = self.config.output_path / "transdecoder_prediction"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.info("🔬 使用TransDecoder识别编码区域 | Using TransDecoder to identify coding regions")
        
        # 复制Trinity组装文件到TransDecoder工作目录 | Copy Trinity assembly to TransDecoder working directory
        trinity_copy = output_dir / "transcripts.fasta"
        if not self._copy_trinity_assembly(trinity_copy):
            return False
        
        # Step 1: 提取长ORF | Extract long ORFs
        if not self._extract_long_orfs(trinity_copy, output_dir):
            return False
        
        # Step 2: 预测编码区域 | Predict coding regions
        if not self._predict_coding_regions(trinity_copy, output_dir):
            return False
        
        self.logger.info("🎉 TransDecoder编码区预测完成 | TransDecoder coding region prediction completed")
        return True
    
    def _copy_trinity_assembly(self, target_file: Path) -> bool:
        """📋 复制Trinity组装文件 | Copy Trinity assembly file"""
        try:
            import shutil
            shutil.copy2(self.config.trinity_assembly, target_file)
            return True
        except Exception as e:
            self.logger.error(f"❌ 复制Trinity组装文件失败 | Failed to copy Trinity assembly: {e}")
            return False
    
    def _extract_long_orfs(self, transcripts_file: Path, output_dir: Path) -> bool:
        """🧬 提取长ORF | Extract long ORFs"""
        self.logger.info("🧬 提取长ORF | Extracting long ORFs")
        
        # 构建TransDecoder.LongOrfs命令 | Build TransDecoder.LongOrfs command
        cmd = (
            f"{self.config.transdecoder_longorfs_path} "
            f"-t {transcripts_file} "
            f"-m {self.config.transdecoder_min_protein_len} "
            f"--genetic_code {self.config.transdecoder_genetic_code}"
        )
        
        # 添加可选参数 | Add optional parameters
        if self.config.transdecoder_complete_orfs_only:
            cmd += " --complete_orfs_only"
        
        # 切换到输出目录 | Change to output directory
        original_working_dir = self.cmd_runner.working_dir
        self.cmd_runner.working_dir = output_dir
        
        success = self.cmd_runner.run(cmd, "🧬 TransDecoder长ORF提取 | TransDecoder long ORF extraction")
        
        # 恢复工作目录 | Restore working directory
        self.cmd_runner.working_dir = original_working_dir
        
        if success:
            # 检查输出文件 | Check output files
            transdecoder_dir = output_dir / f"{transcripts_file.name}.transdecoder_dir"
            longest_orfs = transdecoder_dir / "longest_orfs.pep"
            
            if longest_orfs.exists():
                self.config.longest_orfs = str(longest_orfs)
                self.logger.info(f"✅ 长ORF文件生成 | Long ORFs file generated: {longest_orfs}")
            else:
                self.logger.error("❌ 未找到长ORF输出文件 | Long ORF output file not found")
                return False
        
        return success
    
    def _predict_coding_regions(self, transcripts_file: Path, output_dir: Path) -> bool:
        """🎯 预测编码区域 | Predict coding regions"""
        self.logger.info("🎯 预测编码区域 | Predicting coding regions")
        
        # 构建TransDecoder.Predict命令 | Build TransDecoder.Predict command
        cmd = (
            f"{self.config.transdecoder_predict_path} "
            f"-t {transcripts_file} "
            f"--genetic_code {self.config.transdecoder_genetic_code}"
        )
        
        # 添加可选参数 | Add optional parameters
        if self.config.transdecoder_complete_orfs_only:
            cmd += " --single_best_only"
        
        # 切换到输出目录 | Change to output directory
        original_working_dir = self.cmd_runner.working_dir
        self.cmd_runner.working_dir = output_dir
        
        success = self.cmd_runner.run(cmd, "🎯 TransDecoder编码区预测 | TransDecoder coding region prediction")
        
        # 恢复工作目录 | Restore working directory
        self.cmd_runner.working_dir = original_working_dir
        
        if success:
            # 检查输出文件 | Check output files
            base_name = transcripts_file.name
            output_files = {
                'cds': output_dir / f"{base_name}.transdecoder.cds",
                'proteins': output_dir / f"{base_name}.transdecoder.pep",
                'gff': output_dir / f"{base_name}.transdecoder.gff3",
                'bed': output_dir / f"{base_name}.transdecoder.bed"
            }
            
            self.config.transdecoder_outputs = {}
            for key, file_path in output_files.items():
                if file_path.exists():
                    self.config.transdecoder_outputs[key] = str(file_path)
                    self.logger.info(f"📄 TransDecoder输出 | TransDecoder output ({key}): {file_path}")
                else:
                    self.logger.warning(f"⚠️ TransDecoder输出文件不存在 | TransDecoder output file not found ({key}): {file_path}")
        
        return success
