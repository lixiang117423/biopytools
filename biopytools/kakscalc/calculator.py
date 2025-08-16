"""
🧮 Ka/Ks Calculator计算引擎模块
功能: 调用KaKs_Calculator2.0进行Ka/Ks计算 | Ka/Ks calculation engine using KaKs_Calculator2.0

重要说明 | Important Notes:
- KaKs_Calculator2.0 需要AXT格式输入文件，不是FASTA格式 | Requires AXT format input, not FASTA
- AXT格式: 序号 seq1_name seq2_name 开始位置1 长度1 方向1 开始位置2 长度2 方向2
             sequence1
             sequence2
             空行
- 方法名称必须使用官方名称，如GMYN而不是gamma-MYN | Method names must use official names
"""

import os
import subprocess
import tempfile
from typing import Dict, List, Tuple
from .config import KaKsConfig
from .logger import Logger

class KaKsCalculator:
    """🧮 Ka/Ks计算引擎 | Ka/Ks calculation engine"""
    
    def __init__(self, logger: Logger, kaks_path: str = "KaKs_Calculator"):
        """
        🏗️ 初始化计算引擎 | Initialize calculation engine
        
        Args:
            logger: 日志器实例 | Logger instance
            kaks_path: KaKs_Calculator可执行文件路径 | Path to KaKs_Calculator executable
        """
        self.logger = logger
        self.kaks_path = kaks_path
        self.config = KaKsConfig()
        self._check_kaks_installation()
    
    def _check_kaks_installation(self):
        """🔍 检查KaKs_Calculator是否可用 | Check if KaKs_Calculator is available"""
        try:
            self.logger.info("检查KaKs_Calculator安装 | Checking KaKs_Calculator installation", "🔍")
            result = subprocess.run([self.kaks_path, "-h"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                raise FileNotFoundError()
            self.logger.success("KaKs_Calculator2.0 可用 | KaKs_Calculator2.0 is available")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.logger.error("未找到KaKs_Calculator2.0 | KaKs_Calculator2.0 not found")
            self.logger.info("请安装KaKs_Calculator2.0或使用--kaks-path指定路径 | Please install KaKs_Calculator2.0 or specify path with --kaks-path")
            raise RuntimeError("KaKs_Calculator2.0 not available")
    
    def prepare_input_file(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str], 
                          pairs: List[Tuple[str, str, str]], temp_dir: str) -> str:
        """
        📝 准备KaKs_Calculator输入文件(AXT格式) | Prepare input file for KaKs_Calculator (AXT format)
        
        Args:
            seq1_dict: 第一个FASTA的序列字典 | Sequences from first FASTA
            seq2_dict: 第二个FASTA的序列字典 | Sequences from second FASTA
            pairs: 配对列表 | Pairs list
            temp_dir: 临时目录 | Temporary directory
            
        Returns:
            输入文件路径 | Input file path
        """
        input_file = os.path.join(temp_dir, self.config.OUTPUT_FILES['temp_axt'])
        
        try:
            self.logger.info(f"准备AXT格式输入文件 | Preparing AXT format input file with {len(pairs)} pairs", "📝")
            
            with open(input_file, 'w') as f:
                for idx, (seq1_id, seq2_id, pair_name) in enumerate(pairs):
                    seq1 = seq1_dict[seq1_id]
                    seq2 = seq2_dict[seq2_id]
                    
                    # 确保序列长度一致 | Ensure sequences have same length
                    if len(seq1) != len(seq2):
                        self.logger.warning(f"序列长度不一致 | Sequence length mismatch: {pair_name} ({len(seq1)} vs {len(seq2)})")
                        # 截取到较短长度 | Truncate to shorter length
                        min_len = min(len(seq1), len(seq2))
                        # 确保是3的倍数 | Ensure multiple of 3
                        min_len = (min_len // 3) * 3
                        seq1 = seq1[:min_len]
                        seq2 = seq2[:min_len]
                    
                    # 写入AXT格式 | Write in AXT format
                    # 格式: 序号 seq1_name seq2_name 开始位置1 长度1 方向1 开始位置2 长度2 方向2
                    f.write(f"{idx} {seq1_id} {seq2_id} 1 {len(seq1)} + 1 {len(seq2)} +\n")
                    f.write(f"{seq1}\n")
                    f.write(f"{seq2}\n")
                    f.write("\n")  # 空行分隔 | Empty line separator
            
            # 📊 验证文件大小 | Verify file size
            file_size = os.path.getsize(input_file)
            self.logger.debug(f"AXT输入文件大小: {file_size} bytes | AXT input file size: {file_size} bytes")
            
            self.logger.success(f"AXT格式输入文件准备完成 | AXT format input file prepared: {input_file}")
            return input_file
            
        except Exception as e:
            self.logger.error(f"准备AXT输入文件失败 | Failed to prepare AXT input file: {e}")
            raise
    
    def _verify_axt_file(self, axt_file: str) -> bool:
        """
        🔍 验证AXT文件格式 | Verify AXT file format
        
        Args:
            axt_file: AXT文件路径 | AXT file path
            
        Returns:
            是否有效 | Whether valid
        """
        try:
            with open(axt_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 4:  # 至少需要4行：头部+seq1+seq2+空行
                self.logger.warning("AXT文件行数不足 | AXT file has insufficient lines")
                return False
            
            # 检查第一个条目的格式 | Check first entry format
            first_line = lines[0].strip().split()
            if len(first_line) != 9:
                self.logger.warning(f"AXT头部格式不正确 | Invalid AXT header format: {lines[0]}")
                return False
            
            self.logger.debug(f"AXT文件验证通过 | AXT file validation passed: {len(lines)} lines")
            return True
            
        except Exception as e:
            self.logger.warning(f"AXT文件验证失败 | AXT file validation failed: {e}")
            return False
    
    def run_calculation(self, input_file: str, method: str, temp_dir: str) -> str:
        """
        🚀 运行Ka/Ks计算 | Run Ka/Ks calculation
        
        Args:
            input_file: 输入文件路径 | Input file path
            method: 计算方法 | Calculation method
            temp_dir: 临时目录 | Temporary directory
            
        Returns:
            输出文件路径 | Output file path
        """
        output_file = os.path.join(temp_dir, "kaks_output.txt")
        
        try:
            # 🔍 验证AXT文件格式 | Verify AXT file format
            if not self._verify_axt_file(input_file):
                raise ValueError(f"无效的AXT文件格式 | Invalid AXT file format: {input_file}")
            
            # 🏃‍♂️ 构建命令 | Build command
            cmd = [
                self.kaks_path,
                "-i", input_file,
                "-o", output_file,
                "-m", method
            ]
            
            self.logger.info(f"运行计算 | Running calculation: {' '.join(cmd)}", "🚀")
            self.logger.info(f"使用方法 | Using method: {self.config.get_method_description(method)}", "🧮")
            
            # 🔄 执行计算 | Execute calculation
            result = subprocess.run(cmd, capture_output=True, text=True, 
                                  timeout=self.config.CALCULATION_TIMEOUT)
            
            # 📝 记录详细输出信息 | Log detailed output information
            self.logger.debug(f"KaKs_Calculator返回码 | Return code: {result.returncode}")
            if result.stdout:
                self.logger.debug(f"KaKs_Calculator标准输出 | Stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"KaKs_Calculator错误输出 | Stderr: {result.stderr}")
            
            if result.returncode != 0:
                error_msg = f"KaKs_Calculator执行失败 | KaKs_Calculator failed (code {result.returncode})"
                if result.stderr:
                    error_msg += f": {result.stderr}"
                self.logger.error(error_msg)
                raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
            
            # ✅ 检查输出文件 | Check output file
            if not os.path.exists(output_file):
                raise FileNotFoundError(f"输出文件未生成 | Output file not created: {output_file}")
            
            # 📊 检查文件内容 | Check file content
            with open(output_file, 'r') as f:
                lines = f.readlines()
                if len(lines) <= 1:  # 只有标题行 | Only header line
                    raise ValueError("输出文件为空或无有效数据 | Output file is empty or has no valid data")
            
            self.logger.success(f"计算完成 | Calculation completed: {output_file}")
            
            # 🔍 记录输出信息 | Log output info
            if result.stdout:
                self.logger.debug(f"KaKs_Calculator输出 | KaKs_Calculator output: {result.stdout}")
            
            return output_file
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"计算超时 ({self.config.CALCULATION_TIMEOUT}秒) | Calculation timed out ({self.config.CALCULATION_TIMEOUT}s)")
            raise
        except Exception as e:
            self.logger.error(f"Ka/Ks计算失败 | Ka/Ks calculation failed: {e}")
            raise
