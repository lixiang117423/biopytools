"""
📤 输出格式化模块 | Output Formatting Module
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
from .utils import CommandRunner

class KmerOutputFormatter:
    """📝 K-mer输出格式化器 | K-mer Output Formatter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_unik_to_fasta(self, unik_file: str, output_fasta: str, is_fastq_input: bool = False, 
                             fastq_samples: Dict[str, Tuple[str, str]] = None, 
                             fasta_files: List[str] = None) -> bool:
        """将unik文件转换为FASTA格式 | Convert unik file to FASTA format"""
        self.logger.info(f"📝 转换为FASTA格式 | Converting to FASTA format: {output_fasta}")
        
        try:
            if is_fastq_input:
                return self._convert_fastq_to_fasta(unik_file, output_fasta, fastq_samples)
            else:
                return self._convert_fasta_to_fasta(unik_file, output_fasta, fasta_files)
        except Exception as e:
            self.logger.error(f"❌ 转换FASTA失败 | FASTA conversion failed: {e}")
            return False
    
    def convert_jellyfish_to_fasta(self, jf_file: str, output_fasta: str) -> bool:
        """将Jellyfish二进制文件转换为FASTA格式 | Convert Jellyfish binary file to FASTA format"""
        self.logger.info(f"🐟 将Jellyfish文件转换为FASTA格式 | Converting Jellyfish file to FASTA format: {output_fasta}")
        
        # 使用jellyfish dump命令转换
        dump_cmd = f"{self.config.jellyfish_path} dump {jf_file}"
        
        try:
            result = subprocess.run(dump_cmd, shell=True, capture_output=True, text=True, cwd=self.config.output_dir)
            
            if result.returncode != 0:
                self.logger.error(f"❌ jellyfish dump失败 | jellyfish dump failed: {result.stderr}")
                return False
            
            # 解析dump输出并转换为FASTA格式
            lines = result.stdout.strip().split('\n')
            kmers = []

            for line in lines:
                line = line.strip()
                if line and not line.startswith('>'):  # 跳过以>开头的计数行，只保留序列行
                    kmers.append(line)
            
            # 写入FASTA文件
            with open(output_fasta, 'w') as f:
                for i, kmer in enumerate(kmers, 1):
                    f.write(f">kmer_{i}\n{kmer}\n")
            
            self.logger.info(f"✅ 成功生成FASTA文件 | Successfully generated FASTA file: {len(kmers)} k-mers")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ Jellyfish转换FASTA失败 | Jellyfish FASTA conversion failed: {e}")
            return False
    
    def _convert_fastq_to_fasta(self, unik_file: str, output_fasta: str, fastq_samples: Dict[str, Tuple[str, str]]) -> bool:
        """转换FASTQ输入的k-mer为FASTA | Convert FASTQ input k-mers to FASTA"""
        view_cmd = f"{self.config.unikmer_path} view -t {unik_file}"
        
        result = subprocess.run(view_cmd, shell=True, capture_output=True, text=True, cwd=self.config.output_dir)
        
        if result.returncode != 0:
            self.logger.error(f"❌ unikmer view失败 | unikmer view failed: {result.stderr}")
            return False
        
        lines = result.stdout.strip().split('\n')
        kmers = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
        
        with open(output_fasta, 'w') as f:
            for i, kmer in enumerate(kmers, 1):
                f.write(f">kmer_{i}\n{kmer}\n")
        
        self.logger.info(f"✅ 成功生成FASTA文件 | Successfully generated FASTA file: {len(kmers)} k-mers")
        return True
    
    def _convert_fasta_to_fasta(self, unik_file: str, output_fasta: str, fasta_files: List[str]) -> bool:
        """转换FASTA输入的k-mer为FASTA（带位置信息）| Convert FASTA input k-mers to FASTA (with position info)"""
        # 生成临时BED文件来获取位置信息
        temp_bed = os.path.join(self.config.output_dir, "temp_locate_for_fasta")
        locate_cmd = f"{self.config.unikmer_path} locate {unik_file} -g {' '.join(fasta_files)} -o {temp_bed}"
        
        if not self.cmd_runner.run(locate_cmd, "📍 获取k-mer位置信息用于FASTA生成 | Getting k-mer position info for FASTA generation"):
            self.logger.warning("⚠️ 无法获取位置信息，使用简化格式 | Cannot get position info, using simplified format")
            return self._convert_fastq_to_fasta(unik_file, output_fasta, None)
        
        # 找到实际生成的BED文件
        bed_file = None
        possible_files = [temp_bed, f"{temp_bed}.bed", f"{temp_bed}.txt", f"{temp_bed}.out"]
        
        for possible_file in possible_files:
            if os.path.exists(possible_file):
                bed_file = possible_file
                break
        
        if not bed_file:
            self.logger.warning("⚠️ 找不到unikmer locate输出文件 | Cannot find unikmer locate output file")
            return self._convert_fastq_to_fasta(unik_file, output_fasta, None)
        
        # 读取BED文件并生成FASTA
        fasta_entries = []
        try:
            with open(bed_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            seq_name = parts[0]
                            start_pos = int(parts[1]) + 1  # BED格式是0-based，转为1-based
                            end_pos = int(parts[2])
                            kmer_seq = parts[3]
                            
                            # 生成FASTA ID：>OV12_1_51
                            fasta_id = f">{seq_name}_{start_pos}_{end_pos}"
                            fasta_entries.append((seq_name, start_pos, fasta_id, kmer_seq))
        
            # 按序列名称和位置排序
            fasta_entries.sort(key=lambda x: (x[0], x[1]))
            
            # 写入FASTA文件
            with open(output_fasta, 'w') as f:
                for _, _, fasta_id, kmer_seq in fasta_entries:
                    f.write(f"{fasta_id}\n{kmer_seq}\n")
            
            # 清理临时文件
            if os.path.exists(bed_file):
                os.remove(bed_file)
            
            self.logger.info(f"✅ 成功生成FASTA文件 | Successfully generated FASTA file: {len(fasta_entries)} k-mers")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 解析BED输出失败 | Failed to parse BED output: {e}")
            if bed_file and os.path.exists(bed_file):
                os.remove(bed_file)
            return self._convert_fastq_to_fasta(unik_file, output_fasta, None)
    
    def generate_bed_file(self, unik_file: str, output_bed: str, fasta_files: List[str]) -> bool:
        """生成BED文件（仅用于FASTA输入）| Generate BED file (for FASTA input only)"""
        self.logger.info(f"📋 生成BED文件 | Generating BED file: {output_bed}")
        
        try:
            bed_output = output_bed if output_bed.endswith('.bed') else f"{output_bed}.bed"
            
            base_name = self.config.base_name
            temp_output = os.path.join(self.config.output_dir, base_name + "_locate")
            
            locate_cmd = f"{self.config.unikmer_path} locate {unik_file} -g {' '.join(fasta_files)} -o {temp_output}"
            
            if self.cmd_runner.run(locate_cmd, "📋 生成BED位置信息 | Generating BED position information"):
                # 检查所有可能的输出文件
                possible_files = [
                    temp_output,
                    os.path.join(self.config.output_dir, base_name),
                    f"{temp_output}.bed",
                    f"{temp_output}.txt"
                ]
                
                found_file = None
                for possible_file in possible_files:
                    if os.path.exists(possible_file):
                        found_file = possible_file
                        self.logger.info(f"📍 找到unikmer locate输出文件 | Found unikmer locate output file: {found_file}")
                        break
                
                if found_file:
                    # 对文件进行排序
                    self._sort_bed_file(found_file)
                    
                    # 如果文件名不是目标名称，则重命名
                    if found_file != bed_output:
                        os.rename(found_file, bed_output)
                        self.logger.info(f"📝 文件已重命名 | File renamed: {os.path.basename(found_file)} -> {os.path.basename(bed_output)}")
                    
                    self.logger.info(f"✅ 成功生成并排序BED文件 | Successfully generated and sorted BED file: {bed_output}")
                    return True
                
                self.logger.warning(f"⚠️ 未找到有效的unikmer locate输出，将创建简化版本 | No valid unikmer locate output found")
                return self._create_simple_bed(unik_file, bed_output)
            else:
                self.logger.warning(f"⚠️ unikmer locate命令失败 | unikmer locate command failed")
                return self._create_simple_bed(unik_file, bed_output)
                
        except Exception as e:
            self.logger.warning(f"⚠️ 生成BED文件出错 | Error generating BED file: {e}")
            return self._create_simple_bed(unik_file, bed_output)
    
    def _sort_bed_file(self, bed_file: str) -> bool:
        """对BED文件按位置排序 | Sort BED file by position"""
        try:
            if not os.path.exists(bed_file):
                return False
                
            with open(bed_file, 'r') as f:
                lines = f.readlines()
            
            bed_entries = []
            for line in lines:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            seq_name = parts[0]
                            start_pos = int(parts[1])
                            bed_entries.append((seq_name, start_pos, line.strip()))
                        except ValueError:
                            bed_entries.append(('', float('inf'), line.strip()))
            
            # 按序列名称和起始位置排序
            bed_entries.sort(key=lambda x: (x[0], x[1]))
            
            # 写回文件
            with open(bed_file, 'w') as f:
                for _, _, line in bed_entries:
                    f.write(f"{line}\n")
            
            self.logger.info(f"✅ BED文件排序完成 | BED file sorted: {len(bed_entries)} entries")
            return True
            
        except Exception as e:
            self.logger.warning(f"⚠️ BED文件排序失败 | BED file sorting failed: {e}")
            return False
    
    def _create_simple_bed(self, unik_file: str, output_bed: str) -> bool:
        """创建简化的BED文件 | Create simplified BED file"""
        try:
            bed_output = output_bed if output_bed.endswith('.bed') else f"{output_bed}.bed"
            
            view_cmd = f"{self.config.unikmer_path} view -t {unik_file}"
            result = subprocess.run(view_cmd, shell=True, capture_output=True, text=True, cwd=self.config.output_dir)
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                kmers = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
                
                with open(bed_output, 'w') as f:
                    for i, kmer in enumerate(kmers, 1):
                        f.write(f"seq\t{i}\t{i+self.config.kmer_length-1}\t{kmer}\n")
                
                self.logger.info(f"✅ 创建简化BED文件成功 | Successfully created simplified BED file")
                return True
            else:
                self.logger.error(f"❌ 无法获取k-mer列表 | Cannot get k-mer list")
                return False
                
        except Exception as e:
            self.logger.error(f"❌ 创建简化BED文件失败 | Failed to create simplified BED file: {e}")
            return False
