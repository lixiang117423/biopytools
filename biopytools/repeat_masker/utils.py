"""
重复序列分析工具函数模块 | Repeat Sequence Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, List

class RepeatLogger:
    """重复序列分析日志管理器 | Repeat Sequence Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "repeat_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False

class SequenceValidator:
    """序列文件验证器 | Sequence File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def validate_fasta(self, fasta_file: str) -> bool:
        """验证FASTA文件格式 | Validate FASTA file format"""
        try:
            with open(fasta_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"无效的FASTA格式 | Invalid FASTA format: {fasta_file}")
                    return False
            
            self.logger.info(f"FASTA文件格式验证通过 | FASTA file format validation passed: {fasta_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"读取FASTA文件时出错 | Error reading FASTA file: {e}")
            return False
    
    def get_sequence_stats(self, fasta_file: str) -> Dict:
        """获取序列统计信息 | Get sequence statistics"""
        stats = {
            'num_sequences': 0,
            'total_length': 0,
            'gc_content': 0,
            'n_content': 0
        }
        
        try:
            with open(fasta_file, 'r') as f:
                sequence = ""
                gc_count = 0
                n_count = 0
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if sequence:
                            stats['total_length'] += len(sequence)
                            gc_count += sequence.upper().count('G') + sequence.upper().count('C')
                            n_count += sequence.upper().count('N')
                            sequence = ""
                        stats['num_sequences'] += 1
                    else:
                        sequence += line
                
                # 处理最后一个序列 | Process the last sequence
                if sequence:
                    stats['total_length'] += len(sequence)
                    gc_count += sequence.upper().count('G') + sequence.upper().count('C')
                    n_count += sequence.upper().count('N')
                
                # 计算百分比 | Calculate percentages
                if stats['total_length'] > 0:
                    stats['gc_content'] = (gc_count / stats['total_length']) * 100
                    stats['n_content'] = (n_count / stats['total_length']) * 100
            
            return stats
            
        except Exception as e:
            self.logger.error(f"计算序列统计信息时出错 | Error calculating sequence statistics: {e}")
            return stats

class RepeatResultParser:
    """重复序列结果解析器 | Repeat Result Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_repeatmasker_out(self, out_file: str) -> List[Dict]:
        """解析RepeatMasker .out文件 | Parse RepeatMasker .out file"""
        repeats = []
        
        try:
            with open(out_file, 'r') as f:
                # 跳过头部行 | Skip header lines
                for _ in range(3):
                    next(f)
                
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 15:
                        repeat = {
                            'score': int(parts[0]),
                            'divergence': float(parts[1]),
                            'deletion': float(parts[2]),
                            'insertion': float(parts[3]),
                            'query_seq': parts[4],
                            'query_start': int(parts[5]),
                            'query_end': int(parts[6]),
                            'query_left': parts[7],
                            'strand': parts[8],
                            'repeat_name': parts[9],
                            'repeat_class': parts[10],
                            'repeat_start': parts[11],
                            'repeat_end': parts[12],
                            'repeat_left': parts[13],
                            'id': parts[14] if len(parts) > 14 else ''
                        }
                        repeats.append(repeat)
            
            self.logger.info(f"解析了 {len(repeats)} 个重复序列 | Parsed {len(repeats)} repeat sequences")
            return repeats
            
        except Exception as e:
            self.logger.error(f"解析RepeatMasker结果时出错 | Error parsing RepeatMasker results: {e}")
            return []
    
    def parse_trf_output(self, trf_file: str) -> List[Dict]:
        """解析TRF输出文件 | Parse TRF output file"""
        repeats = []
        
        try:
            with open(trf_file, 'r') as f:
                current_sequence = ""
                
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('Sequence:'):
                        current_sequence = line.split(':')[1].strip()
                    elif line and not line.startswith('@'):
                        parts = line.split()
                        if len(parts) >= 15:
                            repeat = {
                                'sequence': current_sequence,
                                'start': int(parts[0]),
                                'end': int(parts[1]),
                                'period_size': int(parts[2]),
                                'copy_number': float(parts[3]),
                                'consensus_size': int(parts[4]),
                                'percent_matches': int(parts[5]),
                                'percent_indels': int(parts[6]),
                                'score': int(parts[7]),
                                'a_percent': int(parts[8]),
                                'c_percent': int(parts[9]),
                                'g_percent': int(parts[10]),
                                't_percent': int(parts[11]),
                                'entropy': float(parts[12]),
                                'consensus': parts[13],
                                'sequence_match': parts[14] if len(parts) > 14 else ''
                            }
                            repeats.append(repeat)
            
            self.logger.info(f"解析了 {len(repeats)} 个串联重复序列 | Parsed {len(repeats)} tandem repeats")
            return repeats
            
        except Exception as e:
            self.logger.error(f"解析TRF结果时出错 | Error parsing TRF results: {e}")
            return []
    
    def parse_edta_output(self, edta_dir: str) -> List[Dict]:
        """解析EDTA输出文件 | Parse EDTA output files"""
        repeats = []
        
        try:
            # 查找EDTA生成的注释文件 | Find EDTA annotation files
            import glob
            
            # EDTA通常生成以下文件：
            # *.EDTA.TEanno.gff3 - 主要的TE注释文件
            # *.EDTA.TElib.fa - 重复序列库
            # *.EDTA.TEanno.sum - 统计文件
            
            gff_files = glob.glob(os.path.join(edta_dir, "*.EDTA.TEanno.gff3"))
            
            for gff_file in gff_files:
                with open(gff_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('#') or not line:
                            continue
                        
                        parts = line.split('\t')
                        if len(parts) >= 9:
                            # 解析GFF3格式的TE注释
                            attributes = {}
                            for attr in parts[8].split(';'):
                                if '=' in attr:
                                    key, value = attr.split('=', 1)
                                    attributes[key] = value
                            
                            repeat = {
                                'sequence_id': parts[0],
                                'source': parts[1],
                                'feature_type': parts[2],
                                'start': int(parts[3]),
                                'end': int(parts[4]),
                                'score': parts[5] if parts[5] != '.' else None,
                                'strand': parts[6],
                                'phase': parts[7],
                                'attributes': attributes,
                                'te_classification': attributes.get('Classification', 'Unknown'),
                                'te_name': attributes.get('Name', 'Unknown')
                            }
                            repeats.append(repeat)
            
            self.logger.info(f"解析了 {len(repeats)} 个EDTA转座元件 | Parsed {len(repeats)} EDTA transposable elements")
            return repeats
            
        except Exception as e:
            self.logger.error(f"解析EDTA结果时出错 | Error parsing EDTA results: {e}")
            return []
