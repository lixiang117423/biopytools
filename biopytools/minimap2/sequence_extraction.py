"""
序列提取模块 | Sequence Extraction Module
"""

import os
from .utils import CommandRunner

class SequenceExtractor:
    """序列提取器 | Sequence Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def extract_unmapped_sequences(self, unmapped_regions):
        """提取未比对序列 | Extract unmapped sequences"""
        query_genome = self.config.query_genome
        bed_file = self.config.bed_file
        unmapped_fasta = self.config.unmapped_fasta
        seqkit_path = self.config.seqkit_path
        
        self.logger.info(f"使用seqkit提取未比对序列 | Using seqkit to extract unmapped sequences")
        self.logger.info(f"输入基因组文件 | Input genome file: {query_genome}")
        self.logger.info(f"BED文件 | BED file: {bed_file}")
        self.logger.info(f"输出序列文件 | Output sequence file: {unmapped_fasta}")
        
        if not unmapped_regions:
            self.logger.warning("没有找到未比对区间，跳过序列提取 | No unmapped regions found, skipping sequence extraction")
            return True
        
        # 创建临时序列文件 | Create temporary sequence file
        temp_fasta = unmapped_fasta + ".tmp"
        
        # 使用seqkit subseq提取序列 | Use seqkit subseq to extract sequences
        command = f"{seqkit_path} subseq --bed {bed_file} {query_genome} > {temp_fasta}"
        
        success = self.cmd_runner.run(command, "seqkit提取未比对序列 | seqkit extract unmapped sequences")
        
        if not success:
            return False
        
        # 修改序列ID格式并筛选长度 | Modify sequence ID format and filter by length
        return self._process_extracted_sequences(temp_fasta, unmapped_fasta, unmapped_regions)
    
    def _process_extracted_sequences(self, temp_fasta, final_fasta, unmapped_regions):
        """处理提取的序列 | Process extracted sequences"""
        self.logger.info("处理序列ID格式并筛选长度 | Processing sequence ID format and filtering by length")
        
        try:
            # 检查临时文件是否存在且不为空 | Check if temp file exists and is not empty
            if not os.path.exists(temp_fasta):
                self.logger.error(f"临时序列文件不存在 | Temporary sequence file does not exist: {temp_fasta}")
                return False
            
            if os.path.getsize(temp_fasta) == 0:
                self.logger.warning(f"临时序列文件为空 | Temporary sequence file is empty: {temp_fasta}")
                # 创建空的最终文件 | Create empty final file
                with open(final_fasta, 'w') as f:
                    pass
                return True
            
            # 创建区间映射字典 | Create interval mapping dictionary
            region_map = {}
            for region in unmapped_regions:
                key = f"{region['query_name']}_{region['start']}_{region['end']}"
                region_map[key] = region
            
            sequences_written = 0
            sequences_read = 0
            
            self.logger.info(f"开始处理临时文件 | Start processing temp file: {temp_fasta}")
            
            # 创建中间文件 | Create intermediate file
            intermediate_fasta = final_fasta + ".intermediate"
            
            with open(temp_fasta, 'r') as input_file, open(intermediate_fasta, 'w') as output_file:
                current_seq = ""
                current_header = ""
                
                for line_num, line in enumerate(input_file, 1):
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # 处理前一个序列 | Process previous sequence
                        if current_header and current_seq:
                            sequences_read += 1
                            if self._write_sequence_if_valid(output_file, current_header, current_seq, region_map):
                                sequences_written += 1
                            else:
                                self.logger.debug(f"序列被过滤 | Sequence filtered: {current_header}")
                        
                        current_header = line
                        current_seq = ""
                        
                        # 调试：显示前几个header | Debug: show first few headers
                        if sequences_read < 5:
                            self.logger.info(f"读取序列头 | Reading sequence header: {current_header}")
                            
                    else:
                        current_seq += line
                
                # 处理最后一个序列 | Process last sequence
                if current_header and current_seq:
                    sequences_read += 1
                    if self._write_sequence_if_valid(output_file, current_header, current_seq, region_map):
                        sequences_written += 1
                    else:
                        self.logger.debug(f"序列被过滤 | Sequence filtered: {current_header}")
            
            # 使用seqkit格式化序列长度为60bp每行 | Use seqkit to format sequence length to 60bp per line
            if sequences_written > 0:
                self.logger.info("格式化序列长度为60bp每行 | Formatting sequence length to 60bp per line")
                seqkit_path = self.config.seqkit_path
                format_command = f"{seqkit_path} seq -w 60 {intermediate_fasta} > {final_fasta}"
                
                success = self.cmd_runner.run(format_command, "格式化序列长度 | Format sequence length")
                if not success:
                    self.logger.warning("序列格式化失败，使用原始文件 | Sequence formatting failed, using original file")
                    # 如果格式化失败，直接复制中间文件 | If formatting fails, copy intermediate file
                    import shutil
                    shutil.move(intermediate_fasta, final_fasta)
                else:
                    self.logger.info("序列格式化完成 | Sequence formatting completed")
            else:
                # 如果没有序列，创建空文件 | If no sequences, create empty file
                with open(final_fasta, 'w') as f:
                    pass
            
            # 删除临时文件 | Remove temporary files
            if os.path.exists(temp_fasta):
                os.remove(temp_fasta)
            if os.path.exists(intermediate_fasta):
                os.remove(intermediate_fasta)
            
            self.logger.info(f"总共读取 {sequences_read} 个序列 | Total sequences read: {sequences_read}")
            self.logger.info(f"成功提取并处理 {sequences_written} 个序列 | Successfully extracted and processed {sequences_written} sequences")
            self.logger.info(f"最终序列文件 | Final sequence file: {final_fasta}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"处理序列时出错 | Error processing sequences: {e}")
            return False
    
    def _write_sequence_if_valid(self, output_file, header, sequence, region_map):
        """如果序列有效则写入文件 | Write sequence to file if valid"""
        # 调试：显示header格式 | Debug: show header format
        self.logger.debug(f"处理序列头 | Processing header: {header}")
        
        # 解析seqkit输出的header格式 | Parse seqkit output header format
        # seqkit subseq的输出格式：>sequence_name_start-end:.
        # 需要处理末尾的 :. 
        header_clean = header[1:]  # 移除开头的 > | Remove leading >
        
        # 移除末尾的 :. 如果存在 | Remove trailing :. if exists
        if header_clean.endswith(':.'):
            header_clean = header_clean[:-2]
        
        # 解析格式：query_name_start-end | Parse format: query_name_start-end
        if '_' in header_clean:
            # 从右边找最后一个下划线，分离序列名和坐标 | Find last underscore from right to separate name and coordinates
            last_underscore = header_clean.rfind('_')
            query_name = header_clean[:last_underscore]
            coord_part = header_clean[last_underscore + 1:]
            
            if '-' in coord_part:
                try:
                    start_str, end_str = coord_part.split('-')
                    start = int(start_str)
                    end = int(end_str)
                    
                    # 检查序列长度 | Check sequence length
                    seq_length = len(sequence)
                    self.logger.debug(f"序列长度 | Sequence length: {seq_length} for {query_name}:{start}-{end}")
                    
                    if seq_length >= self.config.min_unmapped_length:
                        # 生成新的序列ID | Generate new sequence ID
                        new_id = f"{query_name}:{start}-{end}"
                        output_file.write(f">{new_id}\n")
                        output_file.write(f"{sequence}\n")
                        self.logger.debug(f"写入序列 | Written sequence: {new_id} (length: {seq_length})")
                        return True
                    else:
                        self.logger.debug(f"序列长度不足 | Sequence too short: {seq_length} < {self.config.min_unmapped_length}")
                        
                except ValueError as e:
                    self.logger.debug(f"解析坐标失败 | Failed to parse coordinates: {e}")
        
        # 如果上面的解析失败，记录原始header | If above parsing fails, log original header
        self.logger.debug(f"无法解析header格式 | Cannot parse header format: {header}")
        return False
