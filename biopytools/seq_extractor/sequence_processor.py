"""
序列处理模块 🧬 | Sequence Processing Module
"""

import subprocess
from pathlib import Path
from typing import List, Tuple, Dict
from .utils import CommandRunner, parse_regions_file, format_region_name

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

class SequenceProcessor:
    """序列处理器 🔬 | Sequence Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 检查Biopython是否可用 | Check if Biopython is available
        if not BIOPYTHON_AVAILABLE and self.config.sequence_type == "protein":
            self.logger.warning("Biopython未安装，将尝试使用基础方法处理蛋白质序列 | "
                              "Biopython not installed, will try basic methods for protein sequences")
    
    def extract_sequences(self) -> bool:
        """提取序列 🎯 | Extract sequences"""
        try:
            # 解析区域文件 | Parse regions file
            self.logger.info(f"解析区域文件 📊 | Parsing regions file: {self.config.regions_file}")
            regions = parse_regions_file(self.config.regions_file)
            self.logger.info(f"找到 {len(regions)} 个区域 | Found {len(regions)} regions")
            
            # 检查是否有负链信息 | Check if there are negative strand regions
            negative_strands = [r for r in regions if r[3] == '-']
            if negative_strands:
                self.logger.info(f"检测到 {len(negative_strands)} 个负链区域，将自动进行反向互补 🔄 | "
                               f"Detected {len(negative_strands)} negative strand regions, will auto reverse complement")
            
            if self.config.sequence_type == "dna":
                return self._extract_dna_sequences(regions)
            else:
                return self._extract_protein_sequences(regions)
                
        except Exception as e:
            self.logger.error(f"序列提取失败 ❌ | Sequence extraction failed: {e}")
            return False
    
    def _extract_dna_sequences(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """提取DNA序列 🧬 | Extract DNA sequences"""
        self.logger.info("使用samtools提取DNA序列 | Using samtools to extract DNA sequences")
        
        # 检查是否需要建立索引 | Check if index is needed
        fai_file = f"{self.config.sequence_file}.fai"
        if not Path(fai_file).exists():
            self.logger.info("建立序列索引 🗂️ | Building sequence index")
            index_cmd = f"{self.config.samtools_path} faidx {self.config.sequence_file}"
            if not self.cmd_runner.run(index_cmd, "建立FASTA索引 | Building FASTA index"):
                return False
        
        if self.config.merge_output:
            return self._extract_dna_merged(regions)
        else:
            return self._extract_dna_separate(regions)
    
    def _extract_dna_merged(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """合并提取DNA序列 📦 | Extract DNA sequences merged"""
        # 构建samtools faidx命令 | Build samtools faidx command
        region_strings = []
        for chrom, start, end, strand in regions:
            # samtools faidx使用1-based坐标 | samtools faidx uses 1-based coordinates
            region_strings.append(f"{chrom}:{start}-{end}")
        
        # 分批处理以避免命令行过长 | Process in batches to avoid command line too long
        batch_size = 100
        all_sequences = []
        
        for i in range(0, len(region_strings), batch_size):
            batch = region_strings[i:i + batch_size]
            batch_regions = regions[i:i + batch_size]
            
            self.logger.info(f"处理批次 {i // batch_size + 1}/{(len(region_strings) - 1) // batch_size + 1} "
                           f"| Processing batch {i // batch_size + 1}/{(len(region_strings) - 1) // batch_size + 1}")
            
            temp_output = f"{self.config.output_file}.batch_{i // batch_size + 1}.tmp"
            cmd = f"{self.config.samtools_path} faidx {self.config.sequence_file} {' '.join(batch)} > {temp_output}"
            
            if not self.cmd_runner.run(cmd, f"提取第{i // batch_size + 1}批序列 | Extract batch {i // batch_size + 1} sequences"):
                return False
            
            # 读取临时文件并重新格式化 | Read temp file and reformat
            batch_sequences = self._process_extracted_sequences(temp_output, batch_regions)
            all_sequences.extend(batch_sequences)
            
            # 删除临时文件 | Remove temp file
            Path(temp_output).unlink(missing_ok=True)
        
        # 写入最终输出文件 | Write final output file
        return self._write_final_output(all_sequences)
    
    def _extract_dna_separate(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """分别提取DNA序列 📂 | Extract DNA sequences separately"""
        output_dir = Path(self.config.output_file).parent
        base_name = Path(self.config.output_file).stem
        
        for i, (chrom, start, end, strand) in enumerate(regions, 1):
            region_str = f"{chrom}:{start}-{end}"
            strand_suffix = "_rc" if strand == '-' else ""
            output_file = output_dir / f"{base_name}_{i:04d}_{chrom}_{start}_{end}{strand_suffix}.fasta"
            
            cmd = f"{self.config.samtools_path} faidx {self.config.sequence_file} {region_str} > {output_file}"
            
            if not self.cmd_runner.run(cmd, f"提取区域 {i}/{len(regions)}: {region_str} ({strand}链) | Extract region {i}/{len(regions)}: {region_str} ({strand} strand)"):
                return False
            
            # 如果是负链，需要后处理进行反向互补 | Post-process reverse complement for negative strand
            if strand == '-':
                self._post_process_reverse_complement(output_file)
        
        self.logger.info(f"所有序列已分别保存到 {output_dir} | All sequences saved separately to {output_dir}")
        return True
    
    def _post_process_reverse_complement(self, fasta_file: Path):
        """后处理反向互补 🔄 | Post-process reverse complement"""
        sequences = []
        with open(fasta_file, 'r') as f:
            current_header = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header is not None:
                        seq = ''.join(current_seq)
                        rev_comp_seq = self._reverse_complement(seq)
                        new_header = current_header.replace('>', '>rc_')
                        sequences.append((new_header, rev_comp_seq))
                    
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # 处理最后一个序列 | Process last sequence
            if current_header is not None:
                seq = ''.join(current_seq)
                rev_comp_seq = self._reverse_complement(seq)
                new_header = current_header.replace('>', '>rc_')
                sequences.append((new_header, rev_comp_seq))
        
        # 重写文件 | Rewrite file
        with open(fasta_file, 'w') as f:
            for header, seq in sequences:
                f.write(f"{header}\n")
                for i in range(0, len(seq), self.config.line_width):
                    f.write(f"{seq[i:i + self.config.line_width]}\n")
    
    def _extract_protein_sequences(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """提取蛋白质序列 🧪 | Extract protein sequences"""
        if BIOPYTHON_AVAILABLE:
            return self._extract_protein_with_biopython(regions)
        else:
            return self._extract_protein_basic(regions)
    
    def _extract_protein_with_biopython(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """使用Biopython提取蛋白质序列 🐍 | Extract protein sequences with Biopython"""
        self.logger.info("使用Biopython提取蛋白质序列 | Using Biopython to extract protein sequences")
        
        # 读取序列文件 | Read sequence file
        sequences = {}
        for record in SeqIO.parse(self.config.sequence_file, "fasta"):
            sequences[record.id] = str(record.seq)
        
        extracted_sequences = []
        
        for chrom, start, end, strand in regions:
            if chrom not in sequences:
                self.logger.warning(f"序列 {chrom} 在文件中未找到 | Sequence {chrom} not found in file")
                continue
            
            seq = sequences[chrom]
            
            # 蛋白质序列使用1-based坐标 | Protein sequences use 1-based coordinates
            if start > len(seq) or end > len(seq):
                self.logger.warning(f"区域 {chrom}:{start}-{end} 超出序列长度 {len(seq)} | "
                                  f"Region {chrom}:{start}-{end} exceeds sequence length {len(seq)}")
                continue
            
            # 提取子序列 (转换为0-based索引) | Extract subsequence (convert to 0-based index)
            subseq = seq[start-1:end]
            
            # 格式化头部 | Format header
            if self.config.include_headers:
                header = format_region_name(chrom, start, end, chrom)
                if strand == '-':
                    header = header.replace('>', '>rc_')
            else:
                prefix = "rc_" if strand == '-' else ""
                header = f">{prefix}region_{len(extracted_sequences) + 1}"
            
            extracted_sequences.append((header, subseq))
        
        return self._write_final_output(extracted_sequences)
    
    def _extract_protein_basic(self, regions: List[Tuple[str, int, int, str]]) -> bool:
        """基础方法提取蛋白质序列 📝 | Extract protein sequences with basic method"""
        self.logger.info("使用基础方法提取蛋白质序列 | Using basic method to extract protein sequences")
        
        # 读取FASTA文件 | Read FASTA file
        sequences = self._parse_fasta_basic(self.config.sequence_file)
        
        extracted_sequences = []
        
        for chrom, start, end, strand in regions:
            if chrom not in sequences:
                self.logger.warning(f"序列 {chrom} 在文件中未找到 | Sequence {chrom} not found in file")
                continue
            
            seq = sequences[chrom]
            
            if start > len(seq) or end > len(seq):
                self.logger.warning(f"区域 {chrom}:{start}-{end} 超出序列长度 {len(seq)} | "
                                  f"Region {chrom}:{start}-{end} exceeds sequence length {len(seq)}")
                continue
            
            # 提取子序列 | Extract subsequence
            subseq = seq[start-1:end]
            
            # 格式化头部 | Format header
            if self.config.include_headers:
                header = format_region_name(chrom, start, end, chrom)
                if strand == '-':
                    header = header.replace('>', '>rc_')
            else:
                prefix = "rc_" if strand == '-' else ""
                header = f">{prefix}region_{len(extracted_sequences) + 1}"
            
            extracted_sequences.append((header, subseq))
        
        return self._write_final_output(extracted_sequences)
    
    def _parse_fasta_basic(self, fasta_file: str) -> Dict[str, str]:
        """基础FASTA解析 📖 | Basic FASTA parsing"""
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # 提取序列ID | Extract sequence ID
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        
        # 处理最后一条序列 | Process last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    def _process_extracted_sequences(self, temp_file: str, regions: List[Tuple[str, int, int, str]]) -> List[Tuple[str, str]]:
        """处理提取的序列 ✨ | Process extracted sequences"""
        sequences = []
        current_header = None
        current_seq = []
        region_index = 0
        
        with open(temp_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存前一个序列 | Save previous sequence
                    if current_header is not None and region_index < len(regions):
                        chrom, start, end, strand = regions[region_index]
                        if self.config.include_headers:
                            new_header = format_region_name(chrom, start, end, current_header)
                            if strand == '-':
                                new_header = new_header.replace('>', '>rc_')
                        else:
                            prefix = "rc_" if strand == '-' else ""
                            new_header = f">{prefix}region_{region_index + 1}"
                        
                        seq = ''.join(current_seq)
                        
                        # 处理DNA特殊选项 | Handle DNA special options
                        if self.config.sequence_type == "dna":
                            # 负链自动反向互补 | Auto reverse complement for negative strand
                            if strand == '-' or self.config.reverse_complement:
                                seq = self._reverse_complement(seq)
                            if self.config.translate_dna:
                                seq = self._translate_dna(seq)
                        
                        sequences.append((new_header, seq))
                        region_index += 1
                    
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
        
        # 处理最后一个序列 | Process last sequence
        if current_header is not None and region_index < len(regions):
            chrom, start, end, strand = regions[region_index]
            if self.config.include_headers:
                new_header = format_region_name(chrom, start, end, current_header)
                if strand == '-':
                    new_header = new_header.replace('>', '>rc_')
            else:
                prefix = "rc_" if strand == '-' else ""
                new_header = f">{prefix}region_{region_index + 1}"
            
            seq = ''.join(current_seq)
            
            if self.config.sequence_type == "dna":
                if strand == '-' or self.config.reverse_complement:
                    seq = self._reverse_complement(seq)
                if self.config.translate_dna:
                    seq = self._translate_dna(seq)
            
            sequences.append((new_header, seq))
        
        return sequences
    
    def _reverse_complement(self, seq: str) -> str:
        """反向互补 🔄 | Reverse complement"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
        
        rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
        return rev_comp
    
    def _translate_dna(self, seq: str) -> str:
        """翻译DNA为蛋白质 🔄 | Translate DNA to protein"""
        if BIOPYTHON_AVAILABLE:
            from Bio.Seq import Seq
            return str(Seq(seq).translate())
        else:
            # 基础翻译表 | Basic translation table
            genetic_code = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
            }
            
            protein = []
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3].upper()
                amino_acid = genetic_code.get(codon, 'X')
                protein.append(amino_acid)
            
            return ''.join(protein)
    
    def _write_final_output(self, sequences: List[Tuple[str, str]]) -> bool:
        """写入最终输出文件 💾 | Write final output file"""
        try:
            with open(self.config.output_file, 'w') as f:
                for header, seq in sequences:
                    f.write(f"{header}\n")
                    
                    # 按指定宽度分行 | Break lines by specified width
                    for i in range(0, len(seq), self.config.line_width):
                        f.write(f"{seq[i:i + self.config.line_width]}\n")
            
            self.logger.info(f"✅ 成功提取 {len(sequences)} 条序列到 {self.config.output_file} | "
                           f"Successfully extracted {len(sequences)} sequences to {self.config.output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 写入输出文件失败 | Failed to write output file: {e}")
            return False
