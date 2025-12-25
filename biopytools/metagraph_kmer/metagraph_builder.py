"""
🏗️ MetaGraph库构建模块 | MetaGraph Library Construction Module
"""

import os
from pathlib import Path
from typing import List

class MetaGraphBuilder:
    """MetaGraph库构建器 | MetaGraph Library Builder"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_graph(self, ref_files: List[str], file_format: str) -> bool:
        """
        构建MetaGraph图 | Build MetaGraph graph
        
        Args:
            ref_files: 参考文件列表 | Reference file list
            file_format: 文件格式 ('fasta' or 'fastq')
        
        Returns:
            是否成功 | Success status
        """
        self.logger.info("🏗️  步骤1.1: 构建de Bruijn图 | Step 1.1: Building de Bruijn graph")
        
        # 构建文件列表 | Build file list
        files_str = " ".join(ref_files)
        
        # 构建MetaGraph命令 | Build MetaGraph command
        mode = "canonical" if self.config.canonical else "basic"
        
        # 使用绝对路径
        dbg_file_abs = self.config.dbg_file.resolve()

        cmd = (
            f"{self.config.metagraph_path} build "
            f"-k {self.config.kmer_length} "
            f"--mode {mode} "
            f"-p {self.config.threads} "
            f"-o {dbg_file_abs} "
            f"-v "
            f"{files_str}"
        )
                
        if not self.cmd_runner.run(cmd, "构建MetaGraph de Bruijn图"):
            raise RuntimeError("❌ MetaGraph图构建失败 | MetaGraph graph build failed")
        
        self.logger.info(f"✅ DBG文件已生成 | DBG file generated: {self.config.dbg_file}")
        return True
    
    def annotate_with_coordinates(self, ref_files: List[str], file_format: str) -> bool:
        """
        添加坐标注释 | Add coordinate annotations
        
        Args:
            ref_files: 参考文件列表 | Reference file list
            file_format: 文件格式 | File format
        
        Returns:
            是否成功 | Success status
        """
        if file_format != 'fasta':
            self.logger.warning("⚠️  只有FASTA格式支持坐标注释 | Only FASTA format supports coordinate annotation")
            return False
        
        self.logger.info("🏗️  步骤1.2: 添加坐标注释 | Step 1.2: Adding coordinate annotations")
        
        # 构建文件列表 | Build file list
        files_str = " ".join(ref_files)
        
        # 使用绝对路径
        dbg_file_abs = self.config.dbg_file.resolve()
        output_prefix_abs = (self.config.output_path / self.config.base_name).resolve()

        cmd = (
            f"{self.config.metagraph_path} annotate "
            f"-i {dbg_file_abs} "
            f"--coordinates "
            f"--anno-filename "
            f"-p {self.config.threads} "
            f"-o {output_prefix_abs} "
            f"-v "
            f"{files_str}"
        )
        
        if not self.cmd_runner.run(cmd, "添加MetaGraph坐标注释"):
            raise RuntimeError("❌ 坐标注释失败 | Coordinate annotation failed")
        
        self.logger.info(f"✅ 注释文件已生成 | Annotation file generated: {self.config.anno_file}")
        return True
    
    def query_coordinates(self, ref_files: List[str]) -> Path:
        """
        从MetaGraph查询坐标 | Query coordinates from MetaGraph
        
        Args:
            ref_files: 参考文件列表 | Reference file list
        
        Returns:
            坐标文件路径 | Coordinates file path
        """
        self.logger.info("🔍 步骤2: 从MetaGraph提取k-mer坐标 | Step 2: Extracting k-mer coordinates from MetaGraph")
        
        coords_raw_file = self.config.output_path / "ref_kmers_coords_raw.txt"

        # 使用绝对路径
        coords_raw_file_abs = coords_raw_file.resolve()
        dbg_file_abs = self.config.dbg_file.resolve()
        anno_file_abs = self.config.anno_file.resolve()

        # 构建文件列表 | Build file list
        files_str = " ".join(ref_files)

        cmd = (
            f"{self.config.metagraph_path} query "
            f"-i {dbg_file_abs} "
            f"-a {anno_file_abs} "
            f"--query-mode coords "
            f"--min-kmers-fraction-label 0 "
            f"-v "
            f"{files_str} "
            f"> {coords_raw_file_abs}"
        )
        
        if not self.cmd_runner.run(cmd, "查询k-mer坐标"):
            raise RuntimeError("❌ 坐标查询失败 | Coordinate query failed")
        
        self.logger.info(f"✅ 原始坐标已提取 | Raw coordinates extracted: {coords_raw_file}")
        
        # 解析并格式化坐标 | Parse and format coordinates
        coords_file = self._parse_coordinates(coords_raw_file)
        
        return coords_file
    
    def _parse_coordinates(self, raw_file: Path) -> Path:
        """
        解析MetaGraph坐标输出 | Parse MetaGraph coordinate output
        
        Args:
            raw_file: 原始坐标文件 | Raw coordinates file
        
        Returns:
            格式化的坐标文件 | Formatted coordinates file
        """
        self.logger.info("📝 解析并格式化坐标信息 | Parsing and formatting coordinates")
        
        coords_file = self.config.output_path / "ref_kmers_coords.tsv"
        
        kmer_count = 0
        
        with open(raw_file, 'r') as infile, open(coords_file, 'w') as outfile:
            # 写入表头 | Write header
            outfile.write("kmer\tseq_id\tstart\tend\tstrand\n")
            
            current_seq = None
            
            for line in infile:
                line = line.strip()
                
                # 序列名称行 | Sequence name line
                if line.startswith('>'):
                    current_seq = line[1:]
                    continue
                
                # 坐标行格式: filename-seqid:start-end
                # 例如: reference.fasta-seq1:100-131
                if current_seq and ':' in line and '-' in line:
                    try:
                        # 解析坐标 | Parse coordinates
                        parts = line.split(':')
                        if len(parts) == 2:
                            seq_id = parts[0].split('-')[-1]  # 获取序列ID
                            positions = parts[1].split('-')
                            if len(positions) == 2:
                                start = positions[0]
                                end = positions[1]
                                
                                # 这里暂时无法获取实际的k-mer序列，使用占位符
                                # 实际k-mer需要从原始序列中提取
                                outfile.write(f"PLACEHOLDER\t{seq_id}\t{start}\t{end}\t+\n")
                                kmer_count += 1
                                
                                if kmer_count % 1000000 == 0:
                                    self.logger.info(f"   已解析 | Parsed: {kmer_count:,} coordinates")
                    except Exception as e:
                        self.logger.debug(f"解析坐标行失败 | Failed to parse coordinate line: {line}")
                        continue
        
        self.logger.info(f"✅ 共解析 {kmer_count:,} 个坐标 | Total parsed: {kmer_count:,} coordinates")
        self.logger.info(f"📝 坐标文件已保存 | Coordinates saved: {coords_file}")
        
        return coords_file
    
    def extract_kmers_from_fasta(self, ref_files: List[str]) -> Path:
        """
        直接从FASTA提取k-mer（备选方案，更准确）
        Extract k-mers directly from FASTA (alternative, more accurate)
        """
        self.logger.info("🧬 从FASTA文件直接提取k-mer和坐标 | Extracting k-mers and coordinates directly from FASTA")
        
        from Bio import SeqIO
        import gzip
        
        coords_file = self.config.output_path / "ref_kmers_coords.tsv"
        
        with open(coords_file, 'w') as f:
            f.write("kmer\tseq_id\tstart\tend\tstrand\n")
            
            kmer_count = 0
            
            for fasta_file in ref_files:
                self.logger.info(f"   处理文件 | Processing: {os.path.basename(fasta_file)}")
                
                if fasta_file.endswith('.gz'):
                    handle = gzip.open(fasta_file, 'rt')
                else:
                    handle = open(fasta_file, 'r')
                
                for record in SeqIO.parse(handle, "fasta"):
                    seq_id = record.id
                    sequence = str(record.seq).upper()
                    
                    for i in range(len(sequence) - self.config.kmer_length + 1):
                        kmer = sequence[i:i + self.config.kmer_length]
                        
                        if 'N' in kmer:
                            continue
                        
                        start = i
                        end = i + self.config.kmer_length
                        
                        # Canonical k-mer处理
                        if self.config.canonical:
                            kmer_rc = self._reverse_complement(kmer)
                            if kmer <= kmer_rc:
                                canonical_kmer, strand = kmer, '+'
                            else:
                                canonical_kmer, strand = kmer_rc, '-'
                        else:
                            canonical_kmer, strand = kmer, '+'
                        
                        f.write(f"{canonical_kmer}\t{seq_id}\t{start}\t{end}\t{strand}\n")
                        kmer_count += 1
                        
                        if kmer_count % 1000000 == 0:
                            self.logger.info(f"   已提取 | Extracted: {kmer_count:,} k-mers")
                
                handle.close()
        
        self.logger.info(f"✅ 共提取 {kmer_count:,} 个k-mer | Total extracted: {kmer_count:,} k-mers")
        return coords_file
    
    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """计算反向互补 | Calculate reverse complement"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(seq))
