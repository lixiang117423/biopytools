"""
基因组组装工具函数模块 |Genome Assembly Utility Functions Module
"""

import subprocess
import os
from pathlib import Path

def check_dependencies(logger):
    """检查依赖软件 |Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        result = subprocess.run(['hifiasm', '-V'],
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            logger.info(f"hifiasm 可用|available: {result.stdout.strip()}")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        logger.error("hifiasm 未安装或不在PATH中|not installed or not in PATH")
        raise RuntimeError("缺少必需软件: hifiasm|Missing required software: hifiasm")

def run_command(cmd: str, logger, work_dir: str = None, capture_output: bool = False) -> bool:
    """执行命令 |Execute command"""
    try:
        logger.info(f"执行命令|Executing command: {cmd}")

        if work_dir:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=work_dir,
                capture_output=capture_output,
                text=True,
                check=True
            )
        else:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=capture_output,
                text=True,
                check=True
            )

        if capture_output and result.stdout:
            logger.info(f"输出|Output: {result.stdout}")

        logger.info("命令执行成功|Command executed successfully")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False

def get_fasta_stats(fasta_file: str) -> dict:
    """获取FASTA文件统计信息|Get FASTA file statistics"""
    if not os.path.exists(fasta_file):
        return {}
    
    stats = {}
    
    # 计算序列数和总长度
    with open(fasta_file) as f:
        seq_lengths = []
        current_seq = []
        
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    seq_lengths.append(len(''.join(current_seq)))
                    current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_seq:
            seq_lengths.append(len(''.join(current_seq)))
    
    stats['num_seqs'] = len(seq_lengths)
    stats['total_len'] = sum(seq_lengths)
    
    # 计算最长序列
    stats['max_len'] = max(seq_lengths) if seq_lengths else 0
    
    # 计算N50
    if seq_lengths:
        sorted_lengths = sorted(seq_lengths, reverse=True)
        total_len = sum(sorted_lengths)
        half_len = total_len / 2
        
        cum_len = 0
        for length in sorted_lengths:
            cum_len += length
            if cum_len >= half_len:
                stats['n50'] = length
                break
    
    return stats

def format_time(seconds: int) -> str:
    """格式化时间|Format time"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours}小时 {minutes}分钟 {seconds}秒"


def generate_contig_reads_map_from_gfa(gfa_file: str, output_file: str, logger) -> bool:
    """
    从GFA文件生成contig-reads映射文件|Generate contig-reads mapping file from GFA

    解析GFA文件中的A行(read alignments)，生成contig到reads的映射文件
    Parse A lines (read alignments) in GFA file, generate contig to reads mapping

    Args:
        gfa_file: GFA文件路径|GFA file path
        output_file: 输出映射文件路径|Output mapping file path
        logger: 日志对象|Logger object

    Returns:
        bool: 是否成功|Whether successful
    """
    logger.info("=" * 60)
    logger.info("生成Contig-Reads映射文件|Generating Contig-Reads mapping file")
    logger.info("=" * 60)

    if not os.path.exists(gfa_file):
        logger.warning(f"GFA文件不存在|GFA file not found: {gfa_file}")
        return False

    try:
        # 第一步：解析GFA文件中的A行|Step 1: Parse A lines in GFA file
        logger.info(f"解析GFA文件|Parsing GFA file: {gfa_file}")

        contig_to_reads = {}
        read_count = 0

        with open(gfa_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('A\t'):
                    continue

                # A行格式: A\tcontig_id\tposition\tstrand\tread_name\tread_start\tread_end\ttags
                # A line format: A\tcontig_id\tposition\tstrand\tread_name\tread_start\tread_end\ttags
                parts = line.split('\t')
                if len(parts) >= 5:
                    contig_id = parts[1]
                    read_name = parts[4]

                    if contig_id not in contig_to_reads:
                        contig_to_reads[contig_id] = []

                    # 避免重复添加同一个read到同一个contig
                    # Avoid adding duplicate read to same contig
                    if read_name not in contig_to_reads[contig_id]:
                        contig_to_reads[contig_id].append(read_name)
                        read_count += 1

        logger.info(f"  共处理|Total alignments processed: {read_count:,}")
        logger.info(f"  共发现|Total contigs found: {len(contig_to_reads):,}")

        # 第二步：写入输出文件（按contig ID排序）
        # Step 2: Write to output file (sorted by contig ID)
        logger.info(f"写入映射文件|Writing mapping file: {output_file}")

        # 按数字部分排序|Sort by numeric part
        sorted_contigs = sorted(contig_to_reads.keys(),
                                key=lambda x: int(''.join(filter(str.isdigit, x))) if ''.join(filter(str.isdigit, x)) else 0)

        with open(output_file, 'w') as f_out:
            # 写入header|Write header
            f_out.write("#contig_id\tread_name\n")

            for contig_id in sorted_contigs:
                for read_name in contig_to_reads[contig_id]:
                    f_out.write(f"{contig_id}\t{read_name}\n")

        # 统计信息|Statistics
        logger.info("映射统计|Mapping statistics:")
        logger.info(f"  总contig数|Total contigs: {len(contig_to_reads):,}")

        # 计算每个contig的reads数量|Calculate reads per contig
        reads_per_contig = [len(reads) for reads in contig_to_reads.values()]
        if reads_per_contig:
            logger.info(f"  平均reads数/contig|Avg reads/contig: {sum(reads_per_contig) / len(reads_per_contig):.1f}")
            logger.info(f"  最少reads数|Min reads: {min(reads_per_contig):,}")
            logger.info(f"  最多reads数|Max reads: {max(reads_per_contig):,}")

        logger.info(f"  输出文件|Output file: {output_file}")
        logger.info("=" * 60)
        logger.info("Contig-Reads映射文件生成完成|Contig-Reads mapping file generated successfully")
        logger.info("=" * 60)

        return True

    except Exception as e:
        logger.error(f"生成contig-reads映射失败|Failed to generate contig-reads mapping: {str(e)}")
        return False
