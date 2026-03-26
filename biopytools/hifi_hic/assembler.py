"""
基因组组装执行模块 |Genome Assembly Execution Module
"""

import os
import glob
from pathlib import Path

try:
    from utils import run_command, generate_contig_reads_map_from_gfa
except ImportError:
    from biopytools.genome_assembler.utils import run_command, generate_contig_reads_map_from_gfa

class HifiasmAssembler:
    """Hifiasm组装器|Hifiasm Assembler"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def run_assembly(self) -> bool:
        """运行hifiasm组装|Run hifiasm assembly"""

        # 根据是否有Hi-C数据选择不同的模式|Choose different mode based on Hi-C data availability
        if self.config.has_hic:
            self.logger.info("运行hifiasm组装(HiFi + Hi-C模式)|Running hifiasm assembly (HiFi + Hi-C mode)")
            mode_str = "HiFi+Hi-C"
        else:
            self.logger.info("运行hifiasm组装(仅HiFi模式)|Running hifiasm assembly (HiFi only mode)")
            mode_str = "HiFi only"

        # 构建hifiasm命令|Build hifiasm command
        cmd_parts = [
            "hifiasm",
            f"-o {self.config.prefix}",
            f"-t {self.config.threads}",
            f"--hg-size {self.config.genome_size}",
            f"--n-hap {self.config.n_hap}",
        ]

        # 如果有Hi-C数据，添加Hi-C参数|Add Hi-C parameters if Hi-C data is available
        if self.config.has_hic:
            cmd_parts.extend([
                f"--h1 {self.config.hic_r1}",
                f"--h2 {self.config.hic_r2}",
            ])

        # 添加 purge level 参数|Add purge level parameter
        if self.config.purge_level is not None:
            cmd_parts.append(f"-l {self.config.purge_level}")

        # 添加 hom_cov 参数|Add homozygous coverage parameter
        if self.config.hom_cov is not None:
            cmd_parts.append(f"--hom-cov {self.config.hom_cov}")

        cmd_parts.append(self.config.hifi_data)

        cmd = " ".join(cmd_parts)

        # 执行命令并记录日志|Execute command and log
        log_file = os.path.join(self.config.log_dir, "hifiasm_assembly.log")
        full_cmd = f"{cmd} 2>&1|tee {log_file}"

        success = run_command(
            full_cmd,
            self.logger,
            work_dir=self.config.raw_dir,
            capture_output=False
        )

        if success:
            self.logger.info(f"组装完成({mode_str})|Assembly completed({mode_str})")
            return True
        else:
            self.logger.error(f"组装失败({mode_str})，请检查错误信息|Assembly failed({mode_str}), check error messages")
            return False

    def convert_gfa_to_fasta(self) -> dict:
        """转换GFA为FASTA格式|Convert GFA to FASTA format"""
        self.logger.info("转换GFA为FASTA格式|Converting GFA to FASTA format")

        # 根据是否有Hi-C数据定义需要转换的GFA文件|Define GFA files to convert based on Hi-C data availability
        if self.config.has_hic:
            gfa_files = {
                f"{self.config.prefix}.hic.hap1.p_ctg.gfa": f"{self.config.prefix}.hap1.fa",
                f"{self.config.prefix}.hic.hap2.p_ctg.gfa": f"{self.config.prefix}.hap2.fa",
                f"{self.config.prefix}.hic.p_ctg.gfa": f"{self.config.prefix}.primary.fa",
                f"{self.config.prefix}.hic.a_ctg.gfa": f"{self.config.prefix}.alternate.fa"
            }
        else:
            gfa_files = {
                f"{self.config.prefix}.hap1.p_ctg.gfa": f"{self.config.prefix}.hap1.fa",
                f"{self.config.prefix}.hap2.p_ctg.gfa": f"{self.config.prefix}.hap2.fa",
                f"{self.config.prefix}.p_ctg.gfa": f"{self.config.prefix}.primary.fa",
                f"{self.config.prefix}.a_ctg.gfa": f"{self.config.prefix}.alternate.fa"
            }

        results = {}

        for gfa_file, fasta_name in gfa_files.items():
            gfa_path = os.path.join(self.config.raw_dir, gfa_file)
            fasta_path = os.path.join(self.config.fasta_dir, fasta_name)

            # 如果标准文件名不存在，尝试带有.bp.前缀的文件名
            # If standard filename doesn't exist, try filename with .bp. prefix
            if not os.path.exists(gfa_path):
                # 对于有Hi-C的情况: {prefix}.hic.hap1.p_ctg.gfa -> {prefix}.bp.hic.hap1.p_ctg.gfa
                # 对于无Hi-C的情况: {prefix}.hap1.p_ctg.gfa -> {prefix}.bp.hap1.p_ctg.gfa
                # hifiasm使用purge-dups时会在prefix后添加.bp.
                if self.config.has_hic:
                    # 尝试在 prefix 后添加 .bp.
                    bp_gfa_file = gfa_file.replace(f"{self.config.prefix}.hic.", f"{self.config.prefix}.bp.hic.")
                else:
                    # 尝试在 prefix 后添加 .bp.
                    bp_gfa_file = gfa_file.replace(f"{self.config.prefix}.", f"{self.config.prefix}.bp.")

                bp_gfa_path = os.path.join(self.config.raw_dir, bp_gfa_file)
                if os.path.exists(bp_gfa_path):
                    gfa_path = bp_gfa_path
                    gfa_file = bp_gfa_file
                    self.logger.info(f"使用带.bp.前缀的文件|Using .bp. prefixed file: {gfa_file}")

            if os.path.exists(gfa_path):
                self.logger.info(f"转换|Converting: {gfa_file} -> {fasta_name}")

                # 使用awk转换并用seqkit格式化序列行|Convert using awk and format with seqkit
                cmd = f"""awk '/^S/{{print ">"$2; print $3}}' {gfa_path} | seqkit seq -w 60 - > {fasta_path}"""

                if run_command(cmd, self.logger):
                    stats = self._get_fasta_summary(fasta_path)
                    results[fasta_name] = stats
                    self.logger.info(f"转换成功|Successfully converted: {fasta_name}")
                    self.logger.info(f"序列数|Sequences: {stats['num_seqs']}")
                    self.logger.info(f"总长度|Total length: {stats['total_len']:,} bp")
                else:
                    self.logger.error(f"转换失败|Failed to convert: {gfa_file}")
            else:
                self.logger.warning(f"文件不存在|File not found: {gfa_file}")

        return results
    
    def _get_fasta_summary(self, fasta_file: str) -> dict:
        """获取FASTA文件摘要|Get FASTA file summary"""
        if not os.path.exists(fasta_file):
            return {}

        # 使用shell命令快速统计|Use shell commands for quick statistics
        cmd = f"""grep -c "^>" {fasta_file}"""
        num_seqs = os.popen(cmd).read().strip()

        cmd = f"""awk '/^>/{{next}}{{sum+=length($0)}}END{{print sum}}' {fasta_file}"""
        total_len = os.popen(cmd).read().strip()

        return {
            'num_seqs': int(num_seqs) if num_seqs else 0,
            'total_len': int(total_len) if total_len else 0
        }

    def generate_contig_reads_mapping(self) -> bool:
        """
        生成contig-reads映射文件|Generate contig-reads mapping file

        从GFA文件中解析read alignment信息，生成contig到reads的映射
        Parse read alignment from GFA files, generate contig to reads mapping

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("生成Contig-Reads映射|Generating Contig-Reads mapping")
        self.logger.info("=" * 60)

        # 根据是否有Hi-C数据定义GFA文件|Define GFA files based on Hi-C data availability
        if self.config.has_hic:
            gfa_prefix = f"{self.config.prefix}.hic"
        else:
            gfa_prefix = f"{self.config.prefix}"

        # 定义需要处理的GFA文件|Define GFA files to process
        gfa_files = [
            f"{gfa_prefix}.p_ctg.gfa",   # primary contigs
            f"{gfa_prefix}.hap1.p_ctg.gfa",  # hap1 contigs
            f"{gfa_prefix}.hap2.p_ctg.gfa",  # hap2 contigs
        ]

        success_count = 0

        for gfa_file in gfa_files:
            gfa_path = os.path.join(self.config.raw_dir, gfa_file)

            # 如果标准文件名不存在，尝试带有.bp.前缀的文件名
            # If standard filename doesn't exist, try filename with .bp. prefix
            actual_gfa_file = gfa_file
            if not os.path.exists(gfa_path):
                # 尝试 .bp. 前缀版本
                # Try with .bp. prefix
                if self.config.has_hic:
                    bp_gfa_file = gfa_file.replace(f"{self.config.prefix}.hic.", f"{self.config.prefix}.bp.hic.")
                else:
                    bp_gfa_file = gfa_file.replace(f"{self.config.prefix}.", f"{self.config.prefix}.bp.")

                bp_gfa_path = os.path.join(self.config.raw_dir, bp_gfa_file)
                if os.path.exists(bp_gfa_path):
                    gfa_path = bp_gfa_path
                    actual_gfa_file = bp_gfa_file
                    self.logger.info(f"使用带.bp.前缀的文件|Using .bp. prefixed file: {actual_gfa_file}")

            if not os.path.exists(gfa_path):
                self.logger.info(f"跳过|Skip: {gfa_file} (不存在|not exist)")
                continue

            # 确定输出文件名|Determine output file name
            # 例如: OV53.hic.p_ctg.gfa -> OV53.p_ctg.contig_reads.tsv
            #      OV53.hic.hap1.p_ctg.gfa -> OV53.hap1.p_ctg.contig_reads.tsv
            #      OV53.bp.p_ctg.gfa -> OV53.p_ctg.contig_reads.tsv (去掉.bp.)
            # e.g.: OV53.hic.p_ctg.gfa -> OV53.p_ctg.contig_reads.tsv
            # 去掉所有特殊前缀（.hic., .bp., .bp.hic.）
            base_name = actual_gfa_file
            # 移除.hic.或.bp.hic.
            base_name = base_name.replace(f"{self.config.prefix}.bp.hic.", f"{self.config.prefix}.")
            base_name = base_name.replace(f"{self.config.prefix}.hic.", f"{self.config.prefix}.")
            # 移除剩余的.bp.
            base_name = base_name.replace(f"{self.config.prefix}.bp.", f"{self.config.prefix}.")

            # 去掉.gfa后缀，添加.contig_reads.tsv|Remove .gfa suffix, add .contig_reads.tsv
            output_name = base_name.replace(".gfa", ".contig_reads.tsv")
            output_path = os.path.join(self.config.fasta_dir, output_name)

            self.logger.info(f"处理|Processing: {actual_gfa_file}")

            if generate_contig_reads_map_from_gfa(gfa_path, output_path, self.logger):
                success_count += 1

        self.logger.info(f"成功处理|Successfully processed: {success_count}/{len(gfa_files)} 个文件")
        self.logger.info("=" * 60)

        return success_count > 0
