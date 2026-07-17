"""
基因组组装执行模块 |Genome Assembly Execution Module
"""

import os
import glob
from pathlib import Path

from .utils import run_shell_command, generate_contig_reads_map_from_gfa, get_fasta_stats

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

        # 构建hifiasm命令参数(list形式)|Build hifiasm command args (list)
        cmd_args = [
            "-o", self.config.prefix,
            "-t", str(self.config.threads),
            "--hg-size", self.config.genome_size,
            "--n-hap", str(self.config.n_hap),
        ]

        # 如果有Hi-C数据，添加Hi-C参数|Add Hi-C parameters if Hi-C data available
        if self.config.has_hic:
            cmd_args.extend(["--h1", self.config.hic_r1, "--h2", self.config.hic_r2])

        # 添加 purge level 参数|Add purge level parameter
        if self.config.purge_level is not None:
            cmd_args.extend(["-l", str(self.config.purge_level)])

        # 添加 hom_cov 参数|Add homozygous coverage parameter
        if self.config.hom_cov is not None:
            cmd_args.extend(["--hom-cov", str(self.config.hom_cov)])

        cmd_args.append(self.config.hifi_data)

        # hifiasm|tee:跨env管道(hifiasm conda + tee系统),完整路径+LD_LIBRARY_PATH(§13.2.1)
        # |hifiasm|tee: cross-env pipeline, full path + LD_LIBRARY_PATH (§13.2.1)
        log_file = os.path.join(self.config.log_dir, "hifiasm_assembly.log")
        shell_cmd = f"{self.config.hifiasm_path} {' '.join(cmd_args)} 2>&1 | tee {log_file}"

        success = run_shell_command(
            shell_cmd, [self.config.hifiasm_path],
            self.logger, work_dir=self.config.raw_dir
        )

        if success:
            self.logger.info(f"组装完成({mode_str})|Assembly completed({mode_str})")
            return True
        else:
            self.logger.error(f"组装失败({mode_str})，请检查错误信息|Assembly failed({mode_str}), check error messages")
            return False

    def _build_gfa_to_fasta_map(self) -> dict:
        """
        根据n_hap动态构建GFA到FASTA的文件映射|Build GFA to FASTA mapping dynamically based on n_hap

        Returns:
            dict: GFA文件名到FASTA文件名的映射|GFA filename to FASTA filename mapping
        """
        gfa_files = {}
        prefix = self.config.prefix

        # 根据是否有Hi-C数据构建前缀|Build prefix based on Hi-C data availability
        if self.config.has_hic:
            # 单倍型文件|Haplotype files
            for i in range(1, self.config.n_hap + 1):
                gfa_files[f"{prefix}.hic.hap{i}.p_ctg.gfa"] = f"{prefix}.hap{i}.fa"
            # primary和alternate文件|Primary and alternate files
            gfa_files[f"{prefix}.hic.p_ctg.gfa"] = f"{prefix}.primary.fa"
            gfa_files[f"{prefix}.hic.a_ctg.gfa"] = f"{prefix}.alternate.fa"
        else:
            # 单倍型文件|Haplotype files
            for i in range(1, self.config.n_hap + 1):
                gfa_files[f"{prefix}.hap{i}.p_ctg.gfa"] = f"{prefix}.hap{i}.fa"
            # primary和alternate文件|Primary and alternate files
            gfa_files[f"{prefix}.p_ctg.gfa"] = f"{prefix}.primary.fa"
            gfa_files[f"{prefix}.a_ctg.gfa"] = f"{prefix}.alternate.fa"

        return gfa_files

    def convert_gfa_to_fasta(self) -> dict:
        """转换GFA为FASTA格式|Convert GFA to FASTA format"""
        self.logger.info("转换GFA为FASTA格式|Converting GFA to FASTA format")

        # 根据n_hap动态构建文件映射|Build file mapping dynamically based on n_hap
        gfa_files = self._build_gfa_to_fasta_map()

        results = {}

        for gfa_file, fasta_name in gfa_files.items():
            gfa_path = os.path.join(self.config.raw_dir, gfa_file)
            fasta_path = os.path.join(self.config.fasta_dir, fasta_name)

            # 断点续传:fasta已存在则跳过转换只统计(§10.2)|Resume: skip conversion if fasta exists
            if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
                self.logger.info(f"FASTA已存在，跳过转换|FASTA exists, skipping conversion: {fasta_name}")
                results[fasta_name] = self._get_fasta_summary(fasta_path)
                continue

            # 如果标准文件名不存在，尝试带有.bp.前缀的文件名
            # If standard filename doesn't exist, try filename with .bp. prefix
            if not os.path.exists(gfa_path):
                # hifiasm使用purge-dups时会在prefix后添加.bp.
                # hifiasm adds .bp. after prefix when using purge-dups
                bp_gfa_file = gfa_file.replace(
                    f"{self.config.prefix}.", f"{self.config.prefix}.bp.", 1
                )

                bp_gfa_path = os.path.join(self.config.raw_dir, bp_gfa_file)
                if os.path.exists(bp_gfa_path):
                    gfa_path = bp_gfa_path
                    gfa_file = bp_gfa_file
                    self.logger.info(f"使用带.bp.前缀的文件|Using .bp. prefixed file: {gfa_file}")

            if os.path.exists(gfa_path):
                self.logger.info(f"转换|Converting: {gfa_file} -> {fasta_name}")

                # 使用awk转换并用seqkit格式化序列行(跨env管道:awk系统+seqkit conda)
                # |Convert using awk and format with seqkit (cross-env: awk system + seqkit conda)
                shell_cmd = f"""awk '/^S/{{print ">"$2; print $3}}' {gfa_path} | {self.config.seqkit_path} seq -w 60 - > {fasta_path}"""

                if run_shell_command(shell_cmd, [self.config.seqkit_path], self.logger):
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
        """获取FASTA文件摘要(Python原生统计,替代os.popen)|Get FASTA summary (native Python)"""
        return get_fasta_stats(fasta_file)

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

        # 根据是否有Hi-C数据定义GFA前缀|Define GFA prefix based on Hi-C data availability
        if self.config.has_hic:
            gfa_prefix = f"{self.config.prefix}.hic"
        else:
            gfa_prefix = f"{self.config.prefix}"

        # 根据n_hap动态构建需要处理的GFA文件列表|Build GFA file list dynamically based on n_hap
        gfa_files = [
            f"{gfa_prefix}.p_ctg.gfa",   # primary contigs
        ]
        for i in range(1, self.config.n_hap + 1):
            gfa_files.append(f"{gfa_prefix}.hap{i}.p_ctg.gfa")

        success_count = 0

        for gfa_file in gfa_files:
            gfa_path = os.path.join(self.config.raw_dir, gfa_file)

            # 如果标准文件名不存在，尝试带有.bp.前缀的文件名
            # If standard filename doesn't exist, try filename with .bp. prefix
            actual_gfa_file = gfa_file
            if not os.path.exists(gfa_path):
                # 尝试 .bp. 前缀版本
                # Try with .bp. prefix
                bp_gfa_file = gfa_file.replace(
                    f"{self.config.prefix}.", f"{self.config.prefix}.bp.", 1
                )

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
