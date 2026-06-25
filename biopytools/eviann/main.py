"""
EviAnn基因组注释主程序模块|EviAnn Genome Annotation Main Module
"""

import argparse
import sys
import os
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List

from .config import EviAnnConfig
from .utils import EviAnnLogger, build_conda_command, expand_path


class EviAnnotator:
    """EviAnn基因组注释主类|Main EviAnn Genome Annotation Class"""

    def __init__(self, **kwargs):
        """
        初始化EviAnn注释器|Initialize EviAnn annotator

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = EviAnnConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        logger_path = Path(self.config.genome).parent / f"{Path(self.config.genome).name}.eviann.log"
        self.logger_manager = EviAnnLogger(str(logger_path))
        self.logger = self.logger_manager.get_logger()

        # 处理RNA-seq输入：如果是FASTQ文件，生成EviAnn格式的描述文件
        # Process RNA-seq input: if FASTQ file, generate EviAnn format description file
        self.rnaseq_desc_file = self._prepare_rnaseq_input()

    def _prepare_rnaseq_input(self) -> str:
        """
        准备RNA-seq输入文件|Prepare RNA-seq input file

        支持多种输入方式：
        1. 已有EviAnn格式描述文件
        2. short_reads和long_reads参数（自动生成描述文件）
        3. rnaseq参数（如果是FASTQ文件，自动生成描述文件）

        Returns:
            str: RNA-seq描述文件路径|RNA-seq description file path
        """
        # 优先使用short_reads和long_reads|Priority: short_reads and long_reads
        if self.config.short_reads or self.config.long_reads:
            return self._generate_rnaseq_description()

        # 其次使用rnaseq参数|Next: use rnaseq parameter
        if self.config.rnaseq:
            return self._process_rnaseq_param()

        return ""

    def _generate_rnaseq_description(self) -> str:
        """
        根据short_reads和long_reads生成EviAnn格式的描述文件
        Generate EviAnn format description file from short_reads and long_reads

        Returns:
            str: 描述文件路径|Description file path
        """
        # 生成描述文件路径|Generate description file path
        desc_file = Path(self.config.genome).parent / "rnaseq_inputs.txt"

        lines = []

        # 处理二代转录组|Process short-reads
        if self.config.short_reads:
            short_path = Path(self.config.short_reads)
            if short_path.is_dir():
                # 扫描目录|Scan directory
                for item in sorted(short_path.iterdir()):
                    if item.suffix in ['.fq', '.fastq', '.fq.gz', '.fastq.gz']:
                        lines.append(f"{str(item.absolute())} fastq")
            elif short_path.is_file():
                # 单个文件|Single file
                lines.append(f"{str(short_path.absolute())} fastq")

        # 处理三代转录组|Process long-reads
        if self.config.long_reads:
            long_path = Path(self.config.long_reads)
            if long_path.is_dir():
                # 扫描目录|Scan directory
                for item in sorted(long_path.iterdir()):
                    if item.suffix in ['.fq', '.fastq', '.fq.gz', '.fastq.gz', '.fa', '.fa.gz', '.fasta', '.fasta.gz']:
                        lines.append(f"{str(item.absolute())} isoseq")
            elif long_path.is_file():
                # 单个文件|Single file
                lines.append(f"{str(long_path.absolute())} isoseq")

        if not lines:
            self.logger.warning("未找到有效的RNA-seq文件|No valid RNA-seq files found")
            return ""

        # 写入描述文件|Write description file
        with open(desc_file, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        self.logger.info(f"生成RNA-seq描述文件|Generated RNA-seq description file: {desc_file}")
        self.logger.info(f"  二代转录组|Short-reads: {self.config.short_reads}")
        self.logger.info(f"  三代转录组|Long-reads: {self.config.long_reads}")
        self.logger.info(f"  包含{len(lines)}个文件|Contains {len(lines)} files")

        return str(desc_file)

    def _process_rnaseq_param(self) -> str:
        """
        处理rnaseq参数（向后兼容）
        Process rnaseq parameter (for backward compatibility)

        Returns:
            str: RNA-seq描述文件路径|RNA-seq description file path
        """
        rnaseq_path = Path(self.config.rnaseq)

        # 如果不存在，返回原值|If not exists, return original value
        if not rnaseq_path.exists():
            return self.config.rnaseq

        # 如果是已存在的描述文件或目录，直接返回|If existing description file or directory, return directly
        if rnaseq_path.is_dir() or (rnaseq_path.is_file() and not self._is_fastq_file(rnaseq_path)):
            return self.config.rnaseq

        # 如果是FASTQ文件，生成描述文件|If FASTQ file, generate description file
        if self._is_fastq_file(rnaseq_path):
            # 生成描述文件路径|Generate description file path
            desc_file = Path(self.config.genome).parent / f"{rnaseq_path.stem}.rnaseq.txt"

            # 判断数据类型|Determine data type
            filename_lower = rnaseq_path.name.lower()
            if any(keyword in filename_lower for keyword in ['pacbio', 'isoseq', 'ont', 'nanopore', '三代', 'longread']):
                data_type = "isoseq"
            else:
                data_type = "fastq"

            # 获取绝对路径|Get absolute path
            abs_path = str(rnaseq_path.absolute())

            # 写入描述文件|Write description file
            with open(desc_file, 'w') as f:
                f.write(f"{abs_path} {data_type}\n")

            self.logger.info(f"生成RNA-seq描述文件|Generated RNA-seq description file: {desc_file}")
            self.logger.info(f"  数据类型|Data type: {data_type}")
            self.logger.info(f"  文件路径|File path: {abs_path}")

            return str(desc_file)

        return self.config.rnaseq

    def _is_fastq_file(self, file_path: Path) -> bool:
        """
        判断是否为FASTQ文件|Check if file is FASTQ format

        Args:
            file_path: 文件路径|File path

        Returns:
            bool: 是否为FASTQ文件|Is FASTQ file or not
        """
        suffixes = ''.join(file_path.suffixes).lower()
        return any(ext in suffixes for ext in ['.fq', '.fastq'])

    def build_eviann_command(self) -> List[str]:
        """
        构建EviAnn命令|Build EviAnn command

        Returns:
            完整的命令列表|Full command list
        """
        args = []

        # 必需参数|Required parameters
        args.extend(['-g', self.config.genome])

        # 数据输入参数|Data input parameters
        if self.rnaseq_desc_file:
            args.extend(['-r', self.rnaseq_desc_file])

        if self.config.transcripts:
            args.extend(['-e', self.config.transcripts])

        if self.config.proteins:
            args.extend(['-p', self.config.proteins])

        if self.config.uniprot:
            args.extend(['-s', self.config.uniprot])

        # 其他参数|Other parameters
        if self.config.threads != 12:
            args.extend(['-t', str(self.config.threads)])

        if self.config.max_intron:
            args.extend(['-m', str(self.config.max_intron)])

        if self.config.ploidy != 2:
            args.extend(['-d', str(self.config.ploidy)])

        if self.config.cds_gff:
            args.extend(['-c', self.config.cds_gff])

        if self.config.lncrna_tpm != 1.0:
            args.extend(['--lncrnamintpm', str(self.config.lncrna_tpm)])

        if self.config.partial:
            args.append('--partial')

        if self.config.functional:
            args.append('--functional')

        if self.config.mito_contigs:
            args.extend(['--mito_contigs', self.config.mito_contigs])

        if self.config.extra_gff:
            args.extend(['--extra', self.config.extra_gff])

        if self.config.debug:
            args.append('--debug')

        if self.config.verbose:
            args.append('--verbose')

        # 构建完整命令|Build full command
        eviann_sh = os.path.join(self.config.eviann_path, 'bin', 'eviann.sh')
        cmd = build_conda_command(eviann_sh, args)

        return cmd

    def run_annotation(self):
        """
        运行EviAnn基因组注释|Run EviAnn genome annotation

        Returns:
            bool: 成功返回True，失败返回False|Success returns True, failure returns False
        """
        start_time = datetime.now()

        try:
            self.logger.info("=" * 80)
            self.logger.info("开始EviAnn基因组注释流程|Starting EviAnn Genome Annotation Pipeline")
            self.logger.info("=" * 80)

            # 显示配置信息|Display configuration
            self._print_config()

            # 构建命令|Build command
            cmd = self.build_eviann_command()

            # 记录命令|Log command
            self.logger.info(f"执行|Executing: EviAnn基因组注释|EviAnn genome annotation")
            self.logger.info(f"命令|Command: {' '.join(cmd)}")

            # 执行命令|Execute command
            self.logger.info("开始运行EviAnn|Starting EviAnn...")
            result = subprocess.run(
                cmd,
                check=False,
                cwd=os.path.dirname(self.config.genome)
            )

            end_time = datetime.now()
            duration = int((end_time - start_time).total_seconds())

            if result.returncode == 0:
                self.logger.info("=" * 80)
                self.logger.info("EviAnn基因组注释完成|EviAnn Genome Annotation Completed")
                self.logger.info(f"总耗时|Total duration: {duration} 秒|seconds")
                self.logger.info("=" * 80)

                # 输出结果文件|Output result files
                self._print_output_files()

                return True
            else:
                self.logger.error(f"EviAnn执行失败，返回码|EviAnn execution failed, return code: {result.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"EviAnn流程出错|Error in EviAnn pipeline: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _print_config(self):
        """打印配置信息|Print configuration"""
        self.logger.info("EviAnn配置|EviAnn Configuration:")
        self.logger.info(f"  基因组文件|Genome file: {self.config.genome}")
        self.logger.info(f"  RNA-seq数据|RNA-seq data: {self.config.rnaseq or 'None'}")
        self.logger.info(f"  转录本数据|Transcripts: {self.config.transcripts or 'None'}")
        self.logger.info(f"  蛋白质数据|Proteins: {self.config.proteins or 'None'}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  倍性|Ploidy: {self.config.ploidy}")

    def _print_output_files(self):
        """打印输出文件信息|Print output files"""
        genome_base = Path(self.config.genome).name

        output_files = [
            f"{genome_base}.pseudo_label.gff",
            f"{genome_base}.proteins.fasta",
            f"{genome_base}.transcripts.fasta"
        ]

        self.logger.info("输出文件|Output files:")
        for output_file in output_files:
            output_path = Path(self.config.genome).parent / output_file
            if output_path.exists():
                self.logger.info(f"   {output_file}")
            else:
                self.logger.warning(f"   {output_file} (未找到|not found)")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='EviAnn基因组注释工具|EviAnn Genome Annotation Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-g', '--genome',
                       required=True,
                       help='基因组FASTA文件|Genome FASTA file')

    # 数据输入参数|Data input arguments
    parser.add_argument('--short-reads',
                       help='二代转录组数据（文件或目录）|Short-read RNA-seq data (file or directory)')
    parser.add_argument('--long-reads',
                       help='三代转录组数据（文件或目录）|Long-read RNA-seq data (file or directory)')
    parser.add_argument('-e', '--transcripts',
                       help='转录本FASTA文件|Transcripts FASTA file')
    parser.add_argument('-p', '--proteins',
                       help='蛋白质FASTA文件|Proteins FASTA file')

    # 可选参数|Optional arguments
    parser.add_argument('-s', '--uniprot',
                       help='UniProt-SwissProt FASTA|UniProt-SwissProt FASTA')
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=12,
                       help='线程数|Number of threads')
    parser.add_argument('-m', '--max-intron',
                       type=int,
                       help='最大内含子长度|Maximum intron length')
    parser.add_argument('-d', '--ploidy',
                       type=int,
                       default=2,
                       help='基因组倍性|Genome ploidy')
    parser.add_argument('--lncrna-tpm',
                       type=float,
                       default=1.0,
                       help='lncRNA最小TPM|Minimum TPM for lncRNA')

    # 布尔选项|Boolean options
    parser.add_argument('--partial',
                       action='store_true',
                       help='包含部分CDS|Include partial CDS')
    parser.add_argument('--functional',
                       action='store_true',
                       help='执行功能注释|Perform functional annotation')
    parser.add_argument('--debug',
                       action='store_true',
                       help='调试模式|Debug mode')
    parser.add_argument('--verbose',
                       action='store_true',
                       help='详细输出|Verbose output')

    args = parser.parse_args()

    try:
        # 创建注释器并运行|Create annotator and run
        annotator = EviAnnotator(
            genome=args.genome,
            short_reads=args.short_reads,
            long_reads=args.long_reads,
            transcripts=args.transcripts,
            proteins=args.proteins,
            uniprot=args.uniprot,
            threads=args.threads,
            max_intron=args.max_intron,
            ploidy=args.ploidy,
            lncrna_tpm=args.lncrna_tpm,
            partial=args.partial,
            functional=args.functional,
            debug=args.debug,
            verbose=args.verbose
        )

        success = annotator.run_annotation()
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
