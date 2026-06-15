"""
EGAPx批量生成器主模块|EGAPx Batch Generator Main Module
"""

import os
import sys
import time
import shutil
from .config import EGAPxBatchConfig
from .utils import EGAPxBatchLogger, GenomeSplitter, TemplateProcessor, JobGenerator


class EGAPxBatchGenerator:
    """EGAPx批量运行生成器|EGAPx Batch Processing Generator"""

    def __init__(self, **kwargs):
        """
        初始化生成器|Initialize generator

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = EGAPxBatchConfig(**kwargs)
        self.config.validate()

        # 创建输出目录|Create output directory
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 初始化日志|Initialize logging
        self.logger_manager = EGAPxBatchLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类|Initialize utility classes
        self.splitter = GenomeSplitter(self.logger)
        self.template_processor = TemplateProcessor(self.logger)
        self.job_generator = JobGenerator(self.logger)

        # 处理reads输入|Process reads inputs
        self.processed_short_reads = self._process_reads_input(
            self.config.short_reads, "short"
        ) if self.config.short_reads else ""
        self.processed_long_reads = self._process_reads_input(
            self.config.long_reads, "long"
        ) if self.config.long_reads else ""

    def _extract_sample_name(self, filename: str) -> str:
        """
        从文件名提取样本名，识别双端测序模式|Extract sample name from filename,识别paired-end patterns

        Args:
            filename: 文件名|Filename

        Returns:
            str: 样本名|Sample name
        """
        import re

        # 去除扩展名|Remove extensions
        name = filename.replace('.fq.gz', '').replace('.fastq.gz', '') \
                       .replace('.fq', '').replace('.fastq', '')

        # 识别双端测序模式并去除后缀|Recognize paired-end patterns and remove suffix
        # 匹配模式: _1, _2, _R1, _R2, .1, .2, .R1, .R2 (可能包含.clean等后缀)
        # Match patterns: _1, _2, _R1, _R2, .1, .2, .R1, .R2 (may include .clean suffix)

        # 尝试各种模式|Try various patterns
        patterns = [
            r'_(\d+)(?:\.clean)?$',  # _1.clean, _2.clean
            r'_(R[12])(?:\.clean)?$',  # _R1.clean, _R2.clean
            r'\.(\d+)(?:\.clean)?$',  # .1.clean, .2.clean
            r'\.(R[12])(?:\.clean)?$',  # .R1.clean, .R2.clean
            r'_(\d+)$',  # _1, _2
            r'_(R[12])$',  # _R1, _R2
            r'\.(\d+)$',  # .1, .2
            r'\.(R[12])$',  # .R1, .R2
        ]

        for pattern in patterns:
            match = re.search(pattern, name)
            if match:
                # 去除匹配的后缀|Remove matched suffix
                name = re.sub(pattern, '', name)
                break

        return name

    def _process_reads_input(self, input_path: str, reads_type: str) -> str:
        """
        处理reads输入(目录或文件)|Process reads input (directory or file)

        Args:
            input_path: 输入路径(文件或目录)|Input path (file or directory)
            reads_type: reads类型(short/long)|Reads type (short/long)

        Returns:
            str: EGAPx格式的文件路径（绝对路径）|Path to EGAPx format file (absolute path)
        """
        if not input_path:
            return ""

        # 生成输出文件路径（绝对路径）|Generate output file path (absolute path)
        output_file = os.path.abspath(os.path.join(
            self.config.output_dir,
            f"{reads_type}_reads_list.txt"
        ))

        # 如果是文件，生成单文件列表|If file, generate single-file list
        if os.path.isfile(input_path):
            self.logger.info(f"  {reads_type}_reads: 使用文件|Using file: {input_path}")

            # 获取绝对路径|Get absolute path
            abs_path = os.path.abspath(input_path)

            # 从文件名提取样本名|Extract sample name from filename
            basename = os.path.basename(abs_path)
            sample_name = self._extract_sample_name(basename)

            # 写入EGAPx格式|Write EGAPx format
            with open(output_file, 'w') as f:
                f.write(f"{sample_name}\t{abs_path}\n")

            self.logger.info(f"  生成reads列表文件|Generated reads list file: {output_file}")
            return output_file

        # 如果是目录，生成EGAPx格式文件|If directory, generate EGAPx format file
        if os.path.isdir(input_path):
            self.logger.info(f"  {reads_type}_reads: 扫描目录|Scanning directory: {input_path}")

            # 扫描目录中的fq/fastq文件|Scan for fq/fastq files in directory
            reads_files = []
            for item in os.listdir(input_path):
                if item.endswith(('.fq', '.fq.gz', '.fastq', '.fastq.gz')):
                    reads_files.append(os.path.join(input_path, item))

            if not reads_files:
                self.logger.warning(f"  目录中未找到fq/fastq文件|No fq/fastq files found in directory")
                return ""

            reads_files.sort()

            # 按样本名分组|Group by sample name
            from collections import defaultdict
            sample_groups = defaultdict(list)

            for reads_file in reads_files:
                # 获取绝对路径|Get absolute path
                abs_path = os.path.abspath(reads_file)

                # 从文件名提取样本名|Extract sample name from filename
                basename = os.path.basename(abs_path)
                sample_name = self._extract_sample_name(basename)

                sample_groups[sample_name].append(abs_path)

            # 生成EGAPx格式文件（使用绝对路径）|Generate EGAPx format file (using absolute paths)
            with open(output_file, 'w') as f:
                for sample_name in sorted(sample_groups.keys()):
                    # 同一个样本的多个文件按顺序写入|Write multiple files of same sample in order
                    for file_path in sorted(sample_groups[sample_name]):
                        f.write(f"{sample_name}\t{file_path}\n")

            self.logger.info(f"  生成reads列表文件|Generated reads list file: {output_file}")
            self.logger.info(f"  包含{len(reads_files)}个文件，{len(sample_groups)}个样本|Contains {len(reads_files)} files, {len(sample_groups)} samples")

            return output_file

        return ""

    def generate(self):
        """
        执行批量生成|Perform batch generation

        Returns:
            是否成功|Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("开始EGAPx批量运行配置生成|Starting EGAPx Batch Config Generation")
            self.logger.info("=" * 60)
            self.logger.info(f"配置信息|Configuration:\n{self.config}")

            # 获取输出目录的绝对路径|Get absolute path of output directory
            output_dir_abs = os.path.abspath(self.config.output_dir)

            # 构建序列列表|Build sequence list
            if self.config.split_genome:
                # 步骤1: 拆分基因组|Step 1: Split genome
                self.logger.info("=" * 60)
                self.logger.info("步骤1/3: 拆分基因组文件|Step 1/3: Split genome file")
                self.logger.info("=" * 60)

                split_dir = os.path.join(self.config.output_dir, "temp_split")

                split_files = self.splitter.split_genome(
                    self.config.genome,
                    split_dir,
                    self.config.chr_prefix
                )

                if not split_files:
                    self.logger.error("未找到任何序列|No sequences found")
                    return False

                self.logger.info(f"成功拆分 {len(split_files)} 个序列|Successfully split {len(split_files)} sequences")
            else:
                # 不拆分，直接使用原始基因组文件|No split, use original genome directly
                split_dir = None
                genome_name = os.path.basename(self.config.genome)
                if genome_name.endswith('.fa.gz') or genome_name.endswith('.fasta.gz'):
                    genome_name = genome_name.replace('.fa.gz', '').replace('.fasta.gz', '')
                else:
                    genome_name = genome_name.replace('.fa', '').replace('.fasta', '')
                import subprocess
                chr_size = int(subprocess.check_output(
                    ['awk', '/^>/{next} !/^>/{n+=length} END{print n}', self.config.genome]
                ).strip())
                split_files = [(self.config.genome, chr_size)]
                self.logger.info(f"跳过拆分|Skipped splitting, using whole genome: {genome_name}")

            # 步骤2: 生成配置文件和脚本|Step 2: Generate configs and scripts
            self.logger.info("=" * 60)
            self.logger.info("步骤2/3: 生成配置文件和脚本|Step 2/3: Generate configs and scripts")
            self.logger.info("=" * 60)

            # 任务列表放在输出目录的上一层|Job list in parent directory of output
            job_list_path = os.path.join(os.path.dirname(self.config.output_dir), "all_jobs_submit.list.sh")

            with open(job_list_path, 'w') as job_list:
                for chr_file, chr_size in split_files:
                    chr_name = os.path.basename(chr_file).replace('.fa', '')
                    chr_dir = os.path.join(self.config.output_dir, chr_name)
                    chr_dir_abs = os.path.join(output_dir_abs, chr_name)

                    # 创建目录结构|Create directory structure
                    os.makedirs(chr_dir, exist_ok=True)
                    os.makedirs(os.path.join(chr_dir, 'work'), exist_ok=True)
                    os.makedirs(os.path.join(chr_dir, 'output'), exist_ok=True)

                    # 处理基因组FASTA文件|Handle genome FASTA file
                    chr_fa = os.path.join(chr_dir_abs, f"{chr_name}.fa")
                    chr_fa_relative = os.path.join(chr_dir, f"{chr_name}.fa")
                    if self.config.split_genome:
                        shutil.move(chr_file, chr_fa_relative)
                    else:
                        if chr_file != chr_fa_relative:
                            shutil.copy2(chr_file, chr_fa_relative)

                    # 生成YAML配置|Generate YAML config
                    new_yaml = os.path.join(chr_dir_abs, f"{chr_name}.yaml")
                    # 构建locus标签前缀|Build locus tag prefix
                    if self.config.locus_tag_prefix:
                        new_prefix = f"{self.config.locus_tag_prefix}_{chr_name}"
                    else:
                        new_prefix = chr_name

                    yaml_content = self.job_generator.generate_yaml(
                        genome=chr_fa,
                        locus_tag_prefix=new_prefix,
                        taxid=self.config.taxid,
                        short_reads=self.processed_short_reads,
                        long_reads=self.processed_long_reads,
                        chr_size=chr_size
                    )

                    with open(new_yaml, 'w') as f:
                        f.write(yaml_content)

                    # 生成运行脚本|Generate run script
                    new_script = os.path.join(chr_dir_abs, f"egapx_{chr_name}.sh")

                    script_content = self.job_generator.generate_script(
                        yaml_path=new_yaml,
                        work_dir=os.path.join(chr_dir_abs, 'work'),
                        output_dir=os.path.join(chr_dir_abs, 'output'),
                        report_name=self.config.report_name,
                        chr_dir=chr_dir_abs,
                        egapx_dir=os.path.dirname(output_dir_abs),
                        local_cache=self.config.local_cache
                    )

                    with open(new_script, 'w') as f:
                        f.write(script_content)

                    # 设置可执行权限|Set executable permission
                    os.chmod(new_script, 0o755)

                    # 添加到任务列表|Add to job list
                    job_list.write(f"bash \"{new_script}\"\n")

                    self.logger.info(f"  {chr_name} 配置完成|{chr_name} config completed")

            # 步骤3: 完成统计|Step 3: Completion statistics
            self.logger.info("=" * 60)
            self.logger.info("步骤3/3: 清理临时目录|Step 3/3: Cleanup temporary directory")
            self.logger.info("=" * 60)

            # 清理临时目录|Cleanup temporary directory
            if split_dir and os.path.exists(split_dir):
                shutil.rmtree(split_dir)

            # 完成统计|Completion statistics
            elapsed_time = time.time() - start_time
            self.config.job_count = len(split_files)

            self.logger.info("=" * 60)
            self.logger.info("配置生成完成|Config generation completed")
            self.logger.info("=" * 60)
            self.logger.info(f"统计信息|Statistics:")
            self.logger.info(f"  - 任务数量|Job count: {self.config.job_count}")
            self.logger.info(f"  - 任务列表|Job list: {job_list_path}")
            self.logger.info(f"  - 总耗时|Total time: {elapsed_time:.2f} seconds")

            # 在输出目录的上一层创建EGAPx软链接|Create EGAPx symlinks in parent directory
            self._create_egapx_symlinks_in_parent()

            self._show_usage_guide(job_list_path)

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}", exc_info=True)
            return False

    def _create_egapx_symlinks_in_parent(self):
        """在输出目录的上一层创建EGAPx运行所需的软链接|Create EGAPx runtime symlinks in parent directory"""
        parent_dir = os.path.dirname(self.config.output_dir)

        # 仅链接运行所需目录|Only link directories required for execution
        required_dirs = ['ui', 'nf', 'egapx_config']

        try:
            for item in required_dirs:
                src = os.path.join(self.config.egapx_path, item)
                dst = os.path.join(parent_dir, item)
                if os.path.exists(src) and not os.path.exists(dst):
                    os.symlink(src, dst)
                    self.logger.info(f"创建软链接|Created symlink: {item}")

            # 自定义SIF镜像|Custom SIF image
            if self.config.sif_image:
                config_dir = os.path.join(parent_dir, 'egapx_config')
                if os.path.islink(config_dir):
                    os.remove(config_dir)
                os.makedirs(config_dir, exist_ok=True)
                sif_config = os.path.join(config_dir, 'singularity.config')
                with open(sif_config, 'w') as f:
                    f.write(f"singularity.enabled = true\n")
                    f.write(f"process.container = '{self.config.sif_image}'\n")
                self.logger.info(f"自定义SIF镜像配置|Custom SIF config: {self.config.sif_image}")

            self.logger.info(f"EGAPx软链接已创建|EGAPx symlinks created in: {parent_dir}")

        except Exception as e:
            self.logger.warning(f"创建软链接失败|Failed to create symlinks: {e}")

    def _show_usage_guide(self, job_list_path):
        """显示使用指南|Show usage guide"""
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("执行方式|Execution methods:")
        self.logger.info("=" * 60)
        self.logger.info("顺序执行|Serial execution:")
        self.logger.info(f"  bash {job_list_path}")
        self.logger.info("")
        self.logger.info("使用 GNU parallel|Using GNU parallel:")
        self.logger.info(f"  cat {job_list_path}|parallel -j 4")
        self.logger.info("")
        self.logger.info("使用 xargs|Using xargs:")
        self.logger.info(f"  cat {job_list_path}|xargs -P 4 -I {{}} bash -c {{}}")
        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='EGAPx批量运行配置生成工具|EGAPx Batch Config Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g genome.fa -o output_dir
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('required arguments|必需参数')
    required.add_argument('-g', '--genome', required=True,
                         help='[FILE] 基因组FASTA文件路径|Genome FASTA file path')
    required.add_argument('-o', '--output', required=True,
                         help='[DIR] 输出目录路径|Output directory path')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('optional arguments|可选参数')
    optional.add_argument('-e', '--egapx',
                         default='~/software/EGAPX_v.0.4.1-alpha/egapx',
                         help='[PATH] EGAPx安装路径|EGAPx installation path')
    optional.add_argument('--local-cache',
                         default='~/software/EGAPX_v.0.4.1-alpha/local_cache',
                         help='[PATH] EGAPx本地缓存路径|EGAPx local cache path')
    optional.add_argument('--sif', default='~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif',
                         help='[FILE] Singularity镜像路径|Singularity image path')
    optional.add_argument('--no-split', action='store_true', default=False,
                         help='不按染色体拆分基因组|Do not split genome by chromosome')
    optional.add_argument('--taxid', default='71234',
                         help='[INT] 物种分类ID|Species taxonomy ID')
    optional.add_argument('-p', '--chr-prefix',
                         help='[STR] 染色体前缀过滤|Chromosome prefix filter')
    optional.add_argument('--locus-prefix', default='',
                         help='[STR] locus标签前缀|Locus tag prefix')
    optional.add_argument('--report-name', default='EGAPx',
                         help='[STR] 报告名称|Report name')
    optional.add_argument('--short-reads', default='',
                         help='[FILE] 短读长测序数据文件路径|Short reads file path')
    optional.add_argument('--long-reads', default='',
                         help='[FILE] 长读长测序数据文件路径|Long reads file path')

    args = parser.parse_args()

    # 创建生成器并运行|Create generator and run
    generator = EGAPxBatchGenerator(
        genome=args.genome,
        output_dir=args.output,
        taxid=args.taxid,
        egapx_path=args.egapx,
        local_cache=args.local_cache,
        sif_image=args.sif,
        split_genome=not args.no_split,
        chr_prefix=args.chr_prefix,
        locus_tag_prefix=args.locus_prefix,
        report_name=args.report_name,
        short_reads=args.short_reads,
        long_reads=args.long_reads
    )

    success = generator.generate()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
