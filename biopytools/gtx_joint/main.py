"""
GTX Joint Calling命令生成器主模块 | GTX Joint Calling Command Generator Main Module
"""

import os
import re
import sys
import time
from pathlib import Path
from .config import GTXJointConfig
from .utils import GTXJointLogger, GTXJointScanner, GTXJointWriter


class GTXJointGenerator:
    """GTX Joint Calling命令生成器 | GTX Joint Calling Command Generator"""

    def __init__(self, **kwargs):
        """
        初始化生成器 | Initialize generator

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        # 初始化配置 | Initialize configuration
        self.config = GTXJointConfig(**kwargs)

        # 初始化日志 | Initialize logging
        log_file = os.path.join(kwargs.get('output_dir', '.'), 'gtx_joint.log') if kwargs.get('output_dir') else None
        self.logger_manager = GTXJointLogger(log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类 | Initialize utility classes
        self.scanner = GTXJointScanner(self.logger)
        self.writer = GTXJointWriter(self.logger)

    def generate_commands(self):
        """
        生成命令脚本 | Generate command script

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            # 验证配置 | Validate configuration
            self.config.validate()
            self.logger.info("=" * 60)
            self.logger.info("🧬 开始生成 GTX Joint Calling 命令 | Starting GTX Joint Calling Command Generation")
            self.logger.info("=" * 60)
            self.logger.info(f"📋 配置信息 | Configuration:\n{self.config}")

            # 检查faketime | Check faketime
            self.config.use_faketime = self.scanner.check_faketime()
            if self.config.use_faketime:
                self.logger.info(f"✅ faketime 可用 | faketime is available")
            else:
                self.logger.warning(f"⚠️ faketime 未找到，将不使用 faketime | faketime not found")

            # 创建输出目录 | Create output directory
            os.makedirs(self.config.output_dir, exist_ok=True)
            self.logger.info(f"✅ 输出目录已创建/存在 | Output directory ready: {self.config.output_dir}")

            # 创建临时目录 | Create temporary directory
            os.makedirs(self.config.tmp_dir, exist_ok=True)
            self.logger.info(f"✅ 临时目录已创建/存在 | Temporary directory ready: {self.config.tmp_dir}")

            # 扫描GVCF文件 | Scan GVCF files
            self.config.gvcf_files = self.scanner.scan_gvcf_files(self.config.gvcf_dir)
            if not self.config.gvcf_files:
                self.logger.error("❌ 未找到GVCF文件 | No GVCF files found")
                return False

            # 读取染色体 | Read chromosomes
            reference_fai = f"{self.config.reference}.fai"
            self.config.chromosomes = self.scanner.read_chromosomes(reference_fai)
            if not self.config.chromosomes:
                self.logger.error("❌ 无法读取染色体信息 | Failed to read chromosome information")
                return False

            # 应用染色体过滤 | Apply chromosome filter
            if self.config.chr_pattern:
                pattern = re.compile(self.config.chr_pattern)
                filtered_chrs = [c for c in self.config.chromosomes if pattern.match(c)]
                if len(filtered_chrs) != len(self.config.chromosomes):
                    self.logger.info(f"🔍 应用过滤模式 | Applying filter pattern: {self.config.chr_pattern}")
                    self.logger.info(f"📊 过滤后 | After filtering: {len(filtered_chrs)}/{len(self.config.chromosomes)} chromosomes")
                    self.config.chromosomes = filtered_chrs

            # 备份旧脚本 | Backup old script
            script_path = os.path.join(self.config.output_dir, self.config.script_name)
            self.writer.backup_old_script(script_path)

            # 生成命令 | Generate commands
            if self.config.window_size:
                success = self._generate_by_windows(script_path)
            else:
                success = self._generate_by_chromosomes(script_path)

            if not success:
                return False

            # 设置可执行权限 | Set executable permission
            self.writer.make_executable(script_path)

            # 输出统计信息 | Output statistics
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("🎉 命令生成完成！| Command generation completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"📊 统计信息 | Statistics:")
            self.logger.info(f"  • 输入 GVCF 文件数 | Input GVCF files: {len(self.config.gvcf_files)}")
            self.logger.info(f"  • 处理染色体数 | Processed chromosomes: {len(self.config.chromosomes)}")
            self.logger.info(f"  • 命令脚本文件 | Script file: {script_path}")
            self.logger.info(f"  • 总耗时 | Total time: {elapsed_time:.2f} seconds")
            self._show_usage_guide(script_path)

            return True

        except Exception as e:
            self.logger.error(f"❌ 程序执行出错 | Program execution error: {str(e)}", exc_info=True)
            return False

    def _generate_by_chromosomes(self, script_path: str) -> bool:
        """
        按染色体生成命令 | Generate commands by chromosomes

        Args:
            script_path: 脚本路径 | Script path

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("📝 正在生成命令脚本... | Generating command script...")
        self.logger.info(f"🔄 运行模式 | Mode: 按染色体拆分 | By chromosome")
        self.logger.info("=" * 60)

        processed_count = 0
        total_chrs = len(self.config.chromosomes)

        # 构建VCF参数 | Build VCF arguments
        vcf_args = " ".join(f"-v {f}" for f in self.config.gvcf_files)
        faketime_prefix = f"faketime '{self.config.faketime}' " if self.config.use_faketime else ""

        for chr_name in self.config.chromosomes:
            output_file = os.path.join(self.config.output_dir, f"{chr_name}.joint.vcf.gz")

            # 生成命令 | Generate command
            command = (
                f"{faketime_prefix}{self.config.gtx_exec} joint "
                f"-r {self.config.reference} "
                f"-o {output_file} "
                f"-L {chr_name} "
                f"--tmp-dir {self.config.tmp_dir} "
                f"-t {self.config.threads} "
                f"{vcf_args}"
            )

            if not self.writer.write_command(script_path, command):
                return False

            processed_count += 1

            # 显示进度 | Show progress
            if processed_count % 10 == 0 or processed_count == total_chrs:
                self.logger.info(f"📊 进度 | Progress: {processed_count}/{total_chrs} 条命令已生成")

        return True

    def _generate_by_windows(self, script_path: str) -> bool:
        """
        按区间生成命令 | Generate commands by windows

        Args:
            script_path: 脚本路径 | Script path

        Returns:
            是否成功 | Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("📝 正在生成命令脚本... | Generating command script...")
        self.logger.info(f"🔄 运行模式 | Mode: 按 {self.config.window_size} bp 区间拆分 | By windows")
        self.logger.info("=" * 60)

        # 读取染色体长度 | Read chromosome lengths
        chr_lengths = self._read_chromosome_lengths(f"{self.config.reference}.fai")

        processed_count = 0
        vcf_args = " ".join(f"-v {f}" for f in self.config.gvcf_files)
        faketime_prefix = f"faketime '{self.config.faketime}' " if self.config.use_faketime else ""

        for chr_name in self.config.chromosomes:
            if chr_name not in chr_lengths:
                self.logger.warning(f"⚠️ 无法获取 {chr_name} 的长度，跳过 | Cannot get length for {chr_name}, skipping")
                continue

            chr_len = chr_lengths[chr_name]
            num_windows = (chr_len + self.config.window_size - 1) // self.config.window_size

            for i in range(num_windows):
                start = i * self.config.window_size + 1
                end = min((i + 1) * self.config.window_size, chr_len)

                output_file = os.path.join(
                    self.config.output_dir,
                    f"{chr_name}_{start}-{end}.joint.vcf.gz"
                )
                region = f"{chr_name}:{start}-{end}"

                command = (
                    f"{faketime_prefix}{self.config.gtx_exec} joint "
                    f"-r {self.config.reference} "
                    f"-o {output_file} "
                    f"-L {region} "
                    f"--tmp-dir {self.config.tmp_dir} "
                    f"-t {self.config.threads} "
                    f"{vcf_args}"
                )

                if not self.writer.write_command(script_path, command):
                    return False

                processed_count += 1

                if processed_count % 50 == 0:
                    self.logger.info(f"📊 进度 | Progress: {processed_count} 条命令已生成")

        return True

    def _read_chromosome_lengths(self, fai_file: str) -> dict:
        """
        读取染色体长度 | Read chromosome lengths

        Args:
            fai_file: 参考基因组索引文件 | Reference genome index file

        Returns:
            染色体长度字典 | Dictionary of chromosome lengths
        """
        chr_lengths = {}
        try:
            with open(fai_file, 'r') as f:
                for line in f:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        chr_lengths[parts[0]] = int(parts[1])
        except Exception as e:
            self.logger.error(f"❌ 读取染色体长度失败 | Failed to read chromosome lengths: {e}")
        return chr_lengths

    def _show_usage_guide(self, script_path: str):
        """显示使用指南 | Show usage guide"""
        self.logger.info("=" * 60)
        self.logger.info("💡 执行方式 | Execution methods:")
        self.logger.info("=" * 60)
        self.logger.info("1️⃣ 串行执行 | Serial execution:")
        self.logger.info(f"   bash {script_path}")
        self.logger.info("")
        self.logger.info("2️⃣ GNU Parallel 并行 | Parallel execution:")
        self.logger.info(f"   parallel -j 4 < {script_path}")
        self.logger.info("")

        if self.config.window_size:
            self.logger.info("📌 区间拆分模式提示 | Window mode tip:")
            self.logger.info("   运行完成后需要合并同一染色体的所有区间VCF文件")
            self.logger.info("   推荐使用 bcftools concat 进行合并 | Use bcftools concat to merge:")
            self.logger.info(f"   bcftools concat -o Chr01.merged.vcf.gz {self.config.output_dir}/Chr01_*.joint.vcf.gz")
            self.logger.info("")

        self.logger.info("=" * 60)


def main():
    """主函数 | Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='🧬 GTX Joint Calling命令生成工具 | GTX Joint Calling Command Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 🚀 基本用法 | Basic usage
  %(prog)s -g /path/to/gtx -r genome.fa -i ./gvcf -o ./output

  # 🔄 按区间拆分 | By windows (10M)
  %(prog)s -g ./gtx -r genome.fa -i ./gvcf -o ./output -w 10000000

  # 🔍 只处理主染色体 | Only main chromosomes
  %(prog)s -g ./gtx -r genome.fa -i ./gvcf -o ./output -p "^Chr[0-9]+$"

  # ⚡ 完整参数 | Full parameters
  %(prog)s -g ./gtx -r genome.fa -i ./gvcf -o ./output \\
      -t 24 -T /tmp -s my_script.sh -p "^Chr" -w 10000000
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments | 必需参数')
    required.add_argument('-g', '--gtx', default='/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx',
                         help='📂 GTX可执行文件路径 (默认: /share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx) | GTX executable path')
    required.add_argument('-r', '--ref', required=True,
                         help='🧬 参考基因组文件路径 | Reference genome file path')
    required.add_argument('-i', '--input', required=True,
                         help='📂 GVCF文件所在目录 | GVCF files directory')
    required.add_argument('-o', '--output', required=True,
                         help='📤 输出结果目录 | Output directory')

    # 可选参数 | Optional parameters
    optional = parser.add_argument_group('optional arguments | 可选参数')
    optional.add_argument('-t', '--threads', type=int, default=88,
                         help='🧵 线程数 | Number of threads (default: 88)')
    optional.add_argument('-T', '--tmp-dir', default='./tmp',
                         help='🗑️ 临时目录 | Temporary directory (default: ./tmp)')
    optional.add_argument('-s', '--script', default='./run_gtx_joint.sh',
                         help='📄 输出脚本文件名 | Output script filename (default: run_gtx_joint.sh)')
    optional.add_argument('-f', '--faketime', default='2020-10-20 00:00:00',
                         help='⏰ faketime时间 | faketime time (default: 2020-10-20 00:00:00)')
    optional.add_argument('-p', '--pattern',
                         help='🔍 染色体过滤正则 | Chromosome filter pattern')
    optional.add_argument('-w', '--window', type=int,
                         help='📏 区间大小(bp) | Window size in bp (e.g., 10000000 for 10M)')

    args = parser.parse_args()

    # 创建生成器并运行 | Create generator and run
    generator = GTXJointGenerator(
        gtx_exec=args.gtx,
        reference=args.ref,
        gvcf_dir=args.input,
        output_dir=args.output,
        threads=args.threads,
        tmp_dir=args.tmp_dir,
        script_name=args.script,
        faketime=args.faketime,
        chr_pattern=args.pattern,
        window_size=args.window
    )

    success = generator.generate_commands()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
