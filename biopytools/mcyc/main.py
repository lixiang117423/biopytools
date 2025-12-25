"""
甲烷循环基因丰度分析主程序模块 | Methane Cycle Gene Abundance Analysis Main Module
"""

import os
import sys
import shutil
import subprocess
import argparse
import re
import time
import pandas as pd
from typing import Optional

from .config import MCycConfig
from .calculator import MatrixCalculator


class MCycAnalyzer:
    """甲烷循环分析主类 | Main Methane Cycle Analysis Class"""

    def __init__(self, config: MCycConfig):
        """
        初始化分析器 | Initialize analyzer

        Args:
            config: MCycConfig配置对象 | MCycConfig configuration object
        """
        self.config = config
        self.calculator = MatrixCalculator(config)

    def print_log(self, msg: str, level: str = "INFO"):
        """打印日志 | Print log"""
        icon_map = {
            "INFO": "ℹ️ ", "WARN": "⚠️ ", "ERROR": "❌ ",
            "SUCCESS": "✅ ", "PROCESS": "⏳ ", "TOOL": "🛠️ "
        }
        icon = icon_map.get(level, "ℹ️ ")
        print(f"{icon} [{time.strftime('%H:%M:%S')}] {msg}")

    def check_dependency(self, tool: str) -> bool:
        """
        检查依赖工具是否可用 | Check if dependency tool is available

        Args:
            tool: 工具名称 | Tool name

        Returns:
            bool: 是否可用 | Whether available
        """
        if not shutil.which(tool):
            self.print_log(f"未找到命令 | Tool not found: {tool}，请先安装或加载环境 | Please install or load environment first", "ERROR")
            return False
        return True

    def run_cmd(self, cmd: str) -> bool:
        """
        执行命令 | Execute command

        Args:
            cmd: 命令字符串 | Command string

        Returns:
            bool: 是否执行成功 | Whether execution was successful
        """
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            return True
        except subprocess.CalledProcessError:
            self.print_log(f"命令执行失败 | Command execution failed: {cmd}", "ERROR")
            return False

    def step1_check_dependencies(self) -> bool:
        """
        步骤1: 检查依赖工具 | Step 1: Check dependencies

        Returns:
            bool: 检查是否通过 | Whether check passed
        """
        self.print_log("正在检查依赖工具... | Checking dependencies...", "PROCESS")

        required_tools = ['diamond', 'seqkit', 'perl']
        for tool in required_tools:
            if not self.check_dependency(tool):
                return False

        if not os.path.exists(self.config.input_list):
            self.print_log(f"输入文件不存在 | Input file not found: {self.config.input_list}", "ERROR")
            return False

        self.print_log("依赖检查完成 | Dependency check completed", "SUCCESS")
        return True

    def step2_prepare_database(self) -> bool:
        """
        步骤2: 准备数据库 | Step 2: Prepare database

        Returns:
            bool: 是否准备成功 | Whether preparation was successful
        """
        self.print_log("准备数据库文件... | Preparing database files...", "PROCESS")

        # 检查并构建Diamond索引 | Check and build Diamond index
        if not os.path.exists(self.config.dmnd_file):
            self.print_log("未找到Diamond索引，正在构建... | Diamond index not found, building...", "WARN")
            cmd = f"diamond makedb --in {self.config.fasta_file} --db {self.config.diamond_db_base} --quiet"
            if not self.run_cmd(cmd):
                self.print_log("Diamond索引构建失败 | Diamond index building failed", "ERROR")
                return False

        # 建立数据库软链接 | Create database symlinks
        self.print_log("建立数据库软链接... | Creating database symlinks...", "TOOL")
        try:
            link_files = ["MCycDB_2021.dmnd", "id2gene.map", "MCycDB_2021.faa"]
            for link in link_files:
                if os.path.exists(link):
                    os.remove(link)

            os.symlink(self.config.dmnd_file, "MCycDB_2021.dmnd")
            os.symlink(self.config.map_file, "id2gene.map")
            os.symlink(self.config.fasta_file, "MCycDB_2021.faa")

        except OSError as e:
            self.print_log(f"创建软链接失败 | Failed to create symlinks: {e}", "ERROR")
            return False

        self.print_log("数据库准备完成 | Database preparation completed", "SUCCESS")
        return True

    def step3_prepare_data(self) -> bool:
        """
        步骤3: 准备数据文件 | Step 3: Prepare data files

        Returns:
            bool: 是否准备成功 | Whether preparation was successful
        """
        self.print_log("正在处理Fastq文件... | Processing Fastq files...", "PROCESS")

        # 创建临时目录 | Create temporary directory
        if not os.path.exists(self.config.staging_dir):
            os.makedirs(self.config.staging_dir)

        # 处理输入文件列表 | Process input file list
        processed_samples = []
        with open(self.config.input_list, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 2:
                    self.print_log(f"第{line_num}行格式错误：需要至少两列 | Line {line_num} format error: need at least 2 columns", "WARN")
                    continue

                sample_name = parts[0]
                paths = parts[1]
                target = os.path.join(self.config.staging_dir, f"{sample_name}.fastq.gz")

                if not os.path.exists(target):
                    if ";" in paths:
                        # 处理双端测序 | Handle paired-end sequencing
                        r1, r2 = paths.split(";")
                        cmd = f"cat {r1} {r2} > {target}"
                        if not self.run_cmd(cmd):
                            self.print_log(f"合并文件失败 | Failed to merge files for sample: {sample_name}", "ERROR")
                            return False
                    else:
                        # 创建软链接 | Create symlink
                        try:
                            os.symlink(os.path.abspath(paths), target)
                        except OSError as e:
                            self.print_log(f"创建软链接失败 | Failed to create symlink for sample {sample_name}: {e}", "ERROR")
                            return False

                processed_samples.append(sample_name)

        if not processed_samples:
            self.print_log("没有找到有效的样本 | No valid samples found", "ERROR")
            return False

        self.print_log(f"成功准备 {len(processed_samples)} 个样本 | Successfully prepared {len(processed_samples)} samples", "SUCCESS")
        return True

    def step4_count_reads(self) -> bool:
        """
        步骤4: 统计Reads数量 | Step 4: Count reads

        Returns:
            bool: 是否统计成功 | Whether counting was successful
        """
        self.print_log("正在统计Reads数量... | Counting reads...", "PROCESS")

        # 使用seqkit统计 | Use seqkit to count
        cmd = f"seqkit stats -j {self.config.thread_count} {self.config.staging_dir}/*.fastq.gz -T"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)

        if result.returncode != 0:
            self.print_log(f"SeqKit统计失败 | SeqKit counting failed", "ERROR")
            return False

        # 解析结果 | Parse results
        try:
            lines = result.stdout.strip().split('\n')
            header = lines[0].split('\t')
            file_idx = header.index("file")
            num_seqs_idx = header.index("num_seqs")

            with open(self.config.sample_info_path, 'w') as f:
                f.write("Sample\tReadCount\n")  # 添加表头 | Add header
                for line in lines[1:]:
                    cols = line.split('\t')
                    if len(cols) > max(file_idx, num_seqs_idx):
                        file_path = cols[file_idx]
                        count = cols[num_seqs_idx].replace(",", "")
                        basename = os.path.basename(file_path)
                        sample_name = basename.replace(".fastq.gz", "")
                        f.write(f"{sample_name}\t{count}\n")

            self.print_log("Reads统计完成 | Read counting completed", "SUCCESS")
            return True

        except (ValueError, IndexError) as e:
            self.print_log(f"SeqKit输出解析失败 | SeqKit output parsing failed: {e}", "ERROR")
            return False

    def step5_run_diamond(self) -> bool:
        """
        步骤5: 运行Diamond比对 | Step 5: Run Diamond alignment

        Returns:
            bool: 是否比对成功 | Whether alignment was successful
        """
        self.print_log("配置并运行MCycDB Perl脚本... | Configuring and running MCycDB Perl script...", "PROCESS")

        # 复制并配置Perl脚本 | Copy and configure Perl script
        try:
            shutil.copy(self.config.perl_script_src, self.config.perl_script_dst)

            # 修改Perl脚本配置 | Modify Perl script configuration
            with open(self.config.perl_script_dst, 'r') as f:
                content = f.read()

            content = re.sub(r'my \$db_file = .*?;', 'my $db_file = "MCycDB_2021.faa";', content)
            content = re.sub(r'my \$diamond_db = .*?;', 'my $diamond_db = "MCycDB_2021.dmnd";', content)
            content = re.sub(r'my \$map_file = .*?;', 'my $map_file = "id2gene.map";', content)

            with open(self.config.perl_script_dst, 'w') as f:
                f.write(content)

        except (IOError, OSError) as e:
            self.print_log(f"Perl脚本配置失败 | Perl script configuration failed: {e}", "ERROR")
            return False

        # 检查是否已有结果 | Check if results already exist
        if not self.config.skip_diamond:
            if os.path.exists(self.config.raw_output) and os.path.getsize(self.config.raw_output) > 0:
                self.print_log("检测到已有中间结果，跳过Diamond比对... | Found intermediate results, skipping Diamond alignment...", "WARN")
                return True

            self.print_log("开始运行MCycDB分析 (Diamond)... | Starting MCycDB analysis (Diamond)...", "PROCESS")
            cmd = f"perl {self.config.perl_script_dst} -d {self.config.staging_dir} -m diamond -f fastq.gz -s nucl -si {self.config.sample_info_path} -o {self.config.raw_output}"

            if not self.run_cmd(cmd):
                self.print_log("Diamond比对失败 | Diamond alignment failed", "ERROR")
                return False

            self.print_log("Diamond比对完成 | Diamond alignment completed", "SUCCESS")
        else:
            self.print_log("跳过Diamond比对（配置要求）| Skipping Diamond alignment (configuration requirement)", "INFO")

        return True

    def step6_calculate_matrices(self) -> bool:
        """
        步骤6: 计算丰度矩阵 | Step 6: Calculate abundance matrices

        Returns:
            bool: 是否计算成功 | Whether calculation was successful
        """
        self.print_log("正在进行矩阵计算... | Performing matrix calculations...", "PROCESS")
        return self.calculator.run_all_calculations()

    def step7_cleanup(self) -> bool:
        """
        步骤7: 清理临时文件 | Step 7: Clean up temporary files

        Returns:
            bool: 是否清理成功 | Whether cleanup was successful
        """
        if self.config.keep_temp_files:
            self.print_log("保留临时文件（配置要求）| Keeping temporary files (configuration requirement)", "INFO")
            return True

        self.print_log("清理临时文件... | Cleaning up temporary files...", "TOOL")

        # 删除软链接 | Remove symlinks
        link_files = ["MCycDB_2021.dmnd", "id2gene.map", "MCycDB_2021.faa"]
        for link in link_files:
            if os.path.exists(link):
                try:
                    os.remove(link)
                except OSError:
                    pass

        # 删除临时目录 | Remove temporary directory
        if os.path.exists(self.config.staging_dir):
            try:
                shutil.rmtree(self.config.staging_dir)
            except OSError:
                pass

        # 删除临时文件 | Remove temporary files
        temp_files = [self.config.sample_info_path, self.config.raw_output, self.config.perl_script_dst]
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except OSError:
                    pass

        self.print_log("临时文件清理完成 | Temporary files cleanup completed", "SUCCESS")
        return True

    def run_full_analysis(self) -> bool:
        """
        运行完整的分析流程 | Run complete analysis pipeline

        Returns:
            bool: 是否成功完成 | Whether completed successfully
        """
        print("========================================================")
        self.print_log("MCycDB分析流程启动 | MCycDB analysis pipeline started")
        print("========================================================")

        # 验证配置 | Validate configuration
        try:
            self.config.validate()
        except (ValueError, FileNotFoundError) as e:
            self.print_log(f"配置验证失败 | Configuration validation failed: {e}", "ERROR")
            return False

        # 执行分析步骤 | Execute analysis steps
        steps = [
            (self.step1_check_dependencies, "检查依赖 | Checking dependencies"),
            (self.step2_prepare_database, "准备数据库 | Preparing database"),
            (self.step3_prepare_data, "准备数据 | Preparing data"),
            (self.step4_count_reads, "统计Reads | Counting reads"),
            (self.step5_run_diamond, "Diamond比对 | Diamond alignment"),
            (self.step6_calculate_matrices, "计算矩阵 | Calculating matrices"),
            (self.step7_cleanup, "清理文件 | Cleaning files")
        ]

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.print_log(f"执行步骤 {i} | Executing step {i}: {step_name}", "INFO")
            if not step_func():
                self.print_log(f"步骤 {i} 失败 | Step {i} failed", "ERROR")
                return False
            self.print_log(f"步骤 {i} 完成 | Step {i} completed", "SUCCESS")

        # 输出最终结果 | Output final results
        self.print_output_summary()

        print("========================================================")
        self.print_log("分析全部完成！| Analysis completed successfully!", "SUCCESS")
        return True

    def print_output_summary(self):
        """输出结果摘要 | Output results summary"""
        print("\n📂 输出文件 | Output files:")

        # 获取矩阵信息 | Get matrix information
        gene_count, sample_count = self.calculator.get_matrix_info()
        if gene_count > 0 and sample_count > 0:
            print(f"   矩阵维度 | Matrix dimension: {gene_count} 基因 × {sample_count} 样本")

        # 获取输出文件状态 | Get output file status
        outputs = self.calculator.get_output_summary()
        output_names = {
            "raw_counts": "Matrix_1_RawCounts.txt (原始计数 | Raw Counts)",
            "tpm": "Matrix_2_TPM.txt       (相对丰度 | Relative Abundance)",
            "clr": "Matrix_3_CLR.txt       (CLR变换 | CLR Transformed)"
        }

        for key, name in output_names.items():
            if outputs[key]["exists"]:
                print(f"   ✅ {name}")


def main():
    """命令行入口函数 | Command line entry function"""
    parser = argparse.ArgumentParser(
        description="MCycDB 甲烷循环基因丰度分析工具 | MCycDB Methane Cycle Gene Abundance Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:

  # 基本分析 | Basic analysis
  python -m biopytools.mcyc samples.txt

  # 指定输出目录 | Specify output directory
  python -m biopytools.mcyc samples.txt -o ./results

  # 跳过Diamond比对（如果已有结果）| Skip Diamond alignment (if results already exist)
  python -m biopytools.mcyc samples.txt --skip-diamond

  # 保留临时文件 | Keep temporary files
  python -m biopytools.mcyc samples.txt --keep-temp

  # 高性能分析 | High performance analysis
  python -m biopytools.mcyc samples.txt -t 8
        """
    )

    parser.add_argument('input_list', help='输入文件列表路径 | Input file list path (sample_name\\tfile_path)')
    parser.add_argument('-o', '--output', help='输出目录路径 | Output directory path (default: current directory)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='线程数 | Thread count (default: 4)')
    parser.add_argument('--mcyc-base', help='MCycDB基础路径 | MCycDB base path')
    parser.add_argument('--skip-diamond', action='store_true', help='跳过Diamond比对 | Skip Diamond alignment')
    parser.add_argument('--keep-temp', action='store_true', help='保留临时文件 | Keep temporary files')

    args = parser.parse_args()

    # 创建配置 | Create configuration
    try:
        config = MCycConfig(
            input_list=args.input_list,
            output_dir=args.output,
            mcyc_base=args.mcyc_base,
            thread_count=args.threads,
            skip_diamond=args.skip_diamond,
            keep_temp_files=args.keep_temp
        )
    except Exception as e:
        print(f"❌ 配置创建失败 | Configuration creation failed: {e}")
        sys.exit(1)

    # 创建分析器并运行 | Create analyzer and run
    analyzer = MCycAnalyzer(config)
    success = analyzer.run_full_analysis()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()