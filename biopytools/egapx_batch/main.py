"""
EGAPx批量生成器主模块 | EGAPx Batch Generator Main Module
"""

import os
import sys
import time
import shutil
from .config import EGAPxBatchConfig
from .utils import EGAPxBatchLogger, GenomeSplitter, TemplateProcessor, JobGenerator


class EGAPxBatchGenerator:
    """EGAPx批量运行生成器 | EGAPx Batch Processing Generator"""

    def __init__(self, **kwargs):
        """
        初始化生成器 | Initialize generator

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        # 初始化配置 | Initialize configuration
        self.config = EGAPxBatchConfig(**kwargs)
        self.config.validate()

        # 创建输出目录 | Create output directory
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 初始化日志 | Initialize logging
        self.logger_manager = EGAPxBatchLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类 | Initialize utility classes
        self.splitter = GenomeSplitter(self.logger)
        self.template_processor = TemplateProcessor(self.logger)
        self.job_generator = JobGenerator(self.logger)

    def generate(self):
        """
        执行批量生成 | Perform batch generation

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("开始EGAPx批量运行配置生成 | Starting EGAPx Batch Config Generation")
            self.logger.info("=" * 60)
            self.logger.info(f"配置信息 | Configuration:\n{self.config}")

            # 步骤1: 拆分基因组 | Step 1: Split genome
            self.logger.info("=" * 60)
            self.logger.info("步骤 1/3: 拆分基因组文件 | Step 1/3: Split genome file")
            self.logger.info("=" * 60)

            split_dir = os.path.join(self.config.output_dir, "temp_split")

            # 获取输出目录的绝对路径 | Get absolute path of output directory
            output_dir_abs = os.path.abspath(self.config.output_dir)

            # 使用awk快速拆分 | Use awk to split (much faster)
            split_files = self.splitter.split_genome(
                self.config.genome,
                split_dir,
                self.config.chr_prefix
            )

            if not split_files:
                self.logger.error("未找到任何序列 | No sequences found")
                return False

            self.logger.info(f"成功拆分 {len(split_files)} 个序列 | Successfully split {len(split_files)} sequences")

            # 步骤2: 生成配置文件和脚本 | Step 2: Generate configs and scripts
            self.logger.info("=" * 60)
            self.logger.info("步骤 2/3: 生成配置文件和脚本 | Step 2/3: Generate configs and scripts")
            self.logger.info("=" * 60)

            job_list_path = os.path.join(self.config.output_dir, "all_jobs_submit.list.sh")

            with open(job_list_path, 'w') as job_list:
                for chr_file in split_files:
                    chr_name = os.path.basename(chr_file).replace('.fa', '')
                    chr_dir = os.path.join(self.config.output_dir, chr_name)
                    chr_dir_abs = os.path.join(output_dir_abs, chr_name)

                    # 创建目录结构 | Create directory structure
                    os.makedirs(chr_dir, exist_ok=True)
                    os.makedirs(os.path.join(chr_dir, 'work'), exist_ok=True)
                    os.makedirs(os.path.join(chr_dir, 'output'), exist_ok=True)

                    # 移动染色体文件 | Move chromosome file
                    chr_fa = os.path.join(chr_dir_abs, f"{chr_name}.fa")
                    chr_fa_relative = os.path.join(chr_dir, f"{chr_name}.fa")
                    shutil.move(chr_file, chr_fa_relative)

                    # 创建EGAPx软链接 | Create EGAPx symlinks
                    self._create_egapx_symlinks(chr_dir)

                    # 生成YAML配置 | Generate YAML config
                    new_yaml = os.path.join(chr_dir_abs, f"{chr_name}.yaml")
                    # 构建locus标签前缀 | Build locus tag prefix
                    if self.config.locus_tag_prefix:
                        new_prefix = f"{self.config.locus_tag_prefix}_{chr_name}"
                    else:
                        new_prefix = chr_name

                    yaml_content = self.job_generator.generate_yaml(
                        genome=chr_fa,
                        locus_tag_prefix=new_prefix,
                        short_reads=self.config.short_reads,
                        long_reads=self.config.long_reads
                    )

                    with open(new_yaml, 'w') as f:
                        f.write(yaml_content)

                    # 生成运行脚本 | Generate run script
                    new_script = os.path.join(chr_dir_abs, f"egapx_{chr_name}.sh")
                    new_report = f"{self.config.report_name}_{chr_name}"

                    script_content = self.job_generator.generate_script(
                        yaml_path=new_yaml,
                        work_dir=os.path.join(chr_dir_abs, 'work'),
                        output_dir=os.path.join(chr_dir_abs, 'output'),
                        report_name=new_report,
                        chr_dir=chr_dir_abs
                    )

                    with open(new_script, 'w') as f:
                        f.write(script_content)

                    # 设置可执行权限 | Set executable permission
                    os.chmod(new_script, 0o755)

                    # 添加到任务列表 | Add to job list
                    job_list.write(f"bash \"{new_script}\"\n")

                    self.logger.info(f"  {chr_name} 配置完成 | {chr_name} config completed")

            # 步骤3: 生成并行执行脚本 | Step 3: Generate parallel execution script
            self.logger.info("=" * 60)
            self.logger.info("步骤 3/3: 生成并行执行脚本 | Step 3/3: Generate parallel execution script")
            self.logger.info("=" * 60)

            parallel_script = os.path.join(self.config.output_dir, "run_all_parallel.sh")
            self._generate_parallel_script(parallel_script, job_list_path)

            # 清理临时目录 | Cleanup temporary directory
            shutil.rmtree(split_dir)

            # 完成统计 | Completion statistics
            elapsed_time = time.time() - start_time
            self.config.job_count = len(split_files)

            self.logger.info("=" * 60)
            self.logger.info("配置生成完成！| Config generation completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"统计信息 | Statistics:")
            self.logger.info(f"  - 任务数量 | Job count: {self.config.job_count}")
            self.logger.info(f"  - 任务列表 | Job list: {job_list_path}")
            self.logger.info(f"  - 并行脚本 | Parallel script: {parallel_script}")
            self.logger.info(f"  - 总耗时 | Total time: {elapsed_time:.2f} seconds")

            self._show_usage_guide(job_list_path, parallel_script)

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错 | Program execution error: {str(e)}", exc_info=True)
            return False

    def _create_egapx_symlinks(self, chr_dir):
        """创建EGAPx软链接 | Create EGAPx symlinks"""
        symlink_record = os.path.join(chr_dir, '.egapx_symlinks')

        try:
            # 创建软链接 | Create symlinks
            for item in os.listdir(self.config.egapx_path):
                src = os.path.join(self.config.egapx_path, item)
                dst = os.path.join(chr_dir, item)
                if os.path.exists(src):
                    os.symlink(src, dst)

            # 记录软链接 | Record symlinks
            symlinks = [f for f in os.listdir(chr_dir)
                       if os.path.islink(os.path.join(chr_dir, f))]
            with open(symlink_record, 'w') as f:
                for link in symlinks:
                    f.write(f"{link}\n")

        except Exception as e:
            self.logger.warning(f"创建软链接失败 | Failed to create symlinks: {e}")

    def _generate_parallel_script(self, script_path, job_list_path):
        """生成并行执行脚本 | Generate parallel execution script"""
        content = f"""#!/bin/bash
# EGAPx并行执行脚本 | EGAPx parallel execution script
# 自动生成时间 | Auto-generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
JOB_LIST="$SCRIPT_DIR/all_jobs_submit.list.sh"

if [[ ! -f "$JOB_LIST" ]]; then
    echo "[ERROR] 任务列表不存在 | Job list not found: $JOB_LIST"
    exit 1
fi

JOBS=${{1:-4}}  # 默认并行4个任务 | Default parallel 4 jobs

echo "[INFO] 开始并行执行 | Starting parallel execution, 并行度 | parallelism: $JOBS"
echo "[INFO] 任务列表 | Job list: $JOB_LIST"
echo ""

if command -v parallel &> /dev/null; then
    echo "[INFO] 使用 GNU parallel | Using GNU parallel"
    cat "$JOB_LIST" | parallel -j "$JOBS"
else
    echo "[WARN] GNU parallel 未安装，使用 xargs | GNU parallel not found, using xargs"
    cat "$JOB_LIST" | xargs -P "$JOBS" -I {{}} bash -c {{}}
fi

echo ""
echo "[INFO] 所有任务执行完成 | All jobs completed"
"""

        with open(script_path, 'w') as f:
            f.write(content)

        os.chmod(script_path, 0o755)
        self.logger.info(f"并行脚本已生成 | Parallel script generated: {script_path}")

    def _show_usage_guide(self, job_list_path, parallel_script):
        """显示使用指南 | Show usage guide"""
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("执行方式 | Execution methods:")
        self.logger.info("=" * 60)
        self.logger.info("[1] 顺序执行 | Serial execution:")
        self.logger.info(f"   bash {job_list_path}")
        self.logger.info("")
        self.logger.info("[2] 并行执行（推荐）| Parallel execution (recommended):")
        self.logger.info(f"   bash {parallel_script} 4")
        self.logger.info("")
        self.logger.info("[3] 使用 GNU parallel:")
        self.logger.info(f"   cat {job_list_path} | parallel -j 4")
        self.logger.info("")
        self.logger.info("[4] 使用 xargs:")
        self.logger.info(f"   cat {job_list_path} | xargs -P 4 -I {{}} bash -c {{}}")
        self.logger.info("=" * 60)


def main():
    """主函数 | Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='🧬 EGAPx批量运行配置生成工具 | EGAPx Batch Config Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 🚀 基本用法 (使用内置模板) | Basic usage (with built-in templates)
  %(prog)s \\
      -g genome.fa \\
      -o output_dir

  # 🏷️ 自定义locus标签前缀和报告名 | Custom locus tag prefix and report name
  %(prog)s \\
      -g genome.fa \\
      -o output_dir \\
      --locus-prefix Gene \\
      --report-name MyEGAPx

  # 📋 使用自定义模板 | Use custom templates
  %(prog)s \\
      -g genome.fa \\
      -o output_dir \\
      -y template.yaml \\
      -s template.sh

  # 🔍 只处理特定前缀的染色体 | Only process chromosomes with specific prefix
  %(prog)s \\
      -g genome.fa \\
      -o output_dir \\
      --chr-prefix Chr
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments | 必需参数')
    required.add_argument('-g', '--genome', required=True,
                         help='[FILE] 基因组FASTA文件路径 | Genome FASTA file path')
    required.add_argument('-o', '--output', required=True,
                         help='[DIR] 输出目录路径 | Output directory path')

    # 可选参数 | Optional parameters
    optional = parser.add_argument_group('optional arguments | 可选参数')
    optional.add_argument('-y', '--yaml',
                         default=None,
                         help='[FILE] YAML模板文件路径 (可选，内置默认模板) | YAML template file path (optional, built-in default template)')
    optional.add_argument('-s', '--script',
                         default=None,
                         help='[FILE] Shell脚本模板路径 (可选，内置默认模板) | Shell script template path (optional, built-in default template)')
    optional.add_argument('-e', '--egapx',
                         default='/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx',
                         help='[PATH] EGAPx安装路径 (默认: /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx) | EGAPx installation path')
    optional.add_argument('-p', '--chr-prefix',
                         help='[STR] 染色体前缀过滤 | Chromosome prefix filter')
    optional.add_argument('--locus-prefix', default='',
                         help='[STR] locus标签前缀 (默认: 空) | Locus tag prefix (default: empty)')
    optional.add_argument('--report-name', default='EGAPx',
                         help='[STR] 报告名称 (默认: EGAPx) | Report name (default: EGAPx)')
    optional.add_argument('--short-reads', default='',
                         help='[FILE] 短读长测序数据文件路径 | Short reads file path')
    optional.add_argument('--long-reads', default='',
                         help='[FILE] 长读长测序数据文件路径 | Long reads file path')

    args = parser.parse_args()

    # 创建生成器并运行 | Create generator and run
    generator = EGAPxBatchGenerator(
        genome=args.genome,
        yaml_template=args.yaml,
        script_template=args.script,
        output_dir=args.output,
        egapx_path=args.egapx,
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
