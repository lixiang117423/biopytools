"""
PGGB泛基因组图构建核心模块|PGGB Pangenome Graph Builder Core Module
"""

import os
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime
from typing import Optional, List

from .config import PGGBConfig
from .utils import (
    PGGBLogger,
    build_conda_command,
    check_pggb_dependencies,
    ensure_fasta_index,
    format_number,
    format_size
)


class PGGBRunner:
    """PGGB泛基因组图构建运行器|PGGB Pangenome Graph Builder Runner"""

    def __init__(self, config: PGGBConfig, logger: Optional[PGGBLogger] = None):
        """
        初始化运行器|Initialize runner

        Args:
            config: PGGB配置对象|PGGB configuration object
            logger: 日志对象|Logger object (可选|optional)
        """
        self.config = config
        self.logger_obj = logger
        self.logger = logger.get_logger() if logger else None
        self.start_time = None

    def run(self) -> bool:
        """
        运行PGGB泛基因组图构建|Run PGGB pangenome graph building

        Returns:
            是否成功|Whether successful
        """
        try:
            self.start_time = time.time()

            # 1. 设置日志|Setup logging
            if self.logger_obj is None:
                log_dir = Path(self.config.output_dir) / "99_logs"
                log_dir.mkdir(parents=True, exist_ok=True)
                self.logger_obj = PGGBLogger(log_dir / "pggb.log")
                self.logger = self.logger_obj.get_logger()

            # 2. 打印开始信息|Print start information
            self._print_header()

            # 3. 验证配置|Validate configuration
            self.config.validate()

            # 4. 检查依赖|Check dependencies
            if not check_pggb_dependencies(self.config, self.logger):
                return False

            # 5. 确保FASTA索引存在|Ensure FASTA index exists
            if not ensure_fasta_index(self.config.input_fa, self.logger):
                return False

            # 6. 构建并执行命令|Build and execute command
            success = self._run_pggb()

            # 7. 打印结束信息|Print end information
            self._print_footer(success)

            return success

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Program execution error: {e}", exc_info=True)
            return False

    def _print_header(self):
        """打印头部信息|Print header information"""
        if not self.logger:
            return

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("PGGB 泛基因组图构建|PGGB Pangenome Graph Builder")
        self.logger.info("=" * 60)
        self.logger.info(f"版本|Version: 1.0.0")
        self.logger.info(f"时间|Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  输入文件|Input FASTA: {self.config.input_fa}")

        input_size = Path(self.config.input_fa).stat().st_size
        self.logger.info(f"  输入大小|Input size: {format_size(input_size)}")

        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  conda环境|Conda env: {self.config.conda_env}")
        self.logger.info(f"  段长度|Segment length: {self.config.segment_length}")
        self.logger.info(f"  一致度|Map pct id: {self.config.map_pct_id}")

        if self.config.n_haplotypes > 0:
            self.logger.info(f"  单倍型数|N haplotypes: {self.config.n_haplotypes}")
        else:
            self.logger.info(f"  单倍型数|N haplotypes: auto")

        if self.config.vcf_spec:
            self.logger.info(f"  VCF参考|VCF spec: {self.config.vcf_spec}")

        if self.config.resume:
            self.logger.info(f"  断点续传|Resume: enabled")

        if self.config.input_paf:
            self.logger.info(f"  外部PAF|Input PAF: {self.config.input_paf}")

        self.logger.info("")

    def _build_pggb_command(self) -> List[str]:
        """构建pggb命令|Build pggb command"""
        args = [
            self.config.input_fa,
            '-o', self.config.output_dir,
            '-t', str(self.config.threads),
        ]

        # wfmash参数|wfmash parameters
        args.extend(['-s', str(self.config.segment_length)])
        args.extend(['-p', str(self.config.map_pct_id)])
        args.extend(['-c', str(self.config.n_mappings)])
        args.extend(['-K', str(self.config.mash_kmer)])

        # mash_kmer_thres需要特殊处理(避免科学计数法问题)
        thres = self.config.mash_kmer_thres
        if thres == 0.001:
            args.extend(['-F', '0.001'])
        else:
            args.extend(['-F', str(thres)])

        if self.config.block_length > 0:
            args.extend(['-l', str(self.config.block_length)])

        if self.config.no_splits:
            args.append('-N')

        if self.config.sparse_map:
            args.extend(['-x', self.config.sparse_map])

        if self.config.exclude_delim:
            args.extend(['-Y', self.config.exclude_delim])

        # seqwish参数|seqwish parameters
        args.extend(['-k', str(self.config.min_match_length)])

        if self.config.sparse_factor > 0:
            args.extend(['-f', str(self.config.sparse_factor)])

        args.extend(['-B', self.config.transclose_batch])

        # smoothxg参数|smoothxg parameters
        if self.config.skip_normalization:
            args.append('-X')

        if self.config.n_haplotypes > 0:
            args.extend(['-n', str(self.config.n_haplotypes)])

        if self.config.max_path_jump > 0:
            args.extend(['-j', str(self.config.max_path_jump)])

        if self.config.max_edge_jump > 0:
            args.extend(['-e', str(self.config.max_edge_jump)])

        args.extend(['-G', self.config.target_poa_length])

        if self.config.poa_params:
            args.extend(['-P', self.config.poa_params])

        args.extend(['-O', str(self.config.poa_padding)])
        args.extend(['-d', str(self.config.pad_max_depth)])

        # odgi/vg参数|odgi/vg parameters
        if self.config.skip_viz:
            args.append('-v')

        if self.config.do_stats:
            args.append('-S')

        if self.config.vcf_spec:
            args.extend(['-V', self.config.vcf_spec])

        # 通用参数|General parameters
        if self.config.poa_threads > 0:
            args.extend(['-T', str(self.config.poa_threads)])

        if self.config.resume:
            args.append('-r')

        if self.config.compress:
            args.append('-Z')

        if self.config.keep_temp_files:
            args.append('-A')

        if self.config.input_paf:
            args.extend(['-a', self.config.input_paf])

        if self.config.temp_dir:
            args.extend(['-D', self.config.temp_dir])

        return args

    def _run_pggb(self) -> bool:
        """执行pggb命令|Execute pggb command"""
        pggb_args = self._build_pggb_command()

        # 使用conda run包装|Wrap with conda run
        cmd = build_conda_command('pggb', pggb_args)

        self.logger.info("=" * 60)
        self.logger.info("执行PGGB|Executing PGGB")
        self.logger.info("=" * 60)
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info("")

        try:
            result = subprocess.run(
                cmd,
                cwd=None,
                timeout=None  # PGGB运行时间可能很长，不设超时
            )

            if result.returncode == 0:
                self.logger.info("")
                self.logger.info("PGGB 执行完成|PGGB completed successfully")
                self._summarize_outputs()
                return True
            else:
                self.logger.error(f"PGGB 执行失败|PGGB failed with return code: {result.returncode}")
                return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"PGGB 执行出错|PGGB execution error: {e}")
            return False
        except Exception as e:
            self.logger.error(f"未知错误|Unknown error: {e}")
            return False

    def _summarize_outputs(self):
        """总结输出文件|Summarize output files"""
        output_dir = Path(self.config.output_dir)

        self.logger.info("")
        self.logger.info("输出文件|Output files:")
        self.logger.info("-" * 40)

        # 按类型分类输出文件|Categorize output files by type
        output_patterns = {
            'GFA': ['*.gfa', '*.gfa.gz'],
            'OG': ['*.og', '*.og.gz'],
            'GBZ': ['*.gbz'],
            'PAF': ['*.paf', '*.paf.gz'],
            'VCF': ['*.vcf.gz', '*.vcf'],
            '其他|Other': ['*.log', '*.yml', '*.maf', '*.maf.gz'],
        }

        total_size = 0
        file_count = 0

        for category, patterns in output_patterns.items():
            files = []
            for pattern in patterns:
                files.extend(output_dir.glob(pattern))

            if files:
                for f in sorted(files, key=lambda x: x.stat().st_size, reverse=True):
                    size = f.stat().st_size
                    total_size += size
                    file_count += 1
                    self.logger.info(f"  [{category}] {f.name} ({format_size(size)})")

        self.logger.info("-" * 40)
        self.logger.info(f"  共 {file_count} 个文件|Total {file_count} files, {format_size(total_size)}")

    def _print_footer(self, success: bool):
        """打印结束信息|Print footer information"""
        if not self.logger:
            return

        elapsed = time.time() - self.start_time
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)
        seconds = int(elapsed % 60)

        self.logger.info("")
        self.logger.info("=" * 60)
        if success:
            self.logger.info("PGGB 泛基因组图构建完成|PGGB Pangenome Graph Builder Completed")
        else:
            self.logger.info("PGGB 泛基因组图构建失败|PGGB Pangenome Graph Builder Failed")
        self.logger.info(f"总耗时|Total time: {hours}h {minutes}m {seconds}s")
        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='PGGB泛基因组图构建工具|PGGB Pangenome Graph Builder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('required arguments|必需参数')
    required.add_argument('-i', '--input', required=True,
                         help='[FILE] 输入FASTA文件|Input FASTA file')
    required.add_argument('-o', '--output', required=True,
                         help='[DIR] 输出目录|Output directory')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('optional arguments|可选参数')
    optional.add_argument('-t', '--threads', type=int, default=24,
                         help='[INT] 线程数 (default: 24)')
    optional.add_argument('--conda-env', default='pggb_v.0.7.4',
                         help='[STR] conda环境名 (default: pggb_v.0.7.4)')

    # wfmash参数|wfmash parameters
    wfmash = parser.add_argument_group('wfmash parameters|wfmash参数')
    wfmash.add_argument('-s', '--segment-length', type=int, default=5000,
                       help='[INT] 比对分段长度 (default: 5000)')
    wfmash.add_argument('-l', '--block-length', type=int, default=0,
                       help='[INT] 最小block长度 (default: 0, auto=5*segment-length)')
    wfmash.add_argument('-p', '--map-pct-id', type=int, default=90,
                       help='[INT] 比对一致度 (default: 90)')
    wfmash.add_argument('-c', '--n-mappings', type=int, default=1,
                       help='[INT] 每segment的mapping数 (default: 1)')
    wfmash.add_argument('-K', '--mash-kmer', type=int, default=19,
                       help='[INT] mash kmer大小 (default: 19)')
    wfmash.add_argument('--no-splits', action='store_true',
                       help='禁用序列拆分|Disable sequence splitting')
    wfmash.add_argument('--sparse-map', default='',
                       help='[STR] 稀疏映射比例|Sparse mapping fraction')
    wfmash.add_argument('--input-paf', default='',
                       help='[FILE] 外部PAF文件(跳过wfmash)|External PAF file')

    # smoothxg参数|smoothxg parameters
    smooth = parser.add_argument_group('smoothxg parameters|smoothxg参数')
    smooth.add_argument('-n', '--n-haplotypes', type=int, default=0,
                       help='[INT] 单倍型数 (default: 0, auto-detect)')
    smooth.add_argument('--skip-normalization', action='store_true',
                       help='跳过图归一化|Skip graph normalization')

    # 输出选项|Output options
    output = parser.add_argument_group('output options|输出选项')
    output.add_argument('--vcf-spec', default='',
                       help='[STR] VCF输出参考规范|VCF output reference spec')
    output.add_argument('--stats', action='store_true',
                       help='生成统计信息|Generate statistics')
    output.add_argument('--resume', action='store_true',
                       help='断点续传|Resume from existing outputs')
    output.add_argument('--compress', action='store_true',
                       help='压缩输出|Compress output files')
    output.add_argument('--keep-temp', action='store_true',
                       help='保留临时文件|Keep intermediate files')

    args = parser.parse_args()

    # 创建配置|Create configuration
    config = PGGBConfig(
        input_fa=args.input,
        output_dir=args.output,
        conda_env=args.conda_env,
        threads=args.threads,
        segment_length=args.segment_length,
        block_length=args.block_length,
        map_pct_id=args.map_pct_id,
        n_mappings=args.n_mappings,
        mash_kmer=args.mash_kmer,
        no_splits=args.no_splits,
        sparse_map=args.sparse_map,
        n_haplotypes=args.n_haplotypes,
        skip_normalization=args.skip_normalization,
        vcf_spec=args.vcf_spec,
        do_stats=args.stats,
        resume=args.resume,
        compress=args.compress,
        keep_temp_files=args.keep_temp,
        input_paf=args.input_paf,
    )

    # 创建运行器并执行|Create runner and execute
    runner = PGGBRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
