"""
亚基因组归属主程序模块|Subgenome Assignment Main Module

编排 combine → align → assign → split 完整流程
|Orchestrates the full combine → align → assign → split pipeline
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from .config import SubgenomeAssignConfig
from .utils import SubgenomeLogger, CommandRunner, check_dependencies
from .aligner import ParentAligner
from .assigner import SubgenomeAssigner
from .splitter import FastaSplitter


def parse_parents_arg(parent_list: List[str]) -> Dict[str, List[str]]:
    """
    解析 --parent NAME:hap1,hap2 形式的参数列表
    |Parse --parent NAME:hap1,hap2 style arguments

    Args:
        parent_list: ["Ca:Ca_hap1.fa,Ca_hap2.fa", "Ch:Ch_hap1.fa,Ch_hap2.fa"]

    Returns:
        {"Ca": ["Ca_hap1.fa", "Ca_hap2.fa"], "Ch": ["Ch_hap1.fa", "Ch_hap2.fa"]}
    """
    result = {}
    for spec in parent_list:
        if ':' not in spec:
            raise ValueError(
                f"--parent 参数格式错误|invalid --parent format: {spec}. "
                f"应为|expected NAME:hap1.fa,hap2.fa"
            )
        name, haps_str = spec.split(':', 1)
        name = name.strip()
        haps = [h.strip() for h in haps_str.split(',') if h.strip()]
        if not name:
            raise ValueError(f"--parent 亲本名为空|empty parent name: {spec}")
        if not haps:
            raise ValueError(f"--parent {name} 没有提供 hap 文件|no hap files")
        result[name] = haps
    return result


class SubgenomeAssignRunner:
    """亚基因组归属运行器|Subgenome Assignment Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = SubgenomeAssignConfig(**kwargs)

        # 校验|Validate
        errors = self.config.validate()
        if errors:
            raise ValueError("\n".join(errors))

        # 初始化日志|Initialize logging
        self.logger_manager = SubgenomeLogger(self.config.log_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各组件|Initialize components
        self.aligner = ParentAligner(self.config, self.logger, self.cmd_runner)
        self.assigner = SubgenomeAssigner(self.config, self.logger)
        self.splitter = FastaSplitter(self.config, self.logger, self.cmd_runner)

        # 流程起始时间|Pipeline start time
        self.start_time = datetime.now()

    def run(self) -> bool:
        """运行完整归属流程|Run complete assignment pipeline"""
        self.logger.info("=" * 60)
        self.logger.info(
            f"开始亚基因组归属流程|Starting subgenome assignment pipeline "
            f"({len(self.config.parents)} parents)"
        )
        self.logger.info("=" * 60)

        # 列出亲本配置|List parents
        for name, haps in self.config.parents.items():
            self.logger.info(
                f"  亲本|Parent {name}: {len(haps)} 个 hap|haps"
            )

        try:
            # Step 1: 检查依赖|Check dependencies
            self.logger.info("步骤1/5: 检查依赖|Step 1/5: Checking dependencies")
            if not check_dependencies(self.config, self.logger):
                return False

            # Step 2: 合并亲本 hap FASTA|Combine parental hap FASTAs
            self.logger.info("步骤2/5: 合并亲本 FASTA|Step 2/5: Combining parental FASTAs")
            if not self.aligner.combine_parent_haps():
                return False

            # Step 3: minimap2 比对|minimap2 alignment
            self.logger.info("步骤3/5: 目标基因组 vs 各亲本比对|Step 3/5: Aligning target vs parents")
            paf_paths = self.aligner.align_target_to_all_parents()
            if not paf_paths:
                self.logger.error("比对失败|Alignment failed")
                return False

            # Step 4: 解析 PAF + 归属判定|Parse PAFs + assignment
            self.logger.info("步骤4/5: 归属判定|Step 4/5: Assigning subgenomes")

            # 先确保 .fai 存在|Ensure target .fai exists
            target_fai = self._ensure_fai()
            if not target_fai:
                return False

            rows = self.assigner.assign(target_fai, paf_paths)
            tsv_file = self.config.assignment_dir / "subgenome_assignment.tsv"
            self.assigner.write_tsv(rows, tsv_file)
            self.assigner.summarize(rows)

            # Step 5: 拆分 FASTA|Split FASTA
            self.logger.info("步骤5/5: 拆分 FASTA|Step 5/5: Splitting FASTA")
            split_results = self.splitter.split(rows, self.config.target, self.config.split_dir)

            # 元数据|Metadata
            self._write_pipeline_info(rows, split_results)

            self.logger.info("=" * 60)
            self.logger.info("亚基因组归属流程完成|Subgenome assignment completed")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info("=" * 60)
            return True

        except Exception as e:
            self.logger.error(f"流程异常|Pipeline error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _ensure_fai(self) -> Path:
        """确保目标基因组 .fai 存在|Ensure target genome .fai exists"""
        target = self.config.target
        # .fai 命名规则：.fa → .fai; .fa.gz → .fai（samtools 行为）
        # |.fai naming: .fa -> .fai; .fa.gz -> .fai (samtools behavior)
        if target.endswith('.gz'):
            fai = Path(target[:-3] + '.fai') if not Path(target[:-3] + '.fai').exists() else Path(target[:-3] + '.fai')
            # samtools faidx 对 .gz 文件生成的 .fai 路径规则较复杂
            # 简化：直接在 alignment_dir 里生成 faidx 后查找
        # 一般情况 .fai 是 target + ".fai"
        fai = Path(str(target) + '.fai')
        if not fai.exists():
            self.logger.info(f"构建 .fai 索引|Building .fai: {fai.name}")
            cmd = [self.config.samtools_path, 'faidx', target]
            if not self.cmd_runner.run(cmd, "samtools faidx|target .fai"):
                return None

        # 重新检查|Re-check
        if not fai.exists():
            # 兜底：在工作目录下手动建一个|Fallback: build in alignment_dir
            fai = self.config.alignment_dir / (Path(target).name + '.fai')
            if not fai.exists():
                self.logger.error(f"无法找到/生成 .fai|Cannot find/create .fai for: {target}")
                return None
        return fai

    def _write_pipeline_info(self, rows: list, split_results: dict):
        """生成 software_versions.yml 和 pipeline_params.yaml|Generate metadata"""
        self.logger.info("写入流程元数据|Writing pipeline metadata")
        self._write_software_versions()
        self._write_pipeline_params(rows, split_results)

    def _write_software_versions(self):
        """生成 software_versions.yml|Generate software_versions.yml"""
        end_time = datetime.now()
        runtime = int((end_time - self.start_time).total_seconds())

        versions = {}
        for tool_name, tool_path, version_args in [
            ('minimap2', self.config.minimap2_path, ['--version']),
            ('samtools', self.config.samtools_path, ['--version']),
        ]:
            try:
                result = subprocess.run(
                    [tool_path] + version_args,
                    capture_output=True, text=True, timeout=30,
                )
                first_line = (result.stdout or result.stderr).strip().split('\n')[0]
                versions[tool_name] = {
                    'version': first_line,
                    'path': tool_path,
                }
            except Exception as e:
                versions[tool_name] = {
                    'version': 'unknown', 'path': tool_path, 'error': str(e),
                }

        info = {
            'pipeline': {
                'name': 'biopytools.subgenome_assign',
                'version': '1.0.0',
            },
            'tools': versions,
            'parameters': self.config.to_param_dict(),
            'execution': {
                'start_time': self.start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime,
            },
        }

        out_file = self.config.pipeline_info_dir / "software_versions.yml"
        with open(out_file, 'w') as f:
            self._dump_yaml(info, f, 0)
        self.logger.info(f"软件版本|Software versions: {out_file}")

    def _write_pipeline_params(self, rows: list, split_results: dict):
        """生成 pipeline_params.yaml|Generate pipeline_params.yaml"""
        from collections import Counter
        out_file = self.config.pipeline_info_dir / "pipeline_params.yaml"
        params = self.config.to_param_dict()
        params['result_summary'] = {
            'total_chromosomes': len(rows),
            'per_parent_count': dict(Counter(r['assigned'] for r in rows)),
            'split_fastas': split_results,
        }
        with open(out_file, 'w') as f:
            self._dump_yaml(params, f, 0)
        self.logger.info(f"流程参数|Pipeline params: {out_file}")

    def _dump_yaml(self, data, f, indent: int):
        """简易YAML序列化|Minimal YAML serializer"""
        pad = '  ' * indent
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, (dict, list)):
                    f.write(f"{pad}{k}:\n")
                    self._dump_yaml(v, f, indent + 1)
                else:
                    f.write(f"{pad}{k}: {self._yaml_scalar(v)}\n")
        elif isinstance(data, list):
            for item in data:
                if isinstance(item, (dict, list)):
                    f.write(f"{pad}-\n")
                    self._dump_yaml(item, f, indent + 1)
                else:
                    f.write(f"{pad}- {self._yaml_scalar(item)}\n")

    @staticmethod
    def _yaml_scalar(v) -> str:
        """简单YAML标量格式化|Simple YAML scalar formatter"""
        if v is None:
            return 'null'
        if isinstance(v, bool):
            return 'true' if v else 'false'
        if isinstance(v, (int, float)):
            return str(v)
        s = str(v)
        if any(c in s for c in ':#{}[]&*!|>"\'@\n'):
            return f'"{s}"'
        return s


def main():
    """命令行入口|Command-line entry"""
    parser = argparse.ArgumentParser(
        description='亚基因组归属工具|Subgenome Assignment Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 必需参数|Required
    req = parser.add_argument_group('必需|Required')
    req.add_argument('-i', '--target', required=True,
                     help='目标多倍体基因组 FASTA|Target polyploid genome FASTA')
    req.add_argument('--parent', action='append', required=True, metavar='NAME:hap1.fa,hap2.fa',
                     help='亲本名:hap1,hap2,...（可重复指定多个亲本）'
                          '|Parent name:hap1,hap2,... (can be repeated for multiple parents)')

    # 输出|Output
    out_g = parser.add_argument_group('输出|Output')
    out_g.add_argument('-o', '--output-dir', default='./subgenome_assign_output',
                       help='输出目录|Output directory')

    # 比对参数|Alignment
    aln_g = parser.add_argument_group('比对|Alignment')
    aln_g.add_argument('--preset', default='asm10',
                       choices=['asm5', 'asm10', 'asm20', 'asm25'],
                       help='minimap2 -x 预设（asm5=<1%%差异, asm10=1-5%%, asm20=5-15%%）'
                          '|minimap2 preset')
    aln_g.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    aln_g.add_argument('--minimap2-secondary', action='store_true',
                       help='保留次要比对（默认 --secondary=no）|Keep secondary alignments')

    # 归属判定|Assignment
    asg_g = parser.add_argument_group('归属判定|Assignment')
    asg_g.add_argument('--min-conf', type=float, default=0.65,
                       help='置信度阈值，低于此值标记 LOW_CONFIDENCE'
                          '|Confidence threshold for LOW_CONFIDENCE flag')

    # FASTA 拆分|FASTA splitting
    spl_g = parser.add_argument_group('FASTA 拆分|FASTA splitting')
    spl_g.add_argument('--no-split', action='store_true',
                       help='不输出拆分的 FASTA|Do not output split FASTAs')
    spl_g.add_argument('--no-keep-unassigned', action='store_true',
                       help='不输出未归属染色体的 FASTA|Do not output unassigned FASTA')

    # 工具路径|Tool paths
    tool_g = parser.add_argument_group('工具路径|Tool paths')
    tool_g.add_argument('--minimap2-path',
                        default='~/miniforge3/envs/cphasing/bin/minimap2',
                        help='minimap2 二进制路径|minimap2 binary path')
    tool_g.add_argument('--samtools-path',
                        default='~/.local/bin/samtools',
                        help='samtools 二进制路径|samtools binary path')

    args = parser.parse_args()

    try:
        # 解析 --parent 列表|Parse parent list
        parents = parse_parents_arg(args.parent)

        runner = SubgenomeAssignRunner(
            target=args.target,
            parents=parents,
            output_dir=args.output_dir,
            preset=args.preset,
            threads=args.threads,
            minimap2_secondary=args.minimap2_secondary,
            min_conf=args.min_conf,
            split_fasta=not args.no_split,
            keep_unassigned=not args.no_keep_unassigned,
            minimap2_path=args.minimap2_path,
            samtools_path=args.samtools_path,
        )
        success = runner.run()
        sys.exit(0 if success else 1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
