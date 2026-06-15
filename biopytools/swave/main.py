"""
Swave主程序模块|Swave Main Module
"""

import argparse
import sys
import os
from .config import SwaveConfig
from .utils import SwaveLogger, CommandRunner, check_dependencies
from .sv_caller import SwaveSVCaller
from .pav_extractor import PAVExtractor
import subprocess


class SwaveRunner:
    """Swave运行器|Swave Runner"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # 对于非call命令，可能不需要完整配置
        self.config = None
        self.logger = None
        self.cmd_runner = None

        if kwargs.get('command') == 'call':
            # 移除command参数，避免传入SwaveConfig
            # Remove command kwarg to prevent passing to SwaveConfig
            config_kwargs = {k: v for k, v in kwargs.items() if k != 'command'}
            self.config = SwaveConfig(**config_kwargs)
            self.config.validate()

            # 初始化日志|Initialize logging
            self.logger_manager = SwaveLogger(
                output_dir=self.config.output_dir,
                log_file='swave_pipeline.log',
                log_level="INFO"
            )
            self.logger = self.logger_manager.get_logger()

            # 初始化命令执行器|Initialize command runner
            self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

            # 检查依赖|Check dependencies
            if not check_dependencies(self.config, self.logger):
                self.logger.error("依赖检查失败|Dependency check failed")
                sys.exit(1)

    def run_call(self):
        """运行SV检测|Run SV detection"""
        self.logger.info("开始Swave SV检测流程|Starting Swave SV detection pipeline")

        # 复用SwaveSVCaller的参数构建，避免重复代码
        # Reuse SwaveSVCaller.build_call_args to avoid code duplication
        swave_script = os.path.join(self.config.swave_path, 'Swave.py')
        swave_args = SwaveSVCaller.build_call_args(self.config)

        python_path, env = SwaveSVCaller._get_swave_env()
        args = [python_path, swave_script] + swave_args

        # 执行命令|Execute command
        self.logger.info(f"执行命令|Executing: {' '.join(args)}")

        try:
            result = subprocess.run(
                args,
                shell=False,
                check=False,
                cwd=self.config.swave_path,
                env=env
            )

            if result.returncode == 0:
                self.logger.info("Swave执行成功|Swave execution completed successfully")

                # 自动将图路径编号转换为实际碱基序列
                # Automatically convert graph path IDs to actual sequences
                self.logger.info("开始自动序列转换|Starting automatic sequence conversion")
                sv_caller = SwaveSVCaller(self.config, self.logger, self.cmd_runner)
                if sv_caller.run_convert_seq():
                    self.logger.info("序列转换成功|Sequence conversion completed")
                else:
                    self.logger.warning("序列转换失败，请手动运行 biopytools swave convert_seq|"
                                       "Sequence conversion failed, please run convert_seq manually")

                return True
            else:
                self.logger.error(f"Swave执行失败|Swave execution failed: return code {result.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {str(e)}")
            return False


def main():
    """主函数|Main function"""
    # 创建顶层解析器|Create top-level parser
    parser = argparse.ArgumentParser(
        description='Swave结构变异检测工具|Swave Structural Variant Detection Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='子命令|Sub-commands')

    # ========== call 子命令 ==========
    parser_call = subparsers.add_parser('call', help='检测结构变异|Call structural variants')

    # 必需参数|Required parameters
    parser_call.add_argument('-i', '--assemblies-tsv', required=True,
                            help='样本组装TSV文件|Assemblies TSV file')
    parser_call.add_argument('-r', '--ref-fasta', required=True,
                            help='参考基因组FASTA文件|Reference genome FASTA file')
    parser_call.add_argument('-g', '--gfa-file', required=True,
                            help='泛基因组图GFA文件|Pangenome graph GFA file')
    parser_call.add_argument('-s', '--gfa-source', required=True,
                            choices=['minigraph', 'cactus', 'pggb'],
                            help='GFA文件来源|GFA file source (minigraph/cactus/pggb)')

    # 路径配置|Path configuration
    parser_call.add_argument('--swave-path',
                            default='~/software/swave/Swave-main',
                            help='Swave软件路径|Swave software path')
    parser_call.add_argument('-o', '--output-dir',
                            default='./swave_output',
                            help='输出目录|Output directory')

    # 可选参数|Optional parameters
    parser_call.add_argument('--decomposed-vcf',
                            help='Decomposed VCF文件(cactus/pggb必需)|Decomposed VCF file (required for cactus/pggb)')
    parser_call.add_argument('--output-mode',
                            default='auto',
                            choices=['auto', 'population', 'single'],
                            help='输出模式|Output mode')
    parser_call.add_argument('--spec-samples',
                            nargs='+',
                            help='指定样本|Specify samples')

    # SV检测参数|SV detection parameters
    parser_call.add_argument('--min-sv-size',
                            type=int, default=50,
                            help='最小SV大小|Minimum SV size')
    parser_call.add_argument('--max-sv-size',
                            type=int, default=1000000,
                            help='最大SV大小|Maximum SV size')
    parser_call.add_argument('--max-sv-comps',
                            type=int, default=5,
                            help='最大SV组件数|Maximum number of SV components')

    # 处理选项|Processing options
    parser_call.add_argument('--dup-to-ins',
                            action='store_true',
                            help='将duplication报告为insertion|Report duplications as insertions')
    parser_call.add_argument('--remove-small',
                            action='store_true',
                            help='移除小于min_sv_size的节点|Remove nodes smaller than min_sv_size')
    parser_call.add_argument('--force-reverse',
                            action='store_true',
                            help='强制调用反向映射snarls|Force call reversed mapping snarls')

    # 性能参数|Performance parameters
    parser_call.add_argument('-t', '--threads',
                            type=int, default=12,
                            help='线程数|Number of threads')

    # 外部工具路径|External tool paths
    parser_call.add_argument('--minigraph-path',
                            default='minigraph',
                            help='minigraph工具路径|minigraph tool path')
    parser_call.add_argument('--gfatools-path',
                            default='gfatools',
                            help='gfatools工具路径|gfatools tool path')

    # 高级选项|Advanced options
    parser_call.add_argument('--spec-snarl',
                            help='只调用特定snarl|Only call specific snarl')
    parser_call.add_argument('--spec-path',
                            help='只调用特定path|Only call specific path')

    # ========== convert_seq 子命令 ==========
    parser_convert = subparsers.add_parser('convert_seq',
                                           help='转换图路径为序列|Convert graph paths to sequences')
    parser_convert.add_argument('--vcf-path', required=True,
                               help='VCF文件路径|VCF file path')
    parser_convert.add_argument('--gfa-path', required=True,
                               help='GFA文件路径|GFA file path')
    parser_convert.add_argument('--ref-path', required=True,
                               help='参考基因组FASTA文件|Reference genome FASTA file')
    parser_convert.add_argument('--swave-path',
                               default='~/software/swave/Swave-main',
                               help='Swave软件路径|Swave software path')
    parser_convert.add_argument('--output-path',
                               help='输出路径|Output path')
    parser_convert.add_argument('--force-pangenie',
                               action='store_true',
                               help='强制输出满足pangenie要求的序列|Force output sequences to meet pangenie requirements')

    # ========== convert_Plines 子命令 ==========
    parser_convert_plines = subparsers.add_parser('convert_Plines',
                                                  help='转换VCF为GFA P lines|Convert VCF to GFA P lines')
    parser_convert_plines.add_argument('--gfa-path', required=True,
                                      help='GFA文件路径|GFA file path')
    parser_convert_plines.add_argument('--vcf-path', required=True,
                                      help='VCF文件路径|VCF file path')
    parser_convert_plines.add_argument('--ref-vcf-path',
                                      help='参考VCF文件路径|Reference VCF file path')
    parser_convert_plines.add_argument('--swave-path',
                                      default='~/software/swave/Swave-main',
                                      help='Swave软件路径|Swave software path')
    parser_convert_plines.add_argument('--output-path',
                                      help='输出路径|Output path')
    parser_convert_plines.add_argument('--force-vg',
                                      action='store_true',
                                      help='强制输出满足vg要求的序列|Force output sequences to meet vg requirements')

    # ========== extract_csv 子命令 ==========
    parser_extract_csv = subparsers.add_parser('extract_csv',
                                               help='从VCF提取CSV|Extract CSV from VCF')
    parser_extract_csv.add_argument('--vcf-path', required=True,
                                   help='VCF文件路径|VCF file path')
    parser_extract_csv.add_argument('--swave-path',
                                   default='~/software/swave/Swave-main',
                                   help='Swave软件路径|Swave software path')
    parser_extract_csv.add_argument('--spec-csv',
                                   choices=['INV', 'DUP', 'All'],
                                   help='特定CSV类型|Specific CSV type (INV/DUP/All)')
    parser_extract_csv.add_argument('--output-path',
                                   help='输出路径|Output path')

    # ========== extract_sample 子命令 ==========
    parser_extract_sample = subparsers.add_parser('extract_sample',
                                                 help='提取特定样本的SV|Extract SVs for specific samples')
    parser_extract_sample.add_argument('--vcf-path', required=True,
                                      help='VCF文件路径|VCF file path')
    parser_extract_sample.add_argument('--spec-samples', required=True, nargs='+',
                                      help='指定样本|Specify samples')
    parser_extract_sample.add_argument('--swave-path',
                                      default='~/software/swave/Swave-main',
                                      help='Swave软件路径|Swave software path')
    parser_extract_sample.add_argument('--output-path',
                                      help='输出路径|Output path')

    # ========== pav 子命令 ==========
    parser_pav = subparsers.add_parser('pav',
                                       help='提取PAV矩阵|Extract PAV (Presence/Absence Variation) matrix')
    parser_pav.add_argument('-i', '--vcf-file', required=True,
                            help='swave converted VCF文件|swave converted VCF file')
    parser_pav.add_argument('-o', '--output-file', default='pav_matrix.tsv',
                            help='输出TSV文件|Output TSV file')
    parser_pav.add_argument('--min-ac', type=int, default=1,
                            help='最小等位基因数|Minimum allele count (default: 1)')
    parser_pav.add_argument('--no-strip-prefix', action='store_true',
                            help='保留CHROM中的样本前缀|Keep sample prefix in CHROM')
    parser_pav.add_argument('--svtype',
                            nargs='+',
                            help='仅保留指定SV类型|Only keep specified SV types (e.g., DUP INS DEL)')

    # 解析参数|Parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        if args.command == 'call':
            # 创建运行器并运行|Create runner and run
            runner = SwaveRunner(command='call',
                                assemblies_tsv=args.assemblies_tsv,
                                ref_fasta=args.ref_fasta,
                                gfa_file=args.gfa_file,
                                gfa_source=args.gfa_source,
                                swave_path=args.swave_path,
                                output_dir=args.output_dir,
                                decomposed_vcf=args.decomposed_vcf,
                                output_mode=args.output_mode,
                                spec_samples=args.spec_samples,
                                min_sv_size=args.min_sv_size,
                                max_sv_size=args.max_sv_size,
                                max_sv_comps=args.max_sv_comps,
                                dup_to_ins=args.dup_to_ins,
                                remove_small=args.remove_small,
                                force_reverse=args.force_reverse,
                                threads=args.threads,
                                minigraph_path=args.minigraph_path,
                                gfatools_path=args.gfatools_path,
                                spec_snarl=args.spec_snarl,
                                spec_path=args.spec_path)

            success = runner.run_call()
            sys.exit(0 if success else 1)

        elif args.command == 'pav':
            # PAV矩阵提取|PAV matrix extraction
            extractor = PAVExtractor()
            vcf_file = os.path.expanduser(args.vcf_file)
            output_file = os.path.expanduser(args.output_file)

            try:
                result = extractor.extract(
                    vcf_file=vcf_file,
                    output_file=output_file,
                    min_ac=args.min_ac,
                    strip_prefix=not args.no_strip_prefix,
                    svtype_only=args.svtype
                )
                print(f"PAV矩阵已保存|PAV matrix saved: {result}")
                sys.exit(0)
            except Exception as e:
                print(f"错误|Error: {e}", file=sys.stderr)
                sys.exit(1)

        else:
            # 其他子命令直接调用Swave.py|Other subcommands call Swave.py directly
            swave_path = getattr(args, 'swave_path', '~/software/swave/Swave-main')
            swave_path = os.path.expanduser(swave_path)
            swave_script = os.path.join(swave_path, 'Swave.py')

            python_path, env = SwaveSVCaller._get_swave_env()
            cmd_list = [python_path, swave_script, args.command]

            # 添加参数|Add arguments
            for arg_name, arg_value in vars(args).items():
                if arg_name in ['command', 'swave_path']:
                    continue
                if arg_value is None or arg_value is False:
                    continue
                if arg_value is True:
                    cmd_list.append(f'--{arg_name}')
                elif isinstance(arg_value, list):
                    cmd_list.extend([f'--{arg_name}'] + arg_value)
                else:
                    cmd_list.extend([f'--{arg_name}', str(arg_value)])

            # 执行命令|Execute command
            result = subprocess.run(cmd_list, shell=False, check=False, env=env)
            sys.exit(result.returncode)

    except Exception as e:
        print(f"错误|Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
