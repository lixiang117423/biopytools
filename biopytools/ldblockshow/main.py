"""
连锁不平衡热图分析主程序模块|LD Heatmap Analysis Main Module
"""

import argparse
import shutil
import sys
import os
from typing import List

from .config import LDBlockShowConfig
from .utils import LDBlockShowLogger, CommandRunner


class LDBlockShowAnalyzer:
    """连锁不平衡热图分析主类|Main LD Heatmap Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = LDBlockShowConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = LDBlockShowLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

    def build_command(self) -> List[str]:
        """构建LDBlockShow命令|Build LDBlockShow command

        Returns:
            List[str]: 命令列表|Command list
        """
        cmd: List[str] = []

        # 可执行文件|Executable
        cmd.append(self.config.ldblockshow_path)

        # 输入文件（三选一）|Input file (one of three)
        if self.config.vcf_file:
            cmd.extend(["-InVCF", self.config.vcf_file])
        elif self.config.in_genotype:
            cmd.extend(["-InGenotype", self.config.in_genotype])
        elif self.config.in_plink:
            cmd.extend(["-InPlink", self.config.in_plink])

        # 必需参数|Required parameters
        cmd.extend(["-OutPut", self.config.output_prefix])
        cmd.extend(["-Region", self.config.region])

        # LD度量选择|LD statistic selection
        if self.config.sele_var != 1:  # 1是默认值|1 is default
            cmd.extend(["-SeleVar", str(self.config.sele_var)])

        # 过滤参数|Filter parameters
        if self.config.maf != 0.05:
            cmd.extend(["-MAF", str(self.config.maf)])

        if self.config.miss != 0.25:
            cmd.extend(["-Miss", str(self.config.miss)])

        if self.config.hwe != 0.0:
            cmd.extend(["-HWE", str(self.config.hwe)])

        if self.config.het != 1.0:
            cmd.extend(["-Het", str(self.config.het)])

        if self.config.enable_oth_var:
            cmd.append("-EnableOthVar")

        # Block检测参数|Block detection parameters
        if self.config.block_type != 1:  # 1是默认值|1 is default
            cmd.extend(["-BlockType", str(self.config.block_type)])

        if self.config.block_type == 3 and self.config.block_cut != "0.85:0.90":
            cmd.extend(["-BlockCut", self.config.block_cut])

        if self.config.block_type == 4 and self.config.fix_block:
            cmd.extend(["-FixBlock", self.config.fix_block])

        # 可视化参数|Visualization parameters
        if self.config.in_gwas:
            cmd.extend(["-InGWAS", self.config.in_gwas])

        if self.config.in_gff:
            cmd.extend(["-InGFF", self.config.in_gff])

        if self.config.mer_min_snp_num != 50:
            cmd.extend(["-MerMinSNPNum", str(self.config.mer_min_snp_num)])

        # 输出格式|Output format
        cmd.append("-OutPng")
        if self.config.out_pdf:
            cmd.append("-OutPdf")

        # 其他参数|Other parameters
        if self.config.sub_pop:
            cmd.extend(["-SubPop", self.config.sub_pop])

        if self.config.tag_snp_cut != 0.80:
            cmd.extend(["-TagSNPCut", str(self.config.tag_snp_cut)])

        # LDBlockShow是编译的独立二进制，不需要conda run包装
        # LDBlockShow is a compiled standalone binary, no conda run wrapping needed
        return cmd

    def run_analysis(self):
        """运行连锁不平衡热图分析|Run LD heatmap analysis"""
        try:
            self.logger.info("连锁不平衡热图分析流程开始|LD heatmap analysis pipeline started")
            self.logger.info(f"分析区域|Region: {self.config.region}")
            self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
            self.logger.info(f"LD度量|LD statistic: {self._get_sele_var_name()}")
            self.logger.info(f"Block类型|Block type: {self._get_block_type_name()}")

            # 记录输入文件|Log input file
            if self.config.vcf_file:
                self.logger.info(f"VCF文件|VCF file: {self.config.vcf_file}")
            elif self.config.in_genotype:
                self.logger.info(f"Genotype文件|Genotype file: {self.config.in_genotype}")
            elif self.config.in_plink:
                self.logger.info(f"Plink文件|Plink file: {self.config.in_plink}")

            # 步骤1: LDBlockShow LD计算|Step 1: LDBlockShow LD calculation
            ld_done = self._run_step("LDBlockShow LD计算|LDBlockShow LD calculation",
                                     f"{self.config.output_prefix}.site.gz",
                                     self._run_ldblockshow)

            # 步骤2: ShowLDSVG绘图|Step 2: ShowLDSVG figure generation
            if ld_done:
                svg_key = f"{self.config.output_prefix}.svg"
                self._run_step("ShowLDSVG绘图|ShowLDSVG figure generation",
                              svg_key if self.config.out_pdf else None,
                              self._run_showldsvg)

            # 检查输出文件|Check output files
            self._check_output_files()

            self.logger.info("=" * 60)
            self.logger.info("分析完成！|Analysis completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"输出文件位于|Output files located at: {self.config.output_prefix}")

            # 列出输出文件|List output files
            output_files = [
                f"{self.config.output_prefix}.svg",
                f"{self.config.output_prefix}.png" if self.config.out_png else None,
                f"{self.config.output_prefix}.pdf" if self.config.out_pdf else None,
                f"{self.config.output_prefix}.site.gz",
                f"{self.config.output_prefix}.blocks.gz",
                f"{self.config.output_prefix}.TriangleV.gz"
            ]

            for file in output_files:
                if file and os.path.exists(file):
                    self.logger.info(f"  - {file}")

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}", exc_info=True)
            sys.exit(1)

    def _is_step_completed(self, output_file: str) -> bool:
        """检查步骤是否已完成|Check if step is done"""
        if not output_file:
            return False
        return os.path.exists(output_file) and os.path.getsize(output_file) > 0

    def _run_step(self, step_name: str, output_file: str, run_func) -> bool:
        """带断点续传的步骤执行器|Step runner with checkpoint resume"""
        if output_file and self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: {step_name}")
            return True

        self.logger.info(f"开始步骤|Starting step: {step_name}")
        return run_func()

    def _run_ldblockshow(self) -> bool:
        """执行LDBlockShow LD计算|Execute LDBlockShow LD calculation"""
        cmd = self.build_command()

        self.logger.info("=" * 60)
        self.logger.info("执行LDBlockShow分析|Executing LDBlockShow analysis")
        self.logger.info("=" * 60)

        success = self.cmd_runner.run(cmd, "LDBlockShow分析|LDBlockShow analysis", timeout=86400)

        if not success:
            self.logger.error("LDBlockShow分析失败|LDBlockShow analysis failed")
            sys.exit(1)

        # 检查是否有SNP数据|Check if SNP data exists
        site_gz = f"{self.config.output_prefix}.site.gz"
        if not self._is_step_completed(site_gz):
            self.logger.warning("无SNP位点数据，跳过绘图|No SNP site data, skipping figure generation")
            return False

        return True

    def _run_showldsvg(self):
        """调用ShowLDSVG生成图像|Call ShowLDSVG to generate figure

        LDBlockShow内部通过std::system调用ShowLDSVG时使用tab分隔参数，
        导致shell无法正确解析。此处由wrapper直接调用ShowLDSVG。
        """
        showldsvg_path = os.path.join(
            os.path.dirname(self.config.ldblockshow_path), 'ShowLDSVG'
        )

        if not os.path.exists(showldsvg_path):
            self.logger.warning(f"ShowLDSVG未找到|ShowLDSVG not found: {showldsvg_path}")
            return

        # 构建ShowLDSVG命令|Build ShowLDSVG command
        cmd: List[str] = [showldsvg_path]
        cmd.extend(["-InPreFix", self.config.output_prefix])
        cmd.extend(["-OutPut", f"{self.config.output_prefix}.svg"])

        if self.config.in_gwas:
            cmd.extend(["-InGWAS", self.config.in_gwas])

        if self.config.in_gff:
            cmd.extend(["-InGFF", self.config.in_gff])

        if self.config.mer_min_snp_num != 50:
            cmd.extend(["-MerMinSNPNum", str(self.config.mer_min_snp_num)])

        if self.config.out_png:
            cmd.append("-OutPng")

        if self.config.out_pdf:
            cmd.append("-OutPdf")

        if self.config.cutline:
            cmd.extend(["-Cutline", str(self.config.cutline)])

        if self.config.point_size is not None:
            cmd.extend(["-PointSize", str(self.config.point_size)])

        if self.config.top_site:
            cmd.extend(["-TopSite", self.config.top_site])

        if self.config.no_log_p:
            cmd.append("-NoLogP")

        if self.config.no_gene_name:
            cmd.append("-NoGeneName")

        if self.config.show_num:
            cmd.append("-ShowNum")

        if self.config.spe_snp_name:
            cmd.extend(["-SpeSNPName", self.config.spe_snp_name])

        if self.config.show_gwas_spe_snp:
            cmd.append("-ShowGWASSpeSNP")

        if self.config.resize_h is not None:
            cmd.extend(["-ResizeH", str(self.config.resize_h)])

        if self.config.no_show_ldist is not None:
            cmd.extend(["-NoShowLDist", str(self.config.no_show_ldist)])

        self.logger.info("执行ShowLDSVG绘图|Executing ShowLDSVG figure generation")
        success = self.cmd_runner.run(cmd, "ShowLDSVG绘图|ShowLDSVG figure generation", timeout=3600)

        if not success:
            self.logger.error("ShowLDSVG绘图失败|ShowLDSVG figure generation failed")
            sys.exit(1)

        # PNG兜底转换：ShowLDSVG内置的batik对大SVG可能内存不足
        # PNG fallback: ShowLDSVG's built-in batik may OOM on large SVGs
        if self.config.out_png:
            self._convert_svg_to_png()

    def _convert_svg_to_png(self):
        """SVG转PNG兜底，当ShowLDSVG内置转换失败时使用其他工具|SVG to PNG fallback"""
        svg_file = f"{self.config.output_prefix}.svg"
        png_file = f"{self.config.output_prefix}.png"

        if not os.path.exists(svg_file):
            return
        if os.path.exists(png_file) and os.path.getsize(png_file) > 0:
            return

        self.logger.info("ShowLDSVG内置PNG转换可能失败，尝试兜底转换|ShowLDSVG built-in PNG conversion may have failed, trying fallback")

        # 尝试rsvg-convert|Try rsvg-convert
        rsvg = shutil.which('rsvg-convert')
        if rsvg:
            self.logger.info(f"使用rsvg-convert转换|Converting with rsvg-convert")
            success = self.cmd_runner.run(
                [rsvg, '-z', '3', svg_file, '-o', png_file],
                "SVG转PNG(rsvg-convert)|SVG to PNG (rsvg-convert)",
                timeout=3600
            )
            if success and os.path.exists(png_file) and os.path.getsize(png_file) > 0:
                self.logger.info(f"PNG转换成功|PNG conversion done: {png_file}")
                return

        # 尝试ImageMagick convert|Try ImageMagick convert
        convert = shutil.which('convert')
        if convert:
            self.logger.info(f"使用ImageMagick convert转换|Converting with ImageMagick convert")
            success = self.cmd_runner.run(
                [convert, svg_file, png_file],
                "SVG转PNG(ImageMagick)|SVG to PNG (ImageMagick)",
                timeout=3600
            )
            if success and os.path.exists(png_file) and os.path.getsize(png_file) > 0:
                self.logger.info(f"PNG转换成功|PNG conversion done: {png_file}")
                return

        self.logger.warning(f"PNG转换失败，请手动转换: {svg_file}|PNG conversion failed, please convert manually: {svg_file}")

    def _check_output_files(self):
        """检查输出文件是否存在|Check if output files exist"""
        expected_files = [
            ("SVG图像|SVG image", f"{self.config.output_prefix}.svg"),
            ("SNP位点文件|SNP sites file", f"{self.config.output_prefix}.site.gz"),
            ("Block文件|Block file", f"{self.config.output_prefix}.blocks.gz"),
            ("配对LD值文件|Pairwise LD file", f"{self.config.output_prefix}.TriangleV.gz"),
        ]

        if self.config.out_png:
            expected_files.append(("PNG图像|PNG image", f"{self.config.output_prefix}.png"))

        if self.config.out_pdf:
            expected_files.append(("PDF图像|PDF image", f"{self.config.output_prefix}.pdf"))

        for desc, file_path in expected_files:
            if os.path.exists(file_path):
                self.logger.info(f"输出文件已生成|Output file generated: {file_path}")
            else:
                self.logger.warning(f"输出文件未找到|Output file not found: {file_path}")

    def _get_sele_var_name(self) -> str:
        """获取LD度量名称|Get LD statistic name"""
        names = {
            1: "D'|D'",
            2: "R²|R²",
            3: "R² and D'|R² and D'",
            4: "D' and R²|D' and R²"
        }
        return names.get(self.config.sele_var, "Unknown")

    def _get_block_type_name(self) -> str:
        """获取Block类型名称|Get block type name"""
        names = {
            1: "Gabriel et al. (PLINK)|Gabriel et al. (PLINK)",
            2: "Solid Spine of LD|Solid Spine of LD",
            3: "BlockCut with self-defined cutoff|BlockCut with self-defined cutoff",
            4: "FixBlock by input file|FixBlock by input file",
            5: "No Block|No Block"
        }
        return names.get(self.config.block_type, "Unknown")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="连锁不平衡热图分析工具：基于VCF/Genotype/Plink生成LD热图|LD heatmap analysis tool: Generate LD heatmap from VCF/Genotype/Plink",
        epilog='示例|Examples: %(prog)s -i variants.vcf.gz -o ld_result -r chr1:1000000-2000000'
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')

    required.add_argument("-o", "--output-prefix", required=True,
                         help="输出文件前缀（包含路径）|Output file prefix (including path)")
    required.add_argument("-r", "--region", required=True,
                         help="分析区域，格式chr:start-end|Analysis region, format chr:start-end")

    # 输入文件（三选一）|Input files (one of three)
    input_params = parser.add_argument_group(
        '输入文件（三选一）|Input files (one required)'
    )

    input_params.add_argument("-i", "--vcf-file",
                             help="VCF变异文件路径|VCF variant file path")
    input_params.add_argument("--in-genotype",
                             help="SNP Genotype格式文件路径|SNP Genotype format file path")
    input_params.add_argument("--in-plink",
                             help="Plink文件前缀(bed+bim+fam或ped+map)|Plink file prefix (bed+bim+fam or ped+map)")

    # LD度量选择|LD statistic selection
    ld_stat = parser.add_argument_group('LD度量选择|LD statistic selection')

    ld_stat.add_argument("--sele-var", type=int, choices=[1, 2, 3, 4], default=1,
                        help='选择LD度量统计量|Select LD statistic: 1=D\' (default), 2=R², 3/4=Both')

    # 过滤参数|Filter parameters
    filter_params = parser.add_argument_group('过滤参数|Filter parameters')

    filter_params.add_argument("--maf", type=float, default=0.05,
                             help="最小次要等位基因频率|Minimum minor allele frequency (default: 0.05)")
    filter_params.add_argument("--miss", type=float, default=0.25,
                             help="最大缺失率|Maximum missing ratio (default: 0.25)")
    filter_params.add_argument("--hwe", type=float, default=0.0,
                             help="Hardy-Weinberg平衡P值阈值|Hardy-Weinberg equilibrium P-value threshold (default: 0.0)")
    filter_params.add_argument("--het", type=float, default=1.0,
                             help="最大杂合率|Maximum heterozygosity ratio (default: 1.0)")
    filter_params.add_argument("--enable-oth-var", action="store_true",
                             help="允许indel/SV/CNV变异|Allow bi-indel bi-sv bi-cnv variants")

    # Block检测参数|Block detection parameters
    block_params = parser.add_argument_group('Block检测参数|Block detection parameters')

    block_params.add_argument("--block-type", type=int, choices=[1, 2, 3, 4, 5], default=1,
                             help='Block检测方法|Block detection method: 1=Gabriel(default), 2=SolidSpine, 3=BlockCut, 4=FixBlock, 5=NoBlock')
    block_params.add_argument("--block-cut", default="0.85:0.90",
                             help="BlockType3的cutoff（格式：cutoff:ratio）|Cutoff for BlockType3 (format: cutoff:ratio)")
    block_params.add_argument("--fix-block",
                             help="固定block文件路径|Fixed block file path (for BlockType=4)")

    # 可视化参数|Visualization parameters
    vis_params = parser.add_argument_group('可视化参数|Visualization parameters')

    vis_params.add_argument("--in-gwas",
                          help="GWAS P值文件，制表符/空格分隔，3列: chr pos pvalue，染色体名须与VCF一致|GWAS P-value file, tab/space delimited, 3 columns: chr pos pvalue, chr names must match VCF")
    vis_params.add_argument("--in-gff",
                          help="GFF3注释文件路径|GFF3 annotation file path")
    vis_params.add_argument("--mer-min-snp-num", type=int, default=50,
                          help="合并网格的最小SNP数|Minimum SNP number to merge grids (default: 50)")

    # 输出格式|Output format
    output_format = parser.add_argument_group('输出格式|Output format')

    output_format.add_argument("--no-out-png", action="store_true",
                             help="不输出PNG格式图像|Do not output PNG format image")
    output_format.add_argument("--out-pdf", action="store_true",
                             help="输出PDF格式图像|Output PDF format image")

    # 其他参数|Other parameters
    other_params = parser.add_argument_group('其他参数|Other parameters')

    other_params.add_argument("--sub-pop",
                             help="亚群样本文件路径|Subgroup sample file path")
    other_params.add_argument("--tag-snp-cut", type=float, default=0.80,
                             help="TagSNP的LD cutoff|LD cutoff for TagSNP (default: 0.80)")

    # ShowLDSVG绘图参数|ShowLDSVG drawing parameters
    svg_params = parser.add_argument_group(
        '绘图参数(ShowLDSVG)|Drawing parameters (ShowLDSVG)'
    )

    svg_params.add_argument("--cutline", type=float, default=5.0,
                           help="GWAS P值显著性阈值(-log10)|GWAS P-value significance cutoff (-log10) (default: 5)")
    svg_params.add_argument("--point-size", type=int,
                           help="GWAS散点大小|GWAS point size")
    svg_params.add_argument("--top-site",
                           help="指定GWAS峰值位点(chr:pos)|Specify GWAS peak site (chr:pos)")
    svg_params.add_argument("--no-log-p", action="store_true",
                           help="不对P值取-log10|Do not -log10 transform P-value")
    svg_params.add_argument("--no-gene-name", action="store_true",
                           help="不显示基因名|Do not show gene names")
    svg_params.add_argument("--show-num", action="store_true",
                           help="在热图中显示R²/D'值|Show R²/D' values in heatmap")
    svg_params.add_argument("--spe-snp-name",
                           help="特殊SNP名称文件(chr site Name)|Special SNP name file (chr site Name)")
    svg_params.add_argument("--show-gwas-spe-snp", action="store_true",
                           help="在GWAS图中显示特殊SNP名称|Show special SNP names in GWAS plot")
    svg_params.add_argument("--resize-h", type=int,
                           help="图像高度，宽度按比例自动调整|Image height, width auto-adjusted")
    svg_params.add_argument("--no-show-ldist", type=int,
                           help="超过此距离的SNP对不显示LD|NoShow pairwise LD over this distance")

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = LDBlockShowAnalyzer(
            output_prefix=args.output_prefix,
            region=args.region,
            vcf_file=args.vcf_file,
            in_genotype=args.in_genotype,
            in_plink=args.in_plink,
            sele_var=args.sele_var,
            maf=args.maf,
            miss=args.miss,
            hwe=args.hwe,
            het=args.het,
            enable_oth_var=args.enable_oth_var,
            block_type=args.block_type,
            block_cut=args.block_cut,
            fix_block=args.fix_block,
            in_gwas=args.in_gwas,
            in_gff=args.in_gff,
            mer_min_snp_num=args.mer_min_snp_num,
            out_png=not args.no_out_png,
            out_pdf=args.out_pdf,
            sub_pop=args.sub_pop,
            tag_snp_cut=args.tag_snp_cut,
            cutline=args.cutline,
            point_size=args.point_size,
            top_site=args.top_site,
            no_log_p=args.no_log_p,
            no_gene_name=args.no_gene_name,
            show_num=args.show_num,
            spe_snp_name=args.spe_snp_name,
            show_gwas_spe_snp=args.show_gwas_spe_snp,
            resize_h=args.resize_h,
            no_show_ldist=args.no_show_ldist
        )

        analyzer.run_analysis()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("操作被用户中断|Operation interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
