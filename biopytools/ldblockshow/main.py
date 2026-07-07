"""
连锁不平衡热图分析主程序模块|LD Heatmap Analysis Main Module
"""

import argparse
import glob
import re
import shutil
import sys
import os
from dataclasses import dataclass
from typing import List, Optional

from .config import LDBlockShowConfig
from .utils import LDBlockShowLogger, CommandRunner


@dataclass
class BedRegion:
    """BED 一行 region(0-based 半开)|One BED row (0-based half-open)"""

    chrom: str
    bed_start: int   # 0-based, BED 原始|raw BED start
    bed_end: int      # 0-based 半开, BED 原始|raw BED end (exclusive)
    name: Optional[str] = None

    @property
    def region_str(self) -> str:
        """ldblockshow -Region 用的 1-based 闭区间 chr:start-end
        |1-based inclusive for ldblockshow -Region (aligned to VCF POS)

        BED 0-based 半开 [start,end) → 1-based 闭 [start+1, end]
        |BED 0-based half-open [start,end) → 1-based inclusive [start+1, end]
        """
        return f"{self.chrom}:{self.bed_start + 1}-{self.bed_end}"

    @property
    def label(self) -> str:
        """输出文件命名：name 优先(清洗)，否则 chr_start_end(1-based)
        |Output label: name (sanitized) if present, else 1-based chr_start_end"""
        if self.name:
            return _sanitize_label(self.name)
        return f"{self.chrom}_{self.bed_start + 1}_{self.bed_end}"


# BED name 允许字符：字母数字 . - _，其余替换为 _ |allowed chars; others → _
_LABEL_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _sanitize_label(name: str) -> str:
    """清洗 BED name 用于输出文件名(去空格/斜杠等)|Sanitize BED name for filename use"""
    return _LABEL_SAFE.sub("_", name)


# region 字符串(chr:start-end)中需替换为 _ 的分隔符|separators in region string → _
_REGION_SEP = re.compile(r"[:\-]")


def _region_str_to_label(region_str: str) -> str:
    """单 region 字符串 chr:start-end → 文件名 stem chr_start_end
    |Single-region string chr:start-end → filename stem chr_start_end"""
    return _REGION_SEP.sub("_", region_str)


def _find_rsvg_convert() -> Optional[str]:
    """
    定位 rsvg-convert：PATH > RSVG_CONVERT_PATH 环境变量 > ~/miniforge3/envs/*/bin
    |Locate rsvg-convert: PATH > RSVG_CONVERT_PATH env > ~/miniforge3/envs/*/bin

    ShowLDSVG 内置的 batik 对大 SVG(数百万 polygon)会渲染出空白 PNG，
    rsvg-convert 是可靠的 PNG 渲染器。PATH 里通常没有，故额外检索已装 librsvg 的 conda 环境。
    |ShowLDSVG's built-in batik renders blank PNGs on large SVGs (millions of polygons);
    rsvg-convert is the reliable renderer. It's usually absent from PATH, so also search
    conda envs that already bundle librsvg.
    """
    # 1. PATH
    found = shutil.which('rsvg-convert')
    if found:
        return found
    # 2. 环境变量|env var
    env_path = os.environ.get('RSVG_CONVERT_PATH')
    if env_path and os.path.exists(env_path):
        return env_path
    # 3. 已安装 librsvg 的 conda 环境(复用现有安装)|reuse existing installs in conda envs
    matches = sorted(glob.glob(os.path.expanduser('~/miniforge3/envs/*/bin/rsvg-convert')))
    if matches:
        return matches[0]
    return None


def _parse_bed_regions(bed_path) -> List[BedRegion]:
    """
    解析基因组 BED 文件为 region 列表|Parse a genomic BED file into regions

    - 接受 BED3+ (chrom, start, end[, name], ...)；第4列作 name(可选)
    - 跳过空行、# 注释、track/browser 头|Skip blank/#/track/browser lines
    - 严格校验：整数坐标、end>start，否则 ValueError(带行号) fail fast
    - 至少需 1 个有效 region|Need >=1 valid region

    Args:
        bed_path: BED 文件路径|BED file path

    Returns:
        region 列表|list of BedRegion
    """
    regions: List[BedRegion] = []
    with open(bed_path, "r") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line:
                continue
            first = line.split()[0]
            # 跳过 UCSC 头/注释|Skip UCSC headers/comments
            if first in ("track", "browser") or first.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"BED 第{lineno}行无效|BED line {lineno} invalid: 至少3列|need >=3 columns")
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise ValueError(
                    f"BED 第{lineno}行无效|BED line {lineno} invalid: 坐标非整数|non-integer coords"
                ) from None
            if end <= start:
                raise ValueError(
                    f"BED 第{lineno}行无效|BED line {lineno} invalid: end({end})<=start({start})")
            name = fields[3] if len(fields) >= 4 else None
            regions.append(BedRegion(chrom, start, end, name))

    if not regions:
        raise ValueError("BED 无有效 region|BED has no valid regions")
    return regions


def _per_region_prefix(output_dir: str, region: BedRegion, seen: set) -> str:
    """
    为单个 region 生成输出目录内的不冲突前缀(目录/label)|Non-colliding prefix inside the output dir (dir/label).

    label 重名时追加 _2, _3...|Append _2, _3... on duplicate labels.

    Args:
        output_dir: --output-dir 给的目录|output directory
        region: 当前 region|current region
        seen: 已使用前缀集合(会被修改)|set of used prefixes (mutated)

    Returns:
        该 region 的完整输出前缀(目录/label)|full output prefix (dir/label)
    """
    label = region.label
    candidate = label
    suffix = 2
    full = os.path.join(output_dir, candidate)
    while full in seen:
        candidate = f"{label}_{suffix}"
        suffix += 1
        full = os.path.join(output_dir, candidate)
    seen.add(full)
    return full


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
        if self.config.no_snp_filter:
            # 强制覆盖LDBlockShow内部默认值，确保不过滤任何SNP
            # Override LDBlockShow internal defaults to disable all SNP filtering
            cmd.extend(["-MAF", "0.0", "-Miss", "1.0"])
        else:
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
        """运行连锁不平衡热图分析(单region入口，完成后退出)|Single-region entry: run then exit"""
        try:
            ok = self._run_pipeline()
        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            sys.exit(130)
        sys.exit(0 if ok else 1)

    def _run_pipeline(self) -> bool:
        """
        执行LD热图全流程，返回是否成功(不调用sys.exit)|Run full pipeline, return success (no sys.exit)

        单region run_analysis 与 batch 多region 循环复用本方法。
        |Shared by single-region run_analysis and the batch multi-region loop.
        KeyboardInterrupt 上抛由调用方决定(单region退出 / batch停止)。
        |KeyboardInterrupt propagates: single-region exits / batch stops.
        """
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

            return True

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            raise  # 上抛由调用方处理(单region退出 / batch停止)|propagate
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}", exc_info=True)
            return False

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
        """SVG转PNG：优先 rsvg-convert(batik 对大 SVG 渲染空白)，兜底 ImageMagick
        |SVG to PNG: prefer rsvg-convert (batik renders large SVGs blank), fallback IM

        ShowLDSVG/LDBlockShow 内置的 batik 对含数百万 polygon 的大 SVG 会渲染出**空白但非零**
        的 PNG，不能用 getsize>0 判断好坏。故只要 rsvg-convert 可用就总是重新渲染覆盖。
        |The built-in batik renders a **blank-but-nonzero** PNG on large SVGs (millions of
        polygons), so getsize>0 cannot judge correctness. Always re-render with rsvg-convert
        when available.
        """
        svg_file = f"{self.config.output_prefix}.svg"
        png_file = f"{self.config.output_prefix}.png"

        if not os.path.exists(svg_file):
            return

        rsvg = _find_rsvg_convert()
        if rsvg:
            self.logger.info(
                f"使用 rsvg-convert 渲染 PNG(避免 batik 大 SVG 空白)|"
                f"Render PNG with rsvg-convert (avoids batik blank on large SVG): {rsvg}")
            success = self.cmd_runner.run(
                [rsvg, '-z', '3', svg_file, '-o', png_file],
                "SVG转PNG(rsvg-convert)|SVG to PNG (rsvg-convert)",
                timeout=7200
            )
            if success and os.path.exists(png_file) and os.path.getsize(png_file) > 0:
                self.logger.info(
                    f"PNG渲染完成|PNG rendered: {png_file} "
                    f"({os.path.getsize(png_file) // 1024} KB)")
                return
            self.logger.warning("rsvg-convert 渲染失败，尝试兜底|rsvg-convert failed, trying fallback")

        # 兜底：ImageMagick convert(其 SVG 委托也依赖 rsvg-convert，缺失时可能空白)
        # |Fallback: ImageMagick convert (its SVG delegate also needs rsvg-convert; blank if missing)
        convert = shutil.which('convert')
        if convert:
            self.logger.info("兜底使用 ImageMagick convert|Fallback: ImageMagick convert")
            success = self.cmd_runner.run(
                [convert, svg_file, png_file],
                "SVG转PNG(ImageMagick)|SVG to PNG (ImageMagick)",
                timeout=7200
            )
            if success and os.path.exists(png_file) and os.path.getsize(png_file) > 0:
                self.logger.info(f"PNG转换完成|PNG converted: {png_file}")
                return

        # 都失败：batik 已写的 PNG 很可能空白，明确告警并给出修复建议
        # |All failed: batik's PNG is likely blank; warn clearly with a fix
        self.logger.error(
            "无法生成有效 PNG：未找到 rsvg-convert，而 batik 对大 SVG 会渲染空白。"
            "请安装 librsvg(如 conda install -c conda-forge librsvg)或设置环境变量 "
            "RSVG_CONVERT_PATH 指向 rsvg-convert。SVG 本身正常，可手动转换。"
            "|Cannot produce a valid PNG: rsvg-convert not found, and batik renders large "
            "SVGs blank. Install librsvg (e.g. conda install -c conda-forge librsvg) or set "
            "RSVG_CONVERT_PATH to point to rsvg-convert. The SVG itself is fine; convert manually.")

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

    required.add_argument("-o", "--output-dir", required=True,
                         help="输出目录(自动创建)；每 region 产物落在 目录/<label>.* "
                         "|Output directory (auto-created); per-region outputs land in dir/<label>.*")

    # 分析区域：-r 与 -b 二选一|Region: exactly one of -r or -b
    region_params = parser.add_argument_group('分析区域（-r 与 -b 二选一）|Region (-r XOR -b)')
    region_params.add_argument("-r", "--region",
                              help="单个分析区域，格式chr:start-end|Single region, format chr:start-end")
    region_params.add_argument("-b", "--bed",
                              help="基因组BED文件(每行 chrom start end [name])，等价多个 -r 批量出图"
                              "|Genomic BED (cols: chrom start end [name]), equivalent to multiple -r")

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

    filter_params.add_argument("--enable-snp-filter", action="store_true", default=False,
                             help="启用SNP过滤（默认关闭），使用-MAF/-Miss/-HWE/-Het参数过滤SNP|Enable SNP filtering (default OFF), filter SNPs using -MAF/-Miss/-HWE/-Het parameters")
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

    # -r/--region 与 -b/--bed 恰好一个|exactly one of -r/--region or -b/--bed
    if bool(args.region) == bool(args.bed):
        print("错误|Error: 必须且只能指定 -r/--region 或 -b/--bed 之一"
              "|Must specify exactly one of -r/--region or -b/--bed", file=sys.stderr)
        sys.exit(1)

    # 创建分析器并运行|Create analyzer(s) and run
    try:
        # 创建输出目录|create output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # 共享参数(region/output_prefix 随 region 变)|shared kwargs (region/output_prefix vary per region)
        shared = dict(
            vcf_file=args.vcf_file,
            in_genotype=args.in_genotype,
            in_plink=args.in_plink,
            sele_var=args.sele_var,
            no_snp_filter=not args.enable_snp_filter,
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
            no_show_ldist=args.no_show_ldist,
        )

        if args.bed:
            # BED 多 region 批处理|BED multi-region batch
            if not os.path.exists(args.bed):
                print(f"错误|Error: BED文件不存在|BED file not found: {args.bed}", file=sys.stderr)
                sys.exit(1)
            regions = _parse_bed_regions(args.bed)
            print(f"BED 解析到|BED parsed: {len(regions)} 个 region|regions", file=sys.stderr)

            seen = set()
            succeeded = 0
            failed = 0
            for idx, region in enumerate(regions, 1):
                per_prefix = _per_region_prefix(args.output_dir, region, seen)
                print(f"[{idx}/{len(regions)}] region={region.region_str} -> {per_prefix}",
                      file=sys.stderr)
                analyzer = LDBlockShowAnalyzer(
                    output_prefix=per_prefix, region=region.region_str, **shared)
                try:
                    ok = analyzer._run_pipeline()
                except KeyboardInterrupt:
                    raise  # 停止整个 batch|stop the whole batch
                if ok:
                    succeeded += 1
                else:
                    failed += 1
                    print(f"[{idx}/{len(regions)}] 失败，继续下一个|failed, continuing",
                          file=sys.stderr)
            print(f"BED 批处理完成|BED batch done: {succeeded} 成功|succeeded, "
                  f"{failed} 失败|failed (共|total {len(regions)})", file=sys.stderr)
            sys.exit(0 if failed == 0 else 1)
        else:
            # 单 region：stem 由 region 字符串派生(chr_start_end)|single region: stem from region string
            stem = _region_str_to_label(args.region)
            output_prefix = os.path.join(args.output_dir, stem)
            analyzer = LDBlockShowAnalyzer(
                output_prefix=output_prefix, region=args.region, **shared)
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
