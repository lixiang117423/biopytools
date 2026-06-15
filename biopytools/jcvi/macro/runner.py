"""
宏观共线性可视化运行器|Macrosynteny Visualization Runner

支持两物种和多物种(3+)的karyotype共线性图自动生成
Supports both two-species and multi-species (3+) karyotype collinearity visualization

全流程|Full pipeline:
1-6. 共享管道 (JcviPipeline): 发现样本→提取序列→GFF→BED→比对→过滤→synteny scan
7.   jcvi.compara.synteny screen → .anchors.simple
8.   自动生成 seqids (从BED提取染色体, 过滤小scaffold)
9.   自动生成 layout (自动排版, 自动配色, 相邻物种edges)
10.  jcvi.graphics.karyotype → PDF

支持 --replot: 跳过步骤1-9, 仅重新执行最后绘图步骤
"""

import os
import time
from pathlib import Path
from datetime import datetime
from typing import Optional, List
from collections import defaultdict

from .config import MacroConfig
from ..config import JcviBaseConfig
from ..utils import (
    JcviLogger,
    build_jcvi_command,
    get_jcvi_stem,
)
from ..pipeline import JcviPipeline

# 多物种默认调色板(与matplotlib一致)|Multi-species default color palette
DEFAULT_COLORS = ['r', 'b', 'g', 'm', 'c', 'y', 'k', 'orange',
                  'purple', 'brown', 'pink', 'gray']


class MacroRunner:
    """宏观共线性可视化运行器|Macrosynteny Visualization Runner"""

    VERSION = "2.0.0"

    def __init__(self, config: MacroConfig, logger: Optional[JcviLogger] = None):
        self.config = config
        self.logger_obj = logger
        self.logger = logger.get_logger() if logger else None
        self.start_time = None

    def run(self) -> bool:
        try:
            self.start_time = time.time()

            if self.logger_obj is None:
                log_dir = Path(self.config.output_dir) / "99_logs"
                log_dir.mkdir(parents=True, exist_ok=True)
                self.logger_obj = JcviLogger(log_dir / "macro.log", "MacroSynteny")
                self.logger = self.logger_obj.get_logger()

            self._print_header()
            self.config.validate()

            species_list = self.config.species_list
            n_species = len(species_list)
            mode = "多物种" if n_species >= 3 else "两物种"

            self.logger.info(
                f"模式|Mode: {mode} ({n_species}物种|species): "
                f"{', '.join(species_list)}"
            )

            # 构建物种对列表(相邻物种)
            pair_list = [(species_list[i], species_list[i + 1])
                         for i in range(n_species - 1)]
            # 转为 pipeline 需要的格式: ["A,B"]
            pair_strs = [f"{a},{b}" for a, b in pair_list]

            if self.config.replot:
                self.logger.info(
                    "  --replot 模式, 跳过管道和screen/seqids/layout生成|"
                    "--replot mode, skipping pipeline and screen/seqids/layout"
                )
                seqids_file = Path(self.config.output_dir) / "seqids"
                layout_file = Path(self.config.output_dir) / "layout"
                if not seqids_file.exists() or not layout_file.exists():
                    self.logger.error(
                        "未找到seqids/layout文件, 请先运行不带--replot的命令|"
                        "seqids/layout not found, run without --replot first"
                    )
                    return False
            else:
                # 步骤1-6: 共享管道 (发现样本→提取序列→BED→比对→过滤→synteny)
                self._log_section(
                    f"步骤1-6/10: 共享管道|Steps 1-6/10: Shared pipeline"
                )
                pipeline_config = self._build_pipeline_config(pair_strs)
                pipeline = JcviPipeline(pipeline_config, self.logger_obj)
                pipeline_result = pipeline.run()

                if not pipeline_result.pair_dirs:
                    self.logger.error("共享管道未产出任何pairwise结果|"
                                      "Shared pipeline produced no pairwise results")
                    return False

                # 查找BED文件 (从pipeline输出)
                self.logger.info("")
                bed_map = {}
                for sp in species_list:
                    bed = self._find_bed(sp)
                    if not bed:
                        self.logger.error(
                            f"未找到物种 {sp} 的BED文件|BED not found for species: {sp}"
                        )
                        return False
                    bed_map[sp] = bed

                # 步骤7: Screen anchors → .anchors.simple
                self.logger.info("")
                self._log_section(
                    f"步骤7/10: Screen anchors ({n_species - 1}对|pairs)|"
                    f"Step 7/10: Screen anchors"
                )
                simple_files = {}
                for sp_a, sp_b in pair_list:
                    stem_a = get_jcvi_stem(sp_a)
                    stem_b = get_jcvi_stem(sp_b)
                    pprefix = f"{stem_a}.{stem_b}"

                    pair_dir = pipeline_result.pair_dirs.get((sp_a, sp_b))
                    if not pair_dir:
                        # 尝试反向查找 (pipeline可能生成反方向目录)
                        pair_dir = pipeline_result.pair_dirs.get((sp_b, sp_a))
                    if not pair_dir:
                        self.logger.error(
                            f"未找到 {sp_a} vs {sp_b} 的pairwise目录|"
                            f"Pairwise dir not found: {sp_a} vs {sp_b}"
                        )
                        return False

                    simple = self._screen_anchors(pair_dir, pprefix)
                    if not simple:
                        return False
                    simple_files[(sp_a, sp_b)] = simple

                # 步骤8: 生成 seqids
                self.logger.info("")
                self._log_section("步骤8/10: 生成seqids|Step 8/10: Generate seqids")
                seqids_file = self._generate_seqids(species_list, bed_map)
                self.logger.info(f"  生成|Generated: {seqids_file}")

                # 步骤9: 生成 layout
                self.logger.info("")
                self._log_section("步骤9/10: 生成layout|Step 9/10: Generate layout")
                layout_file = self._generate_layout(
                    species_list, bed_map, simple_files
                )
                self.logger.info(f"  生成|Generated: {layout_file}")

            # 步骤10: 绘制 karyotype
            self.logger.info("")
            self._log_section("步骤10/10: 绘制karyotype|Step 10/10: Plot karyotype")
            pdf_file = self._plot_karyotype(seqids_file, layout_file)

            if pdf_file and pdf_file.exists():
                self.logger.info(f"\n  输出文件|Output: {pdf_file}")
                self.logger.info(
                    "\n  提示|Tip: 可手动编辑seqids和layout文件后使用 --replot 重新绘图"
                )
                self.logger.info(
                    "  Tip: Edit seqids/layout files and re-run with --replot"
                )

            self._print_footer(True)
            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Error: {e}", exc_info=True)
            return False

    # ---- Pipeline 配置构建 ----

    def _build_pipeline_config(self, pair_strs: List[str]) -> JcviBaseConfig:
        """构建共享管道配置|Build shared pipeline config from MacroConfig"""
        return JcviBaseConfig(
            input_dir=self.config.input_dir,
            output_dir=self.config.output_dir,
            conda_env=self.config.conda_env,
            gff_ext=self.config.gff_ext,
            fa_ext=self.config.fa_ext,
            gff_type=self.config.gff_type,
            gff_key=self.config.gff_key,
            cscore=self.config.cscore,
            min_size=self.config.min_size,
            dbtype=self.config.dbtype,
            align_soft=self.config.align_soft,
            no_strip_names=self.config.no_strip_names,
            no_dotplot=self.config.no_dotplot,
            threads=self.config.threads,
            pairs=pair_strs,
        )

    # ---- 文件查找 ----

    def _find_bed(self, name: str) -> Optional[Path]:
        """查找BED文件 (从pipeline输出目录)|Find BED file from pipeline output"""
        output_dir = Path(self.config.output_dir)
        candidates = [
            output_dir / "02_bed" / f"{get_jcvi_stem(name)}.uniq.bed",
            output_dir / "02_bed" / f"{name}.uniq.bed",
            output_dir / "bed" / f"{name}.uniq.bed",
        ]
        for f in candidates:
            if f.exists() and f.stat().st_size > 0:
                return f
        return None

    # ---- 步骤7: Screen anchors ----

    def _screen_anchors(self, pair_dir: Path, pprefix: str) -> Optional[Path]:
        """步骤7: jcvi.compara.synteny screen → .anchors.simple"""
        anchors_file = pair_dir / f"{pprefix}.anchors"
        simple_out = Path(self.config.output_dir) / f"{pprefix}.anchors.simple"

        if simple_out.exists() and simple_out.stat().st_size > 0:
            self.logger.info(f"  跳过已完成|Skipping completed: screen ({pprefix})")
            return simple_out

        if not anchors_file.exists():
            self.logger.error(f"  未找到anchors文件|Anchors not found: {anchors_file}")
            return None

        cmd = build_jcvi_command('jcvi.compara.synteny', [
            'screen',
            f"--minspan={self.config.minspan}",
            '--simple',
            str(anchors_file),
            str(simple_out),
        ], self.config.conda_env)

        self.logger.info(f"  执行|Executing: screen ({pprefix})")
        self.logger.info(f"  命令|Command: {' '.join(cmd)}")
        result = self._run_cmd(cmd, cwd=str(pair_dir))
        if result != 0:
            self.logger.error(f"  screen anchors 失败|screen failed: {pprefix}")
            return None

        return simple_out

    # ---- 步骤8: 生成 seqids ----

    def _generate_seqids(self, species_list: List[str],
                         bed_map: dict) -> Path:
        """步骤8: 从BED文件为所有物种生成seqids|Generate seqids for all species"""
        seqids_file = Path(self.config.output_dir) / "seqids"

        with open(seqids_file, 'w') as f:
            for sp in species_list:
                chroms = self._extract_chroms(bed_map[sp])
                if chroms:
                    f.write(','.join(chroms) + '\n')
                else:
                    self.logger.warning(
                        f"  物种 {sp} 未找到符合条件的染色体|"
                        f"No qualifying chromosomes for species: {sp}"
                    )

        return seqids_file

    def _extract_chroms(self, bed_file: Path) -> list:
        """从BED文件提取染色体列表, 过滤小scaffold|Extract chromosomes from BED"""
        chrom_genes = defaultdict(int)
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    chrom_genes[parts[0]] += 1

        return [c for c, n in sorted(chrom_genes.items())
                if n >= self.config.min_chr_genes]

    # ---- 步骤9: 生成 layout ----

    def _generate_layout(self, species_list: List[str], bed_map: dict,
                         simple_files: dict) -> Path:
        """步骤9: 自动生成layout|Auto-generate layout file"""
        layout_file = Path(self.config.output_dir) / "layout"
        n = len(species_list)

        with open(layout_file, 'w') as f:
            f.write("# y, xstart, xend, rotation, color, label, va, bed, align\n")

            y_positions = self._calc_y_positions(n)

            for i, sp in enumerate(species_list):
                y = y_positions[i]
                color = DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                va = "bottom" if i == n - 1 else "top"
                xstart, xend = 0.18, 0.95

                f.write(f"{y}, {xstart}, {xend}, 0, {color}, "
                        f"{sp}, {va}, {bed_map[sp]}, center\n")

            f.write("# edges\n")
            for i in range(n - 1):
                sp_a = species_list[i]
                sp_b = species_list[i + 1]
                simple = simple_files[(sp_a, sp_b)]
                f.write(f"e, {i}, {i + 1}, {simple}\n")

        return layout_file

    @staticmethod
    def _calc_y_positions(n: int) -> list:
        """计算n个物种的y坐标|Calculate y positions for n species"""
        if n == 2:
            return [0.6, 0.4]

        y_top = 0.88
        y_bottom = 0.12
        step = (y_top - y_bottom) / (n - 1)
        return [round(y_top - i * step, 3) for i in range(n)]

    # ---- 步骤10: 绘图 ----

    def _plot_karyotype(self, seqids_file: Path, layout_file: Path) -> Optional[Path]:
        """步骤10: jcvi.graphics.karyotype"""
        pdf_file = Path(self.config.output_dir) / "karyotype.pdf"

        args = [str(seqids_file), str(layout_file)]

        n_species = len(self.config.species_list)
        if not self.config.figsize:
            if n_species >= 3:
                self.config.figsize = f"{14}x{max(10, n_species * 2)}"
        if self.config.figsize:
            args.append(f"--figsize={self.config.figsize}")

        args.append(f"--shadestyle={self.config.shadestyle}")
        args.append(f"--chrstyle={self.config.chrstyle}")
        args.append("--notex")

        cmd = build_jcvi_command('jcvi.graphics.karyotype', args,
                                 self.config.conda_env)

        self.logger.info(f"  命令|Command: {' '.join(cmd)}")
        result = self._run_cmd(cmd, cwd=str(self.config.output_dir))
        if result != 0:
            self.logger.error("  karyotype绘图失败|karyotype plotting failed")
            return None

        return pdf_file

    # ---- 工具方法 ----

    def _run_cmd(self, cmd, cwd=None):
        import subprocess
        run_env = None
        if cmd and '/envs/' in cmd[0] and cmd[0].endswith('python'):
            env_bin = os.path.dirname(cmd[0])
            run_env = os.environ.copy()
            run_env['PATH'] = env_bin + os.pathsep + run_env.get('PATH', '')
        result = subprocess.run(cmd, cwd=cwd, env=run_env)
        return result.returncode

    def _print_header(self):
        if not self.logger:
            return
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("JCVI 宏观共线性可视化|JCVI Macrosynteny Visualization")
        self.logger.info("=" * 60)
        self.logger.info(f"版本|Version: {self.VERSION}")
        self.logger.info(f"时间|Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")

    def _print_footer(self, success: bool):
        if not self.logger:
            return
        elapsed = time.time() - self.start_time
        hours, remainder = divmod(int(elapsed), 3600)
        minutes, seconds = divmod(remainder, 60)
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(
            "宏观共线性可视化完成|Macrosynteny Visualization Completed"
            if success else "宏观共线性可视化失败|Failed"
        )
        self.logger.info(f"  总耗时|Total time: {hours}h {minutes}m {seconds}s")
        self.logger.info("=" * 60)

    def _log_section(self, title: str):
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(title)
        self.logger.info("=" * 60)
