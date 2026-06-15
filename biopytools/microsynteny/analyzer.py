"""
微观共线性分析核心逻辑|Microsynteny Analysis Core Logic
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple
import itertools

from .config import MicrosyntenyConfig
from .utils import MicrosyntenyLogger, run_jcvi_command, check_file_exists, find_gene_in_bed


class MicrosyntenyAnalyzer:
    """微观共线性分析器|Microsynteny Analyzer"""

    def __init__(self, config: MicrosyntenyConfig, logger):
        """初始化分析器|Initialize analyzer

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # JCVI Python路径|JCVI Python path
        self.jcvi_python = os.path.join(config.jcvi_path, "bin", "python")

    def run_pipeline(self):
        """运行完整分析流程|Run complete analysis pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("开始微观共线性分析流程|Starting microsynteny analysis pipeline")
        self.logger.info("=" * 60)

        # 步骤1: 数据预处理|Step 1: Data preprocessing
        if self.config.step in ['1', None]:
            self.run_step1_preprocess()

        # 步骤2: 共线性分析|Step 2: Synteny analysis
        if self.config.step in ['2', None]:
            self.run_step2_synteny()

        # 步骤3: 区块提取|Step 3: Block extraction
        if self.config.step in ['3', None]:
            self.run_step3_blocks()

        # 步骤4: 绘图|Step 4: Plotting
        if self.config.step in ['4', None]:
            self.run_step4_plot()

        self.logger.info("=" * 60)
        self.logger.info("分析流程完成|Analysis pipeline completed")
        self.logger.info("=" * 60)

    def run_step1_preprocess(self):
        """步骤1: 数据预处理|Step 1: Data preprocessing

        - GFF → BED转换
        - 提取CDS序列
        - 合并BED文件
        """
        step = "1"
        step_name = "数据预处理|Data preprocessing"

        if self.config.is_step_done(step):
            self.logger.info(f"步骤{step}已完成，跳过|Step {step} already completed, skipping")
            return

        self.logger.info(f"步骤{step}: {step_name}|Step {step}: {step_name}")

        species_list = self.config.get_species_list()

        for sp_id in species_list:
            self.logger.info(f"处理物种|Processing species: {sp_id}")

            # 使用存储的文件路径（支持多种后缀）|Use stored file paths (support multiple extensions)
            gff_file = self.config.genome_files[f"{sp_id}_gff"]
            fa_file = self.config.genome_files[f"{sp_id}_fa"]

            # 1.1 GFF → BED|GFF to BED conversion
            bed_file = self.config.preprocess_dir / f"{sp_id}.bed"

            self.logger.info(f"转换GFF到BED|Converting GFF to BED: {sp_id}")

            args = [
                "--type=mRNA",
                "--key=ID",  # 使用ID而不是Name以确保一致性|Use ID instead of Name for consistency
                gff_file,
                "-o", str(bed_file)
            ]

            if not run_jcvi_command(self.jcvi_python, "formats.gff", "bed",
                                    args, self.logger):
                raise Exception(f"GFF→BED转换失败|GFF to BED conversion failed: {sp_id}")

            self.logger.info(f"BED文件已生成|BED file generated: {bed_file}")

            # 1.2 提取CDS序列|Extract CDS sequences
            # JCVI期望的文件名是 {species}.cds|JCVI expects filename {species}.cds
            # JCVI会从{species}.bed推断{species}.cds|JCVI infers {species}.cds from {species}.bed
            cds_file = self.config.preprocess_dir / f"{sp_id}.cds"

            self.logger.info(f"从GFF提取CDS序列|Extracting CDS sequences from GFF: {sp_id}")

            import os

            # 使用gff load命令从GFF和基因组序列提取CDS|Use gff load to extract CDS from GFF and genome
            args_load = [
                gff_file,
                fa_file,
                "--parents", "mRNA",
                "--children", "CDS",
                "-o", str(cds_file)
            ]

            if not run_jcvi_command(self.jcvi_python, "formats.gff", "load",
                                    args_load, self.logger):
                raise Exception(f"CDS提取失败|CDS extraction failed: {sp_id}")

            self.logger.info(f"CDS文件已生成|CDS file generated: {cds_file}")

        # 1.3 合并所有BED文件|Merge all BED files
        all_bed_file = self.config.preprocess_dir / "all_species.bed"

        self.logger.info("合并所有BED文件|Merging all BED files")

        bed_files = [
            str(self.config.preprocess_dir / f"{sp}.bed")
            for sp in species_list
        ]

        args = bed_files + ["-o", str(all_bed_file)]

        if not run_jcvi_command(self.jcvi_python, "formats.bed", "merge",
                                args, self.logger):
            raise Exception("合并BED文件失败|Failed to merge BED files")

        self.logger.info(f"合并BED文件已生成|Merged BED file generated: {all_bed_file}")

        # 标记步骤完成|Mark step as completed
        self.config.mark_step_done(step)
        self.logger.info(f"步骤{step}完成|Step {step} completed")

    def run_step2_synteny(self):
        """步骤2: 共线性分析|Step 2: Synteny analysis

        - 两两物种序列比对
        - 共线性计算
        - 生成.anchors文件
        """
        step = "2"
        step_name = "共线性分析|Synteny analysis"

        if self.config.is_step_done(step):
            self.logger.info(f"步骤{step}已完成，跳过|Step {step} already completed, skipping")
            return

        self.logger.info(f"步骤{step}: {step_name}|Step {step}: {step_name}")

        species_list = self.config.get_species_list()

        # 两两比对|Pairwise comparison
        species_pairs = list(itertools.combinations(species_list, 2))

        self.logger.info(f"共有|Total {len(species_pairs)} 个物种对需要进行共线性分析|species pairs need synteny analysis")

        # 保存原始工作目录|Save original working directory
        import os
        original_cwd = os.getcwd()
        # 获取preprocess目录的绝对路径|Get absolute path of preprocess directory
        preprocess_abs = os.path.abspath(self.config.preprocess_dir)

        for sp1, sp2 in species_pairs:
            self.logger.info(f"分析物种对|Analyzing species pair: {sp1} vs {sp2}")

            # 切换到preprocess目录，JCVI会在当前目录查找.bed文件|Switch to preprocess dir
            os.chdir(preprocess_abs)

            try:
                # JCVI期望不带.bed后缀的文件名|JCVI expects filenames without .bed extension
                # 从bed文件名中去掉.bed后缀|Strip .bed extension from filenames
                bed1_name = sp1  # Ccu
                bed2_name = sp2  # Ceq

                args = [
                    bed1_name,
                    bed2_name,
                    f"--cscore={self.config.cscore}",
                    "--no_strip_names",
                    "--no_dotplot"  # 跳过dotplot以避免zlib版本问题|Skip dotplot to avoid zlib version issues
                ]

                if not run_jcvi_command(self.jcvi_python, "compara.catalog", "ortholog",
                                        args, self.logger):
                    self.logger.warning(f"共线性分析失败|Synteny analysis failed: {sp1} vs {sp2}")
                    continue

            finally:
                # 恢复原始工作目录|Restore original working directory
                os.chdir(original_cwd)

            # 移动生成的文件到synteny目录|Move generated files to synteny dir
            for ext in ['.last', '.last.filtered', '.anchors', '.lifted.anchors']:
                src = os.path.join(preprocess_abs, f"{sp1}.{sp2}{ext}")
                dst = str(self.config.synteny_dir / f"{sp1}.{sp2}{ext}")

                if os.path.exists(src):
                    os.rename(src, dst)
                    self.logger.info(f"移动文件|Moved file: {src} -> {dst}")
                else:
                    self.logger.warning(f"文件不存在，无法移动|File not found, cannot move: {src}")

            self.logger.info(f"物种对共线性分析完成|Species pair synteny analysis completed: {sp1} vs {sp2}")

        # 标记步骤完成|Mark step as completed
        self.config.mark_step_done(step)
        self.logger.info(f"步骤{step}完成|Step {step} completed")

    def run_step3_blocks(self):
        """步骤3: 区块提取|Step 3: Block extraction

        - 提取微观共线性区块
        - 提取目标基因周围的区块
        """
        step = "3"
        step_name = "区块提取|Block extraction"

        if self.config.is_step_done(step):
            self.logger.info(f"步骤{step}已完成，跳过|Step {step} already completed, skipping")
            return

        self.logger.info(f"步骤{step}: {step_name}|Step {step}: {step_name}")

        species_list = self.config.get_species_list()
        species_pairs = list(itertools.combinations(species_list, 2))

        for sp1, sp2 in species_pairs:
            self.logger.info(f"提取区块|Extracting blocks: {sp1} vs {sp2}")

            lifted_anchors = self.config.synteny_dir / f"{sp1}.{sp2}.lifted.anchors"
            bed1 = str(self.config.preprocess_dir / f"{sp1}.bed")

            if not os.path.exists(lifted_anchors):
                self.logger.warning(f"跳过，文件不存在|Skipping, file not found: {lifted_anchors}")
                continue

            # 使用mcscan提取微观区块
            blocks_file = self.config.blocks_dir / f"{sp1}.{sp2}.blocks"

            args = [
                bed1,
                str(lifted_anchors),
                "--iter=1",
                "-o", str(blocks_file)
            ]

            if not run_jcvi_command(self.jcvi_python, "compara.synteny", "mcscan",
                                    args, self.logger):
                self.logger.warning(f"区块提取失败|Block extraction failed: {sp1} vs {sp2}")
                continue

            self.logger.info(f"区块文件已生成|Block file generated: {blocks_file}")

        # 3.2 提取目标基因所在的区块
        # 对于微观共线性，我们保留每个物种对的独立blocks文件
        # 在绘图阶段会为每个物种对生成独立的共线性图

        self.logger.info("整理区块文件|Organizing block files")

        # 创建一个索引文件，记录所有可用的blocks文件
        blocks_index = self.config.blocks_dir / "blocks_index.txt"
        with open(blocks_index, 'w') as f:
            for sp1, sp2 in species_pairs:
                blocks_file = self.config.blocks_dir / f"{sp1}.{sp2}.blocks"
                if os.path.exists(blocks_file):
                    f.write(f"{sp1}\t{sp2}\t{blocks_file}\n")

        self.logger.info(f"区块索引文件已生成|Block index file generated: {blocks_index}")

        # 标记步骤完成|Mark step as completed
        self.config.mark_step_done(step)
        self.logger.info(f"步骤{step}完成|Step {step} completed")

    def run_step4_plot(self):
        """步骤4: 绘图|Step 4: Plotting

        - 生成包含所有物种的共线性图
        """
        step = "4"
        step_name = "绘图|Plotting"

        if self.config.is_step_done(step):
            self.logger.info(f"步骤{step}已完成，跳过|Step {step} already completed, skipping")
            return

        self.logger.info(f"步骤{step}: {step_name}|Step {step}: {step_name}")

        # 使用pyCirclize可视化
        from .pycirclize_plotter import PycirclizePlotter

        plotter = PycirclizePlotter(self.config, self.logger)

        if not plotter.generate_multi_species_plot():
            self.logger.error("共线性图生成失败|Synteny plot generation failed")
            return

        # 标记步骤完成|Mark step as completed
        self.config.mark_step_done(step)
        self.logger.info(f"步骤{step}完成|Step {step} completed")

    def _generate_pair_layout_file(self, layout_file: Path, sp1: str, sp2: str):
        """生成两物种的layout文件|Generate layout file for two species

        Args:
            layout_file: layout文件路径|Layout file path
            sp1: 物种1|Species 1
            sp2: 物种2|Species 2
        """
        # 两物种的简单layout：上下排列
        # 格式: x, y, rotation, ha, va, color, ratio, label
        with open(layout_file, 'w') as f:
            # 第一个物种在上方 - 添加label显示物种名
            f.write(f"0.5, 0.7, 0, center, top, , 1, {sp1}\n")
            # 第二个物种在下方 - 添加label显示物种名
            f.write(f"0.5, 0.3, 0, center, top, , 1, {sp2}\n")
            # 边连接
            f.write("# edges\n")
            f.write("e, 0, 1\n")

        self.logger.debug(f"Pair layout文件已生成|Pair layout file generated: {layout_file}")
