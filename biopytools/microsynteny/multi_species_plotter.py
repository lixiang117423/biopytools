"""
多物种微观共线性绘图工具|Multi-species microsynteny plotter
"""

import os
import subprocess
from pathlib import Path
from typing import List, Set


class MultiSpeciesPlotter:
    """多物种共线性图生成器|Multi-species synteny plot generator"""

    def __init__(self, config, logger, jcvi_python: str):
        """初始化|Initialize

        Args:
            config: MicrosyntenyConfig对象|MicrosyntenyConfig object
            logger: 日志器|Logger
            jcvi_python: JCVI Python路径|JCVI Python path
        """
        self.config = config
        self.logger = logger
        self.jcvi_python = jcvi_python

    def generate_multi_species_plot(self):
        """生成多物种共线性图|Generate multi-species synteny plot"""
        species_list = self.config.get_species_list()

        # 选择第一个物种作为参考物种
        reference_species = species_list[0]
        self.logger.info(f"选择参考物种|Reference species: {reference_species}")

        # 收集所有与参考物种比较的blocks文件
        all_blocks = []
        species_to_plot = [reference_species]

        for sp in species_list[1:]:
            blocks_file = self.config.blocks_dir / f"{reference_species}.{sp}.blocks"
            if os.path.exists(blocks_file):
                all_blocks.append(blocks_file)
                species_to_plot.append(sp)

        # 限制最多3-4个物种以避免图太复杂
        if len(species_to_plot) > 4:
            self.logger.warning(f"物种数量过多({len(species_to_plot)})，只使用前4个|Too many species ({len(species_to_plot)}), using first 4")
            species_to_plot = species_to_plot[:4]
            all_blocks = all_blocks[:3]

        if not all_blocks:
            self.logger.error("没有找到blocks文件|No blocks files found")
            return False

        self.logger.info(f"将绘制|Plotting {len(species_to_plot)}个物种: {', '.join(species_to_plot)}")

        # 合并所有物种的blocks到多列格式
        merged_multi_blocks = self.config.plot_dir / "multi_species.blocks"

        self.logger.info("合并多物种blocks文件|Merging multi-species blocks")

        # 使用JCVI的base join命令合并blocks
        join_cmd = [
            self.jcvi_python,
            "-m", "jcvi.formats.base",
            "join"
        ]

        # 添加所有blocks文件作为参数
        for blocks_file in all_blocks:
            join_cmd.extend(["--noheader", str(blocks_file)])

        # 执行join命令并处理输出
        try:
            result = subprocess.run(
                join_cmd,
                capture_output=True,
                text=True,
                check=True,
                env=os.environ.copy()
            )

            # 处理输出：提取需要的列（不进行列过滤，留到提取目标基因后再过滤）
            # 对于N个物种对（都是reference vs others），join输出会有 2*N-1 列
            # 我们需要：列1（参考物种），列2, 4, 6, 8...（其他物种）
            raw_lines = []
            for line in result.stdout.split('\n'):
                if not line.strip():
                    continue
                parts = line.strip().split('\t')

                # 提取需要的列：第1列 + 偶数列（2, 4, 6, ...）
                selected_parts = []
                for i, part in enumerate(parts):
                    if i == 0:  # 第一列（参考物种）
                        selected_parts.append(part)
                    elif i % 2 == 1:  # 奇数索引（偶数列：2, 4, 6, ...）
                        selected_parts.append(part)

                raw_lines.append(selected_parts)

            # 直接写入合并的blocks文件，不进行列过滤
            if raw_lines:
                with open(merged_multi_blocks, 'w') as outf:
                    for line in raw_lines:
                        outf.write('\t'.join(line) + '\n')
                n_cols = len(raw_lines[0])
                self.logger.info(f"多物种blocks文件已生成，包含{n_cols}列|Multi-species blocks file generated with {n_cols} columns")
            else:
                # 如果没有数据，创建空文件
                with open(merged_multi_blocks, 'w') as outf:
                    pass
                n_cols = 0

        except subprocess.CalledProcessError as e:
            self.logger.error(f"合并blocks失败|Failed to merge blocks: {e.stderr}")
            return False

        # 提取包含目标基因的blocks
        target_genes = set()
        for sp in species_to_plot:
            if sp in self.config.species_genes:
                target_genes.update(self.config.species_genes[sp])

        extracted_blocks_raw = self.config.plot_dir / "multi_species.extracted.raw.blocks"

        # 提取包含目标基因的blocks，但限制数量以避免图太密集
        MAX_BLOCKS = 50

        blocks_extracted = self._extract_blocks_with_targets(
            merged_multi_blocks, extracted_blocks_raw, target_genes, max_blocks=MAX_BLOCKS
        )

        self.logger.info(f"提取包含目标基因的blocks（最多{MAX_BLOCKS}个）|Extracted {blocks_extracted} blocks containing target genes (max {MAX_BLOCKS})")

        # 备选方案：如果没有找到，使用前50个blocks
        if blocks_extracted == 0:
            self.logger.warning(f"未找到包含目标基因的blocks，使用前{MAX_BLOCKS}个|No blocks with target genes found, using first {MAX_BLOCKS}")
            with open(merged_multi_blocks, 'r') as inf, open(extracted_blocks_raw, 'w') as outf:
                for i, line in enumerate(inf):
                    if i >= MAX_BLOCKS:
                        break
                    outf.write(line)

        # 关键修复：在提取目标基因后，过滤空列
        extracted_blocks = self.config.plot_dir / "multi_species.extracted.blocks"

        # 读取提取的blocks，检查哪些列有有效数据
        extracted_lines = []
        with open(extracted_blocks_raw, 'r') as inf:
            for line in inf:
                if line.strip():
                    extracted_lines.append(line.strip().split('\t'))

        if extracted_lines:
            n_cols_extracted = len(extracted_lines[0])

            # 检查每一列是否有有效基因
            valid_cols_extracted = []
            for col_idx in range(n_cols_extracted):
                has_valid_gene = any(
                    line[col_idx] != "."
                    for line in extracted_lines
                )
                if has_valid_gene:
                    valid_cols_extracted.append(col_idx)

            self.logger.info(f"提取的blocks有{n_cols_extracted}列，{len(valid_cols_extracted)}列有有效数据|Extracted blocks has {n_cols_extracted} columns, {len(valid_cols_extracted)} have data")

            # 只保留有效列
            with open(extracted_blocks, 'w') as outf:
                for line in extracted_lines:
                    filtered_line = [line[i] for i in valid_cols_extracted]
                    outf.write('\t'.join(filtered_line) + '\n')

            # 更新species_to_plot列表，移除没有有效数据的物种
            species_to_plot = [species_to_plot[i] for i in valid_cols_extracted]
            self.logger.info(f"最终绘制|Finally plotting {len(species_to_plot)}个物种: {', '.join(species_to_plot)}")

        # 合并所有物种的BED文件
        all_beds = [str(self.config.preprocess_dir / f"{sp}.bed") for sp in species_to_plot]
        merged_bed = self.config.plot_dir / "all_species.merged.bed"

        self.logger.info("合并所有BED文件|Merging all BED files")

        args_merge = all_beds + ["-o", str(merged_bed)]

        from .utils import run_jcvi_command

        if not run_jcvi_command(self.jcvi_python, "formats.bed", "merge",
                                args_merge, self.logger):
            self.logger.error("BED文件合并失败|Failed to merge BED files")
            return False

        # 生成多物种layout文件
        layout_file = self.config.plot_dir / "multi_species.layout"

        self.logger.info(f"生成layout文件|Generating layout file: {layout_file}")
        self._generate_multi_layout_file(layout_file, species_to_plot)

        # 调用JCVI绘图
        self.logger.info("开始绘图|Starting plotting")

        # 生成SVG和PNG格式（PDF格式在此环境下有问题）
        # Generate SVG and PNG formats (PDF has issues in this environment)
        formats_to_try = [('svg', 'SVG'), ('png', 'PNG')]
        success = False
        generated_files = []  # 记录所有生成的文件以便后续重命名

        for format_name, format_label in formats_to_try:
            self.logger.info(f"生成{format_label}格式图|Generating {format_label} format plot")

            # 直接使用JCVI模块命令，使用绝对路径
            plot_cmd = [
                self.jcvi_python, "-m", "jcvi.graphics.synteny",
                str(extracted_blocks.resolve()),
                str(merged_bed.resolve()),
                str(layout_file.resolve()),
                "--genelabelsize=6",
                "--format", format_name
            ]

            try:
                result = subprocess.run(
                    plot_cmd,
                    capture_output=True,
                    text=True,
                    check=False,  # 不检查返回码，有些情况下即使成功也会返回非0
                    env=os.environ.copy(),
                    timeout=300  # 5分钟超时
                )
                # 检查输出文件是否生成
                output_file = extracted_blocks.with_suffix(f'.{format_name}')
                if output_file.exists() and output_file.stat().st_size > 0:
                    generated_files.append((format_name, output_file))
                    self.logger.info(f"{format_label}文件已生成|{format_label} file generated: {output_file.name} ({output_file.stat().st_size / 1024 / 1024:.1f} MB)")
                    success = True
                else:
                    self.logger.warning(f"{format_label}文件未生成或为空|{format_label} file not generated or empty")
                    if result.stdout:
                        self.logger.warning(f"Stdout: {result.stdout[:1000]}")
                    if result.stderr:
                        self.logger.warning(f"Stderr: {result.stderr[:1000]}")
            except subprocess.TimeoutExpired:
                self.logger.warning(f"{format_label}格式生成超时|{format_label} format generation timeout")
                continue

        # 重命名所有生成的文件
        for format_name, src_file in generated_files:
            dst = self.config.plot_dir / f"multi_species_synteny.{format_name}"
            # 如果目标文件已存在，先删除
            if dst.exists():
                os.remove(dst)
            os.rename(src_file, dst)
            self.logger.info(f"图片已保存|Image saved: {dst} ({dst.stat().st_size / 1024 / 1024:.1f} MB)")

        if not success:
            self.logger.error("所有格式都生成失败|All format generation failed")
            return False

        self.logger.info("绘图完成|Plotting completed")
        return True

    def _extract_blocks_with_targets(
        self,
        input_file: Path,
        output_file: Path,
        target_genes: Set[str],
        max_blocks: int = None
    ) -> int:
        """提取包含目标基因的blocks|Extract blocks containing target genes

        Args:
            input_file: 输入blocks文件|Input blocks file
            output_file: 输出blocks文件|Output blocks file
            target_genes: 目标基因集合|Target genes set
            max_blocks: 最大提取数量|Maximum number to extract

        Returns:
            int: 提取的blocks数量|Number of extracted blocks
        """
        blocks_extracted = 0
        with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
            for line in inf:
                genes_in_block = line.strip().split('\t')
                block_has_target = False

                for gene in genes_in_block:
                    if gene in target_genes:
                        block_has_target = True
                        break

                if block_has_target:
                    outf.write(line)
                    blocks_extracted += 1

                    # 如果设置了最大数量限制，达到后停止
                    if max_blocks and blocks_extracted >= max_blocks:
                        break

        return blocks_extracted

    def _generate_multi_layout_file(self, layout_file: Path, species_list: List[str]):
        """生成多物种layout文件|Generate multi-species layout file

        Args:
            layout_file: layout文件路径|Layout file path
            species_list: 物种列表|Species list
        """
        n_species = len(species_list)

        # 计算y坐标位置|Calculate y coordinates
        # 从上到下排列|Arrange from top to bottom
        y_start = 0.9
        y_end = 0.1
        y_step = (y_start - y_end) / (n_species - 1) if n_species > 1 else 0

        with open(layout_file, 'w') as f:
            # 写入物种行 - 使用微观共线性格式
            # 格式: x, y, rotation, ha, va, color, ratio, label
            for i, sp_id in enumerate(species_list):
                y_pos = y_start - (i * y_step)
                f.write(f"0.5, {y_pos:.2f}, 0, center, top, , 1, {sp_id}\n")

            # 写入边连接|Write edges
            # 连接所有相邻物种对
            f.write("# edges\n")
            for i in range(n_species - 1):
                f.write(f"e, {i}, {i+1}\n")

        self.logger.info(f"Layout文件已生成，包含{n_species}个物种|Layout file generated with {n_species} species")
