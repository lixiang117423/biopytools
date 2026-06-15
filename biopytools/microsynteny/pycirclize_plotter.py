"""
pyCirclize共线性可视化|pyCirclize Synteny Visualization
"""

import subprocess
import json
from pathlib import Path
from typing import Dict
import os
from .utils import build_conda_command


class PycirclizePlotter:
    """基于pyCirclize的共线性可视化类|Synteny visualization class using pyCirclize"""

    def __init__(self, config, logger):
        """初始化|Initialize

        Args:
            config: MicrosyntenyConfig对象|MicrosyntenyConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 获取pyCirclize环境的Python路径
        self.pycirclize_python = os.path.join(config.pycirclize_path, "bin", "python")

    def generate_multi_species_plot(self) -> bool:
        """生成多物种共线性图|Generate multi-species synteny plot

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            species_list = self.config.get_species_list()
            self.logger.info(f"开始生成多物种共线性图|Starting multi-species synteny plot generation")
            self.logger.info(f"物种数量|Number of species: {len(species_list)}")
            self.logger.info(f"使用pyCirclize环境|Using pyCirclize environment: {self.config.pycirclize_path}")

            # 1. 读取merged BED文件
            merged_bed = self.config.preprocess_dir / "all_species.bed"
            if not merged_bed.exists():
                self.logger.error(f"BED文件不存在|BED file not found: {merged_bed}")
                return False

            self.logger.info(f"读取基因位置文件|Reading gene positions: {merged_bed.name}")

            # 2. 读取blocks文件（所有两两物种对的blocks）
            blocks_files = list(self.config.blocks_dir.glob("*.blocks"))
            if not blocks_files:
                self.logger.error(f"Blocks文件不存在|No blocks files found in: {self.config.blocks_dir}")
                return False

            self.logger.info(f"找到|Found {len(blocks_files)} 个共线性区块文件|synteny blocks files")

            # 3. 获取目标基因列表
            target_genes_list = []
            for genes in self.config.species_genes.values():
                target_genes_list.extend(genes)
            self.logger.info(f"目标基因数量|Target genes: {len(target_genes_list)}")

            # 4. 准备数据传递给subprocess
            data = {
                'merged_bed': str(merged_bed),
                'blocks_files': [str(f) for f in blocks_files],
                'plot_dir': str(self.config.plot_dir),
                'target_genes': self.config.species_genes,
                'target_genes_list': target_genes_list,
                'extend_genes': self.config.extend_genes
            }

            # 4. 在指定的Python环境中运行pyCirclize脚本
            self.logger.info("使用pyCirclize生成图表|Generating plot with pyCirclize")
            success = self._run_pycirclize_in_env(data)

            if success:
                self.logger.info("共线性图生成完成|Synteny plot generation completed")
            return success

        except Exception as e:
            self.logger.error(f"共线性图生成失败|Synteny plot generation failed: {e}")
            import traceback
            self.logger.debug(f"错误详情|Error details:\n{traceback.format_exc()}")
            return False

    def _run_pycirclize_in_env(self, data: dict) -> bool:
        """在指定的Python环境中运行pyCirclize脚本|Run pyCirclize script in specified Python environment

        Args:
            data: 数据字典|Data dictionary

        Returns:
            bool: 是否成功|Whether successful
        """
        # 创建临时脚本文件
        script_content = '''
import sys
import pandas as pd
from pycirclize import Circos

# 从参数加载数据
import json
data = json.loads(sys.argv[1])

merged_bed_path = data['merged_bed']
blocks_files_paths = data['blocks_files']
plot_dir = data['plot_dir']
target_genes_dict = data['target_genes']
target_genes_list = data['target_genes_list']
extend_genes = data.get('extend_genes', 30)

# 读取BED文件
bed_df = pd.read_csv(merged_bed_path, sep='\\\\t', header=None,
                     names=['chrom', 'start', 'end', 'gene_id', 'score', 'strand'])

print(f"Loaded {len(bed_df)} genes from BED file")

# 获取目标基因集合
target_genes_set = set(target_genes_list)
print(f"Target genes: {len(target_genes_set)}")

# 筛选包含目标基因的染色体，并扩展区域
selected_chroms = set()
for target_gene in target_genes_set:
    # 匹配基因ID（可能带有后缀）
    matching_rows = bed_df[bed_df['gene_id'].str.startswith(target_gene)]
    if len(matching_rows) > 0:
        for _, row in matching_rows.iterrows():
            chrom = row['chrom']
            selected_chroms.add(chrom)
            print(f"Found target gene {target_gene} on {chrom}")

print(f"Selected {len(selected_chroms)} chromosomes containing target genes")

# 筛选BED数据，只保留选中的染色体
filtered_bed = bed_df[bed_df['chrom'].isin(selected_chroms)].copy()
print(f"Filtered to {len(filtered_bed)} genes on selected chromosomes")

# 对每个染色体，提取目标基因周围的基因
extended_regions = {}
for chrom in selected_chroms:
    chrom_genes = filtered_bed[filtered_bed['chrom'] == chrom].reset_index(drop=True)

    # 找出该染色体上的所有目标基因索引
    target_indices = []
    for idx, row in chrom_genes.iterrows():
        gene_id = row['gene_id']
        if any(gene_id.startswith(g) for g in target_genes_set):
            target_indices.append(idx)

    if not target_indices:
        continue

    # 扩展每个目标基因的区域
    extended_indices = set()
    for idx in target_indices:
        start_idx = max(0, idx - extend_genes)
        end_idx = min(len(chrom_genes), idx + extend_genes + 1)
        extended_indices.update(range(start_idx, end_idx))

    # 获取扩展区域的基因
    extended_regions[chrom] = chrom_genes.iloc[list(extended_indices)]

    print(f"Chromosome {chrom}: {len(target_indices)} target genes, extended to {len(extended_indices)} genes")

# 合并所有扩展区域的基因
final_bed = pd.concat(extended_regions.values(), ignore_index=True)
print(f"Final BED contains {len(final_bed)} genes in extended regions")

# 计算扇区大小（使用扩展后的最大位置）
sector_sizes = {}
for chrom, genes in extended_regions.items():
    max_pos = genes['end'].max()
    sector_sizes[chrom] = max_pos

print(f"Created {len(sector_sizes)} sectors with extended regions")

# 提取物种信息并分配颜色
def get_species_from_chrom(chrom_name):
    """从染色体名称提取物种名|Extract species name from chromosome name"""
    # 尝试不同的分隔符
    for sep in ['_', '-']:
        if sep in chrom_name:
            parts = chrom_name.split(sep)
            if len(parts) >= 2:
                # 第一部分通常是物种名
                return parts[0]
    # 如果没有分隔符，使用前2-3个字符
    return chrom_name[:3]

# 为每个物种分配颜色
species_list = set(get_species_from_chrom(chrom) for chrom in extended_regions.keys())
color_map = {
    'Ccu': '#FF6B6B',  # 红色
    'Ceq': '#4ECDC4',  # 青色
    'Cgl': '#45B7D1',  # 蓝色
    'Ptr': '#96CEB4',  # 绿色
    'Ath': '#FFEEAD',  # 黄色
}
# 如果有其他物种，使用默认颜色
default_colors = ['#FF9FF3', '#54A0FF', '#5F27CD', '#00D2D3', '#FF9F43']
for i, species in enumerate(species_list):
    if species not in color_map:
        color_map[species] = default_colors[i % len(default_colors)]

print(f"Species colors: {color_map}")

# 初始化Circos
circos = Circos(sector_sizes, space=5)

# 为每个扇区绘制基因
for sector in circos.sectors:
    if sector.name not in extended_regions:
        continue

    # 获取物种并设置颜色
    species = get_species_from_chrom(sector.name)
    sector_color = color_map.get(species, '#CCCCCC')

    # 添加染色体名称标签
    sector.text(sector.name, size=10, orientation="horizontal", r=110)

    chrom_genes = extended_regions[sector.name]
    track = sector.add_track((90, 100), r_pad_ratio=0.1)
    track.axis(ec="lightgrey", lw=0.5, fc=sector_color, alpha=0.3)

    for _, gene in chrom_genes.iterrows():
        gene_id = gene['gene_id']
        is_target = any(gene_id.startswith(g) for g in target_genes_set)
        color = "salmon" if gene['strand'] == '+' else "skyblue"
        alpha = 1.0 if is_target else 0.5

        try:
            track.genomic_features(
                {'start': gene['start'], 'end': gene['end'],
                 'strand': 1 if gene['strand'] == '+' else -1},
                plotstyle="arrow",
                fc=color,
                alpha=alpha,
                lw=0.5
            )
        except Exception as e:
            pass

        if is_target:
            # 清理基因名称：去掉transient前缀和多余后缀
            gene_name = gene_id.split('-')[0]
            if ':' in gene_name:
                gene_name = gene_name.split(':')[-1]

            try:
                # 使用最简单的参数
                track.text((gene['start'] + gene['end']) / 2, gene_name)
            except Exception as e:
                print(f"Failed to add label for {gene_name}: {e}")
                pass

# 绘制连接线 - 读取所有blocks文件
# 首先收集所有连接，找出有共线性的染色体
gene_pos_index = final_bed.set_index('gene_id')[['chrom', 'start', 'end']].to_dict('index')
chromosomes_with_links = set()

# 第一遍：找出所有有连接的染色体
all_links = []
for blocks_file_path in blocks_files_paths:
    blocks_df = pd.read_csv(blocks_file_path, sep='\\\\t', header=None)

    for _, row in blocks_df.iterrows():
        genes = [g for g in row if pd.notna(g) and g != '.']

        for i in range(len(genes) - 1):
            gene1, gene2 = genes[i], genes[i + 1]

            if gene1 in gene_pos_index and gene2 in gene_pos_index:
                pos1 = gene_pos_index[gene1]
                pos2 = gene_pos_index[gene2]
                all_links.append((pos1, pos2))
                chromosomes_with_links.add(pos1['chrom'])
                chromosomes_with_links.add(pos2['chrom'])

print(f"Found {len(chromosomes_with_links)} chromosomes with synteny links")

# 过滤extended_regions，只保留有连接的染色体
filtered_regions = {chrom: genes for chrom, genes in extended_regions.items()
                   if chrom in chromosomes_with_links}

print(f"Filtered from {len(extended_regions)} to {len(filtered_regions)} chromosomes with links")

# 如果有染色体被过滤掉，需要重新初始化Circos
if len(filtered_regions) < len(extended_regions):
    # 重新计算扇区大小
    sector_sizes = {}
    for chrom, genes in filtered_regions.items():
        max_pos = genes['end'].max()
        sector_sizes[chrom] = max_pos

    print(f"Re-initialized Circos with {len(sector_sizes)} sectors")

    # 重新初始化Circos
    circos = Circos(sector_sizes, space=5)

    # 重新绘制扇区
    for sector in circos.sectors:
        if sector.name not in filtered_regions:
            continue

        # 获取物种并设置颜色
        species = get_species_from_chrom(sector.name)
        sector_color = color_map.get(species, '#CCCCCC')

        # 添加染色体名称标签
        sector.text(sector.name, size=10, orientation="horizontal", r=110)

        chrom_genes = filtered_regions[sector.name]
        track = sector.add_track((90, 100), r_pad_ratio=0.1)
        track.axis(ec="lightgrey", lw=0.5, fc=sector_color, alpha=0.3)

        for _, gene in chrom_genes.iterrows():
            gene_id = gene['gene_id']
            is_target = any(gene_id.startswith(g) for g in target_genes_set)
            color = "salmon" if gene['strand'] == '+' else "skyblue"
            alpha = 1.0 if is_target else 0.5

            try:
                track.genomic_features(
                    {'start': gene['start'], 'end': gene['end'],
                     'strand': 1 if gene['strand'] == '+' else -1},
                    plotstyle="arrow",
                    fc=color,
                    alpha=alpha,
                    lw=0.5
                )
            except Exception as e:
                pass

            if is_target:
                # 清理基因名称：去掉transient前缀和多余后缀
                gene_name = gene_id.split('-')[0]
                if ':' in gene_name:
                    gene_name = gene_name.split(':')[-1]

                try:
                    # 使用最简单的参数
                    track.text((gene['start'] + gene['end']) / 2, gene_name)
                except Exception as e:
                    print(f"Failed to add label for {gene_name}: {e}")
                    pass

# 第二遍：绘制连接线
links_created = 0
for pos1, pos2 in all_links:
    # 确保两个染色体都在过滤后的列表中
    if pos1['chrom'] in filtered_regions and pos2['chrom'] in filtered_regions:
        try:
            circos.link(
                (pos1['chrom'], pos1['start'], pos1['end']),
                (pos2['chrom'], pos2['start'], pos2['end']),
                color="orange", alpha=0.3, linewidth=0.5, ec="black", lw=0.3
            )
            links_created += 1
        except Exception as e:
            pass

print(f"Created {links_created} links total")

# 保存图片
for fmt in ['png', 'svg', 'pdf']:
    try:
        output_file = f"{plot_dir}/multi_species_synteny.{fmt}"
        circos.savefig(output_file, dpi=300)
        print(f"Saved {output_file}")
    except Exception as e:
        print(f"Failed to save {fmt}: {e}")

sys.exit(0)
'''

        script_file = self.config.plot_dir / "plot_script.py"
        with open(script_file, 'w') as f:
            f.write(script_content)

        # 运行脚本（使用conda包装）|Run script with conda wrapping
        # 提取Python命令名用于conda环境检测|Extract Python command name for conda env detection
        python_cmd = os.path.basename(self.pycirclize_python)  # python
        script_args = [str(script_file.absolute()), json.dumps(data)]

        wrapped_cmd = build_conda_command(python_cmd, script_args)

        self.logger.debug(f"运行命令|Running command: {' '.join(wrapped_cmd)}")

        result = subprocess.run(
            wrapped_cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5分钟超时|5 min timeout
        )

        # 清理临时脚本
        try:
            script_file.unlink()
        except:
            pass

        if result.returncode != 0:
            self.logger.error(f"pyCirclize执行失败|pyCirclize execution failed")
            if result.stderr:
                self.logger.error(f"错误信息|Error:\n{result.stderr}")
            return False

        if result.stdout:
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    self.logger.info(f"pyCirclize: {line}")

        # 检查输出文件
        for fmt in ['png', 'svg', 'pdf']:
            output_file = self.config.plot_dir / f"multi_species_synteny.{fmt}"
            if output_file.exists() and output_file.stat().st_size > 0:
                size_mb = output_file.stat().st_size / 1024 / 1024
                self.logger.info(f"图片已保存|Image saved: {output_file.name} ({size_mb:.2f} MB)")

        return True
