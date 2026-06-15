"""输出模块|Output module"""

import pandas as pd
from pathlib import Path
from typing import List
from .detector import InsertionSite


def write_results(results: List[InsertionSite], output_dir: str) -> Path:
    """写入检测结果到表格|Write detection results to table

    Args:
        results: 插入位点列表|List of insertion sites
        output_dir: 输出目录|Output directory

    Returns:
        Path: 输出文件路径|Output file path
    """
    output_path = Path(output_dir) / "04_results"
    output_path.mkdir(parents=True, exist_ok=True)

    # 转换为DataFrame|Convert to DataFrame
    data = []
    for site in results:
        data.append({
            'sample_id': site.sample_id,
            'chromosome': site.chromosome,
            'start': site.start,
            'end': site.end,
            'orientation': site.orientation,
            'insert_position': site.tdna_position,
            'junction': site.junction,
            'support_reads': site.support_reads,
            'score': site.score,
            'group_counts': site.group_counts
        })

    df = pd.DataFrame(data)

    # 按样品和得分排序|Sort by sample and score
    df = df.sort_values(['sample_id', 'score'], ascending=[True, False])

    # 输出TSV文件|Output TSV file
    output_file = output_path / "insertion_sites.tsv"
    df.to_csv(output_file, sep='\t', index=False)

    return output_file


def write_summary(results: List[InsertionSite], output_dir: str):
    """写入汇总统计|Write summary statistics

    Args:
        results: 插入位点列表|List of insertion sites
        output_dir: 输出目录|Output directory
    """
    output_path = Path(output_dir) / "04_results"

    if not results:
        return

    # 统计每个样品的插入位点数|Count insertion sites per sample
    sample_counts = {}
    for site in results:
        sample_counts[site.sample_id] = sample_counts.get(site.sample_id, 0) + 1

    # 统计染色体分布|Count chromosome distribution
    chrom_counts = {}
    for site in results:
        chrom_counts[site.chromosome] = chrom_counts.get(site.chromosome, 0) + 1

    # 写入汇总文件|Write summary file
    summary_file = output_path / "summary.txt"
    with open(summary_file, 'w') as f:
        f.write("插入位点检测汇总|Insertion Site Detection Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"总插入位点数|Total insertion sites: {len(results)}\n\n")

        f.write("按样品统计|By sample:\n")
        for sample, count in sorted(sample_counts.items()):
            f.write(f"  {sample}: {count}\n")
        f.write("\n")

        f.write("按染色体统计|By chromosome:\n")
        for chrom, count in sorted(chrom_counts.items()):
            f.write(f"  {chrom}: {count}\n")
        f.write("\n")

        # 统计得分分布|Score distribution
        scores = [site.score for site in results]
        f.write(f"得分范围|Score range: {min(scores)} - {max(scores)}\n")
        f.write(f"平均得分|Average score: {sum(scores)/len(scores):.1f}\n")
        f.write(f"中位数得分|Median score: {sorted(scores)[len(scores)//2]}\n")
