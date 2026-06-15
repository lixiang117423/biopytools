"""
GWAS2Gene主程序模块|GWAS2Gene Main Module
"""

import sys
import logging
from .config import GWAS2GeneConfig
from .utils import (
    parse_gff, parse_gwas, parse_function_file,
    calc_distance, detect_chromosome_format, normalize_chromosome
)


class GWAS2GeneLogger:
    """GWAS2Gene日志管理器|GWAS2Gene Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """
        设置日志|Setup logging

        Args:
            log_level: 日志级别|Log level
        """
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            import os
            log_path = self.log_file
            log_path_dir = os.path.dirname(log_path)
            if log_path_dir:
                os.makedirs(log_path_dir, exist_ok=True)
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class GWAS2GeneFinder:
    """GWAS候选基因筛选器|GWAS Candidate Gene Finder"""

    def __init__(self, **kwargs):
        """
        初始化GWAS2Gene分析器|Initialize GWAS2Gene analyzer

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = GWAS2GeneConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        import os
        log_file = os.path.splitext(self.config.output_file)[0] + '.log'
        self.logger_manager = GWAS2GeneLogger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

    def run(self):
        """运行GWAS候选基因筛选流程|Run GWAS candidate gene finding pipeline"""
        try:
            self.logger.info("=" * 55)
            self.logger.info("  GWAS Candidate Gene Finder")
            self.logger.info("=" * 55)

            # 1. 解析GFF
            self.logger.info("[1/4] 解析GFF文件|Parsing GFF file...")
            genes, gene_to_transcripts, chrom_bounds = parse_gff(self.config.gff_file)
            gff_chroms = set(chrom_bounds.keys())
            self.logger.info(f"  - 共读取|Total genes: {len(genes)}, 涉及|chromosomes: {len(chrom_bounds)}")

            # 2. 解析GWAS，筛选显著SNP
            self.logger.info(
                f"[2/4] 解析GWAS结果，筛选P < {self.config.threshold:.2e} 的SNP|"
                f"Parsing GWAS, filtering SNPs with P < {self.config.threshold:.2e}..."
            )
            sig_snps = parse_gwas(
                self.config.gwas_file,
                self.config.pval_col,
                self.config.threshold,
                self.logger
            )
            self.logger.info(f"  - 筛选到|Significant SNPs: {len(sig_snps)}")

            if not sig_snps:
                self.logger.info("  没有显著SNP，程序退出|No significant SNPs found, exiting")
                return

            # 检测染色体名称格式|Detect chromosome name format
            gwas_chroms = set(snp['chrom'] for snp in sig_snps)
            gwas_format = detect_chromosome_format(list(gwas_chroms))
            gff_format = detect_chromosome_format(list(gff_chroms))

            self.logger.info(f"  - GWAS染色体格式|GWAS chrom format: {gwas_format}")
            self.logger.info(f"  - GFF染色体格式|GFF chrom format: {gff_format}")

            # 3. 解析功能注释（可选）
            func_map = {}
            if self.config.func_file:
                self.logger.info(f"[3/4] 解析功能注释文件|Parsing function annotation file...")
                func_map = parse_function_file(self.config.func_file, self.logger)
                self.logger.info(f"  - 读取|Annotations loaded: {len(func_map)}")
            else:
                self.logger.info(f"[3/4] 未提供功能注释文件，跳过|No function annotation file, skipping")

            # 4. 区间重叠 & 输出
            self.logger.info(
                f"[4/4] 查找窗口（+/-{self.config.window:,} bp）内的候选基因|"
                f"Finding candidate genes within +/-{self.config.window:,} bp window..."
            )

            total_records = self._write_results(
                sig_snps, genes, gene_to_transcripts,
                chrom_bounds, func_map, gwas_format, gff_format
            )

            self.logger.info(f"  - 共输出|Total records: {total_records}")
            self.logger.info(f"[OK] 结果已保存至|Results saved to: {self.config.output_file}")
            self.logger.info("=" * 55)

        except Exception as e:
            self.logger.error(f"分析流程失败|Pipeline failed: {e}")
            raise

    def _write_results(self, sig_snps, genes, gene_to_transcripts, chrom_bounds, func_map, gwas_format, gff_format):
        """
        写入结果文件|Write results to file

        Args:
            sig_snps: 显著SNP列表|List of significant SNPs
            genes: 基因列表|List of genes
            gene_to_transcripts: 基因到转录本的映射|Gene to transcripts mapping
            chrom_bounds: 染色体边界|Chromosome boundaries
            func_map: 功能注释映射|Function annotation mapping
            gwas_format: GWAS染色体格式|GWAS chromosome format
            gff_format: GFF染色体格式|GFF chromosome format

        Returns:
            int: 总记录数|Total records count
        """
        has_func = bool(func_map)

        out_header = [
            'snp_id', 'snp_chrom', 'snp_pos', 'snp_pval',
            'gene_id', 'gene_name', 'transcript_id',
            'gene_chrom', 'gene_start', 'gene_end', 'gene_strand',
            'distance'
        ]
        if has_func:
            out_header.append('function')

        total_records = 0

        with open(self.config.output_file, 'w') as out:
            out.write('\t'.join(out_header) + '\n')

            for snp in sig_snps:
                snp_chrom = snp['chrom']
                pos = snp['pos']

                # 遍历所有gene，找重叠|Iterate all genes, find overlaps
                for gene in genes:
                    gene_chrom = gene['chrom']

                    # 规范化染色体名称进行比较|Normalize chromosome names for comparison
                    # 统一转换为GFF格式进行比较|Convert both to GFF format for comparison
                    norm_snp_chrom = normalize_chromosome(snp_chrom, gff_format)
                    norm_gene_chrom = normalize_chromosome(gene_chrom, gff_format)

                    if norm_snp_chrom != norm_gene_chrom:
                        continue

                    # 窗口边界，截断到染色体范围|Window bounds, truncated to chromosome range
                    if gene_chrom in chrom_bounds:
                        chrom_min, chrom_max = chrom_bounds[gene_chrom]
                    else:
                        chrom_min, chrom_max = 1, float('inf')

                    win_start = max(pos - self.config.window, chrom_min)
                    win_end = min(pos + self.config.window, chrom_max)

                    g_start = gene['start']
                    g_end = gene['end']

                    # 判断重叠：窗口与基因坐标有交集|Check overlap
                    if g_end < win_start or g_start > win_end:
                        continue

                    gene_id = gene['gene_id']
                    transcripts = gene_to_transcripts.get(gene_id, [])
                    transcript_str = ','.join(transcripts) if transcripts else 'NA'

                    distance = calc_distance(pos, g_start, g_end, gene['strand'])

                    # 功能注释匹配|Function annotation matching
                    func_str = 'NA'
                    if has_func:
                        if gene_id in func_map:
                            func_str = func_map[gene_id]
                        else:
                            for tid in transcripts:
                                if tid in func_map:
                                    func_str = func_map[tid]
                                    break

                    row = [
                        snp['snp_id'],
                        snp_chrom,
                        str(pos),
                        f"{snp['pval']:.2e}",
                        gene_id,
                        gene['gene_name'],
                        transcript_str,
                        gene_chrom,
                        str(g_start),
                        str(g_end),
                        gene['strand'],
                        str(distance),
                    ]
                    if has_func:
                        row.append(func_str)

                    out.write('\t'.join(row) + '\n')
                    total_records += 1

        return total_records


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='GWAS候选基因筛选工具|GWAS Candidate Gene Finder',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--gwas', required=True,
                       help='GWAS结果文件路径|GWAS result file path')
    parser.add_argument('--pval-col', required=True,
                       help='P值所在列名或列号（1-based）|P-value column name or index (1-based)')
    parser.add_argument('--threshold', type=float, default=1e-5,
                       help='P值阈值|P-value threshold (default: 1e-5)')
    parser.add_argument('--window', type=int, default=200000,
                       help='上下游窗口大小|Window size upstream/downstream in bp (default: 200000)')
    parser.add_argument('--gff', required=True,
                       help='GFF3注释文件路径|GFF3 annotation file path')
    parser.add_argument('--output', required=True,
                       help='输出文件路径|Output file path')
    parser.add_argument('--func', default=None,
                       help='功能注释文件路径|Function annotation file path (optional)')

    args = parser.parse_args()

    # 构建配置参数字典|Build configuration kwargs
    kwargs = {
        'gwas_file': args.gwas,
        'pval_col': args.pval_col,
        'threshold': args.threshold,
        'window': args.window,
        'gff_file': args.gff,
        'output_file': args.output,
        'func_file': args.func
    }

    # 创建查找器并运行|Create finder and run
    finder = GWAS2GeneFinder(**kwargs)
    finder.run()


if __name__ == "__main__":
    main()
