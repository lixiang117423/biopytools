"""ATOMM主入口模块|ATOMM Main Entry Module"""

import argparse
import sys
from pathlib import Path
import numpy as np

from .config import ATOMMConfig
from .utils import (
    ATOMMLogger, read_genotype, read_phenotype,
    write_estimate_tsv, write_marginal_host_tsv, write_marginal_pathogen_tsv,
    write_interaction_tsv, write_results_excel, format_number
)
from .kinship import calculate_kinship
from .model import ATOMMModel
from .convert import vcf_to_genotype, phenotype_matrix_to_atomm, detect_vcf_encoding


class ATOMMAnalyzer:
    """ATOMM分析器|ATOMM Analyzer

    Orchestrates the full ATOMM pipeline:
    1. Read input data
    2. Calculate kinship matrices
    3. Estimate heritability (null model)
    4. Run marginal tests (host + pathogen)
    5. Run interaction tests
    6. Write results
    """

    def __init__(self, config: ATOMMConfig, logger=None):
        self.config = config
        self.logger = logger or ATOMMLogger().get_logger()

    def _convert_inputs(self):
        """转换VCF/表型矩阵为ATOMM格式|Convert VCF/phenotype matrix to ATOMM format"""
        logger = self.logger
        output_dir = self.config.output_path

        if self.config.host_vcf:
            encoding = self.config.encoding
            if encoding == 'auto':
                encoding = detect_vcf_encoding(self.config.host_vcf)

            host_out = str(output_dir / "host_genotype.txt")
            logger.info(f"转换宿主VCF→ATOMM格式|Converting host VCF to ATOMM format: {self.config.host_vcf}")
            vcf_to_genotype(
                self.config.host_vcf, host_out,
                maf_threshold=self.config.convert_maf_threshold,
                encoding=encoding
            )
            self.config.host_genotype = host_out

        if self.config.pathogen_vcf:
            encoding = self.config.encoding
            if encoding == 'auto':
                encoding = detect_vcf_encoding(self.config.pathogen_vcf)

            pathogen_out = str(output_dir / "pathogen_genotype.txt")
            logger.info(f"转换病原VCF→ATOMM格式|Converting pathogen VCF to ATOMM format: {self.config.pathogen_vcf}")
            vcf_to_genotype(
                self.config.pathogen_vcf, pathogen_out,
                maf_threshold=self.config.convert_maf_threshold,
                encoding=encoding
            )
            self.config.pathogen_genotype = pathogen_out

        if self.config.phenotype_matrix:
            pheno_out = str(output_dir / "phenotype.txt")
            logger.info(f"转换表型矩阵→ATOMM格式|Converting phenotype matrix to ATOMM format: {self.config.phenotype_matrix}")
            phenotype_matrix_to_atomm(
                self.config.phenotype_matrix, pheno_out,
                missing_value=self.config.missing_value
            )
            self.config.phenotype = pheno_out

    def run(self):
        """运行完整ATOMM分析流程|Run full ATOMM analysis pipeline"""
        logger = self.logger
        logger.info("开始ATOMM分析|Starting ATOMM analysis")

        if self.config.has_convert_inputs:
            self._convert_inputs()

        # Step 1: 读取输入数据
        logger.info("读取宿主基因型数据|Reading host genotype data")
        host_chrom, host_snp, genotype_h = read_genotype(self.config.host_genotype)
        n_h_snps, n_h = genotype_h.shape
        logger.info(f"宿主数据|Host data: {n_h_snps} SNPs, {n_h} individuals")

        logger.info("读取病原基因型数据|Reading pathogen genotype data")
        pathogen_chrom, pathogen_snp, genotype_p = read_genotype(self.config.pathogen_genotype)
        n_p_snps, n_p = genotype_p.shape
        logger.info(f"病原数据|Pathogen data: {n_p_snps} SNPs, {n_p} individuals")

        logger.info("读取表型数据|Reading phenotype data")
        host_ids, pathogen_ids, covariates, phenotypes = read_phenotype(self.config.phenotype)
        n_obs = len(phenotypes)
        n_cov = covariates.shape[1]
        logger.info(f"表型数据|Phenotype data: {n_obs} observations, {n_cov} covariates")

        # Step 2: 计算GRM
        logger.info("计算宿主遗传关系矩阵|Computing host kinship matrix")
        kinship_h, maf_h = calculate_kinship(genotype_h, self.config.maf_threshold)
        logger.info(f"宿主GRM|Host GRM: {kinship_h.shape}, "
                     f"MAF过滤后保留|kept after MAF filter: {int(np.sum((maf_h > self.config.maf_threshold) & (maf_h < 1 - self.config.maf_threshold)))} SNPs")

        logger.info("计算病原遗传关系矩阵|Computing pathogen kinship matrix")
        kinship_p, maf_p = calculate_kinship(genotype_p, self.config.maf_threshold)
        logger.info(f"病原GRM|Pathogen GRM: {kinship_p.shape}, "
                     f"MAF过滤后保留|kept after MAF filter: {int(np.sum((maf_p > self.config.maf_threshold) & (maf_p < 1 - self.config.maf_threshold)))} SNPs")

        # Step 3: 估计遗传力
        logger.info("估计零模型方差分量|Estimating null model variance components")
        model = ATOMMModel(host_ids, pathogen_ids, covariates, phenotypes,
                           kinship_h, kinship_p, logger)
        herit, Sigma = model.estimate_heritability(
            tol=self.config.optimize_tol,
            maxiter=self.config.optimize_maxiter
        )

        output_dir = self.config.output_path

        # Step 4: 边际检验
        # 宿主
        logger.info("运行宿主边际检验|Running host marginal tests")
        if self.config.host_snp_range:
            start, end = self.config.host_snp_range
            h_index = list(range(start, end + 1))
        else:
            h_index = list(range(n_h_snps))

        h_chrom, h_snp_ids, h_stat, h_pval = model.test_marginal_host(
            genotype_h, h_index, chrom_ids_in=host_chrom, snp_ids_in=host_snp)
        write_estimate_tsv(str(output_dir / "estimate.tsv"), herit)
        write_marginal_host_tsv(str(output_dir / "marginal_host.tsv"), h_chrom, h_snp_ids, h_stat, h_pval)
        logger.info(f"宿主边际检验结果已保存|Host marginal results saved: {len(h_index)} SNPs")

        # 病原
        logger.info("运行病原边际检验|Running pathogen marginal tests")
        if self.config.pathogen_snp_range:
            start, end = self.config.pathogen_snp_range
            p_index = list(range(start, end + 1))
        else:
            p_index = list(range(n_p_snps))

        p_chrom, p_snp_ids, p_stat, p_pval = model.test_marginal_pathogen(
            genotype_p, p_index, chrom_ids_in=pathogen_chrom, snp_ids_in=pathogen_snp)
        write_marginal_pathogen_tsv(str(output_dir / "marginal_pathogen.tsv"), p_chrom, p_snp_ids, p_stat, p_pval)
        logger.info(f"病原边际检验结果已保存|Pathogen marginal results saved: {len(p_index)} SNPs")

        # Step 5: 交互检验
        interaction_data = None
        if self.config.interaction_host_range and self.config.interaction_pathogen_range:
            ih_start, ih_end = self.config.interaction_host_range
            ip_start, ip_end = self.config.interaction_pathogen_range
            ih_index = list(range(ih_start, ih_end + 1))
            ip_index = list(range(ip_start, ip_end + 1))
            n_pairs = len(ih_index) * len(ip_index)
            logger.info(f"运行交互检验|Running interaction tests: {n_pairs} pairs")
            hc, hs, pc, ps, i_stat, i_pval = model.test_interaction(
                genotype_h, genotype_p, ih_index, ip_index,
                h_chrom_in=host_chrom, h_snp_in=host_snp,
                p_chrom_in=pathogen_chrom, p_snp_in=pathogen_snp
            )
            interaction_data = (hc, hs, pc, ps, i_stat, i_pval)
            write_interaction_tsv(str(output_dir / "interaction.tsv"), hc, hs, pc, ps, i_stat, i_pval)
            logger.info(f"交互检验结果已保存|Interaction results saved")
        else:
            logger.info("未指定交互检验SNP范围，跳过|No interaction SNP ranges specified, skipping")

        # Step 6: 输出Excel汇总
        excel_path = str(output_dir / "atomm_results.xlsx")
        write_results_excel(
            excel_path,
            herit=herit,
            host_data=(h_chrom, h_snp_ids, h_stat, h_pval),
            pathogen_data=(p_chrom, p_snp_ids, p_stat, p_pval),
            interaction_data=interaction_data if self.config.interaction_host_range else None
        )
        logger.info(f"Excel汇总结果已保存|Excel summary saved: {excel_path}")

        logger.info("ATOMM分析完成|ATOMM analysis completed")
        logger.info(f"结果保存在|Results saved in: {output_dir}")

        return {
            'heritability': herit,
            'host_marginal': (h_chrom, h_snp_ids, h_stat, h_pval),
            'pathogen_marginal': (p_chrom, p_snp_ids, p_stat, p_pval),
        }


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ATOMM: 双物种混合效应模型关联分析|Two-Organism Mixed-Effect Model Association Analysis"
    )

    input_group = parser.add_argument_group('ATOMM格式输入|ATOMM format input')
    input_group.add_argument('-gh', '--host-genotype',
                             help='宿主基因型文件(ATOMM格式)|Host genotype file (ATOMM format)')
    input_group.add_argument('-gp', '--pathogen-genotype',
                             help='病原基因型文件(ATOMM格式)|Pathogen genotype file (ATOMM format)')
    input_group.add_argument('-p', '--phenotype',
                             help='表型文件(ATOMM格式)|Phenotype file (ATOMM format)')

    convert_group = parser.add_argument_group('VCF/矩阵转换输入|VCF/matrix conversion input')
    convert_group.add_argument('--host-vcf',
                               help='宿主VCF文件(自动转换为ATOMM格式)|Host VCF file (auto-convert to ATOMM format)')
    convert_group.add_argument('--pathogen-vcf',
                               help='病原VCF文件(自动转换为ATOMM格式)|Pathogen VCF file (auto-convert to ATOMM format)')
    convert_group.add_argument('--phenotype-matrix',
                               help='交叉感染表型矩阵(行=宿主,列=病原)|Cross-infection phenotype matrix')
    convert_group.add_argument('--encoding', default='auto', choices=['auto', 'haploid', 'dosage'],
                               help='基因型编码方式|Genotype encoding mode (default: auto)')
    convert_group.add_argument('--convert-maf', type=float, default=0.05,
                               help='VCF转换时MAF阈值|MAF threshold for VCF conversion (default: 0.05)')
    convert_group.add_argument('--missing-value', default='NA',
                               help='表型缺失值标记|Missing value marker (default: NA)')

    parser.add_argument('-o', '--output-dir', default='./output',
                        help='输出目录|Output directory')
    parser.add_argument('--maf', type=float, default=0.05,
                        help='MAF过滤阈值|MAF filter threshold (default: 0.05)')

    parser.add_argument('--host-snp-range', type=int, nargs=2, metavar=('START', 'END'),
                        help='宿主边际检验SNP范围(0-based)|Host marginal test SNP range')
    parser.add_argument('--pathogen-snp-range', type=int, nargs=2, metavar=('START', 'END'),
                        help='病原边际检验SNP范围(0-based)|Pathogen marginal test SNP range')
    parser.add_argument('--interaction-host-range', type=int, nargs=2, metavar=('START', 'END'),
                        help='交互检验宿主SNP范围|Interaction test host SNP range')
    parser.add_argument('--interaction-pathogen-range', type=int, nargs=2, metavar=('START', 'END'),
                        help='交互检验病原SNP范围|Interaction test pathogen SNP range')

    parser.add_argument('--tol', type=float, default=1e-6,
                        help='优化容忍度|Optimizer tolerance (default: 1e-6)')
    parser.add_argument('--maxiter', type=int, default=10000,
                        help='优化器最大迭代次数|Max optimizer iterations (default: 10000)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    has_convert = any([args.host_vcf, args.pathogen_vcf, args.phenotype_matrix])

    if not has_convert and not all([args.host_genotype, args.pathogen_genotype, args.phenotype]):
        raise ValueError(
            "必须提供ATOMM格式输入(-gh/-gp/-p)或VCF/矩阵转换输入(--host-vcf/--pathogen-vcf/--phenotype-matrix)|"
            "Must provide either ATOMM format inputs (-gh/-gp/-p) or VCF/matrix conversion inputs"
        )

    try:
        config = ATOMMConfig(
            host_genotype=args.host_genotype or "",
            pathogen_genotype=args.pathogen_genotype or "",
            phenotype=args.phenotype or "",
            output_dir=args.output_dir,
            maf_threshold=args.maf,
            host_snp_range=tuple(args.host_snp_range) if args.host_snp_range else None,
            pathogen_snp_range=tuple(args.pathogen_snp_range) if args.pathogen_snp_range else None,
            interaction_host_range=tuple(args.interaction_host_range) if args.interaction_host_range else None,
            interaction_pathogen_range=tuple(args.interaction_pathogen_range) if args.interaction_pathogen_range else None,
            optimize_tol=args.tol,
            optimize_maxiter=args.maxiter,
            host_vcf=args.host_vcf or "",
            pathogen_vcf=args.pathogen_vcf or "",
            phenotype_matrix=args.phenotype_matrix or "",
            encoding=args.encoding,
            convert_maf_threshold=args.convert_maf,
            missing_value=args.missing_value,
        )
        config.validate()

        logger_manager = ATOMMLogger()
        logger = logger_manager.get_logger()

        analyzer = ATOMMAnalyzer(config, logger)
        analyzer.run()

    except Exception as e:
        print(f"错误|Error: {e}")
        raise


if __name__ == '__main__':
    main()
