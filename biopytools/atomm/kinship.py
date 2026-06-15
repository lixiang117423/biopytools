"""ATOMM遗传关系矩阵计算模块|ATOMM Kinship Matrix Calculation Module"""

import numpy as np


def calculate_kinship(genotype, maf_threshold=0.05):
    """计算经验遗传关系矩阵(GRM)|Calculate empirical Genetic Relatedness Matrix (GRM)

    使用VanRaden方法，与ATOMM MATLAB实现一致
    VanRaden method, consistent with ATOMM MATLAB implementation

    Args:
        genotype: 基因型矩阵 (n_snps x n_individuals), 0/1编码
        maf_threshold: MAF过滤阈值，低于此阈值的SNP被排除|MAF filter threshold

    Returns:
        tuple: (kinship_matrix, maf_array)
            kinship_matrix: GRM (n_individuals x n_individuals)
            maf_array: 每个SNP的MAF
    """
    n_snps, n_ind = genotype.shape

    maf = np.mean(genotype, axis=1)

    mask = (maf > maf_threshold) & (maf < (1 - maf_threshold))
    maf_filtered = maf[mask]
    genotype_filtered = genotype[mask, :]

    n_kept = int(mask.sum())

    matrix = np.zeros((n_ind, n_kept))
    for i in range(n_kept):
        g = genotype_filtered[i, :]
        matrix[:, i] = (g - maf_filtered[i]) / np.sqrt(maf_filtered[i] * (1 - maf_filtered[i]))

    kinship = matrix @ matrix.T / n_kept

    return kinship, maf
