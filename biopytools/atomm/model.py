"""ATOMM统计模型模块|ATOMM Statistical Model Module

包含零模型(Null Model)方差分量估计和三种Score检验:
- 宿主边际检验 (Host Marginal Test)
- 病原边际检验 (Pathogen Marginal Test)
- 宿主-病原交互检验 (Interaction Test)

数学模型 (Wang et al. 2018):
    Y_k = X_k * beta + G_i^h * gamma_h + G_j^p * gamma_p + G_i^h * G_j^p * gamma_hp
          + eta_i^h + eta_j^p + eta_ij^hp + epsilon_k

    Sigma = xi_h * K_h_total + xi_p * K_p_total + xi_hp * K_cross + (1 - sum_xi) * I
    K_cross = K_h_total .* K_p_total  (逐元素乘积|element-wise product)
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2


class ATOMMModel:
    """ATOMM双物种混合效应模型|ATOMM Two-Organism Mixed-Effect Model"""

    def __init__(self, host_ids, pathogen_ids, covariates, phenotypes,
                 kinship_h, kinship_p, logger=None):
        """初始化ATOMM模型|Initialize ATOMM model

        Args:
            host_ids: 宿主ID数组 (n_obs,)，从1开始|Host ID array
            pathogen_ids: 病原ID数组 (n_obs,)，从1开始|Pathogen ID array
            covariates: 协变量矩阵 (n_obs, n_cov)，包含截距|Covariate matrix (with intercept)
            phenotypes: 表型值 (n_obs,)|Phenotype values
            kinship_h: 宿主GRM (n_h x n_h)|Host GRM
            kinship_p: 病原GRM (n_p x n_p)|Pathogen GRM
            logger: 日志器|Logger instance
        """
        self.logger = logger
        self.n = len(phenotypes)
        self.n_h = kinship_h.shape[0]
        self.n_p = kinship_p.shape[0]

        self.X = covariates
        self.Y = phenotypes
        self.n_cov = covariates.shape[1]

        self.kinship_h = kinship_h
        self.kinship_p = kinship_p

        self.obs_host_ids = host_ids
        self.obs_pathogen_ids = pathogen_ids

        host_idx = host_ids - 1
        pathogen_idx = pathogen_ids - 1

        self.K_h_total = kinship_h[np.ix_(host_idx, host_idx)]
        self.K_p_total = kinship_p[np.ix_(pathogen_idx, pathogen_idx)]
        self.K_cross = self.K_h_total * self.K_p_total

        self.herit = None
        self.Sigma = None

    def estimate_heritability(self, tol=1e-6, maxiter=10000):
        """估计零模型方差分量|Estimate variance components under null model

        使用Nelder-Mead优化最小化profiled负对数似然
        Minimize profiled negative log-likelihood using Nelder-Mead

        Args:
            tol: 优化容忍度|Optimizer tolerance
            maxiter: 最大迭代次数|Max iterations

        Returns:
            tuple: (herit, Sigma)
                herit: [xi_h, xi_p, xi_hp, xi_e]
                Sigma: n x n 协方差矩阵|n x n covariance matrix
        """
        def objective(x):
            if np.min(x) < 1e-10:
                return 1e10

            Sigma = x[0] * self.K_h_total + x[1] * self.K_p_total + x[2] * self.K_cross + np.eye(self.n)
            scale = np.sum(x) + 1
            Sigma = Sigma / scale

            try:
                L = np.linalg.cholesky(Sigma)
                log_det = 2.0 * np.sum(np.log(np.diag(L)))
                Sigma_inv_X = np.linalg.solve(Sigma, self.X)
                XtSX = self.X.T @ Sigma_inv_X
                XtSY = self.X.T @ np.linalg.solve(Sigma, self.Y)
                beta = np.linalg.solve(XtSX, XtSY)
                mu = self.X @ beta
                resid = self.Y - mu
                sigma_t = resid @ np.linalg.solve(Sigma, resid) / self.n
                return log_det + self.n * np.log(sigma_t)
            except np.linalg.LinAlgError:
                return 1e10

        x0 = np.array([0.25, 0.25, 0.25])

        result = minimize(objective, x0, method='Nelder-Mead',
                          options={'maxiter': maxiter, 'xatol': tol, 'fatol': tol, 'adaptive': True})

        x_opt = result.x
        herit = x_opt / (np.sum(x_opt) + 1)
        herit = np.append(herit, 1 - np.sum(herit))

        Sigma = (herit[0] * self.K_h_total + herit[1] * self.K_p_total +
                 herit[2] * self.K_cross + herit[3] * np.eye(self.n))

        self.herit = herit
        self.Sigma = Sigma

        if self.logger:
            self.logger.info(
                f"遗传力估计完成|Heritability estimated: "
                f"host={herit[0]:.4f}, pathogen={herit[1]:.4f}, "
                f"interaction={herit[2]:.4f}, noise={herit[3]:.4f}"
            )

        return herit, Sigma

    def _build_Z_matrix(self, obs_ids, n_individuals):
        """从观察值ID构建指示矩阵|Build incidence matrix from observation IDs

        Z[k, i] = 1 if observation k belongs to individual i, else 0

        Args:
            obs_ids: 每个观察对应的个体ID (n_obs,), 从1开始
            n_individuals: 个体总数

        Returns:
            Z: n_obs x n_individuals 指示矩阵
        """
        idx = obs_ids - 1
        Z = np.zeros((len(obs_ids), n_individuals))
        for i in range(n_individuals):
            Z[idx == i, i] = 1.0
        return Z

    def test_marginal_host(self, genotype_h, index, chrom_ids_in=None, snp_ids_in=None):
        """宿主SNP边际关联检验|Marginal association test for host SNPs

        对应ATOMM_Marginal_host.m的Score检验
        Score test corresponding to ATOMM_Marginal_host.m

        H0: gamma_h = 0, gamma_p = gamma_hp = 0

        Args:
            genotype_h: 宿主基因型矩阵 (n_snps x n_individuals)
            index: 要测试的SNP索引列表 (0-based)|SNP indices to test
            chrom_ids_in: 染色体ID列表|Chromosome ID list
            snp_ids_in: SNP ID列表|SNP ID list

        Returns:
            tuple: (chrom_ids, snp_ids, statistics, p_values)
        """
        Sigma_inv = np.linalg.inv(self.Sigma)

        Z_h = self._build_Z_matrix(self.obs_host_ids, self.n_h)

        beta = np.linalg.solve(self.X.T @ Sigma_inv @ self.X,
                               self.X.T @ Sigma_inv @ self.Y)
        mu = self.X @ beta

        Sigma_inv_Z = Sigma_inv @ Z_h
        deno = (self.Y - mu) @ Sigma_inv_Z @ self.kinship_h @ Sigma_inv_Z.T @ (self.Y - mu)
        temY = (self.Y - mu) @ Sigma_inv_Z
        K_inv = np.linalg.pinv(self.kinship_h, rcond=max(self.kinship_h.shape) * np.finfo(self.kinship_h.dtype).eps * np.linalg.norm(self.kinship_h, ord=2))

        chrom_ids = []
        snp_ids = []
        statistics = np.zeros(len(index))
        p_values = np.zeros(len(index))

        for i, idx in enumerate(index):
            g = genotype_h[idx, :]
            chrom_ids.append(chrom_ids_in[idx] if chrom_ids_in else '1')
            snp_ids.append(snp_ids_in[idx] if snp_ids_in else idx + 1)

            sigma_g = g @ K_inv @ g / (self.n_h - 1)
            T = (temY @ g) ** 2 / (deno * sigma_g)
            statistics[i] = T
            p_values[i] = 1.0 - chi2.cdf(T, 1)

        if self.logger:
            self.logger.info(
                f"宿主边际检验完成|Host marginal test completed: "
                f"{len(index)} SNPs tested"
            )

        return chrom_ids, snp_ids, statistics, p_values

    def test_marginal_pathogen(self, genotype_p, index, chrom_ids_in=None, snp_ids_in=None):
        """病原SNP边际关联检验|Marginal association test for pathogen SNPs

        对应ATOMM_Marginal_pathogen.m的Score检验
        Score test corresponding to ATOMM_Marginal_pathogen.m

        H0: gamma_p = 0, gamma_h = gamma_hp = 0

        Args:
            genotype_p: 病原基因型矩阵 (n_snps x n_individuals)
            index: 要测试的SNP索引列表 (0-based)|SNP indices to test
            chrom_ids_in: 染色体ID列表|Chromosome ID list
            snp_ids_in: SNP ID列表|SNP ID list

        Returns:
            tuple: (chrom_ids, snp_ids, statistics, p_values)
        """
        Sigma_inv = np.linalg.inv(self.Sigma)

        Z_p = self._build_Z_matrix(self.obs_pathogen_ids, self.n_p)

        beta = np.linalg.solve(self.X.T @ Sigma_inv @ self.X,
                               self.X.T @ Sigma_inv @ self.Y)
        mu = self.X @ beta

        Sigma_inv_Z = Sigma_inv @ Z_p
        deno = (self.Y - mu) @ Sigma_inv_Z @ self.kinship_p @ Sigma_inv_Z.T @ (self.Y - mu)
        temY = (self.Y - mu) @ Sigma_inv_Z
        K_inv = np.linalg.pinv(self.kinship_p, rcond=max(self.kinship_p.shape) * np.finfo(self.kinship_p.dtype).eps * np.linalg.norm(self.kinship_p, ord=2))

        chrom_ids = []
        snp_ids = []
        statistics = np.zeros(len(index))
        p_values = np.zeros(len(index))

        for i, idx in enumerate(index):
            g = genotype_p[idx, :]
            chrom_ids.append(chrom_ids_in[idx] if chrom_ids_in else '1')
            snp_ids.append(snp_ids_in[idx] if snp_ids_in else idx + 1)

            sigma_g = g @ K_inv @ g / (self.n_p - 1)
            T = (temY @ g) ** 2 / (deno * sigma_g)
            statistics[i] = T
            p_values[i] = 1.0 - chi2.cdf(T, 1)

        if self.logger:
            self.logger.info(
                f"病原边际检验完成|Pathogen marginal test completed: "
                f"{len(index)} SNPs tested"
            )

        return chrom_ids, snp_ids, statistics, p_values

    def test_interaction(self, genotype_h, genotype_p, index_h, index_p,
                         h_chrom_in=None, h_snp_in=None, p_chrom_in=None, p_snp_in=None):
        """宿主-病原交互检验|Host-pathogen interaction test

        对应ATOMM_Interaction.m的Score检验
        Score test corresponding to ATOMM_Interaction.m

        H0: gamma_hp = 0 (但包含 gamma_h 和 gamma_p 主效应)

        Args:
            genotype_h: 宿主基因型矩阵 (n_snps x n_individuals)
            genotype_p: 病原基因型矩阵 (n_snps x n_individuals)
            index_h: 宿主SNP索引列表 (0-based)|Host SNP indices
            index_p: 病原SNP索引列表 (0-based)|Pathogen SNP indices
            h_chrom_in: 宿主染色体ID列表|Host chromosome ID list
            h_snp_in: 宿主SNP ID列表|Host SNP ID list
            p_chrom_in: 病原染色体ID列表|Pathogen chromosome ID list
            p_snp_in: 病原SNP ID列表|Pathogen SNP ID list

        Returns:
            tuple: (host_chrom, host_snp, pathogen_chrom, pathogen_snp, statistics, p_values)
        """
        Sigma_inv = np.linalg.inv(self.Sigma)
        k_norm = np.linalg.norm(self.K_cross, ord=2)
        k_rcond = max(self.K_cross.shape) * np.finfo(self.K_cross.dtype).eps * k_norm
        K_cross_inv = np.linalg.pinv(self.K_cross, rcond=k_rcond)

        n_pairs = len(index_h) * len(index_p)
        host_chrom = []
        host_snp = []
        pathogen_chrom = []
        pathogen_snp = []
        statistics = np.zeros(n_pairs)
        p_values = np.zeros(n_pairs)

        host_idx = self.obs_host_ids - 1
        pathogen_idx = self.obs_pathogen_ids - 1

        ncount = 0
        for i_idx in index_h:
            gh = genotype_h[i_idx, :]
            gh_obs = gh[host_idx]

            for j_idx in index_p:
                gp = genotype_p[j_idx, :]
                gp_obs = gp[pathogen_idx]

                host_chrom.append(h_chrom_in[i_idx] if h_chrom_in else '1')
                host_snp.append(h_snp_in[i_idx] if h_snp_in else i_idx + 1)
                pathogen_chrom.append(p_chrom_in[j_idx] if p_chrom_in else '1')
                pathogen_snp.append(p_snp_in[j_idx] if p_snp_in else j_idx + 1)

                X_new = np.column_stack([self.X, gh_obs, gp_obs])
                g = gh_obs * gp_obs

                beta = np.linalg.solve(X_new.T @ Sigma_inv @ X_new,
                                       X_new.T @ Sigma_inv @ self.Y)
                mu = X_new @ beta
                resid = self.Y - mu

                deno = resid @ self.K_cross @ Sigma_inv @ self.K_cross @ resid
                sigma_g = g @ K_cross_inv @ g / (self.n - 1)
                T = (resid @ Sigma_inv @ g) ** 2 / (deno * sigma_g)
                statistics[ncount] = T
                p_values[ncount] = 1.0 - chi2.cdf(T, 1)

                ncount += 1

        if self.logger:
            self.logger.info(
                f"交互检验完成|Interaction test completed: "
                f"{n_pairs} pairs tested"
            )

        return host_chrom, host_snp, pathogen_chrom, pathogen_snp, statistics, p_values
