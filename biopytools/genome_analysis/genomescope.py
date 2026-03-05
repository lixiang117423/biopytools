"""
GenomeScope核心算法模块|GenomeScope Core Algorithm Module

4峰基因组模型拟合，用于从k-mer直方图估算基因组特征
4-peak genome model fitting for estimating genome characteristics from k-mer histograms
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import nbinom
import warnings
from typing import Optional, Tuple
import logging

warnings.filterwarnings('ignore')


# ============================================================================
# CONSTANTS
# ============================================================================

NUM_ROUNDS = 4                 # Number of optimization rounds
START_SHIFT = 5                # Coverage steps to trim between rounds
TYPICAL_ERROR = 15             # Typical cutoff for sequencing error
MAX_ITERATIONS = 20            # Max iterations for curve fitting
SCORE_CLOSE = 0.20             # Score percentage threshold (20%)
SCORE_HET_FOLD_DIFFERENCE = 10 # Heterozygosity fold difference threshold

# Colors for plots (matching R version)
COLOR_BGCOLOR = '#D3D3D3'      # light grey
COLOR_HIST = '#56B4E9'         # blue
COLOR_4PEAK = 'black'
COLOR_ERRORS = '#D55E00'       # orange
COLOR_KMERPEAK = 'black'
COLOR_RESIDUAL = 'purple'
COLOR_COVTHRES = 'red'

PLOT_SIZE = 1200               # PNG size in pixels
RESOLUTION = 100               # DPI


# ============================================================================
# 4-PEAK MODEL FUNCTIONS
# ============================================================================

def four_peak_model(x, d, r, kmercov, bias, length, k):
    """
    4-peak genome model for k-mer histogram fitting

    The model considers:
    - Homozygous k-mers (1x)
    - Heterozygous k-mers (2x)
    - Repeat k-mers (3x)
    - Heterozygous repeat k-mers (4x)

    Parameters:
    -----------
    x : array
        Coverage values
    d : float
        Heterozygosity rate
    r : float
        Repeat rate
    kmercov : float
        K-mer coverage
    bias : float
        Sequencing bias
    length : float
        Genome length
    k : int
        K-mer size

    Returns:
    --------
    array : Predicted k-mer frequencies
    """
    # Probability that a random k-mer is from a unique region
    p_unique = (1 - r) ** k

    # Component 1: Homozygous unique k-mers (1x coverage)
    comp1 = (2 * (1 - d) * (1 - (1 - r) ** k) +
             2 * d * (1 - (1 - r) ** k) ** 2 +
             2 * d * ((1 - r) ** k) * (1 - ((1 - r) ** k)))

    # Component 2: Heterozygous unique k-mers (2x coverage)
    comp2 = ((1 - d) * ((1 - r) ** k) +
             d * (1 - (1 - r) ** k) ** 2)

    # Component 3: Homozygous repeat k-mers (3x coverage)
    comp3 = 2 * d * ((1 - r) ** k) * (1 - ((1 - r) ** k))

    # Component 4: Heterozygous repeat k-mers (4x coverage)
    comp4 = d * (1 - r) ** (2 * k)

    # Negative binomial distribution for each coverage level
    nb1 = nbinom.pmf(x, kmercov / bias, kmercov / (kmercov + kmercov / bias))
    nb2 = nbinom.pmf(x, kmercov * 2 / bias, kmercov * 2 / (kmercov * 2 + kmercov * 2 / bias))
    nb3 = nbinom.pmf(x, kmercov * 3 / bias, kmercov * 3 / (kmercov * 3 + kmercov * 3 / bias))
    nb4 = nbinom.pmf(x, kmercov * 4 / bias, kmercov * 4 / (kmercov * 4 + kmercov * 4 / bias))

    # Combine components
    y = (comp1 * nb1 * length +
         comp2 * nb2 * length +
         comp3 * nb3 * length +
         comp4 * nb4 * length)

    return y


class ModelWrapper:
    """Wrapper for curve_fit to pass k parameter"""
    def __init__(self, k):
        self.k = k

    def __call__(self, x, d, r, kmercov, bias, length):
        return four_peak_model(x, d, r, kmercov, bias, length, self.k)


# ============================================================================
# MODEL FITTING FUNCTIONS
# ============================================================================

def fit_4peak_model(x, y, k, est_kmercov, est_length, logger=None):
    """
    Fit 4-peak model to k-mer histogram

    Parameters:
    -----------
    x, y : arrays
        Coverage and frequency data
    k : int
        K-mer size
    est_kmercov : float
        Estimated k-mer coverage
    est_length : float
        Estimated genome length
    logger : logging.Logger
        Logger instance

    Returns:
    --------
    popt : array or None
        Optimized parameters [d, r, kmercov, bias, length]
    """
    if logger:
        logger.debug(f"拟合4峰模型|Fitting 4-peak model: kmercov={est_kmercov:.2f}, length={est_length:.0f}")

    # Initial parameter estimates
    p0 = [0.01, 0.5, est_kmercov, 0.5, est_length]  # [d, r, kmercov, bias, length]

    # Parameter bounds
    bounds = ([0, 0, 1, 0.1, est_length * 0.1],
              [1, 1, est_kmercov * 3, 5, est_length * 10])

    try:
        wrapper = ModelWrapper(k)
        popt, pcov = curve_fit(
            wrapper, x, y,
            p0=p0,
            bounds=bounds,
            maxfev=MAX_ITERATIONS * 1000,
            ftol=1e-10,
            xtol=1e-10
        )

        if logger:
            logger.debug(f"模型收敛|Model converged: d={popt[0]:.4f}, r={popt[1]:.4f}, "
                        f"kmercov={popt[2]:.2f}, bias={popt[3]:.2f}, length={popt[4]:.0f}")

        return popt

    except (RuntimeError, ValueError) as e:
        if logger:
            logger.debug(f"模型拟合失败|Model fitting failed: {e}")
        return None


def score_model(kmer_hist, model_params, k):
    """
    Score model by computing residual errors after excluding sequencing errors

    Parameters:
    -----------
    kmer_hist : tuple (x, y)
        Original k-mer histogram
    model_params : array
        Model parameters [d, r, kmercov, bias, length]

    Returns:
    --------
    dict : Scoring metrics
    """
    x, y = kmer_hist

    # Predict with model
    pred = four_peak_model(x, *model_params, k)

    # Compute error rate (kmers unexplained by model through first peak)
    kcov_floor = int(max(model_params[2] - 2 * 0.1, 1))  # Rough estimate
    error_cutoff_idx = np.where(x == kcov_floor)[0]

    if len(error_cutoff_idx) > 0:
        error_cutoff_idx = error_cutoff_idx[0]
    else:
        error_cutoff_idx = int(len(x) * 0.1)

    # Truncate errors
    error_kmers = y[:error_cutoff_idx] - pred[:error_cutoff_idx]
    error_kmers[error_kmers < 0] = 0

    # Total residual error
    total_error = np.sum(error_kmers)
    total_kmers = np.sum(y)

    # Percentage of kmers modeled
    pct_modeled = 100 * (1 - total_error / total_kmers)

    return {
        'all': pct_modeled,
        'total_error': total_error,
        'total_kmers': total_kmers
    }


# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

def load_histogram(histo_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load k-mer histogram file

    Format: coverage frequency (space-separated)

    Parameters:
    -----------
    histo_file : str
        Path to histogram file

    Returns:
    --------
    x, y : arrays
        Coverage and frequency
    """
    # Read histogram file
    data = pd.read_csv(histo_file, sep=' ', header=None, names=['coverage', 'frequency'])

    # Remove last row (often empty or invalid)
    data = data.iloc[:-1]

    # Skip coverage 0 if present
    if data.iloc[0]['coverage'] == 0:
        data = data.iloc[1:]

    x = data['coverage'].values
    y = data['frequency'].values

    return x, y


def estimate_initial_params(x, y, k, read_length, max_cov_index=None):
    """
    Estimate initial parameters for model fitting

    Parameters:
    -----------
    x, y : arrays
        Coverage and frequency
    k : int
        K-mer size
    read_length : int
        Read length
    max_cov_index : int
        Maximum coverage index to consider

    Returns:
    --------
    dict : Initial parameter estimates
    """
    # Truncate if max_cov_index specified
    if max_cov_index is not None and max_cov_index < len(x):
        x = x[:max_cov_index]
        y = y[:max_cov_index]

    # Find peak (excluding error region)
    error_region = x < TYPICAL_ERROR
    main_peak_idx = np.argmax(y[~error_region]) + np.sum(error_region)

    # Estimate k-mer coverage
    est_kmercov = x[main_peak_idx]

    # Total number of k-mers
    total_kmers = np.sum(y)

    # Number of reads
    num_reads = total_kmers / (read_length - k + 1)

    # Estimated genome length
    est_length = int(num_reads * read_length / est_kmercov)

    return {
        'kmercov': est_kmercov,
        'length': est_length,
        'max_cov_index': max_cov_index
    }


# ============================================================================
# RESULTS REPORTING
# ============================================================================

def report_results(x, y, k, model_params, output_dir: str, logger=None):
    """
    Generate results: plots and summary files

    Parameters:
    -----------
    x, y : arrays
        Full k-mer histogram
    k : int
        K-mer size
    model_params : array
        Best fit model parameters [d, r, kmercov, bias, length]
    output_dir : str
        Output directory
    logger : logging.Logger
        Logger instance
    """
    if logger:
        logger.info("生成结果文件|Generating result files...")

    d, r, kmercov, bias, length = model_params

    # Compute model predictions
    x_plot = np.arange(1, max(x) + 1)
    y_pred = four_peak_model(x_plot, *model_params, k)

    # Calculate derived metrics
    unique_length = length * (1 - r) ** k
    repeat_length = length - unique_length
    total_length = length / (1 - r)  # Corrected for repeats

    # Calculate error rate
    error_cutoff = int(max(kmercov - 2 * 0.5, 1))
    error_cutoff_idx = np.where(x == error_cutoff)[0]
    if len(error_cutoff_idx) > 0:
        error_cutoff_idx = error_cutoff_idx[0]
        error_kmers = y[:error_cutoff_idx] - four_peak_model(x[:error_cutoff_idx], *model_params, k)
        error_kmers[error_kmers < 0] = 0
        error_rate = error_kmers.sum() / y.sum()
    else:
        error_rate = 0.01

    # Calculate model fit
    residuals = y - four_peak_model(x, *model_params, k)
    model_fit = 100 * (1 - np.abs(residuals).sum() / y.sum())

    # ----------------------------------------------------------------------
    # Generate plots
    # ----------------------------------------------------------------------
    if logger:
        logger.info("生成图表|Generating plots...")

    # Linear scale plot
    fig, ax = plt.subplots(figsize=(PLOT_SIZE/100, PLOT_SIZE/100), dpi=RESOLUTION)

    ax.set_xlim(0, max(x) * 1.1)
    ax.set_ylim(0, max(y) * 1.1)
    ax.fill_between([0, max(x) * 1.1], 0, max(y) * 1.1, color=COLOR_BGCOLOR)

    ax.plot(x, y, color=COLOR_HIST, linewidth=0.5)
    ax.plot(x_plot, y_pred, color=COLOR_4PEAK, linewidth=2, label='4-peak model')
    ax.axvline(x=error_cutoff, color=COLOR_COVTHRES, linestyle='--', linewidth=1.5, label=f'Error cutoff={error_cutoff}')

    ax.set_xlabel('Coverage', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'GenomeScope Profile\nk={k}', fontsize=14)
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'plot.png'), dpi=RESOLUTION)
    plt.close()

    # Log scale plot
    fig, ax = plt.subplots(figsize=(PLOT_SIZE/100, PLOT_SIZE/100), dpi=RESOLUTION)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-10, max(x) * 10)
    ax.set_ylim(1e-10, max(y) * 10)
    ax.fill_between([1e-10, max(x) * 10], 1e-10, max(y) * 10, color=COLOR_BGCOLOR)

    ax.plot(x, y, color=COLOR_HIST, linewidth=0.5)
    ax.plot(x_plot, y_pred, color=COLOR_4PEAK, linewidth=2)
    ax.axvline(x=error_cutoff, color=COLOR_COVTHRES, linestyle='--', linewidth=1.5)

    ax.set_xlabel('Coverage', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'GenomeScope Profile (log)\nk={k}', fontsize=14)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'plot.log.png'), dpi=RESOLUTION)
    plt.close()

    # ----------------------------------------------------------------------
    # Write summary file
    # ----------------------------------------------------------------------
    if logger:
        logger.info("写入汇总文件|Writing summary file...")

    summary_file = os.path.join(output_dir, 'summary.txt')
    with open(summary_file, 'w') as f:
        f.write("GenomeScope version 2.0 (Python)\n")
        f.write(f"k = {k}\n")
        f.write("\n")
        f.write(f"{'Property':<25} {'Min':>18} {'Max':>18}\n")
        f.write("-" * 65 + "\n")

        f.write(f"{'Heterozygosity':<25} {d*100:>17.3f}% {d*100:>17.3f}%\n")
        f.write(f"{'Genome Haploid Length':<25} {total_length:>17.0f} {total_length:>17.0f}\n")
        f.write(f"{'Genome Repeat Length':<25} {repeat_length:>17.0f} {repeat_length:>17.0f}\n")
        f.write(f"{'Genome Unique Length':<25} {unique_length:>17.0f} {unique_length:>17.0f}\n")
        f.write(f"{'Model Fit':<25} {model_fit:>17.3f}% {model_fit:>17.3f}%\n")
        f.write(f"{'Read Error Rate':<25} {error_rate*100:>17.3f}% {error_rate*100:>17.3f}%\n")

    # ----------------------------------------------------------------------
    # Write model file
    # ----------------------------------------------------------------------
    model_file = os.path.join(output_dir, 'model.txt')
    with open(model_file, 'w') as f:
        f.write(f"kmercov\t{kmercov:.2f}\n")

    if logger:
        logger.info(f"模型收敛|Model converged: het={d:.4f}, kcov={kmercov:.2f}, err={error_rate:.3f}")


# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

def run_genomescope_analysis(histo_file: str, k: int, read_length: int,
                            output_dir: str, max_coverage: Optional[int] = None,
                            logger: Optional[logging.Logger] = None) -> Optional[float]:
    """
    Run GenomeScope analysis

    Parameters:
    -----------
    histo_file : str
        Path to k-mer histogram file (absolute path)
    k : int
        K-mer size
    read_length : int
        Read length
    output_dir : str
        Output directory (absolute path)
    max_coverage : int
        Maximum coverage to model
    logger : logging.Logger
        Logger instance

    Returns:
    -------
    kcov : float or None
        K-mer coverage value
    """
    if k > read_length:
        if logger:
            logger.error(f"K ({k}) 不能大于读长|K cannot be greater than read length ({read_length})")
        return None

    if logger:
        logger.info(f"GenomeScope分析|GenomeScope analyzing: {histo_file}")
        logger.info(f"参数|Parameters: k={k}, readlen={read_length}, outdir={output_dir}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize progress
    progress_file = os.path.join(output_dir, "progress.txt")
    with open(progress_file, 'w') as f:
        f.write("starting\n")

    try:
        # Load histogram
        if logger:
            logger.info("加载k-mer直方图|Loading k-mer histogram...")
        x, y_orig = load_histogram(histo_file)

        # Determine max coverage index
        if max_coverage is not None:
            max_cov_index = np.sum(x <= max_coverage)
        else:
            max_cov_index = len(x)

        # Find error trim point (local minimum in error region)
        error_region = x < TYPICAL_ERROR
        if np.any(error_region):
            start = np.argmin(y_orig[:TYPICAL_ERROR]) + 1
        else:
            start = 1

        # Main optimization loop
        best_model = None
        best_score = float('inf')

        for round_num in range(NUM_ROUNDS):
            if logger:
                logger.info(f"轮次 {round_num + 1}/{NUM_ROUNDS}|Round {round_num + 1}/{NUM_ROUNDS}: 修剪到覆盖度|trimming to coverage {start}")

            # Trim data
            trim_idx = np.where(x >= start)[0]
            if len(trim_idx) == 0:
                break

            trim_start = trim_idx[0]
            x_trim = x[trim_start:]
            y_trim = y_orig[trim_start:]

            # Estimate initial parameters
            init_params = estimate_initial_params(x_trim, y_trim, k, read_length, max_cov_index - trim_start)

            # Fit model
            model_params = fit_4peak_model(
                x_trim, y_trim, k,
                init_params['kmercov'],
                init_params['length'],
                logger
            )

            # Evaluate model
            if model_params is not None:
                score = score_model((x_trim, y_trim), model_params, k)
                if logger:
                    logger.info(f"收敛|Converged. 得分|Score: {score['all']:.2f}%")

                # Update best model
                if score['all'] < best_score:
                    best_model = model_params
                    best_score = score['all']
                    if logger:
                        logger.info("新的最佳模型|New best model!")

                with open(progress_file, 'a') as f:
                    f.write(f"round {round_num}: converged. score: {score['all']:.2f}\n")
            else:
                if logger:
                    logger.warning("未收敛|Failed to converge")
                with open(progress_file, 'a') as f:
                    f.write(f"round {round_num}: unconverged\n")

            # Shift error trim point for next round
            start += START_SHIFT

        # Check if we found a valid model
        if best_model is None:
            if logger:
                logger.error("所有轮次均未找到有效模型|Failed to fit model in all rounds")
            with open(progress_file, 'a') as f:
                f.write("unconverged\n")
            return None

        # Report results with full histogram
        report_results(x, y_orig, k, best_model, output_dir, logger)

        with open(progress_file, 'a') as f:
            f.write("done\n")

        # Return k-mer coverage
        kcov = best_model[2]  # kmercov is at index 2
        return kcov

    except Exception as e:
        if logger:
            logger.error(f"GenomeScope分析失败|GenomeScope analysis failed: {str(e)}", exc_info=True)
        return None
