import os
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.stats as stats
import statsmodels.stats.multitest as multitest


def statistics_calculation(
    final_results_dataframe, save_loc, multiple_correction_method, plot_dpi=600
):
    """Performs the statistical analysis for the pipeline.

    Args:
        final_results_file (str): File path for the final results .csv file
        save_loc (str): Directory to save the results of the statistical analysis
        plot_dpi (int, optional): Resolution of all resulting plots. Defaults to 600.
        multiple_correction_method (str, optional): Multiple correction method for the calculation of final p-values. Defaults to 'fdr_by'.
    """

    # Load 'final' data and extract normalized scores
    final_results_dataframe.drop(["evalue"], axis=1, inplace=True)

    scores = final_results_dataframe["Normalized_score"].to_numpy()

    # Create and save a normalized histogram of the normalized scores
    # With a KDE estimator line
    save_normalized_histogram(scores, save_loc, plot_dpi)

    # Calculate and append robust Z-scores to new column
    final_results_dataframe["Robust_Zscores"] = calculate_robust_z_scores(scores)

    # Fit gumbel curve to final normalized scores
    gumbel_params = fit_gumbel(scores, save_loc, plot_dpi)

    # Calculate corrected p-values and save as results column
    final_results_dataframe["Corrected_pvalues"] = calculate_gumbel_p_values(
        scores, gumbel_params, multiple_correction_method
    )

    # save the final dataframe

    # df_final_results.to_csv(save_loc + "final_results.csv")

    print(f"Finished statistical analysis")

    return final_results_dataframe


def save_normalized_histogram(scores, save_path, plot_dpi=600):
    """Create and save a normalized histogram of the normalized alignment scores, with a KDE estimator line.

    Args:
        scores (numpy.ndarray): A NumPy array of Normalized_score extracted from the final results dataframe.
        save_path (str): Location to save the plot.
        plot_dpi (int, optional): DPI quality of saved plots. Defaults to 600.

    Returns:
            Boolean: True for successfully creating and saving the histogram, False for any failure to do so.


    """

    print("Creating and saving histogram of normalized scores for final results...")
    try:
        hist_plot = sns.histplot(
            scores, stat="density", kde=True, label="Kernel Density Estimator"
        )

        hist_plot.set_title("Histogram plot of Normalized scores")
        hist_plot.set_xlabel("Normalized Alignment Scores")
        hist_plot.legend()

        # Create a custom legend for the KDE line
        kde_line = hist_plot.lines[
            0
        ]  # Typically, the KDE line is the first line object in the plot

        # Create a legend with a line style
        legend_lines = [
            Line2D(
                [0],
                [0],
                color=kde_line.get_color(),
                lw=2,
                label="Kernel Density Estimator",
            )
        ]

        hist_plot.legend(handles=legend_lines, loc="best")

        plot_name = "Histogram_normalized_scores.png"
        plt.savefig(save_path + plot_name, dpi=plot_dpi)

        plt.close()

        return True

    except Exception as ee:
        print(f"Saving histogram of Normalized scores failed with:{ee}")

        return False


def calculate_robust_z_scores(scores):
    """Calculates robust z_scores based on a NumPy array of normalized scores

    Args:
        scores (numpy.ndarray): A NumPy array of Normalized_score extracted from the final results dataframe.

    Returns:
        numpy.ndarray: NumPy array of calculated Z-scores
    """

    print("Calculating robust Z-scores...")

    median = np.median(scores)
    mad = (
        np.median(np.abs(scores - median)) * 1.4826
    )  # Normalize to be similar to normal Z-scores, assumes normality of data

    if mad == 0:
        # Avoid division by zero, return inf in this case

        print(
            "WARNING: Mean Absolute Deviation computed as '0', please check the normalized scores."
        )
        return np.inf

    # Calculation of robust Z-scores
    z_scores = (scores - median) / mad

    print("Finished calculating robust Z-scores")

    return z_scores


def fit_gumbel(scores, save_loc, plot_dpi=600):
    """Fits the gumbel distribution to all of the Normalized scores.
    Also creates and saves plots meant to check the success of the fit.

    Args:
        scores (numpy.ndarray):  A NumPy array of Normalized_score extracted from the final results dataframe.
        save_loc (str): Location to save the plot.
        plot_dpi (int, optional): DPI quality of saved plots. Defaults to 600.

    Returns:
        tuple: Tuple of parameters from the Gumbel fit. Structure: (mu, beta)
    """

    # Fit Gumbel distribution to the scores
    params = stats.gumbel_r.fit(scores)

    # Save plots for checking the assumptions of fit
    plot_and_save_gumbel_fit(scores, params, save_loc, plot_dpi)
    save_qq_plot_for_gumbel_fit(scores, params, save_loc, plot_dpi)

    return params


def calculate_gumbel_p_values(scores, params, multiple_correction_method):
    """Calculates p-values based on the fitted Gumbel distribution
      assuming a one-tailed test, ergo H1: normalized score is abnormally high.
      Also performs multiple test correction of the resulting p-values.

    Args:
        scores (numpy.ndarray): A NumPy array of Normalized_score extracted from the final results dataframe.
        params (tuple): Parameters from the gumbel_fit,  (mu, beta)
        multiple_correction_method (str, optional): Method selected for multiple test correction, using the statsmodels.stats.multitest module. Defaults to 'fdr_by' == Storey's Q test.

    Returns:
        numpy.ndarray: A NumPy array of p-values (or q-values) for the normalized scores, based on the Gumbel fit.
    """

    mu, beta = params
    cdf_values = stats.gumbel_r.cdf(scores, loc=mu, scale=beta)

    # Calculate the p-values (1 - CDF)
    # This assumes a one-tailed test, ergo we are only interested in high scores
    p_values = 1 - cdf_values

    from statsmodels.stats.multitest import multipletests

    _, corrected_p_values, _, _ = multipletests(
        p_values, method=multiple_correction_method
    )

    return corrected_p_values


def plot_and_save_gumbel_fit(scores, params, save_loc, plot_dpi=600):
    """Create and save a normalized histogram of the alignment scores, in comparison to the Gumbel fit. Good for an initial check.

    Args:
        scores (numpy.ndarray):  A NumPy array of Normalized_score extracted from the final results dataframe.
        params (tuple): Parameters from the gumbel_fit,  (mu, beta)
        save_loc (str): Location to save the plot.
        plot_dpi (int, optional): DPI quality of saved plots. Defaults to 600.

    Returns:
        Boolean: True for successfully creating and saving the histogram, False for any failure to do so.
    """

    try:
        # Generate the distribution for gumbel
        mu, beta = params
        x = np.linspace(min(scores), max(scores), 10_000)
        pdf_fitted = stats.gumbel_r.pdf(x, loc=mu, scale=beta)

        # Plot the histogram of your scores
        plt.hist(
            scores, density=True, alpha=0.6, color="blue", label="Histogram of scores"
        )

        # Plot the fitted Gumbel PDF
        plt.plot(x, pdf_fitted, "r-", lw=2, label="Fitted Gumbel PDF")

        # Add labels and title
        plt.title("Fitted Gumbel Distribution to Alignment Scores")
        plt.xlabel("Normalized Alignment Score")
        plt.ylabel("Density")
        plt.legend()

        plot_name = "Gumbel_fit_histogram.png"
        plt.savefig(save_loc + plot_name, dpi=plot_dpi)
        plt.close()

        return True

    except Exception as ee:
        print(f"Saving histogram for Gumbel fit failed with:{ee}")
        return False


def save_qq_plot_for_gumbel_fit(scores, params, save_loc, plot_dpi=600):
    """Generates a Q-Q plot for the Gumbel fit. This is the main plot to
      consult in order to determine if the fit is successful. Note that
      if the pipeline actually produced a positive hit, that will show as a
      deviation in the higher scores.

    Args:
        scores (numpy.ndarray):  A NumPy array of Normalized_score extracted from the final results dataframe.
        params (tuple): Parameters from the gumbel_fit,  (mu, beta)
        save_loc (str): Location to save the plot.
        plot_dpi (int, optional): DPI quality of saved plots. Defaults to 600.

    Returns:
        Boolean: True for successfully creating and saving the histogram, False for any failure to do so.
    """

    try:
        stats.probplot(scores, dist="gumbel_r", sparams=params, plot=plt)
        plt.title("Q-Q Plot for Gumbel Fit against Normalized scores")

        plot_name = "QQ_plot_gumbel_fit.png"
        plt.savefig(save_loc + plot_name, dpi=plot_dpi)
        plt.close()

        return True

    except Exception as ee:

        print(f"Saving histogram for Gumbel fit failed with:{ee}")
        return False
