import os
import logging
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.stats as stats
import statsmodels.stats.multitest as multitest


logger = logging.getLogger(__name__)




def statistics_calculation(final_results_file, save_path, plot_dpi = 600 ):
    """_summary_

    Args:
        final_results_file (_type_): _description_
        save_path (_type_): _description_
        plot_dpi (int, optional): _description_. Defaults to 600.

    Returns:
        _type_: _description_
    """
    
    # Load 'final' data and extract normalized scores
    df_final_results = pd.read_csv(final_results_file, sep = ",")
    df_final_results.drop(['Unnamed: 0', 'evalue'], axis =1, inplace= True)

    scores = df_final_results["Normalized_score"].to_numpy()

    # Create and save a normalized histogram of the normalized score
    # With a KDE estimator line
    save_normalized_histogram(scores, save_path, plot_dpi)

    return scores



def save_normalized_histogram(scores, save_path, plot_dpi = 600):
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
        hist_plot = sns.histplot(scores,
                                stat = "density", kde = True,
                                label="Kernel Density Estimator")

        hist_plot.set_title("Histogram plot of Normalized scores")
        hist_plot.set_xlabel("Normalized Alignment Scores")
        hist_plot.legend()

        # Create a custom legend for the KDE line
        kde_line = hist_plot.lines[0]  # Typically, the KDE line is the first line object in the plot

        # Create a legend with a line style
        legend_lines = [Line2D([0], [0], 
                        color=kde_line.get_color(), 
                        lw=2, label= "Kernel Density Estimator")]
        
        hist_plot.legend(handles=legend_lines, loc='best')

        plt.savefig(save_path, dpi = plot_dpi)

        plt.close()

        return True
    
    except Exception as ee:
        print(f"Saving histogram of Normalized scores failed with:{ee}")

        return False



def p_value_generation(scores,save_loc, plot_dpi = 600 ):

    print("Fitting results and calculating p-values...")

    # Fit Gumbel distribution to the scores
    params = stats.gumbel_r.fit(scores)




    
def plot_and_save_gumbel_fit(scores, params, save_loc, plot_dpi = 600):
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
        plt.hist(scores, density=True, alpha=0.6, color='blue', label='Histogram of scores')

        # Plot the fitted Gumbel PDF
        plt.plot(x, pdf_fitted, 'r-', lw=2, label='Fitted Gumbel PDF')

        # Add labels and title
        plt.title('Fitted Gumbel Distribution to Alignment Scores')
        plt.xlabel('Normalized Alignment Score')
        plt.ylabel('Density')
        plt.legend()

        plt.savefig(save_loc, dpi = plot_dpi)
        plt.close()
        
        return True

    except Exception as ee:
        print(f"Saving histogram for Gumbel fit failed with:{ee}")
        return False


def save_qq_plot_for_gumbel_fit(scores, params, save_loc, plot_dpi = 600):
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
        plt.title('Q-Q Plot for Gumbel Fit against Normalized scores')
        plt.savefig(save_loc, dpi = plot_dpi)
        plt.close()

        return True
    
    except Exception as ee:

        print(f"Saving histogram for Gumbel fit failed with:{ee}")
        return False

    
