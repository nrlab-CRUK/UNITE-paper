import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# read id list from a csv file: /home/nrlab/wang04/ulyses/models/bill/roc_curve_efficient_net_s_tf_0.01_0.03/TP_99spe/repeat_2_fold_1_TP_samples_99spec.csv
import pandas as pd
id_list = pd.read_csv(
    '/home/nrlab/wang04/ulyses/models/bill/roc_curve_efficient_net_s_tf_0.01_0.03/TP_99spe/repeat_2_fold_1_TP_samples_99spec.csv')

# sort the id list and get unique id list
id_list = id_list.sort_values(by='TP')
id_list = id_list.drop_duplicates(subset='TP', keep='first')

# get the


def plot_channels_separately(tensor, label_selected, filename='output.pdf'):
    num_samples = tensor.shape[0]
    num_channels = tensor.shape[3]

    # Assuming A4 paper size in inches
    fig_width_inches = 11.7
    fig_height_inches = 8.3

    with PdfPages(filename) as pdf:
        num_plots = len(label_selected)

        # Calculate the number of pages needed
        num_pages = (num_plots + 8) // 9  # Each page has 9 plots

        for page in range(num_pages):
            fig, axes = plt.subplots(3, 3, figsize=(
                fig_width_inches, fig_height_inches))

            # Ensure axes is a 2D array
            if num_pages == 1:
                axes = [axes]
            else:
                axes = axes.reshape(1, -1)

            for i in range(9):
                plot_index = page * 9 + i
                if plot_index < num_plots:
                    sample_index = label_selected[plot_index]
                    for j in range(num_channels):
                        row = i // 3
                        col = i % 3
                        # Assuming grayscale images, change cmap if necessary
                        axes[0][i].imshow(
                            tensor[sample_index, :, :, j], cmap='gray')

                        # Add sample ID in orange color
                        axes[0][i].text(0.5, 0.5, f"{sample_index + 1}", horizontalalignment='center',
                                        verticalalignment='center', fontsize=10, color='orange')

                        axes[0][i].axis('off')

            plt.subplots_adjust(wspace=0, hspace=0)
            pdf.savefig(fig)
            plt.close()

# Example usage:
# Assuming tensor is your 4D tensor with dimensions 1000x464x201x3
# and label_selected is a list of indices indicating which samples to include
# plot_channels_separately(tensor, label_selected, filename='output.pdf')
