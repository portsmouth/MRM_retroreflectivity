# Prerequisites:
#
#     python3 -m pip install matplotlib pandas
#
# Running this script:
#
#     python3 brdf_plotting.py
#
# Viewing the output:
#
#     open brdf_plots.pdf
#

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import colorsys

leg_locs = ['lower left', 'lower right', 'upper left',
            'lower left', 'lower right', 'upper left',
            'lower left', 'lower right', 'upper left']


def plot_brdf_data():

    # Path pattern to match CSV files
    files = glob.glob('brdf_theta_i_*_phi_*.csv')
    files.sort()

    # Prepare a 3x3 subplot
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    axes = axes.flatten()

    # Helper function to determine the subplot index based on theta_i and phi values
    def get_index(theta_i, phi):
        theta_index = {'45': 0, '80': 1}
        phi_index   = {'0': 0, '45': 2}
        return phi_index[phi] + theta_index[theta_i]

    for file in files:
        data = pd.read_csv(file)

        # Extract theta_i and phi from the filename
        parts = file.replace('.csv', '').split('_')
        theta_i = parts[3]
        phi = parts[5]
        roughness = float(parts[7])  # roughness

        print(f"Processing file: {file} with theta_i: {theta_i}, phi: {phi}, roughness: {roughness}")

        # Convert theta_i to radians and calculate cos(theta_i)
        theta_i_rad = np.radians(float(theta_i))
        cos_theta_i = np.cos(theta_i_rad)

        # Adjust the BRDF values by multiplying by cos(theta_i)
        #data['Lambert'] *= cos_theta_i
        #data['ON'] *= cos_theta_i
        #data['QON'] *= cos_theta_i
        #data['footnoteQON'] *= cos_theta_i
        #data['fujiiQON'] *= cos_theta_i
        #data['FON'] *= cos_theta_i
        #data['EON'] *= cos_theta_i
        #data['GGX'] *= cos_theta_i
        #data['MRM'] *= cos_theta_i

        # Calculate the correct subplot index based on extracted theta_i and phi
        index = get_index(theta_i, phi)

        # Set custom x-axis ticks
        custom_ticks = np.arange(-90, 91, 45)  # Creates ticks from -90 to 90 with a step of 22.5
        axes[index].set_xticks(custom_ticks)

        # Set custom x-axis tick labels with a degree symbol
        custom_tick_labels = [f"{int(tick)}°" for tick in custom_ticks]
        axes[index].set_xticklabels(custom_tick_labels)

        # Vertical lines at specified theta_o values
        for theta_line in [-80, -45, 0, 45, 80]:
            axes[index].axvline(x=theta_line, color='lightgray', linestyle='dotted', linewidth=0.7)

        # Setting y-ticks for every 0.05 and adding horizontal grid lines
        #axes[index].set_yticks(np.arange(0, 0.28, 0.05))

        # Enable grid only for y-axis with specific styling, avoid x-axis grid
        axes[index].grid(which='both', axis='y', color='lightgray', linestyle='dotted', linewidth=0.7)

        # Plot each BRDF curve in the correct subplot with a thinner line width
        if roughness < 0.5:
            linestyle = 'dashed'
            show_legend = False
        else:
            linestyle = 'solid'
            show_legend = True

        axes[index].plot(data['theta_o'], data['Lambert'], label='Lambert' if show_legend else None, color='gray',   linewidth=1.0, linestyle=linestyle)
        axes[index].plot(data['theta_o'], data['EON'],     label='EON'     if show_legend else None, color='orange', linewidth=1.0, linestyle=linestyle)
        axes[index].plot(data['theta_o'], data['GGX'],     label='GGX'     if show_legend else None, color='blue',   linewidth=1.0, linestyle=linestyle)
        axes[index].plot(data['theta_o'], data['MRM'],     label='MRM'     if show_legend else None, color='green',  linewidth=1.0, linestyle=linestyle)

        # Set the title for each subplot
        axes[index].set_title(f"$\\theta_i: {theta_i}^\\circ$, $\\phi_o: {phi}^\\circ$, $r: {roughness}$", fontsize=12)
        axes[index].set_xlabel(f"$\\theta_o$", fontsize=16, fontweight='bold')
        axes[index].set_ylabel('BRDF value', fontsize=16)
        axes[index].set_xlim(-90, 90)

        axes[index].set_yscale('log')
        axes[index].set_ylim(0.00001, 100000)

        # Add legend
        if show_legend:
            axes[index].legend(loc='upper center', fontsize=9)

    # Overall plot title with LaTeX formatting, disabled for the paper version of the plot
    # plt.suptitle('Comparison of diffuse BRDFs with $\\rho = 0.8$, $\\sigma = \\pi/4$, and $r = 1/2$', fontsize=16)

    plt.tight_layout()
    plt.savefig('BRDF_Plots.pdf')
    plt.show()
    plt.close(fig)

if __name__ == '__main__':
    plot_brdf_data()
