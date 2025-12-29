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
import math

#################################

def rad_to_deg(rad):
    return rad * 180.0 / math.pi

# Load space-separated BRDF data from files and plot them

# Load data, skipping comment lines
data = np.loadtxt('../data/retro-3M-jaune.dat', comments='#')

# Extract columns
theta_L = data[:, 0]      # First column
theta_V = data[:, 1]      # Second column
phi_diff = data[:, 2]     # Third column
brdf_values_times_cosine = data[:, 3]  # Fourth column

# Quantize values to differ at most by 1.0e-10
quantization_precision = 1.0e-10
theta_L = np.round(theta_L / quantization_precision) * quantization_precision
theta_V = np.round(theta_V / quantization_precision) * quantization_precision
phi_diff = np.round(phi_diff / quantization_precision) * quantization_precision
brdf_values = np.round(brdf_values_times_cosine / quantization_precision) * quantization_precision

# Filter out rows where abs(phi_diff) > 0.0001
mask = np.abs(phi_diff) <= 0.0001
theta_L = theta_L[mask]
theta_V = theta_V[mask]
phi_diff = phi_diff[mask]
brdf_values_times_cosine = brdf_values_times_cosine[mask]



# Create a dictionary mapping theta_L to list of (theta_V, brdf_value) pairs
brdf_dict_simple = {}
for tL, tV, brdf_times_cos in zip(theta_L, theta_V, brdf_values_times_cosine):
    if tL not in brdf_dict_simple:
        brdf_dict_simple[tL] = []
    print('tL:', tL, 'cos(tL):', np.cos(tL), 'brdf_times_cos:', brdf_times_cos)
    brdf_dict_simple[tL].append((tV, brdf_times_cos / np.cos(tL)))  # Store BRDF value (not multiplied by cosine)


print("BRDF Data Summary:")
print(f"Number of unique theta_L values: {len(brdf_dict_simple)}")

# Path pattern to match CSV files
calculator_files = glob.glob('brdf_mrm_theta_i_*_phi_*.csv')
calculator_files.sort()

calc_brdf_dict = {}

for file in calculator_files:

    calc_data = pd.read_csv(file)
    #print(f"Columns in {file}: {calc_data.columns.tolist()}")
    #print(calc_data)

    # Extract theta_i and phi from the filename
    parts = file.replace('.csv', '').split('_')

    theta_i = parts[4]
    phi = parts[6]
    calc_roughness = float(parts[8])  # roughness, assume constant for all loaded data

    print(f"Processing file: {file} with theta_i: {theta_i}, phi: {phi}, roughness: {calc_roughness}")

    # Convert theta_i to radians and calculate cos(theta_i)
    theta_i_rad = np.radians(float(theta_i))
    cos_theta_i = np.cos(theta_i_rad)

    calc_brdf_dict[float(theta_i)] = [calc_data['theta_o'], calc_data['MRM']]




# Create plots for each unique theta_L value
fig, axes = plt.subplots(nrows=int(np.ceil(len(brdf_dict_simple) / 3)), ncols=3, figsize=(15, 5 * int(np.ceil(len(brdf_dict_simple) / 3))))
axes = axes.flatten() if len(brdf_dict_simple) > 1 else [axes]

for idx, (tL, data_points) in enumerate(sorted(brdf_dict_simple.items())):
    # Sort data points by theta_V
    data_points.sort(key=lambda x: x[0])
    theta_V_vals = [point[0] for point in data_points]
    brdf_vals = [point[1] for point in data_points]

    tL_deg = rad_to_deg(tL)
    tL_deg = np.round(tL_deg / 0.001) * 0.001
    print('tL_deg:', tL_deg)
    # Convert theta_V values to degrees for plotting
    theta_V_vals_deg = [rad_to_deg(tv) for tv in theta_V_vals]

    # Plot BRDF values vs theta_V (in degrees)
    axes[idx].clear()  # Clear previous plot
    axes[idx].plot(theta_V_vals_deg, brdf_vals, 'o-', color='black', alpha=0.25, linewidth=1.5, markersize=0)
    axes[idx].plot(theta_V_vals_deg, brdf_vals, 'o-', color='black', alpha=0.75, linewidth=0, markersize=3)

    axes[idx].set_xlabel('$\\theta_V$ (degrees)', fontsize=16)
    axes[idx].set_ylabel('BRDF value', fontsize=16)
    axes[idx].set_title(f'$\\theta_L = {tL_deg:.2f}^\\circ$', fontsize=16)
    axes[idx].set_xlim(-90, 90)
    axes[idx].set_ylim(0.001, 50)
    axes[idx].set_yscale('log')
    axes[idx].grid(True, alpha=0.3)

    # Set custom x-axis ticks
    custom_ticks = np.arange(-90, 91, 45)  # Creates ticks from -90 to 90 with a step of 22.5
    axes[idx].set_xticks(custom_ticks)

    # Set custom x-axis tick labels with a degree symbol
    custom_tick_labels = [f"{int(tick)}°" for tick in custom_ticks]
    axes[idx].set_xticklabels(custom_tick_labels)

    axes[idx].axvline(x=tL_deg, color='grey', linestyle=':', linewidth=1.5, alpha=0.7)
    axes[idx].axvline(x=-tL_deg, color='grey', linestyle=':', linewidth=1.5, alpha=0.7)

    axes[idx].plot(calc_brdf_dict[tL_deg][0], calc_brdf_dict[tL_deg][1], 'r--', color='green', label='Calculator MRM', linewidth=1.5)

# Hide unused subplots
for idx in range(len(brdf_dict_simple), len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
plt.savefig('brdf_theta_L_plots.pdf')
plt.show()


