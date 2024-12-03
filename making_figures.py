# file containing functions to be used for plottig results/making figures

import numpy as np
import matplotlib.pyplot as plt

# calculate and plot the difference between two wavefields


def calculate_diff(unscaled_field, scaled_field):
    "pass the data from the fields, p.data[0] and p_s.data[0] for example"
    diff_field = (unscaled_field - scaled_field)
    # Calculate maximum values for color scaling
    vmax_unscaled = np.amax(np.abs(unscaled_field))
    vmax_scaled = np.amax(np.abs(scaled_field))
    vmax_diff = np.amax(np.abs(diff_field))

    # Create a figure with 3 subplots in one row
    fig, ax = plt.subplots(1, 4, figsize=(24, 6))  # 1 row, 3 columns

    # Plot unscaled field
    im1 = ax[0].imshow(unscaled_field, cmap='seismic', vmin=-vmax_unscaled, vmax=vmax_unscaled)
    ax[0].set_title("Unscaled Field")
    ax[0].set_xlabel("X (m)")
    ax[0].set_ylabel("Z (m)")
    fig.colorbar(im1, ax=ax[0], shrink=0.6)

    # Plot scaled field
    im2 = ax[1].imshow(scaled_field, cmap='seismic', vmin=-vmax_scaled, vmax=vmax_scaled)
    ax[1].set_title("Scaled Field")
    ax[1].set_xlabel("X (m)")
    ax[1].set_ylabel("Z (m)")
    fig.colorbar(im2, ax=ax[1], shrink=0.6)

    # Plot the true difference field
    # use the vmax and vmin of either the scaled or unscaled field to normalise the diff against this
    im3 = ax[2].imshow(diff_field, cmap='seismic',  vmin=-vmax_diff, vmax=vmax_diff)
    ax[2].set_title("True difference")
    ax[2].set_xlabel("X (m)")
    ax[2].set_ylabel("Z (m)")
    fig.colorbar(im3, ax=ax[2], shrink=0.6)

    # Plot the normalised difference field
    # use the vmax and vmin of either the scaled or unscaled field to normalise the diff against this
    im4 = ax[3].imshow(diff_field, cmap='seismic',  vmin=-vmax_unscaled, vmax=vmax_unscaled)
    ax[3].set_title("Normalised difference")
    ax[3].set_xlabel("X (m)")
    ax[3].set_ylabel("Z (m)")
    fig.colorbar(im4, ax=ax[3], shrink=0.6)

    # Adjust layout for better spacing
    plt.tight_layout()

    # couple of print statements to show the max/min difference
    print("Maximum value of the difference:", np.amax(diff_field))
    print("Minimum value of the difference:", np.amin(diff_field))

    # Show the plot
    plt.show()
