# file containing functions to be used for plottig results/making figures

import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, mpf, sin, pi, fadd, fsub, fmul, fdiv


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

    # calculate absolute percentage change
    diff_range = np.amax(diff_field) - np.amin(diff_field)
    unscaled_range = np.amax(unscaled_field) - np.amin(unscaled_field)
    change = (diff_range/unscaled_range) * 100

    # couple of print statements to show the max/min difference
    print("Maximum value of the difference:", np.amax(diff_field))
    print("Minimum value of the difference:", np.amin(diff_field))
    print("Percentage change: ", change, "%")

    # Show the plot
    plt.show()

# functions for solving wave equation for given set up and plotting results
# Analytical Solution
def analytical_solution(t, x, c, L):
    "Returns analytical solution at a given timestep"
    t_np = float(t)
    return np.cos(2 * np.pi * c * t_np / L) * np.sin(2 * np.pi * x / L)

def solve(prec, L, c, N, T, courant):
    "Solves given system, returning the numerical solution at t0, tN/2 and tN"
    mp.prec = prec
    dx = 2 * L / (N - 1)          # Spatial step size for [-L, L] domain
    courant = courant
    dt = courant * dx / c             # Time step size (CFL condition)
    T = T             # Total simulation time in seconds (period is T=1)
    x = np.linspace(-L, L, N)  # Spatial grid as a NumPy array

    # Total time steps
    Nt = int(T / dt) + 1

    # Initialize the wave field
    u = analytical_solution(0, x, c, L)  # Initial condition u(x, 0) as float
    u_prev = u.copy()                          # Initial velocity is zero: u_t(x, 0) = 0
    u_next = np.zeros(N, dtype=np.float64)     # Placeholder for the next time step

    # Dynamic time indices for plotting and error storage
    plot_times = [0, Nt // 2 - 1, Nt - 1]      # Start, middle, and end of the time range
    u_num_at_t = {}
    errors = []

    # Main time evolution loop
    for n in range(Nt):
        # Store the results for specific time indices
        if n in plot_times:
            u_num_at_t[n] = u.copy()  # Store the wave profile at the specified time indices

        # Compute the next time step using finite difference
        coeff = fdiv(fmul(mpf(c), mpf(dt)), mpf(dx))
        coeff_squared = fdiv(coeff, coeff)

        for i in range(1, N - 1):
            # Compute the finite difference term with mpf arithmetic
            term2 = fadd(fsub(mpf(u[i + 1]), fmul(mpf(2), mpf(u[i]))), mpf(u[i - 1]))
            update = fmul(coeff_squared, term2)  # (c * dt / dx) ** 2 * (term2)
            u_next[i] = float(fadd(fsub(fmul(mpf(2), mpf(u[i])), mpf(u_prev[i])), update))  # Cast result back to float

        # Apply boundary conditions (no need for mpf here)
        u_next[0] = 0.0
        u_next[-1] = 0.0

        # Compute error at this time step (if plot_times include the current time step)
        t = n * dt
        if n in plot_times:
            u_analytic = analytical_solution(t, x, c, L)
            errors.append({
                "time": t,
                "error": max(abs(u[i] - u_analytic[i]) for i in range(N))
                })

        # Update arrays for the next time step
        u_prev = u.copy()
        u = u_next.copy()
    return u_num_at_t


def plot(solutions, L, c, N, T, C):
    "Produces plot of analytical solution above numerical solution"
    dx = 2 * L / (N - 1)          # Spatial step size for [-L, L] domain
    courant = C
    dt = courant * dx / c             # Time step size (CFL condition)
    T = T             # Total simulation time in seconds (period is T=1)
    x = np.linspace(-L, L, N)  # Spatial grid as a NumPy array

    # Total time steps
    Nt = int(T / dt) + 1
    plot_times = [0, Nt // 2 - 1, Nt - 1]
    ## Plotting the analytical solution
    fig1, axs1 = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    for i, idx in enumerate(plot_times):
        t = idx * dt
        # Analytical solution at time t using NumPy variables
        u_analytic = analytical_solution(t, x, c, L)
        axs1[i].plot(x, u_analytic, 'k-', label=f"Analytical t = {float(t)} s")

        # Plot formatting
        axs1[i].set_xlabel("Position (x) [m]")
        axs1[i].set_title(f"Analytical Solution at t = {float(t)} s")
        axs1[i].legend()
        axs1[i].grid(True)

    # Set global y-axis label
    fig1.supylabel("u(x, t)")
    fig1.suptitle("Analytical Solution")
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Plotting the numerical solution
    fig2, axs2 = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    for i, idx in enumerate(plot_times):
        t = idx * dt
        # Retrieve the numerical solution for the mapped time index
        u_num = solutions[idx]
        axs2[i].plot([float(x_i) for x_i in x], [float(u_n) for u_n in u_num], 'r--', label=f"Numerical t = {float(t)} s")

        # Plot formatting
        axs2[i].set_xlabel("Position (x) [m]")
        axs2[i].set_title(f"Numerical Solution at t = {float(t)} s")
        axs2[i].legend()
        axs2[i].grid(True)

    # Set global y-axis label
    fig2.supylabel("u(x, t)")
    fig2.suptitle(f"Numerical Solution (Precision: {mp.prec} bits)")
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.show()

