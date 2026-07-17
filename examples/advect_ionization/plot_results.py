#!/usr/bin/env /usr/bin/python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_profile(filename):
    data = np.loadtxt(filename, comments='#')
    # Column 0: x (m), Column 5: neutral H fraction
    xs = data[:, 0]
    fracs = data[:, 5]

    # Group by unique x coordinates and compute the mean fraction for each x
    unique_xs, indices = np.unique(xs, return_inverse=True)
    mean_fracs = np.zeros_like(unique_xs)
    for i in range(len(unique_xs)):
        mean_fracs[i] = np.mean(fracs[indices == i])

    return unique_xs, mean_fracs

def main():
    # Load simulated profiles
    xs_0, fracs_0 = load_profile('advection_000.txt')
    xs_1, fracs_1 = load_profile('advection_001.txt')

    # Convert x from meters to parsec for nicer presentation
    pc_in_m = 3.08567758e16
    xs_pc_0 = xs_0 / pc_in_m
    xs_pc_1 = xs_1 / pc_in_m

    # Calculate analytical benchmark
    # Initial profile is 0.9 inside [-2.5 pc, 2.5 pc] and 0.1 elsewhere.
    # It advects at vx = 10 km/s for t = 0.1 Myr.
    # Displacement d = vx * t = 10,000 m/s * 0.1e6 * 3.15576e7 s = 3.15576e16 m ~ 1.0227 pc
    # With inflow (open) boundary conditions:
    # Any advected state from x < -5 pc is filled with the background value (0.1).
    d_pc = 1.0227
    analytical_fracs = []
    for x in xs_pc_1:
        x_init = x - d_pc
        if x_init >= -5.0 and x_init <= 5.0:
            # Evaluate initial condition at x_init
            if x_init >= -2.5 and x_init <= 2.5:
                val = 0.9
            else:
                val = 0.1
        else:
            val = 0.1
        analytical_fracs.append(val)

    analytical_fracs = np.array(analytical_fracs)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(xs_pc_0, fracs_0, 'k--', label='Initial Profile (t=0)')
    plt.plot(xs_pc_1, fracs_1, 'ro-', label='Simulated Profile (t=0.1 Myr)')
    plt.plot(xs_pc_1, analytical_fracs, 'b-', label='Analytical Benchmark (t=0.1 Myr)')

    plt.xlabel('x (pc)')
    plt.ylabel('Neutral H Fraction')
    plt.title('Ionization State Advection Test')
    plt.grid(True)
    plt.legend()
    plt.ylim(-0.05, 1.05)

    plt.savefig('advection_comparison.png', dpi=150)
    print("Plot successfully generated and saved to advection_comparison.png")

if __name__ == '__main__':
    main()
