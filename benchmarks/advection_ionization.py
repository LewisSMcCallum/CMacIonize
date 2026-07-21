#!/usr/bin/env python3

"""Check ionization advection against periodic translation and make a plot."""

import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PC_IN_M = 3.08567758e16
BOX_SIZE_PC = 10.0
INITIAL_LEFT_PC = -2.5
INITIAL_RIGHT_PC = 2.5
BACKGROUND_FRACTION = 0.1
TOP_FRACTION = 0.9
VELOCITY_M_PER_S = 1.0e4
END_TIME_S = 0.1e6 * 365.25 * 24.0 * 3600.0

# The discontinuous profile is intentionally diffused by first-order scalar
# advection. These tolerances detect a missing, reversed, or non-periodic shift
# while allowing the expected diffusion at this inexpensive resolution.
MAXIMUM_L1_ERROR = 0.06
MAXIMUM_MEAN_ERROR = 0.01
FRACTION_TOLERANCE = 1.0e-10


def load_profile(filename):
    """Return the x coordinates and transverse-mean neutral fractions."""
    data = np.loadtxt(filename, comments="#")
    x, inverse = np.unique(data[:, 0] / PC_IN_M, return_inverse=True)
    counts = np.bincount(inverse)
    fraction = np.bincount(inverse, weights=data[:, 5]) / counts
    return x, fraction


def analytic_cell_average(x):
    """Cell-average the exactly translated periodic square profile."""
    dx = np.mean(np.diff(x))
    displacement = VELOCITY_M_PER_S * END_TIME_S / PC_IN_M
    translated_left = INITIAL_LEFT_PC + displacement
    translated_right = INITIAL_RIGHT_PC + displacement
    overlap = np.zeros_like(x)

    # Include periodic copies of the translated high-fraction interval.
    for periodic_offset in (-BOX_SIZE_PC, 0.0, BOX_SIZE_PC):
        left = translated_left + periodic_offset
        right = translated_right + periodic_offset
        overlap += np.maximum(
            0.0, np.minimum(x + 0.5 * dx, right) -
            np.maximum(x - 0.5 * dx, left)
        )
    overlap = np.minimum(overlap, dx)
    return BACKGROUND_FRACTION + (
        TOP_FRACTION - BACKGROUND_FRACTION
    ) * overlap / dx


def main():
    x_initial, initial = load_profile("advection_ionization_000.txt")
    x_final, final = load_profile("advection_ionization_001.txt")
    if not np.allclose(x_initial, x_final, rtol=0.0, atol=1.0e-10):
        raise RuntimeError("Initial and final snapshot coordinates differ.")

    analytic = analytic_cell_average(x_final)
    residual = final - analytic
    l1_error = np.mean(np.abs(residual))
    mean_error = abs(np.mean(final) - np.mean(initial))
    in_range = (
        np.min(final) >= BACKGROUND_FRACTION - FRACTION_TOLERANCE and
        np.max(final) <= TOP_FRACTION + FRACTION_TOLERANCE
    )
    passed = (
        l1_error <= MAXIMUM_L1_ERROR and
        mean_error <= MAXIMUM_MEAN_ERROR and in_range
    )

    status = "PASS" if passed else "FAIL"
    summary = (
        f"status: {status}\n"
        f"mean absolute error: {l1_error:.8f} "
        f"(limit {MAXIMUM_L1_ERROR:.8f})\n"
        f"mean fraction change: {mean_error:.8f} "
        f"(limit {MAXIMUM_MEAN_ERROR:.8f})\n"
        f"final fraction range: [{np.min(final):.8f}, {np.max(final):.8f}]\n"
    )
    print(summary, end="")
    with open("advection_ionization_metrics.txt", "w", encoding="utf-8") as output:
        output.write(summary)

    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})
    axes[0].plot(x_initial, initial, "k--", label="initial")
    axes[0].plot(x_final, analytic, "b-", label="analytic cell average")
    axes[0].plot(x_final, final, "ro-", markersize=3,
                 label="CMacIonize")
    axes[0].set_ylabel("neutral H fraction")
    axes[0].set_ylim(0.0, 1.0)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="best")
    axes[0].set_title(f"Ionization-state advection: {status}")

    axes[1].axhline(0.0, color="k", linewidth=0.8)
    axes[1].plot(x_final, residual, "k.-")
    axes[1].set_xlabel("x (pc)")
    axes[1].set_ylabel("error")
    axes[1].grid(True, alpha=0.3)
    axes[1].text(
        0.02, 0.05,
        f"L1={l1_error:.4f}; mean change={mean_error:.4f}",
        transform=axes[1].transAxes,
    )
    figure.tight_layout()
    figure.savefig("advection_ionization.png", dpi=150)
    plt.close(figure)

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
