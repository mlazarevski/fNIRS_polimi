"""TD-NIRS: Compute Differential Pathlength Factor (DPF) at each optode position."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as stats
import os, sys

sys.path.insert(0, os.path.dirname(__file__))
from config import *


def read_dat(path):
    """Read a TD-NIRS .DAT file → (NREPS, NBINS) float64."""
    with open(path, "rb") as fh:
        fh.read(TD_HEADER)
        curves = np.zeros((TD_NREPS, TD_NBINS), dtype=np.float64)
        for i in range(TD_NREPS):
            fh.read(TD_SUBHDR)
            curves[i] = np.frombuffer(fh.read(TD_NBINS * 2),
                                      dtype=np.uint16).astype(np.float64)
    return curves


def mean_time(curve, t_axis, thresh_frac=0.01):
    """Mean time of flight from a DTOF, thresholded at 1% of peak."""
    mask = curve > curve.max() * thresh_frac
    return np.sum(t_axis[mask] * curve[mask]) / np.sum(curve[mask])


def scholkmann_dpf(wl, age=25):
    """Scholkmann et al. 2013 — age & wavelength-dependent DPF."""
    if age > 50:
        age = 50
    return (223.3 + 0.05624 * age**0.8493
            - 5.723e-7 * wl**3 + 0.001245 * wl**2 - 0.9025 * wl)


def main():
    t_bins = np.arange(TD_NBINS) * TD_PS_BIN  # ps

    # --- Read all files ---
    irf = read_dat(IRF_DAT)
    meas = {}
    for f in MEAS_DATS:
        idx = int(os.path.basename(f).replace("DPF_VmLC", "").replace(".DAT", ""))
        meas[idx] = read_dat(f)

    # --- IRF mean times per wavelength ---
    irf_mt = np.zeros(TD_N_WL)
    for wi in range(TD_N_WL):
        avg = irf[wi * TD_N_REP:(wi + 1) * TD_N_REP].mean(axis=0)
        irf_mt[wi] = mean_time(avg, t_bins)

    # --- Tissue transit times ---
    n_pos = len(meas)
    delta_t = np.zeros((n_pos, TD_N_WL))
    pos_ids = sorted(meas.keys())

    for pi, pid in enumerate(pos_ids):
        for wi in range(TD_N_WL):
            avg = meas[pid][wi * TD_N_REP:(wi + 1) * TD_N_REP].mean(axis=0)
            delta_t[pi, wi] = mean_time(avg, t_bins) - irf_mt[wi]

    # --- DPF for a range of assumed probe distances ---
    print("=" * 70)
    print("TD-NIRS DPF ANALYSIS")
    print("=" * 70)
    print("\nIRF mean times (ps):", np.round(irf_mt, 1))
    print("\nTissue transit times (delta_t in ps):")
    print("{:>5s} {:>7s}  {:>7s} {:>7s} {:>7s} {:>7s}".format(
        "Pos", "Region", "670nm", "730nm", "780nm", "830nm"))
    for i, lbl in enumerate(POS_LABELS):
        print("{:>5s} {:>7s}  {:7.1f} {:7.1f} {:7.1f} {:7.1f}".format(
            lbl, pos_region(lbl),
            delta_t[i, 0], delta_t[i, 1], delta_t[i, 2], delta_t[i, 3]))

    print("\n--- DPF for different probe distances ---")
    for d in [10, 15, 20, 25, 30]:
        dpf = V_TISSUE * delta_t / d
        print("d={:2d}mm: mean DPF = [{:.2f}, {:.2f}, {:.2f}, {:.2f}], overall = {:.2f}".format(
            d, *dpf.mean(axis=0), dpf.mean()))

    D_PROBE = 30.0
    dpf_all = V_TISSUE * delta_t / D_PROBE

    # --- Regional breakdown ---
    regions = ["Left", "Central", "Right"]
    region_masks = {r: [pos_region(lbl) == r for lbl in POS_LABELS] for r in regions}

    print("\n--- Regional DPF (d={:.0f}mm) ---".format(D_PROBE))
    for r in regions:
        m = dpf_all[region_masks[r]]
        print("{:>8s}: n={}, mean={:.2f}, std={:.2f}, per wl={}".format(
            r, m.shape[0], m.mean(), m.std(), np.round(m.mean(axis=0), 2)))

    # --- ANOVA ---
    print("\n--- One-way ANOVA: regional DPF differences ---")
    for wi, wl in enumerate(TD_WAVELENGTHS):
        groups = [dpf_all[region_masks[r], wi] for r in regions]
        F, p = stats.f_oneway(*groups)
        print("  {} nm: F={:.3f}, p={:.4f} {}".format(wl, F, p, "*" if p < 0.05 else ""))

    # --- Scholkmann comparison ---
    print("\n--- Scholkmann (2013) reference DPF (age=25) ---")
    for wl in TD_WAVELENGTHS:
        print("  {} nm: {:.2f}".format(wl, scholkmann_dpf(wl)))

    # === PLOTS ===
    # 1. Example DTOFs
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    fig.suptitle("Example DTOFs — Position 1 (averaged over 5 reps)")
    for wi, ax in enumerate(axes.flat):
        avg_m = meas[1][wi * TD_N_REP:(wi + 1) * TD_N_REP].mean(axis=0)
        avg_i = irf[wi * TD_N_REP:(wi + 1) * TD_N_REP].mean(axis=0)
        ax.semilogy(t_bins / 1e3, avg_m, "b-", lw=1.2, label="Tissue DTOF")
        ax.semilogy(t_bins / 1e3, avg_i, "r--", lw=1.0, label="IRF")
        ax.set_title("{} nm".format(TD_WAVELENGTHS[wi]))
        ax.set_ylabel("Counts")
        ax.set_xlim(2.5, 6.5)
        ax.legend(fontsize=9)
    axes[1, 0].set_xlabel("Time (ns)")
    axes[1, 1].set_xlabel("Time (ns)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "td_dtof_example.png"), dpi=150)
    plt.close()

    # 2. DPF by region bar chart
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors_r = {"Left": "steelblue", "Central": "seagreen", "Right": "tomato"}
    x_pos = np.arange(TD_N_WL)
    width = 0.25
    for ri, r in enumerate(regions):
        m = dpf_all[region_masks[r]]
        means = m.mean(axis=0)
        sems = m.std(axis=0) / np.sqrt(m.shape[0])
        axes[0].bar(x_pos + ri * width, means, width, yerr=sems, capsize=4,
                    label=r, color=colors_r[r], alpha=0.85)
    axes[0].set_xticks(x_pos + width)
    axes[0].set_xticklabels(["{} nm".format(w) for w in TD_WAVELENGTHS])
    axes[0].set_ylabel("DPF (d={:.0f}mm)".format(D_PROBE))
    axes[0].set_title("DPF by region and wavelength")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.3)

    wl_range = np.arange(650, 870, 1)
    dpf_sch = np.array([scholkmann_dpf(w) for w in wl_range])
    axes[1].plot(wl_range, dpf_sch, "k-", lw=2, label="Scholkmann (age=25)")
    for r in regions:
        m = dpf_all[region_masks[r]]
        axes[1].scatter(TD_WAVELENGTHS, m.mean(axis=0), s=80, zorder=5,
                        color=colors_r[r], edgecolor="k", label=r)
    axes[1].set_xlabel("Wavelength (nm)")
    axes[1].set_ylabel("DPF")
    axes[1].set_title("Measured DPF vs. Scholkmann (2013)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "td_dpf_regions.png"), dpi=150)
    plt.close()

    # 3. Sensitivity to probe distance
    fig, ax = plt.subplots(figsize=(10, 5))
    d_range = np.arange(10, 35, 1)
    for wi, wl in enumerate(TD_WAVELENGTHS):
        mean_dpf = [V_TISSUE * delta_t[:, wi].mean() / dd for dd in d_range]
        ax.plot(d_range, mean_dpf, "o-", ms=3, label="{} nm".format(wl))
    ax.axhspan(5.0, 6.5, alpha=0.12, color="green", label="Literature range")
    ax.set_xlabel("Assumed TD probe distance (mm)")
    ax.set_ylabel("Mean DPF")
    ax.set_title("DPF sensitivity to probe S-D distance")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "td_dpf_sensitivity.png"), dpi=150)
    plt.close()

    # Save numerical results
    np.savez(os.path.join(OUT_DIR, "td_results.npz"),
             delta_t=delta_t, dpf_all=dpf_all, irf_mt=irf_mt,
             pos_labels=POS_LABELS, d_probe=D_PROBE)

    print("\nPlots saved to", OUT_DIR)


if __name__ == "__main__":
    main()
