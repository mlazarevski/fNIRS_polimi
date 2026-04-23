"""CW-NIRS: Full analysis pipeline — quality, preprocessing, MBLL, block averaging, statistics."""
import numpy as np
import scipy.io
import scipy.signal as sig
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py, json, os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import *
from td_dpf import scholkmann_dpf


# ═══════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════
def load_cw_data():
    nirs  = scipy.io.loadmat(NIRS_FILE, simplify_cells=True)
    probe = scipy.io.loadmat(PROBE_FILE, simplify_cells=True)
    cal   = json.load(open(CAL_FILE))

    d  = nirs["d"]
    t  = nirs["t"].flatten()
    SD = nirs["SD"]
    s  = nirs["s"]

    src_pos = np.array(SD["SrcPos"])
    det_pos = np.array(SD["DetPos"])
    wavelengths = np.array(SD["Lambda"])
    fs = 1.0 / np.mean(np.diff(t))

    index_c = np.array(probe["probeInfo"]["probes"]["index_c"]).astype(int) - 1
    n_pairs = d.shape[1] // 2

    stim_info = {}
    with h5py.File(SNIRF_FILE, "r") as f:
        for k in ["stim1", "stim2", "stim3", "stim4"]:
            if "nirs/" + k in f:
                name = f["nirs/" + k + "/name"][()].decode()
                data = f["nirs/" + k + "/data"][:]
                stim_info[name] = data

    dists_mm = np.array([np.linalg.norm(src_pos[si] - det_pos[di])
                         for si, di in index_c])
    is_short = dists_mm < 15

    return {
        "d": d, "t": t, "s": s, "fs": fs,
        "src_pos": src_pos, "det_pos": det_pos,
        "wavelengths": wavelengths,
        "index_c": index_c, "n_pairs": n_pairs,
        "stim_info": stim_info,
        "dists_mm": dists_mm, "is_short": is_short,
        "cal": cal,
    }


# ═══════════════════════════════════════════════════════════
# SIGNAL QUALITY
# ═══════════════════════════════════════════════════════════
def assess_quality(data):
    d, t, cal = data["d"], data["t"], data["cal"]
    n_pairs, is_short = data["n_pairs"], data["is_short"]
    dists_mm, fs = data["dists_mm"], data["fs"]

    amps = np.array(cal["amplitudes"])
    cv = np.array(cal["coefficient_of_variation"])
    qlev = np.array([q["level"] for q in cal["signal_quality"]])

    bl_mask = t < 35.0
    snr_wl1 = d[bl_mask, :n_pairs].mean(0) / d[bl_mask, :n_pairs].std(0)
    snr_wl2 = d[bl_mask, n_pairs:].mean(0) / d[bl_mask, n_pairs:].std(0)

    print("\n" + "=" * 70)
    print("CW-NIRS SIGNAL QUALITY")
    print("=" * 70)
    print("All channels quality level 3:", all(q == 3 for q in qlev))
    print("Mean amplitude: wl1={:.4f}, wl2={:.4f}".format(amps[:, 0].mean(), amps[:, 1].mean()))
    print("Mean CV       : wl1={:.4f}, wl2={:.4f}".format(cv[:, 0].mean(), cv[:, 1].mean()))
    print("Long channels : {}".format((~is_short).sum()))
    print("Short channels: {} (D8-D15)".format(is_short.sum()))
    print("Min baseline SNR (long): wl1={:.1f}, wl2={:.1f}".format(
        snr_wl1[~is_short].min(), snr_wl2[~is_short].min()))

    # Flag low-quality channels
    low_snr = (snr_wl1 < 10) | (snr_wl2 < 10)
    low_amp = (amps[:, 0] < 0.01) | (amps[:, 1] < 0.01)
    flagged = np.where(low_snr | low_amp)[0]
    if len(flagged) > 0:
        print("Flagged channels (low SNR or amplitude):", flagged + 1)
        for ch in flagged:
            print("  Ch {}: SNR=({:.1f},{:.1f}), amp=({:.4f},{:.4f}), dist={:.1f}mm {}".format(
                ch + 1, snr_wl1[ch], snr_wl2[ch], amps[ch, 0], amps[ch, 1],
                dists_mm[ch], "[SHORT]" if is_short[ch] else ""))
    else:
        print("No channels flagged for low quality.")

    # --- Plot ---
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    axes[0].bar(np.arange(n_pairs) - 0.15, amps[:, 0], 0.3, label="760 nm", color="C0")
    axes[0].bar(np.arange(n_pairs) + 0.15, amps[:, 1], 0.3, label="850 nm", color="C1")
    axes[0].set_xlabel("Channel pair"); axes[0].set_ylabel("Amplitude")
    axes[0].set_title("Calibration amplitudes"); axes[0].legend(fontsize=9)

    colors_ch = ["tomato" if sh else "steelblue" for sh in is_short]
    axes[1].bar(np.arange(n_pairs), dists_mm, color=colors_ch)
    axes[1].axhline(15, color="gray", ls="--", lw=0.8)
    axes[1].set_xlabel("Channel pair"); axes[1].set_ylabel("Distance (mm)")
    axes[1].set_title("S-D distances (red=short)")

    axes[2].bar(np.arange(n_pairs) - 0.15, snr_wl1, 0.3, label="760 nm", color="C0")
    axes[2].bar(np.arange(n_pairs) + 0.15, snr_wl2, 0.3, label="850 nm", color="C1")
    axes[2].set_xlabel("Channel pair"); axes[2].set_ylabel("SNR")
    axes[2].set_title("Baseline SNR"); axes[2].legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_quality.png"), dpi=150)
    plt.close()

    return {"amps": amps, "cv": cv, "qlev": qlev, "snr_wl1": snr_wl1, "snr_wl2": snr_wl2}


# ═══════════════════════════════════════════════════════════
# PREPROCESSING
# ═══════════════════════════════════════════════════════════
def preprocess(data):
    d, t, n_pairs, fs = data["d"], data["t"], data["n_pairs"], data["fs"]
    is_short, stim_info = data["is_short"], data["stim_info"]

    eps = 1e-12
    I0 = d.mean(axis=0, keepdims=True)
    OD_raw = -np.log((d + eps) / (I0 + eps))

    # Motion artifact detection
    dod_dt = np.diff(OD_raw[:, :n_pairs], axis=0)
    spike_thresh = 0.15
    spike_counts = (np.abs(dod_dt) > spike_thresh).sum(axis=0)

    print("\n" + "=" * 70)
    print("PREPROCESSING: ARTIFACT DETECTION")
    print("=" * 70)
    print("Spike threshold: |dOD/dt| > {}".format(spike_thresh))
    print("Spike counts per channel (first 22 long):")
    for ch in range(min(22, n_pairs)):
        tag = " [SHORT]" if is_short[ch] else ""
        if spike_counts[ch] > 0:
            print("  Ch {:2d}: {:3d} spikes{}".format(ch + 1, spike_counts[ch], tag))
    if spike_counts[~is_short].max() == 0:
        print("  No spikes detected in any long channel.")

    # Bandpass filter
    def bandpass(x, fs, lo=0.01, hi=0.20, order=3):
        nyq = 0.5 * fs
        b, a = sig.butter(order, [lo / nyq, hi / nyq], btype="band")
        return sig.filtfilt(b, a, x, axis=0)

    OD_filt = bandpass(OD_raw, fs)
    print("\nBandpass filter applied: 0.01-0.20 Hz, 3rd-order Butterworth")

    # --- Plots ---
    # Raw intensity heatmap
    fig, axes = plt.subplots(2, 1, figsize=(15, 7), sharex=True)
    fig.suptitle("Raw CW intensity — all channels")
    for ax, sl, title in [(axes[0], slice(0, n_pairs), "760 nm"),
                           (axes[1], slice(n_pairs, 2 * n_pairs), "850 nm")]:
        chunk = d[:, sl]
        ax.imshow(chunk.T, aspect="auto", extent=[t[0], t[-1], 1, n_pairs],
                  cmap="viridis", vmin=np.percentile(chunk, 2), vmax=np.percentile(chunk, 98))
        ax.set_title(title); ax.set_ylabel("Channel pair")
        for _, cdata in stim_info.items():
            for onset, dur, *_ in cdata:
                ax.axvspan(onset, onset + dur, color="gold", alpha=0.12)
    axes[1].set_xlabel("Time (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_raw_heatmap.png"), dpi=150)
    plt.close()

    # OD before/after filtering
    ch_ex = 5
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    fig.suptitle("OD before and after bandpass — Channel {}".format(ch_ex + 1))
    axes[0].plot(t, OD_raw[:, ch_ex], lw=0.5, color="gray")
    axes[0].set_ylabel("Raw OD")
    axes[1].plot(t, OD_filt[:, ch_ex], lw=0.8, color="C0")
    for _, cdata in stim_info.items():
        for onset, dur, *_ in cdata:
            axes[1].axvspan(onset, onset + dur, color="gold", alpha=0.10)
    axes[1].set_ylabel("Filtered OD"); axes[1].set_xlabel("Time (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_od_filter.png"), dpi=150)
    plt.close()

    # Filtered OD heatmap
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle("Filtered OD — all channels")
    for ax, sl, title in [(axes[0], slice(0, n_pairs), "760 nm"),
                           (axes[1], slice(n_pairs, 2 * n_pairs), "850 nm")]:
        im = ax.imshow(OD_filt[:, sl].T, aspect="auto",
                       extent=[t[0], t[-1], 1, n_pairs], cmap="RdBu_r",
                       vmin=-0.05, vmax=0.05)
        ax.set_title(title); ax.set_ylabel("Channel pair"); ax.set_xlabel("Time (s)")
        plt.colorbar(im, ax=ax, label="OD")
        for _, cdata in stim_info.items():
            for onset, dur, *_ in cdata:
                ax.axvspan(onset, onset + dur, color="gold", alpha=0.12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_od_heatmap.png"), dpi=150)
    plt.close()

    return OD_filt


# ═══════════════════════════════════════════════════════════
# MBLL
# ═══════════════════════════════════════════════════════════
def compute_hemoglobin(data, OD_filt):
    t, n_pairs = data["t"], data["n_pairs"]
    dists_mm, is_short = data["dists_mm"], data["is_short"]
    index_c, stim_info = data["index_c"], data["stim_info"]
    src_pos, det_pos = data["src_pos"], data["det_pos"]

    DPF_760 = 6.0
    DPF_850 = 6.0

    print("\n" + "=" * 70)
    print("MBLL: HEMOGLOBIN CONCENTRATION CHANGES")
    print("=" * 70)
    print("DPF (fixed, matching MATLAB pipeline): 760nm={:.2f}, 850nm={:.2f}".format(DPF_760, DPF_850))

    HbO = np.zeros((len(t), n_pairs))
    HbR = np.zeros((len(t), n_pairs))

    for ch in range(n_pairs):
        dOD = np.column_stack([OD_filt[:, ch], OD_filt[:, ch + n_pairs]])
        rho = dists_mm[ch] / 10.0
        A = np.array([
            [EXT_COEFFS[0, 0] * rho * DPF_760, EXT_COEFFS[0, 1] * rho * DPF_760],
            [EXT_COEFFS[1, 0] * rho * DPF_850, EXT_COEFFS[1, 1] * rho * DPF_850],
        ])
        conc = dOD @ np.linalg.pinv(A).T
        HbO[:, ch] = conc[:, 0]
        HbR[:, ch] = conc[:, 1]

    HbT = HbO + HbR
    print("HbO/HbR/HbT computed: shape =", HbO.shape)

    # Plot example long channels
    ex_chs = [i for i in range(n_pairs) if not is_short[i]][:6]
    fig, axes = plt.subplots(len(ex_chs), 1, figsize=(15, 2.5 * len(ex_chs)), sharex=True)
    fig.suptitle("Hemoglobin concentration — example long channels")
    for i, ch in enumerate(ex_chs):
        axes[i].plot(t, HbO[:, ch] * 1e6, color="crimson", lw=0.8, label="HbO")
        axes[i].plot(t, HbR[:, ch] * 1e6, color="royalblue", lw=0.8, label="HbR")
        for _, cdata in stim_info.items():
            for onset, dur, *_ in cdata:
                axes[i].axvspan(onset, onset + dur, color="gold", alpha=0.08)
        axes[i].set_ylabel(r"$\Delta$c ($\mu$M)")
        si, di = index_c[ch]
        axes[i].set_title("Ch {} | S{}-D{} ({:.1f}mm)".format(
            ch + 1, si + 1, di + 1, dists_mm[ch]), fontsize=10)
        if i == 0:
            axes[i].legend(fontsize=9)
    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_hemoglobin.png"), dpi=150)
    plt.close()

    return HbO, HbR, HbT


# ═══════════════════════════════════════════════════════════
# BLOCK AVERAGING & STATISTICS
# ═══════════════════════════════════════════════════════════
def block_average_and_stats(data, HbO, HbR):
    t, fs = data["t"], data["fs"]
    n_pairs, is_short = data["n_pairs"], data["is_short"]
    index_c, stim_info = data["index_c"], data["stim_info"]
    src_pos, det_pos = data["src_pos"], data["det_pos"]

    PRE, POST = 5.0, 25.0
    pre_samp  = int(np.round(PRE * fs))
    post_samp = int(np.round(POST * fs))
    epoch_len = pre_samp + post_samp
    t_epoch   = np.linspace(-PRE, POST, epoch_len)

    # ROI assignment
    src_x = src_pos[:, 0]
    ch_roi = []
    for si, _ in index_c:
        if src_x[si] < -15:
            ch_roi.append("Left")
        elif src_x[si] > 15:
            ch_roi.append("Right")
        else:
            ch_roi.append("Central")
    ch_roi = np.array(ch_roi)
    ch_long = ~is_short

    print("\n" + "=" * 70)
    print("BLOCK AVERAGING & STATISTICAL ANALYSIS")
    print("=" * 70)
    print("Epoch window: [{}, {}] s around onset".format(-PRE, POST))
    print("Baseline correction: [-{}, 0] s".format(PRE))

    for roi in ["Left", "Central", "Right"]:
        chs = [i + 1 for i in range(n_pairs) if ch_roi[i] == roi and ch_long[i]]
        print("  {} ROI: {} long channels {}".format(roi, len(chs), chs))

    # Extract epochs
    cond_list = list(COND_MAP.values())
    epochs_hbo, epochs_hbr = {}, {}

    for code, cdata in stim_info.items():
        cond = COND_MAP[code]
        onsets = cdata[:, 0]
        n_trials = len(onsets)
        ep_o = np.full((n_trials, epoch_len, n_pairs), np.nan)
        ep_r = np.full((n_trials, epoch_len, n_pairs), np.nan)

        for ti, onset in enumerate(onsets):
            idx0 = np.argmin(np.abs(t - onset))
            i0, i1 = idx0 - pre_samp, idx0 + post_samp
            if i0 < 0 or i1 > len(t):
                continue
            seg_o = HbO[i0:i1, :]
            seg_r = HbR[i0:i1, :]
            bl_o = seg_o[:pre_samp].mean(axis=0, keepdims=True)
            bl_r = seg_r[:pre_samp].mean(axis=0, keepdims=True)
            ep_o[ti] = seg_o - bl_o
            ep_r[ti] = seg_r - bl_r

        epochs_hbo[cond] = ep_o
        epochs_hbr[cond] = ep_r
        print("  {}: {} trials".format(cond, n_trials))

    # --- Plot 1: Grand average per condition ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle("Block-averaged HRF — Grand average (all long channels)")
    for ci, (cond, ax) in enumerate(zip(cond_list, axes.flat)):
        ep_o = epochs_hbo[cond][:, :, ch_long]
        ep_r = epochs_hbr[cond][:, :, ch_long]
        m_o = np.nanmean(ep_o, axis=(0, 2)) * 1e6
        m_r = np.nanmean(ep_r, axis=(0, 2)) * 1e6
        se_o = np.nanstd(np.nanmean(ep_o, axis=2), axis=0) / np.sqrt(ep_o.shape[0]) * 1e6
        se_r = np.nanstd(np.nanmean(ep_r, axis=2), axis=0) / np.sqrt(ep_r.shape[0]) * 1e6
        ax.plot(t_epoch, m_o, "r-", lw=2, label="HbO")
        ax.fill_between(t_epoch, m_o - se_o, m_o + se_o, color="r", alpha=0.15)
        ax.plot(t_epoch, m_r, "b-", lw=2, label="HbR")
        ax.fill_between(t_epoch, m_r - se_r, m_r + se_r, color="b", alpha=0.15)
        ax.axvline(0, color="k", ls="--", lw=0.8)
        ax.axvline(10, color="k", ls=":", lw=0.8)
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_title(cond); ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$\Delta$c ($\mu$M)"); ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_blockavg_grand.png"), dpi=150)
    plt.close()

    # --- Plot 2: ROI x Condition ---
    rois = ["Left", "Central", "Right"]
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharex=True)
    fig.suptitle("Block-averaged HbO by ROI and condition")
    for ri, roi in enumerate(rois):
        rmask = np.array([(ch_roi[ch] == roi) and ch_long[ch] for ch in range(n_pairs)])
        for ci, cond in enumerate(cond_list):
            ax = axes[ri, ci]
            if rmask.sum() == 0:
                ax.text(0.5, 0.5, "No channels", transform=ax.transAxes, ha="center")
                continue
            m_o = np.nanmean(epochs_hbo[cond][:, :, rmask], axis=(0, 2)) * 1e6
            m_r = np.nanmean(epochs_hbr[cond][:, :, rmask], axis=(0, 2)) * 1e6
            ax.plot(t_epoch, m_o, "r-", lw=1.5)
            ax.plot(t_epoch, m_r, "b-", lw=1.5)
            ax.axvline(0, color="k", ls="--", lw=0.6)
            ax.axvline(10, color="k", ls=":", lw=0.6)
            ax.axhline(0, color="gray", lw=0.4)
            ax.grid(True, alpha=0.15)
            if ri == 0:
                ax.set_title(cond, fontsize=10)
            if ci == 0:
                ax.set_ylabel(roi + r"  $\Delta$c ($\mu$M)")
            if ri == 2:
                ax.set_xlabel("Time (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_blockavg_roi.png"), dpi=150)
    plt.close()

    # --- Plot 3: All 4 conditions overlaid per ROI ---
    cond_colors = {"Peripheral Left": "C0", "Peripheral Right": "C1",
                   "Central Left": "C2", "Central Right": "C3"}
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    fig.suptitle("HbO block averages: all conditions overlaid per ROI")
    for ri, roi in enumerate(rois):
        rmask = np.array([(ch_roi[ch] == roi) and ch_long[ch] for ch in range(n_pairs)])
        ax = axes[ri]
        for cond in cond_list:
            m_o = np.nanmean(epochs_hbo[cond][:, :, rmask], axis=(0, 2)) * 1e6
            ax.plot(t_epoch, m_o, lw=2, color=cond_colors[cond], label=cond)
        ax.axvline(0, color="k", ls="--", lw=0.6)
        ax.axvline(10, color="k", ls=":", lw=0.6)
        ax.axhline(0, color="gray", lw=0.4)
        ax.set_title(roi + " ROI")
        ax.set_xlabel("Time (s)")
        if ri == 0:
            ax.set_ylabel(r"$\Delta$HbO ($\mu$M)")
            ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_blockavg_overlay.png"), dpi=150)
    plt.close()

    # --- Quantitative: activation window ---
    ACT_WIN = (4, 10)
    act_mask = (t_epoch >= ACT_WIN[0]) & (t_epoch <= ACT_WIN[1])

    results = []
    for roi in rois:
        rmask = np.array([(ch_roi[ch] == roi) and ch_long[ch] for ch in range(n_pairs)])
        for cond in cond_list:
            ep = epochs_hbo[cond][:, :, rmask]
            trial_means = np.nanmean(ep[:, act_mask, :], axis=(1, 2)) * 1e6
            for tm in trial_means:
                results.append({"ROI": roi, "Condition": cond, "HbO_uM": tm})

    df = pd.DataFrame(results)

    # --- Plot 4: Bar chart ---
    colors_r = {"Left": "steelblue", "Central": "seagreen", "Right": "tomato"}
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(cond_list))
    width = 0.22
    for ri, roi in enumerate(rois):
        sub = df[df["ROI"] == roi]
        means = [sub[sub["Condition"] == c]["HbO_uM"].mean() for c in cond_list]
        sems = [sub[sub["Condition"] == c]["HbO_uM"].std() /
                np.sqrt(sub[sub["Condition"] == c].shape[0]) for c in cond_list]
        ax.bar(x + ri * width, means, width, yerr=sems, capsize=4,
               label=roi, color=colors_r[roi], alpha=0.85)
    ax.set_xticks(x + width)
    ax.set_xticklabels(cond_list, rotation=15)
    ax.set_ylabel(r"Mean $\Delta$HbO ($\mu$M) [4-10 s window]")
    ax.set_title("Activation-window HbO by condition and ROI")
    ax.legend(); ax.axhline(0, color="gray", lw=0.5)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_activation_bars.png"), dpi=150)
    plt.close()

    # --- Plot 5: Probe topography ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("Activation topography (mean HbO, 4-10 s)")
    for ci, (cond, ax) in enumerate(zip(cond_list, axes.flat)):
        ep = epochs_hbo[cond]
        metric = np.nanmean(ep[:, act_mask, :], axis=(0, 1)) * 1e6

        ax.scatter(src_pos[:, 0], src_pos[:, 1], marker="^", s=120, color="darkred", zorder=5)
        ax.scatter(det_pos[:, 0], det_pos[:, 1], marker="s", s=60, color="navy", zorder=5)
        for i in range(src_pos.shape[0]):
            ax.text(src_pos[i, 0] + 1.5, src_pos[i, 1] + 1.5,
                    "S{}".format(i + 1), fontsize=7, color="darkred")
        for i in range(det_pos.shape[0]):
            ax.text(det_pos[i, 0] + 1.5, det_pos[i, 1] + 1.5,
                    "D{}".format(i + 1), fontsize=6, color="navy")

        midx, midy, vals = [], [], []
        for ch in range(n_pairs):
            if is_short[ch]:
                continue
            si, di = index_c[ch]
            ax.plot([src_pos[si, 0], det_pos[di, 0]],
                    [src_pos[si, 1], det_pos[di, 1]], "k-", alpha=0.15, lw=0.5)
            midx.append((src_pos[si, 0] + det_pos[di, 0]) / 2)
            midy.append((src_pos[si, 1] + det_pos[di, 1]) / 2)
            vals.append(metric[ch])

        midx, midy, vals = np.array(midx), np.array(midy), np.array(vals)
        vmax = max(abs(np.nanmin(vals)), abs(np.nanmax(vals)), 1e-10)
        sc = ax.scatter(midx, midy, c=vals, s=150, cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, edgecolor="k", lw=0.5, zorder=10)
        plt.colorbar(sc, ax=ax, label=r"$\Delta$HbO ($\mu$M)", shrink=0.8)
        ax.set_title(cond); ax.set_xlabel("X (mm)"); ax.set_ylabel("Y (mm)")
        ax.axis("equal")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cw_topography.png"), dpi=150)
    plt.close()

    # --- Statistics ---
    print("\n--- STATISTICAL RESULTS ---")
    pivot = df.pivot_table(values="HbO_uM", index="ROI", columns="Condition", aggfunc="mean")
    print("\nMean HbO (uM) in activation window:")
    print(pivot.round(4).to_string())

    print("\nOne-way ANOVA (Condition effect) per ROI:")
    for roi in rois:
        sub = df[df["ROI"] == roi]
        groups = [sub[sub["Condition"] == c]["HbO_uM"].values for c in cond_list]
        F, p = stats.f_oneway(*groups)
        sig_str = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print("  {}: F={:.3f}, p={:.4f} {}".format(roi, F, p, sig_str))

    print("\nPairwise t-tests:")
    comparisons = [
        ("Peripheral Left", "Central Left", "Periph vs Central (left stim)"),
        ("Peripheral Right", "Central Right", "Periph vs Central (right stim)"),
        ("Peripheral Left", "Peripheral Right", "Left vs Right (peripheral)"),
        ("Central Left", "Central Right", "Left vs Right (central)"),
    ]
    for c1, c2, desc in comparisons:
        v1 = df[df["Condition"] == c1]["HbO_uM"].values
        v2 = df[df["Condition"] == c2]["HbO_uM"].values
        t_stat, p_val = stats.ttest_ind(v1, v2)
        sig_str = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print("  {} vs {}: t={:.3f}, p={:.4f} {}  ({})".format(c1, c2, t_stat, p_val, sig_str, desc))

    # Save
    np.savez(os.path.join(OUT_DIR, "cw_results.npz"),
             t_epoch=t_epoch, HbO=HbO, HbR=HbR)
    df.to_csv(os.path.join(OUT_DIR, "cw_activation_stats.csv"), index=False)

    print("\nPlots and data saved to", OUT_DIR)


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════
def main():
    data = load_cw_data()
    assess_quality(data)
    OD_filt = preprocess(data)
    HbO, HbR, HbT = compute_hemoglobin(data, OD_filt)
    block_average_and_stats(data, HbO, HbR)


if __name__ == "__main__":
    main()
