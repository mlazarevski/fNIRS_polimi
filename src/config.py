"""Shared paths and constants for the Flickering Wedges analysis."""
import os
import numpy as np

BASE    = "/Users/leonardo/GitHub/Brain/nirs/fNIRS_polimi/Flickering Wedges"
CW_DIR  = os.path.join(BASE, "CW - Wedges")
ROI_DIR = os.path.join(BASE, "ROI")
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "results")
os.makedirs(OUT_DIR, exist_ok=True)

NIRS_FILE  = os.path.join(CW_DIR, "2025-04-11_009.nirs")
SNIRF_FILE = os.path.join(CW_DIR, "2025-04-11_009.snirf")
PROBE_FILE = os.path.join(CW_DIR, "2025-04-11_009_probeInfo.mat")
CAL_FILE   = os.path.join(CW_DIR, "2025-04-11_009_calibration.json")
TRI_FILE   = os.path.join(CW_DIR, "2025-04-11_009.tri")

DAT_FILES = sorted([os.path.join(BASE, f) for f in os.listdir(BASE) if f.endswith(".DAT")])
MEAS_DATS = [f for f in DAT_FILES if "VmLC" in os.path.basename(f)]
IRF_DAT   = [f for f in DAT_FILES if "VsLC" in os.path.basename(f)][0]

# TD-NIRS constants
TD_HEADER    = 764
TD_SUBHDR    = 204
TD_NBINS     = 4096
TD_PS_BIN    = 6.0
TD_N_WL      = 4
TD_N_REP     = 5
TD_NREPS     = TD_N_WL * TD_N_REP
TD_WAVELENGTHS = np.array([670, 730, 780, 830])

# CW wavelengths
CW_WAVELENGTHS = np.array([760, 850])

# Stimulus condition mapping (SNIRF label → human name)
COND_MAP = {"1": "Peripheral Left", "2": "Peripheral Right",
            "4": "Central Left",    "8": "Central Right"}

# Optode position labels for the 22 VmLC files
# s5 is co-located with Iz and was used for IRF
POS_LABELS = ["S1", "S2", "S3", "S4", "S6", "S7", "S8",
              "D1", "D2", "D3", "D4", "D5", "D6", "D7",
              "D8", "D9", "D10", "D11", "D12", "D13", "D14", "D15"]

# Optode x-coordinates from digpts.txt (for region assignment)
POS_X = {
    "S1": -56.6, "S2": -54.0, "S3": -30.2, "S4": -0.2,
    "S6": 30.6,  "S7": 57.0,  "S8": 54.7,
    "D1": -56.1, "D2": -38.8, "D3": -28.0, "D4": -0.02,
    "D5": 39.0,  "D6": 56.7,  "D7": 28.7,
    "D8": -62.6, "D9": -58.8, "D10": -37.7, "D11": -7.5,
    "D12": -7.0, "D13": 22.5, "D14": 51.7,  "D15": 48.6,
}

def pos_region(label):
    x = POS_X.get(label, 0)
    if x < -15:
        return "Left"
    elif x > 15:
        return "Right"
    return "Central"

# Speed of light in tissue
N_TISSUE = 1.4
V_TISSUE = 0.3 / N_TISSUE  # mm/ps

# Extinction coefficients (1/cm/M) from Prahl/OMLC
# Rows: [760nm, 850nm]; Columns: [HbO, HbR]
EXT_COEFFS = np.array([
    [1486.5865, 3843.707],
    [2526.391,  1798.643],
])
