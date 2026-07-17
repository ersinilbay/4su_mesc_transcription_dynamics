from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt

# ======================== USER PATHS =========================
REP1_QC = Path("/Users/ersinilbay/PycharmProjects/Master-internship/rep1_properfiltering/adata_qc_raw_with_kinetics_and_states_fix.h5ad")
REP1_UM = Path("/Users/ersinilbay/PycharmProjects/Master-internship/rep1_properfiltering/adata_umap_with_states_fix.h5ad")
REP2_QC = Path("/Users/ersinilbay/PycharmProjects/Master-internship/rep2_fix/adata_qc_raw_with_kinetics_and_states_rep2_fix.h5ad")
REP2_UM = Path("/Users/ersinilbay/PycharmProjects/Master-internship/rep2_fix/adata_umap_with_states_rep2_fix.h5ad")
OUT_DIR = Path("/Users/ersinilbay/PycharmProjects/Master-internship/fanostuff")

# ======================== FIGURE STYLE =======================
FIGSIZE_SCAT = (4.8, 4.8)
FIGSIZE_RESID_SMALL = (5.2, 5.2)
ANNOTATED_SCALE = 1.15
FIGSIZE_RESID_ANN = (
    FIGSIZE_RESID_SMALL[0] * ANNOTATED_SCALE,
    FIGSIZE_RESID_SMALL[1] * ANNOTATED_SCALE,
)

BLUE = "#1f77b4"
RIDGE = "0.25"
RED = "crimson"

STATE_MAP = {"Pluripotent": "pluri", "Intermediate": "inter", "2-cell like": "2C"}
STATE_COL: str | None = None

RESIDUAL_CFG = dict(
    n_label=35,
    min_mean=0.1,
    max_mean=6.0,
    select_rule="zscore",
    z_thr=1.2,
    min_fano=2.0,
    z_thr_display=0.9,
    min_fano_display=1.6,
    top_frac_global=0.10,
    min_fold=1.5,
    exclude_pref=("mt-", "Rpl", "Rps", "Mrps", "Gm"),
    exclude_suf=("Rik", "-ps"),
    whitelist=(),
    min_detect_frac=0.10,
    winsor_q=0.99,
)

RESIDUAL_STATE_OVERRIDES = {
    "pluri": {
        "n_label": 40,
        "z_thr_display": 4.2,
        "min_fano_display": 1.80,
        "z_thr": 4.5,
        "min_fano": 1.80,
        "whitelist": ("Rpl23", "mt-Nd5"),
        "min_mean": 0.07,
        "max_mean": 6.0,
    },
    "2C": {
        "n_label": 40,
        "z_thr_display": 1.1,
        "min_fano_display": 2.0,
    },
}

ENRICH_STRICT = False
PREFERRED_ENRICH_LIBRARIES = ("GO_Biological_Process_2021",)

GO_MIN_OVERLAP = 3
GO_MAX_TERMS = 10
GO_FDR_MAX = 0.05

GO_FIG_WIDTH_IN = 6.0
GO_ROW_HEIGHT_IN = 0.34
GO_LEFT_MARGIN = 0.46
GO_BAR_HEIGHT = 0.65
GO_TERM_MAXCHARS = 55

PRINT_ENRICHR_CATALOG = True


def apply_plot_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "figure.constrained_layout.use": True,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.0,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.9,
            "ytick.minor.width": 0.9,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "ytick.minor.size": 2.0,
            "legend.frameon": True,
            "legend.framealpha": 0.92,
            "legend.facecolor": "white",
            "legend.edgecolor": "0.3",
        }
    )