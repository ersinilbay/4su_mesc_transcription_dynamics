from __future__ import annotations

import numpy as np
import scanpy as sc

from config import STATE_COL, STATE_MAP


def detect_state_column(ad_um: sc.AnnData) -> str | None:
    """
    Try to detect which obs column in the UMAP AnnData contains cell-state labels.
    """
    if STATE_COL and STATE_COL in ad_um.obs:
        return STATE_COL

    candidates_priority = [
        "state",
        "State",
        "umap_state",
        "UMAP_state",
        "cell_state",
        "CellState",
        "annotation",
        "annotations",
        "label",
        "labels",
        "celltype",
        "CellType",
        "cell_type",
        "leiden",
        "louvain",
        "cluster",
        "clusters",
    ]
    for c in candidates_priority:
        if c in ad_um.obs:
            return c

    hits = []
    for c in ad_um.obs.columns:
        cname = str(c).lower()
        if ("state" in cname) or ("annot" in cname) or ("cluster" in cname) or ("celltype" in cname):
            hits.append(c)

    return hits[0] if hits else None


def canonical_state_code(val):
    """
    Convert pretty state names to short codes where possible.
    Example:
        'Pluripotent' -> 'pluri'
        'Intermediate' -> 'inter'
        '2-cell like' -> '2C'
    """
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None

    if val in STATE_MAP:
        return STATE_MAP[val]

    s = str(val)

    # already short code?
    for pretty, code in STATE_MAP.items():
        if s == code:
            return code

    # normalized text match
    s_norm = s.strip().lower().replace(" ", "").replace("-", "")
    for pretty, code in STATE_MAP.items():
        p_norm = pretty.strip().lower().replace(" ", "").replace("-", "")
        if s_norm == p_norm:
            return code

    return s


def sync_state(ad_qc: sc.AnnData, ad_um: sc.AnnData) -> None:
    """
    Copy the state labels from the UMAP object into the QC object
    based on shared cell barcodes.
    """
    col = detect_state_column(ad_um)
    if col is None:
        print("[WARN] Could not detect a state column in UM.")
        preview = [(c, ", ".join(ad_um.obs[c].astype(str).unique()[:5])) for c in ad_um.obs.columns[:20]]
        import pandas as pd
        print(pd.DataFrame(preview, columns=["obs column", "sample values"]).to_string(index=False))
        return

    shared = ad_qc.obs_names.intersection(ad_um.obs_names)
    if len(shared) == 0:
        print("[WARN] No shared cell barcodes between QC and UM.")
        return

    ad_qc.obs.loc[shared, "state"] = ad_um.obs.loc[shared, col].values
    ad_qc.obs["state"] = ad_qc.obs["state"].map(canonical_state_code)

    vc = ad_qc.obs["state"].value_counts(dropna=False)
    print(f"[INFO] Using UM obs['{col}'] as state; normalized counts:\n{vc}")


def mask_for_state(adata: sc.AnnData, state_code: str) -> np.ndarray | None:
    """
    Return boolean mask for cells belonging to one state.
    """
    if "state" not in adata.obs:
        print("[WARN] obs['state'] not present; using ALL cells.")
        return None

    m = (adata.obs["state"] == state_code).values
    if m.sum() == 0:
        print(f"[WARN] No cells in state '{state_code}'.")
        return None

    return m