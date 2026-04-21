from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def preferred_layer(adata: sc.AnnData, names: list[str]) -> np.ndarray | None:
    """
    Return the first matching layer found in AnnData, converted to dense if needed.
    """
    for nm in names:
        if nm in adata.layers:
            X = adata.layers[nm]
            return X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return None


def fraction_new_from_layers(adata: sc.AnnData) -> np.ndarray | None:
    """
    Compute new / total fraction from layer pairs if explicit new-fraction layer is absent.
    """
    new = preferred_layer(adata, ["new", "C", "new_counts", "newrna"])
    tot = preferred_layer(adata, ["total", "T", "total_counts", "oldrna_plus_newrna", "counts"])

    if new is None or tot is None:
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        frac = new / np.maximum(tot, 1e-12)
        frac[~np.isfinite(frac)] = 0.0

    return frac


def new_fraction_matrix(adata_qc: sc.AnnData) -> np.ndarray:
    """
    Return new-fraction / NTR matrix.
    """
    frac = preferred_layer(adata_qc, ["new_frac", "ntr", "NTR"])
    if frac is not None:
        return frac

    frac = fraction_new_from_layers(adata_qc)
    if frac is not None:
        return frac

    raise ValueError("Couldn't find NEW fraction. Provide 'new_frac' (or NEW/TOTAL layers).")


def total_matrix(adata_qc: sc.AnnData) -> np.ndarray:
    """
    Return total count matrix if available, otherwise fall back to X.
    """
    X = preferred_layer(adata_qc, ["total", "T", "total_counts", "counts"])
    if X is not None:
        return X

    X = adata_qc.X
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def safe_stats_over_cells(
    mat: np.ndarray,
    mask_cells: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-gene mean and variance across cells.
    """
    if mask_cells is not None:
        mat = mat[mask_cells, :]

    if mat.size == 0:
        n_genes = mat.shape[1] if mat.ndim == 2 else 0
        return np.full(n_genes, np.nan), np.full(n_genes, np.nan)

    mat = np.asarray(mat, dtype=np.float64)
    mean = np.nanmean(mat, axis=0)
    var = np.nanvar(mat, axis=0, ddof=1)

    return np.asarray(mean).ravel(), np.asarray(var).ravel()


def detect_fraction(TOT: np.ndarray, mask_cells: np.ndarray | None) -> np.ndarray:
    """
    Fraction of cells in which each gene is detected (>0).
    """
    X = TOT[mask_cells, :] if mask_cells is not None else TOT
    if hasattr(X, "toarray"):
        X = X.toarray()

    return (X > 0).mean(axis=0).ravel()


def gene_stats_from_total(adata: sc.AnnData) -> pd.DataFrame:
    """
    Per-gene mean, variance, and Fano factor from TOTAL counts.
    """
    X = total_matrix(adata)
    mean, var = safe_stats_over_cells(X, mask_cells=None)
    fano = var / np.maximum(mean, 1e-20)

    return pd.DataFrame(
        {"mean": mean, "variance": var, "fano": fano},
        index=adata.var_names,
    )


def list_offscale_fano_genes(
    adata: sc.AnnData,
    *,
    fano_min: float = 10.0,
    mean_min: float = 0.10,
    out_csv: Path | None = None,
    tag: str = "",
) -> list[str]:
    """
    List genes with very large Fano at non-tiny mean; optionally save CSV.
    """
    st = gene_stats_from_total(adata)
    keep = (st["fano"] >= fano_min) & (st["mean"] >= mean_min) & np.isfinite(st["fano"])
    df = st.loc[keep].sort_values(["fano", "mean"], ascending=[False, False])

    if out_csv is not None:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv)
        print(f"[INFO] {tag}: saved off-scale Fano genes -> {out_csv} (n={len(df)})")
    else:
        print(f"[INFO] {tag}: off-scale Fano genes n={len(df)}")

    return df.index.tolist()


def fit_expected_logfano(
    log_mean: np.ndarray,
    log_fano: np.ndarray,
    nbins: int = 40,
    smooth: int = 3,
    min_per_bin: int = 30,
):
    """
    Fit a smoothed expected log-Fano curve as a function of log-mean.
    """
    qs = np.linspace(0.0, 1.0, nbins + 1)
    edges = np.quantile(log_mean, qs)
    edges = np.unique(edges)

    if edges.size < 5:
        edges = np.linspace(log_mean.min(), log_mean.max(), 20)

    bins = np.digitize(log_mean, edges[1:-1], right=False)
    med_x, med_y = [], []

    for b in range(bins.min(), bins.max() + 1):
        m = bins == b
        if m.sum() >= min_per_bin:
            med_x.append(np.median(log_mean[m]))
            med_y.append(np.median(log_fano[m]))

    med_x = np.asarray(med_x)
    med_y = np.asarray(med_y)

    if med_x.size == 0:
        med_x = np.array([np.median(log_mean) - 1e-6, np.median(log_mean) + 1e-6])
        med_y = np.array([np.median(log_fano), np.median(log_fano)])

    if smooth > 1 and med_y.size >= smooth:
        k = np.ones(smooth) / smooth
        med_y = np.convolve(med_y, k, mode="same")

    order = np.argsort(med_x)
    med_x, med_y = med_x[order], med_y[order]

    def evaluate(logm: np.ndarray) -> np.ndarray:
        return np.interp(logm, med_x, med_y, left=med_y[0], right=med_y[-1])

    return med_x, med_y, evaluate


def expected_curve_eval(x: np.ndarray, f: np.ndarray, winsor_q: float):
    """
    Evaluate expected Fano curve across genes.
    """
    logx = np.log10(x)
    logf = np.log10(f)

    if 0.9 < winsor_q < 1.0:
        q = np.quantile(logf, winsor_q)
        logf_fit = np.clip(logf, None, q)
    else:
        logf_fit = logf

    _, _, eval_fn = fit_expected_logfano(
        logx,
        logf_fit,
        nbins=40,
        smooth=3,
        min_per_bin=30,
    )
    exp_f = np.power(10.0, eval_fn(logx))
    return eval_fn, exp_f


def residual_prepare_arrays(
    adata_qc: sc.AnnData,
    mask_cells,
    state_tag: str,
    residual_cfg: dict,
    overrides: dict,
):
    """
    Prepare filtered mean/Fano/name arrays for residual-Fano analysis.
    """
    cfg = residual_cfg.copy()
    ov = overrides.get(state_tag, {})
    for k, v in ov.items():
        cfg[k] = v

    TOT = total_matrix(adata_qc)
    mean, var = safe_stats_over_cells(TOT, mask_cells=mask_cells)
    det_all = detect_fraction(TOT, mask_cells)

    m = np.isfinite(mean) & np.isfinite(var) & (mean > 0)
    x = mean[m]
    f = var[m] / np.maximum(mean[m], 1e-20)
    det = det_all[m]
    names = np.array(adata_qc.var_names)[m].astype(str)

    keep_base = (x >= cfg["min_mean"]) & (x <= cfg["max_mean"])
    n_cells_state = mask_cells.sum() if isinstance(mask_cells, np.ndarray) else adata_qc.n_obs
    eff_min_det_frac = max(cfg.get("min_detect_frac", 0.0), 2.0 / max(n_cells_state, 1))
    keep_base &= det >= eff_min_det_frac

    keep_excl = np.ones_like(keep_base, dtype=bool)
    for pref in tuple(cfg.get("exclude_pref", ())):
        keep_excl &= ~np.char.startswith(names, pref)
    for suf in tuple(cfg.get("exclude_suf", ())):
        keep_excl &= ~np.char.endswith(names, suf)

    wl = set(map(str, cfg.get("whitelist", ())))
    if wl:
        keep_excl |= np.isin(names, list(wl))

    keep = keep_base & keep_excl
    x, f, names, det = x[keep], f[keep], names[keep], det[keep]

    return x, f, names, det, cfg


def residual_arrays_and_masks(
    adata_qc: sc.AnnData,
    mask_cells,
    state_tag: str,
    residual_cfg: dict,
    overrides: dict,
):
    """
    Compute expected-Fano residual arrays and strict/display masks.
    """
    x, f, names, det, cfg = residual_prepare_arrays(
        adata_qc,
        mask_cells,
        state_tag,
        residual_cfg,
        overrides,
    )

    _, exp_f = expected_curve_eval(x, f, cfg["winsor_q"])
    resid_ratio = f / np.maximum(exp_f, 1e-20)

    logf = np.log10(f)
    logexp = np.log10(exp_f)
    resid = logf - logexp
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid))) or 1e-9
    z = resid / sigma

    rule = cfg.get("select_rule", "zscore")
    z_thr = float(cfg.get("z_thr", 1.5))
    min_fano = float(cfg.get("min_fano", 2.5))
    z_thr_display = float(cfg.get("z_thr_display", z_thr))
    min_fano_display = float(cfg.get("min_fano_display", min_fano))

    if rule != "zscore":
        pass

    strict_mask = (z >= z_thr) & (f >= min_fano)
    loose_mask = (z >= z_thr_display) & (f >= min_fano_display)

    return x, f, names, exp_f, strict_mask, loose_mask, cfg, z, resid_ratio, resid, sigma


def print_fano_summary(
    genes: list[str],
    stats1: pd.DataFrame,
    stats2: pd.DataFrame,
    thr: float = 1.5,
) -> pd.DataFrame:
    """
    Print and return simple per-gene Fano summary table across two replicates.
    """
    rows = []

    for g in genes:
        def row(st, rep):
            if g in st.index:
                m = float(st.at[g, "mean"])
                v = float(st.at[g, "variance"])
                f = float(st.at[g, "fano"])
                return [g, rep, m, v, f, ("noisy" if f > thr else "near-Poisson")]
            return [g, rep, np.nan, np.nan, np.nan, "NA"]

        rows += [row(stats1, "rep1"), row(stats2, "rep2")]

    df = pd.DataFrame(rows, columns=["gene", "rep", "mean", "variance", "fano", "label"])

    with pd.option_context("display.precision", 3):
        print("\n=== Per-gene Fano summary ===")
        print(df)

    return df