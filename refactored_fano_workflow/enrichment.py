from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import (
    ENRICH_STRICT,
    GO_BAR_HEIGHT,
    GO_FDR_MAX,
    GO_FIG_WIDTH_IN,
    GO_LEFT_MARGIN,
    GO_MAX_TERMS,
    GO_ROW_HEIGHT_IN,
    GO_TERM_MAXCHARS,
    PREFERRED_ENRICH_LIBRARIES,
)
from gene_sets import (
    PLURI_REP1_BOTTOMRIGHT,
    PLURI_REP1_TOPLEFT,
    PLURI_REP2_BOTTOMRIGHT,
    PLURI_REP2_TOPLEFT,
)
from io_utils import ensure_outdir, save_current_fig

try:
    import gseapy as gp
    HAS_GSEAPY = True
except Exception:
    HAS_GSEAPY = False

SUPP_TABLES: list[tuple[str, pd.DataFrame]] = []
GO_GLOBAL_XMAX = 0.0


def present_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in ["Term", "n", "Adjusted P-value", "log10p", "Overlap", "Genes"] if c in df.columns]


def enrichr_available_libraries() -> set[str]:
    if not HAS_GSEAPY:
        return set()
    try:
        return set(gp.get_library_name())
    except Exception:
        try:
            return set(gp.get_libraries())
        except Exception:
            return set()


def validate_requested_libraries(requested: tuple[str, ...]) -> tuple[list[str], list[str]]:
    if not HAS_GSEAPY:
        raise RuntimeError("gseapy is not installed. Run: pip install -U gseapy pandas")
    available = enrichr_available_libraries()
    valid = [nm for nm in requested if nm in available]
    missing = [nm for nm in requested if nm not in available]
    return valid, missing


def run_one_enrichr(gene_list: list[str], universe: list[str] | None, library) -> pd.DataFrame:
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=library,
            background=universe,
            outdir=None,
            cutoff=1.0,
            verbose=False,
        )
        if enr.results is None or enr.results.empty:
            return pd.DataFrame()

        df = enr.results.copy()
        for cand in ["Adjusted P-value", "Adjusted P value", "Adj P-value", "Adjusted Pval", "FDR q-value", "FDR"]:
            if cand in df.columns:
                df.rename(columns={cand: "Adjusted P-value"}, inplace=True)
                break

        if "Genes" not in df.columns:
            for cand in ["Gene_set", "Gene Set", "GeneSet", "genes"]:
                if cand in df.columns:
                    df.rename(columns={cand: "Genes"}, inplace=True)
                    break

        if "Overlap" not in df.columns and "Genes" in df.columns:
            def make_overlap_from_genes(s):
                if pd.isna(s):
                    return "0/0"
                n = len([g for g in str(s).replace(";", ",").split(",") if g.strip()])
                return f"{n}/?"
            df["Overlap"] = df["Genes"].map(make_overlap_from_genes)

        keep_cols = [c for c in ["Term", "Adjusted P-value", "Overlap", "Genes"] if c in df.columns]
        return df[keep_cols]
    except Exception as e:
        print(f"[WARN] Enrichr failed for '{library}': {e}")
        return pd.DataFrame()


def standardize_enrich_df_full(df: pd.DataFrame, fdr_max: float) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    d = df.copy()
    if "Adjusted P-value" not in d.columns:
        for cand in ["FDR q-value", "FDR", "P.adjust", "adjp", "q-value", "qvalue"]:
            if cand in d.columns:
                d.rename(columns={cand: "Adjusted P-value"}, inplace=True)
                break
    if "Adjusted P-value" not in d.columns:
        return pd.DataFrame()

    if "Overlap" in d.columns:
        def n_from_overlap(s):
            try:
                a, _ = str(s).split("/")
                return int(a)
            except Exception:
                return np.nan
        d["n"] = d["Overlap"].map(n_from_overlap)
    elif "Genes" in d.columns:
        def n_from_genes(s):
            if pd.isna(s):
                return 0
            return len([g for g in str(s).replace(";", ",").split(",") if g.strip()])
        d["n"] = d["Genes"].map(n_from_genes)
    else:
        d["n"] = np.nan

    d = d[d["Adjusted P-value"] <= fdr_max].copy()
    if d.empty:
        return d

    d["log10p"] = -np.log10(d["Adjusted P-value"].clip(lower=1e-300))
    cols = present_cols(d)
    d = d[cols].sort_values(["log10p", "n"], ascending=[False, False])
    return d


def collect_supp_table(sheet_name: str, df: pd.DataFrame):
    if df is None or df.empty:
        return
    safe = sheet_name.replace("/", "_").replace(" ", "_")
    if len(safe) > 31:
        safe = safe[:31]
    SUPP_TABLES.append((safe, df.copy()))


def max_log10p_from_df(d: pd.DataFrame, fdr_max: float) -> float:
    if d is None or d.empty:
        return 0.0
    d = d.copy()
    if "Adjusted P-value" not in d.columns:
        for cand in ["FDR q-value", "FDR", "P.adjust", "adjp", "q-value", "qvalue"]:
            if cand in d.columns:
                d.rename(columns={cand: "Adjusted P-value"}, inplace=True)
                break
    if "Adjusted P-value" not in d.columns:
        return 0.0
    d = d[d["Adjusted P-value"] <= fdr_max]
    if d.empty:
        return 0.0
    return float((-np.log10(d["Adjusted P-value"].clip(lower=1e-300))).max())


def plot_go_bar(
    df: pd.DataFrame,
    title: str,
    out_path,
    max_terms: int = GO_MAX_TERMS,
    fdr_max: float = GO_FDR_MAX,
    figsize=None,
    xlim=None,
    bar_height: float = GO_BAR_HEIGHT,
    show_xgrid: bool = False,
    box_all_spines: bool = True,
) -> pd.DataFrame:
    if df.empty:
        print(f"[INFO] No enriched terms for: {title}")
        return pd.DataFrame()

    df = df.copy()

    def olap_n_from_overlap(s):
        try:
            a, _ = str(s).split("/")
            return int(a)
        except Exception:
            return np.nan

    n_vals = df["Overlap"].map(olap_n_from_overlap) if "Overlap" in df.columns else None
    if (n_vals is None) or (np.isnan(n_vals).all()):
        def count_genes(s):
            if pd.isna(s):
                return 0
            return len([g for g in str(s).replace(";", ",").split(",") if g.strip()])
        n_vals = df["Genes"].map(count_genes) if "Genes" in df.columns else pd.Series([0] * len(df), index=df.index)
    df["n"] = n_vals

    if "Adjusted P-value" not in df.columns:
        for cand in ["FDR q-value", "FDR", "P.adjust", "adjp", "q-value", "qvalue"]:
            if cand in df.columns:
                df.rename(columns={cand: "Adjusted P-value"}, inplace=True)
                break
    if "Adjusted P-value" not in df.columns:
        print(f"[INFO] No 'Adjusted P-value' column for: {title}")
        return pd.DataFrame()

    df = df[df["Adjusted P-value"] <= fdr_max]
    if df.empty:
        print(f"[INFO] All terms filtered by FDR for: {title}")
        return pd.DataFrame()

    df["log10p"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
    df = df.sort_values(["log10p", "n"], ascending=[False, False]).head(max_terms).iloc[::-1]

    n = len(df)
    if figsize is None:
        fig_w = GO_FIG_WIDTH_IN
        base_h = max(1.6, n * GO_ROW_HEIGHT_IN + 0.8)
        fig_h = 0.8 if n == 1 else base_h
    else:
        fig_w, fig_h = figsize

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    def shorten(s, k=GO_TERM_MAXCHARS):
        s = str(s)
        return s if len(s) <= k else s[:k - 1] + "…"

    ylbls = [f"{shorten(t)} (n={int(nv)})" for t, nv in zip(df["Term"], df["n"].fillna(0))]
    y = np.arange(n)
    bar_h_local = 0.45 if n == 1 else bar_height

    ax.barh(y, df["log10p"], height=bar_h_local)
    ax.set_yticks(y)
    ax.set_yticklabels(ylbls, fontsize=9)
    ax.set_xlabel(r"$-\log_{10}$ (adj $p$)", fontsize=9)
    ax.set_title(title, fontsize=11)

    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        xmax = float(np.nanmax(df["log10p"])) if len(df) else 1.0
        ax.set_xlim(0.0, xmax * 1.05)

    ax.spines["top"].set_visible(box_all_spines)
    ax.spines["right"].set_visible(box_all_spines)
    ax.grid(axis="x", color="0.92", visible=show_xgrid)
    ax.margins(x=0.02)
    plt.subplots_adjust(left=GO_LEFT_MARGIN, right=0.98, top=0.90, bottom=0.12)

    save_current_fig(out_path.parent, out_path.name)
    plt.show()
    return df


def run_pluri_regime_combo_enrichments(out_dir):
    global GO_GLOBAL_XMAX

    if not HAS_GSEAPY:
        print("[INFO] Skipping pluri-regime GO: gseapy not installed.")
        return

    universe = None
    highfreq_lowsize = sorted(set(PLURI_REP1_BOTTOMRIGHT) | set(PLURI_REP2_BOTTOMRIGHT))
    lowfreq_highsize = sorted(set(PLURI_REP1_TOPLEFT) | set(PLURI_REP2_TOPLEFT))

    valid_libs, missing_libs = validate_requested_libraries(tuple(PREFERRED_ENRICH_LIBRARIES))
    if missing_libs:
        msg = "[WARN] Missing Enrichr libraries on this server:\n" + "\n".join(f"  - {m}" for m in missing_libs)
        if ENRICH_STRICT:
            raise RuntimeError(msg + "\nSet ENRICH_STRICT=False to auto-skip missing libraries.")
        print(msg + "\nProceeding with available libraries only.")

    fdr_max = GO_FDR_MAX
    tag = "pluri_combo"

    for lib in valid_libs:
        df_high_raw = run_one_enrichr(highfreq_lowsize, universe, lib)
        df_low_raw = run_one_enrichr(lowfreq_highsize, universe, lib)

        if (df_high_raw is None or df_high_raw.empty) and (df_low_raw is None or df_low_raw.empty):
            continue

        collect_supp_table(f"{tag}_highFreq_lowSize_{lib}", standardize_enrich_df_full(df_high_raw, fdr_max))
        collect_supp_table(f"{tag}_lowFreq_highSize_{lib}", standardize_enrich_df_full(df_low_raw, fdr_max))

        pair_max = max(max_log10p_from_df(df_high_raw, fdr_max), max_log10p_from_df(df_low_raw, fdr_max))
        GO_GLOBAL_XMAX = max(GO_GLOBAL_XMAX, pair_max)

        lib_title = str(lib)

        state_tag = "highFreq_lowSize"
        title = f"{tag} {state_tag}: GO/Pathway — {lib_title}"
        outp = out_dir / f"{tag}_GO_{state_tag}_{lib_title}.svg"
        if df_high_raw is not None and not df_high_raw.empty:
            df_plot = plot_go_bar(
                df_high_raw,
                title,
                outp,
                max_terms=GO_MAX_TERMS,
                fdr_max=fdr_max,
                xlim=(0.0, GO_GLOBAL_XMAX),
                show_xgrid=False,
                box_all_spines=True,
            )
            if isinstance(df_plot, pd.DataFrame) and not df_plot.empty:
                safe_lib = lib_title.replace(" ", "_")
                csv_path = out_dir / f"{tag}_{state_tag}_{safe_lib}_GO_table.csv"
                cols = present_cols(df_plot)
                df_plot[cols].to_csv(csv_path, index=False)
                print(f"[INFO] Saved genes table (top terms): {csv_path}")

        state_tag = "lowFreq_highSize"
        title = f"{tag} {state_tag}: GO/Pathway — {lib_title}"
        outp = out_dir / f"{tag}_GO_{state_tag}_{lib_title}.svg"
        if df_low_raw is not None and not df_low_raw.empty:
            df_plot = plot_go_bar(
                df_low_raw,
                title,
                outp,
                max_terms=GO_MAX_TERMS,
                fdr_max=fdr_max,
                xlim=(0.0, GO_GLOBAL_XMAX),
                show_xgrid=False,
                box_all_spines=True,
            )
            if isinstance(df_plot, pd.DataFrame) and not df_plot.empty:
                safe_lib = lib_title.replace(" ", "_")
                csv_path = out_dir / f"{tag}_{state_tag}_{safe_lib}_GO_table.csv"
                cols = present_cols(df_plot)
                df_plot[cols].to_csv(csv_path, index=False)
                print(f"[INFO] Saved genes table (top terms): {csv_path}")


def write_supplementary_workbook(out_dir):
    if not SUPP_TABLES:
        return
    xlsx_path = out_dir / "GO_supplementary_tables.xlsx"
    with pd.ExcelWriter(xlsx_path) as xw:
        for sheet_name, df in SUPP_TABLES:
            if df is not None and not df.empty:
                df.to_excel(xw, index=False, sheet_name=sheet_name)
    print(f"[INFO] Supplementary workbook saved: {xlsx_path}")