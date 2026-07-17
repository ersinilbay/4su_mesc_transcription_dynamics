from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import enrichment as enrich_mod
from config import (
    OUT_DIR,
    REP1_QC,
    REP1_UM,
    REP2_QC,
    REP2_UM,
    apply_plot_style,
)
from enrichment import run_pluri_regime_combo_enrichments, write_supplementary_workbook
from gene_sets import (
    GENES_REP1,
    GENES_REP2,
    REP1_GROUP_COLORS,
    REP1_GROUP_MARKERS,
    REP1_PLURI_GROUPS,
    REP2_GROUP_COLORS,
    REP2_GROUP_MARKERS,
    REP2_PLURI_GROUPS,
    SELECTED_LABELS,
)
from io_utils import ensure_outdir
from state_utils import mask_for_state, sync_state
from stats_utils import (
    gene_stats_from_total,
    list_offscale_fano_genes,
    new_fraction_matrix,
    print_fano_summary,
    safe_stats_over_cells,
    total_matrix,
)
from plotting import (
    make_residual_groups_2x2,
    plot_fano_total_residual_outliers,
    plot_fano_total_with_highlight_labels,
    plot_fano_vs_mean,
    plot_global_mean_variance,
    plot_mean_variance,
)


def load_pluri_subset(path_qc: Path, path_um: Path, tag: str) -> sc.AnnData:
    """
    Load QC + UM objects, sync state labels, and return the pluripotent subset.
    """
    if not path_qc.exists():
        raise FileNotFoundError(f"{tag}: missing QC file: {path_qc}")
    if not path_um.exists():
        raise FileNotFoundError(f"{tag}: missing UM file: {path_um}")

    ad_qc = sc.read_h5ad(path_qc)
    ad_um = sc.read_h5ad(path_um)
    sync_state(ad_qc, ad_um)

    m_pluri = mask_for_state(ad_qc, "pluri")
    if m_pluri is None or m_pluri.sum() == 0:
        raise ValueError(f"{tag}: no pluripotent cells found.")

    ad_pluri = ad_qc[m_pluri, :].copy()
    ad_pluri.obs["state"] = "pluri"
    return ad_pluri


def run_one(
    tag: str,
    path_qc: Path,
    path_um: Path,
    highlight_all: list[str],
    label_only: list[str],
    out_dir: Path,
) -> None:
    """
    Reproduce the per-replicate workflow on pluripotent cells only.
    """
    print(f"[INFO] Loading {tag}...")
    ad_qc = load_pluri_subset(path_qc, path_um, tag=tag)

    state_tag = "pluri"
    mask = None  # already subset to pluri

    # Global / pluri-only TOTAL mean-variance
    plot_global_mean_variance(tag, ad_qc, out_dir)

    # NEW layer
    NEW = new_fraction_matrix(ad_qc)
    new_mean, new_var = safe_stats_over_cells(NEW, mask_cells=mask)
    plot_mean_variance(tag, new_mean, new_var, where="NEW", state_tag=state_tag, out_dir=out_dir)
    plot_fano_vs_mean(tag, new_mean, new_var, where="NEW", state_tag=state_tag, out_dir=out_dir)

    # TOTAL layer
    TOT = total_matrix(ad_qc)
    tot_mean, tot_var = safe_stats_over_cells(TOT, mask_cells=mask)
    plot_mean_variance(tag, tot_mean, tot_var, where="TOTAL", state_tag=state_tag, out_dir=out_dir)
    plot_fano_vs_mean(tag, tot_mean, tot_var, where="TOTAL", state_tag=state_tag, out_dir=out_dir)

    # Residual-Fano panels
    if tag == "rep1":
        plot_fano_total_residual_outliers(
            tag,
            ad_qc,
            mask,
            state_tag,
            out_dir,
            overlay_groups=REP1_PLURI_GROUPS,
            overlay_colors=REP1_GROUP_COLORS,
            overlay_markers=REP1_GROUP_MARKERS,
        )
        make_residual_groups_2x2(
            tag,
            ad_qc,
            state_tag,
            out_dir,
            groups=REP1_PLURI_GROUPS,
            colors=REP1_GROUP_COLORS,
            markers=REP1_GROUP_MARKERS,
        )

    elif tag == "rep2":
        plot_fano_total_residual_outliers(
            tag,
            ad_qc,
            mask,
            state_tag,
            out_dir,
            overlay_groups=REP2_PLURI_GROUPS,
            overlay_colors=REP2_GROUP_COLORS,
            overlay_markers=REP2_GROUP_MARKERS,
        )
        make_residual_groups_2x2(
            tag,
            ad_qc,
            state_tag,
            out_dir,
            groups=REP2_PLURI_GROUPS,
            colors=REP2_GROUP_COLORS,
            markers=REP2_GROUP_MARKERS,
        )
        plot_fano_total_residual_outliers(
            tag,
            ad_qc,
            mask,
            state_tag,
            out_dir,
            overlay_groups={"Model-stable (rep2)": GENES_REP2},
            overlay_colors={"Model-stable (rep2)": "#333333"},
            overlay_markers={"Model-stable (rep2)": "o"},
            show_cutoff=True,
            show_strict_only=False,
            dump_tables=False,
        )

        try:
            disp_csv = out_dir / f"{tag}_{state_tag}_residual_fano_display_only.csv"
            disp_tbl = pd.read_csv(disp_csv)
            overlap = sorted(set(GENES_REP2).intersection(set(disp_tbl["gene"])))
            print(
                f"[INFO] {tag}-{state_tag}: {len(overlap)} of {len(GENES_REP2)} pass display cutoff:\n  "
                + ", ".join(overlap)
            )
        except Exception as e:
            print(f"[WARN] Could not summarize GENES_REP2 overlap: {e}")

    else:
        plot_fano_total_residual_outliers(tag, ad_qc, mask, state_tag, out_dir)

    # Model-stable overlay
    plot_fano_total_with_highlight_labels(
        tag,
        ad_qc,
        highlight_all=highlight_all,
        label_only=label_only,
        mask_cells=mask,
        state_tag=state_tag,
        out_dir=out_dir,
    )


def main() -> None:
    np.random.seed(42)
    apply_plot_style()
    ensure_outdir(OUT_DIR)

    # Reset enrichment module globals for a clean run
    enrich_mod.GO_GLOBAL_XMAX = 0.0
    enrich_mod.SUPP_TABLES.clear()

    # Per-replicate workflow
    run_one(
        "rep1",
        REP1_QC,
        REP1_UM,
        highlight_all=GENES_REP1,
        label_only=SELECTED_LABELS.get("Pluripotent_rep1", []),
        out_dir=OUT_DIR,
    )
    run_one(
        "rep2",
        REP2_QC,
        REP2_UM,
        highlight_all=GENES_REP2,
        label_only=SELECTED_LABELS.get("Pluripotent_rep2", []),
        out_dir=OUT_DIR,
    )

    # Combined pluri-regime GO / pathway enrichment
    run_pluri_regime_combo_enrichments(OUT_DIR)

    # Load pluri subsets for simple summary stats / tables only
    ad1_pluri = load_pluri_subset(REP1_QC, REP1_UM, tag="rep1")
    ad2_pluri = load_pluri_subset(REP2_QC, REP2_UM, tag="rep2")

    stats1 = gene_stats_from_total(ad1_pluri)
    stats2 = gene_stats_from_total(ad2_pluri)

    list_offscale_fano_genes(
        ad1_pluri,
        fano_min=10.0,
        mean_min=0.10,
        out_csv=OUT_DIR / "rep1_pluri_offscale_fano.csv",
        tag="rep1",
    )
    list_offscale_fano_genes(
        ad2_pluri,
        fano_min=10.0,
        mean_min=0.10,
        out_csv=OUT_DIR / "rep2_pluri_offscale_fano.csv",
        tag="rep2",
    )

    print_fano_summary(["Tdgf1", "Farsa", "Gpi1"], stats1, stats2)

    write_supplementary_workbook(OUT_DIR)


if __name__ == "__main__":
    main()