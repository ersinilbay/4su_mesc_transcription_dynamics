from __future__ import annotations

import re
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import (
    BLUE,
    FIGSIZE_RESID_ANN,
    FIGSIZE_SCAT,
    RED,
    RESIDUAL_CFG,
    RESIDUAL_STATE_OVERRIDES,
)
from io_utils import save_current_fig
from stats_utils import (
    expected_curve_eval,
    residual_arrays_and_masks,
    residual_prepare_arrays,
    safe_stats_over_cells,
    total_matrix,
)

# optional label repulsion
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except Exception:
    HAS_ADJUSTTEXT = False


def prep_ax(ax, logx: bool = True, logy: bool = True, grid: bool = True) -> None:
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    if grid:
        ax.yaxis.grid(True, which="major", color="#E0E0E0", alpha=0.6, lw=0.8)
        ax.xaxis.grid(False)
    for sp in ax.spines.values():
        sp.set_alpha(0.9)
        sp.set_linewidth(1.2)


def annotate_points(ax, x, y, names, fontsize: float = 5.8):
    texts = []
    for xi, yi, nm in zip(x, y, names):
        t = ax.text(
            xi,
            yi,
            nm,
            fontsize=fontsize,
            color="black",
            alpha=0.98,
            ha="left",
            va="bottom",
            clip_on=False,
            path_effects=[pe.Stroke(linewidth=1.4, foreground="white", alpha=0.95), pe.Normal()],
            bbox=dict(facecolor="white", alpha=0.55, pad=0.30, edgecolor="none"),
        )
        texts.append(t)

    if HAS_ADJUSTTEXT:
        adjust_text(
            texts,
            x=x,
            y=y,
            ax=ax,
            expand_points=(2.6, 2.8),
            expand_text=(1.5, 1.5),
            force_text=(0.9, 0.9),
            force_points=(0.55, 0.55),
            only_move={"points": "y", "text": "xy"},
            autoalign="y",
            avoid_self=True,
            lim=2800,
            precision=0.001,
        )

    for t, xi, yi in zip(texts, x, y):
        tx, ty = t.get_position()
        p1 = ax.transData.transform((xi, yi))
        p2 = ax.transData.transform((tx, ty))
        dx, dy = (p2 - p1)
        if (dx * dx + dy * dy) ** 0.5 >= 10:
            frac = 0.9
            px = xi + (tx - xi) * frac
            py = yi + (ty - yi) * frac
            ax.plot([xi, px], [yi, py], lw=0.8, alpha=0.9, color="0.25", solid_capstyle="round")

    return texts


def plot_mean_variance(
    tag: str,
    mean: np.ndarray,
    var: np.ndarray,
    where: str,
    state_tag: str,
    out_dir: Path,
) -> None:
    m = np.isfinite(mean) & np.isfinite(var) & (mean > 0) & (var > 0)
    x, y = mean[m], var[m]

    plt.figure(figsize=FIGSIZE_SCAT)
    ax = plt.gca()
    prep_ax(ax, True, True)
    ax.scatter(x, y, s=9, color=BLUE, alpha=0.45, edgecolors="none")
    ax.set_title(f"{tag} {where}: mean–variance ({state_tag})")
    ax.set_xlabel("Mean")
    ax.set_ylabel("Variance")
    plt.tight_layout(pad=0.7)
    save_current_fig(out_dir, f"{tag}_{state_tag}_{where}_mean-variance.svg")
    plt.show()


def plot_fano_vs_mean(
    tag: str,
    mean: np.ndarray,
    var: np.ndarray,
    where: str,
    state_tag: str,
    out_dir: Path,
) -> None:
    m = np.isfinite(mean) & np.isfinite(var) & (mean > 0)
    x = mean[m]
    fano = var[m] / np.maximum(x, 1e-20)

    plt.figure(figsize=FIGSIZE_SCAT)
    ax = plt.gca()
    prep_ax(ax, True, True)
    ax.scatter(x, fano, s=9, color=BLUE, alpha=0.45, edgecolors="none")
    ax.axhline(1.0, ls=":", lw=1.1, color=BLUE, alpha=0.9, label="Poisson Fano = 1")
    ax.set_title(f"{tag} {where}: Fano vs mean ({state_tag})")
    ax.set_xlabel("Mean")
    ax.set_ylabel("Fano = Var/Mean")
    ax.legend(loc="upper left", framealpha=0.92)
    plt.tight_layout(pad=0.7)
    save_current_fig(out_dir, f"{tag}_{state_tag}_{where}_fano-vs-mean.svg")
    plt.show()


def plot_global_mean_variance(tag: str, ad_qc: sc.AnnData, out_dir: Path) -> None:
    state_tag = "all"
    if "state" in ad_qc.obs and ad_qc.obs["state"].nunique() == 1:
        state_tag = str(ad_qc.obs["state"].unique()[0])

    TOT = total_matrix(ad_qc)
    mean, var = safe_stats_over_cells(TOT, mask_cells=None)
    plot_mean_variance(tag, mean, var, where="TOTAL", state_tag=state_tag, out_dir=out_dir)


def plot_fano_total_with_highlight_labels(
    tag: str,
    adata_qc: sc.AnnData,
    highlight_all: list[str],
    label_only: list[str],
    mask_cells: np.ndarray | None,
    state_tag: str,
    out_dir: Path,
) -> None:
    TOT = total_matrix(adata_qc)
    mean, var = safe_stats_over_cells(TOT, mask_cells=mask_cells)

    m = np.isfinite(mean) & np.isfinite(var) & (mean > 0)
    x_all = mean[m]
    f_all = var[m] / np.maximum(mean[m], 1e-20)
    names_all = np.array(adata_qc.var_names)[m]
    idx = {g: i for i, g in enumerate(names_all)}

    present = [g for g in highlight_all if g in idx]
    bx = np.array([x_all[idx[g]] for g in present]) if present else np.array([])
    by = np.array([f_all[idx[g]] for g in present]) if present else np.array([])

    label_present = [g for g in label_only if g in idx]
    extra_names = []
    if state_tag == "pluri" and len(present) > 0:
        scores = np.log10(np.maximum(bx, 1e-12)) + np.log10(np.maximum(by, 1e-12))
        topk = np.argsort(scores)[-3:]
        extra_names = [present[i] for i in topk if present[i] not in label_present]

    all_label_names = label_present + [g for g in extra_names if g not in label_present]
    if len(all_label_names) > 0:
        lx = np.array([x_all[idx[g]] for g in all_label_names])
        ly = np.array([f_all[idx[g]] for g in all_label_names])
    else:
        lx, ly = np.array([]), np.array([])

    plt.figure(figsize=FIGSIZE_RESID_ANN)
    ax = plt.gca()
    prep_ax(ax, True, True)
    ax.scatter(x_all, f_all, s=16, color=BLUE, alpha=0.33, edgecolors="none")
    ax.axhline(1.0, ls=":", lw=1.1, color=BLUE, alpha=0.9, label="Poisson Fano = 1")

    if len(present) > 0:
        ax.scatter(
            bx,
            by,
            s=40,
            color=RED,
            alpha=0.96,
            edgecolors="black",
            linewidths=0.45,
            label="Model-stable list",
        )

    if len(all_label_names) > 0:
        annotate_points(ax, lx, ly, all_label_names, fontsize=8.5 if len(all_label_names) <= 22 else 8.0)

    ax.set_title(f"{tag} TOTAL: Fano vs mean ({state_tag}; model-stable)")
    ax.set_xlabel("Mean (TOTAL)")
    ax.set_ylabel("Fano = Var/Mean")
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98), borderaxespad=0.0)

    plt.tight_layout(pad=0.9)
    save_current_fig(out_dir, f"{tag}_{state_tag}_TOTAL_fano-vs-mean_model-stable.svg")
    plt.show()


def plot_fano_total_residual_outliers(
    tag: str,
    adata_qc: sc.AnnData,
    mask_cells: np.ndarray | None,
    state_tag: str,
    out_dir: Path,
    overlay_groups: dict[str, list[str]] | None = None,
    overlay_colors: dict[str, str] | None = None,
    overlay_markers: dict[str, str] | None = None,
    *,
    show_cutoff: bool = True,
    show_strict_only: bool = False,
    dump_tables: bool = True,
):
    x, f, names, det, cfg = residual_prepare_arrays(
        adata_qc,
        mask_cells,
        state_tag,
        RESIDUAL_CFG,
        RESIDUAL_STATE_OVERRIDES,
    )
    if len(f) == 0:
        print(f"[WARN] {tag}-{state_tag}: no genes pass residual filters.")
        return []

    _, exp_f = expected_curve_eval(x, f, cfg["winsor_q"])
    resid_ratio = f / np.maximum(exp_f, 1e-20)

    logf = np.log10(f)
    logexp = np.log10(exp_f)
    resid = logf - logexp
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid))) or 1e-9

    rule = cfg.get("select_rule", "zscore")
    z_thr = float(cfg.get("z_thr", 1.2))
    min_fano = float(cfg.get("min_fano", 2.0))
    z_thr_display = float(cfg.get("z_thr_display", z_thr))
    min_fano_display = float(cfg.get("min_fano_display", min_fano))
    top_frac_global = cfg.get("top_frac_global", None)
    min_fold = cfg.get("min_fold", None)

    if rule == "zscore":
        z = resid / sigma
        strict_mask = (z >= z_thr) & (f >= min_fano)
        loose_mask = (z >= z_thr_display) & (f >= min_fano_display)
        score = z
    elif rule == "fold":
        thr = float(min_fold if isinstance(min_fold, (int, float)) and np.isfinite(min_fold) else 1.5)
        strict_mask = (resid_ratio >= thr) & (f >= min_fano)
        loose_mask = (resid_ratio >= max(1.3, thr - 0.2)) & (f >= min_fano_display)
        score = resid_ratio
    elif rule == "fano_global":
        qthr = float(np.quantile(f, 1.0 - float(top_frac_global or 0.10)))
        strict_mask = f >= qthr
        loose_mask = f >= np.quantile(f, 1.0 - float(top_frac_global or 0.10) * 1.3)
        score = f
    else:
        raise ValueError(f"Unknown select_rule: {rule}")

    if not np.any(loose_mask):
        print(f"[WARN] {tag}-{state_tag}: no display outliers under current settings.")
        return []

    display_mask = strict_mask if show_strict_only else loose_mask
    f_cut = exp_f * (10.0 ** (z_thr_display * sigma))

    plt.figure(figsize=FIGSIZE_RESID_ANN)
    ax = plt.gca()
    prep_ax(ax, True, True)
    ax.scatter(x, f, s=12, color=BLUE, alpha=0.28, edgecolors="none")

    if show_cutoff:
        order = np.argsort(x)
        ax.plot(x[order], f_cut[order], color="0.35", lw=1.2, ls="--", label="_nolegend_")

    title_extra = ""
    suffix = ""
    do_label = True
    if overlay_groups and len(overlay_groups) == 1:
        gname = next(iter(overlay_groups))
        safe = re.sub(r"[^A-Za-z0-9]+", "_", gname).strip("_")
        title_extra = f" — {gname}"
        suffix = f"_{safe}"
        do_label = True

    x_o, f_o, names_o = x[display_mask], f[display_mask], names[display_mask]
    score_o = score[display_mask]

    order = np.argsort(score_o)[::-1]
    sel_idx_sub = order[: len(order)]
    sel_x, sel_f, sel_n = x_o[sel_idx_sub], f_o[sel_idx_sub], names_o[sel_idx_sub]

    h_groups = {}
    is_group = np.zeros_like(sel_idx_sub, dtype=bool)
    if overlay_groups:
        name_to_i_full = {n.lower(): i for i, n in enumerate(names_o)}
        for grp, glist in overlay_groups.items():
            ids = [name_to_i_full[g.lower()] for g in glist if g.lower() in name_to_i_full]
            if not ids:
                continue
            gx, gy = x_o[ids], f_o[ids]
            h_groups[grp] = ax.scatter(
                gx,
                gy,
                s=58,
                marker=(overlay_markers or {}).get(grp, "o"),
                color=(overlay_colors or {}).get(grp, "#2ca02c"),
                edgecolors="black",
                linewidths=0.55,
                alpha=0.98,
                label=f"{grp} (n={len(ids)})",
            )
            sel_names_set = {n for n in names_o[ids]}
            is_group |= np.array([n in sel_names_set for n in sel_n])

    h_other = ax.scatter(
        sel_x[~is_group],
        sel_f[~is_group],
        s=24,
        color=RED,
        alpha=0.96,
        edgecolors="black",
        linewidths=0.45,
        label=f"Other high-Fano outliers (n={(~is_group).sum()})",
    )

    if do_label:
        overlay_set = set()
        if overlay_groups:
            for glist in overlay_groups.values():
                overlay_set.update(g.lower() for g in glist)

        wl = set(map(str.lower, cfg.get("whitelist", ())))
        allow = np.array([(nm.lower() in overlay_set) or (nm.lower() in wl) for nm in sel_n], dtype=bool)

        if not allow.any():
            k = min(int(cfg.get("n_label", 35)), len(sel_n))
            allow = np.zeros_like(allow, dtype=bool)
            allow[:k] = True

        annotate_points(
            ax,
            sel_x[allow],
            sel_f[allow],
            sel_n[allow],
            fontsize=(5.8 if allow.sum() <= 22 else 5.6),
        )

    ax.set_title(f"{tag} TOTAL: Residual-Fano outliers ({state_tag}){title_extra}")
    ax.set_xlabel("Mean (TOTAL)")
    ax.set_ylabel("Fano = Var/Mean")

    legend_items = list(h_groups.values()) + [h_other]
    ax.legend(
        legend_items,
        [h.get_label() for h in legend_items],
        loc="upper left",
        framealpha=0.92,
        borderaxespad=0.2,
        ncol=1,
    )

    plt.tight_layout(pad=0.9)
    save_current_fig(out_dir, f"{tag}_{state_tag}_TOTAL_fano-vs-mean_residual-outliers{suffix}.svg")
    plt.show()

    if dump_tables:
        z = resid / sigma if rule == "zscore" else score
        disp_mask = (z >= z_thr_display) & (f >= min_fano_display) if rule == "zscore" else display_mask
        core_mask = (z >= z_thr) & (f >= min_fano) if rule == "zscore" else strict_mask

        table = pd.DataFrame(
            {
                "gene": names,
                "mean": x,
                "fano": f,
                "expected_fano": exp_f,
                "z": (resid / sigma) if rule == "zscore" else np.nan,
                "resid_log10": resid,
                "above_display_cut": disp_mask,
                "strict_core": core_mask,
            }
        ).sort_values(["above_display_cut", "strict_core", "z", "fano"], ascending=[False, False, False, False])

        to_print = table.loc[table["above_display_cut"], ["gene", "mean", "fano", "z"]]
        print(
            f"\n=== genes above display cutoff (sorted by z) — display={int(disp_mask.sum())}, strict={int(core_mask.sum())} ==="
        )
        with pd.option_context("display.max_rows", 600, "display.width", 140):
            print(to_print.to_string(index=False))

        out_all = out_dir / f"{tag}_{state_tag}_residual_fano_table.csv"
        out_disp = out_dir / f"{tag}_{state_tag}_residual_fano_display_only.csv"
        table.to_csv(out_all, index=False)
        to_print.to_csv(out_disp, index=False)
        print(f"[INFO] Saved: {out_all}")
        print(f"[INFO] Saved: {out_disp}")

    return list(names[strict_mask])


def make_residual_groups_2x2(
    tag: str,
    adata_qc: sc.AnnData,
    state_tag: str,
    out_dir: Path,
    groups: dict[str, list[str]],
    colors: dict[str, str],
    markers: dict[str, str],
) -> None:
    x, f, names, exp_f, strict_mask, loose_mask, cfg, *_ = residual_arrays_and_masks(
        adata_qc,
        mask_cells=None,
        state_tag=state_tag,
        residual_cfg=RESIDUAL_CFG,
        overrides=RESIDUAL_STATE_OVERRIDES,
    )
    x_all, f_all, names_all = x, f, names
    display_mask = loose_mask

    _, exp_f_full = expected_curve_eval(x_all, f_all, cfg["winsor_q"])
    resid = np.log10(f_all) - np.log10(exp_f_full)
    sigma = 1.4826 * np.median(np.abs(resid - np.median(resid))) or 1e-9
    z_thr_display = float(cfg.get("z_thr_display", cfg.get("z_thr", 1.2)))
    f_cut = exp_f_full * (10.0 ** (z_thr_display * sigma))

    def panel(ax, group_name: str):
        prep_ax(ax, True, True)
        ax.scatter(x_all, f_all, s=12, color=BLUE, alpha=0.28, edgecolors="none")
        order = np.argsort(x_all)
        ax.plot(x_all[order], f_cut[order], color="0.35", lw=1.2, ls="--", label="_nolegend_")

        x_disp = x_all[display_mask]
        f_disp = f_all[display_mask]
        n_disp = names_all[display_mask]
        ax.scatter(
            x_disp,
            f_disp,
            s=18,
            color="lightcoral",
            alpha=0.85,
            edgecolors="black",
            linewidths=0.35,
            label="Other high-Fano outliers",
        )

        want = set(groups.get(group_name, []))
        sel = np.array([g in want for g in n_disp])

        n_in_group = 0
        if sel.any():
            gx, gy, gn = x_disp[sel], f_disp[sel], n_disp[sel]
            ax.scatter(
                gx,
                gy,
                s=62,
                color=colors.get(group_name, "black"),
                marker=markers.get(group_name, "o"),
                edgecolors="black",
                linewidths=0.55,
                alpha=0.98,
            )
            annotate_points(ax, gx, gy, gn, fontsize=5.8 if sel.sum() <= 22 else 5.6)
            n_in_group = int(sel.sum())

        ax.set_title(f"{group_name} (n={n_in_group})", fontsize=10)
        ax.set_xlabel("Mean (TOTAL)")
        ax.set_ylabel("Fano = Var/Mean")

    fig, axs = plt.subplots(2, 2, figsize=(9.2, 7.2))
    order_names = [
        "Pluripotency / early-embryo signaling",
        "Chromatin / genome regulation",
        "RNA biology",
        "Other regulators",
    ]
    for ax, nm in zip(axs.ravel(), order_names):
        panel(ax, nm)

    fig.suptitle(f"{tag} TOTAL: Residual-Fano outliers ({state_tag}) — grouped", y=0.995, fontsize=11)
    plt.tight_layout(pad=1.0)
    outp = out_dir / f"{tag}_{state_tag}_TOTAL_residual_groups_2x2.svg"
    fig.savefig(outp, format="svg", bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)
    print(f"[INFO] Saved: {outp}")