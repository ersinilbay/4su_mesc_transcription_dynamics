from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt


def ensure_outdir(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)


def save_current_fig(out_dir: Path, filename: str) -> None:
    ensure_outdir(out_dir)
    try:
        plt.gcf().canvas.draw()
        plt.savefig(out_dir / filename, format="svg", bbox_inches="tight")
    except Exception as e:
        print(f"[WARN] SVG save failed for {filename}: {e} — trying PNG fallback.")
        plt.savefig(out_dir / Path(filename).with_suffix(".png"), dpi=300)