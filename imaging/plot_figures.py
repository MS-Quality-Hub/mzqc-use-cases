"""Draw the imaging QC figure from one or more mzQC files.

One row per mzQC file: total ion current (MS:4000228), base peak fraction
(MS:4000230) and reference m/z error (MS:4000232) maps. Everything is read back
from the parsed mzQC, never from in-memory metrics, so the figure shows what the
file actually contains.

Usage
-----
    python plot_figures.py out/A.mzQC out/B.mzQC --out figures/imaging_qc
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Drawn at print size: one double-column figure, 180 mm wide, so point sizes
# below are the point sizes on the page.
sns.set_theme(context="paper", style="ticks")
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 5,
        "ytick.labelsize": 5,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)
FIG_WIDTH_IN = 7.09  # 180 mm
ROW_HEIGHT_IN = 1.85


def load_report(path: Path) -> dict[str, Any]:
    """The label, reference m/z and the three matrix metrics of an mzQC file."""
    run = json.loads(path.read_text())["mzQC"]["runQualities"][0]
    metrics = {m["accession"]: m["value"] for m in run["qualityMetrics"]}
    params = {p["accession"]: p["value"] for p in run["metadata"]["cvParameters"]}

    def as_array(accession: str) -> np.ndarray:
        return np.array(
            [[np.nan if v is None else v for v in row] for row in metrics[accession]],
            dtype=float,
        )

    return {
        "label": run["metadata"]["label"],
        "ref_mz": params["MS:4000235"],
        "tic": as_array("MS:4000228"),
        "bpf": as_array("MS:4000230"),
        "mzerr": as_array("MS:4000232"),
    }


def panels(report: dict[str, Any]) -> list[dict[str, Any]]:
    """How each matrix metric is drawn: colour map and limits."""
    tic, bpf, mzerr = report["tic"], report["bpf"], report["mzerr"]
    lim = np.nanpercentile(np.abs(mzerr), 98)
    return [
        dict(img=tic, title="total ion current", cmap="viridis", cbar_label="TIC",
             vmin=0, vmax=np.nanpercentile(tic, 99)),
        dict(img=bpf, title="base peak fraction", cmap="magma", cbar_label="base peak / TIC",
             vmin=np.nanpercentile(bpf, 1), vmax=np.nanpercentile(bpf, 99)),
        dict(img=mzerr, title=f"reference m/z error\nat {report['ref_mz']}", cmap="coolwarm",
             cbar_label="mass error (ppm)", vmin=-lim, vmax=lim),
    ]


def draw_map(ax: plt.Axes, img: np.ndarray, title: str, cmap: str, cbar_label: str,
             vmin: float, vmax: float) -> None:
    """One matrix metric as a spatial map; null (not acquired) cells stay blank."""
    im = ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set(title=title, xlabel="x (pixel)", ylabel="y (pixel)")
    # A divider keeps the colourbar exactly as tall as an aspect-locked image.
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.08)
    ax.figure.colorbar(im, cax=cax).set_label(cbar_label)


def plot(reports: list[dict[str, Any]], out_stem: Path) -> None:
    """TIC, base peak fraction and m/z error maps, one row per report, as PNG and PDF.

    Row heights follow each pixel grid's aspect ratio so every map is the same width.
    """
    fig, axes = plt.subplots(
        len(reports), 3,
        figsize=(FIG_WIDTH_IN, ROW_HEIGHT_IN * len(reports)),
        gridspec_kw={"height_ratios": [r["tic"].shape[0] / r["tic"].shape[1] for r in reports]},
        squeeze=False,
    )
    for row, report in enumerate(reports):
        for col, panel in enumerate(panels(report)):
            ax = axes[row, col]
            ax.text(-0.28, 1.02, chr(ord("a") + row * 3 + col), transform=ax.transAxes,
                    fontsize=8, fontweight="bold", va="bottom", ha="right")
            draw_map(ax, **panel)
    fig.tight_layout()
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        fig.savefig(out_stem.with_suffix(suffix))
        print(f"wrote {out_stem.with_suffix(suffix)}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mzqc", type=Path, nargs="+", help="mzQC file(s), one figure row each")
    parser.add_argument("--out", type=Path, default=Path("figures/imaging_qc"),
                        help="output path without suffix (default: figures/imaging_qc)")
    args = parser.parse_args()
    plot([load_report(p) for p in args.mzqc], args.out)


if __name__ == "__main__":
    main()
