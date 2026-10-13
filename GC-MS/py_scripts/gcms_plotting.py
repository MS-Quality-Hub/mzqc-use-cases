#!/usr/bin/env python3
"""
Chromatogram plotting for GC/MS QC results.

Self-contained (matplotlib only; it never reads or analyzes an mzML file) and
usable on its own with the results of gcms_analysis:

    from gcms_analysis import analyze_mzml, load_expected_compounds
    from gcms_plotting import plot_chromatogram

    analysis = analyze_mzml('sample.mzML', load_expected_compounds('compounds.csv'))
    plot_chromatogram(analysis.chromatogram, analysis.smoothed_chromatogram,
                      analysis.noise_levels, analysis.observed_peaks,
                      analysis.matches, 'sample.png')

The peak dicts it consumes are plain mappings ('rt'/'intensity'/'sn'/... , see
gcms_analysis.PeakInfo and MatchResult); the chromatograms are pyOpenMS
MSChromatogram objects, which are only iterated over, never modified.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from pyopenms import MSChromatogram
    from gcms_analysis import MatchResult, PeakInfo

def add_chrom_to_plot(chromatogram: MSChromatogram, color: str = 'black', label: str = 'Chromatogram') -> None:
    """Helper function to add a chromatogram to the plot."""
    rt_data: list[float] = []
    intensity_data: list[float] = []
    for i in range(chromatogram.size()):
        rt_data.append(chromatogram[i].getRT())
        intensity_data.append(chromatogram[i].getIntensity())
    plt.plot(rt_data, intensity_data, color=color, linewidth=1.0, label=label)

def _points_to_data(ax: Axes, fig: Figure, points: float, axis: str = 'y') -> float:
    """Convert a length in points (1/72 inch) to a delta in data coordinates."""
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((0, 1) if axis == 'y' else (1, 0))
    idx = 1 if axis == 'y' else 0
    data_per_pixel = 1.0 / (p1[idx] - p0[idx])
    pixel_per_point = fig.dpi / 72.0
    return abs(points * pixel_per_point * data_per_pixel)

def _place_expected_peak_labels(
    ax: Axes,
    fig: Figure,
    entries: list[tuple[float, float, str, str]],
) -> float | None:
    """
    Add each (x, y_anchor, label, color) in `entries` as a text artist above its
    anchor, then nudge labels upward (using their actual rendered bounding boxes)
    so that labels for closely-eluting or similarly tall peaks don't overlap.
    Using real text extents (rather than a heuristic proximity threshold) means
    only labels that would truly collide get pushed, regardless of how their
    neighbours' individual peak heights compare.

    Returns the highest y (data coords) reached by any placed label's top edge,
    or None if `entries` is empty.
    """
    gap_pt = 8
    pad_pt = 3
    fig.canvas.draw()
    gap = _points_to_data(ax, fig, gap_pt, 'y')

    texts = []
    for x, y, label, color in entries:
        texts.append(ax.text(x, y + gap, label, ha='center', va='bottom',
                              fontsize=8, color=color, zorder=6))
    fig.canvas.draw()

    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    pad_x = _points_to_data(ax, fig, pad_pt, 'x')
    pad_y = _points_to_data(ax, fig, pad_pt, 'y')

    placed: list[tuple[float, float, float, float]] = []  # finalized (x0, x1, y0, y1) boxes, in data coords
    order = sorted(range(len(texts)), key=lambda i: entries[i][0])
    for i in order:
        t = texts[i]
        bbox = t.get_window_extent(renderer=renderer)
        (x0, y0) = inv.transform((bbox.x0, bbox.y0))
        (x1, y1) = inv.transform((bbox.x1, bbox.y1))
        xlo, xhi = min(x0, x1), max(x0, x1)
        ylo, yhi = min(y0, y1), max(y0, y1)

        shift = 0.0
        moved = True
        while moved:
            moved = False
            cur_ylo, cur_yhi = ylo + shift, yhi + shift
            needed = shift
            for (px0, px1, py0, py1) in placed:
                if xlo - pad_x < px1 and xhi + pad_x > px0 and cur_ylo - pad_y < py1 and cur_yhi + pad_y > py0:
                    needed = max(needed, (py1 - ylo) + pad_y)
            if needed > shift:
                shift = needed
                moved = True

        if shift:
            orig_x, orig_y = t.get_position()
            t.set_position((orig_x, orig_y + shift))
        placed.append((xlo, xhi, ylo + shift, yhi + shift))

    return max((yhi for (_, _, _, yhi) in placed), default=None)

EXPECTED_PEAK_COLORS: list[str] = ['tab:red', 'tab:blue', 'tab:green', 'tab:purple']

def plot_chromatogram(
    chromatogram: MSChromatogram,
    smoothed_chrom: MSChromatogram,
    noise_levels: list[float],
    observed_peaks: list[PeakInfo],
    results: list[MatchResult],
    output_filename: str | Path,
    log_y: bool = False,
) -> None:
    """
    Plot the chromatogram and the peaks matched to expected compounds.

    - Chromatogram: black line
    - Smoothed chromatogram: green line
    - Noise level (from the same S/N estimator used for peak picking): grey line
    - Matched peaks: triangles at the peak apex, labeled with compound name,
      delta RT, and S/N, placed above the peak (nudged to avoid overlap between
      neighbouring labels). Each compound gets its own color, cycled across
      EXPECTED_PEAK_COLORS in RT order.
    - Unmatched picked peaks (no corresponding expected compound): small black
      circles at the peak apex, unlabeled.

    If log_y, the y-axis (intensity) is log-scaled; non-positive intensities are
    then simply not drawn (matplotlib masks them), which is expected for GC/MS
    baseline points.
    """
    # Create the plot
    fig = plt.figure(figsize=(12, 6))
    ax = plt.gca()
    if log_y:
        ax.set_yscale('log')

    # Plot chromatogram
    #add_chrom_to_plot(chromatogram, color='black', label='Chromatogram (TIC)')
    add_chrom_to_plot(smoothed_chrom, color='green', label='Smoothed Chromatogram')
    rt_data = [smoothed_chrom[i].getRT() for i in range(smoothed_chrom.size())]
    plt.plot(rt_data, noise_levels, color='grey', linewidth=1.0, label='Noise Level', zorder=1)

    # Plot picked peaks that weren't matched to any expected compound
    matched_rts = {row['observed_rt'] for row in results if row['observed_rt'] is not None}
    unmatched_peaks = [p for p in observed_peaks if p['rt'] not in matched_rts]
    if unmatched_peaks:
        plt.scatter([p['rt'] for p in unmatched_peaks], [p['intensity'] for p in unmatched_peaks],
                    color='black', s=20, marker='o', label='Unmatched Picked Peaks', zorder=4)

    # Plot only the peaks matched to an expected compound, anchored at the peak
    # apex (or the expected RT at baseline if unmatched), colored per compound
    if results:
        y_max = max((p['intensity'] for p in observed_peaks), default=1.0)
        rows = sorted(results, key=lambda r: r['expected_rt'])
        colors = [EXPECTED_PEAK_COLORS[i % len(EXPECTED_PEAK_COLORS)] for i in range(len(rows))]
        anchors = [
            (row['observed_rt'], row['intensity']) if row['observed_rt'] is not None
            else (row['expected_rt'], y_max * 0.02)
            for row in rows
        ]

        for (x, y), color in zip(anchors, colors):
            plt.scatter([x], [y], color=color, s=60, marker='^', zorder=5)

        entries = []
        for row, (x, y), color in zip(rows, anchors, colors):
            if row['delta_rt'] is not None:
                label = f"{row['compound_name']}\nΔRT: {row['delta_rt']:.2f}s, S/N: {row['sn']:.1f}"
            else:
                label = f"{row['compound_name']}\nno match"
            entries.append((x, y, label, color))
        label_top = _place_expected_peak_labels(ax, fig, entries)

        # Labels can get pushed above the tallest peak by collision avoidance;
        # extend the y-axis so they stay inside the plot instead of floating
        # above the frame or getting clipped. On a log axis the margin has to
        # be multiplicative (an additive one derived from ymin makes no sense
        # once ymin is near-zero-but-positive).
        if label_top is not None:
            ymin, ymax = ax.get_ylim()
            if label_top > ymax:
                new_ymax = label_top * 1.15 if log_y else label_top + (label_top - ymin) * 0.03
                ax.set_ylim(ymin, new_ymax)

    # Labels and formatting
    plt.xlabel('Retention Time (seconds)')
    plt.ylabel('Intensity')
    plt.title('Chromatogram with Matched Expected Peaks' + (' (log scale)' if log_y else ''))
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save as PNG
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_filename}")
